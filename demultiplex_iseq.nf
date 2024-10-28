#!/usr/bin/env nextflow

/*
 * Perform demultiplexing of an Illumina iSeq run that has UMIs. 
 * Requires bcl2fastq
 * Example usage:
 * nextflow run demultiplex_iseq.nf --run_dir ~/path/{<YYYYMMDD>_<Instrument ID>_<Run Number>_<Flow Cell ID>}/ --output_dir ~/path/output_dir --project_id seq_project_id
 * Author: A. Verhinina
 * Date: 28 Oct 2024
 */

// Stage 1: Use RunInfo xml to find UMI Length. 

process processFindUmiLength {
    input:
    val run_dir

    output:
    path 'umi_length.txt'

    script:
    """
    python -c "
import os
from bs4 import BeautifulSoup
import logging

def find_umi_length(run_dir):
    # Find the run XML settings and calculate UMI length.
    # Barcodes are usually 8bp long,
    # Example -  if the barcode setting is 19, in reality it is 8bp of the barcode + 11bp of the UMI.
    # NextSeq or NovaSeq may have another formatting, 
    # but iSeq always has this formatting because it does not understand UMIs natively.

    with open(f'{run_dir}/RunInfo.xml', 'r') as f:
        data = f.read()

    runinfo_data = BeautifulSoup(data, 'xml') # Use XML reader to parse the file
    index_read_1 = runinfo_data.find('Read', {'IsIndexedRead':'Y', 'Number':'2'})
    index_read_2 = runinfo_data.find('Read', {'IsIndexedRead':'Y', 'Number':'3'})
    i7_len = index_read_1.get('NumCycles')
    i5_len = index_read_2.get('NumCycles')
    umi_len = int(i7_len) - int(i5_len)

    logging.info(f'UMI length is {umi_len}')

    return umi_len

umi_length = find_umi_length('${run_dir}')
with open('umi_length.txt', 'w') as f:
    f.write(str(umi_length))
    "
    """
}

// Stage 2: Duplicate SampleSheet with original settings for record keeping.

process processCopySampleSheet {
    input:
    val run_dir

    output:
    path 'SampleSheet_original.csv'

    script:
    """
    python -c "
import shutil
import logging
import os

logging.basicConfig(level=logging.INFO)

def copy_sample_sheet(run_dir):
    try:
        src = os.path.join(run_dir, 'SampleSheet.csv')
        dst = 'SampleSheet_original.csv'
        shutil.copy(src, dst)
        logging.info('Saved the original SampleSheet for the record.')

    except Exception as e:
        logging.error(f'Error copying SampleSheet: {e}')
        raise

copy_sample_sheet('${run_dir}')
    "
    """
}

// Stage 3: Modify SampleSheet to strip Ns. 

// Because iSeq does not natively parse UMIs, users use Ns to extend indexing cycles and trick iSeq into reading UMIs.
// These Ns are used to set Illumina to sequence the UMIs, thus UMI length will be equal to count of Ns.
// In the sample sheet the barcodes will look like ATGCATGCNNNNNNNNNNN, where ATGCATGC is i7. 
// Bcl2fastq does not understand Ns, however, so we will strip them from the SampleSheet. 
// Thanks Illumina for this headache.

process processModifySampleSheet {
    input:
    val run_dir
    path 'umi_length.txt'

    output:
    path 'SampleSheet.csv'

    script:
    """
    echo "The value of run_dir is: ${run_dir}"

    python -c "
import os
import shutil
import csv
import logging

logging.basicConfig(level=logging.INFO)

def strip_trailing_n(value, n):
    #Strips a set count of Ns from a string, if this string ends with Ns.
    
    if isinstance(value, str) and value.endswith('N' * n):
        return value[:-n]
    return value

def modify_sample_sheet(original_run_dir, umi_len):
    logging.info(f'Modified the original SampleSheet.csv and moved it back to the {original_run_dir}.')
    
    # The code below takes the SampleSheet csv, 
    # skips the header, 
    # goes to the data section, 
    # strips Ns from barcode sequence column.
    
    original_file = os.path.join(original_run_dir, 'SampleSheet.csv')
    working_copy = 'SampleSheet.csv'

    shutil.copy(original_file, working_copy)

    temp_file = 'Temp_SampleSheet.csv'

    with open(working_copy, 'r') as infile, open(temp_file, 'w', newline='') as outfile:

        reader = csv.reader(infile)
        writer = csv.writer(outfile)

        inside_data_section = False
        index_i7_col = None

        for row in reader:
            # Write the file header unchanged
            if not inside_data_section:
                writer.writerow(row)

                # accomodate iSeq and NextSeq samplesheet formatting
                if '[Data]' in row or '[BCLConvert_Data]' in row:
                    inside_data_section = True
            else:
                if index_i7_col is None and ('index' in row or 'Index' in row):
                    index_i7_col = row.index('index') if 'index' in row else row.index('Index')
                elif index_i7_col is not None:
                    row[index_i7_col] = strip_trailing_n(row[index_i7_col], umi_len)
                writer.writerow(row)

    os.replace(temp_file, working_copy)
    shutil.copy(working_copy, original_file)
    
    return None

    
# get the UMI count
with open('umi_length.txt', 'r') as f:
    umi_len = int(f.read().strip())

logging.info(f'Umi length is {umi_len}.')

# strip the UMIs from the barcodes
modify_sample_sheet('${run_dir}', umi_len)
    "
    """
}

// Stage 4: Run bcl2fastq

// Conversion is standard, except for a few parts.
// 1) Please note that we use base mark setting to allow UMI reading.
// To make sure UMIs are not masked, we need to use --mask-short-adapter-reads 0
// 2) This process creates 5 files. 2 for barcodes (I1-I2), two for reads (R1-R3) and one for UMIs (R2)
// This is the only way how bcl2fastq can write UMIs out of the box
// Alternative would be to use bcl-convert, but standard iSeqs are not formatted to understand bcl-convert SampleSheet

// If bcl2fastq kills itself, you can set threads using native flags

process processBcl2Fastq {

    publishDir "${output_dir}", mode: 'copy'

    cpus params.cpu_count

    input:
    path run_dir
    val output_dir

    output:
    path "demux_output", emit: demux_fastq

    // container 'bcl2fastq:latest' // Potentiall can use a Docker container to run bcl2fastq. Comment it out to use.

    // Please note that 8bp lenth of a barcode is hardcoded below.
    // If non-8bp barcode length is used, use-bases-mask flags will need to be adjusted.

    script:
    """
    bcl2fastq --runfolder-dir ${run_dir} --mask-short-adapter-reads 0 -r ${task.cpus} -p ${task.cpus} -w ${task.cpus} \
              --output-dir demux_output \
              --use-bases-mask Y*,I8Y*,I8,Y* \
              --ignore-missing-bcls \
              --ignore-missing-filter \
              --ignore-missing-positions \
              --ignore-missing-controls \
              --create-fastq-for-index-reads
    """
}

// Stage 5: Optional - concatenate UMI with i7 barcode using seqkit

// Since we created 5 files instead of 4 in the stage above, we can re-generate 4 standard files I1-I2, R1-R2.
// A workaround that allows 4 files instead of 5 is i7 barcode concatendated with UMI sequence.
// Here we use seqkit, a custom script, and FASTQ headers to achieve this.
// This is a messy workaround, so you can skip if you want to implement your own solution.
// seqkit can be installed with conda and should be located in the same environment as bcl2fastq


process processSeqkit {
    input:
    val project_id
    path output_dir
    path "umi_concat.py"  // Ensure this script is available for execution

    
    output:
    path "seqkit_output/*", emit: seqkit_fastq 


    publishDir "${output_dir}", mode: 'copy'

    script:
    """
    python ./umi_concat.py --in_dir ${output_dir}/demux_output/${project_id} --out_dir seqkit_output
    """

}
// Define the input channels for concatSamples
workflow {

    // Define the params, which you can set in the command line or nextflow.config file
    // For example you can set them right here in the script:
    // params.run_dir = '/path/run_name'
    // params.output_dir = 'path/demux_outdir'
    // params.project_id = "project_001" - this will be the project name as used on the sequencer

    // Use params to define key directories and the specific project_id
    def run_dir = file(params.run_dir)
    def output_dir = file(params.output_dir)
    def project_id = params.project_id

    params.cpu_count = 8 // change here if you need to adjust threading

    // Path to the Python script in the current directory
    umi_concat_script = file("umi_concat.py")
    
    // First process to find the UMI length
    umi_length = processFindUmiLength(params.run_dir)

    // Second process to copy the SampleSheet (it only needs run_dir)
    processCopySampleSheet(params.run_dir)

    // Third process to modify the SampleSheet, it takes both run_dir and the UMI length
    processModifySampleSheet(params.run_dir, umi_length)

    // Fourth process to run bcl2fastq
    processBcl2Fastq(params.run_dir, params.output_dir)

    // Call concat script
    // processSeqkit(params.project_id, params.output_dir, umi_concat_script)


}


# Concatenate sequence of UMI with i7
# Requires seqkit (can be installed with conda)
# See argparse below for oprions
# Author: A. Vershinina
# Date; 28 Oct 2024

import os
import glob
import subprocess
import logging
import shutil
import argparse
import gzip
import tempfile

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(funcName)s - %(message)s',
                    level=logging.INFO)

def find_samplenames(in_dir) ->list:
    """ find samples"""
    onlyfiles = [f for f in os.listdir(in_dir) if os.path.isfile(os.path.join(in_dir, f))]
    # Extracting unique sample names
    sample_names = set()

    for filename in onlyfiles:
        # Split the filename by underscore and take the first part
        sample_name = '_'.join(filename.split('_')[0:-3])
        sample_names.add(sample_name)

    # Convert the set to a list and print the unique sample names
    unique_sample_names = list(sample_names)

    return unique_sample_names

def concat_with_seqkit(i1_path: str, r2_path: str, output_path: str):
    """ seqkit command """

    logging.info(f"Concatenating {i1_path} and {r2_path} into {output_path}")

    try:
        # Command to concatenate i1 and r2 using seqkit
        cmd = [
            "seqkit", "concat", 
            i1_path, r2_path,
            "-o", output_path
        ]

        # Run the command and capture the output
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check for errors
        if result.returncode != 0:
            logging.error("Error: %s", result.stderr)
        else:
            logging.info("Concatenated sequences saved to %s", output_path)

    except Exception as e:
        logging.error(f"An error occurred: {e}")
    
    return None

def correct_read_id(fq_path: str):
    """ 
    Since we rename R3 into R2, we need to fix the counter in fastq headers.
    This function replaces the read ID from 3: to 2:

    Additionally, since we attach UMI to barcode, we modify the read header.
    It will look like "1:N:0:BARCODE+BARCODE|2:N:0:BARCODE+BARCODE".
    This function returns the header back to 1:N:0:BARCODE+BARCODE.
    """
    basename = os.path.basename(fq_path)

    # Create a temporary directory to store the modified content
    with tempfile.TemporaryDirectory() as tmpdirname:
        temp_path = os.path.join(tmpdirname, basename)
    
        with gzip.open(fq_path, 'rt') as infile, gzip.open(temp_path, 'wt') as outfile:
            for i, line in enumerate(infile):
                if i % 4 == 0:  # This is a header line (every 4th line)
                    if '_R2_' in basename:
                        # Replace the first occurrence of ' 3:' with ' 2:' in the read ID for R2 files
                        line = line.replace(' 3:', ' 2:', 1)
                    
                    elif '_I1_' in basename:
                        # Strip the read ID after the pipe for index files (I1)
                        line = line.split('|')[0] + '\n'

                # Write the (possibly modified) line to the output
                outfile.write(line)

        # Overwrite the original file with the modified temp file
        shutil.move(temp_path, fq_path)

    return None




# def replace_fq_filename(filename: str) -> str:
#     """
#     Example function to replace parts of the filename.
#     Modify based on the actual renaming rules you need.
#     """
#     return filename.replace("_R3_", "_renamed_R3_")


def concat_all_samples(fq_dir: str, out_dir: str):
    """ Running concatenation and handling outputput files"""
    logging.info(f'Working in dir: {fq_dir}')

    # Find sample names
    samples = find_samplenames(fq_dir)
    logging.info(f'Samples: {samples}')

    # Check if the output directory exists
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        logging.info(f"Directory '{out_dir}' created successfully.")
    else:
        logging.info(f"Directory '{out_dir}' already exists.")

    for s in samples:
        if 'Undetermined' in s:
            logging.info(f"Skipping sample: {s}")
            continue

        logging.info(f"Processing sample: {s}")

        # Glob paths for i1, r1, r2, and r3 files
        try:
            i1_paths = glob.glob(os.path.join(fq_dir, f"{s}*I1*.fastq.gz"))  # Path to i1 file
            i2_paths = glob.glob(os.path.join(fq_dir, f"{s}*I2*.fastq.gz"))  # Path to i2 file
            r1_paths = glob.glob(os.path.join(fq_dir, f"{s}*R1*.fastq.gz"))  # Path to r1 file
            r3_paths = glob.glob(os.path.join(fq_dir, f"{s}*R3*.fastq.gz"))  # Path to r3 file (to rename)

            if not i1_paths or not i2_paths or not r1_paths or not r3_paths:
                logging.error(f"Missing I1, I2, R1, or R3 files for sample {s}. Skipping.")
                continue

            i1_path = i1_paths[0]
            i2_path = i2_paths[0]
            r1_path = r1_paths[0]
            r3_path = r3_paths[0]

        except IndexError as e:
            logging.error(f"Error finding files for sample {s}: {e}")
            continue

        # Log the paths
        for p in [i1_path, i2_path, r1_path, r3_path]:
            logging.info(f"Found path: {p}")

        # Output path for concatenated i1 + r2 file
        concat_i1 = os.path.join(out_dir, f"{s}_L001_I1_001.fastq.gz")
        logging.info(f"Output concatenated file path: {concat_i1}")

        # Concatenate i1 and r2
        r2_paths = glob.glob(os.path.join(fq_dir, f"{s}*R2*.fastq.gz"))  # Path to r2 file
        if not r2_paths:
            logging.error(f"Missing R2 file for sample {s}. Skipping.")
            continue
        r2_path = r2_paths[0]

        # Concatenate i1 and r2 using seqkit or a similar tool
        concat_with_seqkit(i1_path, r2_path, concat_i1)

        # Rename r3 to r2 in the output directory
        renamed_r3_path = os.path.join(out_dir, f"{s}_L001_R2_001.fastq.gz")

        # Copy I2, R1, and renamed R3 (as R2) to the output directory
        try:
            shutil.copy(i2_path, os.path.join(out_dir, os.path.basename(i2_path)))  # Copy I2 unchanged
            shutil.copy(r1_path, os.path.join(out_dir, os.path.basename(r1_path)))  # Copy R1 unchanged
            shutil.copy(r3_path, renamed_r3_path)  # Copy R3 as R2

            for f in [renamed_r3_path, concat_i1]:
                correct_read_id(f) # Correct fq header after renaming, correct I1 header after concatenation

        except Exception as e:
            logging.error(f"Error during file operations for sample {s}: {e}")
            continue  # Proceed to the next sample
    
    return None

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Concatenate UMI with I1 sequence. \
                                     Usage: python umi_concat.py --in_dir demux_output --out_dir test_out")
    
    # Define the command-line arguments
    parser.add_argument(
        '--in_dir',
        type=str,
        required=True,
        help="Input directory containing the fastq files: I1,I2, R1,R2,R3. R2 is the file with UMI sequence"
    )
    parser.add_argument(
        '--out_dir',
        type=str,
        required=True,
        help="Output directory to store the intermediate files. These can be removed after the concatenation."
    )

    # Parse the arguments
    args = parser.parse_args()

    # Extract arguments
    in_dir = args.in_dir
    out_dir = args.out_dir

    # Ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)

    # Call the function to concatenate the samples
    concat_all_samples(in_dir, out_dir)
    

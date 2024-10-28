# demultiplex_iseq_w_umi
Nextflow pipeline to convert bcl to fastq and accomodate UMIs.

* Before running, you can use requirement files to install necessary packages.
* bcl2fastq is the main one that is required.
  
# Example usage:
```{bash}
nextflow run demultiplex_iseq.nf --run_dir ~/path/{<YYYYMMDD>_<Instrument ID>_<Run Number>_<Flow Cell ID>}/ --output_dir ~/path/output_dir --project_id seq_project_id
```

# Workflow

Author: A. Verhinina
Date: 28 Oct 2024

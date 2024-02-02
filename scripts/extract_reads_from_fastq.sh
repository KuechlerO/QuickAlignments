#!/bin/bash
 
# Input parameters
fastq_file1=${snakemake_input[0]}
fastq_file2=${snakemake_input[1]}
read_ids_file=${snakemake_input[2]}

output_fastq_file1=${snakemake_output[0]}
output_fastq_file1=${snakemake_output[1]}
 
# Extract reads
seqtk subseq $fastq_file1 $read_ids_file > $output_fastq_file1
seqtk subseq $fastq_file2 $read_ids_file > $output_fastq_file2

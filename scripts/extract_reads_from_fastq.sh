#!/bin/bash
 
# Input parameters
fastq_file1=$1
fastq_file2=$2
read_ids_file=$3

output_fastq_file1=$4
output_fastq_file1=$5
 
# Extract reads
seqtk subseq $fastq_file1 $read_ids_file > $output_fastq_file1
seqtk subseq $fastq_file2 $read_ids_file > $output_fastq_file2

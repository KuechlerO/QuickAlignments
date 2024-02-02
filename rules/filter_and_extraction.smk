# Filter from a BAM-file all reads of interest
rule filter_reads_from_bam_file:
    input:
        input_bam_file = 'input/bam_files/{sample}.sorted.bam',
    output:
        text_report_file="output/filtering/{sample}/read_analysis.txt",
        filtered_reads_ids_gene_region_file="output/filtering/{sample}/reads_in_gene_region.txt",
        exon3_clipped_reads_ids_file="output/filtering/{sample}/exon3_clipped_reads.txt",
        exon2_clipped_reads_ids_file="output/filtering/{sample}/exon2_clipped_reads.txt",
        exon3_exon2_clipped_reads_ids_file="output/filtering/{sample}/exon3_exon2_clipped_reads.txt",
    conda:
        "../envs/all.yaml",
    script:
        "../scripts/filter_reads_from_bam_file.py" 


rule extract_reads_from_fastq_file:
    input:
        input_fastq_file_read1 = 'input/fastq_files/{sample}_1.fastq',
        input_fastq_file_read2 = 'input/fastq_files/{sample}_2.fastq',
        reads_of_interest = 'output/filtering/{sample}/{region}.txt',
    output:
        extracted_reads_file_read1="output/filtering/{sample}/extracted_reads/{sample}_1.{region}.fastq",
        extracted_reads_file_read2="output/filtering/{sample}/extracted_reads/{sample}_2.{region}.fastq",
    conda:
        "../envs/all.yaml",
    shell:
        """
        fastq_file1={input.input_fastq_file_read1}
        fastq_file2={input.input_fastq_file_read2}
        read_ids_file={input.reads_of_interest}

        output_fastq_file1={output.extracted_reads_file_read1}
        output_fastq_file1={output.extracted_reads_file_read2}
        
        # Extract reads
        seqtk subseq $fastq_file1 $read_ids_file > $output_fastq_file1
        seqtk subseq $fastq_file2 $read_ids_file > $output_fastq_file2
        """



# Filter from a BAM-file all reads of interest
rule analyse_novel_mappings:
    input:
        input_bam_file = 'output/{mapper}/{sample}.{region}.sorted.bam',
        input_bai_file = "output/{mapper}/{sample}.{region}.sorted.bam.bai",
    output:
        text_report_file="output/analysis/{mapper}/{sample}.{region}/read_analysis.txt",
        csv_report_file="output/analysis/{mapper}/{sample}.{region}/read_analysis.csv",
        filtered_reads_ids_gene_region_file="output/analysis/{mapper}/{sample}.{region}/reads_in_gene_region.txt",
        exon3_clipped_reads_ids_file="output/analysis/{mapper}/{sample}.{region}/exon3_clipped_reads.txt",
        exon2_clipped_reads_ids_file="output/analysis/{mapper}/{sample}.{region}/exon2_clipped_reads.txt",
        exon3_exon2_clipped_reads_ids_file="output/analysis/{mapper}/{sample}.{region}/exon3_exon2_clipped_reads.txt",
        exon3_reads_not_mapped_file="output/extracted_reads/{mapper}/{sample}.{region}/exon3_not_mapped_reads.txt",
        exon2_reads_not_mapped_file="output/extracted_reads/{mapper}/{sample}.{region}/exon2_not_mapped_reads.txt"
    conda:
        "../envs/all.yaml",
    script:
        "../scripts/filter_reads_from_bam_file.py"

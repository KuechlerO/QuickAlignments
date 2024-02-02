configfile: "config.yaml"

# -------------- Load rules ----------------
include:    "rules/filter_and_extraction.smk"
include:    "rules/hisat2.smk"
include:    "rules/samtools.smk"
include:    "rules/star.smk"


rule all:
    input:
        # Star alignments
        expand("output/star/pe/{sample}.{region}.bam",
            sample=config["input_samples"], region=["reads_in_gene_region", "exon3_clipped_reads",
                "exon2_clipped_reads", "exon3_exon2_clipped_reads"]),
        # Hisat2 alignments
        expand("output/hisat2/mapped/{sample}.{region}.bam",
            sample=config["input_samples"], region=["reads_in_gene_region", "exon3_clipped_reads", 
                "exon2_clipped_reads", "exon3_exon2_clipped_reads"]),
    default_target: True    # Makes this rule the default rule

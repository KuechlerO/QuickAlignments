configfile: "config.yaml"

# -------------- Load rules ----------------
include:    "rules/filter_and_extraction.smk"
include:    "rules/hisat2.smk"
include:    "rules/samtools.smk"
include:    "rules/star.smk"
include:    "rules/analyze_mappings.smk"


rule merge_all_csv_outputs:
    input:
        expand("output/extracted_reads/{sample}/read_analysis.csv", sample=config["input_samples"]),
        expand("output/analysis/{tool}/{sample}.{region}/read_analysis.csv",
            tool=["hisat2/mapped", "star/pe", "star_snakesplice/pe", "star_snakesplice2/pe", "star_snakesplice3/pe",
                "star_snakesplice4/pe", "star_snakesplice5/pe", "star_snakesplice6/pe", "star_snakesplice7/pe", "star_snakesplice8/pe",
                "star_snakesplice9/pe", "star_snakesplice10/pe", "star_snakesplice11/pe", "star_snakesplice12/pe",
                "star_snakesplice13/pe"],
            sample=config["input_samples"], 
            region=["reads_in_gene_region", "exon3_clipped_reads", 
                "exon2_clipped_reads", "exon3_exon2_clipped_reads"]),
    output:
        summary_csv_file = "output/summary/summary_table.csv"
    run:
        import pandas as pd

        df = pd.concat([pd.read_csv(f) for f in input ])
        df.to_csv(output.summary_csv_file, index=False)


rule all:
    input:
        # Extract soft clipping info for original data
        expand("output/extracted_reads/{sample}/read_analysis.csv", sample=config["input_samples"]),
        # Extract reads
        expand("output/extracted_reads/{sample}/{sample}_R1.{region}.fastq",
            sample=config["input_samples"], region=["reads_in_gene_region", "exon3_clipped_reads",
                    "exon2_clipped_reads", "exon3_exon2_clipped_reads"]),

        # Hisat2 alignments
        expand("output/hisat2/mapped/{sample}.{region}.sorted.bam.bai",
            sample=config["input_samples"], region=["reads_in_gene_region", "exon3_clipped_reads", 
                "exon2_clipped_reads", "exon3_exon2_clipped_reads"]),
        
        # Star alignments
        expand("output/star/pe/{sample}.{region}.sorted.bam.bai",
            sample=config["input_samples"], region=["reads_in_gene_region", "exon3_clipped_reads",
                "exon2_clipped_reads", "exon3_exon2_clipped_reads"]),
        
        # Analyze new mappings
        expand("output/analysis/{tool}/{sample}.{region}/read_analysis.csv",
            tool=["hisat2/mapped", "star/pe"], sample=config["input_samples"], 
            region=["reads_in_gene_region", "exon3_clipped_reads", 
                "exon2_clipped_reads", "exon3_exon2_clipped_reads"]),

        
        # Star SnakeSplice alignments
        expand("output/{tool}/{sample}.{region}.sorted.bam.bai",
            tool=["star_snakesplice/pe", "star_snakesplice2/pe", "star_snakesplice3/pe",
                "star_snakesplice4/pe", "star_snakesplice5/pe", "star_snakesplice6/pe", "star_snakesplice7/pe", "star_snakesplice8/pe",
                "star_snakesplice9/pe", "star_snakesplice10/pe", "star_snakesplice11/pe", "star_snakesplice12/pe", 
                "star_snakesplice13/pe"], 
            sample=config["input_samples"], region=["reads_in_gene_region", "exon3_clipped_reads",
                "exon2_clipped_reads", "exon3_exon2_clipped_reads"]),
        
        # Summary file
        rules.merge_all_csv_outputs.output,

    default_target: True    # Makes this rule the default rule

rule hisat2_index:
    input:
        fasta = config["reference_genome"],
    output:
        directory("output/hisat2/index")
    params:
        prefix = "output/hisat2/index/ref"
    log:
        "logs/hisat2_index.log"
    threads: 2
    wrapper:
        "v3.12.1/bio/hisat2/index"


rule hisat2_align:
    input:
        reads=["output/extracted_reads/{sample}/{sample}_R1.{region}.fastq",
            "output/extracted_reads/{sample}/{sample}_R2.{region}.fastq"],
        idx="output/hisat2/index/",
    output:
        "output/hisat2/mapped/{sample}.{region}.bam",
    log:
        "logs/hisat2_align_{sample}_{region}.log",
    params:
        extra="",
    threads: 2
    wrapper:
        "v3.12.1/bio/hisat2/align"

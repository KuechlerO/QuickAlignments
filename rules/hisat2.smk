rule hisat2_index:
    input:
        fasta=config["reference_genome"],
    output:
        directory("output/hisat2/hisat2_index"),
    params:
        prefix="output/hisat2/hisat2_index/ref"
    log:
        "logs/hisat2_index.log"
    threads: 2
    wrapper:
        "v3.3.6/bio/hisat2/index"
 
 
rule hisat2_align:
    input:
        reads=["output/filtering/{sample}/extracted_reads/{sample}_1.{region}.fastq",
         "output/filtering/{sample}/extracted_reads/{sample}_2.{region}.fastq"],
        idx="output/hisat2/hisat2_index/",
    output:
        "output/hisat2/mapped/{sample}.{region}.bam",
    log:
        "logs/hisat2_align_{sample}.{region}.log",
    params:
        extra="",
    threads: 2
    wrapper:
        "v3.3.6/bio/hisat2/align"

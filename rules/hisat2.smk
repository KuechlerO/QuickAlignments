rule hisat2_index:
    input:
        fasta="Homo_sapiens.GRCh37.dna.chromosome.2.fa"
    output:
        directory("hisat2_index_hg19_chr2"),
    params:
        prefix="hisat2_index_hg19_chr2/ref"
    log:
        "logs/hisat2_index.log"
    threads: 2
    wrapper:
        "v3.3.6/bio/hisat2/index"
 
 
rule hisat2_align:
    input:
        reads=["{sample}_1.fastq.gz_filtered.fastq", "{sample}_2.fastq.gz_filtered.fastq"],
        idx="hisat2_index_hg19_chr2/",
    output:
        "output/hisat2/mapped/{sample}.bam",
    log:
        "logs/hisat2_align_{sample}.log",
    params:
        extra="",
    threads: 2
    wrapper:
        "v3.3.6/bio/hisat2/align"

rule samtools_sort:
    input:
        "output/star/pe/{sample}.{region}.bam",
    output:
        "output/star/pe/{sample}.{region}.sorted.bam",
    log:
        "logs/sorting_{sample}.{region}.log",
    params:
        extra="-m 4G",
    threads: 8
    wrapper:
        "v3.3.6/bio/samtools/sort"
 
 
rule samtools_index:
    input:
       "output/star/pe/{sample}.{region}.sorted.bam",
    output:
        "output/star/pe/{sample}.{region}.sorted.bam.bai",
    log:
        "logs/samtools_index/{sample}.{region}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v3.3.6/bio/samtools/index"

rule samtools_sort:
    input:
        "output/{mapper}/{sample}.{region}.bam",
    output:
        "output/{mapper}/{sample}.{region}.sorted.bam",
    log:
        "logs/sorting_{mapper}.{sample}.{region}.log",
    params:
        extra="-m 4G",
    threads: 8
    wrapper:
        "v3.3.6/bio/samtools/sort"
 
 
rule samtools_index:
    input:
       "output/{mapper}/{sample}.{region}.sorted.bam",
    output:
        "output/{mapper}/{sample}.{region}.sorted.bam.bai",
    log:
        "logs/samtools_index/{mapper}.{sample}.{region}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v3.3.6/bio/samtools/index"

rule star_index:
    input:
        fasta=config["reference_genome"],
    output:
        directory("output/star/star_index"),
    message:
        "Testing STAR index"
    threads: 1
    params:
        extra="",
    log:
        "logs/star_index.log",
    wrapper:
        "v3.3.6/bio/star/index"
 
 
rule star_pe_multi:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=["output/filtering/{sample}/extracted_reads/{sample}_1.{region}.fastq"],
        # paired end reads needs to be ordered so each item in the two lists match
        fq2=["output/filtering/{sample}/extracted_reads/{sample}_2.{region}.fastq"],  #optional
        # path to STAR reference genome index
        idx="output/star/star_index",
    output:
        # see STAR manual for additional output files
        aln="output/star/pe/{sample}.{region}.bam",
        log="logs/pe/{sample}.{region}/Log.out",
        sj="output/star/pe/{sample}.{region}/SJ.out.tab",
        unmapped=["output/star/pe/{sample}.{region}/unmapped.1.fastq.gz",
            "output/star/pe/{sample}.{region}/unmapped.2.fastq.gz"],
    log:
        "logs/pe/{sample}.{region}.log",
    params:
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate",
    threads: 8
    wrapper:
        "v3.3.6/bio/star/align"

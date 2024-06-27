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
        fq1=["output/extracted_reads/{sample}/{sample}_R1.{region}.fastq"],
        # paired end reads needs to be ordered so each item in the two lists match
        fq2=["output/extracted_reads/{sample}/{sample}_R2.{region}.fastq"],  #optional
        # path to STAR reference genome index
        idx="output/star/star_index",
    output:
        # see STAR manual for additional output files
        aln="output/star/pe/{sample}.{region}.bam",
        log="logs/star/pe/{sample}.{region}/Log.out",
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


rule star_snakesplice_index:
    input:
        # Reference genome
        fasta=config["reference_genome"],
        gtf=config["reference_genome_annotation"],
    output:
        directory("output/star_snakesplice/star_index"),
    threads:
        1
    params:
        sjdbOverhang="99",
        extra=""
    log:
        "logs/star_index_genome.log"
    wrapper:
        "v3.3.3/bio/star/index"


rule star_snakesplice_align:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=["output/extracted_reads/{sample}/{sample}_R1.{region}.fastq"],
        fq2=["output/extracted_reads/{sample}/{sample}_R2.{region}.fastq"],
        idx=rules.star_snakesplice_index.output,
        # GTF file is passed in extra
        gtf_file=config["reference_genome_annotation"],
    output:
        # see STAR manual for additional output files
        aln=os.path.join("output/star_snakesplice/pe/{sample}.{region}.bam"),		# output only BAM, not SAM -> Use samtools view to convert to SAM
        log="logs/star_snakesplice/pe/{sample}.{region}/Log.out",
        sj="output/star_snakesplice/pe/{sample}.{region}/SJ.out.tab",
        unmapped=["output/star_snakesplice/pe/{sample}.{region}/unmapped.1.fastq.gz",
            "output/star_snakesplice/pe/{sample}.{region}/unmapped.2.fastq.gz"],
    log:
        "logs/star_snakesplice/{sample}.{region}_alignment.log"
    params:
        # path to STAR reference genome index
        index=rules.star_snakesplice_index.output,
        # optional parameters
        extra=lambda wildcards, input: config["star_extra_settings"] + " --sjdbGTFfile " + input.gtf_file,
        # define output directory (not required for wrapper...)
        output_dir=lambda wildcards, output: os.path.dirname(output.aln)  
    threads: 1
    wrapper:
        "v3.3.3/bio/star/align"


rule star_snakesplice_align2:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=["output/extracted_reads/{sample}/{sample}_R1.{region}.fastq"],
        fq2=["output/extracted_reads/{sample}/{sample}_R2.{region}.fastq"],
        idx=rules.star_snakesplice_index.output,
        # GTF file is passed in extra
        gtf_file=config["reference_genome_annotation"],
    output:
        # see STAR manual for additional output files
        aln=os.path.join("output/star_snakesplice2/pe/{sample}.{region}.bam"),		# output only BAM, not SAM -> Use samtools view to convert to SAM
        log="logs/star_snakesplice2/pe/{sample}.{region}/Log.out",
        sj="output/star_snakesplice2/pe/{sample}.{region}/SJ.out.tab",
        unmapped=["output/star_snakesplice2/pe/{sample}.{region}/unmapped.1.fastq.gz",
            "output/star_snakesplice2/pe/{sample}.{region}/unmapped.2.fastq.gz"],
    log:
        "logs/star_snakesplice2/{sample}.{region}_alignment.log"
    params:
        # path to STAR reference genome index
        index=rules.star_snakesplice_index.output,
        # optional parameters
        extra=lambda wildcards, input: config["star_extra_settings2"] + " --sjdbGTFfile " + input.gtf_file,
        # define output directory (not required for wrapper...)
        output_dir=lambda wildcards, output: os.path.dirname(output.aln)  
    threads: 1
    wrapper:
        "v3.3.3/bio/star/align"

rule star_snakesplice_align3:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=["output/extracted_reads/{sample}/{sample}_R1.{region}.fastq"],
        fq2=["output/extracted_reads/{sample}/{sample}_R2.{region}.fastq"],
        idx=rules.star_snakesplice_index.output,
        # GTF file is passed in extra
        gtf_file=config["reference_genome_annotation"],
    output:
        # see STAR manual for additional output files
        aln=os.path.join("output/star_snakesplice3/pe/{sample}.{region}.bam"),		# output only BAM, not SAM -> Use samtools view to convert to SAM
        log="logs/star_snakesplice3/pe/{sample}.{region}/Log.out",
        sj="output/star_snakesplice3/pe/{sample}.{region}/SJ.out.tab",
        unmapped=["output/star_snakesplice3/pe/{sample}.{region}/unmapped.1.fastq.gz",
            "output/star_snakesplice3/pe/{sample}.{region}/unmapped.2.fastq.gz"],
    log:
        "logs/star_snakesplice3/{sample}.{region}_alignment.log"
    params:
        # path to STAR reference genome index
        index=rules.star_snakesplice_index.output,
        # optional parameters
        extra=lambda wildcards, input: config["star_extra_settings3"] + " --sjdbGTFfile " + input.gtf_file,
        # define output directory (not required for wrapper...)
        output_dir=lambda wildcards, output: os.path.dirname(output.aln)  
    threads: 1
    wrapper:
        "v3.3.3/bio/star/align"


rule star_snakesplice_align4:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=["output/extracted_reads/{sample}/{sample}_R1.{region}.fastq"],
        fq2=["output/extracted_reads/{sample}/{sample}_R2.{region}.fastq"],
        idx=rules.star_snakesplice_index.output,
        # GTF file is passed in extra
        gtf_file=config["reference_genome_annotation"],
    output:
        # see STAR manual for additional output files
        aln=os.path.join("output/star_snakesplice4/pe/{sample}.{region}.bam"),		# output only BAM, not SAM -> Use samtools view to convert to SAM
        log="logs/star_snakesplice4/pe/{sample}.{region}/Log.out",
        sj="output/star_snakesplice4/pe/{sample}.{region}/SJ.out.tab",
        unmapped=["output/star_snakesplice4/pe/{sample}.{region}/unmapped.1.fastq.gz",
            "output/star_snakesplice4/pe/{sample}.{region}/unmapped.2.fastq.gz"],
    log:
        "logs/star_snakesplice4/{sample}.{region}_alignment.log"
    params:
        # path to STAR reference genome index
        index=rules.star_snakesplice_index.output,
        # optional parameters
        extra=lambda wildcards, input: config["star_extra_settings4"] + " --sjdbGTFfile " + input.gtf_file,
        # define output directory (not required for wrapper...)
        output_dir=lambda wildcards, output: os.path.dirname(output.aln)  
    threads: 1
    wrapper:
        "v3.3.3/bio/star/align"

rule star_snakesplice_align5:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=["output/extracted_reads/{sample}/{sample}_R1.{region}.fastq"],
        fq2=["output/extracted_reads/{sample}/{sample}_R2.{region}.fastq"],
        idx=rules.star_snakesplice_index.output,
    output:
        # see STAR manual for additional output files
        aln=os.path.join("output/star_snakesplice5/pe/{sample}.{region}.bam"),		# output only BAM, not SAM -> Use samtools view to convert to SAM
        log="logs/star_snakesplice5/pe/{sample}.{region}/Log.out",
        sj="output/star_snakesplice5/pe/{sample}.{region}/SJ.out.tab",
        unmapped=["output/star_snakesplice5/pe/{sample}.{region}/unmapped.1.fastq.gz",
            "output/star_snakesplice5/pe/{sample}.{region}/unmapped.2.fastq.gz"],
    log:
        "logs/star_snakesplice5/{sample}.{region}_alignment.log"
    params:
        # path to STAR reference genome index
        index=rules.star_snakesplice_index.output,
        # optional parameters
        extra=lambda wildcards, input: config["star_extra_settings"],
        # define output directory (not required for wrapper...)
        output_dir=lambda wildcards, output: os.path.dirname(output.aln)  
    threads: 1
    wrapper:
        "v3.3.3/bio/star/align"

rule star_snakesplice_align6:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=["output/extracted_reads/{sample}/{sample}_R1.{region}.fastq"],
        # paired end reads needs to be ordered so each item in the two lists match
        fq2=["output/extracted_reads/{sample}/{sample}_R2.{region}.fastq"],  #optional
        # path to STAR reference genome index
        idx=rules.star_snakesplice_index.output,
    output:
        # see STAR manual for additional output files
        aln=os.path.join("output/star_snakesplice6/pe/{sample}.{region}.bam"),		# output only BAM, not SAM -> Use samtools view to convert to SAM
        log="logs/star_snakesplice6/pe/{sample}.{region}/Log.out",
        sj="output/star_snakesplice6/pe/{sample}.{region}/SJ.out.tab",
        unmapped=["output/star_snakesplice6/pe/{sample}.{region}/unmapped.1.fastq.gz",
            "output/star_snakesplice6/pe/{sample}.{region}/unmapped.2.fastq.gz"],
    log:
        "logs/star_snakesplice6/{sample}.{region}_alignment.log"
    params:
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate",
    threads: 8
    wrapper:
        "v3.3.6/bio/star/align"

rule star_snakesplice_align7:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=["output/extracted_reads/{sample}/{sample}_R1.{region}.fastq"],
        # paired end reads needs to be ordered so each item in the two lists match
        fq2=["output/extracted_reads/{sample}/{sample}_R2.{region}.fastq"],  #optional
        # path to STAR reference genome index
        idx=rules.star_snakesplice_index.output,
    output:
        # see STAR manual for additional output files
        aln=os.path.join("output/star_snakesplice7/pe/{sample}.{region}.bam"),		# output only BAM, not SAM -> Use samtools view to convert to SAM
        log="logs/star_snakesplice7/pe/{sample}.{region}/Log.out",
        sj="output/star_snakesplice7/pe/{sample}.{region}/SJ.out.tab",
        unmapped=["output/star_snakesplice7/pe/{sample}.{region}/unmapped.1.fastq.gz",
            "output/star_snakesplice7/pe/{sample}.{region}/unmapped.2.fastq.gz"],
    log:
        "logs/star_snakesplice7/{sample}.{region}_alignment.log"
    params:
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate",
    threads: 8
    wrapper:
        "v3.3.3/bio/star/align"

rule star_snakesplice_align8:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=["output/extracted_reads/{sample}/{sample}_R1.{region}.fastq"],
        # paired end reads needs to be ordered so each item in the two lists match
        fq2=["output/extracted_reads/{sample}/{sample}_R2.{region}.fastq"],  #optional
        # path to STAR reference genome index
        idx=rules.star_snakesplice_index.output,
    output:
        # see STAR manual for additional output files
        aln=os.path.join("output/star_snakesplice8/pe/{sample}.{region}.bam"),		# output only BAM, not SAM -> Use samtools view to convert to SAM
        log="logs/star_snakesplice8/pe/{sample}.{region}/Log.out",
        sj="output/star_snakesplice8/pe/{sample}.{region}/SJ.out.tab",
        unmapped=["output/star_snakesplice8/pe/{sample}.{region}/unmapped.1.fastq.gz",
            "output/star_snakesplice8/pe/{sample}.{region}/unmapped.2.fastq.gz"],
    log:
        "logs/star_snakesplice8/{sample}.{region}_alignment.log"
    params:
        # optional parameters
        extra=lambda wildcards, input: config["star_extra_settings"],
    threads: 8
    wrapper:
        "v3.3.3/bio/star/align"


rule star_snakesplice_align9:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=["output/extracted_reads/{sample}/{sample}_R1.{region}.fastq"],
        # paired end reads needs to be ordered so each item in the two lists match
        fq2=["output/extracted_reads/{sample}/{sample}_R2.{region}.fastq"],  #optional
        # path to STAR reference genome index
        idx=rules.star_snakesplice_index.output,
    output:
        # see STAR manual for additional output files
        aln=os.path.join("output/star_snakesplice9/pe/{sample}.{region}.bam"),		# output only BAM, not SAM -> Use samtools view to convert to SAM
        log="logs/star_snakesplice9/pe/{sample}.{region}/Log.out",
        sj="output/star_snakesplice9/pe/{sample}.{region}/SJ.out.tab",
        unmapped=["output/star_snakesplice9/pe/{sample}.{region}/unmapped.1.fastq.gz",
            "output/star_snakesplice9/pe/{sample}.{region}/unmapped.2.fastq.gz"],
    log:
        "logs/star_snakesplice9/{sample}.{region}_alignment.log"
    params:
        # optional parameters
        extra=lambda wildcards, input: config["star_extra_settings9"],
    threads: 8
    wrapper:
        "v3.3.3/bio/star/align"


rule star_snakesplice_align10:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=["output/extracted_reads/{sample}/{sample}_R1.{region}.fastq"],
        # paired end reads needs to be ordered so each item in the two lists match
        fq2=["output/extracted_reads/{sample}/{sample}_R2.{region}.fastq"],  #optional
        # path to STAR reference genome index
        idx=rules.star_snakesplice_index.output,
    output:
        # see STAR manual for additional output files
        aln=os.path.join("output/star_snakesplice10/pe/{sample}.{region}.bam"),		# output only BAM, not SAM -> Use samtools view to convert to SAM
        log="logs/star_snakesplice10/pe/{sample}.{region}/Log.out",
        sj="output/star_snakesplice10/pe/{sample}.{region}/SJ.out.tab",
        unmapped=["output/star_snakesplice10/pe/{sample}.{region}/unmapped.1.fastq.gz",
            "output/star_snakesplice10/pe/{sample}.{region}/unmapped.2.fastq.gz"],
    log:
        "logs/star_snakesplice10/{sample}.{region}_alignment.log"
    params:
        # optional parameters
        extra=lambda wildcards, input: config["star_extra_settings10"],
    threads: 8
    wrapper:
        "v3.3.3/bio/star/align"

rule star_snakesplice_align11:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=["output/extracted_reads/{sample}/{sample}_R1.{region}.fastq"],
        # paired end reads needs to be ordered so each item in the two lists match
        fq2=["output/extracted_reads/{sample}/{sample}_R2.{region}.fastq"],  #optional
        # path to STAR reference genome index
        idx=rules.star_snakesplice_index.output,
    output:
        # see STAR manual for additional output files
        aln=os.path.join("output/star_snakesplice11/pe/{sample}.{region}.bam"),		# output only BAM, not SAM -> Use samtools view to convert to SAM
        log="logs/star_snakesplice11/pe/{sample}.{region}/Log.out",
        sj="output/star_snakesplice11/pe/{sample}.{region}/SJ.out.tab",
        unmapped=["output/star_snakesplice11/pe/{sample}.{region}/unmapped.1.fastq.gz",
            "output/star_snakesplice11/pe/{sample}.{region}/unmapped.2.fastq.gz"],
    log:
        "logs/star_snakesplice11/{sample}.{region}_alignment.log"
    params:
        # optional parameters
        extra=lambda wildcards, input: config["star_extra_settings11"],
    threads: 8
    wrapper:
        "v3.3.3/bio/star/align"

rule star_snakesplice_align12:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=["output/extracted_reads/{sample}/{sample}_R1.{region}.fastq"],
        # paired end reads needs to be ordered so each item in the two lists match
        fq2=["output/extracted_reads/{sample}/{sample}_R2.{region}.fastq"],  #optional
        # path to STAR reference genome index
        idx=rules.star_snakesplice_index.output,
    output:
        # see STAR manual for additional output files
        aln=os.path.join("output/star_snakesplice12/pe/{sample}.{region}.bam"),		# output only BAM, not SAM -> Use samtools view to convert to SAM
        log="logs/star_snakesplice12/pe/{sample}.{region}/Log.out",
        sj="output/star_snakesplice12/pe/{sample}.{region}/SJ.out.tab",
        unmapped=["output/star_snakesplice12/pe/{sample}.{region}/unmapped.1.fastq.gz",
            "output/star_snakesplice12/pe/{sample}.{region}/unmapped.2.fastq.gz"],
    log:
        "logs/star_snakesplice12/{sample}.{region}_alignment.log"
    params:
        # optional parameters
        extra=lambda wildcards, input: config["star_extra_settings12"],
    threads: 8
    wrapper:
        "v3.3.3/bio/star/align"


rule star_snakesplice_align13:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1=["output/extracted_reads/{sample}/{sample}_R1.{region}.fastq"],
        # paired end reads needs to be ordered so each item in the two lists match
        fq2=["output/extracted_reads/{sample}/{sample}_R2.{region}.fastq"],  #optional
        # path to STAR reference genome index
        idx=rules.star_snakesplice_index.output,
        gtf_file=config["reference_genome_annotation"],
    output:
        # see STAR manual for additional output files
        aln=os.path.join("output/star_snakesplice13/pe/{sample}.{region}.bam"),		# output only BAM, not SAM -> Use samtools view to convert to SAM
        log="logs/star_snakesplice13/pe/{sample}.{region}/Log.out",
        sj="output/star_snakesplice13/pe/{sample}.{region}/SJ.out.tab",
        unmapped=["output/star_snakesplice13/pe/{sample}.{region}/unmapped.1.fastq.gz",
            "output/star_snakesplice13/pe/{sample}.{region}/unmapped.2.fastq.gz"],
    log:
        "logs/star_snakesplice13/{sample}.{region}_alignment.log"
    params:
        # optional parameters
        extra=lambda wildcards, input: config["star_extra_settings13"] + " --sjdbGTFfile " + input.gtf_file,
    threads: 8
    wrapper:
        "v3.3.3/bio/star/align"



import pysam
import sys
import os
 
# Attention: Reads mapped to the reverse strand are stored reverse complemented in the BAM file.
# Hence: No need to search for query sequence AND reverse complement
# It suffices to only search the query sequence
 
# Read bam file as CLI argument
# input_filename = os.path.basename(sys.argv[1])
# bamfile = pysam.AlignmentFile(sys.argv[1], "rb")
input_filename = snakemake.input["input_bam_file"]
bam_file = pysam.AlignmentFile(input_filename, "rb")

text_report_file = snakemake.output["text_report_file"]
filtered_reads_ids_gene_region_file = snakemake.output["filtered_reads_ids_gene_region_file"]
exon3_clipped_reads_ids_file = snakemake.output["exon3_clipped_reads_ids_file"]
exon2_clipped_reads_ids_file = snakemake.output["exon2_clipped_reads_ids_file"]
exon3_exon2_clipped_reads_ids_file = snakemake.output["exon3_exon2_clipped_reads_ids_file"]
 
 
# Filter region: hg19 - SUPT7L Exon3 & Exon 2
# chr2: 27883851 - 27885150
 
# ------ 0. Filter read IDs of reads that mapped to SUPT7L ------
# region of interest: chr2: 27873676 - 27886481
filtered_read_ids = []
for read in bam_file.fetch('2', 27873676, 27886481):
    filtered_read_ids.append(read.query_name)
print("Exporting read IDs of reads that mapped to SUPT7L")
with open(filtered_reads_ids_gene_region_file, "w") as f:
    for read_id in filtered_read_ids:
        f.write(read_id + "\n")
 
 
# ----- 1. Filter insertion & substitution reads at exon3 & exon2 -----
read_count_exon2_exon3 = 0
exon2_exon3_insertion_reads = []
exon2_exon3_insertion_reads_with_motif = []
exon2_exon3_substitution_reads = []
exon2_exon3_substitution_reads_with_motif = []
for read in bam_file.fetch('2', 27883851, 27885150):
    # Total number of reads
    read_count_exon2_exon3 += 1
 
    # 1. Filter reads with insertion of 4 bases (CTGG)
    if "I" in read.cigarstring and read.cigarstring.count("I") == 1:
        exon2_exon3_insertion_reads.append(read)
 
        # 2. Check if CTGG is the insertion as duplication
        if "CTGG" in read.query_sequence and read.query_sequence.count("CTGGCTGG") == 1:
            exon2_exon3_insertion_reads_with_motif.append(read)
 
 
    # 2. Filter reads with substitution of 1 base (G -> A)
    if read.has_tag('NM') and read.get_tag('NM') == 1:
        exon2_exon3_substitution_reads.append(read)
        # Check if the substitution is G -> A
        # ccacCggagt instead of ccacCggagt
        if 'CCACCGGAGT' in read.query_sequence:
            exon2_exon3_substitution_reads_with_motif.append(read)
 
 
# ----- 2. Filter clipped reads --------
 
# ------ 2.1 Exon3: Clipped sequences -----
# Right-clipped (F1R2) :
# TTTGCAGATTCATTGTCCATATGTCTCTTCAAG
exon3_clipped_reads_ids = []
exon3_reads_softclipped = []
exon3_reads_softclipped_with_motif10 = []
exon3_reads_softclipped_with_motif5 = []
exon3_reads_softclipped_with_motif3 = []
for read in bam_file.fetch('2', 27883851, 27884253):
    # Check if clipped sequence is equal to
    if "S" in read.cigarstring and read.cigarstring.count("S") == 1:
        exon3_clipped_reads_ids.append(read.query_name)
        exon3_reads_softclipped.append(read)
 
        query_alignment_start = read.query_alignment_start
        query_alignment_end = read.query_alignment_end
        clipped_sequence1 = read.query_sequence[0:query_alignment_start]
        clipped_sequence2 = read.query_sequence[query_alignment_end:]
 
        if "TTTGCAGATT" in clipped_sequence1 or "TTTGCAGATT" in clipped_sequence2:
            exon3_reads_softclipped_with_motif10.append(read)
        if "TTTGC" in clipped_sequence1 or "TTTGC" in clipped_sequence2:
            exon3_reads_softclipped_with_motif5.append(read)
        if "TTT" in clipped_sequence1 or "TTT" in clipped_sequence2:
            exon3_reads_softclipped_with_motif3.append(read)
 
# Export read IDs of clipped reads
with open(exon3_clipped_reads_ids_file, "w") as f:
    for read_id in exon3_clipped_reads_ids:
        f.write(read_id + "\n")
 
 
# ------ 2.2 Exon2: Clipped sequences -----
# Right-clipped - Forward
# CGGGGGGCTTCGGCTTGTTGGCTGAGGGTTGGTGCAGGGGTGGGTCATGGACTTCCACCAGACGGAACT
exon2_clipped_reads_ids = []
exon2_reads_softclipped = []
exon2_reads_softclipped_with_motif10 = []
exon2_reads_softclipped_with_motif5 = []
exon2_reads_softclipped_with_motif3 = []
 
for read in bam_file.fetch('2', 27885044, 27885150):
    # Check if clipped sequence is equal to
 
    if "S" in read.cigarstring and read.cigarstring.count("S") == 1:
        exon2_reads_softclipped.append(read)
        exon2_clipped_reads_ids.append(read.query_name)
 
        # Check if clipped sequence is equal to "AGACGGAACT"
        query_alignment_start = read.query_alignment_start
        query_alignment_end = read.query_alignment_end
        clipped_sequence1 = read.query_sequence[0:query_alignment_start]
        clipped_sequence2 = read.query_sequence[query_alignment_end:]
 
 
        if "AGACGGAACT" in clipped_sequence1 or "AGACGGAACT" in clipped_sequence2:
            exon2_reads_softclipped_with_motif10.append(read)
        if "AGACG" in clipped_sequence1 or "AGACG" in clipped_sequence2:
            exon2_reads_softclipped_with_motif5.append(read)
        if "AGA" in clipped_sequence1 or "AGA" in clipped_sequence2:
            exon2_reads_softclipped_with_motif3.append(read)
 
# Export read IDs of clipped reads
with open(exon2_clipped_reads_ids_file, "w") as f:
    for read_id in exon2_clipped_reads_ids:
        f.write(read_id + "\n")
 

# Both exon2 and exon3 clipped reads together
with open(exon3_exon2_clipped_reads_ids_file, "w") as f:
    output_list = exon3_clipped_reads_ids + exon2_clipped_reads_ids
    # Remove duplicates
    output_list = list(set(output_list))
    for read_id in output_list:
        f.write(read_id + "\n")

 
# Export results to text file
with open(text_report_file, "w") as f:
    f.write("Analysis of " + input_filename + "\n")
    f.write("-----------------------------------------------------------------\n")
    f.write("Total number of reads mapped to exon3 & exon2: " + str(read_count_exon2_exon3) + "\n")
    f.write("Number of reads with insertions: " + str(len(exon2_exon3_insertion_reads)) + "\n")
    f.write("Number of reads with insertion of 4 bases (CTGG) as duplication: " + str(len(exon2_exon3_insertion_reads_with_motif)) + "\n")
    f.write("Number of reads with substitution of 1 base: " + str(len(exon2_exon3_substitution_reads)) + "\n")
    f.write("Number of reads with substitution of 1 base (G -> A): " + str(len(exon2_exon3_substitution_reads_with_motif)) + "\n")
    f.write("Number of reads with softclipped sequence in exon3: " + str(len(exon3_reads_softclipped)) + "\n")
    f.write("Number of reads with softclipped sequence in exon3 with queried motif (10 bases): " + str(len(exon3_reads_softclipped_with_motif10)) + "\n")
    f.write("Number of reads with softclipped sequence in exon3 with queried motif (5 bases): " + str(len(exon3_reads_softclipped_with_motif5)) + "\n")
    f.write("Number of reads with softclipped sequence in exon3 with queried motif (3 bases): " + str(len(exon3_reads_softclipped_with_motif3)) + "\n")
    f.write("Number of reads with softclipped sequence in exon2: " + str(len(exon2_reads_softclipped)) + "\n")
    f.write("Number of reads with softclipped sequence in exon2 with queried motif (10 bases): " + str(len(exon2_reads_softclipped_with_motif10)) + "\n")
    f.write("Number of reads with softclipped sequence in exon2 with queried motif (5 bases): " + str(len(exon2_reads_softclipped_with_motif5)) + "\n")
    f.write("Number of reads with softclipped sequence in exon2 with queried motif (3 bases): " + str(len(exon2_reads_softclipped_with_motif3)) + "\n")


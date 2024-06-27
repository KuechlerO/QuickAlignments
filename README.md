# QuickAlignments
This project is designed to identify the mapping problem that occurs with STAR and our SUPT7L samples.

## The problem:
Many reads are not soft clipped instead of being split.
Why is this the case?
Is there a way to avoid this problem?


## The approach
Region of interest:     SUPT7L 
Observation:            Soft-clipping at exon 3 instead of splitting the reads

TODO: Capture all values in a table to make it comparable

1. Extract reads from region of interest that are soft clipped
    A. Reads aligned to SUPT7L ('2', 27873676, 27886481)
    B. Extract all reads from Exon3-Exon2 ('2', 27883851, 27885150)
    C. 4-base (CTGG) insertion in Exon3
    D. Check if insertion is duplication
    E. Soft-clipping at Exon3 ('2', 27883851, 27884253):
        E.1. TTTGCAGATT in clipped sequence?
        E.2. TTTGC in clipped sequence?
        E.3. TTT in clipped sequence?
    F. Soft-clipping at Exon2 ('2', 27885044, 27885150)
        F.1. AGACGGAACT in clipped sequence?
        F.2. AGACG in clipped sequence?
        F.3. AGA in clipped sequence?
2. Realign with same settings to validate
3. 

## TODO Steps:
1. Extract reads of interest
2. Get statistics for reads of interest
3. Apply changes to algorithm
4. Collect statistics and compare


## Results so far...
The problem is with the `outSAMstrandField intronMotif`, which was used to make the output Cufflinks/Cuffdiff compatible.

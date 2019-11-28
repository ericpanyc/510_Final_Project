#!/bin/bash
GTF=/scratch/yachenpa/TRGN510/rsem_GRCh38.p2.gtf
CPUS=1
MAPQ=0

featureCounts -p -Q $MAPQ -T $CPUS -a $GTF -o all_counts.txt *bam 

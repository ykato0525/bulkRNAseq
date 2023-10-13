#!/bin/bash

fastqc --nogroup -o ./MCF7_10Gy96h ~/RNAseq_analyzer/data/fastq/MCF7_10Gy_96h_1.fastq.gz

java -jar ~/RNAseq_analyzer/tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -trimlog log.txt -phred33 ~/RNAseq_analyzer/data/fastq/MCF7_10Gy_96h_1.fastq.gz ~/RNAseq_analyzer/data/fastq/MCF7_10Gy_96h_2.fastq.gz paired_output_1.fq unpaired_output_1.fq paired_output_2.fq unpaired_output_2.fq ILLUMINACLIP:/home/ykato/RNAseq_analyzer/tools/Trimmomatic-0.39/adapters/Truseq3-PE-2.fa:2:30:10 MINLEN:50

bowtie2 -p 4 --un-conc MCF710Gy96h_noribo.fq -x ~/RNAseq_analyzer/data/Homo_sapiens/UCSC/hg38/Sequence/AbundantSequences/humRibosomal -1 paired_output_1.fq -2 paired_output_2.fq -S MCF7_10Gy96h_ribosomal_hit.sam


mkdir STAR
cd STAR
STAR --genomeDir ~/RNAseq_analyzer/data/hg38_noalt/STAR_index/ --readFilesIn ../MCF710Gy96h_noribo.1.fq ../MCF710Gy96h_noribo.2.fq --outSAMtype BAM SortedByCoordinate --sjdbGTFfile ~/RNAseq_analyzer/data/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf --outFileNamePrefix STAR_MCF710Gy96h --runThreadN 4 --bamRemoveDuplicatesType UniqueIdentical
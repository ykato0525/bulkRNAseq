#!/bin/bash


# クオリティチェック
mkdir QC/target_of_QC
cd QC

fastqc --nogroup -o ./MCF7_10Gy96h ~/RNAseq_analyzer/data/fastq/MCF7_10Gy_96h_2.fastq.gz
fastqc --nogroup -o ./MDAQMB468NT ~/RNAseq_analyzer/data/fastq/468_NT_2.fastq.gz

java -jar ~/RNAseq_analyzer/tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -trimlog log.txt -phred33 ~/RNAseq_analyzer/data/fastq/MCF7_10Gy_96h_1.fastq.gz ~/RNAseq_analyzer/data/fastq/MCF7_10Gy_96h_2.fastq.gz paired_output_1.fq unpaired_output_1.fq paired_output_2.fq unpaired_output_2.fq ILLUMINACLIP:/home/ykato/RNAseq_analyzer/tools/Trimmomatic-0.39/adapters/Truseq3-PE-2.fa:2:30:10 MINLEN:50
cd 
bowtie2 -p 4 --un-conc MCF710Gy96h_noribo.fq -x ~/RNAseq_analyzer/data/Homo_sapiens/UCSC/hg38/Sequence/AbundantSequences/humRibosomal -1 paired_output_1.fq -2 paired_output_2.fq -S MCF7_10Gy96h_ribosomal_hit.sam

mkdir STAR_MCF710Gy96h
cd STAR_MCF710Gy96h

STAR --genomeDir ~/RNAseq_analyzer/data/hg38_noalt/STAR_index/ --readFilesIn ../MCF710Gy96h_noribo.1.fq ../MCF710Gy96h_noribo.2.fq --outSAMtype BAM SortedByCoordinate --sjdbGTFfile ~/RNAseq_analyzer/data/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf --outFileNamePrefix STAR_MCF710Gy96h --runThreadN 4 --bamRemoveDuplicatesType UniqueIdentical

more STAR_MCF710Gy96hLog.final.out

samtools view -bh -q 255 ./STAR_MCF710Gy96h/STAR_MCF710Gy96hAligned.sortedByCoord.out.bam > STAR_MCF710Gy96h.uniq.rmdup.bam
samtools index STAR_MCF710Gy96h.uniq.rmdup.bam


featureCounts -p -t exon -g gene_id -a ~/RNAseq_analyzer/data/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf -o MCF7_illumina_counts.txt STAR_MCF710Gy96h.uniq.rmdup.bam

#!/bin/bash

##Generate Genome Index (Genome Downloaded from Ensembel Jan 3 2018 v 2.37)
/home/megan/src/STAR-2.7.6a/source/STAR --runThreadN 7 --runMode genomeGenerate --genomeDir /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies --genomeFastaFiles /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_SCAFFOLDS_CpipJ2.fa --sjdbGTFfile /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gff3 --sjdbOverhang 150 --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 13.5



##Run Mapping Job with ENCODE standard options

#M2 - Parous CAL1
/home/megan/src/STAR-2.7.6a/source/STAR --runThreadN 7 --genomeDir /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies --readFilesCommand gunzip -c --readFilesIn /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M2-1_S5_R1_001_paired.fastq.gz /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M2-1_S5_R2_001_paired.fastq.gz --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

mv Aligned.out.sam ./DGE_GenAligned_SamFiles/M2-1_S5.sam
cat Log.final.out >> Log_summary

/home/megan/src/STAR-2.7.6a/source/STAR --runThreadN 7 --genomeDir /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies --readFilesCommand gunzip -c --readFilesIn /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M2-2_S6_R1_001_paired.fastq.gz /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M2-2_S6_R2_001_paired.fastq.gz --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

mv Aligned.out.sam ./M2-2_S6.sam
cat Log.final.out >> Log_summary

/home/megan/src/STAR-2.7.6a/source/STAR --runThreadN 7 --genomeDir /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies --readFilesCommand gunzip -c --readFilesIn /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M2-3_S7_R1_001_paired.fastq.gz /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M2-3_S7_R2_001_paired.fastq.gz --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

mv Aligned.out.sam ./M2-3_S7.sam
cat Log.final.out >> Log_summary

/home/megan/src/STAR-2.7.6a/source/STAR --runThreadN 7 --genomeDir /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies --readFilesCommand gunzip -c --readFilesIn /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M2-4_S8_R1_001_paired.fastq.gz /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M2-4_S8_R2_001_paired.fastq.gz --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

mv Aligned.out.sam ./M2-4_S8.sam
cat Log.final.out >> Log_summary


#M1 - gravid CAL1
/home/megan/src/STAR-2.7.6a/source/STAR --runThreadN 7 --genomeDir /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies --readFilesCommand gunzip -c --readFilesIn /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M1-1_S1_R1_001_paired.fastq.gz /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M1-1_S1_R2_001_paired.fastq.gz --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

mv Aligned.out.sam ./M1-1_S1.sam
cat Log.final.out >> Log_summary

/home/megan/src/STAR-2.7.6a/source/STAR --runThreadN 7 --genomeDir /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies --readFilesCommand gunzip -c --readFilesIn /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M1-2_S2_R1_001_paired.fastq.gz /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M1-2_S2_R2_001_paired.fastq.gz --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

mv Aligned.out.sam ./M1-2_S2.sam
cat Log.final.out >> Log_summary

/home/megan/src/STAR-2.7.6a/source/STAR --runThreadN 7 --genomeDir /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies --readFilesCommand gunzip -c --readFilesIn /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M1-3_S3_R1_001_paired.fastq.gz /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M1-3_S3_R2_001_paired.fastq.gz --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

mv Aligned.out.sam ./M1-3_S3.sam
cat Log.final.out >> Log_summary

/home/megan/src/STAR-2.7.6a/source/STAR --runThreadN 7 --genomeDir /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies --readFilesCommand gunzip -c --readFilesIn /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M1-4_S4_R1_001_paired.fastq.gz /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M1-4_S4_R2_001_paired.fastq.gz --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

mv Aligned.out.sam ./M1-4_S4.sam
cat Log.final.out >> Log_summary

#M4 Pipiens
/home/megan/src/STAR-2.7.6a/source/STAR --runThreadN 7 --genomeDir /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies --readFilesCommand gunzip -c --readFilesIn /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M4-1_S13_R1_001_paired.fastq.gz /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M4-1_S13_R2_001_paired.fastq.gz --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

mv Aligned.out.sam ./M4-1_S13.sam
cat Log.final.out >> Log_summary

/home/megan/src/STAR-2.7.6a/source/STAR --runThreadN 7 --genomeDir /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies --readFilesCommand gunzip -c --readFilesIn /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M4-2_S14_R1_001_paired.fastq.gz /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M4-2_S14_R2_001_paired.fastq.gz --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

mv Aligned.out.sam ./M4-2_S14.sam
cat Log.final.out >> Log_summary

/home/megan/src/STAR-2.7.6a/source/STAR --runThreadN 7 --genomeDir /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies --readFilesCommand gunzip -c --readFilesIn /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M4-3_S15_R1_001_paired.fastq.gz /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M4-3_S15_R2_001_paired.fastq.gz --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

mv Aligned.out.sam ./M4-3_S15.sam
cat Log.final.out >> Log_summary

/home/megan/src/STAR-2.7.6a/source/STAR --runThreadN 7 --genomeDir /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies --readFilesCommand gunzip -c --readFilesIn /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M4-4_S16_R1_001_paired.fastq.gz /media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/M4-4_S16_R2_001_paired.fastq.gz --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000

mv Aligned.out.sam ./M4-4_S16.sam
cat Log.final.out >> Log_summary

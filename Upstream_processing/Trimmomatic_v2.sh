#!/bin/bash
#This is the script that I used to trim adapter sequences from my reads

cd /media/megan/EE9C16C89C168AEB/Noreuil

#M1 - gravid CAL1
java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 M1-1_S1_R1_001.fastq.gz M1-1_S1_R2_001.fastq.gz M1-1_S1_R1_001_paired.fastq.gz M1-1_S1_R1_001_unpaired.fastq.gz M1-1_S1_R2_001_paired.fastq.gz M1-1_S1_R2_001_unpaired.fastq.gz ILLUMINACLIP:/home/megan/src/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 M1-2_S2_R1_001.fastq.gz M1-2_S2_R2_001.fastq.gz M1-2_S2_R1_001_paired.fastq.gz M1-2_S2_R1_001_unpaired.fastq.gz M1-2_S2_R2_001_paired.fastq.gz M1-2_S2_R2_001_unpaired.fastq.gz ILLUMINACLIP:/home/megan/src/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 M1-3_S3_R1_001.fastq.gz M1-3_S3_R2_001.fastq.gz M1-3_S3_R1_001_paired.fastq.gz M1-3_S3_R1_001_unpaired.fastq.gz M1-3_S3_R2_001_paired.fastq.gz M1-3_S3_R2_001_unpaired.fastq.gz ILLUMINACLIP:/home/megan/src/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 M1-4_S4_R1_001.fastq.gz M1-4_S4_R2_001.fastq.gz M1-4_S4_R1_001_paired.fastq.gz M1-4_S4_R1_001_unpaired.fastq.gz M1-4_S4_R2_001_paired.fastq.gz M1-4_S4_R2_001_unpaired.fastq.gz ILLUMINACLIP:/home/megan/src/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50

#M2 - Parous CAL1
java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 M2-1_S5_R1_001.fastq.gz M2-1_S5_R2_001.fastq.gz M2-1_S5_R1_001_paired.fastq.gz M2-1_S5_R1_001_unpaired.fastq.gz M2-1_S5_R2_001_paired.fastq.gz M2-1_S5_R2_001_unpaired.fastq.gz ILLUMINACLIP:/home/megan/src/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 M2-2_S6_R1_001.fastq.gz M2-2_S6_R2_001.fastq.gz M2-2_S6_R1_001_paired.fastq.gz M2-2_S6_R1_001_unpaired.fastq.gz M2-2_S6_R2_001_paired.fastq.gz M2-2_S6_R2_001_unpaired.fastq.gz ILLUMINACLIP:/home/megan/src/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 M2-3_S7_R1_001.fastq.gz M2-3_S7_R2_001.fastq.gz M2-3_S7_R1_001_paired.fastq.gz M2-3_S7_R1_001_unpaired.fastq.gz M2-3_S7_R2_001_paired.fastq.gz M2-3_S7_R2_001_unpaired.fastq.gz ILLUMINACLIP:/home/megan/src/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 M2-4_S8_R1_001.fastq.gz M2-4_S8_R2_001.fastq.gz M2-4_S8_R1_001_paired.fastq.gz M2-4_S8_R1_001_unpaired.fastq.gz M2-4_S8_R2_001_paired.fastq.gz M2-4_S8_R2_001_unpaired.fastq.gz ILLUMINACLIP:/home/megan/src/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50

#M3 - CAL1 Males
java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 M3-1_S9_R1_001.fastq.gz M3-1_S9_R2_001.fastq.gz M3-1_S9_R1_001_paired.fastq.gz M3-1_S9_R1_001_unpaired.fastq.gz M3-1_S9_R2_001_paired.fastq.gz M3-1_S9_R2_001_unpaired.fastq.gz ILLUMINACLIP:/home/megan/src/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 M3-2_S10_R1_001.fastq.gz M3-2_S10_R2_001.fastq.gz M3-2_S10_R1_001_paired.fastq.gz M3-2_S10_R1_001_unpaired.fastq.gz M3-2_S10_R2_001_paired.fastq.gz M3-2_S10_R2_001_unpaired.fastq.gz ILLUMINACLIP:/home/megan/src/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 M3-3_S11_R1_001.fastq.gz M3-3_S11_R2_001.fastq.gz M3-3_S11_R1_001_paired.fastq.gz M3-3_S11_R1_001_unpaired.fastq.gz M3-3_S11_R2_001_paired.fastq.gz M3-3_S11_R2_001_unpaired.fastq.gz ILLUMINACLIP:/home/megan/src/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4  M3-4_S12_R1_001.fastq.gz M3-4_S12_R2_001.fastq.gz M3-4_S12_R1_001_paired.fastq.gz M3-4_S12_R1_001_unpaired.fastq.gz M3-4_S12_R2_001_paired.fastq.gz M3-4_S12_R2_001_unpaired.fastq.gz ILLUMINACLIP:/home/megan/src/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50

#M4 - IL2 nulliparous
java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 M4-1_S13_R1_001.fastq.gz M4-1_S13_R2_001.fastq.gz M4-1_S13_R1_001_paired.fastq.gz M4-1_S13_R1_001_unpaired.fastq.gz M4-1_S13_R2_001_paired.fastq.gz M4-1_S13_R2_001_unpaired.fastq.gz ILLUMINACLIP:/home/megan/src/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 M4-2_S14_R1_001.fastq.gz M4-2_S14_R2_001.fastq.gz M4-2_S14_R1_001_paired.fastq.gz M4-2_S14_R1_001_unpaired.fastq.gz M4-2_S14_R2_001_paired.fastq.gz M4-2_S14_R2_001_unpaired.fastq.gz ILLUMINACLIP:/home/megan/src/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 M4-3_S15_R1_001.fastq.gz M4-3_S15_R2_001.fastq.gz M4-3_S15_R1_001_paired.fastq.gz M4-3_S15_R1_001_unpaired.fastq.gz M4-3_S15_R2_001_paired.fastq.gz M4-3_S15_R2_001_unpaired.fastq.gz ILLUMINACLIP:/home/megan/src/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50

java -jar /home/megan/src/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 M4-4_S16_R1_001.fastq.gz M4-4_S16_R2_001.fastq.gz M4-4_S16_R1_001_paired.fastq.gz M4-4_S16_R1_001_unpaired.fastq.gz M4-4_S16_R2_001_paired.fastq.gz M4-4_S16_R2_001_unpaired.fastq.gz ILLUMINACLIP:/home/megan/src/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10:8:true SLIDINGWINDOW:5:20 MINLEN:50

#moving trimmed pairs and trimmed unpaired files to separate folders.

mkdir trimmed_pairs trimmed_unpaired

mv *_paired.fastq.gz ./trimmed_pairs
mv *_unpaired.fastq.gz ./trimmed_unpaired



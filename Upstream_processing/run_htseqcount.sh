#!/bin/bash

#Script to count uniquely aligned reads
#M.Fritz 10-19-20

#gff to gtf annotation conversion
gffread /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gff3 -T -o /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gtf


#HTSeq version 0.11.1 - our libraries stranded.

cd /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/DGE_GenAligned_SamFiles

mkdir highQual_gene highQual_exon Str_Yes

#highQual_exon
htseq-count -f bam -a 20 -r pos -s reverse -t exon M1-1_S1.bam  /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gtf > ./highQual_exon/M1-1_S1_htseq
htseq-count -f bam -a 20 -r pos -s reverse -t exon M1-2_S2.bam  /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gtf > ./highQual_exon/M1-2_S2_htseq
htseq-count -f bam -a 20 -r pos -s reverse -t exon  M1-3_S3.bam  /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gtf > ./highQual_exon/M1-3_S3_htseq
htseq-count -f bam -a 20 -r pos -s reverse -t exon  M1-4_S4.bam  /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gtf > ./highQual_exon/M1-4_S4_htseq
htseq-count -f bam -a 20 -r pos -s reverse -t exon  M2-1_S5.bam  /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gtf > ./highQual_exon/M2-1_S5_htseq
htseq-count -f bam -a 20 -r pos -s reverse -t exon  M2-2_S6.bam  /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gtf > ./highQual_exon/M2-2_S6_htseq
htseq-count -f bam -a 20 -r pos -s reverse -t exon  M2-3_S7.bam  /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gtf > ./highQual_exon/M2-3_S7_htseq
htseq-count -f bam -a 20 -r pos -s reverse -t exon  M2-4_S8.bam  /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gtf > ./highQual_exon/M2-4_S8_htseq
htseq-count -f bam -a 20 -r pos -s reverse -t exon  M4-1_S13.bam /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gtf > ./highQual_exon/M4-1_S13_htseq
htseq-count -f bam -a 20 -r pos -s reverse -t exon  M4-2_S14.bam /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gtf > ./highQual_exon/M4-2_S14_htseq
htseq-count -f bam -a 20 -r pos -s reverse -t exon  M4-3_S15.bam /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gtf > ./highQual_exon/M4-3_S15_htseq
htseq-count -f bam -a 20 -r pos -s reverse -t exon  M4-4_S16.bam /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gtf > ./highQual_exon/M4-4_S16_htseq

#highQual_gene example - wasn't sure about the function of --idattr so I tested it.  Here is an example command.
htseq-count -f bam -a 20 -r pos -s reverse --idattr=gene_id  M2-1_S5.bam  /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gtf > ./highQual_gene/M2-1_S5_htseq

#DiffStrandedness_gene - quick check.
htseq-count -f bam -a 20 -r pos -s yes --idattr=gene_id  M2-1_S5.bam  /media/fritzlab/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gtf > ./Str_Yes/M2-1_S5_htseq

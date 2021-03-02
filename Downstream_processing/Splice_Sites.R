#Script to analyze splice variants
#02222021 MF

library(reshape2); library(dplyr)

#diff exp dataset
diff_exp_CALparous_v_Pip <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/DGE_GenAligned_SamFiles/highQual_exon/CALparousF_v_PipEvanF_LFCS_padj05.txt", header = F)
diff_exp_CALgravid_v_par <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/DGE_GenAligned_SamFiles/highQual_exon/CALgravidF_v_CALparousF_LFCS_padj05.txt", header = F)
Hi_Exp_Calgrav <- subset(diff_exp_CALgravid_v_par, V2 > 0)
Lo_Exp_Calgrav <- subset(diff_exp_CALgravid_v_par, V2 < 0)
Hi_exp_CAL <- subset(diff_exp_CALparous_v_Pip, V2 < 0)
Hi_exp_Pip <- subset(diff_exp_CALparous_v_Pip, V2 > 0)

#loading my datasets - M1s or BG gravid
M1_1 <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/SJout_files/M1-1_S1_SJ.out.tab", header = F)#checking against M2s to get a sense of differences
M1_2 <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/SJout_files/M1-2_S2_SJ.out.tab", header = F)
M1_3 <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/SJout_files/M1-3_S3_SJ.out.tab", header = F)
M1_4 <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/SJout_files/M1-4_S4_SJ.out.tab", header = F)

#loading my datasets - M2s or BG parous
M2_1 <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/SJout_files/M2-1_S5_SJ.out.tab", header = F)
M2_2 <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/SJout_files/M2-2_S6_SJ.out.tab", header = F)
M2_3 <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/SJout_files/M2-3_S7_SJ.out.tab", header = F)
M2_4 <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/SJout_files/M2-4_S8_SJ.out.tab", header = F)

#loading my datasets - M4s or AG
M4_1 <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/SJout_files/M4-1_S13_SJ.out.tab", header = F)
M4_2 <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/SJout_files/M4-2_S14_SJ.out.tab", header = F)
M4_3 <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/SJout_files/M4-3_S15_SJ.out.tab", header = F)
M4_4 <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/SJout_files/M4-4_S16_SJ.out.tab", header = F)

#Adding dataset names to columns
my_list <- list("M1_1" = M1_1, "M1_2" = M1_2, "M1_3" = M1_3, "M1_4" = M1_4, "M2_1" = M2_1, "M2_2" = M2_2, "M2_3" = M2_3, "M2_4" = M2_4, "M4_1" = M4_1, "M4_2" = M4_2, "M4_3" = M4_3, "M4_4" = M4_4)
my_list <- Map(cbind, my_list, new_clumn = names(my_list))

lapply(names(my_list), function(x) assign(x, my_list[[x]], envir = .GlobalEnv))#splits list of dfs back out


#reading in gtf to get gene names associated with alt splice sites
gtf <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/mosquito_genome_assemblies/Culex-quinquefasciatus-Johannesburg_BASEFEATURES_CpipJ2.4.gtf", sep="\t", stringsAsFactors=F)
ids <- sapply(strsplit(gtf$V9,";"),.subset,1)
gtf$ID <- sapply(strsplit(ids, " "),.subset,2)#gives VB-formatted transcript IDs

#Gathering POSition information for joins
my_cols <- c("V1", "V2", "V3")

M1_1$Pos <- do.call(paste, c(M1_1[my_cols], sep = "_"))
M1_2$Pos <- do.call(paste, c(M1_2[my_cols], sep = "_"))
M1_3$Pos <- do.call(paste, c(M1_3[my_cols], sep = "_"))
M1_4$Pos <- do.call(paste, c(M1_4[my_cols], sep = "_"))

M2_1$Pos <- do.call(paste, c(M2_1[my_cols], sep = "_"))
M2_2$Pos <- do.call(paste, c(M2_2[my_cols], sep = "_"))
M2_3$Pos <- do.call(paste, c(M2_3[my_cols], sep = "_"))
M2_4$Pos <- do.call(paste, c(M2_4[my_cols], sep = "_"))

M4_1$Pos <- do.call(paste, c(M4_1[my_cols], sep = "_"))
M4_2$Pos <- do.call(paste, c(M4_2[my_cols], sep = "_"))
M4_3$Pos <- do.call(paste, c(M4_3[my_cols], sep = "_"))
M4_4$Pos <- do.call(paste, c(M4_4[my_cols], sep = "_"))

#shared high confidence splice sites for each pop, b/c present in all reps.
M1s_intersect <- Reduce(intersect, list(M1_1$Pos, M1_2$Pos, M1_3$Pos, M1_4$Pos))
M2s_intersect <- Reduce(intersect, list(M2_1$Pos, M2_2$Pos, M2_3$Pos, M2_4$Pos))
M4s_intersect <- Reduce(intersect, list(M4_1$Pos, M4_2$Pos, M4_3$Pos, M4_4$Pos))

#should I bother to look at M1s vs. M2s?
M1vM2 <- setdiff(M1s_intersect, M2s_intersect) #there are potentially 7100 diffs, so I should atleast 


#splice site is present in at least one of reps of the population - for comparison against the intersects above - helps if cov. depth is low for subset of samples.
M1s_union <- Reduce(union, list(M1_1$Pos, M1_2$Pos, M1_3$Pos, M1_4$Pos))
M2s_union <- Reduce(union, list(M2_1$Pos, M2_2$Pos, M2_3$Pos, M2_4$Pos))
M4s_union <- Reduce(union, list(M4_1$Pos, M4_2$Pos, M4_3$Pos, M4_4$Pos))


#shared high confidence splice sites between pops
Tot_splice_intersect1 <- intersect(M1s_intersect, M2s_intersect)
Tot_splice_intersect2 <- intersect(M2s_intersect, M4s_intersect)

#diff between pops
Tot_splice_uniqM1s <- data.frame(setdiff(M1s_intersect, M2s_union)) 
names(Tot_splice_uniqM1s)[1] <- "Pos"

Tot_splice_uniqM2_1s <- data.frame(setdiff(M2s_intersect, M1s_union))
names(Tot_splice_uniqM2_1s)[1] <- "Pos"

Tot_splice_uniqM2_2s <- data.frame(setdiff(M2s_intersect, M4s_union))
names(Tot_splice_uniqM2_2s)[1] <- "Pos"

Tot_splice_uniqM4s <- data.frame(setdiff(M4s_intersect, M2s_union))#more unique splice sites for AG
names(Tot_splice_uniqM4s)[1] <- "Pos"

nrow(Tot_splice_uniqM1s)#but this drops it to 177....getting less interesting
nrow(Tot_splice_uniqM2_1s)
nrow(Tot_splice_uniqM2_2s)
nrow(Tot_splice_uniqM4s)


#merging diffs back with full dataset to get strand (V4) information.
#M1s
M1_1_merged <- merge(M1_1, Tot_splice_uniqM1s, by = "Pos")
M1_2_merged <- merge(M1_2, Tot_splice_uniqM1s, by = "Pos")
M1_3_merged <- merge(M1_3, Tot_splice_uniqM1s, by = "Pos")
M1_4_merged <- merge(M1_4, Tot_splice_uniqM1s, by = "Pos")

table(M1_1_merged$V4) #0s present for strandedness seem to represent non-annotated regions based on visual inspection, could also be DNA contamination so filtering out.

M1s_bound <- rbind(M1_1_merged, M1_2_merged, M1_3_merged, M1_4_merged)

#M2s - diff by comparison
M2_1_merged_compBGp <- merge(M2_1, Tot_splice_uniqM2_1s, by = "Pos")
M2_2_merged_compBGp <- merge(M2_2, Tot_splice_uniqM2_1s, by = "Pos")
M2_3_merged_compBGp <- merge(M2_3, Tot_splice_uniqM2_1s, by = "Pos")
M2_4_merged_compBGp <- merge(M2_4, Tot_splice_uniqM2_1s, by = "Pos")

table(M2_1_merged_compBGp$V4) 

M2_1_merged_compAG2 <- merge(M2_1, Tot_splice_uniqM2_2s, by = "Pos")
M2_2_merged_compAG2 <- merge(M2_2, Tot_splice_uniqM2_2s, by = "Pos")
M2_3_merged_compAG2 <- merge(M2_3, Tot_splice_uniqM2_2s, by = "Pos")
M2_4_merged_compAG2 <- merge(M2_4, Tot_splice_uniqM2_2s, by = "Pos")

table(M2_1_merged_compAG2$V4)

M2s_bound_compBGp <- rbind(M2_1_merged_compBGp, M2_2_merged_compBGp, M2_3_merged_compBGp, M2_4_merged_compBGp)
M2s_bound_compAG2 <- rbind(M2_1_merged_compAG2, M2_2_merged_compAG2, M2_3_merged_compAG2, M2_4_merged_compAG2)

#M4s
M4_1_merged <- merge(M4_1, Tot_splice_uniqM4s, by = "Pos")
M4_2_merged <- merge(M4_2, Tot_splice_uniqM4s, by = "Pos")
M4_3_merged <- merge(M4_3, Tot_splice_uniqM4s, by = "Pos")
M4_4_merged <- merge(M4_4, Tot_splice_uniqM4s, by = "Pos")

table(M4_1_merged$V4)

M4s_bound <- rbind(M4_1_merged, M4_2_merged, M4_3_merged, M4_4_merged)


#Strand information 

#M1s - stranded
M1s_bound_stranded<- subset(M1s_bound, V4 != 0)

#M2s - stranded
M2s_bound_stranded_compBGp <- subset(M2s_bound_compBGp, V4 != 0)
M2s_bound_stranded_compAG2<- subset(M2s_bound_compAG2, V4 != 0)

#M4s - stranded
M4s_bound_stranded <- subset(M4s_bound, V4 != 0)

#Setting up to filter on avg RC - must have at least 100
Norm_counts <- read.csv("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/DGE_GenAligned_SamFiles/highQual_exon/Normalized_ReadCounts_AllTreats.csv", header = T)
head(Norm_counts)

Avg_BGgrav <- rowMeans(Norm_counts[,c(2:5)])
Avg_BGgrav_genes <- data.frame(cbind(as.character(Norm_counts$X), as.numeric(Avg_BGgrav)))
Avg_BGgrav_genes$X1 <- paste0(Avg_BGgrav_genes$X1,"-","RA")
sub_Avg_BGgrav_genes <- subset(Avg_BGgrav_genes, Avg_BGgrav > 100)

Avg_BG1 <- rowMeans(Norm_counts[,c(6:9)])
Avg_BG1_genes <- data.frame(cbind(as.character(Norm_counts$X), as.numeric(Avg_BG1)))
Avg_BG1_genes$X1 <- paste0(Avg_BG1_genes$X1,"-","RA")
sub_Avg_BG1_genes <- subset(Avg_BG1_genes, Avg_BG1 > 100)

Avg_AG2 <- rowMeans(Norm_counts[,c(10:13)])
Avg_AG2_genes <- data.frame(cbind(as.character(Norm_counts$X), as.numeric(Avg_AG2)))
Avg_AG2_genes$X1 <- paste0(Avg_AG2_genes$X1,"-","RA")
sub_Avg_AG2_genes <- subset(Avg_AG2_genes, Avg_AG2 > 100)


Final_Avg_genes1 <- merge(sub_Avg_BG1_genes, sub_Avg_BGgrav_genes, by = "X1")#total number of genes that met filtering by read count criteria and could have had splice variants
Final_Avg_genes2 <- merge(sub_Avg_BG1_genes, sub_Avg_AG2_genes, by = "X1")#total number of genes that met filtering by read count criteria and could have had splice variants


#long to wide format conversion - sanity check, match nrow with table results from above.
M1s_stranded_wide <- reshape(M1s_bound_stranded, 
                             timevar = "new_clumn",
                             idvar = c("Pos", "V1", "V2", "V3", "V4"),
                             direction = "wide")

M2s_stranded_wide_compBGp <- reshape(M2s_bound_stranded_compBGp, 
                                     timevar = "new_clumn",
                                     idvar = c("Pos", "V1", "V2", "V3", "V4"),
                                     direction = "wide")

M2s_stranded_wide_compAG2 <- reshape(M2s_bound_stranded_compAG2, 
                             timevar = "new_clumn",
                             idvar = c("Pos", "V1", "V2", "V3", "V4"),
                             direction = "wide")

M4s_stranded_wide <- reshape(M4s_bound_stranded, 
                         timevar = "new_clumn",
                         idvar = c("Pos", "V1", "V2", "V3", "V4"),
                         direction = "wide")



#Adding columns to wide_formatted datasets and gtf for merges
gtf$intron_start <- gtf$V5 + 1

my_gtf_start <- c("V1", "intron_start")#gets combined contig and intron start
gtf$intron_start <- do.call(paste, c(gtf[my_gtf_start], sep = "_"))

my_df_start <- c("V1", "V2")#gets combined contig and intron start from dfs

M1s_stranded_wide$intron_start <- do.call(paste, c(M1s_stranded_wide[my_df_start], sep = "_"))
M2s_stranded_wide_compBGp$intron_start <- do.call(paste, c(M2s_stranded_wide_compBGp[my_df_start], sep = "_"))
M2s_stranded_wide_compAG2$intron_start <- do.call(paste, c(M2s_stranded_wide_compAG2[my_df_start], sep = "_"))
M4s_stranded_wide$intron_start <- do.call(paste, c(M4s_stranded_wide[my_df_start], sep = "_"))


#merging datasets to add genes - notice that only a portion of the unique splice variants seem to be associated with annotated genes.
#BG1 gravid to parous comp
M1s_stranded_wide_withGenes <- merge(gtf, M1s_stranded_wide, by = "intron_start")
M1s_stranded_wide_filtRC <- merge(M1s_stranded_wide_withGenes, Final_Avg_genes1, by.x = "ID", by.y = "X1")#here's the actual filter by avg RC

length(unique(M1s_stranded_wide_filtRC$Pos))#num BG1-gravid specific high quality splice variants
unique(M1s_stranded_wide_filtRC$ID)#genes containing the BG1-gravid specific splice variants

M2s_stranded_wide_withGenes_compBGp <- merge(gtf, M2s_stranded_wide_compBGp, by = "intron_start")
M2s_stranded_wide_filtRC_compBGp <- merge(M2s_stranded_wide_withGenes_compBGp, Final_Avg_genes1, by.x = "ID", by.y = "X1")#here's the actual filter by avg RC

length(unique(M2s_stranded_wide_filtRC_compBGp$Pos))#num BG1 par specific high quality splice variants
unique(M2s_stranded_wide_filtRC_compBGp$ID)#genes containing the BG1 par specific splice variants

#BG1 parous AG2 comp
M2s_stranded_wide_withGenes_compAG2 <- merge(gtf, M2s_stranded_wide_compAG2, by = "intron_start")
M2s_stranded_wide_filtRC_compAG2 <- merge(M2s_stranded_wide_withGenes_compAG2, Final_Avg_genes2, by.x = "ID", by.y = "X1")#here's the actual filter by avg RC

length(unique(M2s_stranded_wide_filtRC_compAG2$Pos))#num BG1 specific high quality splice variants
unique(M2s_stranded_wide_filtRC_compAG2$ID)#genes containing the BG1-specific splice variants

M4s_stranded_wide_withGenes <- merge(gtf, M4s_stranded_wide, by = "intron_start")
M4s_stranded_wide_filtRC <- merge(M4s_stranded_wide_withGenes, Final_Avg_genes2, by.x = "ID", by.y = "X1")

length(unique(M4s_stranded_wide_filtRC$Pos))#num AG2 specific high quality splice variants
unique(M4s_stranded_wide_filtRC$ID)#genes containing the AG2-specific splice variants


#writing Data S2 output
M1s_DataS2 <- M1s_stranded_wide_filtRC[,c(1,12,13,14,15,16,19,24,29,34)]
M1s_DataS2 <- cbind(Samp = "BG1grav", M1s_DataS2)
colnames(M1s_DataS2)<- c("Samp","Gene_ID","Junc_Name", "Scaffold", "Junc_Start", "Junc_Stop", "Strand", "UniqReadsMappingJunc_Samp1","UniqReadsMappingJunc_Samp2", "UniqReadsMappingJunc_Samp3", "UniqReadsMappingJunc_Samp4")
M1s_DataS2_uniq <- M1s_DataS2 %>% distinct(Samp, Junc_Name, .keep_all = TRUE)

M2s_compgrav_DataS2 <- M2s_stranded_wide_filtRC_compBGp[,c(1,12,13,14,15,16,19,24,29,34)]
M2s_compgrav_DataS2 <- cbind(Samp = "BG1par_compgrav", M2s_compgrav_DataS2)
colnames(M2s_compgrav_DataS2)<- c("Samp","Gene_ID","Junc_Name", "Scaffold", "Junc_Start", "Junc_Stop", "Strand", "UniqReadsMappingJunc_Samp1","UniqReadsMappingJunc_Samp2", "UniqReadsMappingJunc_Samp3", "UniqReadsMappingJunc_Samp4")
M2s_compgrav_DataS2_uniq <- M2s_compgrav_DataS2 %>% distinct(Samp, Junc_Name, .keep_all = TRUE)

M2s_compAG_DataS2 <- M2s_stranded_wide_filtRC_compAG2[,c(1,12,13,14,15,16,19,24,29,34)]
M2s_compAG_DataS2 <- cbind(Samp = "BG1par_compAG", M2s_compAG_DataS2)
colnames(M2s_compAG_DataS2)<- c("Samp","Gene_ID","Junc_Name", "Scaffold", "Junc_Start", "Junc_Stop", "Strand", "UniqReadsMappingJunc_Samp1","UniqReadsMappingJunc_Samp2", "UniqReadsMappingJunc_Samp3", "UniqReadsMappingJunc_Samp4")
M2s_compAG_DataS2_uniq <- M2s_compAG_DataS2 %>% distinct(Samp, Junc_Name, .keep_all = TRUE)

M4s_DataS2 <- M4s_stranded_wide_filtRC[,c(1,12,13,14,15,16,19,24,29,34)]
M4s_DataS2 <- cbind(Samp = "AG2", M4s_DataS2)
colnames(M4s_DataS2)<- c("Samp","Gene_ID","Junc_Name","Scaffold", "Junc_Start", "Junc_Stop", "Strand", "UniqReadsMappingJunc_Samp1","UniqReadsMappingJunc_Samp2", "UniqReadsMappingJunc_Samp3", "UniqReadsMappingJunc_Samp4")
M4s_DataS2_uniq <- M4s_DataS2 %>% distinct(Samp, Junc_Name, .keep_all = TRUE)

For_DataS2 <- rbind(M1s_DataS2_uniq, M2s_compgrav_DataS2_uniq, M2s_compAG_DataS2_uniq, M4s_DataS2_uniq)
write.table(For_DataS2, file = "Data_S2_SpliceJunc.txt", row.names = F)


#merging splice variants with diff_Exp datasets
M1s_stranded_diffExp_diffSplice_HiCalgrav <- intersect(M1s_stranded_wide_filtRC$ID, Hi_Exp_Calgrav$V1)
M1s_stranded_diffExp_diffSplice_HiCalgrav
length(M1s_stranded_diffExp_diffSplice_HiCalgrav)

M1s_stranded_diffExp_diffSplice_LoCalgrav <- intersect(M1s_stranded_wide_filtRC$ID, Lo_Exp_Calgrav$V1)
M1s_stranded_diffExp_diffSplice_LoCalgrav
length(M1s_stranded_diffExp_diffSplice_LoCalgrav)

M2s_stranded_diffExp_diffSplice_HiCal <- intersect(M2s_stranded_wide_filtRC$ID, Hi_exp_CAL$V1)
M2s_stranded_diffExp_diffSplice_HiCal
length(M2s_stranded_diffExp_diffSplice_HiCal)

M2s_stranded_diffExp_diffSplice_HiPip <- intersect(M2s_stranded_wide_filtRC$ID, Hi_exp_Pip$V1)
M2s_stranded_diffExp_diffSplice_HiPip
length(M2s_stranded_diffExp_diffSplice_HiPip)

M4s_stranded_diffExp_diffSplice_HiCal <- intersect(M4s_stranded_wide_filtRC$ID, Hi_exp_CAL$V1)
M4s_stranded_diffExp_diffSplice_HiCal
length(M4s_stranded_diffExp_diffSplice_HiCal)

M4s_stranded_diffExp_diffSplice_HiPip <- intersect(M4s_stranded_wide_filtRC$ID, Hi_exp_Pip$V1)
M4s_stranded_diffExp_diffSplice_HiPip
length(M4s_stranded_diffExp_diffSplice_HiPip)

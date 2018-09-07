Culex\_GO\_Analysis\_Level\_3
================
Megan Fritz
September 7, 2018

Plots for Gene Ontology Classification.
---------------------------------------

##### Loading in data sets and libraries.

``` r
library(rlist)

setwd("~/Desktop/CulexRNAseq")

#Getting the top 9 MF categories
DGC_MF_grVpa <- read.table("./data/Culex_DGE_Blast2GO_output/GO_data_analyses/gravid_v_parous/DGC_MF_grVpa.txt", nrows = 9, sep = "\t", header = T)
DGC_MF_paVma <- read.table("./data/Culex_DGE_Blast2GO_output/GO_data_analyses/parous_v_males/DGC_MF_paVma.txt", nrows = 9, sep = "\t", header = T)
DGC_MF_paVpip <- read.table("./data/Culex_DGE_Blast2GO_output/GO_data_analyses/parous_v_pipiens/DGC_MF_paVpip.txt", nrows = 9, sep = "\t", header = T)
DGC_MF_grVpip <- read.table("./data/Culex_DGE_Blast2GO_output/GO_data_analyses/gravid_v_pipiens/DGC_MF_grVpip.txt", nrows = 9, sep = "\t", header = T)

data_MF <- list(DGC_MF_grVpa$X, DGC_MF_paVma$X, DGC_MF_paVpip$X, DGC_MF_grVpip$X)
tot <- c(625,1744,3174,3327)

#Getting the top 9 BP categories
DGC_BP_grVpa <- read.table("./data/Culex_DGE_Blast2GO_output/GO_data_analyses/gravid_v_parous/DGC_BP_grVpa.txt", nrows = 9, sep = "\t", header = T)
DGC_BP_paVma <- read.table("./data/Culex_DGE_Blast2GO_output/GO_data_analyses/parous_v_males/DGC_BP_paVma.txt", nrows = 9, sep = "\t", header = T)
DGC_BP_paVpip <- read.table("./data/Culex_DGE_Blast2GO_output/GO_data_analyses/parous_v_pipiens/DGC_BP_paVpip.txt", nrows = 9, sep = "\t", header = T)
DGC_BP_grVpip <- read.table("./data/Culex_DGE_Blast2GO_output/GO_data_analyses/gravid_v_pipiens/DGC_BP_grVpip.txt", nrows = 9, sep = "\t", header = T)

data_BP <- list(DGC_BP_grVpa$X, DGC_BP_paVma$X, DGC_BP_paVpip$X, DGC_BP_grVpip$X)


#Getting the top 9 CC categories
DGC_CC_grVpa <- read.table("./data/Culex_DGE_Blast2GO_output/GO_data_analyses/gravid_v_parous/DGC_CC_grVpa.txt", nrows = 9, sep = "\t", header = T)
DGC_CC_paVma <- read.table("./data/Culex_DGE_Blast2GO_output/GO_data_analyses/parous_v_males/DGC_CC_paVma.txt", nrows = 9, sep = "\t", header = T)
DGC_CC_paVpip <- read.table("./data/Culex_DGE_Blast2GO_output/GO_data_analyses/parous_v_pipiens/DGC_CC_paVpip.txt", nrows = 9, sep = "\t", header = T)
DGC_CC_grVpip <- read.table("./data/Culex_DGE_Blast2GO_output/GO_data_analyses/gravid_v_pipiens/DGC_CC_grVpip.txt", nrows = 9, sep = "\t", header = T)

data_CC <- list(DGC_CC_grVpa$X, DGC_CC_paVma$X, DGC_CC_paVpip$X, DGC_CC_grVpip$X)
```

##### Getting MF GO percentages by comparison

``` r
inv_misc_cat_MF <- lapply(data_MF, sum)
misc_cat_MF <- tot - (unlist(inv_misc_cat_MF))

data2_MF <- as.data.frame(data_MF, col.names = c("DGC_MF_grVpa", "DGC_MF_paVma", "DGC_MF_paVpip", "DGC_MF_grVpip"))

Updata_MF <- rbind(data2_MF, misc_cat_MF)

Percentages_MF <- data.frame((t(apply(Updata_MF, 1, "/", tot))*100))
Rounded_Percentages_MF <- round(Percentages_MF, 1)
```

##### Getting BP GO percentages by comparison

``` r
inv_misc_cat_BP <- lapply(data_BP, sum)
misc_cat_BP <- tot - (unlist(inv_misc_cat_BP))

data2_BP <- as.data.frame(data_BP, col.names = c("DGC_BP_grVpa", "DGC_BP_paVma", "DGC_BP_paVpip", "DGC_BP_grVpip"))

Updata_BP <- rbind(data2_BP, misc_cat_BP)

Percentages_BP <- data.frame((t(apply(Updata_BP, 1, "/", tot))*100))
Rounded_Percentages_BP <- round(Percentages_BP, 1)
```

##### Getting CC GO percentages by comparison

``` r
inv_misc_cat_CC <- lapply(data_CC, sum)
misc_cat_CC <- tot - (unlist(inv_misc_cat_CC))

data2_CC <- as.data.frame(data_CC, col.names = c("DGC_CC_grVpa", "DGC_CC_paVma", "DGC_CC_paVpip", "DGC_CC_grVpip"))

Updata_CC <- rbind(data2_CC, misc_cat_CC)

Percentages_CC <- data.frame((t(apply(Updata_CC, 1, "/", tot))*100))
Rounded_Percentages_CC <- round(Percentages_CC, 1)
```

##### First Making Pie Charts for Molecular Function

###### Here, I used the following abbreviations for my pie chart labels:

###### "structural constituent of ribosome" = "struc const ribosome"

###### "serine-type endopeptidase activity", "ser-type endopep act"

``` r
lbls <- c(as.character(DGC_MF_grVpa$GO), "Miscellaneous")
lbls <- gsub("structural constituent of ribosome", "struc const ribosome", lbls)
lbls <- paste(lbls, Rounded_Percentages_MF$DGC_MF_grVpa) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels 
grVpa_pie_MF <- pie(Rounded_Percentages_MF$DGC_MF_grVpa, labels =lbls, main = "Top Molecular Function GO Level 3 Categories for Gravid CAL1 versus Parous CAL1 Comparison", radius = 0.9, col = rainbow(length(lbls)), cex = 0.75, cex.main = 0.75, tck=.2, clockwise = T)
```

![](Culex_GO_Figs_09072018_files/figure-markdown_github/MF%20GO%20Figs-1.png)

``` r
lbls <- c(as.character(DGC_MF_paVma$GO), "Miscellaneous")
lbls <- gsub("structural constituent of ribosome", "struc const ribosome", lbls)
lbls <- paste(lbls, Rounded_Percentages_MF$DGC_MF_paVma) 
lbls <- paste(lbls,"%",sep="") 
paVma_pie_MF <- pie(Rounded_Percentages_MF$DGC_MF_paVma, labels =lbls, main = "Top Molecular Function GO Level 3 Categories for Parous CAL1 versus Male CAL1 Comparison", radius = 0.9, col = rainbow(length(lbls)), cex = 0.75, cex.main = 0.75, tck=.2, clockwise = T)
```

![](Culex_GO_Figs_09072018_files/figure-markdown_github/MF%20GO%20Figs-2.png)

``` r
lbls <- c(as.character(DGC_MF_paVpip$GO), "Miscellaneous")
lbls <- gsub("structural constituent of ribosome", "struc const ribosome", lbls)
lbls <- gsub("serine-type endopeptidase activity", "ser-type endopep act", lbls)
lbls <- paste(lbls, Rounded_Percentages_MF$DGC_MF_paVpip) 
lbls <- paste(lbls,"%",sep="") 
paVpip_pie_MF <- pie(Rounded_Percentages_MF$DGC_MF_paVpip, labels =lbls, main = "Top Molecular Function GO Level 3 Categories for Parous CAL1 versus IL2 Comparison", radius = 0.9, col = rainbow(length(lbls)), cex = 0.75, cex.main = 0.75, tck=.2, clockwise = T)
```

![](Culex_GO_Figs_09072018_files/figure-markdown_github/MF%20GO%20Figs-3.png)

``` r
lbls <- c(as.character(DGC_MF_grVpip$GO), "Miscellaneous")
lbls <- gsub("structural constituent of ribosome", "struc const ribosome", lbls)
lbls <- gsub("serine-type endopeptidase activity", "ser-type endopep act", lbls)
lbls <- paste(lbls, Rounded_Percentages_MF$DGC_MF_grVpip) 
lbls <- paste(lbls,"%",sep="") 
grVpip_pie_MF <- pie(Rounded_Percentages_MF$DGC_MF_grVpip, labels =lbls, main = "Top Molecular Function GO Level 3 Categories for Gravid CAL1 versus IL2 Comparison", radius = 0.9, col = rainbow(length(lbls)), cex = 0.75, cex.main = 0.75, tck=.2, clockwise = T)
```

![](Culex_GO_Figs_09072018_files/figure-markdown_github/MF%20GO%20Figs-4.png)

##### Biological Process GO Comparison Charts

###### Here, I used the following abbreviations for my pie chart labels:

###### "regulation of transcription, DNA-templated" = "reg transcript, DNA-temp"

###### "carbohydrate metabolic process", "carb metabol proc"

``` r
lbls <- c(as.character(DGC_BP_grVpa$GO), "Miscellaneous")
lbls <- gsub("regulation of transcription, DNA-templated", "reg transcript, DNA-temp", lbls)
lbls <- paste(lbls, Rounded_Percentages_BP$DGC_BP_grVpa) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels 
grVpa_pie_BP <- pie(Rounded_Percentages_BP$DGC_BP_grVpa, labels =lbls, main = "Top Biological Process GO Level 3 Categories for Gravid CAL1 versus Parous CAL1 Comparison", radius = 0.9, cex = 0.75, cex.main = 0.75, tck=.2, clockwise = T)
```

![](Culex_GO_Figs_09072018_files/figure-markdown_github/BP%20GO%20Figs-1.png)

``` r
lbls <- c(as.character(DGC_BP_paVma$GO), "Miscellaneous")
lbls <- gsub("regulation of transcription, DNA-templated", "reg transcript, DNA-temp", lbls)
lbls <- gsub("carbohydrate metabolic process", "carb metabol proc", lbls)
lbls <- paste(lbls, Rounded_Percentages_BP$DGC_BP_paVma) 
lbls <- paste(lbls,"%",sep="") 
paVma_pie_BP <- pie(Rounded_Percentages_BP$DGC_BP_paVma, labels =lbls, main = "Top Biological Process GO Level 3 Categories for Parous CAL1 versus Male CAL1 Comparison", radius = 0.9, cex = 0.75, cex.main = 0.75, tck=.2, clockwise = T)
```

![](Culex_GO_Figs_09072018_files/figure-markdown_github/BP%20GO%20Figs-2.png)

``` r
lbls <- c(as.character(DGC_BP_paVpip$GO), "Miscellaneous")
lbls <- gsub("regulation of transcription, DNA-templated", "reg transcript, DNA-temp", lbls)
lbls <- gsub("carbohydrate metabolic process", "carb metabol proc", lbls)
lbls <- paste(lbls, Rounded_Percentages_BP$DGC_BP_paVpip) 
lbls <- paste(lbls,"%",sep="") 
paVpip_pie_BP <- pie(Rounded_Percentages_BP$DGC_BP_paVpip, labels =lbls, main = "Top Biological Process GO Level 3 Categories for Parous CAL1 versus IL2 Comparison", radius = 0.9,  cex = 0.75, cex.main = 0.75, tck=.2, clockwise = T)
```

![](Culex_GO_Figs_09072018_files/figure-markdown_github/BP%20GO%20Figs-3.png)

``` r
lbls <- c(as.character(DGC_BP_grVpip$GO), "Miscellaneous")
lbls <- gsub("regulation of transcription, DNA-templated", "reg transcript, DNA-temp", lbls)
lbls <- gsub("carbohydrate metabolic process", "carb metabol proc", lbls)
lbls <- paste(lbls, Rounded_Percentages_BP$DGC_BP_grVpip) 
lbls <- paste(lbls,"%",sep="") 
grVpip_pie_BP <- pie(Rounded_Percentages_BP$DGC_BP_grVpip, labels =lbls, main = "Top Biological Process GO Level 3 Categories for Gravid CAL1 versus IL2 Comparison", radius = 0.9, cex = 0.75, cex.main = 0.75, tck=.2, clockwise = T)
```

![](Culex_GO_Figs_09072018_files/figure-markdown_github/BP%20GO%20Figs-4.png)

##### Cellular Component GO Comparison Charts

###### Here, I used the following abbreviations for my pie chart labels:

###### "regulation of transcription, DNA-templated" = "reg transcript, DNA-temp"

###### "carbohydrate metabolic process", "carb metabol proc"

``` r
col_schema <- c("blue3", "blue2", "blue", "deepskyblue", "lightblue", "orange", "yellow3", "yellow2", "yellow", "lightyellow")

lbls <- c(as.character(DGC_CC_grVpa$GO), "Miscellaneous")
lbls <- gsub("regulation of transcription, DNA-templated", "reg transcript, DNA-temp", lbls)
lbls <- paste(lbls, Rounded_Percentages_CC$DGC_CC_grVpa) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels 
grVpa_pie_CC <- pie(Rounded_Percentages_CC$DGC_CC_grVpa, labels =lbls, main = "Top Cellular Component GO Level 3 Categories for Gravid CAL1 versus Parous CAL1 Comparison", radius = 1, cex = 0.75, cex.main = 0.75, tck=.2, init.angle = 30, col= col_schema)
```

![](Culex_GO_Figs_09072018_files/figure-markdown_github/CC%20GO%20Figs-1.png)

``` r
lbls <- c(as.character(DGC_CC_paVma$GO), "Miscellaneous")
lbls <- gsub("regulation of transcription, DNA-templated", "reg transcript, DNA-temp", lbls)
lbls <- gsub("carbohydrate metabolic process", "carb metabol proc", lbls)
lbls <- paste(lbls, Rounded_Percentages_CC$DGC_CC_paVma) 
lbls <- paste(lbls,"%",sep="") 
paVma_pie_CC <- pie(Rounded_Percentages_CC$DGC_CC_paVma, labels =lbls, main = "Top Cellular Component GO Level 3 Categories for Parous CAL1 versus Male CAL1 Comparison", radius = 1, cex = 0.75, cex.main = 0.75, tck=.2, init.angle = 30, col= col_schema)
```

![](Culex_GO_Figs_09072018_files/figure-markdown_github/CC%20GO%20Figs-2.png)

``` r
lbls <- c(as.character(DGC_CC_paVpip$GO), "Miscellaneous")
lbls <- gsub("regulation of transcription, DNA-templated", "reg transcript, DNA-temp", lbls)
lbls <- gsub("carbohydrate metabolic process", "carb metabol proc", lbls)
lbls <- paste(lbls, Rounded_Percentages_CC$DGC_CC_paVpip) 
lbls <- paste(lbls,"%",sep="") 
paVpip_pie_CC <- pie(Rounded_Percentages_CC$DGC_CC_paVpip, labels =lbls, main = "Top Cellular Component GO Level 3 Categories for Parous CAL1 versus IL2 Comparison", radius = 1,  cex = 0.75, cex.main = 0.75, tck=.2, init.angle = 30, col= col_schema)
```

![](Culex_GO_Figs_09072018_files/figure-markdown_github/CC%20GO%20Figs-3.png)

``` r
lbls <- c(as.character(DGC_CC_grVpip$GO), "Miscellaneous")
lbls <- gsub("regulation of transcription, DNA-templated", "reg transcript, DNA-temp", lbls)
lbls <- gsub("carbohydrate metabolic process", "carb metabol proc", lbls)
lbls <- paste(lbls, Rounded_Percentages_CC$DGC_CC_grVpip) 
lbls <- paste(lbls,"%",sep="") 
grVpip_pie_CC <- pie(Rounded_Percentages_CC$DGC_CC_grVpip, labels =lbls, main = "Top Cellular Component GO Level 3 Categories for Gravid CAL1 versus IL2 Comparison", radius = 1, cex = 0.75, cex.main = 0.75, tck=.2, init.angle = 30, col= col_schema)
```

![](Culex_GO_Figs_09072018_files/figure-markdown_github/CC%20GO%20Figs-4.png)

Culex\_GO\_Analysis\_Level\_3
================
Megan Fritz
September 7, 2018

Starting with Plots for Molecular Function Gene Ontology Classification.
------------------------------------------------------------------------

##### Loading in data sets and libraries.

``` r
library(rlist)

setwd("~/Desktop/CulexRNAseq")

#Getting the top 9 categories
DGC_MF_grVpa <- read.table("./data/Culex_DGE_Blast2GO_output/GO_data_analyses/gravid_v_parous/DGC_MF_grVpa.txt", nrows = 9, sep = "\t", header = T)
DGC_MF_paVma <- read.table("./data/Culex_DGE_Blast2GO_output/GO_data_analyses/parous_v_males/DGC_MF_paVma.txt", nrows = 9, sep = "\t", header = T)
DGC_MF_paVpip <- read.table("./data/Culex_DGE_Blast2GO_output/GO_data_analyses/parous_v_pipiens/DGC_MF_paVpip.txt", nrows = 9, sep = "\t", header = T)
DGC_MF_grVpip <- read.table("./data/Culex_DGE_Blast2GO_output/GO_data_analyses/gravid_v_pipiens/DGC_MF_grVpip.txt", nrows = 9, sep = "\t", header = T)

data <- list(DGC_MF_grVpa$X, DGC_MF_paVma$X, DGC_MF_paVpip$X, DGC_MF_grVpip$X)
tot <- c(625,1744,3174,3327)
```

##### Getting GO percentages by comparison

``` r
inv_misc_cat <- lapply(data, sum)
misc_cat <- tot - (unlist(inv_misc_cat))

data2 <- as.data.frame(data, col.names = c("DGC_MF_grVpa", "DGC_MF_paVma", "DGC_MF_paVpip", "DGC_MF_grVpip"))

Updata <- rbind(data2, misc_cat)

Percentages <- data.frame((t(apply(Updata, 1, "/", tot))*100))
Rounded_Percentages <- round(Percentages, 1)
```

##### Making Pie Charts for Molecular Function

``` r
lbls <- c(as.character(DGC_MF_grVpa$GO), "Miscellaneous")
lbls <- gsub("structural constituent of ribosome", "struc const ribosome", lbls)
lbls <- paste(lbls, Rounded_Percentages$DGC_MF_grVpa) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels 
grVpa_pie <- pie(Rounded_Percentages$DGC_MF_grVpa, labels =lbls, main = "Top Molecular Function GO Level 3 Categories for Gravid CAL1 versus Parous CAL1 Comparison", radius = 0.9, col = rainbow(length(lbls)), cex = 0.75, cex.main = 0.75, tck=.2, clockwise = T)
```

![](Culex_GO_Figs_09072018_files/figure-markdown_github/MF%20GO%20Figs-1.png)

``` r
lbls <- c(as.character(DGC_MF_paVma$GO), "Miscellaneous")
lbls <- gsub("structural constituent of ribosome", "struc const ribosome", lbls)
lbls <- paste(lbls, Rounded_Percentages$DGC_MF_paVma) 
lbls <- paste(lbls,"%",sep="") 
paVma_pie <- pie(Rounded_Percentages$DGC_MF_paVma, labels =lbls, main = "Top Molecular Function GO Level 3 Categories for Parous CAL1 versus Male CAL1 Comparison", radius = 0.9, col = rainbow(length(lbls)), cex = 0.75, cex.main = 0.75, tck=.2, clockwise = T)
```

![](Culex_GO_Figs_09072018_files/figure-markdown_github/MF%20GO%20Figs-2.png)

``` r
lbls <- c(as.character(DGC_MF_paVpip$GO), "Miscellaneous")
lbls <- gsub("structural constituent of ribosome", "struc const ribosome", lbls)
lbls <- gsub("serine-type endopeptidase activity", "ser-type endopep act", lbls)
lbls <- paste(lbls, Rounded_Percentages$DGC_MF_paVpip) 
lbls <- paste(lbls,"%",sep="") 
paVpip_pie <- pie(Rounded_Percentages$DGC_MF_paVpip, labels =lbls, main = "Top Molecular Function GO Level 3 Categories for Parous CAL1 versus IL2 Comparison", radius = 0.9, col = rainbow(length(lbls)), cex = 0.75, cex.main = 0.75, tck=.2, clockwise = T)
```

![](Culex_GO_Figs_09072018_files/figure-markdown_github/MF%20GO%20Figs-3.png)

``` r
lbls <- c(as.character(DGC_MF_grVpip$GO), "Miscellaneous")
lbls <- gsub("structural constituent of ribosome", "struc const ribosome", lbls)
lbls <- gsub("serine-type endopeptidase activity", "ser-type endopep act", lbls)
lbls <- paste(lbls, Rounded_Percentages$DGC_MF_grVpip) 
lbls <- paste(lbls,"%",sep="") 
grVpip_pie <- pie(Rounded_Percentages$DGC_MF_grVpip, labels =lbls, main = "Top Molecular Function GO Level 3 Categories for Gravid CAL1 versus IL2 Comparison", radius = 0.9, col = rainbow(length(lbls)), cex = 0.75, cex.main = 0.75, tck=.2, clockwise = T)
```

![](Culex_GO_Figs_09072018_files/figure-markdown_github/MF%20GO%20Figs-4.png)
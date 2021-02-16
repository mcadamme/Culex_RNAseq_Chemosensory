Plot showing library strandedness
================
Megan Fritz
written Jan. 8, 2021

## Library strandedness and htseq-count

The Illumina TruSeq Stranded mRNA Sample Prep Kit preserves directional information associated with first-strand synthesis. Proper RNA-seq analysis with a stranded dataset requires careful attention to parameter settings in the read counting step of the bioinformatic analysis (Srinivasan et al. 2020).

Here, we are giving an example from our own dataset, using the criteria of Srinivasan et al. 2020 to show that we should be using the --stranded=reverse parameter setting in htseq-count.

``` r
strand_rev <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/DGE_GenAligned_SamFiles/highQual_exon/M2-1_S5_htseq", header = F) #ran with htseq -s reverse
strand_yes <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/DGE_GenAligned_SamFiles/Str_Yes/M2-1_S5_htseq", header = F) #ran with htseq -s yes


#merging dataset for plot
merged <- merge(strand_rev, strand_yes, by = "V1")
```

``` r
#looking at the ambiguous and no_feature counts
head(merged, n=2) #strand_rev has higher ambiguous alignments
```

    ##                       V1    V2.x    V2.y
    ## 1 __alignment_not_unique 2822659 2822659
    ## 2            __ambiguous  188782     254

``` r
tail(merged, n=3) #but strand_yes has higher no_feature counts.  Right away, this indicates strand_rev is better.
```

    ##                    V1    V2.x     V2.y
    ## 19796    __no_feature 4636071 17259522
    ## 19797   __not_aligned       0        0
    ## 19798 __too_low_aQual       0        0

``` r
merged_red <- merged[c(-1,-2, -19796, -19797, -19798),]
str(merged_red)
```

    ## 'data.frame':    19793 obs. of  3 variables:
    ##  $ V1  : Factor w/ 19798 levels "__alignment_not_unique",..: 3 4 5 6 7 8 9 10 11 12 ...
    ##  $ V2.x: int  0 146 0 160 215 3164 0 0 1 13 ...
    ##  $ V2.y: int  0 4 0 0 1 0 0 4 0 36 ...

``` r
#getting sums of counts to compare them - we see that there are many more reads counted toward features using strand_rev
colSums(merged_red[,c(2,3)])
```

    ##     V2.x     V2.y 
    ## 12590760   155837

Srinivasan et al. 2020 state a plot of stranded ‘Reverse’ against ‘Yes’ read counts may be useful if cDNA library construction in unknown (but we knew ours). Stranded ‘Reverse’ would be the correct option if the majority of the genes occur in a steep vertical line, while stranded ‘Yes’ would be correct if the majority of the genes were oriented in a horizontal line.

``` r
plot(merged_red$V2.y, merged_red$V2.x, ylim = c(0,5000), xlim = c(0,5000), ylab = "stranded_rev", xlab = "stranded_yes", title(main= "Feature counts by stranded parameter"))#axes truncated for better visualization
```

![](strandedness_and_htseq_files/figure-markdown_github/generating%20plot%20to%20look%20at%20counts%20per%20feature-1.png)

``` r
#Even the max count value is higher for stranded_rev
max(merged_red$V2.x)
```

    ## [1] 1482649

``` r
max(merged_red$V2.y)
```

    ## [1] 29486

# Choosing gene versus exon counts

``` r
strand_revG <- read.table("/media/megan/EE9C16C89C168AEB/Noreuil/trimmed_pairs/DGE_GenAligned_SamFiles/highQual_gene/M2-1_S5_htseq", header = F) #ran with htseq -s reverse

#merging dataset for decision
merged <- merge(strand_rev, strand_revG, by = "V1")

merged$diff <- merged$V2.x-merged$V2.y
head(merged)
```

    ##                       V1    V2.x    V2.y diff
    ## 1 __alignment_not_unique 2822659 2822659    0
    ## 2            __ambiguous  188782  188782    0
    ## 3             CPIJ000001       0       0    0
    ## 4             CPIJ000002     146     146    0
    ## 5             CPIJ000003       0       0    0
    ## 6             CPIJ000004     160     160    0

``` r
sum(merged$diff)
```

    ## [1] 0

``` r
#no diff here, so using exon dataset.
```

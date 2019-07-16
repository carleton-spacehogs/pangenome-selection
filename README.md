Generating Figures for Moulana et al.
================
Alief Moulana
7/10/2019

``` r
knitr::opts_chunk$set(echo = TRUE)


#load all necessary packages
library(reshape2)
library(ggplot2)
library(tidyverse)
```

    ## ── Attaching packages ──────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ tibble  2.1.3     ✔ purrr   0.3.2
    ## ✔ tidyr   0.8.3     ✔ dplyr   0.8.3
    ## ✔ readr   1.3.1     ✔ stringr 1.4.0
    ## ✔ tibble  2.1.3     ✔ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(tidyr)
library(dplyr)
library(stringr)
library(factoextra)
```

    ## Welcome! Related Books: `Practical Guide To Cluster Analysis in R` at https://goo.gl/13EFCZ

``` r
#load all the data needed

pangenome <- as.matrix(read.table('txt_files/similarity_matrix.txt', 
                                  quote = "",sep='\t',header=FALSE)) 
#load the similarity matrix data

cluster<- read.table('txt_files/Sulfurovum_gene_clusters_summary_cleaned.txt',sep='\t',quote="",
                     header=TRUE) 
#load the summarized sulfurovum pangenome profile from anvi'o
```

## Taxonomy distribution and pangenomic profile

We performed the common metagenomic pipeline of assembly, mapping,
binning, and annotation (see Methods). Then, we generated Supplementary
Figure 1 on a spread sheet using the taxonomy data we obtained from
“phylosift”. Because Sulfurovum was the most abundant genus in both
Mid Cayman Rise and Axial genomes, we used Sulfurovum MAGs in the
subsequent pangenome analyses following “anvi’o” pangenome workflow.
Then we captured Figure 1A using command `anvi-display-pan` with genomes
storage and pangenome databases generated from the previous steps.
Moreover, we obtained the summary of the pangenome profile (the cluster
file) by running the command `anvi-summarize` for the PAN DB. We also
analyzed how the pattern of total and core genome accumulation as more
MAGs considered. We visualize the result from the chunk below
(Supplementary Figure 2).

``` r
cluster.basic.count <- unique(data.frame(Group=cluster$gene_cluster_id,
                                  Number=cluster$num_genomes_gene_cluster_has_hits,
                                  Genome=cluster$genome_name)) 
#create df, remove duplicates
cluster.basic.count$Group <- as.character(cluster.basic.count$Group)

#reshape df
cluster.basic.reshaped <- reshape(cluster.basic.count,v.names="Number",
                                  timevar="Genome",idvar=c("Group"),direction="wide")
cluster.basic.matrix <- as.matrix(cluster.basic.reshaped[,2:23]) #turn into matrix
cluster.basic.matrix[] <- c(TRUE,FALSE)[(is.na(cluster.basic.matrix) )+ 1] #binary

total.each.genome <- apply(cluster.basic.matrix,2,sum) #numb of gene groups per genome

#create a function returning whether or not the sum of a row up until i equals i
count.sum <- function(row,i){
  sum(row[1:i]) == i
}

#create a vector that count the number of core genes taking into account i genomes
core.now <- c()
for(j in 1:22){
  core.now[j] <- sum(apply(cluster.basic.matrix,1,count.sum,j))
}

#create a function to perform the analysis below
count.sum.2 <- function(numb,x){
  var <- sum(apply(as.matrix(cluster.basic.matrix[,1:x-1]),1,sum) ==0 & 
               cluster.basic.matrix[,x] == 1)
  return(numb+var)
}

#create a vector that expands number of genes as more genomes sequenced
total.now <- c()
total.now[1] <- total.each.genome[1]
count <- total.now[1]
for(k in 2:22){
  total.now[k] <- count.sum.2(count,k)
  count <- total.now[k]
}
#merge the core now and total now
count.groups <- data.frame(Number=c(1:22),Core=core.now,Total=total.now)
ggplot(count.groups,aes(x=Number))+
  geom_line(aes(y=Core,colour="Number of Core Genes"))+
  geom_line(aes(y=Total/5,colour="Number of Total Genes"))+
  scale_y_continuous(sec.axis = sec_axis(~.*5, name = "Number of Total Genes"))+
  xlab("Number of MAGs Recovered") +
  ylab("Number of Core Genes") +
  labs(colour = "Parameter")+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        axis.title=element_text(size=14),
        legend.position = c(0.3, 0.85))
```

![](pangenome-selection_files/figure-gfm/pangenome.count-1.png)<!-- -->

Moreover, we also analyze the similarity in gene content among MAGs. We
first run `similarity_matrix.py` which takes the cluster file as an
input and produces the file `similarity_matrix.txt`. The program
produces a 22-by-22 matrix, where each entry \(ij\) represents the
proportion of genes in the \(i\)th MAG that is also contained by the
\(j\)th MAG. We then visualize the matrix and find the distance among
the genomes in the following chunk.
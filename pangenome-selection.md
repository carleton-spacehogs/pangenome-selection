Processing Data and Generating Figures for Moulana et al.
================
Alief Moulana

``` r
knitr::opts_chunk$set(echo = TRUE)


#load all necessary packages
library(reshape2)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(stringr)
library(factoextra)
library(RColorBrewer)
library(tm)

#load all the data needed

pangenome <- as.matrix(read.table('aux_files/similarity_matrix.txt', 
                                  quote = "",sep='\t',header=FALSE)) 
#load the similarity matrix data

cluster<- read.table('aux_files/Sulfurovum_gene_clusters_summary_cleaned.txt',sep='\t',quote="",header=TRUE) 
#load the summarized sulfurovum pangenome profile from anvi'o
```

## Taxonomy distribution and pangenomic profile

We performed the common metagenomic pipeline of assembly, mapping,
binning, and annotation (see Methods; refer to citations). We then
analyzed the taxonomic composition in each region using Excel
(Supplementary Figure 1; from `Taxonomy_data.xlsx`). Because Sulfurovum
was the most abundant genus in both Mid Cayman Rise and Axial genomes,
we used Sulfurovum MAGs in the subsequent pangenome analyses following
“anvi’o” pangenome workflow. Then we captured Figure 1A using command
`anvi-display-pan` with genomes storage and pangenome databases
generated from the previous steps. Moreover, we obtained the summary of
the pangenome profile (the cluster file
`Sulfurovum_gene_clusters_summary_cleaned.txt`) by running the command
`anvi-summarize` for the PAN DB. We also analyzed how the pattern of
total and core genome accumulation as more MAGs considered. We visualize
the result from the chunk below (Supplementary Figure 2).

``` r
cluster.basic.count <- unique(data.frame(Group=cluster$gene_cluster_id,
                                  Number=cluster$num_genomes_gene_cluster_has_hits,
                                  Genome=cluster$genome_name)) 
# create df, remove duplicates
cluster.basic.count$Group <- as.character(cluster.basic.count$Group)

#reshape df
cluster.basic.reshaped <- reshape(cluster.basic.count,v.names="Number",
                                  timevar="Genome",idvar=c("Group"),direction="wide")
cluster.basic.matrix <- as.matrix(cluster.basic.reshaped[,2:23]) #turn into matrix
cluster.basic.matrix[] <- c(TRUE,FALSE)[(is.na(cluster.basic.matrix) )+ 1] #binary

total.each.genome <- apply(cluster.basic.matrix,2,sum) #numb of gene groups per genome

# create a function returning whether or not the sum of a row up until i equals i
count.sum <- function(row,i){
  sum(row[1:i]) == i
}

# create a vector that count the number of core genes taking into account i genomes
core.now <- c()
for(j in 1:22){
  core.now[j] <- sum(apply(cluster.basic.matrix,1,count.sum,j))
}

# create a function to perform the analysis below
count.sum.2 <- function(numb,x){
  var <- sum(apply(as.matrix(cluster.basic.matrix[,1:x-1]),1,sum) ==0 & 
               cluster.basic.matrix[,x] == 1)
  return(numb+var)
}

# create a vector that expands number of genes as more genomes sequenced
total.now <- c()
total.now[1] <- total.each.genome[1]
count <- total.now[1]
for(k in 2:22){
  total.now[k] <- count.sum.2(count,k)
  count <- total.now[k]
}
# merge the core now and total now
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

<img src="pangenome-selection_files/figure-gfm/pangenome.count-1.png" width="60%" height="60%" />

Moreover, we also analyze the similarity in gene content among MAGs. We
first run `similarity_matrix.py` which takes the cluster file as an
input and produces the file `similarity_matrix.txt`. The program
produces a 22-by-22 matrix, where each entry \(ij\) represents the
proportion of genes in the \(i\)th MAG that is also contained by the
\(j\)th MAG. We then visualize the matrix and find the distance among
the genomes in the following chunk.

\<\<\<\<\<\<\< HEAD

``` r
# create genome names
twenty_two <- as.character(c(1:22))
genome_numbering <- ifelse(nchar(twenty_two) == 1,
                           paste(as.character(0),twenty_two,sep=""),
                           as.character(twenty_two))
genome_names <- paste("Sulfurovum",genome_numbering,sep="_")

# then manipulate the data for the heatmap where each square represents the percentage of genes in genome 1 contained in genome 2
heat_map.data <- data.frame(pangenome)
colnames(heat_map.data) <- genome_names
heat_map.data$genome <- genome_names

# melt the data from a large matrix to pairwise rows
heat_map.data.melt <- melt(heat_map.data, id=c("genome"))
colnames(heat_map.data.melt)<-c("MAG_1", "MAG_2", "Contained")

# cluster the data based on the genome content overlap
a <- hclust(as.dist(1-heat_map.data[, -23]))
a$labels <- genome_numbering
library(factoextra)
dend_plot <- fviz_dend(a)
dend_plot # this gets the dendrogram
```

<img src="pangenome-selection_files/figure-gfm/similarity-figure-1.png" width="60%" height="60%" />

``` r
order <- a$order

# order the MAGs based on the clustering
heat_map.data.melt$MAG_1 <- factor(heat_map.data.melt$MAG_1, levels=ifelse(order < 10, paste0("Sulfurovum_0", order), paste0("Sulfurovum_", order)))
heat_map.data.melt$MAG_2 <- factor(heat_map.data.melt$MAG_2, levels=ifelse(order < 10, paste0("Sulfurovum_0", order), paste0("Sulfurovum_", order)))

# Fig 1B
# plot the heatmap
ggplot(data = heat_map.data.melt, aes(x = MAG_1, y = MAG_2)) +
  geom_tile(aes(fill = Contained))+
  scale_fill_gradientn(colours = c("white", "lightblue", "darkblue"), values = c(0,0.5,1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

<img src="pangenome-selection_files/figure-gfm/similarity-figure-2.png" width="60%" height="60%" />

Because our pangenome is created upon metagenomes, we need to be careful
when declaring a gene to be in the core genome. We simulate the
probability that a gene found in \(n\) genomes to be actually a core
genome, exlcluding genes found in 22 genomes.

``` r
# summary file created based on `anvi-summarize` output
summary_data <- read.table('aux_files/sulfurovum_summary.txt',
                           sep='\t',header=TRUE,quote="",na.strings = "",fill = TRUE )

# function that takes in a number i and outputs the probability that a gene is a missing gene in the ith genome
gene.chance<-function(i){
  data.selected <- summary_data[i,] # look at the ith row of data
  mean.complete <- 100*data.selected$number_genes/data.selected$completion # approximate total number of genes
  sd.complete <- -(mean.complete-data.selected$number_genes)/(qnorm(data.selected$redundancy/100)) # approximate standard deviation
  missing.genes <- rnorm(1,mean.complete,sd.complete) - data.selected$number_genes # total number of genes supossedly under normal distribution
  return(missing.genes/(10263-data.selected$number_genes)) # probability a gene outside this genome is a missing gene in this genome
}

# run the simulation (pretty long)
collection <- data.frame(MAGs=as.factor(c()),
                         Probability=as.numeric(c()))
N<-100
for (i in 1:21){
  for (j in 1:N){
    sampled <- as.matrix(sample(1:22, i, replace=F)) # sample i genomes without replacement
    missing <- apply(sampled,1,gene.chance) # row by row apply function above
    collection <- rbind(collection,c(22-i,prod(missing))) # ML by calculate the product
  }
}

colnames(collection) = c("MAGs","Probability")
collection$MAGs <- as.factor(collection$MAGs)


# Supplement figure for probability
ggplot(collection,mapping=aes(x=MAGs,y=Probability))+
  scale_y_log10()+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        axis.title=element_text(size=14),
        )+
  labs(x="Number of MAGs",
       y="Probability of Core")+
  geom_hline(yintercept=0.05,color="red")
```

<img src="pangenome-selection_files/figure-gfm/simulation-contained-1.png" width="60%" height="60%" />

## The distinct functions between high- and low-frequency genes

We then analyzed the gene annotations and studied the annotation
distribution across gene frequency. To get the count table for category
count across number of MAGs, which is in the output file
`COG_category.txt`, run `get_functions.py`.

``` r
# data and manipulation
COG <- na.omit(read.table('aux_files/Sulfurovum_function_list.txt',
                          sep='\t',header=TRUE))
COG <- arrange(COG,Category)
Categories <- read.table('aux_files/COG_explanation.txt',sep='\t',header=FALSE)
colnames(COG) <- c('Category',as.character(c(1:22)))
COG_melt <- melt(COG, id=c("Category"))

# some color palette functions
colourCount = length(unique(COG$Category))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

load("permutation_test.RData") #get the permuation test result
colnames(p.df) <- c(1:22)
p.df$categories <- factor(COG$Category, levels=COG$Category[order])
p.df$explanation <- factor(Categories$V2, levels=Categories$V2[order])
p.df.melt <- melt(p.df, id=c("categories","explanation"))
#p.values significy the frequency such that the observed values are greater than or equal to the simulated values

#get distance matrix between categories
n.categories <- 25
n.genomes <- 22

distance.matrix <- matrix(0,ncol=n.categories,nrow=n.categories)

for(i in 1:n.categories){
  for(j in 1:n.categories){
    current <- 0
    for(k in 1:n.genomes){
      current <- current + (p.matrix[i,k]-p.matrix[j,k])^2
    }
    distance.matrix[i,j]<-current
  }
} #calculate the sum of squared distance in p-values between categories
colnames(distance.matrix) <- COG$Category
rownames(distance.matrix) <- COG$Category
p.cluster <- hclust(as.dist(distance.matrix)) #hierarchical clustering
library(factoextra)
p.dendrogram <- fviz_dend(p.cluster,show_labels = TRUE)
p.dendrogram #get the dendrogram
```

<img src="pangenome-selection_files/figure-gfm/unnamed-chunk-1-1.png" width="60%" height="60%" />

``` r
category.vector <- as.vector(COG$Category) 
explanation.vector <- as.vector(p.df$explanation) #explanation vector
order <- p.cluster$order
p.df.melt$explanation <- factor(p.df$explanation, levels=explanation.vector[order])
COG_melt$Category <- factor(COG_melt$Category, levels=COG_melt$Category[order])

# Fig 2A that shows the trend of COG categories across gene frequency
ggplot(COG_melt,aes(x=variable, y=value, fill=Category))+
  geom_bar( stat="identity", position="fill")+
  xlab("Number of MAGs") +
  ylab("Gene Category Count") +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        axis.title=element_text(size=14),
        legend.position = "right")
```

<img src="pangenome-selection_files/figure-gfm/unnamed-chunk-1-2.png" width="60%" height="60%" />

P-value figure as the following.

``` r
# for p_value (Supplement)
ggplot(data = p.df.melt, aes(x = variable, y = categories)) +
  geom_tile(data=subset(p.df.melt,value < 0.001),aes(fill = (value*10000)-20))+
  geom_tile(data=subset(p.df.melt,value > 0.999),aes(fill = ((value-1)*10000)+20))+
  geom_tile(data=subset(p.df.melt,value >= 0.001 & value <= 0.999),aes(fill = (value-0.5)/100))+
  scale_fill_gradient2(limits=c(-20, 20),guide = "colourbar",
                        low = "darkred", mid= "white", high = "darkblue",midpoint=0)+
  xlab("Number of MAGs (Frequency)") +
  guides(fill=guide_legend(title="P(Observed > Expected)"))+
  theme(plot.title = element_text(hjust = -0.4),
        plot.margin = rep(grid::unit(0.75,"in"),4))+
  theme(axis.text.x = element_text(size=12),
        axis.text.y=element_text(size=12),
        legend.text=element_text(size=12),
        legend.position="none",
        axis.title=element_text(size=14),
        title=element_text(size=14),
        axis.title.y=element_blank())
```

<img src="pangenome-selection_files/figure-gfm/unnamed-chunk-2-1.png" width="60%" />

We then analyzed specific categories

``` r
# melting and data manipulation
COG_prop <-  data.frame(prop.table(as.matrix(COG[,-1]), margin = 2))
COG_prop$Group <- COG$Category
COG_prop$Explanation <- Categories$V2
COG_prop$Group <- factor(COG_prop$Group, levels=COG_prop$Group[order])
COG_prop$Explanation <- factor(COG_prop$Explanation,
                               levels=COG_prop$Explanation[order])

COG_prop.melt <- melt(COG_prop,id=c("Group","Explanation"))
COG_prop.melt$variable <- rep(1:22,each=25)

# here we make sure that the color matches the color in Fig 2A
Color_select <- data.frame(Group= levels(COG_prop$Group), 
                          Color= getPalette(colourCount))
concerned.group <- c('J','E','M','P')
concerned.color <- as.vector(Color_select[Color_select$Group %in% concerned.group,]$Color)

# Fig 2B 
ggplot(COG_prop.melt[COG_prop.melt$Group %in% concerned.group,],aes(x=variable,y=value))+
  geom_bar(stat="identity")+
  theme_bw()+
  xlab("Number of MAGs")+
  ylab("Proportion")+
  theme(axis.text=element_text(size=12),
        legend.position="none",
        axis.title=element_text(size=14),
        strip.text = element_text(size=10))+
  facet_wrap(~Explanation,ncol=2)+
  aes(fill=Group)+
  scale_fill_manual(values=concerned.color) 
```

<img src="pangenome-selection_files/figure-gfm/COG-1.png" width="60%" height="60%" />

We then calculated the enrichment of a gene group in one region vs. the
other by assuming binomial distribution for each gene group. First, run
`parse_vent_function.py` to calculate the occurrence of each gene
cluster in either of the vent regions, which then gives an output of
`Sulfurovum_vent_count_per_cluster.txt`. Then, run
`compare_mcr_axial.py 12 8` to find the binomial CDF for gene clusters
found in Axial genomes only for cluster that are found in at least all
but one MAGs of either Axial or MCR, which generates
`Sulfurovum_Axial_enriched_compressed.txt`. We process the data
below.

``` r
enrichment<-na.omit(read.table('aux_files/Sulfurovum_Axial_enriched_compressed.txt',
                               sep='\t',header=TRUE,quote="",na.strings = "")) # load the data

# the jittery plot in Figure 3
ggplot(data = enrichment, mapping = aes(x=Category, y=Enrichment)) +
  geom_boxplot(alpha = 0,na.rm = TRUE) +
  geom_jitter(alpha = 0.3, aes(colour = Category),na.rm=TRUE,size=3)+
  geom_hline(yintercept=0.05, linetype=2, color = "red", size=1)+
  geom_hline(yintercept=0.95, linetype=2, color = "blue", size=1)+
  theme_bw()+
  ylab("Binomial CDF for Proportion in Axial")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none")
```

<img src="pangenome-selection_files/figure-gfm/mcr_vs_axial-1.png" width="60%" height="60%" />

``` r
# now only look at the P category to see the difference from others
enrichment.P <- subset(enrichment,Category == 'P') # extract this into Table 1 (for Top 15 with lowest CDF score)
enrichment.Non <- subset(enrichment,Category != 'P')
t.test(enrichment.P$Enrichment,enrichment.Non$Enrichment,alternative="less")
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  enrichment.P$Enrichment and enrichment.Non$Enrichment
    ## t = -4.4981, df = 45.083, p-value = 2.388e-05
    ## alternative hypothesis: true difference in means is less than 0
    ## 95 percent confidence interval:
    ##        -Inf -0.1070381
    ## sample estimates:
    ## mean of x mean of y 
    ## 0.2721944 0.4430042

## Signatures of non-neutral evolution in the pangenome

In this section, we first calculated the pN/pS ratios of each gene (not
gene cluster) in the pangenome. To do so, we used
`anvi-script-calculated-pn-ps-ratio` as described in Methods (we ran
`get-variability-profile-scv.py` and `get-pn-ps.sh`). Note that, we did
this step for each sample, so we need to create a summary file for all
samples using `summarize_pn_pn_sulfurovum.py` (have to be run with all
the processed data), giving us
`Sulfurovum_pn_ps_summary.txt`.

``` r
dat<-read.table('aux_files/Sulfurovum_pn_ps_summary.txt',sep='\t',header=TRUE,na.strings="Inf")
dat_counts<-na.omit(dat)
dat_counts<-dat_counts[is.finite(dat_counts$pN_pS),]
dat_counts$num_genomes <- factor(dat_counts$num_genomes)
dat_counts$type <- factor(ifelse(dat_counts$num_genomes == 22, "N=22", ifelse(dat_counts$num_genomes %in% 1:19,"0<N<20",NA))) #arbitrady core and accessory naming

# also look across categories by doing this:
load("dat_category.RData")

# gives Fig 4A
ggplot(data = dat_counts, mapping = aes(x=num_genomes, y=pN_pS)) +
  geom_boxplot(alpha = 0.3,na.rm = TRUE,outlier.size = 1,aes(colour=num_genomes)) +
  #geom_jitter(alpha = 0.3, aes(colour = num_genomes),na.rm=TRUE,size=1)+
  geom_hline(yintercept=1, linetype=2, color = "rosybrown", size=1)+
  xlab("Number of MAGs") +
  ylab("pN/pS") +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none")
```

<img src="pangenome-selection_files/figure-gfm/pn_ps-1.png" width="60%" height="60%" />

``` r
# we then find the Wilcoxon pairwise p-value for each Number of MAGs, then we also find the mean.
allthe.p <- c()
mean.p <- c()
sd.p <- c()
for(i in 1:22){
  dat_counts$type01 <- ifelse(dat_counts$num_genomes==i,'yes','no')
  mean.p[i] <- mean(dat_counts[dat_counts$num_genomes==i,]$pN_pS)
  sd.p[i] <- sd(dat_counts[dat_counts$num_genomes==i,]$pN_pS)
  wiwi <- wilcox.test(pN_pS ~ type01, data = dat_counts)
  allthe.p[i] <- wiwi$p.value
}
p_value.pN_pS <- data.frame(Number= as.factor(c(1:22)),P= allthe.p, Mean = mean.p,
                            Sd = sd.p)

# gives the supplementary figure for the p-value
ggplot(data=p_value.pN_pS,aes(x=Number,y=P))+
  geom_point(size=3)+
  xlab("Number of MAGs") +
  ylab("Wilcoxon P-value") +
  scale_y_log10()+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none")
```

<img src="pangenome-selection_files/figure-gfm/pn_ps-2.png" width="60%" height="60%" />

``` r
# gives the inset in Figure 4A for mean
ggplot(data=p_value.pN_pS,aes(x=Number,y=Mean))+
  geom_point(size=3)+
  xlab("Number of MAGs") +
  ylab("Mean of pN/pS") +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=24),
        legend.position = "none")
```

<img src="pangenome-selection_files/figure-gfm/pn_ps-3.png" width="60%" height="60%" />

``` r
# gives the supplementary for standard deviation
ggplot(data=p_value.pN_pS,aes(x=Number,y=Sd))+
  geom_point(size=3)+
  xlab("Number of MAGs") +
  ylab("Standard Deviation\nof pN/pS") +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none")
```

<img src="pangenome-selection_files/figure-gfm/pn_ps-4.png" width="60%" height="60%" />

``` r
# supplementary for genome name
ggplot(data = dat_counts, mapping = aes(x=gsub("Sulfurovum_","",genome_name), 
                                        y=pN_pS)) +
  geom_boxplot(alpha = 0.3,na.rm = TRUE,outlier.size = 1,aes(colour=genome_name)) +
  geom_hline(yintercept=1, linetype=2, color = "rosybrown", size=1)+
  xlab("MAG") +
  ylab("pN/pS") +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none")
```

<img src="pangenome-selection_files/figure-gfm/pn_ps-5.png" width="60%" height="60%" />

``` r
# gives Fig 4B
ggplot(data = na.omit(dat_counts), mapping = aes(x=type, y=pN_pS)) +
  geom_boxplot(alpha = 0.3,na.rm = TRUE,outlier.size = 1,aes(colour=type)) +
  geom_hline(yintercept=1, linetype=2, color = "rosybrown", size=1)+
  xlab("Number of MAGs") +
  ylab("pN/pS") +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none")
```

<img src="pangenome-selection_files/figure-gfm/pn_ps-6.png" width="60%" height="60%" />

``` r
# now the supplementary figure for data across categories
ggplot(data = dat_category, mapping = aes(x=category, y=pN_pS)) +
  geom_boxplot(alpha = 0.3,na.rm = TRUE, aes(colour = category),outlier.size=1) +
  #geom_jitter(alpha = 0.3, aes(colour = category),na.rm=TRUE,size=1)+
  geom_hline(yintercept=1, linetype=2, color = "rosybrown", size=1)+
  xlab("Category") +
  ylab("pN/pS") +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position = "none")
```

<img src="pangenome-selection_files/figure-gfm/pn_ps-7.png" width="60%" height="60%" />

Then, we load the data from the selective sweep study. We obtained this
data from running `pre_poisson_test.py` which utilizes the SNV profile
files generated from anvi’o script. For each sample, we got a Poisson
CDF score file, which we later collected using `collect_poisson.py`
which outputs `Sulfurovum_poisson_sweep.txt`.

``` r
# load, merge, and manipulate data
data<-read.table('aux_files/Sulfurovum_poisson_sweep.txt',sep='\t',header=TRUE)
data$numb <- factor(data$numb)
data <- na.omit(data)
data$pv_val_norm <- -log10(data$pv_val)
data <- arrange(data,unique_id)
data$pN_pS <- dat$pN_pS
data$order <- seq(1:length(data[,1]))
data$genome_name <- cluster$genome_name
data$cluster_id <- cluster$gene_cluster_id
data$category <- cluster$COG_CATEGORY
data$func <- cluster$COG_FUNCTION

# first, we still process the pN_pS data for the phosphate related genes
phosphate.clusters <- c("GC_00001177",
                       "GC_00001253",
                       "GC_00001293",
                       "GC_00001150",
                       "GC_00001267") #the interesting phosphate gene clusters
data.P <- na.omit(data[data$cluster_id %in% phosphate.clusters,])
data.P$unique_id <- as.factor(data.P$unique_id)

ggplot(data.P,mapping = aes(x=unique_id, y=pN_pS, fill=func)) +
  geom_bar(stat='identity') +
  theme_bw()+
  xlab("Gene") +
  ylab("pN/pS Ratio") +
  theme(axis.text=element_text(size=12),
        axis.text.x=element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        axis.title=element_text(size=14),
        legend.position = c(0.55,0.8))
```

<img src="pangenome-selection_files/figure-gfm/phosphate_pn_ps_and_sweep-1.png" width="60%" height="60%" />

Then, we performed the actual selective sweep analyses. The data for
total number of SNVs in each MAG is directly extracted from the SNV
profiled files provided by anvi’o.

``` r
# first, we consider the SNV density in each MAG
SNV_file <- read.csv("aux_files/normalized_snv.csv")
SNV_file$ID <- gsub('Sulfurovum_', '', as.vector(SNV_file$genome_id))

# Fig 5A
ggplot(SNV_file,mapping=aes(x=ID,y=norm_SNV,fill=vent_field))+
  geom_col()+
  xlab("Sulfurovum MAG") +
  ylab("Normalized SNV Count \n (#SNVs/kbp)") +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_blank(),
        legend.position=c(0.15,0.68))
```

<img src="pangenome-selection_files/figure-gfm/sweep-1.png" width="60%" height="60%" />

``` r
# now analyze the gene-specific sweeps
numb.numb <- nlevels(data$numb)
proportion <- rep(0,numb.numb)

# p-values < 1e-10
for (i in 1:numb.numb){
  proportion[i] = sum(data[data$pv_val_norm>10,]$numb==i,na.rm=TRUE)/
    sum(cluster$num_genomes_gene_cluster_has_hits==i,na.rm=TRUE)
}

# the plot for proportion of genes with p-values < 1e-10 (Supplementary)
ggplot(mapping=aes(x=1:22,y=proportion))+
  geom_point(size=3)+
  xlab("Number of MAGs") +
  ylab("Proportion of Genes \n with P < 1e-10") +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position="none")
```

<img src="pangenome-selection_files/figure-gfm/sweep-2.png" width="60%" height="60%" />

``` r
# plot for P-values vs. SNV 
data$SNV.p <- ifelse(data$SNV==0,0.5,data$SNV)
data_input <- data[data$numb %in% c(1,4,10,15,18,22),]
data_input$label <- paste("No. of MAGs = ",as.character(data_input$numb), sep="")
data_input$label <- as.factor(data_input$label)
data_input$label <- ordered(data_input$label,levels=c(
  "No. of MAGs = 1", "No. of MAGs = 4","No. of MAGs = 10",
  "No. of MAGs = 15", "No. of MAGs = 18", "No. of MAGs = 22"
)) # change labeling 

# Fig 5B
ggplot(data = data_input, mapping = aes(x=SNV.p, y=pv_val)) +
  scale_x_log10(limits=c(0.5,2000),breaks=c(1,10,100,1000))+
  scale_y_log10()+
  geom_point(alpha = 0.2,na.rm=TRUE)+
  geom_hline(yintercept=1e-10, linetype=2, color = "rosybrown", size=1)+
  xlab("Number of SNVs") +
  ylab("P-value") +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        axis.title=element_text(size=14),
        strip.text = element_text(size=12),
        legend.position="none")+
  facet_wrap(~label,ncol=3)
```

<img src="pangenome-selection_files/figure-gfm/sweep-3.png" width="60%" height="60%" />

``` r
# plot for all possible number of MAGs (Supplementary)
data$label <- paste("No. of MAGs = ",as.character(data$numb), sep="")
data$label <- as.factor(data$label)
data$label <- ordered(data$label,levels=paste("No. of MAGs = ",
                                                               as.character(c(1:22)),
                                                               sep=""))
ggplot(data = data, mapping = aes(x=SNV.p, y=pv_val)) +
  scale_x_log10(limits=c(0.5,2000),breaks=c(1,10,100,1000))+
  scale_y_log10()+
  geom_point(alpha = 0.2,na.rm=TRUE)+
  geom_hline(yintercept=1e-10, linetype=2, color = "rosybrown", size=1)+
  xlab("Number of SNVs") +
  ylab("P-value") +
  theme_bw()+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=14),
        strip.text = element_text(size=12),
        legend.position="none")+
  facet_wrap(~label,ncol=4)
```

<img src="pangenome-selection_files/figure-gfm/sweep-4.png" width="60%" height="60%" />

``` r
# then plot only for #SNV = 0
data.0 <- data[data$SNV==0,]

# Fig 5C
ggplot(data = data.0, mapping = aes(x=numb, y=pv_val)) +
  geom_boxplot(alpha = 0.3,na.rm = TRUE,outlier.size = 3,aes(colour=numb)) +
  scale_y_log10()+
  #geom_boxplot(alpha = 0,na.rm = TRUE) +
  #geom_jitter(alpha = 0.3, aes(colour = numb),na.rm=TRUE,size=0.5)+
  geom_hline(yintercept=1e-10, linetype=2, color = "rosybrown", size=1)+
  xlab("Number of MAGs") +
  ylab("P-value") +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position="none")
```

<img src="pangenome-selection_files/figure-gfm/sweep-5.png" width="60%" height="60%" />

Finally, we also studied whether there is a relationship between the
number of MAGs a cluster of some is found in and the gene’s contig
coverage. We obtained coverage data by running `get_contig_coverage.py`
which utilized the contig databases and, conveniently, an anvi’o script
`anvi-export-splits-and-coverages`. As a result, we have
`All_Sulfurovum_self_to_self_coverage.txt` which only includes
self-to-self
coverage.

``` r
coverage<-read.table('aux_files/All_Sulfurovum_self_to_self_coverage.txt',sep='\t',header=TRUE)
coverage<-dplyr::arrange(coverage,unique_id)
coverage$number <- as.factor(cluster$num_genomes_gene_cluster_has_hits)
ggplot(data=coverage,mapping=aes(x=number,y=coverage)) +
  geom_boxplot(alpha = 0.3,na.rm = TRUE,outlier.size = 3,aes(colour=number)) +
  scale_y_log10()+
  xlab("Number of MAGs") +
  ylab("Average Coverage\nof Contig") +
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.position="none")
```

<img src="pangenome-selection_files/figure-gfm/coverage_problem-1.png" width="60%" height="60%" />

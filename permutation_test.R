### this is R script to calculate the p-values of categories proportion over gene frequency

# data and manipulation
COG <- read.table('txt_files/COG_category.txt',sep='\t',header=TRUE)
Categories <- read.table('txt_files/COG_explanation.txt',sep='\t',header=FALSE)
colnames(COG) <- c('Category',as.character(c(1:22)))
COG_prop <-  data.frame(prop.table(as.matrix(COG[,-1]), margin = 2))
COG_prop$Group <- COG$Category
COG_melt <- melt(COG, id=c("Category"))

# now being simulation for p-value using permuation test
x<-16
N<-10000
n.categories <- 24
n.genomes <- 22
per.mag.sum <- colSums(COG[,-1])
per.category.sum <- rowSums(COG[,-1])
sim.cat.pool <- c()
sim.mag.pool <- c()
for(i in 1:24){
  sim.cat.pool <- c(sim.cat.pool,rep(i,per.category.sum[i]))
}
for(i in 1:22){
  sim.mag.pool <- c(sim.mag.pool,rep(i,per.mag.sum[i]))
}
COG.matrix<- as.matrix(COG[,-1])
total.genes<-sum(COG.matrix)
sim.matrix<- matrix(0,nrow=24,ncol=22)
count.less.matrix <- matrix(0,nrow=24,ncol=22)

# permutation test
set.seed(8907)
for(i in 1:N){
  sim.matrix<- matrix(0,nrow=24,ncol=22)
  sim.mag <- sample(sim.mag.pool)
  for(i in 1:total.genes){
    sim.matrix[sim.cat.pool[i],sim.mag[i]] <- 
      sim.matrix[sim.cat.pool[i],sim.mag[i]]+1
  }
  less.matrix <- sim.matrix < COG.matrix
  count.less.matrix <- count.less.matrix + less.matrix
}
p.matrix<-as.matrix(count.less.matrix/N)
p.df <- as.data.frame(count.less.matrix/N)

save(p.matrix,p.df,file="permutation_test.RData")


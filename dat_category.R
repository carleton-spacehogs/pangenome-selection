# this calculates the pN_pS data by categories

dat<-read.table('aux_files/Sulfurovum_pn_ps_summary.txt',sep='\t',header=TRUE,na.strings="Inf")
dat_counts<-na.omit(dat)
dat_counts<-dat_counts[is.finite(dat_counts$pN_pS),]
dat_counts$num_genomes <- factor(dat_counts$num_genomes)
dat_counts$type <- factor(ifelse(dat_counts$num_genomes == 22, "N=22", ifelse(dat_counts$num_genomes %in% 1:19,"0<N<20",NA))) #arbitrady core and accessory naming

separate_cat <- function(row){ # create at function
  split<-unlist(strsplit(as.character(row$category), "[|]"))
  if(length(split) != 0){
    matrix_id <- matrix(NA,ncol=8,nrow=length(split))
    for(i in 1:length(split)){
      row$category <- split[i]
      entry <- unname(row)
      matrix_id[i,] <- as.matrix(entry)
    }
    return(matrix_id)
  }
  else{
    return(as.matrix(unname(row)))
  }
}

#create the dataframe
dat_category <- matrix(ncol=8)
for(i in 1:nrow(dat_counts)){
  dat_category<-rbind(dat_category,separate_cat(dat_counts[i,]))
}
dat_category <- na.omit(as.data.frame(dat_category[dat_category$category!="",]))
colnames(dat_category)<-colnames(dat_counts[,1:8])
dat_category$pN_pS<-as.numeric(as.character(dat_category$pN_pS))

#now save the data frame
save(dat_category,file="dat_category.RData")

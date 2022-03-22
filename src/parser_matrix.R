#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
inputfile <- args[1]
chrlist <- args[2]
outputname <- args[3]

library(stringr)
data <- read.table(inputfile,header = F,stringsAsFactors = F)
header <- read.table(chrlist,header = F,stringsAsFactors = F)

colnames(data) <- header$V1

Nrow <- nrow(data)
for(i in 1:ncol(data)){
  name <- colnames(data)[i]
  name2 <- unlist(strsplit(name,'_'))[1]
  tmp_cluster <- read.table(paste0(name2,'_reannotation/',name,'.cluster'),header = F,stringsAsFactors = F)
  tmp_cluster_num <- unique(tmp_cluster$V5)
  have_cluster_num <- data[,i]
  have_cluster_num <- have_cluster_num[have_cluster_num!="."]
  have_cluster_num <- unlist(strsplit(as.character(have_cluster_num),','))
  diff <- setdiff(tmp_cluster_num, have_cluster_num)
  if(length(diff)>0){
    for(num in 1:length(diff)){
      data[Nrow+num,name] <- diff[num]
    }
    Nrow <- nrow(data) 
  }
}

data[is.na(data)] <- '.'

# filter pseduogene

speclist <- header$V1
for(spec in speclist){
  name2 <- unlist(strsplit(spec,'_'))[1]
  assign(paste0(spec,'.cluster'),
         read.table(paste0(name2,'_reannotation/',spec,'.cluster'),header=F,stringsAsFactors = F),
         envir = .GlobalEnv)
}

data2 <- data

remain_c <- c()

for(i in 1:nrow(data2)){
  remain <- F
  for(j in 1:ncol(data2)){
    spec <- colnames(data)[j]
    if(data[i,j]!="."){
      #
      cluster_num <- unlist(strsplit(data[i,j],','))
      assign('cluster',get(paste0(spec,'.cluster')))
      cluster_tmp <- cluster[which(cluster[,5]%in%cluster_num),]
      pseduo_gum <- length(which(str_detect(cluster_tmp[,4],'pseudogene')))
      total_num <- nrow(cluster_tmp)
      data2[i,j] <- paste0(pseduo_gum,'/',total_num)
      if(pseduo_gum/total_num<1){
        remain <- T
      }
    }
  }
  remain_c <- c(remain_c,remain)
}

data3 <- data[remain_c,]

write.table(data,paste0(outputname,'_allCollinearityMatrix.txt'),row.names = F,sep = '\t',quote = F)
write.table(data3,paste0(outputname,'_CodingGeneCollinearityMatrix.txt'),row.names = F,sep = '\t',quote = F)

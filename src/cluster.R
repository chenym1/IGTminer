#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
spec <- args[1]
total_bed_file <- args[2]
bed_file <- args[3]
distance <- args[4]
outfile <- args[5]

# cluster
system(paste0("bedtools cluster -i ",bed_file," -d ",distance," | gawk -vOFS='\t' '{print $1,$2,$3,$4,$7,$6}' > ",outfile))

data <- read.table(outfile,header = F,stringsAsFactors = F)

bed <- read.table(total_bed_file,header = F,stringsAsFactors = F)

chromosome_list <- unique(bed$V1)

# call list for merging cluster
merge_result <- c()

for(chromosome in chromosome_list){
  sub_data <- data[which(data$V1==chromosome),]
  if(nrow(sub_data)>0){
    cluster_list <- unique(sub_data$V5)
    if(length(cluster_list)>1){
      for(i in 1:(length(cluster_list)-1)){
        #
        cluster_tmp_1 <- sub_data[which(sub_data$V5==cluster_list[i]),]
        cluster_tmp_1 <- cluster_tmp_1[order(cluster_tmp_1$V2),]
        end_index_1 <- which(bed$V4==cluster_tmp_1[nrow(cluster_tmp_1),4])
        end_location_1 <- cluster_tmp_1[nrow(cluster_tmp_1),3]
        #
        cluster_tmp_2 <- sub_data[which(sub_data$V5==cluster_list[i+1]),]
        cluster_tmp_2 <- cluster_tmp_2[order(cluster_tmp_2$V2),]
        start_index_2 <- which(bed$V4==cluster_tmp_2[1,4])
        start_location_2 <- cluster_tmp_2[1,2]
        # merge ( number of gene in internal region <= 1)
        if(start_index_2-end_index_1 == 1){
          merge_result <- c(merge_result,paste0(cluster_list[i],',',cluster_list[i+1]))
        }
      } 
    }
  }
}

# merge cluster
if(length(merge_result)>0){
  to_merge_list <- unique(unlist(strsplit(merge_result,',')))
  tmp_group <- matrix(rep(0,(length(to_merge_list)*length(to_merge_list))),ncol = length(to_merge_list))
  dimnames(tmp_group) <- list(`1`=to_merge_list,`2`=to_merge_list)
  for(i in 1:length(merge_result)){
    tmp_out <- unlist(strsplit(merge_result[i],','))
    tmp_group[tmp_out[1],tmp_out[2]] <- 1
    tmp_group[tmp_out[2],tmp_out[1]] <- 1
  }
  library(MCL)
  out <- MCL::mcl(x = tmp_group, addLoops = T,max.iter = 1000,inflation = 1.5)
  out_group <- data.frame(cluster=rownames(tmp_group),group=out$Cluster)
  
  for(i in unique(out_group$group)){
    to_merge_list_tp <- as.character(out_group[which(out_group$group==i),1])
    data[which(data$V5%in%to_merge_list_tp),5] <- paste(to_merge_list_tp,collapse = '_')
  } 
}

data2 <- data
cluster_list <- unique(data2$V5)
for(i in 1:length(cluster_list)){
  data2[which(data2$V5==cluster_list[i]),5] <- i
}

write.table(data2,outfile,quote = F,row.names = F,col.names = F,sep = '\t')




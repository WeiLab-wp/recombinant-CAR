rm(list = ls())
library(pheatmap)

my_data <-read.csv("../Data/1.Normalization.fpkm.mean.csv",row.names = 1)
subset.fpkm<- subset(my_data,my_data[,1]>1|my_data[,2]>1|my_data[,3]>1|
                             my_data[,4]>1|my_data[,5]>1|my_data[,6]>1|
                             my_data[,7]>1|my_data[,8]>1|my_data[,9]>1|
                             my_data[,10]>1)

write.csv(subset.fpkm,"../Data/2_2.screen_1.fpkm.csv")

result <- subset.fpkm[,c(1,3,5,7,9)]
for(i in 1:5){
  result[,i] <- (subset.fpkm[,i*2-1]) / (subset.fpkm[,i*2])}

result <- log2(result)
write.csv(result,"../Data/4.screen_1.log_fold_change.csv")

pheatmap(result,clustering_method = "ward.D", cluster_cols = T,cluster_rows = T,treeheight_row = 50,border= F,
         color = colorRampPalette(c("blue3","royalblue","white","red" ,"red3"))(120),annotation_names_row = F,show_rownames = F)



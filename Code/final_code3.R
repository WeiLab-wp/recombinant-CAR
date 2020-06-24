#############NFAT-NKB factor nalysys#############
rm(list = ls())
library(pheatmap)
mydata <- read.csv("../sw.txt",sep = "")
rownames(mydata) <- mydata$time.h.
mydata <- mydata[,-1]

MAX1 <- max(mydata[,1:4])
MIN1 <- min(mydata[,1:4])
mydata_p1 <- (mydata[,1:4]-MIN1)/(MAX1-MIN1)

MAX2 <- max(mydata[,5:8])
MIN2 <- min(mydata[,5:8])
mydata_p2 <- (mydata[,5:8]-MIN2)/(MAX2-MIN2)
mydata_1 <- cbind(mydata_p1,mydata_p2)

x <- mydata_1

pheatmap(x,clustering_method = "ward.D", cluster_cols = F,treeheight_row = 50,border= F,
         color = colorRampPalette(c("blue","white", "red"))(120),
         gaps_col = 4, cutree_rows = 3)
###############################################
x1 <- mydata_p1
pheatmap(x1,clustering_method = "ward.D", cluster_cols = F,cluster_rows = F,treeheight_row = 50,border= F,
         color = colorRampPalette(c("lightgoldenrodyellow", "olivedrab2","green4"))(120))

x2 <- mydata_p2
pheatmap(x2,clustering_method = "ward.D", cluster_cols = F,cluster_rows = F,treeheight_row = 50,border= F,
         color = colorRampPalette(c("linen", "orangered1","red"))(120))



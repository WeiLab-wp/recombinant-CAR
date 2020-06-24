rm(list = ls())
all.fpkm<-read.csv("../Data/0.all.fpkm.csv")
rownames(all.fpkm)<- all.fpkm$gene_id
fpkm <- all.fpkm[,c(12,4,6,8,10,14,16,18,20,22,24,26,28,30,32,34,36,38,40)]

                ###################Normalization#####################
GAPDH <- fpkm[(rownames(fpkm)%in% 'ENSG00000111640'),]
for (i in 1:19) {
  fpkm[,i]=fpkm[,i]/GAPDH[,i]*1000
}

write.csv(fpkm,"../Data/1.Normalization.fpkm.csv")

fpkm.mean<- fpkm[,c(1,3,5,7,9,11,13,15,17,19)]
for(i in 1:9) {
  fpkm.mean[,i+1] <- (fpkm[,2*i]+fpkm[,2*i+1])/2}
fpkm.mean <- fpkm.mean[,c(1,4,2,3,5:10)]
fpkm.mean <- fpkm.mean+0.00000001

write.csv(fpkm.mean,"../Data/1.Normalization.fpkm.mean.csv")
               
                  #####################screening#########################
subset.fpkm<- subset(fpkm.mean,fpkm.mean[,1]>5|fpkm.mean[,2]>5|fpkm.mean[,3]>5|
                               fpkm.mean[,4]>5|fpkm.mean[,5]>5|fpkm.mean[,6]>5|
                               fpkm.mean[,7]>5|fpkm.mean[,8]>5|fpkm.mean[,9]>5|
                               fpkm.mean[,10]>5)

write.csv(subset.fpkm,"./Data/2.screen_5.fpkm.csv")

                  ###################gene ID_name translation##########################
library("biomaRt", lib.loc="~/R/win-library/3.5")
my_mar<- useMart("ENSEMBL_MART_ENSEMBL")
my_dataset <- useDataset("hsapiens_gene_ensembl", mart = my_mar)
ensembl <- rownames(subset.fpkm)
symb <-getBM(values = ensembl, attributes = c("hgnc_symbol","ensembl_gene_id"),
             filters = "ensembl_gene_id",  mart = my_dataset)
symb <- cbind(subset(subset.fpkm,rownames(subset.fpkm)%in% symb$ensembl_gene_id),symb$hgnc_symbol)
write.csv(symb,"../Data/3.with_symb.csv")

                  ###################recommanded genes##########################

recommanded <-read.csv("../Data/recommanded.txt")
recommanded.fpkm <- subset(subset.fpkm,rownames(subset.fpkm)%in% recommanded$symbl)

#sti
sti.exp <- recommanded.fpkm[,c(9,3,7,5)]
normalization.sti.exp <- sti.exp
for(i in 1:37){
  Norma <- max(sti.exp[i,])-min(sti.exp[i,])
  normalization.sti.exp [i,]<- (sti.exp[i,]-min(sti.exp[i,]))/Norma
}
write.csv(normalization.sti.exp,'recommanded.sti.exp.csv')
#unsti
unsti.exp <- recommanded.fpkm[,c(10,4,8,6)]
normalization.unsti.exp <- unsti.exp
for(i in 1:37){
  Norma <- max(unsti.exp[i,])-min(unsti.exp[i,])
  normalization.unsti.exp [i,]<- (unsti.exp[i,]-min(unsti.exp[i,]))/Norma
}
write.csv(normalization.unsti.exp,'recommanded.unsti.exp.csv')

pheatmap(normalization.sti.exp,clustering_method = "ward.D", cluster_cols = F,cluster_rows = F,treeheight_row = 50,border= F,
         color = colorRampPalette(c("blue3","royalblue","white","red" ,"red3"))(120),annotation_names_row = F,show_rownames=T)

pheatmap(normalization.unsti.exp,clustering_method = "ward.D", cluster_cols = F,cluster_rows = T,treeheight_row = 50,border= F,
         color = colorRampPalette(c("blue3","royalblue","white","red" ,"red3"))(120),annotation_names_row = F,show_rownames=T)
                  
               ####################################

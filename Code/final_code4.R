rm(list = ls())

my_data <-read.csv("../Data/1.Normalization.fpkm.mean.csv",row.names = 1)

log2_fold_change <- function(a,b){
  s1 <- a[,c(b,b+1)]
  s2 <- subset(s1,s1[,1]>1|s1[,2]>1) ## screen:fpkm>1
  s2[,3] <- log2(s2[,1]/s2[,2]) ## log_fold_change
  s2 <- subset(s2,s2[,3]>1)  ## upregulated
  names(s2)[3] <- 'log2_foldchange'
  c <- names(s2)[1]
  write.csv(s2,file = paste(paste('../Data/up_',c,sep = ''),'.csv',sep = ''))
}
for (i in c(3,5,7,9)) {
  log2_fold_change(my_data,i)
}

### DEG by GENEWIZ
L.CD28 <-read.csv("../Data/up_S.102.sti.1_FPKM.csv")
L.BB <-read.csv("../Data/up_S.73.sti.1_FPKM.csv")
L.Third <-read.csv("../Data/up_S.103.sti.1_FPKM.csv")
L.NM <-read.csv("../Data/up_S.124.sti.1_FPKM.csv")

# ID conversion
library(org.Hs.eg.db)
idfound <- mappedRkeys(org.Hs.egREFSEQ)
ENSEMBL <- toTable(org.Hs.egENSEMBL)
SYMBOL <- toTable(org.Hs.egSYMBOL)

m <- match(ENSEMBL$gene_id, SYMBOL$gene_id)
ENSEMBL$symbol <- SYMBOL$symbol[m]


L.CD28.symbol <- subset(L.CD28, L.CD28$X %in% ENSEMBL$ensembl_id)
L.BB.symbol <- subset(L.BB, L.BB$X %in% ENSEMBL$ensembl_id)
L.Third.symbol <- subset(L.Third, L.Third$X %in% ENSEMBL$ensembl_id)
L.NM.symbol <- subset(L.NM, L.NM$X %in% ENSEMBL$ensembl_id)

m.CD28 <- match(L.CD28.symbol$X, ENSEMBL$ensembl_id)
L.CD28.symbol$symbol <- ENSEMBL$symbol[m.CD28]

m.BB <- match(L.BB.symbol$X, ENSEMBL$ensembl_id)
L.BB.symbol$symbol <- ENSEMBL$symbol[m.BB]

m.Third <- match(L.Third.symbol$X, ENSEMBL$ensembl_id)
L.Third.symbol$symbol <- ENSEMBL$symbol[m.Third]

m.NM <- match(L.NM.symbol$X, ENSEMBL$ensembl_id)
L.NM.symbol$symbol <- ENSEMBL$symbol[m.NM]


vennlist.4CAR <- list( CD28 = na.omit(unique(L.CD28.symbol$symbol)),
                       BB = na.omit(unique(L.BB.symbol$symbol)),
                       Third =na.omit(unique(L.Third.symbol$symbol)),
                       NM =na.omit(unique(L.NM.symbol$symbol)))
venn.diagram(vennlist.4CAR,
             fill=c("#EE9A00","#CD3278","grey","SpringGreen4"), alpha=c(0.5,0.5,0.5,0.5),cat.cex = 2,cat.pos = 0, sub.cex=4,
             col = "transparent",
             cex=2,  fontfamily = "serif",
             filename = "../Data/L.BB_28_Third_NM.png")


area.4CAR <- calculate.overlap(vennlist.4CAR)
only.28 <- data.frame(gene = area.4CAR[["a9"]])
only.BB <- data.frame(gene = area.4CAR[["a14"]])
only.Third <- data.frame(gene = area.4CAR[["a1"]])
only.NM <- data.frame(gene = area.4CAR[["a3"]])

colnames(L.CD28.symbol) = colnames(L.BB.symbol) = colnames(L.Third.symbol) = colnames(L.NM.symbol) <- c('id', 'unsti', 'sti', 'log2FC','symbol')
All.car <- rbind(L.CD28.symbol, L.BB.symbol, L.Third.symbol, L.NM.symbol)

symbols =data.frame(symbols = unique(All.car$symbol))
data_car <- matrix(0,nrow(symbols),4)
colnames(data_car) <- c('CD28', 'BB', 'Third','NM' )

for (i in 1:nrow(symbols)){
  sybi =  symbols[i,1] 
  # CD28 
  if (sybi %in% L.CD28.symbol$symbol == TRUE){
    datai = subset(L.CD28.symbol$log2FC, L.CD28.symbol$symbol == sybi)
    data_car[i,1] = datai
  }
  else {
    data_car[i,1] = 0 
  }
  
  # BB
  if (sybi %in% L.BB.symbol$symbol == TRUE){
    datai = subset(L.BB.symbol$log2FC, L.BB.symbol$symbol == sybi)
    data_car[i,2] = datai
  }
  else {
    data_car[i,2] = 0 
  }  
  
  # Third
  if (sybi %in% L.Third.symbol$symbol == TRUE){
    datai = subset(L.Third.symbol$log2FC, L.Third.symbol$symbol == sybi)
    data_car[i,3] = datai
  }
  else {
    data_car[i,3] = 0 
  }    
  
  # NM
  if (sybi %in% L.NM.symbol$symbol == TRUE){
    datai = subset(L.NM.symbol$log2FC, L.NM.symbol$symbol == sybi)
    data_car[i,4] = datai
  }
  else {
    data_car[i,4] = 0 
  }   
}


data_car_symbols <- cbind(symbols,data_car)
data_S <- matrix(0,nrow(data_car_symbols),1)

for (i in 1:nrow(data_car_symbols)){
  data_S[i] <- which.max(data_car_symbols[i,2:5])}

max.CD28 <- subset(data_car_symbols, data_S[,1] == 1)
max.CD28 <- max.CD28[order(-max.CD28$CD28),]

max.BB <- subset(data_car_symbols, data_S[,1] == 2)
max.BB <- max.BB[order(-max.BB$BB),]

max.3rd <- subset(data_car_symbols, data_S[,1] == 3)
max.3rd <- max.3rd[order(-max.3rd$Third),]

max.NM <- subset(data_car_symbols, data_S[,1] == 4)
max.NM <- max.NM[order(-max.NM$NM),]

UP_max_ordered <- rbind(max.CD28, max.BB, max.3rd, max.NM)
# UP_max_ordered <- rbind(max.3rd, max.NM)


library(pheatmap)
library(RColorBrewer)

groups =data.frame(`CAR Type` = c("CD19-CD28", "CD19-41BB", "CD19-3rd", "CD19-NM"))

# Color <- colorRampPalette(c("#FFF0F5","#DB7093","#CD5C5C"))(50)
Color <- colorRampPalette(c("#F8F8FF","#DB7093","#CD5C5C"))(50)

pdf(file =paste("../Data/L.heatmap.4CAR.overlap2.pdf",sep=""),width = 12,height = 22)
pheatmap(UP_max_ordered[,2:5],
         display_numbers = F,cluster_cols = F,border=FALSE,
         treeheight_row = 30,
         cellwidth = 40, cellheight =1,fontsize=9,
         cluster_rows = F,
         color = Color,
         # scale = "row",
         scale = "none",
         legend = TRUE,
         # gaps_row = c(1061,1446),
         # cutree_row = 3,
         fontsize_row = 0.5, fontsize_col = 10,
         labels_row = UP_max_ordered$gene)
dev.off()
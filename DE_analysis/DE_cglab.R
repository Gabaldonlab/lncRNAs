library(DESeq2)
library("RColorBrewer")
colors <-  brewer.pal(8, "Set1")
library(reshape2)
library(ggplot2)
library(ggrepel)
library(tibble)   ### for inserting coulmn in dataframe
library(EDASeq)
library(VennDiagram)
library(openxlsx)
library(RUVSeq)
library(GO.db)
library(clusterProfiler)




##### setting wd to DE_analysis folder
local<-"~/path/to/DE_analysis/"
setwd(local)
################### C. GLABRATA DATA ####################


cglab_reordered<-as.matrix(read.table("./cglab_data.txt",check.names = F))


### Merge technical replicates
yeast<-cglab_reordered


yeast<-as.data.frame(yeast)
yeast<-add_column(yeast, "24+reseq" = yeast$`24`+yeast$`24reseq`, .after = "23")
yeast$`24reseq`<-NULL
yeast$`24`<-NULL


################### for PCA in all_spp_comparision

CGLAB_data<-yeast

########## DESeq2 analysis yeast #######

norm_counts <- t(t(CGLAB_data)/colSums(CGLAB_data))*100


colData_yeast<-data.frame(time=factor(x=c(rep("0",2), rep("3",3), rep("12", 3), rep("24",3), rep("24c",2)),levels=c("0","3","12","24", "24c")))
pre_dds_cglab <- DESeqDataSetFromMatrix(yeast, colData=colData_yeast, ~time)
dds_cglab <- DESeq(pre_dds_cglab)


##### all compared with 0 time point########
res_cglab_0_3 <- results(dds_cglab, cooksCutoff=FALSE, contrast = c("time","3","0"))
write.table(res_cglab_0_3, file="./cglab/res_cglab_0_vs_3.txt", sep = "\t", quote = FALSE)

res_cglab_0_12 <- results(dds_cglab, cooksCutoff=FALSE, contrast = c("time","12","0"))
write.table(res_cglab_0_12, file="./cglab/res_cglab_0_vs_12.txt", sep = "\t", quote = FALSE)

res_cglab_0_24 <- results(dds_cglab, cooksCutoff=FALSE, contrast = c("time","24","0"))
write.table(res_cglab_0_24, file="./cglab/res_cglab_0_vs_24.txt", sep = "\t", quote = FALSE)


res_cglab_0_24c <- results(dds_cglab, cooksCutoff=FALSE, contrast = c("time","24c","0"))
write.table(res_cglab_0_24c, file="./cglab/res_cglab_0_vs_24c.txt", sep = "\t", quote = FALSE)


########################################


all_data_cglab<-cbind.data.frame(res_cglab_0_3$log2FoldChange,res_cglab_0_3$padj,
                                 res_cglab_0_12$log2FoldChange,res_cglab_0_12$padj,
                                 res_cglab_0_24$log2FoldChange,res_cglab_0_24$padj)


rownames(all_data_cglab)<-rownames(res_cglab_0_12)
colnames(all_data_cglab)<-c("3", "padj3", "12", "padj12", "24", "padj24") #, "24c", "padj24c")


all_data_cglab<-all_data_cglab[complete.cases(all_data_cglab),]

lncRNA_data<-as.data.frame(subset(all_data_cglab, grepl("MSTRG", row.names(all_data_cglab))))



lncRNA_data$`3`[lncRNA_data$padj3 > 0.01] <- 0   #### if padj > 0.01 then LFC is 0 (if not significant then no LFC)
lncRNA_data$`12`[lncRNA_data$padj12 > 0.01] <- 0
lncRNA_data$`24`[lncRNA_data$padj24 > 0.01] <- 0


LFC_cglab<-cbind.data.frame("0" = 0, lncRNA_data$`3`, lncRNA_data$`12`, lncRNA_data$`24`)
colnames(LFC_cglab)<-c("0","3","12", "24")
rownames(LFC_cglab)<-rownames(lncRNA_data)


####################plotting interpolated values############################################

LFC<-LFC_cglab
dim(LFC)
################## FILTERING !!!

LFC<-LFC[apply(LFC, 1, function(x) !all(x==0)),]
dim(LFC)
############ Plotting of interpolated LFC #################3
time<-c(0,3,12,24)
time_interpol<-approx(time)$y


LFC_no_NA<-LFC[complete.cases(LFC),]
LFC_no_NA_t<-t(LFC_no_NA)

### Interpolation
interpol_LFC<-apply(LFC_no_NA_t, 2, approx)

### rbind only y values from interpol_LFC

combined_interplolated<-do.call(rbind, lapply(interpol_LFC, '[[','y'))

### parsing and plotting
colnames(combined_interplolated)<-time_interpol

LFC_plotting_inter<-melt(combined_interplolated,value.name = "value", varnames=c('Var1', 'Var2'))


ggplot(data=LFC_plotting_inter, aes(x=LFC_plotting_inter$Var2, y=LFC_plotting_inter$value, group=LFC_plotting_inter$Var1, colour=LFC_plotting_inter$value)) +
  geom_line() + xlab("Time (h)")+ylab("")+
  scale_colour_gradient2(low = "royalblue3", mid = "grey", high = "red",space = "Lab",na.value = "grey50", guide = "colourbar")+theme_bw()+labs(colour = "L2FC\ngradient")+
  theme(axis.text = element_text(size = 18),axis.title = element_text(size=18),
        legend.position = "none")+labs(colour = "L2FC\ngradient")+
  scale_y_continuous(breaks=seq(-10,10,2))+expand_limits(y=c(-10,10))+scale_x_continuous(breaks=c(0,3, 12, 24))+
  ggsave("./cglab/cglab_lncRNA_dyn.pdf", width = 6, height = 4,dpi = 600)


#### Numbers of DE lncRNAs
#### 3h
dim(LFC[LFC$`3`>1.5 & grepl("\\|u\\|", row.names(LFC)),])
dim(LFC[LFC$`3`>1.5 & grepl("\\|x\\|", row.names(LFC)),])

dim(LFC[LFC$`3`< -1.5 & grepl("\\|u\\|", row.names(LFC)),])
dim(LFC[LFC$`3`< -1.5 & grepl("\\|x\\|", row.names(LFC)),])

#### 12h
dim(LFC[LFC$`12`>1.5 & grepl("\\|u\\|", row.names(LFC)),])
dim(LFC[LFC$`12`>1.5 & grepl("\\|x\\|", row.names(LFC)),])

dim(LFC[LFC$`12`< -1.5 & grepl("\\|u\\|", row.names(LFC)),])
dim(LFC[LFC$`12`< -1.5 & grepl("\\|x\\|", row.names(LFC)),])

#### 24h
dim(LFC[LFC$`24`>1.5 & grepl("\\|u\\|", row.names(LFC)),])
dim(LFC[LFC$`24`>1.5 & grepl("\\|x\\|", row.names(LFC)),])

dim(LFC[LFC$`24`< -1.5 & grepl("\\|u\\|", row.names(LFC)),])
dim(LFC[LFC$`24`< -1.5 & grepl("\\|x\\|", row.names(LFC)),])




########## test LRT############

##### cglab


colData_yeast_lrt<-data.frame(time=factor(x=c(rep("0",2), rep("3",3), rep("12", 3), rep("24",3)),levels=c("0","3","12","24")))
yeast_lrt<-yeast[,-c(dim(yeast)[2]-1,dim(yeast)[2])] ### remove last two columns

pre_dds_glab_lrt <- DESeqDataSetFromMatrix(yeast_lrt, colData=colData_yeast_lrt, ~time)

dds_lrt_cglab <- DESeq(pre_dds_glab_lrt, test="LRT", reduced=~1)
res_lrt_cglab <- results(dds_lrt_cglab)


non_DE_cglab<-res_lrt_cglab[!is.na(res_lrt_cglab$padj) & res_lrt_cglab$padj>0.5 & res_lrt_cglab$baseMean>100,]


write.table(res_lrt_cglab, file="./cglab/res_lrt_cglab.txt", sep = "\t", quote = FALSE)


### 24 vs 24c comparison

lncRNA_0_24<-rownames(res_cglab_0_24[!is.na(res_cglab_0_24$padj) & res_cglab_0_24$padj<0.01,])
lncRNA_0_24<-lncRNA_0_24[grepl("^MSTRG", lncRNA_0_24)]

lncRNA_0_24c<-rownames(res_cglab_0_24c[!is.na(res_cglab_0_24c$padj) & res_cglab_0_24c$padj<0.01,])
lncRNA_0_24c<-lncRNA_0_24c[grepl("^MSTRG", lncRNA_0_24c)]


overlap<-as.data.frame(Reduce(intersect,list(lncRNA_0_24,lncRNA_0_24c)))

# 
# venn.diagram(list(" "=lncRNA_0_24,
#                   " "=lncRNA_0_24c),
#              fill=c("green3", "grey"), height = 5, width = 5, resolution = 600,
#              main = " ", cat.cex=0.4,cat.pos=c(30,320), cex=1.8, main.cex = 0.7, 
#              filename="./cglab/cglab_24_24c.svg",  
#              imagetype = "svg",inverted=TRUE)

pdf("./cglab/cglab_lncRNA_24_24c.pdf")
venn_plot<-draw.pairwise.venn(area1=length(lncRNA_0_24),
                              area2=length(lncRNA_0_24c),
                              cross.area = length(intersect(lncRNA_0_24,lncRNA_0_24c)),
                              fill=c("green3", "grey"),cex=3,inverted=TRUE )
dev.off()


only_infection_cglab<-setdiff(lncRNA_0_24, lncRNA_0_24c)
infection_specific_lncRNAs<-res_cglab_0_24[only_infection_cglab,]
write.table(infection_specific_lncRNAs,"./cglab/infection_lncRNAs_cglab_24_24c.txt", quote = FALSE, sep = "\t")


### Creating a new xlsx file

res_cglab_lrt<-res_lrt_cglab

wb <- createWorkbook("cglab_data")


header<-c("gene", "baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
fungal_data_files<-c("res_cglab_0_3",
                     "res_cglab_0_12",
                     "res_cglab_0_24",
                     "res_cglab_0_24c",
                     "res_cglab_lrt",
                     "infection_specific_lncRNAs")



descriptions<-c("C. glabrata at 3 hpi compared to 0 h control",
                "C. glabrata at 12 hpi compared to 0 h control",
                "C. glabrata at 24 hpi compared to 0 h control",
                "C. glabrata at 24 h control compared to 0 h control",
                "LRT test for C. glabrata samples (differential expression at any time point)",
                "Infection-specific DE lncRNAs in C. glabrata (DE exclusively during infection)")

sheets_and_descriptions<-cbind.data.frame("Data_sheets"=fungal_data_files,"Description"=descriptions)

addWorksheet(wb,"Sheets_and_descriptions")
writeData(wb, "Sheets_and_descriptions", sheets_and_descriptions,keepNA = T)




for (file in fungal_data_files){
  assign("fung",as.data.frame(eval(parse(text=paste0(file)))))
  fung <- cbind("genes"=rownames(fung), data.frame(fung, row.names=NULL))
  addWorksheet(wb, file)
  writeData(wb, file, fung,keepNA = T)
}

saveWorkbook(wb, "./Suppl_dataset_S6.xlsx", overwrite = TRUE)


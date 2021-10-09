library(openxlsx)
library(DESeq2)
library("RColorBrewer")
colors <-  brewer.pal(8, "Set1")
library(reshape2)
library(ggplot2)
library(ggrepel)
library(tibble)   ### for inserting coulmn in dataframe
library(EDASeq)
library(RUVSeq)
library(VennDiagram)
library(GO.db)
library(clusterProfiler)


##### setting wd to DE_analysis folder
local<-"~/Dropbox/lncRNA/for_github/DE_analysis/"
setwd(local)
################### C. ALBICANS DATA ####################


calb_with_e_and_D<-as.matrix(read.table("./calb_data.txt", check.names = F))


### remove e and D columns
calb <- subset(calb_with_e_and_D, select = -c(e4,e5,e6,e7,e8,e9,e1,e2,e3))


####################################################

### Merging technical replicates
yeast<-calb

yeast<-as.data.frame(yeast)
yeast<-add_column(yeast, "A0C1+reseq" = yeast$A0C1+yeast$A0C1reseq, .after = "2")
yeast$A0C1<-NULL
yeast$A0C1reseq<-NULL


yeast<-add_column(yeast, "A0C2+reseq" = yeast$A0C2+yeast$A0C2reseq, .after = "A0C1+reseq")
yeast$A0C2<-NULL
yeast$A0C2reseq<-NULL


yeast<-add_column(yeast, "3+reseq" = yeast$`3`+yeast$`3reseq`, .after = "A0C2+reseq")
yeast$`3reseq`<-NULL
yeast$`3`<-NULL


yeast<-add_column(yeast, "A24C1+reseq" = yeast$A24C1+yeast$A24C1reseq, .after = "18")
yeast$A24C1<-NULL
yeast$A24C1reseq<-NULL


yeast<-add_column(yeast, "A24C2+reseq" = yeast$A24C2+yeast$A24C2reseq, .after = "A24C1+reseq")
yeast$A24C2<-NULL
yeast$A24C2reseq<-NULL


########## DESeq2 analysis #######

colData_yeast<-data.frame(time=factor(x=c(rep("0",4), rep("1.5",3), rep("3",3), rep("3c",2), rep("12", 3), rep("24",4), rep("24c",4)), levels=c("0","1.5","3","3c","12","24", "24c")))


pre_dds_calb <- DESeqDataSetFromMatrix(yeast, colData=colData_yeast, ~time)
dds_calb <- DESeq(pre_dds_calb)


##### all compared with 0 time point########
res_calb_0_1.5 <- results(dds_calb, cooksCutoff=FALSE, contrast = c("time","1.5","0"))
write.table(res_calb_0_1.5, file="./calb/res_calb_0_vs_1.5.txt", sep = "\t", quote = FALSE)


res_calb_0_3 <- results(dds_calb, cooksCutoff=FALSE, contrast = c("time","3","0"))
write.table(res_calb_0_3, file="./calb/res_calb_0_vs_3.txt", sep = "\t", quote = FALSE)


res_calb_0_3c <- results(dds_calb, cooksCutoff=FALSE, contrast = c("time","3c","0"))
write.table(res_calb_0_3c, file="./calb/res_calb_0_vs_3c.txt", sep = "\t", quote = FALSE)


res_calb_0_12 <- results(dds_calb, cooksCutoff=FALSE, contrast = c("time","12","0"))
write.table(res_calb_0_12, file="./calb/res_calb_0_vs_12.txt", sep = "\t", quote = FALSE)


res_calb_0_24 <- results(dds_calb, cooksCutoff=FALSE, contrast = c("time","24","0"))
write.table(res_calb_0_24, file="./calb/res_calb_0_vs_24.txt", sep = "\t", quote = FALSE)


res_calb_0_24c <- results(dds_calb, cooksCutoff=FALSE, contrast = c("time","24c","0"))
write.table(res_calb_0_24c, file="./calb/res_calb_0_vs_24c.txt", sep = "\t", quote = FALSE)



# ############## 4 time-points ##########################
# 
# 
# all_data<-cbind.data.frame(res_calb_0_1.5$log2FoldChange,res_calb_0_1.5$padj,
#                             res_calb_0_3$log2FoldChange,res_calb_0_3$padj,
#                             res_calb_0_12$log2FoldChange,res_calb_0_12$padj,
#                             res_calb_0_24$log2FoldChange,res_calb_0_24$padj)
# 
# 
# rownames(all_data)<-rownames(res_calb_0_12)
# colnames(all_data)<-c("1.5","padj1.5","3", "padj3", "12", "padj12", "24", "padj24") 
# 
# lncRNA_data<-as.data.frame(subset(all_data, grepl("MSTRG", row.names(all_data))))
# 
# lncRNA_data<-lncRNA_data[complete.cases(lncRNA_data),]
# 
# lncRNA_data$`1.5`[lncRNA_data$padj1.5 > 0.01] <- 0
# lncRNA_data$`3`[lncRNA_data$padj3 > 0.01] <- 0
# lncRNA_data$`12`[lncRNA_data$padj12 > 0.01] <- 0
# lncRNA_data$`24`[lncRNA_data$padj24 > 0.01] <- 0
# 
# 
# LFC<-cbind.data.frame("0" = 0,lncRNA_data$`1.5`, lncRNA_data$`3`, lncRNA_data$`12`, lncRNA_data$`24`)
# colnames(LFC)<-c("0","1.5","3","12", "24")
# rownames(LFC)<-rownames(lncRNA_data)
# 
# ################## FILTERING !!!
# LFC<-LFC[apply(LFC, 1, function(x) !all(x==0)),]    ### discard gene with no changes in expression
# 
# dim(LFC)
# ############ Plotting of interpolated LFC #################3
# time<-c(0,1.5,3,12,24)
# time_interpol<-approx(time)$y
# 
# 
# LFC_no_NA<-LFC[complete.cases(LFC),]
# LFC_no_NA_t<-t(LFC_no_NA)
# 
# ### Interpolation
# interpol_LFC<-apply(LFC_no_NA_t, 2, approx)
# 
# ### rbind only y values from interpol_LFC
# 
# combined_interplolated<-do.call(rbind, lapply(interpol_LFC, '[[','y'))
# 
# ### parsing and plotting
# colnames(combined_interplolated)<-time_interpol
# 
# LFC_plotting_inter<-melt(combined_interplolated,value.name = "value", varnames=c('Var1', 'Var2'))
# 
# png('./calb/calb_lncRNA_dyn_4time_points.png', units="in", width=10, heigh=7, res=400)
# ggplot(data=LFC_plotting_inter, aes(x=LFC_plotting_inter$Var2, y=LFC_plotting_inter$value, group=LFC_plotting_inter$Var1, colour=LFC_plotting_inter$value)) +
#   geom_line() + xlab("Time (h)")+ylab("Log2 fol change") +  ggtitle("Dynamics of lncRNAs in C. albicans (|L2FC|>0)")+
#   scale_colour_gradient2(low = "blue", mid = "grey", high = "red",space = "Lab",na.value = "grey50", guide = "colourbar")+theme_bw()+labs(colour = "L2FC\ngradient")+
#   theme(axis.text = element_text(size = 18),axis.title = element_text(size=18),
#         legend.text = element_text(size=14),legend.title = element_text(size=14))+labs(colour = "L2FC\ngradient")+
#   scale_y_continuous(breaks=seq(-10,10,2))+expand_limits(y=c(-10,10))+scale_x_continuous(breaks=c(0,1.5 ,3, 12, 24)) ### this is for changing ticks
# dev.off()




all_data<-cbind.data.frame(res_calb_0_3$log2FoldChange,res_calb_0_3$padj,
                           res_calb_0_12$log2FoldChange,res_calb_0_12$padj,
                           res_calb_0_24$log2FoldChange,res_calb_0_24$padj)


rownames(all_data)<-rownames(res_calb_0_12)
colnames(all_data)<-c("3", "padj3", "12", "padj12", "24", "padj24") #, "24c", "padj24c")

lncRNA_data<-as.data.frame(subset(all_data, grepl("MSTRG", row.names(all_data))))
lncRNA_data<-lncRNA_data[complete.cases(lncRNA_data),]


lncRNA_data$`3`[lncRNA_data$padj3 > 0.01] <- 0
lncRNA_data$`12`[lncRNA_data$padj12 > 0.01] <- 0
lncRNA_data$`24`[lncRNA_data$padj24 > 0.01] <- 0


LFC<-cbind.data.frame("0" = 0, lncRNA_data$`3`, lncRNA_data$`12`, lncRNA_data$`24`)
colnames(LFC)<-c("0","3","12", "24")
rownames(LFC)<-rownames(lncRNA_data)


################## FILTERING !!!
LFC<-LFC[apply(LFC, 1, function(x) !all(x==0)),]    ### discard gene with no changes in expression
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
  geom_line() + xlab("Time (h)")+ylab("Log2(fold change)") +  
  scale_colour_gradient2(low = "royalblue3", mid = "grey", high = "red",space = "Lab",na.value = "grey50", guide = "colourbar")+theme_bw()+labs(colour = "L2FC\ngradient")+
  theme(axis.text = element_text(size = 18),axis.title = element_text(size=18),
        legend.position = "none")+
  labs(colour = "L2FC\ngradient")+
  scale_y_continuous(breaks=seq(-10,10,2))+expand_limits(y=c(-10,10))+scale_x_continuous(breaks=c(0,3, 12, 24))+
  ggsave("./calb/calb_lncRNA_dyn.pdf", width = 6, height = 4,dpi = 600)



#### Number of DE lncRNAs
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

colData_yeast_lrt<-data.frame(time=factor(x=c(rep("0",4), rep("1.5",3), rep("3",3), rep("12", 3), rep("24",4)), levels=c("0","1.5","3","12","24")))

yeast_lrt<-yeast[,-c(11,12,20,21,22,23)] ### remove 3h controls and 4 last columns (24 h controls)

pre_dds_calb_lrt <- DESeqDataSetFromMatrix(yeast_lrt, colData=colData_yeast_lrt, ~time)

dds_lrt_calb <- DESeq(pre_dds_calb_lrt, test="LRT", reduced=~1)
res_lrt_calb <- results(dds_lrt_calb)

non_DE_calb<-res_lrt_calb[!is.na(res_lrt_calb$padj) & res_lrt_calb$padj>0.5 & res_lrt_calb$baseMean>100,]

write.table(res_lrt_calb, file="./calb/res_lrt_calb.txt", sep = "\t", quote = FALSE)


##########################################

### 24 vs 24c comparison
lncRNA_0_24<-rownames(res_calb_0_24[!is.na(res_calb_0_24$padj) & res_calb_0_24$padj<0.01,])
lncRNA_0_24<-lncRNA_0_24[grepl("^MSTRG", lncRNA_0_24)]


lncRNA_0_24c<-rownames(res_calb_0_24c[!is.na(res_calb_0_24c$padj) & res_calb_0_24c$padj<0.01,])
lncRNA_0_24c<-lncRNA_0_24c[grepl("^MSTRG", lncRNA_0_24c)]


overlap<-as.data.frame(Reduce(intersect,list(lncRNA_0_24,lncRNA_0_24c)))




pdf("./calb/calb_lncRNA_24_24c.pdf")
venn_plot<-draw.pairwise.venn(area1=length(lncRNA_0_24),
                              area2=length(lncRNA_0_24c),
                              cross.area = length(intersect(lncRNA_0_24,lncRNA_0_24c)),
                              fill=c("green3", "grey"),cex=3 )
dev.off()


only_infection_lncRNA<-setdiff(lncRNA_0_24, lncRNA_0_24c)
infection_specific_lncRNAs<-data.frame(res_calb_0_24[only_infection_lncRNA,])
infection_specific_lncRNAs[!is.na(infection_specific_lncRNAs$padj) & infection_specific_lncRNAs$padj<0.01 & abs(infection_specific_lncRNAs$log2FoldChange)>1.5 & infection_specific_lncRNAs$baseMean>100,]


write.table(infection_specific_lncRNAs,"./calb/infection_lncRNAs_calb_24_24c.txt", quote = FALSE, sep = "\t")


# ### ECE1 CALB
# 
# calb_with_e_and_D<-data.frame(calb_with_e_and_D)
# calb_ece1<-cbind(yeast[,1:4],"e4"=calb_with_e_and_D$e4,"e5"=calb_with_e_and_D$e5,"e6"=calb_with_e_and_D$e6,
#                  yeast[,5:19],"e7"=calb_with_e_and_D$e7,"e8"=calb_with_e_and_D$e8,"e9"=calb_with_e_and_D$e9,
#                  yeast[,20:23],"e1"=calb_with_e_and_D$e1,"e2"=calb_with_e_and_D$e2,"e3"=calb_with_e_and_D$e3)
# colData_yeast<-data.frame(time=factor(x=c(rep("0",4),rep("0ece1",3), rep("1.5",3), rep("3",3),rep("3c",2), rep("12", 3), rep("24",4),rep("24ece1",3) ,rep("24c",4),rep("24ece1_c",3)), 
#                                       levels=c("0","0ece1","1.5","3","3c","12","24","24ece1","24c","24ece1_c")))
# 
# 
# #### Removing batch effect
# 
# calb_D_set <- newSeqExpressionSet(as.matrix(calb_ece1),
#                                    phenoData = cbind(colData_yeast,row.names=colnames(calb_ece1)))
# 
# calb_D_removed_batch_effect <- RUVg(calb_D_set, row.names(non_DE_calb), k=3)
# 
# 
# ### DESeq2
# pre_dds_calb_D_no_batch <- DESeqDataSetFromMatrix(counts(calb_D_removed_batch_effect), colData=pData(calb_D_removed_batch_effect), ~W_1+time)
# dds_calb_D_no_batch <- DESeq(pre_dds_calb_D_no_batch)
# 
# 
# cabl_ece1_norm<-vst(dds_calb_D_no_batch)
# plot_calb_ece1<-plotPCA(cabl_ece1_norm, intgroup="time", ntop=5000)
# plot_calb_ece1+geom_text_repel(aes(label=colnames(cabl_ece1_norm)), show.legend = FALSE)+
#   scale_colour_discrete(name  ="Time points")+ theme_bw()+labs(title = "PCA plot of calb samples")
# 
# 
# res_calb_0_vs_0e <- results(dds_calb_D_no_batch, cooksCutoff=FALSE, contrast = c("time","0ece1","0"))
# write.table(res_calb_0_vs_0e, file="./calb/res_calb_0_vs_0ece1_no_batch.txt", sep = "\t", quote = FALSE)
# 
# 
# res_calb_0e_vs_24e <- results(dds_calb_D_no_batch, cooksCutoff=FALSE, contrast = c("time","24ece1","0ece1"))
# write.table(res_calb_0e_vs_24e, file="./calb/res_calb_0ece1_vs_24ece1_no_batch.txt", sep = "\t", quote = FALSE)
# 
# 
# res_calb_24e.c_vs_24e <- results(dds_calb_D_no_batch, cooksCutoff=FALSE, contrast = c("time","24ece1","24ece1_c"))
# write.table(res_calb_24e.c_vs_24e, file="./calb/res_calb_24ece1_c_vs_24ece1_no_batch.txt", sep = "\t", quote = FALSE)
# 
# 
# res_calb_24_vs_24e <- results(dds_calb_D_no_batch, cooksCutoff=FALSE, contrast = c("time","24ece1","24"))
# write.table(res_calb_24_vs_24e, file="./calb/res_calb_24_vs_24ece1_no_batch.txt", sep = "\t", quote = FALSE)


### Creating a new xlsx file


res_calb_lrt<-res_lrt_calb


wb <- createWorkbook("calb_data")
# 
# descriptions<-c("C. albicans at 1.5 hpi compared to 0 h control",
#                 "C. albicans at 3 hpi compared to 0 h control",
#                 "C. albicans at 3 h control compared to 0 h control",
#                 "C. albicans at 12 hpi compared to 0 h control",
#                 "C. albicans at 24 hpi compared to 0 h control",
#                 "C. albicans at 24 h control compared to 0 h control",
#                 "LRT test for C. albicans samples (differential expression at any time point)",
#                 "Infection-specific DE lncRNAs in C. albicans (DE exclusively during infection)",
#                 "C. albicans ece1 mutant at 0 h compared to 0 h control",
#                 "C. albicans ece1 mutant at 24 hpi compared to 0 h of C. albicans ece1 mutant",
#                 "C. albicans ece1 mutant at 24 hpi compared to 24 h control of C. albicans ece1 mutant",
#                 "C. albicans ece1 mutant at 24 hpi compared to C. albicans at 24 hpi"
# )
# 
# fungal_data_files<-c("res_calb_0_1.5",
#                      "res_calb_0_3",
#                      "res_calb_0_3c",
#                      "res_calb_0_12",
#                      "res_calb_0_24",
#                      "res_calb_0_24c",
#                      "res_calb_lrt",
#                      "infection_specific_lncRNAs",
#                      "res_calb_0_vs_0e",
#                      "res_calb_0e_vs_24e",
#                      "res_calb_24e.c_vs_24e",
#                      "res_calb_24_vs_24e")




descriptions<-c("C. albicans at 1.5 hpi compared to 0 h control",
                "C. albicans at 3 hpi compared to 0 h control",
                "C. albicans at 3 h control compared to 0 h control",
                "C. albicans at 12 hpi compared to 0 h control",
                "C. albicans at 24 hpi compared to 0 h control",
                "C. albicans at 24 h control compared to 0 h control",
                "LRT test for C. albicans samples (differential expression at any time point)",
                "Infection-specific DE lncRNAs in C. albicans (DE exclusively during infection)")



fungal_data_files<-c("res_calb_0_1.5",
               "res_calb_0_3",
               "res_calb_0_3c",
               "res_calb_0_12",
               "res_calb_0_24",
               "res_calb_0_24c",
               "res_calb_lrt",
               "infection_specific_lncRNAs")

sheets_and_descriptions<-cbind.data.frame("Data_sheets"=fungal_data_files,"Description"=descriptions)


addWorksheet(wb,"Sheets_and_descriptions")
writeData(wb, "Sheets_and_descriptions", sheets_and_descriptions,keepNA = T)


for (file in fungal_data_files){
  assign("fung",as.data.frame(eval(parse(text=paste0(file)))))
  fung <- cbind("genes"=rownames(fung), data.frame(fung, row.names=NULL))
  addWorksheet(wb, file)
  writeData(wb, file, fung,keepNA = T)
}

saveWorkbook(wb, "./Suppl_dataset_S3.xlsx", overwrite = TRUE)


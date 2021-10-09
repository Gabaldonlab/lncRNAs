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

################### C. TROPICALIS DATA ####################
ctrop_reordered<-as.matrix(read.table("./ctrop_data.txt",check.names = F))


#### discard putative tRNAs
to_discard <- read.table("../ploting_and_network_analysis/ctrop/ctrop_to_discard.txt")
to_discard <- as.character(to_discard$V1)
ctrop_reordered <- ctrop_reordered[!rownames(ctrop_reordered) %in% to_discard,]

yeast<-ctrop_reordered

### Merging technical repliates


yeast<-as.data.frame(yeast)
yeast<-add_column(yeast, "45C+reseq" = yeast$`45C`+yeast$`45Creseq`, .after = "45")
yeast$`45C`<-NULL
yeast$`45Creseq`<-NULL


yeast<-as.data.frame(yeast)
yeast<-add_column(yeast, "52+reseq" = yeast$`52`+yeast$`52reseq`+yeast$`52reseq2`, .after = "51")
yeast$`52`<-NULL
yeast$`52reseq`<-NULL
yeast$`52reseq2`<-NULL


yeast<-as.data.frame(yeast)
yeast<-add_column(yeast, "53+reseq" = yeast$`53`+yeast$`53reseq`, .after = "52+reseq")
yeast$`53reseq`<-NULL
yeast$`53`<-NULL


ctrop_data<-yeast

########## DESeq2 analysis yeast #######

norm_counts <- t(t(ctrop_data)/colSums(ctrop_data))*100

colData_yeast<-data.frame(time=factor(x=c(rep("0",3), rep("3",3), rep("12", 3), rep("24",3), rep("24c",2)),levels=c("0","3","12","24", "24c")))
pre_dds_ctrop <- DESeqDataSetFromMatrix(yeast, colData=colData_yeast, ~time)
dds_ctrop <- DESeq(pre_dds_ctrop)




##### all compared with 0 time point########
res_ctrop_0_3 <- results(dds_ctrop, cooksCutoff=FALSE, contrast = c("time","3","0"))
write.table(res_ctrop_0_3, file="./ctrop/res_ctrop_0_vs_3.txt", sep = "\t", quote = FALSE)

res_ctrop_0_12 <- results(dds_ctrop, cooksCutoff=FALSE, contrast = c("time","12","0"))
write.table(res_ctrop_0_12, file="./ctrop/res_ctrop_0_vs_12.txt", sep = "\t", quote = FALSE)

res_ctrop_0_24 <- results(dds_ctrop, cooksCutoff=FALSE, contrast = c("time","24","0"))
write.table(res_ctrop_0_24, file="./ctrop/res_ctrop_0_vs_24.txt", sep = "\t", quote = FALSE)


res_ctrop_0_24c <- results(dds_ctrop, cooksCutoff=FALSE, contrast = c("time","24c","0"))
write.table(res_ctrop_0_24c, file="./ctrop/res_ctrop_0_vs_24c.txt", sep = "\t", quote = FALSE)


########################################

all_data_ctrop<-cbind.data.frame(res_ctrop_0_3$log2FoldChange,res_ctrop_0_3$padj,
                                 res_ctrop_0_12$log2FoldChange,res_ctrop_0_12$padj,
                                 res_ctrop_0_24$log2FoldChange,res_ctrop_0_24$padj)


rownames(all_data_ctrop)<-rownames(res_ctrop_0_12)
colnames(all_data_ctrop)<-c("3", "padj3", "12", "padj12", "24", "padj24") #, "24c", "padj24c")


all_data_ctrop<-all_data_ctrop[complete.cases(all_data_ctrop),]

lncRNA_data<-as.data.frame(subset(all_data_ctrop, grepl("MSTRG", row.names(all_data_ctrop))))


lncRNA_data$`3`[lncRNA_data$padj3 > 0.01] <- 0   #### if padj > 0.01 then LFC is 0 (if not significant then no LFC)
lncRNA_data$`12`[lncRNA_data$padj12 > 0.01] <- 0
lncRNA_data$`24`[lncRNA_data$padj24 > 0.01] <- 0


LFC_ctrop<-cbind.data.frame("0" = 0, lncRNA_data$`3`, lncRNA_data$`12`, lncRNA_data$`24`)
colnames(LFC_ctrop)<-c("0","3","12", "24")
rownames(LFC_ctrop)<-rownames(lncRNA_data)
  
  
  ####################plotting interpolated values############################################
  
  LFC<-LFC_ctrop
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
    ggsave("./ctrop/ctrop_lncRNA_dyn.pdf", width = 6, height = 4,dpi = 600)
  
  
  
  
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
  
  ##### ctrop
  
  
  colData_yeast_lrt<-data.frame(time=factor(x=c(rep("0",3), rep("3",3), rep("12", 3), rep("24",3)),levels=c("0","3","12","24")))
  yeast_lrt<-yeast[,-c(dim(yeast)[2]-1,dim(yeast)[2])] ### remove last two columns
  
  pre_dds_ctrop_lrt <- DESeqDataSetFromMatrix(yeast_lrt, colData=colData_yeast_lrt, ~time)
  
  dds_lrt_ctrop <- DESeq(pre_dds_ctrop_lrt, test="LRT", reduced=~1)
  res_lrt_ctrop <- results(dds_lrt_ctrop)
  
  
  non_DE_ctrop<-res_lrt_ctrop[!is.na(res_lrt_ctrop$padj) & res_lrt_ctrop$padj>0.5 & res_lrt_ctrop$baseMean>100,]
  
  
  #write.table(res_lrt_ctrop, file="./ctrop/res_lrt_ctrop.txt", sep = "\t", quote = FALSE)
  
  
  ### 24 vs 24c comparison
  
  lncRNA_0_24<-rownames(res_ctrop_0_24[!is.na(res_ctrop_0_24$padj) & res_ctrop_0_24$padj<0.01,])
  lncRNA_0_24<-lncRNA_0_24[grepl("^MSTRG", lncRNA_0_24)]
  
  lncRNA_0_24c<-rownames(res_ctrop_0_24c[!is.na(res_ctrop_0_24c$padj) & res_ctrop_0_24c$padj<0.01,])
  lncRNA_0_24c<-lncRNA_0_24c[grepl("^MSTRG", lncRNA_0_24c)]
  
  overlap<-as.data.frame(Reduce(intersect,list(lncRNA_0_24,lncRNA_0_24c)))
  
  
  # venn.diagram(list(" "=lncRNA_0_24,
  #                   " "=lncRNA_0_24c),
  #              fill=c("green3", "grey"), height = 5, width = 5, resolution = 600,
  #              main = " ", cat.cex=0.4, cex=1.8,cat.pos=c(30,320), main.cex = 0.7, 
  #              filename="./ctrop/ctrop_24_24c.svg",  
  #              imagetype = "svg",inverted=TRUE)
  # 
  pdf("./ctrop/ctrop_lncRNA_24_24c.pdf")
  venn_plot<-draw.pairwise.venn(area1=length(lncRNA_0_24),
                                area2=length(lncRNA_0_24c),
                                cross.area = length(intersect(lncRNA_0_24,lncRNA_0_24c)),
                                fill=c("green3", "grey"),cex=3,inverted=TRUE )
  dev.off()
  
  
  
  
  
  only_infection_ctrop<-setdiff(lncRNA_0_24, lncRNA_0_24c)
  infection_specific_lncRNAs<-res_ctrop_0_24[only_infection_ctrop,]
  write.table(infection_specific_lncRNAs,"./ctrop/infection_lncRNAs_ctrop_24_24c.txt", quote = FALSE, sep = "\t")
  
  
  
  ### Creating a new xlsx file
  res_ctrop_lrt<-res_lrt_ctrop
  wb <- createWorkbook("ctrop_data")
  
  
  fungal_data_files<-c("res_ctrop_0_3",
                       "res_ctrop_0_12",
                       "res_ctrop_0_24",
                       "res_ctrop_0_24c",
                       "res_ctrop_lrt",
                       "infection_specific_lncRNAs")
  
  
  descriptions<-c("C. tropicalis at 3 hpi compared to 0 h control",
                  "C. tropicalis at 12 hpi compared to 0 h control",
                  "C. tropicalis at 24 hpi compared to 0 h control",
                  "C. tropicalis at 24 h control compared to 0 h control",
                  "LRT test for C. tropicalis samples (differential expression at any time point)",
                  "Infection-specific DE lncRNAs in C. tropicalis (DE exclusively during infection)")
  
  sheets_and_descriptions<-cbind.data.frame("Data_sheets"=fungal_data_files,"Description"=descriptions)
  
  addWorksheet(wb,"Sheets_and_descriptions")
  writeData(wb, "Sheets_and_descriptions", sheets_and_descriptions,keepNA = T)
  
  
  
  for (file in fungal_data_files){
    assign("fung",as.data.frame(eval(parse(text=paste0(file)))))
    fung <- cbind("genes"=rownames(fung), data.frame(fung, row.names=NULL))
    addWorksheet(wb, file)
    writeData(wb, file, fung,keepNA = T)
  }
  
  saveWorkbook(wb, "./Suppl_dataset_S4.xlsx", overwrite = TRUE)
  
  
  
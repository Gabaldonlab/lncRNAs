library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)
library(WGCNA)
disableWGCNAThreads()
library(ggrepel)
library(GO.db)
library(clusterProfiler)
library(patchwork)
library(ggpubr)
library(VennDiagram)
library(ggh4x)
library(tidyr)
library(taRifx)
library(openxlsx)
library(colorspace)


### specify here the path to the ploting_and_network_analysis. For example, we specified it for our file system.
### note that for ggsave commands we specify a directory of convenience, so you will need to change them (can be easily done with search/replace)

wd<-"~/Dropbox/lncRNA/for_github/ploting_and_network_analysis/"


setwd(wd)
#############################
#### Import lncRNA data #####
#############################


for (spp in c("calb","ctrop","cpar","cglab","caur")){
  
  bed_initial<-read.table(paste0(spp,"/",spp,".bed"),header = F)
  colnames(bed_initial)<-c("Chr","Start","End","TransID","Strand")
  
  prot_cod_ids <- read.table(paste0(spp,"/",spp,"_prot_cod_ids.txt"),header = F)
  prot_cod_ids <-as.character(prot_cod_ids$V1)
  assign(paste0(spp,"_prot_cod_ids"),prot_cod_ids)
  dim(bed_initial)
  
  only_lncRNA<-bed_initial[grepl("MSTRG",bed_initial$TransID),]
  only_known_feature<-bed_initial[!grepl("MSTRG",bed_initial$TransID),]
  only_known_feature_separated<-cbind(only_known_feature,data.frame(do.call('rbind',strsplit(as.character(bed_initial[!grepl("MSTRG",bed_initial$TransID),]$TransID),"|",fixed=TRUE)))[1])
  
  
  only_prot_cod<-only_known_feature_separated[as.character(only_known_feature_separated$X1) %in% prot_cod_ids,]
  
  bed<-rbind(only_lncRNA,only_prot_cod[,c(1:5)])
  
  bed$Length<-abs(bed$End-bed$Start)+1
  bed$Spp<-spp
  bed$Type<-ifelse(grepl("\\|x",bed$TransID), "x", ifelse(grepl("\\|u",bed$TransID),"u","pc"))
  
  
  if (spp=="ctrop"){
    
    to_discard <- read.table("ctrop/ctrop_to_discard.txt")
    to_discard <- as.character(to_discard$V1)
    bed <- bed[!bed$TransID %in% to_discard,]
    
  }
  
  else if(spp=="caur"){
    to_discard <- read.table("caur/caur_to_discard.txt")
    to_discard <- as.character(to_discard$V1)
    
    bed <- bed[!bed$TransID %in% to_discard,]
  }
  
  assign(spp,bed)
  
  lncRNA_bed<-bed[grep("^MSTR", bed$TransID), ]
  assign(paste0(spp,"_lncRNA"),lncRNA_bed)
  
  
  lncRNA_bed_splitted<-data.frame(do.call('rbind', strsplit(as.character(lncRNA_bed$TransID),'|',fixed=TRUE)))
  colnames(lncRNA_bed_splitted)<-c("ID","Length","Class_code","Gene","Chr","Spp")
  lncRNA_bed_splitted$ID <- paste0(lncRNA_bed_splitted$ID,"_",spp)
  
  assign(paste0(spp,"_lncRNA_splitted"),lncRNA_bed_splitted)
}

all_data<-rbind(calb,cglab,cpar,ctrop,caur)
all_data_lncRNA<-rbind.data.frame(calb_lncRNA_splitted,cglab_lncRNA_splitted,cpar_lncRNA_splitted,ctrop_lncRNA_splitted, caur_lncRNA_splitted)



#################################
#### Plot numbers of lncRNAs ####
#################################


num_u_x<-all_data_lncRNA %>%
  dplyr::count("Class_code"=all_data_lncRNA$Class_code, "Spp"=all_data_lncRNA$Spp)

num_u_x$Spp <- factor(num_u_x$Spp, levels = c("calb","ctrop","cpar","caur","cglab"))


num_u_x$Species <- ifelse(num_u_x$Spp=="calb","C. albicans",ifelse(num_u_x$Spp=="ctrop","C. tropicalis",ifelse(num_u_x$Spp=="cpar","C. parapsilosis",ifelse(num_u_x$Spp=="caur","C. auris","C. glabrata"))))
num_u_x$Species <- factor(num_u_x$Species, levels = c("C. albicans", "C. tropicalis","C. parapsilosis", "C. auris", "C. glabrata"))

num_u_x$Class_code_2<-ifelse(num_u_x$Class_code=="u","i","a")
num_u_x$Class_code_2<-factor(num_u_x$Class_code_2,levels = c("i","a"))



number_plot<-ggplot(num_u_x, aes(x=Class_code_2, y=n, fill=Class_code_2)) + facet_grid(Species ~ ., switch = "y")+
  geom_bar(position="stack", stat="identity")+
  geom_text(aes(label=n),size = 5, position = position_stack(vjust = 0.5))+theme_bw()+
  scale_fill_manual(values = c("grey","steelblue2"))+guides(fill=guide_legend(title="Class code"))+
  ylab("Number of lncRNAs")+xlab("")+scale_x_discrete(position = "top")+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size=16),
        strip.text = element_text(size=16, face = "italic"),
        legend.position = "none")+coord_flip()

#################################################
#### Plot length of lncRNAs and coding genes ####
#################################################

all_data$Spp <- factor(all_data$Spp, levels = c("calb","ctrop","cpar","caur","cglab"))
all_data$Type_2<-ifelse(all_data$Type=="u","i",ifelse(all_data$Type=="pc","pc","a"))
all_data$Type_2<-factor(all_data$Type_2,levels = c("pc","i","a"))



all_data$Species <- ifelse(all_data$Spp=="calb","C. albicans",ifelse(all_data$Spp=="ctrop","C. tropicalis",ifelse(all_data$Spp=="cpar","C. parapsilosis",ifelse(all_data$Spp=="caur","C. auris","C. glabrata"))))
all_data$Species <- factor(all_data$Species, levels = c("C. albicans", "C. tropicalis","C. parapsilosis", "C. auris", "C. glabrata"))



length_plot<-ggplot(all_data, aes(x=Type_2, y=log2(Length),fill=Type_2))+geom_boxplot()+facet_grid(. ~ Species)+
  theme_bw()+ylab("\n\nLog2(Transcript length)")+guides(fill=guide_legend(title="Class code"))+
  theme(axis.text = element_text(size=16),axis.title = element_text(size=13),
        #axis.text.x=element_blank(),
        legend.text = element_text(size=16),legend.title = element_text(size=16),
        strip.text = element_text(size=16, face = "italic"),
        axis.title.x=element_blank())+
  scale_fill_manual(values = c("mediumorchid3","grey","steelblue2"))+theme(legend.position = "none")

#### p-values for length
length_pval<-function(spp){
  u<<-all_data[all_data$Spp==spp & all_data$Type == "u",]
  x<<-all_data[all_data$Spp==spp & all_data$Type == "x",]
  p<<-all_data[all_data$Spp==spp & all_data$Type == "pc",]
  
  u_vs_p<-wilcox.test(u$Length, p$Length)
  x_vs_p<-wilcox.test(x$Length, p$Length)
  x_vs_u<-t.test(x$Length, u$Length)
  
  print(paste0("For ",spp," u_vs_p is "))
  print(u_vs_p$p.value)
  print(paste0("For ",spp," x_vs_p is "))
  print(x_vs_p$p.value)
  print(paste0("For ",spp," x_vs_u is "))
  print(x_vs_u$p.value)
}

length_pval("calb")
length_pval("cglab")
length_pval("cpar")
length_pval("ctrop")
length_pval("caur")


#########################
#### Plot GC-content ####
#########################


GC_content_all<-NULL
for (spp in c("calb","ctrop","cpar","caur","cglab")){
  gc_content<-read.table(paste0(spp,"/",spp,"_gc.tsv"),header = F)
  
  gc_content[] <- lapply(gc_content, function(x) gsub("-T", "", x))
  
  gc_lncRNAs <- gc_content[grepl("MSTRG",gc_content$V1),]
  gc_known_features <- gc_content[!grepl("MSTRG",gc_content$V1),]
  
  gc_only_prot_cod <- gc_known_features[gc_known_features$V1 %in% eval(parse(text=paste0(spp,"_prot_cod_ids"))),]
  
  
  gc_content_both <- rbind(gc_lncRNAs,gc_only_prot_cod)
  
  
  if (spp=="ctrop"){
    
    to_discard <- read.table("ctrop/ctrop_to_discard.txt")
    to_discard <- as.character(to_discard$V1)
    gc_content_both <- gc_content_both[!gc_content_both$V1 %in% to_discard,]
    
  }
  
  else if(spp=="caur"){
    to_discard <- read.table("caur/caur_to_discard.txt")
    to_discard <- as.character(to_discard$V1)
    
    gc_content_both <- gc_content_both[!gc_content_both$V1 %in% to_discard,]
  }
  
  gc_intergenic<-read.table(paste0(spp,"/",spp,"_intergenic_gc.tsv"),header = F)
  
  GC_content_all<-rbind(GC_content_all,gc_content_both, gc_intergenic)
  
  
}

colnames(GC_content_all)<-c("gene","Class_code","GC","spp")
GC_content_all$GC <- as.numeric(GC_content_all$GC)

GC_content_all$spp <- factor(GC_content_all$spp, levels = c("calb","ctrop","cpar","caur","cglab"))
GC_content_all$spp<-factor(GC_content_all$spp,levels = c("calb","ctrop","cpar","caur","cglab"))


GC_content_all$Class_code_2<-ifelse(GC_content_all$Class_code=="u","i",ifelse(GC_content_all$Class_code=="pc","pc",
                                                                              ifelse(GC_content_all$Class_code=="x","a","inter")))
GC_content_all$Class_code_2<-factor(GC_content_all$Class_code_2,levels = c("pc","i","a","inter"))


GC_content_all$Species <- ifelse(GC_content_all$spp=="calb","C. albicans",ifelse(GC_content_all$spp=="ctrop","C. tropicalis",ifelse(GC_content_all$spp=="cpar","C. parapsilosis",ifelse(GC_content_all$spp=="caur","C. auris","C. glabrata"))))
GC_content_all$Species <- factor(GC_content_all$Species, levels = c("C. albicans", "C. tropicalis","C. parapsilosis", "C. auris", "C. glabrata"))




GC_plot<-ggplot(GC_content_all, aes(x=Class_code_2, y=log2(GC),fill=Class_code_2))+geom_boxplot()+facet_grid(. ~ Species)+
  theme_bw()+ylab("\n\nLog2(GC content)")+guides(fill=guide_legend(title="Class code"))+
  theme(axis.text = element_text(size = 16),axis.title = element_text(size=13),
        axis.title.x = element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        legend.text = element_text(size=16),legend.title = element_text(size=16),
        strip.text = element_text(size=16, face = "italic"))+scale_fill_manual(values = c("mediumorchid3","grey","steelblue2","red"))+
  ylim(-3,-0.5)+theme(legend.position = "none")


### pvalue of gc content

gc_pval<-function(spp){
  u<<-GC_content_all[GC_content_all$spp==spp & GC_content_all$Class_code == "u",]
  x<<-GC_content_all[GC_content_all$spp==spp & GC_content_all$Class_code == "x",]
  p<<-GC_content_all[GC_content_all$spp==spp & GC_content_all$Class_code == "pc",]
  inter<<-GC_content_all[GC_content_all$spp==spp & GC_content_all$Class_code == "inter",]
  
  u_vs_p<-wilcox.test(u$GC, p$GC)
  x_vs_p<-wilcox.test(x$GC, p$GC)
  x_vs_u<-wilcox.test(x$GC, u$GC)
  
  x_vs_inter <-wilcox.test(x$GC, inter$GC)
  i_vs_inter <-wilcox.test(u$GC, inter$GC)
  
  print(paste0("For ",spp," u_vs_p is "))
  print(u_vs_p$p.value)
  print(paste0("For ",spp," x_vs_p is "))
  print(x_vs_p$p.value)
  print(paste0("For ",spp," x_vs_u is "))
  print(x_vs_u$p.value)
  
  print(paste0("For ",spp," x_vs_inter is "))
  print(x_vs_inter$p.value)
  
  print(paste0("For ",spp," i_vs_inter is "))
  print(i_vs_inter$p.value)
  
  
}

gc_pval("calb")
gc_pval("cglab")
gc_pval("cpar")
gc_pval("ctrop")
gc_pval("caur")


###########################################
#### Plot expression and make the PCAs #### ####
###########################################


analyze_exp <- function(spp){
  
  expression<-as.matrix(read.table(sprintf("./%s/%s_expression.txt",spp,spp),check.names = FALSE))
  
  length<-data.frame(do.call('rbind', strsplit(as.character(row.names(expression)),'|',fixed=TRUE)))
  norm_exp_len <- expression/as.numeric(as.character(length[,2]))
  TPM_initial <- t(t(norm_exp_len) * 1e6 / colSums(norm_exp_len))
  
  
  ### for networks
  
  ids <- colnames(expression)
  colnames(TPM_initial)<-ids
  assign(paste0(spp,"_TPM"),TPM_initial, envir = .GlobalEnv)
  
  ####  
  
  TPM_initial<-cbind.data.frame("ID"=rownames(TPM_initial),TPM_initial)
  
  
  
  
  only_lncRNA<-TPM_initial[grepl("MSTRG",TPM_initial$ID),]
  only_known_feature<-TPM_initial[!grepl("MSTRG",TPM_initial$ID),]
  only_known_feature_separated<-cbind(data.frame(do.call('rbind',strsplit(as.character(only_known_feature$ID),"|",fixed=TRUE)))[1], only_known_feature)
  
  
  only_prot_cod<-only_known_feature_separated[as.character(only_known_feature_separated$X1) %in% eval(parse(text=paste0(spp,"_prot_cod_ids"))),]
  
  TPM<-rbind.data.frame(only_lncRNA[,2:ncol(only_lncRNA)], only_prot_cod[,3:ncol(only_prot_cod)])
  
  
  
  TPM_yeast<-TPM  ### for co-expression analysis
  ids <- colnames(expression)
  colnames(TPM_yeast)<-ids

  
  
  expression_lncrna<-expression[grepl("MSTR",rownames(expression)),]
  TPM_yeast_lncrna<-TPM_yeast[grepl("MSTR",rownames(TPM_yeast)),]
  
  
  TPM_mean<-as.data.frame(cbind("Mean"=rowMeans(TPM), "Spp"=spp))
  
  TPM_mean$Type<-ifelse(grepl("\\|x",row.names(TPM_mean)), "x", ifelse(grepl("\\|u",row.names(TPM_mean)),"u","pc"))
  
  expression_long<-melt(TPM_mean, id.vars=c("Spp","Type"))
  
  assign(paste0(spp,"_expression_long"),expression_long,envir = .GlobalEnv)
  
  #### makes pca plot for all species
  
  TPM_log<-log2(TPM_yeast+0.01)
  pca<-prcomp(t(TPM_log))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  
  
  metadata<-read.table(sprintf("./%s/%s_SraRunTable.txt",spp,spp),sep = "\t", fill=TRUE,header = T)
  metadata_filt<-metadata[metadata$Run %in% str_replace_all(colnames(TPM_yeast[,grepl("RR",colnames(TPM_yeast))]),"_counts",""),]
  
  
  
  N_samples_bigRNASeq<-length(colnames(TPM_yeast[,!grepl("RR",colnames(TPM_yeast))]))
  study<-c(c(rep("S dataset",N_samples_bigRNASeq)),as.character(metadata_filt$SRA_Study))
  read_lenght<-c(c(rep(75,N_samples_bigRNASeq),metadata_filt$AvgSpotLen))

  
  df_plotting<-cbind.data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],study,read_lenght)
  
  if (spp=="calb"){
    
    N_samples_bigRNASeq<-length(colnames(TPM_yeast[,!grepl("RR",colnames(TPM_yeast))]))
    N_samples_others<-length(colnames(TPM_yeast[,grepl("RR",colnames(TPM_yeast))]))
    study<-c(c(rep("S dataset",N_samples_bigRNASeq)),c(rep("B dataset",N_samples_others)))
    
    df_plotting<-cbind.data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],study,read_lenght)
    df_plotting$study <- factor(df_plotting$study, levels = c("S dataset", "B dataset"))
    pca_plot_calb<-ggplot(df_plotting[order(df_plotting$study,decreasing = T),])+
      geom_point(aes(PC1,PC2,colour=study),size=2)+
      xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) + 
      theme_bw()+
      labs(colour="Dataset")+
      theme(axis.text = element_text(size = 14),axis.title = element_text(size=14),
            legend.text=element_text(size=10))
    

    assign(paste0(spp,"_PCA"),pca_plot_calb,envir = .GlobalEnv)
    
  }
  else if (spp == "caur"){
    pca_plot_caur<-ggplot(df_plotting, aes(PC1,PC2,label=rownames(df_plotting)))+geom_point(aes(colour=as.factor(read_lenght),shape=study), size=4)+
      xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) + theme_bw()+labs(colour="Read length", shape="Project")+
      guides(color = guide_legend(order=1), shape = guide_legend(order=2))+
      theme(axis.text = element_text(size = 14),axis.title = element_text(size=14),
            legend.text=element_text(size=10))+
      guides(color = guide_legend(override.aes = list(size = 4)))+
      guides(shape = guide_legend(override.aes = list(size = 4)))

    assign(paste0(spp,"_PCA"),pca_plot_caur,envir = .GlobalEnv)
  }
  
  else{
    pca_plot<-ggplot(df_plotting, aes(PC1,PC2,label=rownames(df_plotting)))+geom_point(aes(colour=study,shape=as.factor(read_lenght)), size=4)+
      xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) + theme_bw()+labs(colour="Project", shape="Read length")+
      guides(color = guide_legend(order=1), shape = guide_legend(order=2))+
      theme(axis.text = element_text(size = 14),axis.title = element_text(size=14),
            legend.text=element_text(size=10))+
      guides(color = guide_legend(override.aes = list(size = 4)))+
      guides(shape = guide_legend(override.aes = list(size = 4)))

    assign(paste0(spp,"_PCA"),pca_plot,envir = .GlobalEnv)
    
  }
}

analyze_exp("calb")
analyze_exp("ctrop")
analyze_exp("cpar")
analyze_exp("caur")
analyze_exp("cglab")

#### PCA plots
calb_PCA+ctrop_PCA+cpar_PCA+caur_PCA+cglab_PCA+plot_layout(ncol = 2)+plot_annotation(tag_levels = "a") +
  ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/pca_plots_composite.pdf", units="in", width=14, heigh=16, dpi=600)

calb_PCA+ctrop_PCA+cpar_PCA+caur_PCA+cglab_PCA+plot_layout(ncol = 2)+plot_annotation(tag_levels = "a") +
  ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/pca_plots_composite.png", units="in", width=14, heigh=16, dpi=600)
####




all_data_expression_long<-rbind(calb_expression_long,cglab_expression_long,cpar_expression_long,ctrop_expression_long, caur_expression_long)

all_data_expression_long$Spp <- factor(all_data_expression_long$Spp, levels = c("calb","ctrop","cpar","caur","cglab"))

all_data_expression_long$Species <- ifelse(all_data_expression_long$Spp=="calb","C. albicans",ifelse(all_data_expression_long$Spp=="ctrop","C. tropicalis",ifelse(all_data_expression_long$Spp=="cpar","C. parapsilosis",ifelse(all_data_expression_long$Spp=="caur","C. auris","C. glabrata"))))
all_data_expression_long$Species <- factor(all_data_expression_long$Species, levels = c("C. albicans", "C. tropicalis","C. parapsilosis", "C. auris", "C. glabrata"))



all_data_expression_long$Type_2<-ifelse(all_data_expression_long$Type=="u","i",ifelse(all_data_expression_long$Type=="pc","pc","a"))
all_data_expression_long$Type_2<-factor(all_data_expression_long$Type_2,levels = c("pc","i","a"))

expression_plot<-ggplot(all_data_expression_long, aes(x=Type_2, y=log2(as.numeric(as.character(value))+0.01),fill=Type_2))+geom_boxplot()+facet_grid(. ~ Species)+
  theme_bw()+ylab("\n\nLog2(mean TPM+0.01)")+guides(fill=guide_legend(title="Class code"))+
  theme(axis.text = element_text(size = 16),axis.title = element_text(size=13),
        legend.text = element_text(size=16),legend.title = element_text(size=16),
        legend.position = "none",
        strip.text = element_text(size=16, face = "italic"),
        axis.title.x=element_blank())+scale_fill_manual(values = c("mediumorchid3","grey","steelblue2"))

#### p-values expression

expression_pval<-function(spp){
  u<<-all_data_expression_long[all_data_expression_long$Spp==spp & all_data_expression_long$Type == "u",]
  x<<-all_data_expression_long[all_data_expression_long$Spp==spp & all_data_expression_long$Type == "x",]
  p<<-all_data_expression_long[all_data_expression_long$Spp==spp & all_data_expression_long$Type == "pc",]
  
  u_vs_p<-wilcox.test(as.numeric(as.character(u$value)), as.numeric(as.character(p$value)))
  x_vs_p<-wilcox.test(as.numeric(as.character(x$value)), as.numeric(as.character(p$value)))
  x_vs_u<-wilcox.test(as.numeric(as.character(x$value)), as.numeric(as.character(u$value)))
  
  print(paste0("For ",spp," u_vs_p is "))
  print(u_vs_p$p.value)
  print(paste0("For ",spp," x_vs_p is "))
  print(x_vs_p$p.value)
  print(paste0("For ",spp," x_vs_u is "))
  print(x_vs_u$p.value)
}


expression_pval("calb")
expression_pval("cglab")
expression_pval("cpar")
expression_pval("ctrop")

                                                      ########################
                                                      #### Plot SNPs data ####
                                                      ########################
#Due to very large size, the data for this plot is not availabe at github directory. Request it by email from grant.hovhannisyan@gmail.com

for (spp in c("calb","ctrop","cpar","cglab","caur")){
  
  snps_genic_and_lncrna<-read.table(sprintf("~/users/tg/current/hhovhannisyan/jena_project/lncRNAs_bigRNAseq_and_public/SNP_analysis/%s_snp_counts_genic_and_lncRNA.txt",spp))
  
  snps_genic_and_lncrna$V8<-snps_genic_and_lncrna$V7/(snps_genic_and_lncrna$V3-snps_genic_and_lncrna$V2+1)
  snps_genic_and_lncrna$V9<-ifelse(grepl("\\|x",snps_genic_and_lncrna$V4), "a", ifelse(grepl("\\|u",snps_genic_and_lncrna$V4),"i","pc"))
  
  colnames(snps_genic_and_lncrna)<-c("Chromosome","Start","End","ID","Strand","Num","N_variants","N_varians_per_length","Type")
  
  only_lncRNA<-snps_genic_and_lncrna[grepl("MSTRG",snps_genic_and_lncrna$ID),]
  only_known_feature<-snps_genic_and_lncrna[!grepl("MSTRG",snps_genic_and_lncrna$ID),]
  only_known_feature_separated<-cbind(data.frame(do.call('rbind',strsplit(as.character(only_known_feature$ID),"|",fixed=TRUE)))[1], only_known_feature)
  only_prot_cod<-only_known_feature_separated[as.character(only_known_feature_separated$X1) %in% eval(parse(text=paste0(spp,"_prot_cod_ids"))),]
  snps_pc_and_lncRNAs<-rbind.data.frame(only_lncRNA,only_prot_cod[,2:ncol(only_prot_cod)])
  
  
  intergenic<-read.table(sprintf("~/users/tg/current/hhovhannisyan/jena_project/lncRNAs_bigRNAseq_and_public/SNP_analysis/%s_snp_counts_intergenic.txt",spp))
  
  intergenic<-cbind(intergenic[1:3],paste0("intergenic_",intergenic$V4),"NA",intergenic$V4,intergenic$V5,intergenic$V5/(intergenic$V3-intergenic$V2),"inter")
  colnames(intergenic)<-c("Chromosome","Start","End","ID","Strand","Num","N_variants","N_varians_per_length","Type")
  
  snps<-rbind(snps_pc_and_lncRNAs,intergenic)
  snps$Spp<-spp
  
  
  if (spp=="ctrop"){
    
    to_discard <- read.table("~/users/tg/current/hhovhannisyan/jena_project/lncRNAs_bigRNAseq_and_public/revision/blast/ctrop_to_discard.txt")
    to_discard <- as.character(to_discard$V1)
    snps <- snps[!snps$ID %in% to_discard,]
    
  }
  
  else if(spp=="caur"){
    to_discard <- read.table("~/users/tg/current/hhovhannisyan/jena_project/lncRNAs_bigRNAseq_and_public/revision/blast/caur_to_discard.txt")
    to_discard <- as.character(to_discard$V1)
    
    snps <- snps[!snps$ID %in% to_discard,]
  }
  
  
  #snps$Type<-factor(allsnps$Type,levels = c("pc","i","a","inter"))
  
  assign(paste0(spp,"_snps"),snps)
}

all_snps<-rbind(calb_snps,cglab_snps, cpar_snps,caur_snps,ctrop_snps)
all_snps$Type<-factor(all_snps$Type,levels = c("pc","i","a","inter"))
all_snps$Spp<-factor(all_snps$Spp,level=c("calb","ctrop","cpar","caur","cglab"))


all_snps$Species <- ifelse(all_snps$Spp=="calb","C. albicans",ifelse(all_snps$Spp=="ctrop","C. tropicalis",ifelse(all_snps$Spp=="cpar","C. parapsilosis",ifelse(all_snps$Spp=="caur","C. auris","C. glabrata"))))
all_snps$Species <- factor(all_snps$Species, levels = c("C. albicans", "C. tropicalis","C. parapsilosis", "C. auris", "C. glabrata"))
snps_to_plot<-all_snps[all_snps$N_varians_per_length>0,]

snps_plot<-ggplot(snp_data, aes(x=Type, y=log2(N_varians_per_length),fill=Type))+geom_boxplot()+facet_grid(. ~ Species)+
  theme_bw()+ylab("\n\nLog2(N varians/Region length)")+guides(fill=guide_legend(title="Class code"))+
  theme(axis.text = element_text(size=16),axis.title = element_text(size=13),
        legend.text = element_text(size=16),legend.title = element_text(size=16),
        legend.position = "none",
        strip.text = element_text(size=16, face = "italic"),
        axis.title.x=element_blank())+
  scale_fill_manual(values = c("mediumorchid3","grey","steelblue2","red"))+theme(legend.position = "none")

### p-value snps
snps_pval<-function(spp){
  i<<-all_snps[all_snps$Spp==spp & all_snps$Type == "i",]
  a<<-all_snps[all_snps$Spp==spp & all_snps$Type == "a",]
  p<<-all_snps[all_snps$Spp==spp & all_snps$Type == "pc",]
  inter<<-all_snps[all_snps$Spp==spp & all_snps$Type == "inter",]
  
  i_vs_p<-wilcox.test(as.numeric(as.character(i$N_varians_per_length)), as.numeric(as.character(p$N_varians_per_length)))
  a_vs_p<-wilcox.test(as.numeric(as.character(a$N_varians_per_length)), as.numeric(as.character(p$N_varians_per_length)))
  a_vs_i<-wilcox.test(as.numeric(as.character(a$N_varians_per_length)), as.numeric(as.character(i$N_varians_per_length)))
  
  i_vs_inter<-wilcox.test(as.numeric(as.character(i$N_varians_per_length)), as.numeric(as.character(inter$N_varians_per_length)))
  a_vs_inter<-wilcox.test(as.numeric(as.character(a$N_varians_per_length)), as.numeric(as.character(inter$N_varians_per_length)))
  p_vs_inter<-wilcox.test(as.numeric(as.character(p$N_varians_per_length)), as.numeric(as.character(inter$N_varians_per_length)))
  
  print(paste0("For ",spp," i_vs_p is "))
  print(i_vs_p$p.value)
  print(paste0("For ",spp," a_vs_p is "))
  print(a_vs_p$p.value)
  print(paste0("For ",spp," a_vs_i is "))
  print(a_vs_i$p.value)
  
  
  print(paste0("For ",spp," i_vs_inter is "))
  print(i_vs_inter$p.value)
  print(paste0("For ",spp," a_vs_inter is "))
  print(a_vs_inter$p.value)
  print(paste0("For ",spp," p_vs_inter is "))
  print(a_vs_inter$p.value)
  
}


snps_pval("calb")
snps_pval("ctrop")
snps_pval("cglab")
snps_pval("cpar")
snps_pval("caur")

#### combine all plots

length_data<-cbind(all_data[,c(6,10,9)],"Measure"="Length")
colnames(length_data)<-c("value","spp","class_code","measure")


gc_data <- cbind(GC_content_all[,c(3,6,5)], "Measure"= "GC content")
colnames(gc_data)<-c("value","spp","class_code","measure")
gc_data<-gc_data[gc_data$val>0.1,]


expr_data<-cbind(all_data_expression_long[,c(4,5,6)],"Measure"="Expression")
colnames(expr_data)<-c("value","spp","class_code","measure")


variant_data<-cbind(snps_to_plot[,c(8,11,9)],"Measure"="Variants")
colnames(variant_data)<-c("value","spp","class_code","measure")

all_data_to_plot<-rbind.data.frame(length_data,gc_data,expr_data,variant_data)

all_data_to_plot$class_code<-as.character(all_data_to_plot$class_code)

all_data_to_plot[all_data_to_plot == "inter"] <- "ir"
all_data_to_plot$class_code <- factor(all_data_to_plot$class_code, c("pc","i","a","ir"))

all_data_to_plot$measure <- factor(all_data_to_plot$measure, c("Length","Expression","GC content","Variants"))


all_data_plot <- ggplot(all_data_to_plot, aes(x=class_code, y=log2(as.numeric(value)),fill=class_code))+geom_boxplot()+facet_grid(measure ~ spp,scales = "free")+
  theme_bw()+ylab("")+
  guides(fill=guide_legend(title="Class code"))+
  theme(axis.text = element_text(size=16),axis.title = element_text(size=13),
        legend.text = element_text(size=16),legend.title = element_text(size=16),
        legend.position = "none",
        strip.text.x = element_text(size=16, face = "italic"),
        strip.text.y = element_text(size=16),
        axis.title.x=element_blank())+
  scale_fill_manual(values = c("mediumorchid3","grey","steelblue2","red"))+theme(legend.position = "none")


number_plot + all_data_plot+plot_layout(widths = c(1,3)) + 
  ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/composite.pdf",device = "pdf", units="in", width=13, heigh=9, dpi=600)

number_plot + all_data_plot+plot_layout(widths = c(1,3)) + 
  ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/composite.png",device = "png", units="in", width=13, heigh=9, dpi=600)







                                  #####################################################################
                                  #### Network analysis - Choose a set of soft-thresholding powers ####
                                  #####################################################################
library(psych)

choose_soft_threshold_and_plot_results<-function(TPM_table,spp,TPM_filt,percent_sample_filt,SRP_to_exclude){
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  
  SRR_SRP<-read.table(sprintf("./%s/SRR_SRP_correspondence_%s.txt",spp,spp), sep = "\t", header = T)
  
  SRR_to_exclude<<-as.character(SRR_SRP[SRR_SRP$SRA_Study %in% SRP_to_exclude,]$Run)
  TPM_table_excl<-TPM_table[,setdiff(colnames(TPM_table),SRR_to_exclude),drop = FALSE]
  
  #TPM_table_excl<-TPM_table_excl[!grepl("MSTRG",rownames(TPM_table_excl)),]
  print(colnames(TPM_table_excl))
  
  TPM_table_filt<-TPM_table_excl[rowSums(TPM_table_excl >= TPM_filt) >= percent_sample_filt*ncol(TPM_table_excl),]
  TPM_table_filt<-TPM_table_filt[!grepl("\\|x\\|",rownames(TPM_table_filt)),]
  sft_spp = pickSoftThreshold(t(log2(TPM_table_filt+1)), powerVector = powers, verbose = 5,networkType = "unsigned")
  
  # Plot the results:
  # Scale-free topology fit index as a function of the soft-thresholding power
  scale_plot<-ggplot(sft_spp$fitIndices,aes(sft_spp$fitIndices[,1], -sign(sft_spp$fitIndices[,3])*sft_spp$fitIndices[,2]))+
    geom_point(size=2)+geom_text(label=sft_spp$fitIndices$Power,hjust=0.5, vjust=-1, size=6)+geom_hline(yintercept = 0.8, color="red")+
    xlab("Soft Threshold (power)")+ylab("Scale Free Topology Model Fit,signed R^2")+ggtitle(sprintf("Scale Independence (%s)\nTPM>=%s in %s%% of samples",spp,TPM_filt,percent_sample_filt*100))+theme_bw()+
    theme(plot.title = element_text(hjust = 0.5,size = 20),
          axis.text = element_text(size = 18),axis.title = element_text(size=18))+
    scale_y_continuous(breaks=seq(-1, 1.1, 0.2))+expand_limits(y=c(-1,1.1))
  
  
  connectivity_plot<-ggplot(sft_spp$fitIndices,aes(sft_spp$fitIndices[,1],  sft_spp$fitIndices[,5]))+
    geom_point(size=2)+geom_text(label=sft_spp$fitIndices$Power,hjust=0.5, vjust=-1, size=6)+
    xlab("Soft Threshold (power)")+ylab("Mean Connectivity")+ggtitle(sprintf("Mean Connectivity (%s)\nTPM>=%s in %s%% of samples",spp,TPM_filt,percent_sample_filt*100))+theme_bw()+
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.text = element_text(size = 18),axis.title = element_text(size=18))+
    scale_y_continuous(breaks=seq(0, 3000, 400))+expand_limits(y=c(0, 3000))
  
  ggarrange(scale_plot, connectivity_plot,ncol = 2, nrow = 1)+
    ggsave(sprintf("./%s/coexpression/%s_soft_threshold_TPM_more_or_equal_to_%s_in_%s_prcnt_of_samples.png",spp,spp,TPM_filt,percent_sample_filt*100),units="in", width=12, heigh=6, dpi=300 )
  
}

choose_soft_threshold_and_plot_results(calb_TPM,"calb",0.1,0.8,c())
choose_soft_threshold_and_plot_results(ctrop_TPM,"ctrop",0.1,0.8,c("SRP083839","SRP099169"))
choose_soft_threshold_and_plot_results(cpar_TPM,"cpar",0.1,0.8, c("SRP151798","SRP041812"))
choose_soft_threshold_and_plot_results(cglab_TPM,"cglab",0.1,0.8, c("SRP065276"))
choose_soft_threshold_and_plot_results(caur_TPM,"caur",0.1,0.8, c())


                                    #################################################################
                                    #### Building network and modules, and run GO, PFAM and KEGG ####
                                    #################################################################

build_network_and_modules<-function(TPM_table,spp,soft_threshold,TPM_filt,percent_sample_filt, SRP_to_exclude){
  
  
  SRR_SRP<-read.table(sprintf("./%s/SRR_SRP_correspondence_%s.txt",spp,spp), sep = "\t", header = T)
  
  SRR_to_exclude<-as.character(SRR_SRP[SRR_SRP$SRA_Study %in% SRP_to_exclude,]$Run)
  TPM_table_excl<-TPM_table[,setdiff(colnames(TPM_table),SRR_to_exclude),drop = FALSE]
  
  
  ### Filtering
  TPM_table_filt<-TPM_table_excl[rowSums(TPM_table_excl >= TPM_filt) >= percent_sample_filt*ncol(TPM_table_excl),]
  TPM_table_filt<-TPM_table_filt[!grepl("\\|x\\|",rownames(TPM_table_filt)),]
  TPM_table_filt_out<<-TPM_table_filt
  
  
  #### calculate adjacency, TOM, cluster tree and plot
  adjacency = adjacency(t(log2(TPM_table_filt+1)), power = soft_threshold)
  TOM = TOMsimilarity(adjacency)
  
  dissTOM = 1-TOM
  geneTree = hclust(as.dist(dissTOM), method = "average")
  
  ###### Find modules
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2,
                              pamRespectsDendro = FALSE,minClusterSize = 30)
  table(dynamicMods)
  
  #Plot module assignment
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  
  
  ### Merge similar modules
  
  # Calculate eigengenes
  MEList = moduleEigengenes(t(log2(TPM_table_filt+1)), colors = dynamicColors,subHubs=TRUE,excludeGrey = T)
  MEs = MEList$eigengenes
  
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs)
  
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average")
  
  
  #Plot the cut line into the dendrogram
  MEDissThres = 0.25
  
  # Call an automatic merging function
  merge = mergeCloseModules(t(log2(TPM_table_filt+1)), dynamicColors, cutHeight = MEDissThres, verbose = 3)
  
  # The merged module colors
  mergedColors = merge$colors
  
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs
  pdf(sprintf("./%s/coexpression/%s_modules_TPM_more_or_equal_to_%s_in_%s_prcnt_of_samples.pdf",spp,spp,TPM_filt,percent_sample_filt*100),width = 12,height = 5)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Gene\nclustering", sprintf("Eigengene\nclustering")),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05, main=sprintf("Cluster Dendrogram (%s)", spp))
  dev.off()
  
  # Rename to moduleColors
  moduleColors = mergedColors
  moduleColors<<-moduleColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", standardColors(50))
  moduleLabels = match(moduleColors, colorOrder)-1
  MEs = mergedMEs
  
  ### Selecting hub genes
  
  top_hubs<-chooseTopHubInEachModule(t(log2(TPM_table_filt+1)), colorh=moduleColors, power = soft_threshold ,type = "unsigned")
  
  assign(paste0(spp,"_top_hubs"),top_hubs,envir = .GlobalEnv)
  
  
  ### Connectivity
  
  iK=intramodularConnectivity.fromExpr(t(log2(TPM_table_filt+1)), power=soft_threshold, colors=moduleColors)
  row.names(iK)=row.names(TPM_table_filt)
  
  iK$Type<-ifelse(grepl("\\|u",row.names(iK)), "u","pc")
  
  assign(paste0(spp,"_connectivity"),iK,envir = .GlobalEnv)
  
  ggplot(iK, aes(x=Type, y=log2(kWithin+0.01), fill=Type)) +
    geom_boxplot()+
    theme_bw()+ylab("Log2(intermodular connectivity)")+guides(fill=guide_legend(title="Class code"))+
    ggtitle(spp)+
    theme(axis.text.y = element_text(size = 22),axis.title = element_text(size=22),
          axis.text.x=element_blank(),
          legend.text = element_text(size=22),legend.title = element_text(size=22),
          strip.text = element_text(size=22),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank())+scale_fill_manual(values = c("mediumorchid3","grey"))+
    ggsave(sprintf("./%s/coexpression/%s_connectivity_TPM_more_or_equal_to_%s_in_%s_prcnt_of_samples.png",spp,spp,TPM_filt,percent_sample_filt*100), width = 10, height = 7,dpi = 300)
  
  if (spp=="calb"){
    
    
    ### calb darkred
    
    darkred_ids<-names(as.data.frame((t(log2(TPM_table_filt+1)))))[moduleColors=="darkred" ]
    module_matr<-is.finite(match(moduleColors,c("darkred")))
    darkred_matr<-TOM[module_matr,module_matr]
    dimnames(darkred_matr)<-list(darkred_ids,darkred_ids)
    
    
    cyt_calb = exportNetworkToCytoscape(darkred_matr,
                                        edgeFile = "./calb/coexpression/edges_darkred_calb.txt",
                                        nodeFile = "./calb/coexpression/node_darkred_calb.txt",
                                        weighted = TRUE,
                                        nodeNames = darkred_ids,
                                        nodeAttr = moduleColors[module_matr],
                                        threshold = 0)
    
    
  }
  
  else if (spp=="ctrop"){
    
    ### ctrop coral1
    coral1_ids<-names(as.data.frame((t(log2(TPM_table_filt+1)))))[moduleColors=="coral1" ]
    coral1_ids_out<<-names(as.data.frame((t(log2(TPM_table_filt+1)))))[moduleColors=="coral1" ]
    module_matr<-is.finite(match(moduleColors,c("coral1")))
    coral1_matr<-TOM[module_matr,module_matr]
    dimnames(coral1_matr)<-list(coral1_ids,coral1_ids)
    
    
    cyt_ctrop = exportNetworkToCytoscape(coral1_matr,
                                         edgeFile = "./ctrop/coexpression/edge_coral1_ctrop_new.txt",
                                         nodeFile = "./ctrop/coexpression/node_coral1_ctrop_new.txt",
                                         weighted = TRUE,
                                         nodeNames = coral1_ids,
                                         nodeAttr = moduleColors[module_matr],
                                         threshold = 0)
    
    
  }
  else{
    print("no cytoscape network for this species")
  }
  
  
  ### count number of lncRNA in modules
  module_stats<-NULL
  all_lncRNAs_in_modules<-NULL
  all_lncRNAs_in_modules_id_module<-NULL
  for (module in unique(moduleColors)){
    if(!(module=="grey")){
      #print(module)
      
      ## names is the same as colnames.
      gene_ids<-names(as.data.frame((t(log2(TPM_table_filt+1)))))[moduleColors==module]
      
      
      
      
      assign(paste0("gene_ids_",module),gene_ids,envir = .GlobalEnv)
      
      module_stats<-rbind(module_stats,
                          c(module,
                            as.numeric(length(gene_ids)),
                            as.numeric(length(grep("MSTRG", gene_ids))),
                            round(length(grep("MSTRG", gene_ids))/length(gene_ids)*100,2),spp))
      all_lncRNAs_in_modules<-c(all_lncRNAs_in_modules,gene_ids[grepl("MSTRG", gene_ids)])
      
      if (length(gene_ids[grepl("MSTRG", gene_ids)])>0){
        all_lncRNAs_in_modules_id_module<-rbind.data.frame(all_lncRNAs_in_modules_id_module,cbind.data.frame(gene_ids[grepl("MSTRG", gene_ids)],module))
      }
      
      #print(paste(module,length(gene_ids),length(grep("MSTRG", gene_ids))))
      
    }
  }
  all_lncRNAs_in_modules_unique<-unique(all_lncRNAs_in_modules)
  assign(paste0(spp,"_all_lncRNAs_in_modules"),all_lncRNAs_in_modules_unique, envir = .GlobalEnv)
  assign(paste0(spp,"_all_lncRNAs_in_modules_id_module"),all_lncRNAs_in_modules_id_module, envir = .GlobalEnv)


  ######### GO terms  PFAM  and KEGG ############
 
  all_go_terms<-as.data.frame(GOTERM)

  all_go_terms<-all_go_terms[!(all_go_terms$go_id %in% c("GO:0003674","GO:0008150","GO:0005575")),]


  MF_terms<-all_go_terms$go_id[all_go_terms$Ontology=="MF"]
  BP_terms<-all_go_terms$go_id[all_go_terms$Ontology=="BP"]
  CC_terms<-all_go_terms$go_id[all_go_terms$Ontology=="CC"]


  ### Association of GOID with description
  goterms <- Term(GOTERM)
  a<-as.data.frame(goterms)
  go_names<<-cbind(row.names(a),a)


  ### Yeast GOIDS
  yeast_go<-read.table(sprintf("./%s/%s_go.txt",spp,spp))


  MF_universe<-yeast_go[yeast_go$V1 %in% MF_terms,]
  BP_universe<-yeast_go[yeast_go$V1 %in% BP_terms,]
  CC_universe<<-yeast_go[yeast_go$V1 %in% CC_terms,]

  go_module_kegg_module_table <- NULL
   for (module in unique(moduleColors)){
    if (!(module=="grey")){

      gene_ids<-names(as.data.frame((t(log2(TPM_table_filt+1)))))[moduleColors==module]
      gene_ids_no_lncRNA<-sapply(strsplit(gene_ids[!grepl("^MSTRG",gene_ids)],"\\|"), "[", 1 )

      gene_ids_lncRNA<-sapply(gene_ids[grepl("^MSTRG",gene_ids)], "[", 1 )

      print(module)

      for (ontology in c("MF","CC","BP")){

        ego<-enricher(gene_ids_no_lncRNA, pvalueCutoff = 0.05,
                      pAdjustMethod = "BH", universe = eval(parse(text=paste0("as.character(",ontology,"_universe$V2)"))),minGSSize = 2,
                      maxGSSize = NA, TERM2GENE = eval(parse(text=paste0(ontology,"_universe"))),TERM2NAME = go_names)
        if (!is.null(ego)){
          p<-dotplot(ego,showCategory = 5)
          p+ggplot2::ggsave(sprintf("./%s/coexpression/go_terms/%s_module_%s_%s.png",spp,spp,module,ontology),
                            units="in", width=10, height =7, dpi=600)
          write.table(ego,sprintf("./%s/coexpression/go_terms/%s_module_%s_%s.txt",spp,spp,module,ontology), sep = "\t")

          assign(sprintf("%s_module_%s_%s.txt",spp,module,ontology),as.data.frame(ego), envir = .GlobalEnv)

          ### just output the ids of lncrnas in each module
          write.table(gene_ids_lncRNA,sprintf("./%s/coexpression/go_terms/%s_module_%s_%s_lncrna_ids.txt",spp,spp,module,ontology),quote = F)
        }
        ### this will be plotted on barplot
        
        if (length(gene_ids_lncRNA) > 0){
          
          lncRNAs <- toString(as.character(gene_ids_lncRNA))
        }
        else{
          lncRNAs <- "NA"
        }
        
        if (ontology=="BP"){
          if (dim(as.data.frame(ego))[1]==0){
            go_module<-cbind.data.frame("No enrichment",module,toString(as.character(gene_ids_no_lncRNA)) ,lncRNAs)
          }
          else if (dim(ego)[1]<3){
            go_module<-cbind.data.frame(toString(as.character(ego[,2])),module,toString(as.character(gene_ids_no_lncRNA)),lncRNAs)
          }
          else{
            go_module<-cbind.data.frame(toString(as.character(ego[c(1:5),2])),module,toString(as.character(gene_ids_no_lncRNA)),lncRNAs)
            
            #print(module_go_terms)
          }
        }
      }
 
  #print(go_module)

      if (spp=="calb"){
        kegg_spp = "cal"
      }
      else if (spp == "ctrop"){
        kegg_spp = "ctp"
      }
      else if (spp=="caur"){
        kegg_spp = "caur"
      }
      else if (spp == "cglab"){
        kegg_spp = "cgr"
      }
      else {
        print("No KEGG for CPAR")
      }
      
      
      
      if (spp=="calb"){
        gene_ids_no_lncRNA_for_kegg <- gsub("^C","CAALFM_C",gsub("_","",gsub("_A","A",gene_ids_no_lncRNA)))
      }
      else if (spp=="caur"){
        caur_orth = read.table("./caur/caur_strain_orthologs.tsv",header = F)
        gene_ids_no_lncRNA_for_kegg <- as.character(caur_orth[caur_orth$V1 %in% gene_ids_no_lncRNA,][,2])
      }
      else{
        gene_ids_no_lncRNA_for_kegg <- gene_ids_no_lncRNA
      }

      
             
      if (spp!="cpar"){
      ego_kegg <- enrichKEGG(gene = gene_ids_no_lncRNA_for_kegg,
                       organism = kegg_spp,
                       pvalueCutoff = 0.05)

      if (dim(as.data.frame(ego_kegg))[1]==0){
        kegg_module<-cbind.data.frame("No enrichment",module)
      }
      else if (dim(as.data.frame(ego_kegg))[1]<3){
        kegg_module<-cbind.data.frame(toString(as.character(ego_kegg[,2])),module)
      }
      else{
        kegg_module<-cbind.data.frame(toString(as.character(ego_kegg[c(1:5),2])),module)

      }
      }

      else{
        kegg_module<-cbind.data.frame("NA",module)
      }
      
      #### PFAM
    
      
      if (spp=="calb"){
        pfam_name = "C_albicans_SC5314_iprscan.out"
      }
      else if (spp == "ctrop"){
        pfam_name = "C_tropicalis_MYA-3404_iprscan.out"
      }
      else if (spp=="caur"){
        pfam_name = "C_auris_B8441_iprscan.out"
      }
      else if (spp == "cglab"){
        pfam_name = "C_glabrata_CBS138_iprscan.out"
      }
      else {
        pfam_name = "C_parapsilosis_CDC317_iprscan.out"
      }
      
      pfam <- read.table(sprintf("./%s/%s",spp,pfam_name),header = F,fill = T, sep = "\t", quote = "\"")
      universe_pfam <- pfam[,c(12,1)]
      universe_pfam <- unique(universe_pfam[universe_pfam$V12!="NULL",])
      
      
      pfam_name <- pfam[,c(12,13)]
      pfam_name <- unique(pfam_name[pfam_name$V12!="NULL",])

      ego_pfam<-enricher(gene_ids_no_lncRNA, pvalueCutoff = 0.05,
                    pAdjustMethod = "BH", universe = as.character(universe_pfam$V1),minGSSize = 2,
                    maxGSSize = 10000, TERM2GENE = universe_pfam,TERM2NAME = pfam_name)
      
      if (dim(as.data.frame(ego_pfam))[1]==0){
        pfam_module<-cbind.data.frame("No enrichment",module)
      }
      else if (dim(as.data.frame(ego_pfam))[1]<3){
        pfam_module<-cbind.data.frame(toString(as.character(ego_pfam[,2])),module)
      }
      else{
        pfam_module<-cbind.data.frame(toString(as.character(ego_pfam[c(1:5),2])),module)
      }
      
      
      
      colnames(go_module) <- c("go","module", "genes","lncRNAs")
      go_module_out<<-go_module
      
      colnames(kegg_module) <- c("kegg","module")
      colnames(pfam_module) <- c("pfam","module")
      #print(cbind.data.frame(go_module,kegg_module))
      #print(kegg_module)
      #colnames(go_module_kegg_module) <-  c("go","module","kegg","module2")
      go_module_kegg_module_table <- rbind.data.frame(go_module_kegg_module_table,cbind.data.frame(go_module,pfam_module,kegg_module))
     }
   }
  ############

   colnames(go_module_kegg_module_table)<-c("GO terms","Module","Gene ids","lncRNA ids","PFAM","Module2","KEGG","Module3")
   colnames(module_stats)<-c("Module","N_genes","N_lncRNAs","prcnt_lncRNAs","Spp")
   
   module_stats<-as.data.frame(module_stats)
   

   module_stats_sort<-merge(module_stats,go_module_kegg_module_table, by=c("Module"))
  
   module_stats_sort_final <<- module_stats_sort[order(as.numeric(as.character(module_stats_sort$N_genes)),decreasing = F),]
   print(module_stats_sort)
   
   module_stats_sort_final$Module<-factor(module_stats_sort_final$Module, levels=as.character(module_stats_sort_final$Module))
   
   
   assign(paste0(spp,"_module_stats_sort_final"),module_stats_sort_final, envir = .GlobalEnv)
   
   ### plots without GO terms
   ggplot(data=module_stats_sort_final, aes(x=Module,y=as.numeric(as.character(N_genes)),width=0.9, fill=Module,
                                            label = paste0(module_stats_sort_final$N_lncRNAs,"(",module_stats_sort_final$prcnt_lncRNAs,"%)"))) +
     geom_bar(stat="identity")+scale_fill_manual(values = as.character(module_stats_sort_final$Module))+
     geom_text(hjust=0, size=5)+theme_bw()+coord_flip()+
     ylab("Number of genes in modules")+xlab(sprintf("Modules (%s)",spp))+
     theme(legend.position = "none",
           axis.text = element_text(size = 18),axis.title = element_text(size=18),
           axis.text.y = element_blank(),
           axis.ticks.y = element_blank(),
           axis.text.x = element_text(angle = 90,vjust=0.5, hjust=1),
           plot.title = element_text(hjust = 0.5))+
     scale_y_continuous(breaks=seq(0, max(as.numeric(as.character(module_stats_sort_final$N_genes)))+300, 200))+
     expand_limits(y=c(0, max(as.numeric(as.character(module_stats_sort_final$N_genes)))+300))+
     scale_x_discrete(labels=as.character(module_stats_sort_final$Module))+
     ggsave(sprintf("./%s/coexpression/%s_N_genes_in_modules_TPM_more_or_equal_to_%s_in_%s_prcnt_of_samples.pdf",spp,spp,TPM_filt,percent_sample_filt*100), width = 12, height = 6,dpi = 300)
   
   
   assign(paste0("TOM_",spp),TOM,envir = .GlobalEnv)
   

  }
ctrop_top_hubs
cpar_top_hubs
caur_top_hubs
cglab_top_hubs


ctrop_TPM[coral1_ids_out,]

build_network_and_modules(calb_TPM,"calb",12,0.1,0.8,c(""))
build_network_and_modules(ctrop_TPM,"ctrop",16,0.1,0.8,c("SRP083839","SRP099169"))
build_network_and_modules(cpar_TPM,"cpar",14,0.1,0.8,c("SRP151798","SRP041812"))
build_network_and_modules(cglab_TPM,"cglab",12,0.1,0.8, c("SRP065276"))
build_network_and_modules(caur_TPM,"caur",12,0.1,0.8, c(""))


                                      ##########################################################
                                      #### Make a supplemenatry file with GO, PFAM and KEGG ####
                                      ##########################################################




all_go_terms_for_modules<-rbind.data.frame(calb_module_stats_sort_final,
                                           ctrop_module_stats_sort_final,
                                           cpar_module_stats_sort_final,
                                           caur_module_stats_sort_final,
                                           cglab_module_stats_sort_final)


ctrop_module_stats_sort_final

#all_go_terms_for_modules<-b
all_go_terms_for_modules<-all_go_terms_for_modules[,c(1,2,3,4,5,6,9,11,7,8)]
colnames(all_go_terms_for_modules) <- c("Module name","Total number of genes","Number of lncRNAs","% of lncRNAs","Species","GO term enrichment","PFAM domain enrichment","KEGG pathway enrichment","IDs of known features","lncRNA IDs")

all_go_terms_for_modules$Species <- ifelse(all_go_terms_for_modules$Species=="calb","C. albicans",ifelse(all_go_terms_for_modules$Species=="ctrop","C. tropicalis",ifelse(all_go_terms_for_modules$Species=="cpar","C. parapsilosis",ifelse(all_go_terms_for_modules$Species=="caur","C. auris","C. glabrata"))))


wb_S2 <- createWorkbook("Supp_table_S7.xlsx")
addWorksheet(wb_S2, "Supp_table_S7")
writeData(wb_S2, "Supp_table_S7", all_go_terms_for_modules,keepNA = T)
saveWorkbook(wb_S2, "./Supp_table_S7.xlsx", overwrite = TRUE)


                                                  ##############################
                                                  #### Analyze connectivity ####
                                                  ##############################


connectivity_all<-rbind.data.frame(cbind.data.frame(calb_connectivity,"spp"="calb"),
                                   cbind.data.frame(ctrop_connectivity,"spp"="ctrop"),
                                   cbind.data.frame(cglab_connectivity,"spp"="cglab"),
                                   cbind.data.frame(caur_connectivity,"spp"="caur"),
                                   cbind.data.frame(cpar_connectivity,"spp"="cpar"))


connectivity_all$Type<-ifelse(grepl("MSTR", rownames(connectivity_all)),"i","pc")
connectivity_all$Type <- factor(connectivity_all$Type, levels = c("pc","i"))
connectivity_all$spp<-factor(connectivity_all$spp, levels = c("calb","ctrop","cpar","caur","cglab"))

connectivity_all$Species <- ifelse(connectivity_all$spp=="calb","C. albicans",ifelse(connectivity_all$spp=="ctrop","C. tropicalis",ifelse(connectivity_all$spp=="cpar","C. parapsilosis",ifelse(connectivity_all$spp=="caur","C. auris","C. glabrata"))))
connectivity_all$Species <- factor(connectivity_all$Species, levels = c("C. albicans", "C. tropicalis","C. parapsilosis", "C. auris", "C. glabrata"))



ggplot(connectivity_all, aes(x=Type, y=log2(kTotal+0.01),fill=Type))+geom_boxplot()+facet_grid(. ~ Species)+
  theme_bw()+ylab("Log2(Total connectivity+0.01)")+guides(fill=guide_legend(title="Class code"))+
  theme(axis.text = element_text(size = 20),axis.title = element_text(size=16),
        legend.text = element_text(size=22),legend.title = element_text(size=22),
        #legend.position = "none",
        strip.text = element_text(size=22, face = "italic"),
        axis.title.x=element_blank(),
        legend.position = "none")+scale_fill_manual(values = c("mediumorchid3","grey"))+
  ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/connectivity_all.png", width = 12, height = 8,dpi = 300)
  #ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/connectivity_all.pdf", width = 12, height = 8,dpi = 300)



#### Connectivity p-values

connectivity_pval<-function(df,spp){
  
  df$Type<-ifelse(grepl("MSTR", rownames(df)),"i","pc")
  u<<-df[df$Type == "i",]
  #x<<-df[df$Type == "x",]
  p<<-df[df$Type == "pc",]
  #a[grepl("AB",a)]
  u_vs_p<-t.test(as.numeric(as.character(u$kWithin)), as.numeric(as.character(p$kWithin)))
  #x_vs_p<-t.test(as.numeric(as.character(x$kWithin)), as.numeric(as.character(p$kWithin)))
  #x_vs_u<-t.test(as.numeric(as.character(x$kWithin)), as.numeric(as.character(u$kWithin)))
  spp<-"calb"
  print(paste0("For ",spp," u_vs_p is "))
  print(u_vs_p)
  #print(paste0("For ",spp," x_vs_p is "))
  #print(x_vs_p)
  #print(paste0("For ",spp," x_vs_u is "))
  #print(x_vs_u)
}

connectivity_pval(calb_connectivity,"calb")
connectivity_pval(ctrop_connectivity,"ctrop")
connectivity_pval(cpar_connectivity,"cpar")
connectivity_pval(caur_connectivity,"caur")
connectivity_pval(cglab_connectivity,"cglab")


######### connectivity of inf-specific lncRNAs
connectivity_plot_all<-NULL
for (spp in c("calb","ctrop","cpar","cglab")){
  
  expression<-read.table(sprintf("./%s/infection_lncRNAs_%s_24_24c.txt",spp,spp))
  inf_spec_genes<-rownames(expression)
  assign(paste0(spp,"_inf_spec_genes_data"),expression)
  assign(paste0(spp,"_inf_spec_genes"),inf_spec_genes)
  
  connectivity<-eval(parse(text=paste0(spp,"_connectivity")))
  
  lncRNA_connectivity<-connectivity[grepl("MSTRG",rownames(connectivity)),]
  
  lncRNA_connectivity<-lncRNA_connectivity[complete.cases(lncRNA_connectivity),]
  
  inf_spec_connectivity<-lncRNA_connectivity[rownames(lncRNA_connectivity) %in% as.character(inf_spec_genes),]
  other_connectivity<-lncRNA_connectivity[setdiff(rownames(lncRNA_connectivity),rownames(inf_spec_connectivity)),]
  connectivity_plot<-rbind(cbind(inf_spec_connectivity,"status"="Inf-specific","spp"=spp),
                           cbind(other_connectivity,"status"="Not Inf-specific","spp"=spp))
  
  print(spp)
  wilcox<-wilcox.test(connectivity_plot$kWithin[connectivity_plot$status=="Inf-specific"],
                      connectivity_plot$kWithin[connectivity_plot$status=="Not Inf-specific"])$p.value
  print(wilcox)
  
  connectivity_plot_all<-rbind.data.frame(connectivity_plot_all,connectivity_plot)
  
}



connectivity_plot_all$status<-factor(connectivity_plot_all$status, levels = c("Inf-specific","Not Inf-specific"))
connectivity_plot_all$spp<-factor(connectivity_plot_all$spp, levels = c("calb","ctrop","cpar","cglab"))

ggplot(connectivity_plot_all, aes(x=status, y=log2(kTotal+0.01),fill=status))+geom_boxplot()+facet_grid(. ~ spp)+
  theme_bw()+ylab("Log2(Total connectivity+0.01)")+guides(fill=guide_legend(title="lncRNA type"))+
  theme(axis.text.y = element_text(size = 22),axis.title = element_text(size=22),
        axis.text.x=element_blank(),
        legend.text = element_text(size=22),legend.title = element_text(size=22),
        #legend.position = "none",
        strip.text = element_text(size=22),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank())+scale_fill_manual(values = c("tan1","slateblue1"))

wilcox.test(connectivity_plot_all$kWithin[connectivity_plot_all$status=="Inf-specific" & connectivity_plot_all$spp=="cglab"],
            connectivity_plot_all$kWithin[connectivity_plot_all$status=="Not Inf-specific" & connectivity_plot_all$spp=="cglab"])


### check if inf-spec lncRNA are in families

all_lncRNA_families<-read.table("./intergenic.fam")

colnames(all_lncRNA_families)<-c("FAM","Gene")



calb_famlily_genes<-intersect(calb_inf_spec_genes,as.character(all_lncRNA_families$Gene))
calb_families<-as.character(all_lncRNA_families[all_lncRNA_families$Gene %in% calb_famlily_genes,][,1])
all_lncRNA_families[all_lncRNA_families$FAM %in% calb_families ,]

ctrop_famlily_genes<-intersect(ctrop_inf_spec_genes,as.character(all_lncRNA_families$Gene))
ctrop_families<-as.character(all_lncRNA_families[all_lncRNA_families$Gene %in% ctrop_famlily_genes,][,1])
all_lncRNA_families[all_lncRNA_families$FAM %in% ctrop_families ,]

cpar_famlily_genes<-intersect(cpar_inf_spec_genes,as.character(all_lncRNA_families$Gene))
cpar_families<-as.character(all_lncRNA_families[all_lncRNA_families$Gene %in% cpar_famlily_genes,][,1])
all_lncRNA_families[all_lncRNA_families$FAM %in% cpar_families ,]

cglab_famlily_genes<-intersect(cglab_inf_spec_genes,as.character(all_lncRNA_families$Gene))
cglab_families<-as.character(all_lncRNA_families[all_lncRNA_families$Gene %in% cglab_famlily_genes,][,1])
all_lncRNA_families[all_lncRNA_families$FAM %in% cglab_families ,]



#### shared families of inf-specific lncRNAs

venn.diagram(list("calb"=calb_families,
                  "cpar"=cpar_families,
                  "ctrop"=ctrop_families,
                  "cglab"=cglab_families),
             fill=c("red", "blue", "green", "orange"), height = 5, width = 5, resolution = 600,
             filename= "./lncRNA_families_inf_spec_venn.svg",
             cat.cex=0.8, cex=2, main.cex = 1,
             imagetype = "svg")

venn_plot<-draw.quintuple.venn(area1=length(calb_families), 
                               area2=length(cpar_families), 
                               area3=length(ctrop_families), 
                               area4=length(calb_families),
                               category = c("CALB", "CPAR", "CTROP","CAUR" ,"CGLAB"),cex=2, 
                               n12=overlapalbpar, n13=overlapalbtrop, n14=overlapalbaur, n15=overlapalbglab, n23=overlappartrop, n24=overlapparaur, n25=overlapparglab, n34=overlaptropaur, n35=overlaptropglab, n45=overlapaurglab, n123=overlapalbpartrop, n124=overlapalbparaur, n125=overlapalbparglab, n134=overlapalbtropaur, n135=overlapalbtropglab, n145=overlapalbaurglab, n234=overlappartropaur, n235=overlappartropglab, n245=overlapparaurglab,n345=overlaptropaurglab, n1234=overlapalbpartropaur, n1235=overlapalbpartropglab, n1245=overlapalbparaurglab, n1345=overlapalbtropaurglab, n2345=overlappartropaurglab, n12345=overlapalbpartropaurglab,
                               fill = c("red", "blue", "green", "orange"))
pdf("genes.pdf")
grid.draw(venn_plot)
dev.off()


intersect(calb_families,ctrop_families)
intersect(calb_families,cpar_families)
intersect(ctrop_families,cpar_families)



### are inf_spec lncRNAs in modules
percent_lncRNAs_in_modules<- function(lncRNAs_in_modules,inf_spec_genes,total,spp,id_module){
  N_i <- length(inf_spec_genes[grepl("\\u",inf_spec_genes)])
  
  print(N_i)
  percent<-length(lncRNAs_in_modules[lncRNAs_in_modules %in% inf_spec_genes])/N_i*100
  percent_total<-length(lncRNAs_in_modules[lncRNAs_in_modules %in% inf_spec_genes])/total*100
  
  
  sprintf("%s%% of %s inf-spec lncRNAs are in network modules and  %s%% of total intergenic lncrnas are infection-specific and involved in modules",round(percent,2), spp, round(percent_total,2))

}

percent_lncRNAs_in_modules(calb_all_lncRNAs_in_modules,calb_inf_spec_genes,1459,"calb",calb_all_lncRNAs_in_modules_id_module) # 7216
percent_lncRNAs_in_modules(ctrop_all_lncRNAs_in_modules,ctrop_inf_spec_genes,1568,"ctrop",ctrop_all_lncRNAs_in_modules_id_module) # 2568
percent_lncRNAs_in_modules(cpar_all_lncRNAs_in_modules,cpar_inf_spec_genes,1499,"cpar",cpar_all_lncRNAs_in_modules_id_module) # 4534
percent_lncRNAs_in_modules(cglab_all_lncRNAs_in_modules,cglab_inf_spec_genes,449,"cglab",cglab_all_lncRNAs_in_modules_id_module) #1403



### genes of fam147 in modules

#all_lncRNA_families[all_lncRNA_families$FAM == "fam147" ,]



                                                   #############################
                                                   #### Saturation plots   #####
                                                   #############################


calb_satur<-read.table("calb/saturation_data_calb.txt",header = F)
cglab_satur<-read.table("cglab/saturation_data_cglab.txt",header = F)
cpar_satur<-read.table("cpar/saturation_data_cpar.txt",header = F)
ctrop_satur<-read.table("ctrop/saturation_data_ctrop.txt",header = F)
caur_satur<-read.table("caur/saturation_data_caur.txt",header = F)

all_saturation<-rbind.data.frame(calb_satur,cglab_satur,cpar_satur,ctrop_satur, caur_satur)

all_saturation$sampling<-ifelse(grepl("subseq",all_saturation$V1), "Subsequent", "Random")

colnames(all_saturation)<-c("sample_id","Species","N_samples","a","i","Sampling_strategy")


all_saturation <- gather(all_saturation, Type, N_lncRNAs, i:a, factor_key=TRUE)
all_saturation$Species <- factor(all_saturation$Species,levels = c("calb","ctrop","cpar","caur","cglab"))
all_saturation[all_saturation[,4]=="Random",]

all_saturation$Species <- ifelse(all_saturation$Species=="calb","C. albicans",ifelse(all_saturation$Species=="ctrop","C. tropicalis",ifelse(all_saturation$Species=="cpar","C. parapsilosis",ifelse(all_saturation$Species=="caur","C. auris","C. glabrata"))))
all_saturation$Species <- factor(all_saturation$Species, levels = c("C. albicans", "C. tropicalis","C. parapsilosis", "C. auris", "C. glabrata"))

all_saturation$Type <- ifelse(all_saturation$Type=="i","Intergenic","Antisense")
all_saturation$Type <- factor(all_saturation$Type,levels = c("Intergenic","Antisense"))


ggplot(all_saturation[all_saturation[,4]=="Random",],aes(x=N_samples,y=N_lncRNAs))+geom_line()+geom_point()+
  facet_nested(Type~Species, scales = "free")+
  theme(panel.spacing = unit(0,"line"))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 12, face = "italic"),
        strip.text.y = element_text(size = 12), 
        axis.text = element_text(size=12))+ylab("Number of lncRNAs")
  ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/saturation_plots.png", width = 12, height = 6,dpi = 600)
  #ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/saturation_plots.pdf", width = 12, height = 6,dpi = 600)



                                      #########################################################
                                      #### Plot distribution of lncRNAs across chromosomes ####
                                      #########################################################

library(chromPlot)
library(stringr)

rename<-function(spp,df){
  if (spp=="calb"){
    df$Chrom<-gsub('_C_albicans_SC5314', '', df$Chrom)
    df$Chrom<-gsub('Ca22chr', '', df$Chrom)  }
  else if (spp=="cglab"){
    df$Chrom<-gsub('_C_glabrata_CBS138','',df$Chrom)
    df$Chrom<-gsub('Chr','',df$Chrom)}
  else if (spp=="caur"){
    df$Chrom<-gsub('_C_auris_B8441','',df$Chrom)
    df$Chrom<-gsub('PEKT0200000','',df$Chrom)
    df$Chrom<-gsub('PEKT020000','',df$Chrom)}
  else if (spp=="cpar"){
    df$Chrom<-gsub('_C_parapsilosis_CDC317','',df$Chrom)
    df$Chrom<-gsub('Contig00','',df$Chrom)}
  else {
    df$Chrom<-gsub('_C_tropicalis_MYA-3404','',df$Chrom)
    df$Chrom<-gsub('Supercontig_3.','',df$Chrom)
    df<-df[df[,1] %in% c("1","2","3","4","5","6","7","8","9","10","11","12"),]
      }
  return(df)
}


plot_chromosomes<-function(spp,n_window,n_cols){

  chr<-read.table(sprintf("./%s/%s_chrNameLength.txt",spp,spp))
  chr<-chr[,c(1,3,4,5)]
  colnames(chr)<-c("Chrom","Start","End","Name")
  chr<<-remove.factors(chr)
  print(chr)
  
  chr<-rename(spp,chr)

  ### bed file
  bed<-read.table(sprintf("./%s/%s_lncRNAs.bed",spp,spp))
  
  if (spp=="ctrop"){
    
    to_discard <- read.table("./ctrop/ctrop_to_discard.txt")
    to_discard <- as.character(to_discard$V1)
    bed <- bed[!bed$V4 %in% to_discard,]
    
  }
  
  else if(spp=="caur"){
    to_discard <- read.table("./caur/caur_to_discard.txt")
    to_discard <- as.character(to_discard$V1)
    
    bed <- bed[!bed$V4 %in% to_discard,]
  }
  
  
  
  
  
  bed<-bed[,c(1,2,3,4)]
  colnames(bed)<-c("Chrom","Start","End","Name")
  
  
  bed<-rename(spp,bed)
  bed<<-remove.factors(bed)

  
  bed_u<-bed[grepl("\\|u",bed$Name),]
  bed_x<-bed[grepl("\\|x",bed$Name),]
  
  ### plot
  
  #,width = 6,height = 3,units = "in"
  pdf(sprintf("./%s/chrom_%s_wind_%s.pdf",spp,spp,n_window),width = 6,height = 3)
  chromPlot(gaps  = chr, bands  = cbind(chr,"Colors"="black"),
            annot1 = bed_x,annot2 = bed_u,colAnnot1 = "steelblue2",colAnnot2 = "grey",
            bin = n_window,chrSide=c(-1,1,-1,1,1,1,1,1), figCols=n_cols)
  dev.off()
}
plot_chromosomes("calb",5e4,9)
plot_chromosomes("cglab",5e4,14)
plot_chromosomes("cpar",5e4,9)
plot_chromosomes("ctrop",5e4,12)
plot_chromosomes("caur",5e4,15)


                                            ##########################################
                                            #### Plot lncRNA across chr distance #####
                                            ##########################################



all_distance <- NULL
for (spp in c("calb","ctrop","cpar","caur","cglab")){
  distance <- read.table(sprintf("./%s/%s_distance.bed",spp,spp))
  distance$V6 <- spp
  all_distance <- rbind.data.frame(all_distance,distance)
}

colnames(all_distance)<- c("chr","start","end","N_lncRNA","dist","spp")
all_distance$spp <- factor(all_distance$spp, levels = c("calb","ctrop","cpar","caur","cglab"))

filterred<-all_distance[all_distance$N_lncRNA>0 & all_distance$dist>0,]

filterred$dist <- as.character(filterred$dist)
#Then turn it back into a factor with the levels in the correct order
filterred$dist <- factor(filterred$dist, levels=unique(filterred$dist))


filterred$Species <- ifelse(filterred$spp=="calb","C. albicans",ifelse(filterred$spp=="ctrop","C. tropicalis",ifelse(filterred$spp=="cpar","C. parapsilosis",ifelse(filterred$spp=="caur","C. auris","C. glabrata"))))
filterred$Species <- factor(filterred$Species, levels = c("C. albicans", "C. tropicalis","C. parapsilosis", "C. auris", "C. glabrata"))


ggplot(filterred,aes(x=dist,y=N_lncRNA))+geom_boxplot()+facet_grid(Species~.,scales = "free")+  geom_line()+
  scale_y_continuous(breaks = seq(min(all_distance$N_lncRNA),max(all_distance$N_lncRNA),5))+
  # scale_x_continuous(breaks = seq(0,max(all_distance$dist),5e4))+
  xlab("Distance from the closest telomere")+ylab("Number of lncRNAs")+
  theme_bw()+theme(axis.text.x = element_text(angle = 75,vjust = 0.5),
                   strip.text.y = element_text(face = "italic"))+
  stat_summary(fun=mean, geom="line", aes(group=1,color="red"))  + theme(legend.position = "none")+
  ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/lncRNA_chr_distance.png",device = "png",height = 8,width = 6,dpi = 300,units = "in")
#ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/lncRNA_chr_distance.pdf",device = "pdf",height = 8,width = 6,dpi = 300,units = "in")



                                          ##############################################
                                          #### C. glabrata plots with known ncRNAs #####
                                          ##############################################


for (spp in c("cglab")){
  
  bed_initial<-read.table(paste0(spp,"/",spp,".bed"),header = F)
  colnames(bed_initial)<-c("Chr","Start","End","TransID","Strand")
  
  ## prot cod ids
  prot_cod_ids <- read.table(paste0(spp,"/",spp,"_prot_cod_ids.txt"),header = F)
  prot_cod_ids <-as.character(prot_cod_ids$V1)
  
  ## known ncRNA ids
  ncRNA <- read.table("./cglab/known_ncRNAs.txt")
  cglab_ncRNA_ids <- as.character(ncRNA$V1)
  
  
  
  only_lncRNA<-bed_initial[grepl("MSTRG",bed_initial$TransID),]
  only_lncRNA$Type <- ifelse(grepl("\\|x",only_lncRNA$TransID), "a", "i")
  
  
  only_known_feature<-bed_initial[!grepl("MSTRG",bed_initial$TransID),]
  only_known_feature_separated<-cbind(only_known_feature,data.frame(do.call('rbind',strsplit(as.character(bed_initial[!grepl("MSTRG",bed_initial$TransID),]$TransID),"|",fixed=TRUE)))[1])
  
  
  only_prot_cod<-only_known_feature_separated[as.character(only_known_feature_separated$X1) %in% prot_cod_ids,]
  
  only_known_ncRNAs <- only_known_feature_separated[as.character(only_known_feature_separated$X1) %in% cglab_ncRNA_ids,]
  only_known_ncRNAs_longer_200 <- only_known_ncRNAs[only_known_ncRNAs$End-only_known_ncRNAs$Start>200,]
  ncRNA_ids_longer_200 <- as.character(only_known_ncRNAs_longer_200$X1)
  
  
  bed<-rbind(only_lncRNA,cbind(only_prot_cod[,c(1:5)],"Type"="pc"),cbind(only_known_ncRNAs_longer_200[,c(1:5)],"Type"="known_ncRNAs"))
  
  bed$Length<-abs(bed$End-bed$Start)+1
  bed$Spp<-spp

}


bed$Type<-factor(bed$Type,levels = c("pc","i","a","known_ncRNAs"))

cglab_length_plot<-ggplot(bed, aes(x=Type, y=log2(Length),fill=Type))+geom_boxplot()+facet_grid(. ~ Spp)+
  theme_bw()+ylab("\n\nLog2(Transcript length)")+guides(fill=guide_legend(title="Class code"))+
  theme(axis.text = element_text(size=16),axis.title = element_text(size=13),
        #axis.text.x=element_blank(),
        legend.text = element_text(size=16),legend.title = element_text(size=16),
        strip.text = element_text(size=16),
        axis.title.x=element_blank())+
  scale_fill_manual(values = c("mediumorchid3","grey","steelblue2","green"))+theme(legend.position = "none")


### GC content

GC_content_cglab<-NULL
for (spp in c("cglab")){
  gc_content<-read.table(paste0(spp,"/",spp,"_gc.tsv"),header = F)
  
  gc_content[] <- lapply(gc_content, function(x) gsub("-T", "", x))
  
  gc_lncRNAs <- gc_content[grepl("MSTRG",gc_content$V1),]
  gc_lncRNAs$V2 <- ifelse(grepl("\\|x",gc_lncRNAs$V1), "a", "i")
  
  gc_known_features <- gc_content[!grepl("MSTRG",gc_content$V1),]
  gc_only_prot_cod <- gc_known_features[gc_known_features$V1 %in% prot_cod_ids,]
  gc_only_known_ncRNAs <- gc_known_features[as.character(gc_known_features$V1) %in% ncRNA_ids_longer_200,]
  gc_only_known_ncRNAs$V2 <- "known_ncRNAs"
  
  
  gc_content_both <- rbind(gc_lncRNAs,gc_only_prot_cod,gc_only_known_ncRNAs)
  
  gc_intergenic<-read.table(paste0(spp,"/",spp,"_intergenic_gc.tsv"),header = F)
  
  GC_content_cglab<-rbind(GC_content_cglab,gc_content_both, gc_intergenic)
  
  
}

colnames(GC_content_cglab)<-c("gene","Class_code","GC","spp")
GC_content_cglab$GC <- as.numeric(GC_content_cglab$GC)


GC_content_cglab$Class_code<-factor(GC_content_cglab$Class_code,levels = c("pc","i","a","inter","known_ncRNAs"))

GC_plot<-ggplot(GC_content_cglab, aes(x=Class_code, y=log2(GC),fill=Class_code))+geom_boxplot()+facet_grid(. ~ spp)+
  theme_bw()+ylab("\n\nLog2(GC content)")+guides(fill=guide_legend(title="Class code"))+
  theme(axis.text = element_text(size = 16),axis.title = element_text(size=13),
        axis.title.x = element_blank(),
        legend.text = element_text(size=16),legend.title = element_text(size=16),
        strip.text = element_text(size=16))+scale_fill_manual(values = c("mediumorchid3","grey","steelblue2","red","green"))+
  ylim(-4,1)+theme(legend.position = "none")



#### expression cglab


for (spp in c("cglab")){

  
  expression<-as.matrix(read.table(sprintf("./%s/%s_expression.txt",spp,spp),check.names = FALSE))
  
  
  length<-data.frame(do.call('rbind', strsplit(as.character(row.names(expression)),'|',fixed=TRUE)))
  norm_exp_len <- expression/as.numeric(as.character(length[,2]))
  TPM_initial <- t(t(norm_exp_len) * 1e6 / colSums(norm_exp_len))
  
  TPM_mean<-as.data.frame(cbind("Mean"=rowMeans(TPM_initial), "Spp"=spp))
  
  only_lncRNA<-TPM_mean[grepl("MSTRG",rownames(TPM_mean)),]
  only_lncRNA$Type <- ifelse(grepl("\\|x",rownames(only_lncRNA)), "a", "i")
  
  
  only_known_feature<-TPM_mean[!grepl("MSTRG",rownames(TPM_mean)),]
  only_known_feature <- cbind.data.frame(data.frame(do.call('rbind',strsplit(rownames(only_known_feature),"|",fixed=TRUE)))[1],only_known_feature)
  
  prot_cod <- cbind(only_known_feature[only_known_feature$X1 %in% prot_cod_ids,],"Type"="pc")
  known_ncRNAs <- cbind(only_known_feature[only_known_feature$X1 %in% ncRNA_ids_longer_200,],"Type"="known_ncRNAs")
  
  expression<-rbind(only_lncRNA,prot_cod[,2:4],known_ncRNAs[,2:4])
}


expression$Type<-factor(expression$Type,levels = c("pc","i","a","inter","known_ncRNAs"))


expression_plot<-ggplot(expression, aes(x=Type, y=log2(as.numeric(as.character(Mean))+0.01),fill=Type))+geom_boxplot()+facet_grid(. ~ Spp)+
  theme_bw()+ylab("\n\nLog2(mean TPM+0.01)")+guides(fill=guide_legend(title="Class code"))+
  theme(axis.text = element_text(size = 16),axis.title = element_text(size=13),
        legend.text = element_text(size=16),legend.title = element_text(size=16),
        legend.position = "none",
        strip.text = element_text(size=16),
        axis.title.x=element_blank())+scale_fill_manual(values = c("mediumorchid3","grey","steelblue2","green"))

#### SNPs data (available upon request)

for (spp in c("cglab")){
  
  snps_genic_and_lncrna<-read.table("~/users/tg/current/hhovhannisyan/jena_project/lncRNAs_bigRNAseq_and_public/SNP_analysis/cglab_snp_counts_genic_and_lncRNA.txt")
  
  snps_genic_and_lncrna$V8<-snps_genic_and_lncrna$V7/(snps_genic_and_lncrna$V3-snps_genic_and_lncrna$V2+1)
  #snps_genic_and_lncrna$V9<-ifelse(grepl("\\|x",snps_genic_and_lncrna$V4), "a", ifelse(grepl("\\|u",snps_genic_and_lncrna$V4),"i","pc"))
  
  #colnames(snps_genic_and_lncrna)<-c("Chromosome","Start","End","ID","Strand","Num","N_variants","N_varians_per_length","Type")
  
  only_lncRNA<-snps_genic_and_lncrna[grepl("MSTRG",snps_genic_and_lncrna$V4),]
  only_lncRNA$V9<-ifelse(grepl("\\|x",only_lncRNA$V4), "a", "i")
  
  
  only_known_feature<-snps_genic_and_lncrna[!grepl("MSTRG",snps_genic_and_lncrna$V4),]
  only_known_feature_separated<-cbind(data.frame(do.call('rbind',strsplit(as.character(only_known_feature$V4),"|",fixed=TRUE)))[1], only_known_feature)
  
  
  only_prot_cod<-cbind(only_known_feature_separated[as.character(only_known_feature_separated$X1) %in% prot_cod_ids,],"V9"="pc")
  known_ncRNAs<-cbind(only_known_feature_separated[as.character(only_known_feature_separated$X1) %in% ncRNA_ids_longer_200,],"V9"="known_ncRNAs")
  
  
  snps_pc_and_lncRNAs<-rbind.data.frame(only_lncRNA,only_prot_cod[,2:ncol(only_prot_cod)],known_ncRNAs[,2:ncol(known_ncRNAs)])
  colnames(snps_pc_and_lncRNAs)<-c("Chromosome","Start","End","ID","Strand","Num","N_variants","N_varians_per_length","Type")
  
  intergenic<-read.table(sprintf("~/users/tg/current/hhovhannisyan/jena_project/lncRNAs_bigRNAseq_and_public/SNP_analysis/%s_snp_counts_intergenic.txt",spp))
  
  intergenic<-cbind(intergenic[1:3],paste0("intergenic_",intergenic$V4),"NA",intergenic$V4,intergenic$V5,intergenic$V5/(intergenic$V3-intergenic$V2),"inter")
  colnames(intergenic)<-c("Chromosome","Start","End","ID","Strand","Num","N_variants","N_varians_per_length","Type")
  
  snps<-rbind(snps_pc_and_lncRNAs,intergenic)
  snps$Spp<-"cglab"
  
}



snps$Type<-factor(snps$Type,levels = c("pc","i","a","inter","known_ncRNAs"))

#snps_to_plot<-all_snps[all_snps$N_varians_per_length>0,]
snps_plot<-ggplot(snps, aes(x=Type, y=log2(N_varians_per_length),fill=Type))+geom_boxplot()+facet_grid(. ~ Spp)+
  theme_bw()+ylab("\n\nLog2(N varians/Region length)")+guides(fill=guide_legend(title="Class code"))+
  theme(axis.text = element_text(size=16),axis.title = element_text(size=13),
        legend.text = element_text(size=16),legend.title = element_text(size=16),
        legend.position = "none",
        strip.text = element_text(size=16),
        axis.title.x=element_blank())+
  scale_fill_manual(values = c("mediumorchid3","grey","steelblue2","green","red"))+theme(legend.position = "none")


len <- cbind(bed[,c(6:8)], "Measure"="Length")
colnames(len) <-  c("type","value","spp","measure")

expr <- cbind(expression[,c(3,1,2)],"Measurue"="Expression")
colnames(expr) <- c("type","value","spp","measure")

gc <- cbind(GC_content_cglab[,c(2,3,4)],"Measure"="GC content")
colnames(gc) <- c("type","value","spp","measure")

snips <- cbind(snps[,c(9,8,10)],"Measure"="Variants")
colnames(snips) <- c("type","value","spp","measure")


cglab_all <- rbind(len,expr,gc,snips)
cglab_all$measure <- factor(cglab_all$measure, c("Length","Expression","GC content","Variants"))

cglab_all$Species <- "C. glabrata"

ggplot(cglab_all, aes(x=type, y=log2(as.numeric(value)),fill=type))+geom_boxplot()+facet_grid(measure ~ Species,scales = "free")+
  theme_bw()+ylab("")+
  guides(fill=guide_legend(title="Class code"))+
  theme(axis.text = element_text(size=16),axis.title = element_text(size=13),
        legend.text = element_text(size=16),legend.title = element_text(size=16),
        legend.position = "none",
        strip.text.x = element_text(size=16, face = "italic"),
        strip.text.y = element_text(size=16),
        axis.title.x=element_blank())+
  scale_fill_manual(values = c("mediumorchid3","grey","steelblue2","red","green"))+theme(legend.position = "none")+
  #ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/composite_cglab.pdf",device = "pdf",width = 12, height = 14,dpi = 300)
  ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/composite_cglab.png",width = 12, height = 14,dpi = 300)




                                          #############################################
                                          #### Fold change of inf-specific lncRNAs ####
                                          #############################################

expression_inf_spec<-NULL
for (spp in c("ctrop","calb","cpar","cglab")){
  
  expression<-read.table(sprintf("./%s/infection_lncRNAs_%s_24_24c.txt",spp,spp))
  if (spp=="ctrop"){
    spp_name = "C. tropicalis"
  }
  else if (spp == "calb"){
    spp_name = "C. albicans"
  }
  else if (spp == "cpar"){
    spp_name = "C. parapsilosis"
  }
  else{
    spp_name = "C. glabrata"
  }
  expression_inf_spec<-rbind.data.frame(expression_inf_spec,cbind.data.frame(expression$log2FoldChange,spp_name))
  
}
colnames(expression_inf_spec)<-c("FC","Species")

expression_inf_spec$reg <- ifelse(expression_inf_spec$FC>0,"Up-regulation","Down-regulation")
expression_inf_spec$reg <- factor(expression_inf_spec$reg, levels = c("Up-regulation","Down-regulation"))
expression_inf_spec$Species <- factor(expression_inf_spec$Species, levels = c("C. albicans","C. tropicalis","C. parapsilosis","C. glabrata"))

summary(expression_inf_spec[expression_inf_spec$FC< 0,])

ggplot(expression_inf_spec, aes(x=Species, y=FC,fill=Species))+geom_boxplot()+facet_grid(reg ~ .,scales = "free")+ylab("Log2(fold-change)")+
  theme_bw()+ theme(axis.text = element_text(size=12),
                    axis.title = element_text(size=12),
                    legend.position = "none",
                    strip.text = element_text(size=12),
                    axis.text.x = element_text(face = "italic"))+
  ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/inf_spec_fold_change.png", units="in", width=6, heigh=6, dpi=600)
  #ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/inf_spec_fold_change.pdf", units="in", width=6, heigh=6, dpi=600)




                                                 ################################
                                                 #### Trinity vs Stringtie   ####
                                                 ################################




trinity_vs_stringtie_classcodes <- function(spp,N_samples){

  trinity_classcodes <- as.data.frame(read.table(sprintf("../mapping_assembly_lncRNA_prediction/trinity/%s/all_class_codes_trinity.txt",spp), header = F))
  trinity_classcodes$Mean <- trinity_classcodes$V1/N_samples
  trinity_classcodes$Species <- spp
  trinity_classcodes$Software <- "Trinity"
  
  stringtie_classcodes <- as.data.frame(read.table(sprintf("../mapping_assembly_lncRNA_prediction/trinity/%s/all_class_codes_stringtie.txt",spp), header = F))
  stringtie_classcodes$Mean <- stringtie_classcodes$V1/N_samples
  stringtie_classcodes$Species <- spp
  stringtie_classcodes$Software <- "Stringtie"
  
  merged <- rbind.data.frame(trinity_classcodes,stringtie_classcodes)
  if (spp=="ctrop"){
    merged$Species_to_plot = "C. tropicalis"
  }
  else if (spp == "calb"){
    merged$Species_to_plot = "C. albicans"
  }
  else if (spp == "cpar"){
    merged$Species_to_plot = "C. parapsilosis"
  }
  else if (spp == "caur"){
    merged$Species_to_plot = "C. auris"
  }
  else{
    merged$Species_to_plot = "C. glabrata"
  }
  
  colnames(merged) <- c("raw_classcode","Class code","Mean class codes","spp","Software","Species")
  return(merged)
}

class_codes_calb <- trinity_vs_stringtie_classcodes("calb",703)
class_codes_ctrop <- trinity_vs_stringtie_classcodes("ctrop",53)
class_codes_cpar <- trinity_vs_stringtie_classcodes("cpar",86)
class_codes_caur <- trinity_vs_stringtie_classcodes("caur",61)
class_codes_cglab <- trinity_vs_stringtie_classcodes("cglab",51)

class_codes_all_spp <- rbind.data.frame(class_codes_calb,class_codes_ctrop,class_codes_cpar,class_codes_caur,class_codes_cglab)
class_codes_all_spp$Species <- factor(class_codes_all_spp$Species, levels = c("C. albicans","C. tropicalis","C. parapsilosis","C. auris","C. glabrata"))


class_codes <- ggplot(class_codes_all_spp, aes(x=class_codes_all_spp$Software, y=class_codes_all_spp$`Mean class codes`, fill=class_codes_all_spp$`Class code`)) +
  geom_bar(stat="identity",colour="black") +
  facet_nested(~Species)+theme_bw()+
  ylab("Mean number of class codes across samples")+
  xlab("Software")+
  guides(fill=guide_legend(title="Class code"))+
  theme(strip.text = element_text(face = "italic",size = 12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))#+
#ggsave("~/Dropbox (Personal)/lncRNA/Nat_Com_submission/revision/plots/class_codes.png", units="in", width=10, heigh=6, dpi=600, device = "png")


#### sensitivity and precision

trinity_vs_stringtie_sensitivity <- function(spp){

  trinity_sens <- as.data.frame(read.table(sprintf("../mapping_assembly_lncRNA_prediction/trinity/%s/sensitivity_precision_trinity.txt",spp), header = F))
  trinity_sens$Species <- spp
  trinity_sens$Software <- "Trinity"
  
  stringtie_sens <- as.data.frame(read.table(sprintf("../mapping_assembly_lncRNA_prediction/trinity/%s/sensitivity_precision_stringtie.txt",spp), header = F))
  stringtie_sens$Species <- spp
  stringtie_sens$Software <- "Stringtie"
  
  merged <- rbind.data.frame(trinity_sens,stringtie_sens)
  
  if (spp=="ctrop"){
    merged$Species_to_plot = "C. tropicalis"
  }
  else if (spp == "calb"){
    merged$Species_to_plot = "C. albicans"
  }
  else if (spp == "cpar"){
    merged$Species_to_plot = "C. parapsilosis"
  }
  else if (spp == "caur"){
    merged$Species_to_plot = "C. auris"
  }
  else{
    merged$Species_to_plot = "C. glabrata"
  }
  merged_final <- merged[,c(3,5,7,8,9)]
  colnames(merged_final) <- c("Sensitivity","Precision","spp","Software","Species")
  return(merged_final)
}

sense_calb <- trinity_vs_stringtie_sensitivity("calb")
sense_ctrop <- trinity_vs_stringtie_sensitivity("ctrop")
sense_cpar <- trinity_vs_stringtie_sensitivity("cpar")
sense_caur <- trinity_vs_stringtie_sensitivity("caur")
sense_cglab <- trinity_vs_stringtie_sensitivity("cglab")



sense_all_spp <- rbind.data.frame(sense_calb,sense_ctrop,sense_cpar,sense_caur,sense_cglab)
sense_all_spp$Species <- factor(sense_all_spp$Species, levels = c("C. albicans","C. tropicalis","C. parapsilosis","C. auris","C. glabrata"))


sensitivity <- ggplot(sense_all_spp, aes(x=sense_all_spp$Software, y=sense_all_spp$Sensitivity, fill=sense_all_spp$Software)) +
  geom_boxplot() +  facet_nested(~Species)+theme_bw()+
  ylab("Sensitivity (in %)")+
  xlab("Software")+
  guides(fill=guide_legend(title="Sensitivity"))+
  theme(strip.text = element_text(face = "italic",size = 12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position = "none")#+
#ggsave("~/Dropbox (Personal)/lncRNA/Nat_Com_submission/revision/plots/class_codes.png", units="in", width=10, heigh=6, dpi=600, device = "png")

precision <- ggplot(sense_all_spp, aes(x=sense_all_spp$Software, y=sense_all_spp$Precision, fill=sense_all_spp$Software)) +
  geom_boxplot() +  facet_nested(~Species)+theme_bw()+
  ylab("Specificity (in %)")+
  xlab("Software")+
  guides(fill=guide_legend(title="Precision"))+
  theme(strip.text = element_text(face = "italic",size = 12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position = "none") #+
#ggsave("~/Dropbox (Personal)/lncRNA/Nat_Com_submission/revision/plots/class_codes.png", units="in", width=10, heigh=6, dpi=600, device = "png")
(sensitivity/precision)+plot_annotation(tag_levels = 'a',tag_suffix = "")+ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/sensitivity_specificity.png", units="in", width=8, heigh=6, dpi=600, device = "png")



### Unique for stringite and trinity

uniq_transcripts <- read.table("../mapping_assembly_lncRNA_prediction/trinity/unique_transcirpts_trinity_vs_stringtie.txt")
uniq_transcripts$prcnt <- uniq_transcripts$V1/uniq_transcripts$V2*100
uniq_transcripts$Species <- ifelse(uniq_transcripts$V4=="calb","C. albicans",ifelse(uniq_transcripts$V4=="ctrop","C. tropicalis",ifelse(uniq_transcripts$V4=="cpar","C. parapsilosis",ifelse(uniq_transcripts$V4=="caur","C. auris","C. glabrata"))))
colnames(uniq_transcripts) <- c("Unique","Total","Software","spp","Percentage","Species")
uniq_transcripts$Species <- factor(uniq_transcripts$Species, levels = c("C. albicans","C. tropicalis","C. parapsilosis","C. auris","C. glabrata"))


total <- ggplot(uniq_transcripts, aes(x=uniq_transcripts$Software, y=uniq_transcripts$Total, fill=uniq_transcripts$Software)) +
  geom_boxplot() +  facet_nested(~Species)+theme_bw()+
  ylab("Total number of transcipts")+
  xlab("Software")+
  guides(fill=guide_legend(title="Total numer of transcirpts"))+
  theme(strip.text = element_text(face = "italic",size = 12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position = "none")

unique <- ggplot(uniq_transcripts, aes(x=uniq_transcripts$Software, y=uniq_transcripts$Unique, fill=uniq_transcripts$Software)) +
  geom_boxplot() +  facet_nested(~Species)+theme_bw()+
  ylab("Number of unique transcripts")+
  xlab("Software")+
  #guides(fill=guide_legend(title="Percent (%) of unique transcripts"))+
  theme(strip.text = element_text(face = "italic",size = 12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position = "none")+coord_cartesian(ylim = c(0, 5000))

unique_percent <- ggplot(uniq_transcripts, aes(x=uniq_transcripts$Software, y=uniq_transcripts$Percentage, fill=uniq_transcripts$Software)) +
  geom_boxplot() +  facet_nested(~Species)+theme_bw()+
  ylab("Percent (%) of unique transcripts")+
  xlab("Software")+
  #guides(fill=guide_legend(title="Percent (%) of unique transcripts"))+
  theme(strip.text = element_text(face = "italic",size = 12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position = "none")+coord_cartesian(ylim = c(0, 100))+
  stat_summary(fun.y=mean, colour="blue", geom="point", size=2,show_guide = FALSE)+
  stat_summary(fun.y=mean, colour="black", geom="text", show_guide = FALSE,
               vjust=-0.7, aes( label=round(..y.., digits=1)))


(total+unique+class_codes+unique_percent)+plot_annotation(tag_levels = 'a',tag_suffix = "")+plot_layout(ncol = 1,heights = c(1,1,2,1))+ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/total_unique_classcode_all_features.png", units="in", width=8, heigh=11, dpi=600, device = "png")


###### u and x data

uniq_transcripts_ux <- read.table("../mapping_assembly_lncRNA_prediction/trinity/unique_transcirpts_trinity_vs_stringtie_ux.txt")
uniq_transcripts_ux$prcnt <- uniq_transcripts_ux$V1/uniq_transcripts_ux$V2*100
uniq_transcripts_ux$Species <- ifelse(uniq_transcripts_ux$V4=="calb","C. albicans",ifelse(uniq_transcripts_ux$V4=="ctrop","C. tropicalis",ifelse(uniq_transcripts_ux$V4=="cpar","C. parapsilosis",ifelse(uniq_transcripts_ux$V4=="caur","C. auris","C. glabrata"))))
colnames(uniq_transcripts_ux) <- c("Unique","Total","Software","spp","Class_code","Percentage","Species")

uniq_transcripts_ux$Species <- factor(uniq_transcripts_ux$Species, levels = c("C. albicans","C. tropicalis","C. parapsilosis","C. auris","C. glabrata"))
uniq_transcripts_ux$Type <- ifelse(uniq_transcripts_ux$Class_code=="u","Intergenic","Antisense")
uniq_transcripts_ux$Type <- factor(uniq_transcripts_ux$Type, levels = c("Intergenic","Antisense"))

total_ux <- ggplot(uniq_transcripts_ux, aes(x=uniq_transcripts_ux$Software, y=uniq_transcripts_ux$Total, fill=uniq_transcripts_ux$Software)) +
  geom_boxplot() +  facet_nested(Type~Species,scales = "free")+theme_bw()+
  ylab("Total number of transcirpts")+
  xlab("Software")+
  guides(fill=guide_legend(title="Total numer of transcirpts"))+
  theme(strip.text.x = element_text(face = "italic",size = 12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position = "none")

unique_ux <- ggplot(uniq_transcripts_ux, aes(x=uniq_transcripts_ux$Software, y=uniq_transcripts_ux$Unique, fill=uniq_transcripts_ux$Software)) +
  geom_boxplot() +  facet_nested(Type~Species, scales = "free")+theme_bw()+
  ylab("Number of unique transcripts")+
  xlab("Software")+
  
  theme(strip.text.x = element_text(face = "italic",size = 12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position = "none")

unique_ux_percent <- ggplot(uniq_transcripts_ux, aes(x=uniq_transcripts_ux$Software, y=uniq_transcripts_ux$Percentage, fill=uniq_transcripts_ux$Software)) +
  geom_boxplot() +  facet_nested(Type~Species, scales = "free")+theme_bw()+
  ylab("Percent (%s) of unique transcripts")+
  xlab("Software")+
  theme(strip.text.x = element_text(face = "italic",size = 12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.position = "none")+  stat_summary(fun.y=mean, colour="blue", geom="point", size=2,show_guide = FALSE)+
  stat_summary(fun.y=mean, colour="black", geom="text", show_guide = FALSE,
               vjust=-0.7, aes( label=round(..y.., digits=1)))



(total_ux/unique_ux/unique_ux_percent)+plot_annotation(tag_levels = 'a',tag_suffix = "")+ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/total_vs_unique_ux.png", units="in", width=8, heigh=11, dpi=600, device = "png")



                                                  ##############################
                                                  #### Plot repeat matching ####
                                                  ##############################


repeat_data <- read.table("../repeat_calling/repeat_overlap_stats.txt", header = T, sep = "\t")
repeat_data$Percent <- (repeat_data$N/repeat_data$Total)*100
repeat_data$species <- ifelse(repeat_data$Species=="calb","C. albicans",ifelse(repeat_data$Species=="ctrop","C. tropicalis",ifelse(repeat_data$Species=="cpar","C. parapsilosis",ifelse(repeat_data$Species=="caur","C. auris","C. glabrata"))))
repeat_data$species <- factor(repeat_data$species, levels = c("C. albicans","C. tropicalis","C. parapsilosis","C. auris","C. glabrata"))
repeat_data$lncRNA_type <- ifelse(repeat_data$Class_code=="u","Intergenic",ifelse(repeat_data$Class_code=="x","Antisense","Protein coding"))
repeat_data$lncRNA_type <- factor(repeat_data$lncRNA_type, levels = c("Intergenic","Antisense","Protein coding"))
repeat_data$Type <- factor(repeat_data$Type, levels = c("Simple repeat","Low complexity repeat","LINE/SINE","LTR","Unknown","No overlap"))


ggplot(repeat_data, aes(x=repeat_data$Type, y=repeat_data$Percent, fill=repeat_data$Type)) +
  geom_bar(stat = "identity") +  facet_nested(lncRNA_type~species)+theme_bw()+ylim(0,100)+
  ylab("Percent (%) of lncRNAs overalaping repeats")+geom_text(aes(label=paste0(round(repeat_data$Percent,1),"\n","(",repeat_data$N,")")), position=position_dodge(width=0.9), vjust=-0.25, size = 3)+
  xlab("Types of repeats")+
  theme(strip.text.x = element_text(face = "italic",size = 12),
        strip.text.y = element_text(size = 12),
        axis.text = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  guides(fill=guide_legend(title=""))+scale_fill_brewer(palette="Dark2")+
  theme(legend.position="bottom")+
  ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/repeats.png", units="in", width=11, heigh=9, dpi=600)
  ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/repeats.pdf", units="in", width=11, heigh=9, dpi=600)



                    ##########################################################################
                    ### Plotting Expression of inf-specific co-expressed genes and lncRNAs ###
                    ##########################################################################

make_line_plot_ctrop<-function(id){
coral1_edges <- read.csv("../ploting_and_network_analysis/ctrop/coexpression/edge_coral1_ctrop_new.txt",sep = "\t", stringsAsFactors = F)
inf_edges<-coral1_edges[(coral1_edges$fromNode %in% ctrop_inf_spec_genes) | (coral1_edges$toNode %in% ctrop_inf_spec_genes),]
inf_edges_filtered<-inf_edges[as.numeric(inf_edges$weight)>0.1,]
a<-inf_edges_filtered[order(inf_edges_filtered$weight, decreasing =T ),]

b<-coral1_edges[order(as.numeric(coral1_edges$weight), decreasing = F),]


class(coral1_edges$weight)
#colnames(ctrop_TPM)
to_plot <- t(ctrop_TPM[unlist(a[id,c(1,2)]),])
to_plot_long<-melt(to_plot)
to_plot_long<-cbind(to_plot_long,data.frame(do.call('rbind', strsplit(as.character(to_plot_long$Var2),'|',fixed=TRUE)))[,1])
colnames(to_plot_long)<-c("Sample ID","Transcripts", "log2(TPM)","Transcript ID")

spearman <- cor(to_plot_long[grepl("MSTRG",to_plot_long$`Transcript ID`),][,3],
                to_plot_long[!grepl("MSTRG",to_plot_long$`Transcript ID`),][,3], method = "spearman")

c<-ctrop_inf_spec_protcod <-  read.table("../ploting_and_network_analysis/ctrop/coexpression/ctrop_inf_spec_protcod.txt", stringsAsFactors = F)

id <-ggplot(to_plot_long, 
            aes(x=`Sample ID`,y=log2(`log2(TPM)`),group=`Transcript ID`, colour=`Transcript ID`))+
  geom_line(size=1.5)+
  geom_point(size=2)+
  theme_bw()+
  ylab("log2(TPM)")+
  theme(axis.text.x = element_text(angle = 90))+
  scale_color_manual(values = c("mediumorchid3","grey"))+
  annotate("text",label=paste0("Spearman's rho = ",round(spearman,2)),
           x=10,y=max(log2(to_plot_long$`log2(TPM)`)))
  #ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/ctrop_coexp.pdf", units="in", width=10, heigh=4, dpi=600)+
 # ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/ctrop_coexp.png", units="in", width=10, heigh=4, dpi=600)
return(id)
}

plot_1<-make_line_plot_ctrop(1)
plot_2<-make_line_plot_ctrop(2)
plot_3<-make_line_plot_ctrop(12)

### calb
  edges <- read.csv("../ploting_and_network_analysis/calb/coexpression/edges_darkred_calb.txt",sep = "\t", stringsAsFactors = F)
  inf_edges<-edges[(edges$fromNode %in% calb_inf_spec_genes) | (edges$toNode %in% calb_inf_spec_genes),]
  inf_edges_filtered<-inf_edges[as.numeric(inf_edges$weight)>0.0001,]
  a_calb<-inf_edges_filtered[order(inf_edges_filtered$weight, decreasing =T ),]
  
  b_calb<-edges[order(as.numeric(edges$weight), decreasing = F),]
  
  
  #class(coral1_edges$weight)
  
  to_plot <- t(calb_TPM[unlist(a_calb[1,c(1,2)]),])
  to_plot_long<-melt(to_plot)
  to_plot_long<-cbind(to_plot_long,data.frame(do.call('rbind', strsplit(as.character(to_plot_long$Var2),'|',fixed=TRUE)))[,1])
  colnames(to_plot_long)<-c("Sample ID","Transcripts", "log2(TPM)","Transcript ID")
  
  
  spearman <- cor(to_plot_long[grepl("MSTRG",to_plot_long$`Transcript ID`),][,3],
      to_plot_long[!grepl("MSTRG",to_plot_long$`Transcript ID`),][,3], method = "spearman")
  
  
  
  to_plot_long<-to_plot_long[(grepl("concat",to_plot_long$`Sample ID`)) |(grepl("calb",to_plot_long$`Sample ID`)) ,]
  plot_4<-ggplot(to_plot_long, aes(x=`Sample ID`,y=log2(`log2(TPM)`),group=`Transcript ID`, colour=`Transcript ID`))+
    geom_line(size=1.5)+
    geom_point(size=2)+
    #geom_violin()+
    theme_bw()+
    ylab("log2(TPM)")+
    theme(axis.text.x = element_text(angle = 90))+
    scale_color_manual(values = c("mediumorchid3","grey"))+
    annotate("text",label=paste0("Spearman's rho = ",round(spearman,2)),
            x=10,y=max(log2(to_plot_long$`log2(TPM)`)))#+
  #ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/ctrop_coexp.pdf", units="in", width=10, heigh=4, dpi=600)+
  # ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/ctrop_coexp.png", units="in", width=10, heigh=4, dpi=600)

  ((plot_4/plot_1/plot_1)|(plot_1/plot_2/plot_3))+
    ggsave("~/Dropbox/lncRNA/Nat_Com_submission/revision/plots/coexpression_lines.pdf", units="in", width=16, heigh=14, dpi=600)
  
  
  
### Generate mean stats for length, expression, GC and SNPs ###
length_mean <- NULL
for (spp in c("calb","ctrop","cpar","caur","cglab")){
  length_mean <- rbind.data.frame(length_mean, cbind("length",spp,as.data.frame(t(aggregate( all_data[all_data$Spp == spp,][,6], list(all_data[all_data$Spp == spp,][,9]), mean)))[2,],"NA"))
}
colnames(length_mean)<- c("category","spp","pc","i","a","inter")



expression_mean <- NULL
for (spp in c("calb","ctrop","cpar","caur","cglab")){
  expression_mean <- rbind.data.frame(expression_mean, cbind("expression",spp,as.data.frame(t(aggregate( as.numeric(all_data_expression_long[all_data_expression_long$Spp == spp,][,4]), list(all_data_expression_long[all_data_expression_long$Spp == spp,][,6]), mean)))[2,],"NA"))
}
colnames(expression_mean)<- c("category","spp","pc","i","a","inter")


gc_mean <- NULL
for (spp in c("calb","ctrop","cpar","caur","cglab")){
  gc_mean <- rbind.data.frame(gc_mean, cbind("GC",spp,as.data.frame(t(aggregate( as.numeric(GC_content_all[GC_content_all$spp == spp,][,3]), list(GC_content_all[GC_content_all$spp == spp,][,5]), mean)))[2,]))
}
colnames(gc_mean)<- c("category","spp","pc","i","a","inter")


variants_mean <- NULL
for (spp in c("calb","ctrop","cpar","caur","cglab")){
  variants_mean <- rbind.data.frame(variants_mean, cbind("variants",spp,as.data.frame(t(aggregate(as.numeric(all_snps[all_snps$Spp == spp,][,8]), list(all_snps[all_snps$Spp == spp,][,9]), mean)))[2,]))
}

colnames(variants_mean)<- c("category","spp","pc","i","a","inter")


all_mean_values <- rbind.data.frame(length_mean,expression_mean,gc_mean,variants_mean)
write.table(all_mean_values,"./mean_stats.txt", quote = F,sep = "\t")



                                                ##############################
                                                ####          END         ####
                                                ##############################

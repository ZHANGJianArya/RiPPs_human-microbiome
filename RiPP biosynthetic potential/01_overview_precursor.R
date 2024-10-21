library(ggplot2)
library(reshape2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(ggradar)
library(ggsci)
library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)
library(VennDiagram)
library(vegan)
library(ggthemes)
rm(list = ls())
load("RiPP_precursor_identified_from_human_microbiome.Rdata")
DeepRiPP_stat <- as.data.frame(table(DeepRiPP_1$class))
colnames(DeepRiPP_stat) <- c("class","counts")
TrRiPP_stat <- as.data.frame(table(TrRiPP_1$class))
colnames(TrRiPP_stat) <- c("class","counts")

# Load ggplot2
# Compute the position of labels
##DeepRiPP#####
if(T){
  BB = DeepRiPP_stat
  BB$class <- factor(BB$class, levels = c(
    "rSAM_MODIFIED_RiPP",
    "LANTHIPEPTIDE",
    "AUTO_INDUCING_PEPTIDE",
    "LASSO_PEPTIDE",
    "BACTERIAL_HEAD_TO_TAIL_CYCLIZED",
    "GRASPETIDE",
    "CYANOBACTIN",
    "THIOPEPTIDE",
    "OTHER"
  ))
  BB <- BB %>% 
    arrange(desc(class)) %>%
    mutate(prop = counts / sum(BB$counts) *100)
  # Basic piechart
  pdf(file = "Human microiome_RiPP_DeepRiPP_stat_1015.pdf",width = 8,height =5)
  colourCount = length(unique(BB$class))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  ggplot(BB, aes(x="", y=prop, fill=class)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 10),
          legend.position="right") +
    scale_fill_lancet()
  #scale_fill_manual(values = getPalette(colourCount))
  
  dev.off()  
}
####TrRiPP####
if(T){
  TT = TrRiPP_stat
  TT$class <- factor(TT$class,levels = c(
    "rSAM_MODIFIED_RiPP",
    "LANTHIPEPTIDE",
    "AUTO_INDUCING_PEPTIDE",
    "LASSO_PEPTIDE",
    "BACTERIAL_HEAD_TO_TAIL_CYCLIZED",
    "GRASPETIDE",
    "CYANOBACTIN",
    "THIOPEPTIDE",
    "OTHER"
  ))
  TT <- TT %>% 
    arrange(desc(class)) %>%
    mutate(prop = counts / sum(TT$counts) *100) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  # Basic piechart
  pdf(file = "Human microbiome_RiPP_TrRiPP_stat_150.pdf",width = 8,height =5)
  ggplot(TT, aes(x="", y=prop, fill=class)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() + 
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 10),
          legend.position="right") +
    scale_fill_lancet()
  
  dev.off()  
}

####Venn diagram######

RiPP_Tr_pre <- TrRiPP_1[1]
colnames(RiPP_Tr_pre) <- c("RiPP_pre")
Deep_RiPP_pre <- DeepRiPP_1[1]
colnames(Deep_RiPP_pre) <- c("RiPP_pre")

venn.plot <- venn.diagram(x = list(DeepRiPP=DeepRiPP_1$Meta, TrRiPP = TrRiPP_1$Meta),
                          lwd = 4,
                          filename = NULL,
                          fill = c("#6b72b1", "#f21826"),
                          alpha =0.7,
                          label.col = "white",
                          cex =4 ,
                          at.pos = c(-20, 14)) 
pdf(file = "Venn.pdf",height = 7, width = 7)
par(mar=c(4,4,4,4)+0.1,xpd=TRUE)
grid.draw(venn.plot)
dev.off()

##################################################################
########Supplementary figure 3 The proportion of precursor########
##################################################################
df <- as.data.frame(table(All_pre_taxa_class_1$category))


colnames(df) <-c("Category","counts")
df$precursor <- "precursor"
df$Category <- factor(df$Category,levels=c("DeepRiPP_TrRiPP_antiSMASH","DeepRiPP_TrRiPP","DeepRiPP_antiSMASH","TrRiPP_antiSMASH","DeepRiPP","TrRiPP"))

df <- df %>% 
  arrange(desc(Category)) %>%
  mutate(prop = counts / sum(df$counts) *100)

# Stacked + percent
pdf(file="All_predicted_precursor_category_231125.pdf",width = 6,height = 2.5)
ggplot(df, aes(fill=Category, y=counts,x=precursor)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_npg()+
  theme_bw(base_size = 15)+
  theme(legend.position = "top",legend.title = element_blank())+
  theme(axis.title = element_blank())+
  coord_flip()
dev.off()

########################################
#########taxa########
########################################
aa=All_pre_taxa_class_1[,c(1,3,4,11)]

phylum <- function(x){
  res <- str_split(x,";")[[1]]
  paste(str_replace(res[2],"p__",""))
}
aa$phylum <- unlist(lapply(aa$Lineage,phylum))
pre_p <- as.data.frame(table(aa$phylum))
colnames(pre_p) <- c("phylum","precursor_counts")

gg <- aa[!duplicated(aa$Genome),]
gg_1 <- as.data.frame(table(gg$phylum))
colnames(gg_1) <- c("phylum","genome_counts")

pre_stat <- merge(pre_p,gg_1,by="phylum")
pre_stat$average <- pre_stat$precursor_counts / pre_stat$genome_counts

df_pre <- pre_stat[,c(1,4)]

library(tidyverse)
df_pre <- df_pre %>% 
  add_row(phylum = "Halobacteriota", average=0)

####heatmap  ######
library(dplyr)
library(ggplot2)
df_pre$phylum <- factor(df_pre$phylum,levels =c("Eremiobacterota","Elusimicrobiota","Planctomycetota","Verrucomicrobiota","Fibrobacterota","Halobacteriota","Deinococcota","Desulfobacterota","Desulfobacterota_I","Fusobacteriota","Chlamydiota","Proteobacteria","Cyanobacteria","Thermoplasmatota","Campylobacterota","Spirochaetota","Bacteroidota","Synergistota","Actinobacteriota","Firmicutes_C","Methanobacteriota","Myxococcota","Chloroflexota","Firmicutes_G","Patescibacteria","Bdellovibrionota","Firmicutes","Firmicutes_A","Firmicutes_B"))

df_pre$phylum <- fct_rev(df_pre$phylum)
pdf(file = "All_pre_taxa_stat_231125.pdf",width = 2.2,height =12)

ggplot(df_pre, aes(x=0,fill=average,y=phylum))+ 
  geom_raster()+ 
  geom_text(aes(label = sprintf("%.0f",average)), position = "identity")+
  scale_fill_gradient2(low = "#80a9d7", high = "#f3ad61")+
  theme_bw(base_size = 14)+
  theme(axis.text.x=element_blank())+
  theme(axis.title = element_blank())+
  theme(legend.position = "top",legend.title = element_blank())

dev.off()

###proportion of genome related with RiPP
df_p <- data.frame(genome =c("RiPP_BGC","Other_g","precursor","OtherPre_g"),
                   counts=c(239927,66554,292041,14440))
df_BGC <- df_p[c(1,2), ]
df_pre_g <- df_p[c(3,4), ]
##BGC stat##
# Compute the position of labels
df_BGC <- df_BGC %>% 
  arrange(desc(genome)) %>%
  mutate(prop = counts / sum(df_BGC$counts) *100) 
# Basic piechart

pdf(file = "RiPP_BGC_genome_proportion_stat.pdf",width = 2,height =2)
ggplot(df_BGC, aes(x="", y=prop, fill=genome)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position="right") +
  scale_fill_lancet()

dev.off() 

##Pre stat##
# Compute the position of labels
df_pre_g <- df_pre_g %>% 
  arrange(desc(genome)) %>%
  mutate(prop = counts / sum(df_pre_g$counts) *100) 


pdf(file = "RiPP_pre_genome_proportion_stat.pdf",width = 2,height =2)
ggplot(df_pre_g, aes(x="", y=prop, fill=genome)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position="right") +
  scale_fill_lancet()

dev.off() 

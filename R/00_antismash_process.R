library(dplyr)
library(tidyverse)
library(plyr)
library(data.table)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggsci)
library(openxlsx)
###############################
####Human microbiome BGC ######
###############################
load("../Human_microbiome_BGC_GTDB.Rdata")
##proportion of BGC class######
HM_BGC_count <- as.data.frame(table(HM_BGCs_g$Class_bigscape_rule))
colnames(HM_BGC_count) <- c("group","value")

BB = HM_BGC_count
##BGC stat##
# Compute the position of labels
BB <- BB %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(BB$value) *100)
  #mutate(ypos = cumsum(prop)- 0.5*prop )
# Basic piechart
BB$group <- factor(BB$group,levels = c("NRPS",
                                       "RiPPs",
                                       "Terpene",
                                       "PKS_other",
                                       "PKS/NRPS Hybrids",
                                       "PKS_I",
                                       "Saccharides",
                                       "Others"))
pdf(file = "HM_BGC_stat.pdf",width = 5,height =5)
ggplot(BB, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.title=element_blank(),legend.position="right") +
  scale_fill_lancet()
#geom_text(aes(y = ypos), color = "white", size=6) +
#scale_color_futurama()
dev.off()  

######## Stacked + percent
BB$group <- factor(BB$group,levels = c("RiPPs",
                                       "NRPS",
                                       "Terpene",
                                       "PKS_other",
                                       "PKS/NRPS Hybrids",
                                       "PKS_I",
                                       "Saccharides",
                                       "Others"))
pdf(file="All_predicted_BGC_category_stacked_bar.pdf",width = 6,height = 2.5)
ggplot(BB, aes(fill=group, y=value,x="")) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_npg()+
  theme_bw(base_size = 15)+
  theme(legend.position = "top",legend.title = element_blank())+
  theme(axis.title = element_blank())+
  coord_flip()
dev.off()



######BGC counts per genome#########
taxa <- function(x){
  res <- str_split(x,";")[[1]]
  paste(str_replace(res[2],"p__",""))
}
aa=HM_BGCs_g
aa$phylum <- unlist(lapply(aa$Lineage,taxa))
phylum_BGC_stat <- table(aa$phylum,aa$Class_bigscape_rule)
write.table(phylum_BGC_stat,file="phylum_BGC_stat_GTDB.tsv",sep="\t",quote=F,row.names = T)
map_1 <- read.csv("phylum_BGC_stat_GTDB.tsv",sep="\t")
map_2 <- rownames_to_column(map_1)
##genome per phylum
bb <- aa[!duplicated(aa$Genome), ]
bb_n <- as.data.frame(table(bb$phylum))
colnames(bb_n) <- c("phylum","Genome_counts")

cc <- merge(map_2,bb_n, by.x = "rowname",by.y="phylum")
cc_1 <- tibble::column_to_rownames(cc,var = "rowname")

df2 <- cc_1/cc_1$Genome_counts
df2 <- df2[-9]

###split by kingdom
kingdom <- function(x){
  res <- str_split(x,";")[[1]]
  paste(str_replace(res[1],"d__",""))
}

aa$kingdom <- unlist(lapply(aa$Lineage,kingdom))
aa_1 <- aa[!duplicated(aa$phylum), ]
aa_2 <- aa_1[,c(14,15)]

###heatmap figure S1B #####
if(T){
  library(ComplexHeatmap)
  library(circlize)
  
  map <- as.matrix(df2)
  pdf(file="SMBGCs from Human microbiome_GTDB.pdf",width = 7,height = 12)
  col_fun = colorRamp2(c(0,0.0000000001, max(map)), c("white", "blue", "red"))
  #col_fun = colorRamp2(c(0,0.0000000001,max(map)), c("white","#ccacff", "blue"))
  
    htmp <- Heatmap(map,
                  col=col_fun,
                  row_names_gp = gpar(fontsize = 14),
                  row_names_side = c("left"),
                  column_names_side = c("bottom"),
                  show_row_dend = F,
                  show_column_dend = F,
                  row_names_max_width = unit(6, "cm"),
                  column_names_gp = gpar(fontsize = 14),
                  cluster_rows =T,cluster_columns =T,
                  rect_gp = gpar(col = "black", lwd = 1),
                  #row_order = sort(rownames(map)),
                  heatmap_legend_param=list(title = "BGC counts per genome", legend_direction = c("horizontal"),legend_width = unit(4, "cm")))
  draw(htmp, heatmap_legend_side="top", annotation_legend_side="right",
       legend_grouping = "original")
  
  dev.off()
}

####BGC number in each BGC classes####
BGC_classes <- HM_BGCs_g[ ,c(1,6,7,8)]
BGC_classes$Cluster_number <- as.numeric(BGC_classes$Cluster_number)
BGC_classes[which(BGC_classes$Cluster_number > 1), "Type"]="Combination"
class_type <- as.data.frame(table(BGC_classes$Class_bigscape_rule,BGC_classes$Type))
class_type <- class_type[class_type$Freq >0, ]
colnames(class_type) <- c("BGC_class", "BGC_type","Counts")
write.table(class_type,"class_type_stat_GTDB.tsv",sep="\t",quote=F,row.names = F)

####heatmap figure S1C ######
library(dplyr)
library(ggplot2)
df <- read.delim2("class_type_stat_GTDB.tsv")

df$BGC_type <- factor(df$BGC_type, levels = c(
"bottromycin","RiPPs_Combination","cyanobactin","cyclic-lactone-autoinducer","epipeptide","glycocin","lanthipeptide-class-i","lanthipeptide-class-ii","lanthipeptide-class-iii","lanthipeptide-class-iv","lanthipeptide-class-v","LAP","lassopeptide","linaridin","microviridin","proteusin","ranthipeptide","RaS-RiPP","redox-cofactor","RiPP-like","RRE-containing","sactipeptide","thioamitides","thiopeptide","PKS/NRPS Hybrids","PKS_other_Combination","hglE-KS","PKS-like","T2PKS","T3PKS","transAT-PKS","transAT-PKS-like","T1PKS","acyl_amino_acids","arylpolyene","betalactone","blactam","butyrolactone","CDPS","Others_Combination","ectoine","furan","hserlactone","ladderane","NAGGN","nucleoside","other","phenazine","phosphonate","prodigiosin","resorcinol","siderophore","tropodithietic-acid","NRPS_Combination","NAPAA","NRPS","NRPS-like","thioamide-NRP","terpene","amglyccycl","oligosaccharide"))
df$BGC_type <- fct_rev(df$BGC_type)
pdf(file="BGC_class_type_count.pdf", width =4.5,height = 12)
ggplot(df, aes(x=0,fill=log10(Counts),y=BGC_type))+ 
  geom_raster()+ 
  geom_text(aes(label = Counts), position = "identity")+
  scale_fill_gradient2(low = "#80a9d7", high = "blue")+
  theme_bw(base_size = 14)+
  theme(axis.text.x=element_blank())
dev.off()







genome_BGC_counts <- as.data.frame(table(Bacteria_BGC$Genome))
colnames(genome_BGC_counts) <- c("Genome","BGC_counts")
genome_BGC_counts$BGC_counts <- as.numeric(genome_BGC_counts$BGC_counts)
mean(genome_BGC_counts$BGC_counts)

##Phylogenetic trends ###
uni_genome <- Bacteria_BGC[!duplicated(Bacteria_BGC$Genome), ]
aa <- Bacteria_BGC
aa$Phylum <- str_replace(aa$Phylum,"p__","")
phylum_BGC_stat <- table(aa$Phylum,aa$Class_bigscape_rule)
write.table(phylum_BGC_stat,file="phylum_BGC_stat.tsv",sep="\t",quote=F,row.names = T)

####### supplementary data 1 #######
SD <- HM_BGCs_g[-12]
write.table(SD, file="Supplementary Data 1 Secondary metabolites BGCs detected in human microbiome.tsv",sep="\t",quote=F,row.names = F)

write.xlsx(SD,"Supplementary Data 1 Secondary metabolites BGCs detected in human microbiome.xlsx")

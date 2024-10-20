library(readxl)
library(tidyr)
library(mlr)
library(tidyverse)
library(openxlsx)
library(parallelMap)
library(parallel)
library(ggthemes)
library(ggpubr)
library(ggsci)
library(circlize)
library(webshot)
library(stringr)
library(reshape2)

###################################################################
#################### activity of single RiPP ######################
###################################################################
load("RiPP_family_novelty_bioactivity.Rdata")
df <- DeepBGC_output
df$product_activity <- df$product_activity %>% replace_na("unknown")
df$product_class <- df$product_class %>% replace_na("No_confident_class")

df_pro <- table(df$sequence_id,df$product_class)
write.table(df_pro,file="df_product_class.tsv",sep="\t",quote = F)
df_pro <- read.delim2("df_product_class.tsv")

df_fun <- table(df$sequence_id,df$product_activity)
write.table(df_fun,file="df_function.tsv",sep="\t",quote = F)
df_fun <- read.delim2("df_function.tsv")

# re-define product class ### ##GC: genomic context ###
df_pro$GC_class <- "novel"
df_pro$hybrid <- df_pro$NRP +df_pro$NRP.Polyketide + df_pro$Other + df_pro$Polyketide + df_pro$Polyketide.Terpene +df_pro$Saccharide + df_pro$Saccharide.Terpene + df_pro$Terpene
df_pro[which(df_pro$hybrid >= 1),"GC_class"]= "Hybrids"
df_pro[which(df_pro$RiPP >= 1),"GC_class"]= "RiPP"
table(df_pro$GC_class)

### reclass activity ###
df_fun$activity <- apply(df_fun, 1, function(x) paste(colnames(df_fun)[which(x >= 1)], collapse = "_"))

df_fun_stat <- as.data.frame(table(df_fun$activity))

df_fun$final_activity <- "multiple"
df_fun[which(df_fun$activity=="antibacterial"),"final_activity"]="antibacterial"
df_fun[which(df_fun$activity=="antibacterial_unknown"),"final_activity"]="antibacterial"
df_fun[which(df_fun$activity=="inhibitor"),"final_activity"]="inhibitor"
df_fun[which(df_fun$activity=="inhibitor_unknown"),"final_activity"]="inhibitor"
df_fun[which(df_fun$activity=="cytotoxic"),"final_activity"]="cytotoxic"
df_fun[which(df_fun$activity=="cytotoxic_unknown"),"final_activity"]="cytotoxic"
df_fun[which(df_fun$activity=="unknown"),"final_activity"]="unknown"

df_fun_stat2 <- as.data.frame(table(df_fun$final_activity))

df_fun_2 <- df_fun[9]
df_pro_2 <- df_pro[11]

df_pro_fun <- merge(df_pro_2,df_fun_2,by="row.names")

#################################################################
### Genomic context : antismash + deepbgc (single precursor) ####
#################################################################

pre_anti <- Prioritized_pre_domain_class[,c(1,10)]
antismash_RiPP <- pre_anti[str_detect(pre_anti$source,"antismash"),]

deepbgc_RiPP <- df_pro_fun[df_pro_fun$GC_class=="RiPP", ]
deepbgc_novel <- df_pro_fun[df_pro_fun$GC_class=="novel", ]
deepbgc_Hybrids <- df_pro_fun[df_pro_fun$GC_class=="Hybrids", ]

## + deepbgc
pre_GC_class <- pre_anti
pre_GC_class$Class <- "Others"
pre_GC_class[which(pre_GC_class$RiPP_pre %in% deepbgc_RiPP$Row.names),"Class"]="deepbgc_RiPP"
pre_GC_class[which(pre_GC_class$RiPP_pre %in% deepbgc_novel$Row.names),"Class"]="deepbgc_Novel"
pre_GC_class[which(pre_GC_class$RiPP_pre %in% deepbgc_Hybrids$Row.names),"Class"]="Others"
pre_GC_class[which(pre_GC_class$RiPP_pre %in% antismash_RiPP$RiPP_pre),"Class"]="antismash_RiPP"

pre_GC_class_stat <- as.data.frame(table(pre_GC_class$Class))

pre_GC_class_stat$Var1 <- factor(pre_GC_class_stat$Var1,levels=c("antismash_RiPP","deepbgc_RiPP",
                                                                 "deepbgc_Novel",
                                                                 "Others"))
pre_GC_class_stat$x <- "x"
pdf(file="genomic_context_antismash_deepbgc_20240529.pdf",width=3,height=4)

ggplot(pre_GC_class_stat, aes(fill=Var1, y=Freq, x=x)) + 
  geom_bar(position="fill", stat="identity")+
  theme_bw()+
  scale_fill_lancet()
dev.off()
pre_GC_class_stat$pro <- pre_GC_class_stat$Freq/sum(pre_GC_class_stat$Freq) *100
############################################################
### Genomic context : antismash + deepbgc (RiPP family) ####
############################################################
#####################################################
######### RiPP family: class and function ###########
#####################################################
####### family function ######

pre_family <- Prioritized_pre_domain_class[,c(1,11)]
pre_family$family_id <- str_replace_all(pre_family$family_id,"bacteria_family","BF_")
pre_family$family_id <- str_replace_all(pre_family$family_id,"archaea_family","AF_")

df_pro_fun[which(df_pro_fun$Row.names %in% antismash_RiPP$RiPP_pre),"GC_class"]="RiPP"

#df_pro_fun_1 <- df_pro_fun[!df_pro_fun$final_activity=="unknown", ]

df_fun_all <- merge(pre_family,df_pro_fun,by.x="RiPP_pre",by.y = "Row.names",all = T)

df_fun_all <- df_fun_all[-3]
df_fun_all$final_activity <- df_fun_all$final_activity %>% replace_na("unknown")
table(df_fun_all$final_activity)

##### family function stat #####
fm_fun <- table(df_fun_all$family_id,df_fun_all$final_activity) %>% as.matrix()

df_fm_fun <- data.frame(family_id=rep(rownames(fm_fun), times=ncol(fm_fun)),
                        activity=rep(colnames(fm_fun), each=nrow(fm_fun)),
                        value=as.vector(fm_fun),
                        stringsAsFactors = FALSE)
## long to wide 
df_fm_fun_w <- dcast(df_fm_fun, family_id ~ activity)

x <- df_fm_fun_w
df_fm_fun_prop <- cbind(id = x[, 1], x[, -1]/rowSums(x[, -1])*100)
rownames(df_fm_fun_prop) <- df_fm_fun_prop$id
df_fm_fun_prop <- df_fm_fun_prop[-1]

df_fm_fun_prop$activity <- apply(df_fm_fun_prop, 1, function(x) paste(colnames(df_fm_fun_prop)[which(x >= 75)], collapse = "_"))

df_fm_fun_prop[which(df_fm_fun_prop$activity==""),"activity"]="unknown"

table(df_fm_fun_prop$activity)

df_family_function <- df_fm_fun_prop[6]

###### class ####
GC_class <- pre_anti
GC_class$Class <- "Others"
GC_class[which(GC_class$RiPP_pre %in% deepbgc_RiPP$Row.names),"Class"]="RiPP"
GC_class[which(GC_class$RiPP_pre %in% deepbgc_novel$Row.names),"Class"]="Novel"
GC_class[which(GC_class$RiPP_pre %in% deepbgc_Hybrids$Row.names),"Class"]="Hybrids"
GC_class[which(GC_class$RiPP_pre %in% antismash_RiPP$RiPP_pre),"Class"]="RiPP"

GC_class_stat <- as.data.frame(table(GC_class$Class))

df_gcclass_all <- merge(pre_family,GC_class,by="RiPP_pre")

##### family function stat #####
fm_gcclass <- table(df_gcclass_all$family_id,df_gcclass_all$Class) %>% as.matrix()

df_fm_gcclass <- data.frame(family_id=rep(rownames(fm_gcclass), times=ncol(fm_gcclass)),
                            fm_gcclass=rep(colnames(fm_gcclass), each=nrow(fm_gcclass)),
                            value=as.vector(fm_gcclass),
                            stringsAsFactors = FALSE)

df_fm_gcclass_w <- dcast(df_fm_gcclass, family_id ~ fm_gcclass)

gc <- df_fm_gcclass_w
df_fm_gcclass_prop <- cbind(id = gc[, 1], gc[, -1]/rowSums(gc[, -1])*100)
rownames(df_fm_gcclass_prop) <- df_fm_gcclass_prop$id
df_fm_gcclass_prop <- df_fm_gcclass_prop[-1]

df_fm_gcclass_prop$gc_class <- apply(df_fm_gcclass_prop, 1, function(x) paste(colnames(df_fm_gcclass_prop)[which(x >= 50)], collapse = "_"))

table(df_fm_gcclass_prop$gc_class)

df_fm_gcclass_prop$GC_Novel <- "Others"
df_fm_gcclass_prop[which(df_fm_gcclass_prop$gc_class=="Hybrids"),"GC_Novel"]="Others"
df_fm_gcclass_prop[which(df_fm_gcclass_prop$gc_class=="Hybrids_Others"),"GC_Novel"]="Others"
df_fm_gcclass_prop[which(df_fm_gcclass_prop$gc_class=="Others"),"GC_Novel"]="Others"

df_fm_gcclass_prop[which(df_fm_gcclass_prop$gc_class=="Hybrids_RiPP"),"GC_Novel"]="RiPP"
df_fm_gcclass_prop[which(df_fm_gcclass_prop$gc_class=="Novel_RiPP"),"GC_Novel"]="RiPP"
df_fm_gcclass_prop[which(df_fm_gcclass_prop$gc_class=="Others_RiPP"),"GC_Novel"]="RiPP"
df_fm_gcclass_prop[which(df_fm_gcclass_prop$gc_class=="RiPP"),"GC_Novel"]="RiPP"

df_fm_gcclass_prop[which(df_fm_gcclass_prop$gc_class=="Hybrids_Novel"),"GC_Novel"]="Novel"
df_fm_gcclass_prop[which(df_fm_gcclass_prop$gc_class=="Novel_Others"),"GC_Novel"]="Novel"
df_fm_gcclass_prop[which(df_fm_gcclass_prop$gc_class=="Novel"),"GC_Novel"]="Novel"


table(df_fm_gcclass_prop$GC_Novel)

pre_n <- RiPP_family_domain_class[,c(1,8)]

pre_n$pre_novel <- "Novel"
pre_n[which(!is.na(pre_n$max_domain)),"pre_novel"]="RiPP"
table(pre_n$pre_novel)

pre_n <- pre_n[-2]

pre_GC_novel <- merge(pre_n,df_fm_gcclass_prop,by.x="family_id",by.y="row.names")
pre_GC_novel <- pre_GC_novel[,c(1,2,8)]
pre_GC_novel$merge_novel <- paste(pre_GC_novel$pre_novel,pre_GC_novel$GC_Novel,sep="_")
table(pre_GC_novel$merge_novel)

pre_GC_novel$final_type <- "Others"
pre_GC_novel[which(pre_GC_novel$merge_novel=="RiPP_RiPP"),"final_type"]="classic"
pre_GC_novel[which(pre_GC_novel$merge_novel=="RiPP_Hybrids"),"final_type"]="Novel"
pre_GC_novel[which(pre_GC_novel$merge_novel=="Novel_RiPP"),"final_type"]="Novel"
pre_GC_novel[which(pre_GC_novel$merge_novel=="RiPP_Novel"),"final_type"]="Novel"
pre_GC_novel[which(pre_GC_novel$merge_novel=="Novel_Novel"),"final_type"]="Novel"
pre_GC_novel[which(pre_GC_novel$merge_novel=="Novel_Hybrids"),"final_type"]="Novel"

table(pre_GC_novel$final_type)

df_family_type <- pre_GC_novel[,c(1,5)]

df_family_type_function <- merge(df_family_type,df_family_function,by.x="family_id",by.y="row.names")

df_family_info <- merge(RiPP_family_domain_class,df_family_type_function,by="family_id")
write.xlsx(df_family_info,"Supplementary Data 5 Overview of RiPP precursor families_20241020.xlsx")

##### family_class_function stat highlight function#####
fm_type_class <- table(df_family_type_function$activity,df_family_type_function$final_type) %>% as.matrix()

df_fm_cor <- data.frame(from=rep(rownames(fm_type_class), times=ncol(fm_type_class)),
                        to=rep(colnames(fm_type_class), each=nrow(fm_type_class)),
                        value=as.vector(fm_type_class),
                        stringsAsFactors = FALSE)

df_fm_cor <- data.frame(from=rep(rownames(fm_type_class), times=ncol(fm_type_class)),
                        to=rep(colnames(fm_type_class), each=nrow(fm_type_class)),
                        value=as.vector(fm_type_class),
                        stringsAsFactors = FALSE)

orders <- c("Others","Novel","classic","inhibitor","cytotoxic","antibacterial","multiple","unknown")

circos.clear()
grid.col = c(Others = "#1e4586", Novel = "#5d3bbb", classic = "#64b350", 
             multiple = "#f2b296", inhibitor ="#340000",
             cytotoxic = "#005ba4", unknown = "grey",
             antibacterial = "#d93621")
grid.col2 = c(classic = "#050505", Novel = "#f3323c", Others = "#4967ab", 
              multiple = "#f2b296", inhibitor ="#340000",
              cytotoxic = "#005ba4", unknown = "grey",
              antibacterial = "#5d3bbb")

## adjust gap between elements
pdf(file="final_type_function_chord_diagram_final_240529.pdf",width=6,height=6)
chordDiagram(df_fm_cor, 
             order = orders, 
             grid.col = grid.col, transparency = 0.2,
             scale = TRUE, 
             annotationTrackHeight = c(0.05, 0.05))

dev.off()



df_1 <- as.data.frame(table(df_family_type_function$final_type))
df_1$pro <- df_1$Freq/sum(df_1$Freq) *100

df_2 <- as.data.frame(table(df_family_type_function$activity))
df_2$pro <- df_2$Freq/sum(df_2$Freq) *100

###### family novelty (Figure SI) #####
fm_fun_stat <- as.data.frame(table(pre_GC_novel$merge_novel))

fm_fun_stat$Var1 <- factor(fm_fun_stat$Var1,
                           levels=c("Novel_Others","RiPP_Others",
                                    "Novel_Novel","RiPP_Novel","Novel_RiPP",
                                    "RiPP_RiPP"))

fm_fun_stat$Var1 <- fct_rev(fm_fun_stat$Var1)
fm_fun_stat$x <- "x"

pdf(file="novelty of RiPPs_20240529_abs.pdf",width=3,height=4)

ggplot(fm_fun_stat, aes(fill=Var1, y=Freq, x=x)) + 
  geom_bar(stat="identity")+
  theme_bw()+
  theme(legend.position = "right")+
  scale_fill_lancet()
dev.off()

fm_fun_stat$pro <- fm_fun_stat$Freq/sum(fm_fun_stat$Freq) *100





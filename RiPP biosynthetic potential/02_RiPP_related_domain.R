library(dplyr)
library(tidyverse)
library(plyr)
library(data.table)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggsci)
library(stringi)
library(strex)
library(VennDiagram)
library(vegan)

########################################################################
##################### Biosynthetic domain Figure 1C#####################
########################################################################
rm(list=ls())
load("RiPP_related_domain_all.Rdata")
##RiPP related domain##
Prioritized_pre_PTM_1 <- Prioritized_PTM_all[,c(14,16)]
colnames(Prioritized_pre_PTM_1) <- c("RiPP_pre","RiPP_related_domain")
colnames(RiPP_precursor_domain) <- c("RiPP_pre","RiPP_related_domain")
Pre_RiPP_related_domain <- rbind(Prioritized_pre_PTM_1,RiPP_precursor_domain)

df <- All_pre_taxa_class[,c(1,10,14)]

df$Class <- "Others"
df[which(df$RiPP_pre %in% BE3141_pre$precursor),"Class"] = "OtherBE_domain"
df[which(df$RiPP_pre %in% Other_BGC_pre_Overlap$RiPP_pre),"Class"] = "Others"
df[which(df$RiPP_pre %in% Pre_RiPP_related_domain$RiPP_pre),"Class"] = "RiPP_related"
df[which(str_detect(df$category,"antiSMASH")),"Class"] = "RiPP_related"

table(df$Class)

# plot chord diagram
df_chord <- table(df$Class, df$prediction) %>% as.matrix()
df_cor <- data.frame(from=rep(rownames(df_chord), times=ncol(df_chord)),
                     to=rep(colnames(df_chord), each=nrow(df_chord)),
                     value=as.vector(df_chord),
                     stringsAsFactors = FALSE)
# scale the BGC counts with log10 for better vasulization
df_cor_log <- log10(df_chord+1)
orders <- c("DeepRiPP_TrRiPP","TrRiPP","DeepRiPP",
            "Others","OtherBE_domain","RiPP_related")
circos.clear()
grid.col = c(DeepRiPP_TrRiPP = "#d93621", TrRiPP = "#1e4586", DeepRiPP = "#64b350", 
             Others = "grey", OtherBE_domain ="#340000",
             RiPP_related_domain = "#5d3bbb")#Other_BGC = "#f2b296",#005ba4
#RiPP_BGC = "#5d3bbb")
grid.col2 = c(DeepRiPP_TrRiPP = "#050505", TrRiPP = "#f3323c", DeepRiPP = "#4967ab", 
              Others = "grey", OtherBE_domain ="#340000",
              RiPP_related_domain = "#5d3bbb") #Other_BGC = "#f2b296",
#RiPP_BGC = "#5d3bbb")

## adjust gap between elements
pdf(file="chord_diagram_final_precursor_genomic_context_domain_240529.pdf",width=6,height=6)

circos.par(gap.after=c(rep(2,ncol(df_cor_log)-1),12,rep(2,nrow(df_cor_log)-1),10),start.degree=95, clock.wise=FALSE)
chordDiagram(t(df_chord), order = orders,grid.col = grid.col, transparency = 0.2, scale = TRUE)

dev.off()

df_1 <- as.data.frame(table(df$prediction))
df_1$pro <- df_1$Freq/sum(df_1$Freq) *100

df_2 <- as.data.frame(table(df$Class))
df_2$pro <- df_2$Freq/sum(df_2$Freq) *100




#########################################################################
####### RiPP precursor domain_Supplementary Figure 3 ####################
#########################################################################
###venn
aa <- Prioritized_pre_domain
table(aa$prediction)
aa_DeepRiPP <- aa[str_detect(aa$prediction,"DeepRiPP"),]
aa_TrRiPP <- aa[str_detect(aa$prediction,"TrRiPP"),]

venn.plot <- venn.diagram(x = list(DeepRiPP=aa_DeepRiPP$RiPP_pre, TrRiPP = aa_TrRiPP$RiPP_pre),
                          lwd = 4,
                          filename = NULL,
                          fill = c("#90c98a", "#e99376"),
                          alpha =0.7,
                          label.col = "white",
                          cex =4 ,
                          at.pos = c(-20, 14)) 
pdf(file = "Venn_domain_2tool_predicted.pdf",height = 7, width = 7)
par(mar=c(4,4,4,4)+0.1,xpd=TRUE)
grid.draw(venn.plot)
dev.off()
###bar plot###

TrRiPP_stat <- as.data.frame(table(aa_TrRiPP$precursor_related_domain))
colnames(TrRiPP_stat) <- c("Class","Counts")
TrRiPP_stat$tool <- "TrRiPP"

DeepRiPP_stat <- as.data.frame(table(aa_DeepRiPP$precursor_related_domain))
colnames(DeepRiPP_stat) <- c("Class","Counts")
DeepRiPP_stat$tool <- "DeepRiPP"

#merge all domain containing pre#
stat_1 <- rbind(TrRiPP_stat,DeepRiPP_stat)

library(ggplot2)
library(ggsci)
if(T){
  pdf(file = "RiPP_pre_domain_bytool_1.pdf")
  ggplot(data=stat_1, aes(x=Class, y=Counts, fill=tool)) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=Counts), vjust=1, color="black",
              position = position_dodge(1), size=3.5)+
    scale_fill_d3()+
    #scale_fill_brewer(palette="Paired")+
    theme(axis.title=element_text(size = 16),axis.text=element_text(size = 16))+
    theme(legend.title=element_blank())+
    theme_bw(base_size = 14)+
    ylab("Counts") +
    xlab("RiPP precursor domain") +
    coord_flip()
  dev.off()
}

#########################################################################
####### RiPP biosynthetic enzyme domain_Supplementary Figure 3 ##########
#########################################################################
###PTM
KRiPP_stat <-table(Prioritized_PTM_all$precursor,Prioritized_PTM_all$responsible_for_RiPP_class)
###TrRiPP
write.table(KRiPP_stat, file ="Known_RiPP_intersection_PTM80.tsv",sep ="\t", row.names=T,quote = F)
KRiPP_PTM <- read.delim('Known_RiPP_intersection_PTM80.tsv', sep = '\t', check.names = FALSE)
KRiPP_PTM$PTM_no =  apply(KRiPP_PTM, 1, function(i) sum(i > 0))
table(KRiPP_PTM$PTM_no)
KRiPP_PTM_SF_1 <- KRiPP_PTM[KRiPP_PTM$PTM_no ==1, ]
KRiPP_PTM_SF <- KRiPP_PTM_SF_1[,c(1:12)]
KRiPP_PTM_SF$PTM <- apply(KRiPP_PTM_SF, 1, function(x) paste0(names(KRiPP_PTM_SF)[x == max(x)],collapse = '_'))

KRiPP_PTM_PMCS <- KRiPP_PTM[KRiPP_PTM$PTM_no >1, ]
KRiPP_PTM_PMCS <- KRiPP_PTM_PMCS[ ,c(1:12)]
KRiPP_PTM_PMCS$PTM <- "hybrid_PTM"
KRiPP_PTM_assign <- rbind(KRiPP_PTM_SF,KRiPP_PTM_PMCS)
table(KRiPP_PTM_assign$PTM)


##statistic###
DeepRiPP_P_PTM <- KRiPP_PTM_assign[rownames(KRiPP_PTM_assign) %in% DeepRiPP_pre$RiPP_pre, ]

DeepRiPP_stat <- as.data.frame(table(DeepRiPP_P_PTM$PTM))
colnames(DeepRiPP_stat) <- c("Class","Counts")
DeepRiPP_stat$tool <- "DeepRiPP"


TrRiPP_P_PTM <- KRiPP_PTM_assign[rownames(KRiPP_PTM_assign) %in% TrRiPP_pre$RiPP_pre, ]

TrRiPP_stat <- as.data.frame(table(TrRiPP_P_PTM$PTM))
colnames(TrRiPP_stat) <- c("Class","Counts")
TrRiPP_stat$tool <- "TrRiPP"

overlap <- TrRiPP_P_PTM[rownames(TrRiPP_P_PTM) %in% rownames(DeepRiPP_P_PTM), ]
##plot
stat_1 <- rbind(TrRiPP_stat,DeepRiPP_stat)

library(ggplot2)
library(ggsci)
pdf(file = "Prioritized_RiPP_PTM_domain_bytool_1.pdf")
ggplot(data=stat_1, aes(x=Class, y=Counts, fill=tool)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=Counts), vjust=1, color="black",
            position = position_dodge(1), size=3.5)+
  scale_fill_d3()+
  #scale_fill_brewer(palette="Paired")+
  theme(axis.title=element_text(size = 16),axis.text=element_text(size = 16))+
  theme(legend.title=element_blank())+
  theme_bw(base_size = 14)+
  ylab("Enzyme-precursor co-occurrence frequency") +
  xlab("RiPP-related PTM enzyme domain") +
  coord_flip()
dev.off()





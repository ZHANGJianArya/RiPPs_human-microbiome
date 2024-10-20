library(ggplot2)
library(reshape2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(ggradar)
library(ggsci)
rm(list = ls())
#####abundance#####
#abd1 = read.table("../oral/PRJNA396840/deseq2/diff_expr_results_caries.csv",header = T,sep = ',')
#abd2 = read.table("../PRJEB21446_BV/deseq2/diff_expr_results_BV.csv",header = T,sep = ',')
abd3 = read.table("../PRJEB53891_CRC_MT/Human_microbiome/deseq2/deseq2/diff_expr_results_CRC.csv",header = T,sep = ',')
abd4 = read.table("../PRJNA289586_T1D/deseq2/diff_expr_results_T1D.csv",header = T,sep = ',')
abd5 = read.table("../PRJNA727609_LC/Human_microbiome/deseq2/deseq2/diff_expr_results_LC_HCV.csv",header = T,sep = ',')
abd6 = read.table("../IBD/deseq2/diff_expr_results_CD.csv",header = T,sep = ',')
abd7 = read.table("../IBD/deseq2/diff_expr_results_UC.csv",header = T,sep = ',')

#abd1$disease = "caries"
#abd2$disease = "BV"
abd3$disease = "CRC"
abd4$disease = "T1D"
abd5$disease = "LC_HCV"
abd6$disease = "CD"
abd7$disease = "UC"

mydata_abd = rbind(abd3,abd4,abd5,abd6,abd7)
colnames(mydata_abd)[1] <- c("family")

###########################################################
######################## figure 5A&B ######################
###########################################################

###down in control###
mydata_abd_down <- mydata_abd[mydata_abd$group=="DOWN", ]
mydata_abd_up <- mydata_abd[mydata_abd$group=="UP", ]

down <- as.data.frame(table(mydata_abd_down$disease))
down$class <- "down"
up <- as.data.frame(table(mydata_abd_up$disease))
up$class <- "up"

differ <- rbind(down,up)

differ$Var1 <- as.factor(differ$Var1)
table(differ$Var1)
##plot figure 2c

pdf("fig6b_all_differential_families_gut_0530.pdf")
differ$Var1 <- factor(differ$Var1,levels = c("CRC","T1D","CD","UC","LC_HCV"))
ggplot(differ,aes(x=Var1,y=Freq,group=class,fill=class)) +
  geom_bar(stat = 'identity')+
  geom_text(aes(label=Freq), vjust= -0.2, color="black",
            position = position_dodge(1), size=3.5)+
  facet_grid(class~.,scales = "free")+
  scale_fill_manual(values =  c("#d9544d","#004988"), label = c("control","case"))+
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y = element_blank())
dev.off()

############Upset plot######################
###higher in control###
b=mydata_abd_down[ ,c(1,3,9)]
b[which(b$log2FoldChange <= -1),"class1"] = 1
data.plot = dcast(b[,c(1,3,4)],disease~family,mean,fill = 0)

rownames(data.plot) = data.plot[,1]
data.plot = data.plot[,-1]

##ordered by type
data1 = data.frame(t(data.plot))
data5=data1

##UpSet plot
uu_f <- data5[, c(1:5)]
m = make_comb_mat(uu_f)
pdf(file="Multidisease_upset_both_higher_in control_gut.pdf",width = 6,height = 3)
#UpSet(m)
UpSet(m, 
      pt_size = unit(3, "mm"), 
      lwd = 2,
      top_annotation = upset_top_annotation(m, add_numbers = TRUE),
      right_annotation = upset_right_annotation(m, add_numbers = TRUE))
dev.off()

####lower in control######
l=mydata_abd_up[ ,c(1,3,9)]
l[which(l$log2FoldChange >= 1),"class1"] = 1
data.plot_l = dcast(l[,c(1,3,4)],disease~family,mean,fill = 0)

rownames(data.plot_l) = data.plot_l[,1]
data.plot_l = data.plot_l[,-1]

##ordered by type
data1_l = data.frame(t(data.plot_l))
data5_l=data1_l

##UpSet plot
uu_f_l <- data5_l[ ,c(1:5)]
m_l = make_comb_mat(uu_f_l)
pdf(file="Multidisease_upset_both_higher in disease_gut.pdf",width = 5,height = 3)
#UpSet(m)
UpSet(m_l, 
      pt_size = unit(3, "mm"), 
      lwd = 2,
      top_annotation = upset_top_annotation(m_l, add_numbers = TRUE),
      right_annotation = upset_right_annotation(m_l, add_numbers = TRUE))
dev.off()


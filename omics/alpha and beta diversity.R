#################################################################
#################### Alpha-diversity ####################
#################################################################
# Load the required packages
library(vegan)
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(patchwork)
library(agricolae)
library(FSA)
library(rcompanion)
library(ggsci)
rm(list = ls())
#prevalence
load("countdata_metadata_IBD.Rdata")

countdata_3 <- countdata_3[ ,colSums(countdata_3)>0]
otu_raw <- countdata_3
data_otu <- t(otu_raw)

data_grp_1 <- data_grp[match(rownames(data_otu), data_grp$id), ]
data_shannon <- diversity(data_otu, index = "shannon")      # calculate Shannon index using vegan package
data_alphadiv <- cbind(data_grp_1, data_shannon) # combine all indices in one data table

my_comparison <- list(c("nonIBD","CD"),c("nonIBD","UC"))
data_alphadiv$condition <- factor(data_alphadiv$condition,levels = c("CD","nonIBD","UC"))
##shannon
pdf(file="Alpha-diversty_shannon.pdf",width = 3,height = 3)
IBD_MT <- ggplot(data_alphadiv, aes(x=condition, y=data_shannon)) +
  geom_boxplot(fill=c("green","blue","red")) +
  stat_compare_means(comparisons = my_comparison, method="wilcox.test")+
  labs(x= ' ', y= 'Shannon') +
  theme_bw(base_size = 12)+
  theme(legend.position = "none")+
  scale_fill_lancet()
IBD_MT
dev.off()

#################################################################
#################### Beta-diversity ####################
#################################################################
library(openxlsx)
library(ade4)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(vegan)
library(ggsci)
library(ggpubr)
library(reshape2)
library(dplyr)
library(tidyr)

load("betadiversity_source_data.Rdata")
countdata_2 <- countdata_1[ ,colSums(countdata_1) != 0]
countdata_3 <- countdata_2[rowSums(countdata_2>0)>ncol(countdata_2)*0.05, ]

countdata <- as.data.frame(countdata_3)

ct <- as.data.frame(t(countdata))
library(tibble)
ct_1 <- tibble::rownames_to_column(ct,var="id")
ct_2 <-select(ct_1,c(1))
#group_1<-metadata
group <- merge(group,ct_2,by="id")
#write.table(group_1, file="mymeta_group.csv",sep = "\t",quote=F,row.names = F)
#save(countdata_1,group,file="betadiversity_source_data.Rdata")

data <- countdata[,as.vector(group$id)]
data <- data[,colSums(data)>0]
distance <- vegdist(t(data), method = 'bray')
distance

pcoa <- cmdscale(distance, k = (nrow(t(data)) - 1), eig = TRUE)

# plot using 'ordiplot', a in-house function in vegan
ordiplot(scores(pcoa)[ ,c(1, 2)], type = 't')
# check info
summary(pcoa)
# check the coordinate value
pcoa$eig
point <- data.frame(pcoa$point)
# save the coordinate info
#write.csv(point, 'pcoa.sample.csv')

# compute first two coordinate value
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)

# extract first two coordinate value
sample_condition <- data.frame({pcoa$point})[1:2]
sample_condition$id <- rownames(sample_condition)
names(sample_condition)[1:2] <- c('PCoA1', 'PCoA2')

sample_condition <- merge(sample_condition, group, by = 'id', all.x = FALSE)

#write.csv(sample_condition, 'sample_condition.csv', quote = F)
### PERMANOVA for conditions
dt <- data.frame(t(data))
#write.table(data,file="data_f.tsv",quote = F)

group_1 <- read.table('mymeta_group.csv',header=T,sep="\t")
group_1 <- group_1[group_1$id %in% sample_condition$id, ]
adonis_result  <- adonis(dt ~ condition, group_1, distance = 'bray', permutations = 999)
#adonis_result[is.na(adonis_result)] <- 0
# save result
otuput <- data.frame(adonis_result$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
#write.table(otuput, file = 'PERMANOVA.result_all_condition.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')
# pairwise
group_name <- unique(group_1$condition)

adonis_result_two <- NULL
for (i in 1:(length(group_name) - 1)) {
  for (j in (i + 1):length(group_name)) {
    group_ij <- subset(group_1, condition %in% c(group_name[i], group_name[j]))
    otu_ij <- dt[group_ij$id, ]
    adonis_result_otu_ij <- adonis(otu_ij~condition, group_ij, permutations = 999, distance = 'bray')   
    adonis_result_two <- rbind(adonis_result_two, c(paste(group_name[i], group_name[j], sep = '/'), 'Bray-Curtis', unlist(data.frame(adonis_result_otu_ij$aov.tab, check.names = FALSE)[1, ])))
  }
}
adonis_result_two <- data.frame(adonis_result_two, stringsAsFactors = FALSE)
names(adonis_result_two) <- c('group', 'distance', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')

# marker significance
for (i in 1:nrow(adonis_result_two)) {
  if (adonis_result_two[i, 'Pr (>F)'] <= 0.001) adonis_result_two[i, 'Sig'] <- '***'
  else if (adonis_result_two[i, 'Pr (>F)'] <= 0.01) adonis_result_two[i, 'Sig'] <- '**'
  else if (adonis_result_two[i, 'Pr (>F)'] <= 0.05) adonis_result_two[i, 'Sig'] <- '*'
}

# save result
sample_condition$condition <- factor(sample_condition$condition, levels = c("UC","nonIBD","CD"))
write.table(adonis_result_two, 'PERMANOVA.result_pairwise_999_all_conditions.txt', row.names = FALSE, sep = '\t', quote = FALSE, na = '')

pdf(file="Pcoa_MT_t4p.pdf",width = 5,height = 6)
ggscatter(sample_condition, x= "PCoA1", y = "PCoA2",
          color = "condition",
          #shape = "condition",
          ellipse = TRUE,  
          mean.point = TRUE, star.plot = T,   
          ellipse.level = 0.95,  
          ggtheme = theme_minimal()) +
  labs(title="IBD_HMP2_metagenome",
       subtitle = "CD vs nonIBD:P-value= 0.005
UC vs nonIBD:P-value= 0.006",
       x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'),
       y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%'))+
  theme_classic()+
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  #geom_point(size = 4)+ 
  theme(panel.grid = element_line(color = 'black', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank())+
  theme(axis.title = element_text(size = 18, colour="black"),
        axis.text = element_text(size = 16, colour = "black"),
        legend.text = element_text(size = 12),legend.position="bottom")+scale_color_lancet()

dev.off()

















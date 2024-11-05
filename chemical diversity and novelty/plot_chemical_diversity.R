##Multidimensional Scaling
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
# read the the data of distances between 9 RiPP subfamilies
dist.RiPP <- read.csv("representaive_pairwise_dist.tsv",sep = "\t")
dist.RiPP
row.names(dist.RiPP) <- dist.RiPP[, 1]
dist.RiPP <- dist.RiPP[, -1]
dist.RiPP
# read with-in family data
nodedata_1 <-read.csv("within_ripps_Tanimoto.tsv",sep = "\t", row.names = 1,  encoding="UTF-8")

precursor_used <- read.delim2("Clevaged_RiPP_pre_class_chemical_space.tsv")
class_no <- as.data.frame(table(precursor_used$prediction))

nodedata_2 <- merge(nodedata_1,class_no,by.x = "row.names",by.y="Var1")
row.names(nodedata_2) <- nodedata_2[, 1]
nodedata <- nodedata_2[, -1]
colnames(nodedata)[2] <- "Unique_sequence"


#run Multidimensional Scaling (MDS) with function cmdscale(), and get x and y coordinates
MDS <- cmdscale(dist.RiPP, eig = TRUE, k = 2)
site <- MDS$point
site <- data.frame(site)
site$name <- rownames(site)
merged=merge(site,nodedata,by="row.names",all.x=TRUE)
#set Continuous Color Range
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-0.5, 0.5))
fun_color_range <- colorRampPalette(c("#329337", "#dc2e2a")) 
my_colors <- fun_color_range(1000)
##visualize
pdf(file="MDS_cleavaged_231124.pdf",width = 8,height = 6)
p4 <- ggplot(data = merged, aes(X1, X2)) +
  scale_x_continuous(name= 'MDS1',limits=c(-0.5, 0.35))+ 
  scale_y_continuous(name= 'MDS2',limits=c(-0.35, 0.3))+
  geom_point(data = merged, aes(X1,X2, size= Unique_sequence, color= within_Tanimoto), show.legend = TRUE)+
  scale_size_continuous(range = c(5, 15))+
  geom_text_repel(data = merged, aes(label=name))+
  theme_linedraw(base_size = 15)+
  theme(legend.background = element_blank()) +
  theme(legend.key = element_blank())+
  theme(legend.key.size=unit(0.8,'cm'))+
  theme(legend.position = "top",
        legend.box = "horizontal")+
  theme(legend.title=element_blank())+
  #guides(colour=guide_legend(title="Intra-family Tanimoto coefficient"))+
  #sc
  scale_colour_gradientn(colors = my_colors, limits=c(0.5, 0.7))

p4
dev.off()
#save Fig
ggsave('Charting the chemical space of predicted RiPPs.tiff', p4, dpi = 600, width = 10, height = 8)

pdf(file="MDS_2.pdf",width = 8,height = 6)

p4 <- ggplot(data = merged, aes(X1, X2)) +
  scale_x_continuous(name= 'MDS1',limits=c(-0.35, 0.45))+ 
  scale_y_continuous(name= 'MDS2',limits=c(-0.2, 0.4))+
  geom_point(data = merged, aes(X1,X2, size= Unique_sequence, color= within_Tanimoto), show.legend = TRUE)+
  scale_size_continuous(range = c(5, 15))+
  geom_text_repel(data = merged, aes(label=name))+
  theme_linedraw(base_size = 15)+
  theme(legend.background = element_blank()) +
  theme(legend.key = element_blank())+
  theme(legend.key.size=unit(0.8,'cm'))+
  theme(legend.position = "bottom",
        legend.box = "horizontal")+
  theme(legend.title=element_blank())+
  #sc
  scale_colour_gradientn(colors = my_colors, limits=c(0.3, 0.45))

p4
dev.off()







####chordDiagram#####
library(dplyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(circlize)
library(viridis)
load("Prioritized_family_domain_class_230922.Rdata")

data <- RiPP_family_domain_class
data[is.na(data)] <- "NA"
table(data$max_domain)
data[which(data$max_domain == "Auto_inducing_peptide"), "max_domain"]="AgrD"
data[which(data$max_domain == "rSAM_modified_RiPP"), "max_domain"]="rSAM_modified_RiPP"
data[which(data$max_domain == "Graspetide"), "max_domain"]="1"
data[which(data$max_domain == "Lanthipeptide"), "max_domain"]="2"
data[which(data$max_domain == "Lasso_peptide"), "max_domain"]="3"
data[which(data$max_domain == "Thiopeptide"), "max_domain"]="4"
data[which(data$max_domain == "RiPP_like"), "max_domain"]="5"
data[which(data$max_domain == "LAP"), "max_domain"]="6"
data[which(data$max_domain == "Other_Known_RiPP"), "max_domain"]="7"

table(data$max_domain)
table(data$class)
stat_d <- as.data.frame(table(data$family_class,data$max_domain))
write.table(stat_d,file="RiPP_precursor_families_novelty.tsv",sep="\t",quote=F,row.names = F)

data_long <- stat_d

# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 1, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

# color palette
mycolor <- viridis(20, alpha = 1, begin = 0, end = 1, option = "D")
mycolor <- mycolor[sample(1:20)]

# Base plot
pdf(file="CD_pre_domain_MIBIG_random_20230922_8.pdf",width=12,height=12)
circos.par(start.degree = 90)
chordDiagram(
  x = data_long, 
  grid.col = mycolor,
  order = c("AUTO_INDUCING_PEPTIDE","BACTERIAL_HEAD_TO_TAIL_CYCLIZED",
            "LANTHIPEPTIDE","CYANOBACTIN","LASSO_PEPTIDE","GRASPETIDE","rSAM_MODIFIED_RiPP","THIOPEPTIDE","OTHER",
            "NA","1","2","3","4","5","6","7","rSAM_modified_RiPP","AgrD","MIBIG"),
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.05,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE,
  big.gap = 10)
#small.gap=1)
#abline(v = 0, lty = 2, col = "#00000080")
# Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 4, 
      labels = sector.index, 
      facing = "bending", 
      adj = c(0.3, 0),
      cex = 1
      #adj = c(0.5, 1)
    )
    
    # Add graduation on axis
    circos.axis(
      h = "top", 
      #major.at = seq(from = 5, to = xlim[100], by = ifelse(test = xlim[100]>200, yes = 200, no = 1)), 
      #major.at = seq(from = 0, to = xlim[5], by = ifelse(test = xlim[5]>10, yes = 10, no = 1)), 
      minor.ticks = 10, 
      major.tick.percentage = 0.5,
      labels.niceFacing = T)
  }
)
dev.off()
table(data$family_class)
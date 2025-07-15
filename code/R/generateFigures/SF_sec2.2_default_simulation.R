library(ggplot2)
library(ggpubr)
rm(list = ls())
path <- "C:/Users/Rudan/Documents/GitHub/KineticGP/"
# path <- "/home/mpimp-golm.mpg.de/xu2004/KineticGP/"
setwd(path)
library(png)
library(grid)
library(gridExtra)
library(ggbeeswarm)
library(ggrepel)
library(cowplot)


plot1 <- readPNG('figures/tempFigures/SF2_initial_simACa22and23_original.png')

plot2 <- readPNG('figures/tempFigures/SF2_initial_simGsCa22and23_original.png')

plot3 <- readPNG('figures/tempFigures/SF2_initial_simAQ22and23_original.png')

tmp <- arrangeGrob(rasterGrob(plot1),rasterGrob(plot2),rasterGrob(plot3),ncol = 3)

p <- as_ggplot(tmp) 


ggsave('figures/SFsec2.2_final.png',p, width=8,height=5)

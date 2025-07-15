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

source("code/R/generateFigures/MF2a_Sensitivity.R")
source("code/R/generateFigures/MF2b_parameterization_chi2.R")

fig2a <- plotmf2a()
fig2b <- plotmf2b()

fig_ab <- ggarrange(fig2a,fig2b,labels = c('a', 'b'),widths = c(1,1.2))
          
fig_ab

plot3 <- readPNG('figures/tempFigures/MF2b_new.png')

tmp <- arrangeGrob(fig_ab,rasterGrob(plot3), 
                   layout_matrix = rbind(c(1, 1), c(2, 2)),
                   heights = c(1, 1),  # Adjust relative row heights (smaller value reduces space between rows)
                   widths = c(1, 1))

p <- as_ggplot(tmp) +                                # transform to a ggplot
  draw_plot_label(label = c("", "c"), size = 14,
                  x = c(0, 0), y = c(1.3,  0.5)) # Add labels

p

ggsave('figures/MF2_final_new.png',p, width=8,height=10)

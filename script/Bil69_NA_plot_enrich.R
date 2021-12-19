#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
require(cowplot)

PlotCompareFit_Rep <- function(Bil69_data, graphname){
    textsize <- 8
    p <- ggplot() +
           geom_rect(data=NULL,aes(xmin=log10(2),xmax=Inf,ymin=log10(2),ymax=Inf), color=NA, fill=alpha('grey60', 0.5)) +
           geom_point(data=Bil69_data, aes(x=log10(Rep1Enrich), y=log10(Rep2Enrich)), pch=16, size=0.6, color='black', alpha=0.5) +
           #scale_color_manual(values=c('grey30'),drop=FALSE) +
           theme_cowplot(12) +
           theme(axis.title=element_text(size=textsize,face="bold"),
                 axis.text=element_text(size=textsize,face="bold"),
                 legend.title=element_blank(),
                 legend.key.size=unit(0.1,'in'),
                 legend.spacing.x=unit(0.03, 'in'),
                 legend.text=element_text(size=textsize,face="bold"),
                 legend.position='top') +
           labs(x=expression(bold(log['10']~enrich~'(Rep 1)')),y=expression(bold(log['10']~enrich~'(Rep 2)')))
    ggsave(graphname,p,height=2,width=2,dpi=600, bg='white')
    } 

Bil69_data <- read_tsv('result/Bil69_MultiMutLib_filtered.tsv')
PlotCompareFit_Rep(Bil69_data,'graph/Bil69_mutlib_rep_compare.png')

#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(qualpalr)
require(cowplot)

PlotCompareFit_Rep <- function(SD93_data, graphname){
    textsize <- 7
    palette <- qualpal(n = 2, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex
    p <- ggplot() +
           geom_rect(data=NULL,aes(xmin=log10(0.001),xmax=0,ymin=log10(0.001),ymax=0), color=NA, fill=alpha('grey60', 0.5)) +
           geom_point(data=SD93_data, aes(x=log10(Rep1Freq), y=log10(Rep2Freq), color=class), pch=16,size=0.4,alpha=0.75) +
           scale_color_manual(values=c('brown','grey30'),drop=FALSE) +
           theme_cowplot(12) +
           theme(axis.title=element_text(size=textsize,face="bold"),
                 axis.text=element_text(size=textsize,face="bold"),
                 legend.title=element_blank(),
                 legend.text=element_text(size=textsize,face="bold"),
                 legend.position='right',
                 legend.justification='center',
                 legend.spacing.x=unit(0.03, 'in'),
                 legend.key.size = unit(0.3,"line")) +
           labs(x=expression(bold(log['10']~freq~'(Rep 1)')),y=expression(bold(log['10']~freq~'(Rep 2)')))
    ggsave(graphname,p,height=1.5,width=2.5,dpi=600,bg='white')
    } 

coloring <- function(denovo){
  if (grepl('Y336N', denovo)){
    return ('with Y336N')
    }
  else {
    return ('without Y336N')
    }
  }

SD93_data <- read_tsv('result/SD93_MultiMutLib_filtered.tsv') %>%
               mutate(class=mapply(coloring, denovo))
PlotCompareFit_Rep(SD93_data,'graph/SD93_mutlib_rep_compare.png')

rm(list=ls(all=TRUE))
###################################################
## Basic
###################################################

CURRENT_WORKING_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
CURRENT_WORKING_DIR
setwd(CURRENT_WORKING_DIR)

###################################################
## Install package
###################################################
{
  library(xlsx)
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(survminer)
  library(survival)
  library(ggpmisc)
  library(ggsignif)
}
###################################################
## Graph option
###################################################
{
  pltsetb <- theme(legend.position="none") + 
    theme( plot.background = element_blank() ) +
    theme( text=element_text("Helvetica") ) +
    theme( axis.text = element_text(colour = "black", size = 15)) +
    theme( axis.ticks = element_line(colour = "black", size = 1) ) +
    theme( axis.ticks.length = unit(0.2,"cm") ) + 
    theme( axis.line.x = element_line(colour = "black", size = 1) ) + 
    theme( axis.line.y = element_line(colour = "black", size = 1) ) +
    theme( axis.title = element_text(colour = "black", size = 15) ) +
    theme( panel.background=element_blank() ) +
    theme( panel.grid.major=element_blank() ) +
    theme( panel.grid.minor=element_blank() )
}

###################################################
## PLOT
###################################################

g1df <- read.csv("../../Data/g1_feature_count.csv") %>% 
  mutate(freq=count/50,
         feature_name=substr(feature_name, 3, nchar(feature_name))) %>% 
  filter(freq>=0.25)
g3df <- read.csv("../../Data/g3_feature_count.csv") %>% 
  mutate(freq=count/50,
         feature_name=substr(feature_name, 3, nchar(feature_name))) %>% 
  filter(freq>=0.25)

plt <- ggplot(g1df, aes(y=reorder(feature_name,count),x=freq))+
  geom_bar(stat = "identity",fill="#009966")+
  labs(x="Frequency",y="")+
  pltsetb
ggsave("output/Fig3C.png",plt,width = 6,height = 6)

plt <- ggplot(g3df, aes(y=reorder(feature_name,count),x=freq))+
  geom_bar(stat = "identity",fill="#CC3366")+
  labs(x="Frequency",y="")+
  pltsetb
ggsave("output/Fig3D.png",plt,width = 6,height = 6)

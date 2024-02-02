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
  ##グラフのオプション
  caka <- "#CC3366"; cakao <- "#CC336650";
  cao <- "#3366CC"; caoo <- "#3366CC50";
  cmidori <- "#009966"; cmidorio <- "#00996650";
  corange <- "#FF6600"; corangeo <- "#FFCC33";
  cmurasaki <- "#6600FF"; cmurasakio <- "#CC66FF"
  cpink <- "deeppink3"; cpinko <- "deeppink1"
  ccha <- "#BB4513"; cchao <- "#CD661D"
  paka <- "#ff7f7f"; pakao <- "#ff9e9e"
  pao <- "#8484ff"; paoo <- "#a3a3ff"
  pmidori <- "#7fff7f"; pmidorio <- "#9eff9e"
  
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

###Fig1C###
df.umap <- read.csv("../../Data/clustering_umap.csv")

plt <- ggplot(df,aes(x=X0,y=X1,color=cluster))+
  geom_point(size=4,alpha=.8)+
  scale_color_manual(values = c("G1"="#00CC00",
                                "G2"="#0000FF",
                                "G3"="#FF0000"))+
  labs(x="UMAP1",y="UMAP2")+
  pltsetb+
  #theme( axis.ticks = element_blank() )+
  theme( axis.text = element_blank())
ggsave("output/Fig1B.png",width = 4,height = 4)

###Fig1D###

load("../../Data/phase1.rda")
load("../../Data/phase2.rda")
load("../../Data/phase3.rda")
load("../../Data/phase4.rda")

ColFun <- colorRamp2(breaks = c(0,0.5,1), 
                     colors = c("#cc6600","#cc0066","#00cccc"),
                     transparency = 0,
                     space = "RGB")
cols <- ColFun(seq(0, 1, length = 7))
ColFun2 <- colorRamp2(breaks = c(0,0.5,1), 
                      colors = c("#cc660050","#cc006650","#00cccc50"),
                      transparency = 0,
                      space = "RGB")
coltrans <- ColFun2(seq(0, 1, length = 7))
# plot
lwd <- 2.5
alpha_ribbon <- 0.1
alpha_line <- 0.75
clnames <- c("G1","G3","G2")
{
  #G1
  i=1
  plt1 <- ggplot()+
    geom_ribbon(data = L1.df[[i]], aes(x=time, ymin=log10(lower), ymax=log10(upper)), fill=coltrans[1], alpha=alpha_ribbon)+
    geom_line(data = L1.df[[i]], aes(x=time, y=log10(R1)), lwd=lwd, color = cols[1],alpha=0.75)+
    geom_ribbon(data = L2.df[[i]], aes(x=time, ymin=log10(lower), ymax=log10(upper)), fill=coltrans[4], alpha=alpha_ribbon)+
    geom_line(data = L2.df[[i]], aes(x=time, y=log10(R1)), lwd=lwd, color = cols[4],alpha=0.75)+
    geom_ribbon(data = L3.df[[i]], aes(x=time, ymin=log10(lower), ymax=log10(upper)), fill=coltrans[5], alpha=alpha_ribbon)+
    geom_line(data = L3.df[[i]], aes(x=time, y=log10(R1)), lwd=lwd, color = cols[5],alpha=0.75)+
    geom_ribbon(data = L4.df[[i]], aes(x=time, ymin=log10(lower), ymax=log10(upper)), fill=coltrans[6], alpha=alpha_ribbon)+
    geom_line(data = L4.df[[i]], aes(x=time, y=log10(R1)), lwd=lwd, color = cols[6],alpha=0.75)+
    labs(x="Time after lesion onset (Days)",y="Log10 Lesion count")+
    scale_y_continuous(limits = c(0,3.2),breaks = seq(0,3,1))+
    pltsetb
  #G3
  i=2
  plt2 <- ggplot()+
    geom_ribbon(data = L1.df[[i]], aes(x=time, ymin=log10(lower), ymax=log10(upper)), fill=coltrans[1], alpha=alpha_ribbon)+
    geom_line(data = L1.df[[i]], aes(x=time, y=log10(R2)), lwd=lwd, color = cols[1],alpha=0.75)+
    geom_ribbon(data = L2.df[[i]], aes(x=time, ymin=log10(lower), ymax=log10(upper)), fill=coltrans[4], alpha=alpha_ribbon)+
    geom_line(data = L2.df[[i]], aes(x=time, y=log10(R2)), lwd=lwd, color = cols[4],alpha=0.75)+
    geom_ribbon(data = L3.df[[i]], aes(x=time, ymin=log10(lower), ymax=log10(upper)), fill=coltrans[5], alpha=alpha_ribbon)+
    geom_line(data = L3.df[[i]], aes(x=time, y=log10(R2)), lwd=lwd, color = cols[5],alpha=0.75)+
    geom_ribbon(data = L4.df[[i]], aes(x=time, ymin=log10(lower), ymax=log10(upper)), fill=coltrans[6], alpha=alpha_ribbon)+
    geom_line(data = L4.df[[i]], aes(x=time, y=log10(R2)), lwd=lwd, color = cols[6],alpha=0.75)+
    labs(x="Time after lesion onset (Days)",y="Log10 Lesion count")+
    scale_x_continuous(breaks=seq(0,40,by=10), limits=c(0,39.99))+
    scale_y_continuous(limits = c(0,3.2),breaks = seq(0,3,1))+
    pltsetb
  #G2
  i=3
  plt3 <- ggplot()+
    geom_ribbon(data = L1.df[[i]], aes(x=time, ymin=log10(lower), ymax=log10(upper)), fill=coltrans[1], alpha=alpha_ribbon)+
    geom_line(data = L1.df[[i]], aes(x=time, y=log10(R3)), lwd=lwd, color = cols[1],alpha=0.75)+
    geom_ribbon(data = L2.df[[i]], aes(x=time, ymin=log10(lower), ymax=log10(upper)), fill=coltrans[4], alpha=alpha_ribbon)+
    geom_line(data = L2.df[[i]], aes(x=time, y=log10(R3)), lwd=lwd, color = cols[4],alpha=0.75)+
    geom_ribbon(data = L3.df[[i]], aes(x=time, ymin=log10(lower), ymax=log10(upper)), fill=coltrans[5], alpha=alpha_ribbon)+
    geom_line(data = L3.df[[i]], aes(x=time, y=log10(R3)), lwd=lwd, color = cols[5],alpha=0.75)+
    geom_ribbon(data = L4.df[[i]], aes(x=time, ymin=log10(lower), ymax=log10(upper)), fill=coltrans[6], alpha=alpha_ribbon)+
    geom_line(data = L4.df[[i]], aes(x=time, y=log10(R3)), lwd=lwd, color = cols[6],alpha=0.75)+
    labs(x="Time after lesion onset (Days)",y="Log10 Lesion count")+
    scale_x_continuous(breaks=seq(0,40,by=10), limits=c(0,39.99))+
    scale_y_continuous(limits = c(0,3.2),breaks = seq(0,3,1))+
    pltsetb
}
all <- plt1+plt3+plt2+plot_layout(ncol = 3)
ggsave("output/Fig1C.png",all,width =15,height = 4.5)

###Fig1D###
df <- read.csv("../../Data/lesion_count.csv")
plot.df <- df %>% 
  pivot_longer(., cols = colnames(df)[2:5],names_to = 'type', values_to = 'val')

#各状態のAUC
plt <- ggplot(plot.df,aes(x=type,y=log10(val),fill=cluster))+
  geom_boxplot()+
  scale_fill_manual(values = c("G1"="#00CC00",
                               "G3"="#FF0000",
                               "G2"="#0000FF"))+
  labs(x="Lesion Stage",y="Log10 Lesion count")+
  scale_x_discrete(labels=c("p1auc"="Stage1",
                            "p2auc"="Stage2",
                            "p3auc"="Stage3",
                            "p4auc"="Stage4"))+
  scale_y_continuous(limits = c(1.5,4.5),breaks = seq(2,4,1))+
  pltsetb
ggsave("output/fig1D.png",width = 6,height = 4)

###Fig1E###
df <- read.csv("../../Data/enginnering_parameter.csv")

#Total Lesion
t.res <- pairwise.wilcox.test(x=df$total.lesion,g=df$cluster,paired=F,p.adjust.method="bonferroni")
pval <- data.frame(pname=c("AB","BC","CA"),
                   p.value=c(t.res$p.value[1],t.res$p.value[4],t.res$p.value[2])) %>% 
  mutate(pmark=ifelse(p.value>0.05,"N.S.",
                      ifelse(p.value>0.01,"*",
                             ifelse(p.value>0.001,"**",
                                    ifelse(p.value>0.0001,"***","****")))))


plt <- ggplot(df,aes(x=cluster,y=total.lesion,fill=cluster))+
  geom_boxplot()+
  scale_fill_manual(values = c("G1"="#00CC00",
                               "G3"="#FF0000",
                               "G2"="#0000FF"))+
  labs(x="",y="Log10 Total lesion count")+
  scale_y_continuous(limits = c(2.3,5.4),breaks = seq(3,5,1))+
  geom_signif(comparisons = list(c("G1","G2"),c("G2","G3"),c("G1","G3")),
              step_increase = 0.1,annotations = pval$pmark,
              size = .5,textsize = 5,color="black")+
  pltsetb
ggsave("output/Fig1E-1.png",plt,width=3,height=4)

#Lesion duration
t.res <- pairwise.wilcox.test(x=df$lesion.duration,g=test$cluster,paired=F,p.adjust.method="bonferroni")
pval <- data.frame(pname=c("AB","BC","CA"),
                   p.value=c(t.res$p.value[1],t.res$p.value[4],t.res$p.value[2])) %>% 
  mutate(pmark=ifelse(p.value>0.05,"N.S.",
                      ifelse(p.value>0.01,"*",
                             ifelse(p.value>0.001,"**",
                                    ifelse(p.value>0.0001,"***","****")))))

plt <- ggplot(df,aes(x=cluster,y=lesion.duration,fill=cluster))+
  geom_boxplot()+
  scale_fill_manual(values = c("G1"="#00CC00",
                               "G3"="#FF0000",
                               "G2"="#0000FF"))+
  labs(x="",y="Lesion duration (Days)")+
  scale_y_continuous(limits = c(16,78),breaks = seq(20,60,20),labels = expression(20,40,60))+
  geom_signif(comparisons = list(c("G1","G2"),c("G2","G3"),c("G1","G3")),
              step_increase = 0.1,annotations = pval$pmark,
              size = .5,textsize = 5,color="black")+
  pltsetb
ggsave("output/Fig1E-2.png",plt,width=3,height=4)


###Fig1F###

#reconstruct 
df <- read.csv("../../Data/enginnering_parameter.csv")[,1:2] %>% 
  mutate(lesion.surv=lesion.duration,
         event=1) %>% 
  dplyr::select("lesion.surv","event","cluster")
df2 <- df %>% mutate(cluster="Total")

dz <- rbind(df,df2)

survobj <- Surv(time = dz$lesion.surv,
                event = dz$event)
linelistsurv_fit <- survfit(survobj ~ cluster,data=dz)

pltemp2 <- ggsurvplot(linelistsurv_fit,
                      conf.int = T,
                      palette=c("#00CC00","#0000FF","#FF0000","black"),
                      risk.table = "abs_pct",
                      surv.scale="percent",
                      legend.labs=c("G1","G2","G3","G4"),
                      legend.title="Group",
                      xlab="Onset of lesion (Days)")
plt <- pltemp2$plot+
  pltsetb
ggsave("output/Fig1F.png",plt,width=6,height=4)


rm(list=ls(all=TRUE))
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
  library(FME)
  library(patchwork)
  library(ggsignif)
  library(ggpmisc)
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
## Plotting
###################################################

### Fig2B ###
df <- read.csv("../../Data/virus_indpar.csv")

Tmin <- 0.0
Tmax <- 100
step_size <- 0.01
stime <- seq(Tmin,Tmax,step_size)

Covfun <- function(pars){
  
  V0 <- as.numeric(pars["V0"])
  delta <- as.numeric(pars["delta"])
  
  derivs<-function(time, y, pars){
    with(as.list(c(pars, y)),{
      dV <- -delta*V
      
      return(list(c(dV)))
    })
  }
  y<-c(V=V0)
  
  times<-c(seq(Tmin,Tmax,step_size))
  out<-lsoda(y=y, parms=pars, times=times, func=derivs, rtol=0.00004, atol=0.00000000000001)
  out2<-cbind(time=out[,1],aV=((log10(out[,2]))))
  as.data.frame(out2)
}
alldf <- c()
for (i in 1:nrow(df)){
  ipar <- df[i,]
  cl <- ipar$col
  
  sR_pars <- c(V0 = ipar$V0, delta = ipar$delta)
  fitted <- Covfun(sR_pars)
  kari <- data.frame(id=rep(i,nrow(fitted)),
                     time=fitted$time,
                     aV=fitted$aV,
                     class=rep(df$cluster[i],nrow(fitted))) 
  alldf <- rbind(alldf,kari)
  
}

alldf$class <- factor(alldf$class, levels = c("0", "1", "2"))
pall <- ggplot(alldf, aes(x = time, y = aV, group = id, color = factor(class))) +
  geom_line(lwd = 1, alpha = 0.3) +
  scale_x_continuous(breaks = seq(0, 100, by = 25), limits = c(0, 100)) +
  scale_y_continuous(breaks = seq(-2, 8, by = 2),
                     labels = expression(10^(-2), 10^0, 10^2, 10^4, 10^6, 10^8),
                     limits = c(-0.5, 8.5)) +
  scale_color_manual(values = c("0" = "#00CC00", "1" = "#FF0000", "2" = "#0000FF")) +
  labs(x="Time after lesion onset (Days)",y="Viral load (genomes/mL)")+
  pltsetb
ggsave("output/Fig2B.png",pall,dpi=300,height = 4,width = 6)

### Fig2C ###
df <- read.csv("../../Data/enginnering_parameter.csv")

pt <- 4
# VL at onset Lesion
value <- df$onset.VL
class <- df$cluster
t.res <- pairwise.wilcox.test(x=value,g=class,paired=F,p.adjust.method="bonferroni")
pval <- data.frame(pname=c("AB","BC","CA"),
                   p.value=c(t.res$p.value[1],t.res$p.value[4],t.res$p.value[2])) %>% 
  mutate(pmark=ifelse(p.value>0.05,"N.S.",
                      ifelse(p.value>0.01,"*",
                             ifelse(p.value>0.001,"**",
                                    ifelse(p.value>0.0001,"***","****")))))

plt1 <- ggplot(df,aes(x=cluster,y=onset.VL,fill=cluster))+
  geom_boxplot()+
  labs(y="Viral load at the onset\n(log10 genomes/mL)",x="")+
  scale_fill_manual(values = c("G1"="#00CC00",
                               "G3"="#FF0000",
                               "G2"="#0000FF"))+
  scale_y_continuous(limits = c(3,7),breaks = seq(3,6,1))+
  geom_signif(comparisons = list(c("G1","G2"),c("G2","G3"),c("G1","G3")),
              step_increase = 0.1,annotations = pval$pmark,
              size = .5,textsize = pt,color="black")+
  pltsetb

# VL clearance
value <- df$delta
class <- df$cluster
t.res <- pairwise.wilcox.test(x=value,g=class,paired=F,p.adjust.method="bonferroni")
pval <- data.frame(pname=c("AB","BC","CA"),
                   p.value=c(t.res$p.value[1],t.res$p.value[4],t.res$p.value[2])) %>% 
  mutate(pmark=ifelse(p.value>0.05,"N.S.",
                      ifelse(p.value>0.01,"*",
                             ifelse(p.value>0.001,"**",
                                    ifelse(p.value>0.0001,"***","****")))))

plt2 <- ggplot(df,aes(x=cluster,y=delta,fill=cluster))+
  geom_boxplot()+
  labs(y="Clearance rate (/day)",x="")+
  scale_fill_manual(values = c("G1"="#00CC00",
                               "G3"="#FF0000",
                               "G2"="#0000FF"))+
  scale_y_continuous(limits = c(0.146,0.159),breaks = seq(0.146,0.158,0.004))+
  geom_signif(comparisons = list(c("G1","G2"),c("G2","G3"),c("G1","G3")),
              step_increase = 0.1,annotations = pval$pmark,
              size = .5,textsize = pt,color="black")+
  pltsetb

# VL AUC
value <- df$AUC
class <- df$cluster
t.res <- pairwise.wilcox.test(x=value,g=class,paired=F,p.adjust.method="bonferroni")
pval <- data.frame(pname=c("AB","BC","CA"),
                   p.value=c(t.res$p.value[1],t.res$p.value[4],t.res$p.value[2])) %>% 
  mutate(pmark=ifelse(p.value>0.05,"N.S.",
                      ifelse(p.value>0.01,"*",
                             ifelse(p.value>0.001,"**",
                                    ifelse(p.value>0.0001,"***","****")))))
plt3 <- ggplot(df,aes(y=AUC,x=cluster,fill=cluster))+
  geom_boxplot()+
  labs(y="Total viral load\n(log10 genomes/mL)",x="")+
  scale_fill_manual(values = c("G1"="#00CC00",
                               "G3"="#FF0000",
                               "G2"="#0000FF"))+
  scale_y_continuous(limits = c(30,380),breaks = seq(0,300,100))+
  geom_signif(comparisons = list(c("G1","G2"),c("G2","G3"),c("G1","G3")),
              step_increase = 0.1,annotations = pval$pmark,
              size = .5,textsize = pt,color="black")+
  pltsetb

# VL at disappear lesion
value <- df$lesion.disapper.VL
class <- df$cluster
t.res <- pairwise.wilcox.test(x=value,g=class,paired=F,p.adjust.method="bonferroni")
pval <- data.frame(pname=c("AB","BC","CA"),
                   p.value=c(t.res$p.value[1],t.res$p.value[4],t.res$p.value[2])) %>% 
  mutate(pmark=ifelse(p.value>0.05,"N.S.",
                      ifelse(p.value>0.01,"*",
                             ifelse(p.value>0.001,"**",
                                    ifelse(p.value>0.0001,"***","****")))))
plt4 <- ggplot(df,aes(x=cluster,y=lesion.disapper.VL,fill=cluster))+
  geom_boxplot()+
  labs(y="Viral load at lesion diappearance\n(log10 genomes/mL)",x="")+
  scale_fill_manual(values = c("G1"="#00CC00",
                               "G3"="#FF0000",
                               "G2"="#0000FF"))+
  scale_y_continuous(limits = c(0,5),breaks = seq(0,4,1))+
  geom_signif(comparisons = list(c("G1","G2"),c("G2","G3"),c("G1","G3")),
              step_increase = 0.1,annotations = pval$pmark,
              size = .5,textsize = pt,color="black")+
  pltsetb

all <- plt1+plt2+plt3+plt4+plot_layout(ncol = 2)
ggsave("output/Fig2C.png",all,width=8,height=8)

### Fig2D ###
df <- read.csv("../../Data/enginnering_parameter.csv")

size_n <- 2
alpha_n <- .5
t.size <- 4

plt1 <- ggplot(df,aes(y=lesion.duration,x=lesion.disapper.VL))+
  geom_point(size=size_n,alpha=alpha_n) +
  geom_smooth(method = "lm",formula = y~x,lwd=2,color="black",fill="gray")+
  stat_correlation(use_label(c("R","P")), method="pearson",size=t.size)+
  labs(y="Lesion duration (Days)",x="Viral load at lesion disappearance\n(log10 genomes/mL)")+
  scale_y_continuous(limits = c(16,70),breaks = seq(20,60,20))+
  scale_x_continuous(limits = c(0.5,4.05),breaks = seq(-1,4,1))+
  pltsetb

plt2 <- ggplot(df,aes(x=onset.VL,y=lesion.duration))+
  geom_point(size=size_n,alpha=alpha_n) +
  geom_smooth(method = "lm",formula = y~x,lwd=2,color="black",fill="gray")+
  stat_correlation(use_label(c("R","P")), method="pearson",size=t.size)+
  labs(x="Viral load at lesion onset\n(log10 genomes/mL)",y="Lesion duration (Days)")+
  scale_x_continuous(limits = c(3.5,6.5))+
  scale_y_continuous(limits = c(16,70),breaks = seq(20,60,20))+
  pltsetb

plt3 <- ggplot(df,aes(x=onset.VL,y=total.lesion))+
  geom_point(size=size_n,alpha=alpha_n) +
  geom_smooth(method = "lm",formula = y~x,lwd=2,color="black",fill="gray")+
  stat_correlation(use_label(c("R","P")), method="pearson",size=t.size)+
  labs(x="Viral load at lesion onset\n(log10 genomes/mL)",y="Total lesion count\n(log10 lesion count)")+
  scale_x_continuous(limits = c(3.5,6.5))+
  scale_y_continuous(limits = c(2.2,5))+
  pltsetb

all <- plt1+plt2+plt3+plot_layout(nrow=1)
ggsave("output/Fig2D.png",all,width = 12,height = 4.5)

### Fig2H ###
df <- read.csv("../../Data/enginnering_parameter.csv")

plt <- ggplot(df,aes(x=onset.VL,y=factor(cluster,levels = c("G3","G2","G1")),color=cluster))+
  geom_vline(xintercept = 4.6,lwd=.75,lty=2,alpha=.5)+
  geom_beeswarm(size=1.5,cex=2.5)+
  labs(x="Viral load at lesion onset (log10 genomes/mL)",y="")+
  scale_color_manual(values=c("G1"="#009966","G2"="#3366CC","G3"="#CC3366"))+
  pltsetb
ggsave("output/Fig2G.png",plt,width = 5,height = 2)

### Fig2H ###
df <- read.csv("../../Data/enginnering_parameter.csv")%>% 
  mutate(V_class=ifelse(onset.VL>4.6,"High","Low"))

pltz <- ggplot(df,aes(x=factor(V_class,levels = c("High","Low")),y=lesion.duration,fill=V_class))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("High","Low")),
              map_signif_level = T,step_increase = 0.3,test = "wilcox.test",
              annotations = "****",size = .5,textsize = 5,color="black")+
  scale_fill_manual(values = c("Low"="deeppink3","High"="#6600FF"))+
  labs(x="",y="Lesion duration (Days)")+
  coord_flip()+
  pltsetb
pltz
ggsave("output/Fig2H.png",pltz,width = 5,height = 2)

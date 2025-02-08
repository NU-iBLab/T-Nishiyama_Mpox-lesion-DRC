
rm(list=ls(all=TRUE))
CURRENT_WORKING_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
CURRENT_WORKING_DIR
setwd(CURRENT_WORKING_DIR)

###################################################
## Install package
###################################################
{
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
  library(FME)
  library(patchwork)
  library(ggsignif)
  library(ggpmisc)
  library(ggbeeswarm)
  library(pROC)
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
df <- read.csv("../../Data/VL_profile.csv")[,1:4] %>% 
  mutate(cluster=ifelse(cluster==0,0,1))

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
for (i in 1:length(df$id)){
  ipar <- df[i,]
  cl <- ipar$col
  
  sR_pars <- c(V0 = ipar$V0, delta = ipar$delta)
  fitted <- Covfun(sR_pars)
  kari <- data.frame(id=rep(df[i,"id"],nrow(fitted)),
                     time=fitted$time,
                     aV=fitted$aV,
                     class=rep(df$cluster[i],nrow(fitted))) 
  alldf <- rbind(alldf,kari)
  
}

alldf$class <- factor(alldf$class, levels = c("0", "1"))
pall <- ggplot(alldf, aes(x = time, y = aV, group = id, color = factor(class))) +
  geom_line(lwd = 1, alpha = 0.3) +
  scale_x_continuous(breaks = seq(0, 100, by = 25), limits = c(0, 100)) +
  scale_y_continuous(breaks = seq(-2, 8, by = 2),
                     labels = expression(10^(-2), 10^0, 10^2, 10^4, 10^6, 10^8),
                     limits = c(-0.5, 8.5)) +
  scale_color_manual(values = c("0" = "#0000FF", "1" = "#FF0000")) +
  labs(x="Time after lesion onset (Days)",y="Viral load (genomes/mL)")+
  pltsetb
ggsave("output/Fig2B.png",pall,dpi=300,height = 4,width = 6)

### Fig2C ###
df <- read.csv("../../Data/enginnering_parameter.csv")

pt <- 4
# VL at onset Lesion
t.res <- wilcox.test(df$onset.VL~df$cluster,paired=F,p.adjust.method="bonferroni")
pval <- data.frame(pname=c("AB"),
                   p.value=c(t.res$p.value[1])) %>% 
  mutate(pmark=ifelse(p.value>0.05,"N.S.",
                      ifelse(p.value>0.01,"*",
                             ifelse(p.value>0.001,"**",
                                    ifelse(p.value>0.0001,"***","****")))))

plt1 <- ggplot(df,aes(x=cluster,y=onset.VL,fill=cluster))+
  geom_boxplot()+
  labs(y="Viral load at lesion onset\n(log10 genomes/mL)",x="")+
  scale_fill_manual(values = c("G1"="#0000FF",
                               "G2"="#FF0000"))+
  scale_y_continuous(limits = c(3,6.5),breaks = seq(3,6,1))+
  geom_signif(comparisons = list(c("G1","G2")),
              step_increase = 0.1,annotations = pval$pmark,
              size = .5,textsize = pt,color="black")+
  pltsetb

# VL clearance
t.res <- wilcox.test(df$delta~df$cluster,paired=F,p.adjust.method="bonferroni")
pval <- data.frame(pname=c("AB"),
                   p.value=c(t.res$p.value[1])) %>% 
  mutate(pmark=ifelse(p.value>0.05,"N.S.",
                      ifelse(p.value>0.01,"*",
                             ifelse(p.value>0.001,"**",
                                    ifelse(p.value>0.0001,"***","****")))))

plt2 <- ggplot(df,aes(x=cluster,y=delta,fill=cluster))+
  geom_boxplot()+
  labs(y="Clearance rate (/day)",x="")+
  scale_fill_manual(values = c("G1"="#0000FF",
                               "G2"="#FF0000"))+
  #scale_y_continuous(limits = c(0.146,0.157),breaks = seq(0.146,0.156,0.002))+
  geom_signif(comparisons = list(c("G1","G2")),
              step_increase = 0.1,annotations = pval$pmark,
              size = .5,textsize = pt,color="black")+
  pltsetb
plt2
# VL AUC
t.res <- wilcox.test(df$AUC~df$cluster,paired=F,p.adjust.method="bonferroni")
pval <- data.frame(pname=c("AB"),
                   p.value=c(t.res$p.value[1])) %>% 
  mutate(pmark=ifelse(p.value>0.05,"N.S.",
                      ifelse(p.value>0.01,"*",
                             ifelse(p.value>0.001,"**",
                                    ifelse(p.value>0.0001,"***","****")))))
plt3 <- ggplot(df,aes(y=AUC,x=cluster,fill=cluster))+
  geom_boxplot()+
  labs(y="AUC of virus dynamics",x="")+
  scale_fill_manual(values = c("G1"="#0000FF",
                               "G2"="#FF0000"))+
  #scale_y_continuous(limits = c(30,320),breaks = seq(0,300,100))+
  geom_signif(comparisons = list(c("G1","G2")),
              step_increase = 0.1,annotations = pval$pmark,
              size = .5,textsize = pt,color="black")+
  pltsetb
plt3
# VL at disappear lesion
t.res <- wilcox.test(df$lesion.disapper.VL~df$cluster,paired=F,p.adjust.method="bonferroni")
pval <- data.frame(pname=c("AB"),
                   p.value=c(t.res$p.value[1])) %>% 
  mutate(pmark=ifelse(p.value>0.05,"N.S.",
                      ifelse(p.value>0.01,"*",
                             ifelse(p.value>0.001,"**",
                                    ifelse(p.value>0.0001,"***","****")))))
plt4 <- ggplot(df,aes(x=cluster,y=lesion.disapper.VL,fill=cluster))+
  geom_boxplot()+
  labs(y="Viral load at lesion diappearance\n(log10 genomes/mL)",x="")+
  scale_fill_manual(values = c("G1"="#0000FF",
                               "G2"="#FF0000"))+
  scale_y_continuous(limits = c(0,5),breaks = seq(0,4,1))+
  geom_signif(comparisons = list(c("G1","G2")),
              step_increase = 0.1,annotations = pval$pmark,
              size = .5,textsize = pt,color="black")+
  pltsetb
plt4
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
  scale_x_continuous(limits = c(0,4.2),breaks = seq(0,4,1))+
  pltsetb
plt1
plt2 <- ggplot(df,aes(x=onset.VL,y=lesion.duration))+
  geom_point(size=size_n,alpha=alpha_n) +
  geom_smooth(method = "lm",formula = y~x,lwd=2,color="black",fill="gray")+
  stat_correlation(use_label(c("R","P")), method="pearson",size=t.size)+
  labs(x="Viral load at lesion onset\n(log10 genomes/mL)",y="Lesion duration (Days)")+
  scale_x_continuous(limits = c(3.5,6.5))+
  scale_y_continuous(limits = c(16,70),breaks = seq(20,60,20))+
  pltsetb
plt2
plt3 <- ggplot(df,aes(x=onset.VL,y=total.lesion))+
  geom_point(size=size_n,alpha=alpha_n) +
  geom_smooth(method = "lm",formula = y~x,lwd=2,color="black",fill="gray")+
  stat_correlation(use_label(c("R","P")), method="pearson",size=t.size)+
  labs(x="Viral load at lesion onset\n(log10 genomes/mL)",y="Log10 Total AUC of lesion dynamics\n for each stage")+
  scale_x_continuous(limits = c(3.5,6.5))+
  scale_y_continuous(limits = c(2.2,5))+
  pltsetb
plt3
all <- plt1+plt2+plt3+plot_layout(nrow=1)
ggsave("output/Fig2D.png",all,width = 12,height = 4.5)

### Fig2F ###
df <- read.csv("../../Data/enginnering_parameter.csv") %>% 
  mutate(V_class=ifelse(onset.VL>4.61, "High", "Low"))

pltz <- ggplot(df,aes(x=factor(V_class,levels = c("High","Low")),y=lesion.duration,fill=V_class))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("High","Low")),
              map_signif_level = T,step_increase = 0.3,test = "wilcox.test",
              annotations = "****",size = .5,textsize = 5,color="black")+
  scale_fill_manual(values = c("Low"="deeppink3","High"="#6600FF"))+
  labs(x="Group of viral load at lesion onset",y="Lesion duration (Days)")+
  pltsetb
ggsave("output/Fig2F.png",pltz,width = 5,height = 5)

mean(df[df$V_class=="High","lesion.duration"])
mean(df[df$V_class=="Low","lesion.duration"])






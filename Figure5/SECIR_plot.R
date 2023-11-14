#' 2022-02-10
#' Script to produce a figure similar to figure 4 of 
#' https://arxiv.org/abs/2111.09803
#' author: Jacques Hermes <jacques.hermes@fdm.uni-freiburg.de>

# Sets current working directory to the location of this script
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))
library(numDeriv)
library(ggplot2)
library(ggthemes)
library(magrittr)
library(tidyr)
library(stringr)
library(dMod)
library(ggpubr)
library(lmf)
library(pracma)
library(fBasics)
library(cowplot)

load(file = "Result_SECIR.RData")
load(file = "Meas_SECIR.RData")
load(file = "Data_SECIR.RData")

cols <- c("True value" = dMod:::dMod_colors[4],
          "Incidence-based s=4 days, tau= 4 (RKI)" = dMod:::dMod_colors[2],
          "AKS method" = dMod:::dMod_colors[1],
          "Different values of s"= "lightgrey")

results$place <- "AKS method"

NumObs <- length(Meas$dI)

RKI2 <- rep(NA,NumObs)
for (i in 6:NumObs){
RKI2[i-1]<- sum(Meas$dI[i-0:3])/sum(Meas$dI[i-2:5])
}

RKI3 <- rep(NA,NumObs)
for (i in 7:NumObs){
RKI3[i-1]<- sum(Meas$dI[i-0:3])/sum(Meas$dI[i-3:6])
}

RKI4 <- rep(NA,NumObs)
for (i in 8:NumObs){
RKI4[i-1]<- sum(Meas$dI[i-0:3])/sum(Meas$dI[i-4:7])
}

RKI5 <- rep(NA,NumObs)
for (i in 9:NumObs){
RKI5[i-1]<- sum(Meas$dI[i-0:3])/sum(Meas$dI[i-5:8])
}

RKI6 <- rep(NA,NumObs)
for (i in 10:NumObs){
RKI6[i]<- sum(Meas$dI[i-0:3])/sum(Meas$dI[i-6:9])
}
RKI7 <- rep(NA,NumObs)
for (i in 11:NumObs){
RKI7[i]<- sum(Meas$dI[i-0:3])/sum(Meas$dI[i-7:10])
}

RKI8 <- rep(NA,NumObs)
for (i in 12:NumObs){
RKI8[i]<- sum(Meas$dI[i-0:3])/sum(Meas$dI[i-8:11])
}
RKI9 <- rep(NA,NumObs)
for (i in 13:NumObs){
RKI9[i]<- sum(Meas$dI[i-0:3])/sum(Meas$dI[i-9:12])
}

RKI10 <- rep(NA,NumObs)
for (i in 14:NumObs){
RKI10[i]<- sum(Meas$dI[i-0:3])/sum(Meas$dI[i-10:13])
}

RKI12 <- rep(NA,NumObs)
for (i in 19:NumObs){
RKI12[i]<- sum(Meas$dI[i-0:3])/sum(Meas$dI[i-12:15])
}

RKI <- data.frame(RKI2 = RKI2,RKI3 = RKI3, RKI4=RKI4,RKI5 = RKI5, RKI6=RKI6,RKI7 = RKI7, RKI8=RKI8,RKI9 = RKI9, RKI10=RKI10,RKI12=RKI12, time=1:NumObs)

#plot R0
p11 <- ggplot(results, aes(x=time, y=exp(Rt),col=place))+
ylab(expression(R[t], parse = TRUE))+xlab("Time ")+
geom_hline(yintercept = 1, linetype="dashed", colour="black")+
geom_line(data=RKI, aes(y=RKI2, x=time, col ="Different values of s"))+
geom_line(data=RKI, aes(y=RKI3, x=time, col ="Different values of s"))+
geom_line(data=RKI, aes(y=RKI4, x=time, col ="Incidence-based s=4 days, tau= 4 (RKI)"))+
geom_line(data=RKI, aes(y=RKI5, x=time, col ="Different values of s"))+
geom_line(data=RKI, aes(y=RKI6, x=time, col ="Different values of s"))+      
geom_line(data=RKI, aes(y=RKI7, x=time, col ="Different values of s"))+
geom_line(data=RKI, aes(y=RKI8, x=time, col ="Different values of s"))+
geom_line(data=RKI, aes(y=RKI9, x=time, col ="Different values of s"))+
geom_line(data=RKI, aes(y=RKI10, x=time, col ="Different values of s"))+
geom_line(data=data, aes(x=time,y=R_t, col ="True value"))+
geom_line()+
scale_y_log10()+
scale_color_manual(values=cols)+theme_dMod(base_size = 10)+ theme(legend.text = element_text(size=10),legend.title = element_blank(),legend.position = c(0.65, 0.8),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+geom_ribbon(data=results,aes(ymin=exp(Rt-sRt),ymax=exp(Rt+sRt)), color = NA,alpha=0.1)

ggsave(paste0("SECIR_figure.pdf"), p11,height = 4.5 , width = 5)
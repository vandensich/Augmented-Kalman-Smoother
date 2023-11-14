#' 2022-02-10
#' Script to produce a figure similar to figure 5 of 
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

load( file = "Result_Germany.RData")
load(file = "Data_multi_7_SIR_aktuell_August.RData")

i<-1
Parameter <- data.frame(Rt = exp(results$Rt))
Parameter$gamma <- exp(results$gamma)
Parameter$theta <- exp(results$theta)
Parameter$time <- results$time
Parameter$I <- results$I
Parameter$R <- results$R
Parameter$D <- results$D
Parameter$dI <- results$dI
Parameter$dD <- results$dD
Parameter$dR <- results$dR

save(Parameter, file = "Parameter_Bestfit_Germany.RData")

results$place <- "AKS method"

R0_rki <- 0

#Data taken from 
#https://github.com/robert-koch-institut/SARS-CoV-2-Nowcasting_und_-R-Schaetzung
#Whihch is the officially published R-value for Germany
R0_rki <- read.csv("RKI.csv")

RKI <- data.frame(beta= rep(NA, length(R0_rki$R)))
RKI$R <- R0_rki$R
RKI$time <- as.Date(0:(length(RKI$R)-1), "2020-03-06")

results$RKI4 <- rep(NA,length(results$dI))
for (i in 37:length(results$dI)){
results$RKI4[i]<- sum(data$dI[i-0:6])/sum(data$dI[i-4:10])
}

RKI$col <- "R RKI 7 days" 

#plot R0
p11 <- ggplot(results, aes(x=as.Date(time, "2020-01-22"), y=exp(Rt),col=place))+geom_line()+
geom_line(data=RKI,aes(y=R, x=time, colour="RKI 7 days"))+
geom_vline(xintercept = as.Date("2020-03-17"), col="purple")+
annotate("text", x = as.Date("2020-02-20"), y = 1, label = "1. Lockdown \n (start)", color="purple")+
geom_vline(xintercept = as.Date("2020-05-6"), col="blue")+
annotate("text", x = as.Date("2020-04-10"), y = 15, label = "1. Lockdown \n (stop)", color="blue")+
geom_vline(xintercept = as.Date("2020-06-16"), col="green")+
annotate("text", x = as.Date("2020-07-12"), y = 15, label = "Corona \n warn app", color="green")+
geom_vline(xintercept = as.Date("2020-10-28"), col="orange")+
annotate("text", x = as.Date("2020-10-1"), y = 5, label = "Lockdown \n light (start)", color="orange")+
geom_vline(xintercept = as.Date("2020-12-13"), col="red")+
annotate("text", x = as.Date("2021-01-10"), y = 5, label = "2. Lockdown \n (start)", color="red")+
geom_vline(xintercept = as.Date("2021-02-15"), col="seagreen")+
annotate("text", x = as.Date("2021-03-12"), y = 10, label = "Alpha variant\n spreads", color="seagreen")+
geom_vline(xintercept = as.Date("2021-06-21"), col="turquoise")+
annotate("text", x = as.Date("2021-05-15"), y = 3, label = "Delta variant\n dominant", color="turquoise")+
geom_line(data = results, aes( y = RKI4, colour ="Incidence method s=4 days "))+
ylab("Value")+ggtitle(expression(paste("Effective Reproduction Number ",R[t])))+xlab("Time ")+ theme_dMod()+scale_x_date(date_breaks = "2 months", date_labels = "%b %y")+
scale_color_dMod()+scale_y_log10()+ theme(legend.title = element_blank(),legend.position = "bottom",plot.title = element_text(hjust = 0.5))+geom_ribbon(data=results,aes(ymin=exp(Rt-sRt),ymax=exp(Rt+sRt)), color = NA,alpha=0.1)

#Total plot
blob <- data.frame()
blob <-results
error1 <- subset(blob,select=c("sI", "sR","sD","sdI","sdR","sdD", "sgamma","stheta","time"))
colnames(error1) <- c("I", "R","D","v[I]", "v[R]", "v[D]", "gamma","theta","time")
smooth1 <- subset(blob, select=c("I", "R","D","dI", "dR", "dD", "gamma","theta","time"))
colnames(smooth1) <- c("I", "R","D","v[I]", "v[R]", "v[D]", "gamma","theta","time")
mysmooth <- gather(smooth1,key="name", value="value", -c("time"))
mysmooth$name <- factor(mysmooth$name, levels=c("I","R","D","v[I]", "v[R]", "v[D]","gamma","theta"))
mysmooth$iteration <- "AKS"
myerror <- gather(error1, key="name", value = "error", -c("time"))
mysmooth$error <- myerror$error
colnames(data)<- c("v[I]", "v[R]", "v[D]","time")
mymeas <- cbind(gather(data,key="name", value="value", -c("time")))
mymeas$name <- factor(mymeas$name, levels=c("v[I]", "v[R]", "v[D]"))
mymeas$iteration <-"Recorded data"

base_size <- 10

cols <- c("Recorded data" = dMod:::dMod_colors[4],
          "AKS" = dMod:::dMod_colors[1])

my_labeller <- as_labeller(c(gamma="gamma",
                             theta="theta",
                             default = label_parsed))

plot_fluxes <- function(smooth,meas) {
  ggplot() + 
    geom_line(data = subset(smooth,name%in%c("v[I]", "v[R]", "v[D]")),
              aes(x=as.Date(time, "2020-01-22"), y=value/1e4, col=iteration)) +
    geom_ribbon(data = subset(smooth,name%in%c("v[I]", "v[R]", "v[D]")),
                aes(x=as.Date(time, "2020-01-22"), ymin = (value-error)/1e4, ymax = (value+error)/1e4),color=NA,alpha=0.3)+
    facet_wrap(~name, scales="free_y", ncol=1, labeller = label_parsed) + 
    ylab(expression(Value~"["~10^4~"]"))+xlab("")+ theme_dMod(base_size = base_size)+
    geom_point(data=meas,
               aes(x=as.Date(time, "2020-01-22"), y=value/1e4, col=iteration),
               size = 0.01, shape=20 
               
    )+
    scale_color_manual(values=cols)+scale_y_continuous(n.breaks = 3)+scale_x_date(date_breaks = "4 months", date_labels = "%b %y")+
    theme(legend.title = element_blank(), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

myfluxes <-plot_fluxes( mysmooth, mymeas)

plot_params <- function(smooth){
  plot_data <- subset(smooth,name%in%c("gamma","theta"))
  
  ggplot() + 
    geom_line(data = plot_data,aes(x=as.Date(time, "2020-01-22"), y=exp(value), col=iteration)) +
    geom_ribbon(data = plot_data,aes(x=as.Date(time, "2020-01-22"), ymin = exp(value-error), ymax = exp(value+error)),color=NA,alpha=0.3)+
    facet_wrap(~name, scales="free_y", ncol=1, labeller = label_parsed) + 
    ylab("Value")+xlab("")+ theme_dMod(base_size = base_size)+
    scale_color_manual(values=cols)+scale_y_log10(n.breaks = 5)+scale_x_date(date_breaks = "4 months", date_labels = "%b %y")+
    theme(legend.title = element_blank(), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

myparams <- plot_params(mysmooth)

pdf(file = "Germany_plot.pdf",width= 12, height=11, useDingbats = FALSE)
P <- ggdraw() + 
  draw_plot(myfluxes, x = 0.0, y = .6, width = .5, height = .4) +
  draw_plot(myparams, x = 0.5, y = .6, width = .5, height = .4)+
  draw_plot(p11, x = 0.0, y = 0, width = 1, height = .6)+
  draw_plot_label(c("A", "B", "C"), x=c(0,0.5,0), y=c(1,1,.6), size=15)

print(P)
dev.off()
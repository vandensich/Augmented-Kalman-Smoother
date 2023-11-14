#' 2022-02-10
#' Script to produce a figure similar to figure 3 of 
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

load(file = "Result_simulations_SEIR.RData")
load( file = "Error_simulations_SEIR.RData")
load( file = "Data_simulations_SEIR.RData")
load( file = "Meas_simulations_SEIR.RData")

mysmooth <- list()
mymeas <- list()
mydata <- list()

for (i in 1:length(dat)){
  mysmooth[[i]] <- gather(result[[i]],key="name", value="value", -c("time", "iteration"))
  for (h in 1:length(mysmooth[[i]]$name)){
    if (mysmooth[[i]][h,]$name=="Rt"){
      mysmooth[[i]][h,]$name<-"R[t]"
    }
    if (mysmooth[[i]][h,]$name=="dI"){
      mysmooth[[i]][h,]$name<-"v[I]"
    }
    if (mysmooth[[i]][h,]$name=="dR"){
      mysmooth[[i]][h,]$name<-"v[R]"
    }
    if (mysmooth[[i]][h,]$name=="dD"){
      mysmooth[[i]][h,]$name<-"v[D]"
    }
  }
  myerror <- gather(error[[i]], key="name", value = "error")
  mysmooth[[i]]$error <- myerror$error
  mymeas[[i]] <- cbind(gather(meas[[i]],key="name", value="value", -c("time")), iteration="Sim. data")
  for (h in 1:length(mymeas[[i]]$name)){
    if (mymeas[[i]][h,]$name=="dI"){
      mymeas[[i]][h,]$name<-"v[I]"
    }
    if (mymeas[[i]][h,]$name=="dR"){
      mymeas[[i]][h,]$name<-"v[R]"
    }
    if (mymeas[[i]][h,]$name=="dD"){
      mymeas[[i]][h,]$name<-"v[D]"
    }
  }
  mydata[[i]] <- cbind(gather(dat[[i]], key="name", value="value",-c("time")), iteration="True value")
  for (h in 1:length(mydata[[i]]$name)){
    if (mydata[[i]][h,]$name=="R0"){
      mydata[[i]][h,]$name<-"R[t]"
    }
    if (mydata[[i]][h,]$name=="dI"){
      mydata[[i]][h,]$name<-"v[I]"
    }
    if (mydata[[i]][h,]$name=="dR"){
      mydata[[i]][h,]$name<-"v[R]"
    }
    if (mydata[[i]][h,]$name=="dD"){
      mydata[[i]][h,]$name<-"v[D]"
    }
  }
  mysmooth[[i]]$name <- factor(mysmooth[[i]]$name, levels=c("I","R","D","v[I]", "v[R]", "v[D]","R[t]","gamma","theta"))
  mymeas[[i]]$name <- factor(mymeas[[i]]$name, levels=c("v[I]", "v[R]", "v[D]"))
  mydata[[i]]$name <- factor(mydata[[i]]$name, levels=c("I","R","D","v[I]", "v[R]", "v[D]","R[t]","gamma","theta"))
  
}

base_size <- 10
cols <- c("True value" = dMod:::dMod_colors[2],
          "Sim. data" = dMod:::dMod_colors[4],
          "AKS" = dMod:::dMod_colors[1])

my_labeller <- as_labeller(c(gamma="gamma",
                             theta="theta",
                             'R[t]'="R[t]"),
                           default = label_parsed)

plot_fluxes <- function(smooth,meas) {
  ggplot() + 
    geom_line(data = subset(smooth,name%in%c("v[I]", "v[R]", "v[D]")),
              aes(x=time, y=value/1e5, col=iteration)) +
    geom_ribbon(data = subset(smooth,name%in%c("v[I]", "v[R]", "v[D]")),
                aes(x=time, ymin = (value-error)/1e5, ymax = (value+error)/1e5),color=NA,alpha=0.3)+
    facet_wrap(~name, scales="free_y", ncol=1, labeller = label_parsed) + 
    ylab(expression(Value~"["~10^5~"]"))+xlab("Time")+ theme_dMod(base_size = base_size)+
    geom_point(data=meas,
               aes(x=time, y=value/1e5, col=iteration),
               size = 0.01, shape=20 
    )+
    scale_color_manual(values=cols)+scale_y_continuous(n.breaks = 3)+
    theme(legend.title = element_blank(), legend.position = "bottom",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

myfluxes <- mapply(plot_fluxes, mysmooth, mymeas, SIMPLIFY = FALSE)

plot_params <- function(smooth, data){
  plot_data <- subset(smooth,name%in%c("R[t]","gamma","theta"))
  blankdata <- data.frame(name=c("R[t]","R[t]","gamma","gamma","theta","theta"),
                          time=c(0,300,0,300,0,300),
                          value=c(0.1,10,0.01,0.1,0.001,0.01))
  ggplot() + 
    geom_line(data = plot_data,aes(x=time, y=exp(value), col=iteration)) +
    geom_ribbon(data = plot_data,aes(x=time, ymin = exp(value-error), ymax = exp(value+error)),color=NA,alpha=0.3)+
    facet_wrap(~name, scales="free_y", ncol=1, labeller = my_labeller) + 
    ylab("Value")+xlab("Time")+ theme_dMod(base_size = base_size)+
    geom_line(data=subset(data,name%in%c("R[t]","gamma","theta")),
              aes(y=value, x=time, col=iteration),linetype="solid")+
    geom_blank(data=blankdata, aes(x=time, y=value)) +
    scale_color_manual(values=cols)+scale_y_log10(n.breaks = 5)+
    theme(legend.title = element_blank(), legend.position = "bottom",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

mybeta <- mapply(plot_params, mysmooth, mydata, SIMPLIFY = F)

plot_states <- function(smooth,data){
  ggplot() + 
    geom_line(data = subset(smooth,name%in%c("I", "R", "D")),
              aes(x=time, y=value/1e6, col=iteration)) +
    geom_ribbon(data = subset(smooth,name%in%c( "I", "R","D")),
                aes(x=time, ymin = (value-error)/1e6, ymax = (value+error)/1e6),color=NA ,alpha=0.3)+
    geom_line(data=subset(data, name%in%c( "I","R","D")),
              aes(y=value/1e6, x=time, col=iteration),linetype="longdash")+
    facet_wrap(~name, scales="free_y", ncol=1) + 
    ylab(expression(Value~"["~10^6~"]"))+xlab("Time")+ theme_dMod(base_size = base_size)+
    scale_color_manual(values=cols)+ scale_y_continuous(n.breaks = 4.5)+
    theme(legend.title = element_blank(), legend.position = "bottom",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

mystates <- mapply( plot_states,mysmooth,mydata, SIMPLIFY = F)

plot_R <- function(smooth, data){
  ggplot() + 
    geom_line(data = subset(smooth,name%in%c("R[t]")), 
              aes(x=time, y=exp(value), col=iteration)) +
    geom_ribbon(data = subset(smooth,name%in%c("R[t]")), 
                aes(x=time, ymin = exp(value-error), ymax = exp(value+error),color=NA),alpha=0.3)+
    #ggtitle("Infected flux")+
    facet_wrap(~name, scales="free_y", nrow=1, labeller = label_parsed) + 
    ylab("")+xlab("Time")+ theme_dMod(base_size = 8)+
    geom_line(data=subset(data, name%in%c("R[t]")),
              aes(y=value, x=time, col=iteration),linetype="solid")+
    scale_color_dMod() +
    theme(legend.title = element_blank(), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

myR <- mapply(plot_R, mysmooth, mydata, SIMPLIFY = F)

pdf(file = "SEIRD_plot.pdf",width= 10, height=9, useDingbats = FALSE)
P <- ggdraw() + 
  draw_plot(mybeta[[1]]+ theme(legend.position = "none"), x = 0.03, y = .5, width = .24, height = .48) +
  draw_plot(myfluxes[[1]]+ theme(legend.position = "none"), x = 0.03, y = 0.05, width = .24, height = .45) +
  draw_plot(mybeta[[2]]+ theme(legend.position = "none"), x = 0.27, y = .5, width = .24, height = .48) +
  draw_plot(myfluxes[[2]]+ theme(legend.position = "none"), x = 0.27, y = 0.05, width = .24, height = .45) +
  draw_plot(mybeta[[4]]+ theme(legend.position = "none"), x = 0.51, y = .5, width = .24, height = .48) +
  draw_plot(myfluxes[[4]]+ theme(legend.position = "none"), x = 0.51, y = 0.05, width = .24, height = .45) +
  draw_plot(mybeta[[5]]+ theme(legend.position = "none"), x = 0.75, y = .5, width = .24, height = .48) +
  draw_plot(myfluxes[[5]]+ theme(legend.position = "none"), x = 0.75, y = 0.05, width = .24, height = .45) +
  draw_plot_label(c("A","B","C","D"), x=c(0.03,.27,0.51,0.75), y=c(1,1,1,1), size=12)+
  draw_label(expression(paste(delta^-1,"= 4 days")), x=0.11, y=.99, size=10,fontface = "plain")+
  draw_label(expression(paste(delta^-1,"= 7 days")), x=0.35, y=.99, size=10,fontface = "plain")+
  draw_label(expression(paste(delta^-1,"= 14 days")), x=0.59, y=.99, size=10,fontface = "plain")+
  draw_label(expression(paste(delta^-1,"= 21 days")), x=0.83, y=.99, size=10,fontface = "plain")+
  draw_label("Parameters", x = 0.015, y= .75, angle = 90, size = 19, fontface = "bold")+
  draw_label("Observations", x = 0.015, y= .25, angle = 90, size = 19, fontface = "bold")
print(P)
dev.off()
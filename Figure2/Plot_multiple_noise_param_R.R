#' 2022-02-10
#' Script to produce a figure similar to figure 2 of 
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

load(file = "Result_simulations_SIRD_noise.RData")
load( file = "Error_simulations_SIRD_noise.RData")
load( file = "Meas_simulations_SIRD_noise.RData")
load(file = "SIR_flux_sim_pred_Jacques_all.RData")

mysmooth <- list()
mymeas <- list()
mydata <- list()
dat <- list()

dat[[1]]<- as.data.frame(prediction[[6]])
dat[[1]]$gamma <- 0.02083
dat[[1]]$theta <- 0.00416
dat[[1]]$"R[t]" <- dat[[1]]$beta/(dat[[1]]$gamma+dat[[1]]$theta)* dat[[1]]$S/8e7

for (i in 1:length(meas)){
  mysmooth[[i]] <- gather(result[[i]],key="name", value="value", -c("time"))
  mysmooth[[i]]$iteration <- "AKS"
  
  myerror <- gather(error[[i]], key="name", value = "error", -c("time"))
  mysmooth[[i]]$error <- myerror$error
  mymeas[[i]] <- cbind(gather(meas[[i]],key="name", value="value", -c("time")))
  mymeas[[i]]$iteration <-"Sim. data"
  
  mydata[[i]] <- cbind(gather(dat[[1]], key="name", value="value",-c("time")))
  mydata[[i]]$iteration <-"True value"
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
    blankdata <- data.frame(name=c("v[I]","v[I]","v[R]","v[R]","v[D]","v[D]"),
                            time=c(0,300,0,300,0,300),
                            value=c(0,6,0,5,0,0.7))
    ggplot() + 
      geom_line(data = subset(smooth,name%in%c("v[I]", "v[R]", "v[D]")),
                aes(x=time, y=value/1e5, col=iteration)) +
      geom_ribbon(data = subset(smooth,name%in%c("v[I]", "v[R]", "v[D]")),
                  aes(x=time, ymin = (value-error)/1e5, ymax = (value+error)/1e5),color=NA,alpha=0.3)+
      geom_blank(data=blankdata, aes(x=time, y=value)) +
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

myparams <- mapply(plot_params, mysmooth, mydata, SIMPLIFY = F)

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

mystates <- mapply(plot_states,mysmooth,mydata, SIMPLIFY = F)

plot_R <- function(smooth, data){
  blankdata <- data.frame(name=c("R[t]", "R[t]"),
                          time=c(0,300),
                          value=c(0.1,10))
  ggplot() + 
    geom_line(data = subset(smooth,name%in%c("R[t]")), 
              aes(x=time, y=exp(value), col=iteration)) +
    geom_ribbon(data = subset(smooth,name%in%c("R[t]")), 
                aes(x=time, ymin = exp(value-error), ymax = exp(value+error),color=NA),alpha=0.3)+
    facet_wrap(~name, scales="free_y", nrow=1, labeller = label_parsed) + 
    ylab("")+xlab("Time")+ theme_dMod(base_size = 10)+
    geom_line(data=subset(data, name%in%c("R[t]")),
              aes(y=value, x=time, col=iteration),linetype="solid")+
    scale_color_dMod() +scale_y_log10(n.breaks = 5)+
    geom_blank(data=blankdata, aes(x=time, y=value)) +
    theme(legend.title = element_blank(), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}
myR <- mapply(plot_R, mysmooth, mydata, SIMPLIFY = F)

pdf(file = "SIRD_noise_plot.pdf",width= 10, height=9, useDingbats = FALSE)
P <- ggdraw() + 
  draw_plot(myparams[[2]]+ theme(legend.position = "none"), x = 0.03, y = .5, width = .24, height = .48) +
  draw_plot(myfluxes[[2]]+ theme(legend.position = "none"), x = 0.03, y = 0.05, width = .24, height = .45) +
  draw_plot(myparams[[3]]+ theme(legend.position = "none"), x = 0.27, y = .5, width = .24, height = .48) +
  draw_plot(myfluxes[[3]]+ theme(legend.position = "none"), x = 0.27, y = 0.05, width = .24, height = .45) +
  draw_plot(myparams[[4]]+ theme(legend.position = "none"), x = 0.51, y = .5, width = .24, height = .48) +
  draw_plot(myfluxes[[4]]+ theme(legend.position = "none"), x = 0.51, y = 0.05, width = .24, height = .45) +
  draw_plot(myparams[[5]]+ theme(legend.position = "none"), x = 0.75, y = .5, width = .24, height = .48) +
  draw_plot(myfluxes[[5]]+ theme(legend.position = "none"), x = 0.75, y = 0.05, width = .24, height = .45) +
  draw_plot_label(c("A","B","C","D"), x=c(0.03,.27,0.51,0.75), y=c(1,1,1,1), size=12)+
  draw_label("5% rel. Error", x=0.11, y=.99, size=8,fontface = "plain")+
  draw_label("10% rel. Error", x=0.35, y=.99, size=8,fontface = "plain")+
  draw_label("20% rel. Error", x=0.59, y=.99, size=8,fontface = "plain")+
  draw_label("30% rel. Error", x=0.83, y=.99, size=8,fontface = "plain")+
  draw_label("Parameters", x = 0.015, y= .75, angle = 90, size = 19, fontface = "bold")+
  draw_label("Observations", x = 0.015, y= .25, angle = 90, size = 19, fontface = "bold")
print(P)
dev.off()
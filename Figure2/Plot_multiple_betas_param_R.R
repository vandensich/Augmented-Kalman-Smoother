#' 2022-02-10
#' Script to produce a figure similar to figure 1 of 
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

load(file = "Result_simulations_SIRD.RData")
load( file = "Error_simulations_SIRD.RData")
load( file = "Meas_simulations_SIRD.RData")
load(file = "SIR_flux_sim_pred_Jacques_all_gamma=002.RData")

mysmooth <- list()
mymeas <- list()
mydata <- list()
dat <- list()

for (i in 1:length(prediction)){
  dat[[i]]<- as.data.frame(prediction[[i]])
  dat[[i]]$gamma <- 0.02
  dat[[i]]$theta <- 0.004
  dat[[i]]$"R[t]" <- dat[[i]]$beta/(dat[[i]]$gamma+dat[[i]]$theta)* dat[[i]]$S/8e7
}

load("SIR_flux_sim_pred_Jacques_all_parameters_gamma=002.RData")
dat[[i+1]]<- as.data.frame(prediction[[1]])
dat[[i+1]]$gamma <- 0.02
dat[[i+1]]$"R[t]" <- dat[[i+1]]$beta/(dat[[i+1]]$gamma+dat[[i+1]]$theta)* dat[[i+1]]$S/8e7

Parameter <- data.frame(Rt = result[[7]]$"R[t]")
Parameter$gamma <- result[[7]]$gamma
Parameter$theta <- result[[7]]$theta
Parameter$time <- result[[7]]$time
Parameter$I <- result[[7]]$I
Parameter$R <- result[[7]]$R
Parameter$D <- result[[7]]$D
Parameter$dI <- result[[7]]$'v[I]'
Parameter$dD <- result[[7]]$'v[D]'
Parameter$dR <- result[[7]]$'v[R]'

for (i in 1:length(meas)){
  mysmooth[[i]] <- gather(result[[i]],key="name", value="value", -c("time"))
  mysmooth[[i]]$iteration <- "AKS"
  
  myerror <- gather(error[[i]], key="name", value = "error", -c("time"))
  mysmooth[[i]]$error <- myerror$error
  mymeas[[i]] <- cbind(gather(meas[[i]],key="name", value="value", -c("time")))
  mymeas[[i]]$iteration <-"Sim. data"
  
  mydata[[i]] <- cbind(gather(dat[[i]], key="name", value="value",-c("time")))
  mydata[[i]]$iteration <-"True value"
  mysmooth[[i]]$name <- factor(mysmooth[[i]]$name, levels=c("S","I","R","D","v[I]", "v[R]", "v[D]","R[t]","beta","gamma","theta"))
  mymeas[[i]]$name <- factor(mymeas[[i]]$name, levels=c("v[I]", "v[R]", "v[D]"))
  mydata[[i]]$name <- factor(mydata[[i]]$name, levels=c("S","I","R","D","v[I]", "v[R]", "v[D]","R[t]","beta","gamma","theta"))
  
}

base_size <- 8

cols <- c("True value" = dMod:::dMod_colors[2],
          "Sim. data" = dMod:::dMod_colors[4],
          "AKS" = dMod:::dMod_colors[1])

my_labeller <- as_labeller(c(gamma="gamma",
                             theta="theta",
                             'R[t]'="R[t]",
                             beta = "beta"),
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

plot_params_mit_beta <- function(smooth, data){
  plot_data <- subset(smooth,name%in%c("R[t]","gamma","theta"))
  blankdata <- data.frame(name=c("R[t]","R[t]","gamma","gamma","theta","theta", "beta", "beta"),
                          time=c(0,300,0,300,0,300,0,300),
                          value=c(0.1,10,0.01,0.1,0.001,0.01,0.01,0.1))
  ggplot() + 
    geom_line(data = plot_data,aes(x=time, y=exp(value), col=iteration)) +
    geom_ribbon(data = plot_data,aes(x=time, ymin = exp(value-error), ymax = exp(value+error)),color=NA,alpha=0.3)+
    facet_wrap(~name, scales="free_y", ncol=1, labeller = my_labeller) + 
    ylab("Value")+xlab("Time")+ theme_dMod(base_size = base_size)+
    geom_line(data=subset(data,name%in%c("R[t]","gamma","theta","beta")),
              aes(y=value, x=time, col=iteration),linetype="solid")+
    geom_blank(data=blankdata, aes(x=time, y=value)) +
    scale_color_manual(values=cols)+scale_y_log10(n.breaks = 5)+
    theme(legend.title = element_blank(), legend.position = "bottom",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}
myparams_mitbeta <- mapply(plot_params_mit_beta, mysmooth, mydata, SIMPLIFY = F)

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
    geom_line(data=subset(data, name%in%c("S","I","R","D")),
              aes(y=value/1e6, x=time, col=iteration),linetype="longdash",size=.8)+
    facet_wrap(~name, scales="free_y", ncol=1) + 
    ylab(expression(Value~"["~10^6~"]"))+xlab("Time")+ theme_dMod(base_size = base_size)+
    scale_color_manual(values=cols)+ scale_y_continuous(n.breaks = 4)+
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
                aes(x=time, ymin = exp(value-error), ymax = exp(value+error)),color=NA,alpha=0.3)+
    facet_wrap(~name, scales="free_y", nrow=1, labeller = label_parsed) + 
    ylab("Value")+xlab("Time")+ theme_dMod(base_size = base_size)+
    geom_line(data=subset(data, name%in%c("R[t]")),
              aes(y=value, x=time, col=iteration),linetype="solid")+
    scale_color_dMod() +scale_y_log10(n.breaks = 5)+
    geom_blank(data=blankdata, aes(x=time, y=value)) +
    theme(legend.title = element_blank(), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}
myR <- mapply(plot_R, mysmooth, mydata, SIMPLIFY = F)

pdf(file = "SIRD_flux_plot.pdf",width= 18/2.54, height=20/2.54, useDingbats = FALSE)
P <- ggdraw() + 
  draw_plot(myfluxes[[7]] + theme(legend.position = "none"), x = 0.666, y = .5, width = .333, height = .48) +
  draw_plot(mystates[[7]] + theme(legend.position = "none"), x = 0.333, y = .5, width = .333, height = .48) +
  draw_plot(myparams_mitbeta[[7]] + theme(legend.position = "none"), x = 0, y = .5, width = .333, height = .48) +
  draw_plot(myparams[[6]] + facet_wrap(~name, scales="free_y", ncol=3, labeller = label_parsed) + theme(legend.position = "none"),
            x = 0, y = 0.3, width = 1, height = .2) + 
  draw_plot(myparams[[4]] + facet_wrap(~name, scales="free_y", ncol=3, labeller = label_parsed) + theme(legend.position = "none"),
            x = 0, y = 0.1, width = 1, height = .2)+
  draw_plot_label(c("A", "B", "C"), x=c(0,0,0), y=c(1,.5,.3), size=12)+
  draw_plot_label(c("Parameters","States", "Observations"), x=c(1/12,5.5/12,8.7/12), y=c(1,1,1), size = 10)
print(P)
dev.off()





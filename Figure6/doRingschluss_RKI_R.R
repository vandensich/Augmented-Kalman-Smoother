#' 2022-02-10
#' Script of the analysis employed to produce the results shown in figure 6 of 
#' https://arxiv.org/abs/2111.09803
#' author: Jacques Hermes <jacques.hermes@fdm.uni-freiburg.de>

# Sets current working directory to the location of this script
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

# Loading Libraries
library(deSolve)
library(trust)
library(parallel)
library(ggplot2)
library(ggthemes)
library(cOde)
library(dMod)
library(cowplot)
library(magrittr)
library(tidyr)
library(fBasics)

# Model Definition - Equations

model_name <- "corona_Ringschluss"
flist <- NULL %>% 
  addReaction("", "I", "Rt*(gamma+theta)*I+dummy") %>%
  addReaction("I", "R", "gamma*I") %>%
  addReaction("I", "D", "theta*I") 

# Model Definition - Observables

observables <- eqnvec(dI = "Rt*(gamma+theta)*I",
                      dR = "gamma *I",
                      dD = "theta*I")

# Model Generation

modelCorona <- odemodel(flist, forcings = c("Rt","theta", "gamma"),
                        fixed=NULL, modelname = paste0("odemodel_", model_name),
                        jacobian = "inz.lsodes", compile = TRUE)
myt <- seq(0,500, 1)

# Parameters obtained from Analysis and Plot generation of Figure 5 

load("Parameter_Bestfit_Germany.RData")
myinitime <- 62
Parameter$time <- Parameter$time - myinitime

myforcings <- subset(gather(Parameter, key="name", value="value", -time), name%in%c("gamma", "theta", "Rt"))

myforcings2 <- gather(Parameter, key="name", value="value", -time)

Parameter$Rt <- rep(NA,length(Parameter$dI))
for (i in 37:length(Parameter$dI)){
  Parameter$Rt[i]<- sum(Parameter$dI[i-0:6])/sum(Parameter$dI[i-4:10])
}

myforcings3 <- subset(gather(Parameter, key="name", value="value", -time), name%in%c("gamma", "theta", "Rt"))

g <- Y(observables, states = c("Rt","I","theta","gamma"),parameters = NULL, compile=TRUE, modelname=paste0("g_",model_name))

# Define parameters, initial parameter trafo and steady-state transformations

constraints <- resolveRecurrence(c(
  dI = as.character(subset(myforcings2, name=="dI" & time==myinitime)$value),
  dR = as.character(subset(myforcings2, name=="dR" & time==myinitime)$value),
  dD = as.character(subset(myforcings2, name=="dD" & time==myinitime)$value),
  I = as.character(subset(myforcings2, name=="I" & time==myinitime)$value),
  R = as.character(subset(myforcings2, name=="R" & time==myinitime)$value),
  D = as.character(subset(myforcings2, name=="D" & time==myinitime)$value)))

innerpars <- unique(c(getParameters(modelCorona), getParameters(g)))
names(innerpars) <- innerpars

trafo <- replaceSymbols(names(observables), observables, innerpars)
trafo <- replaceSymbols(names(constraints), constraints, innerpars)

# Generate condition.grid
covtable <- data.frame(
                       row.names = c("AKS","RKI"))


# Parameter transformation
trafoL <- branch(trafo, table=covtable) 

# Specify prediction functions
tolerances <- 1e-12
p0 <- x <- NULL
for (C in names(trafoL)) {
  p0 <- p0 + P(trafoL[[C]], condition = C)
}

x <- Xs(modelCorona, myforcings, condition = "AKS",
        optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances))+
  Xs(modelCorona, myforcings3, condition = "RKI",
     optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances))

# Generate objective function and initial parameter set
outerpars <- attr(p0, "parameters")
pouter <- structure(rnorm(length(outerpars)), names = outerpars)
prior <- rep(0, length(outerpars)); names(prior) <- outerpars

pouter["dummy"] <- 0

# Loading data which was taken from 
# https://github.com/CSSEGISandData/COVID-19 
# and averaged over 7 days
load("Data_multi_7_SIR_aktuell_August.RData")


data<-subset(data, time >=myinitime)
colnames(data)<-c("v[I]", "v[R]", "v[D]","time")
data <-subset(gather(data, key="name", value="value", -time), name%in%c("v[I]", "v[R]", "v[D]"))


prediction <- (g*x*p0)(myt,pouter)
colnames(prediction[[1]])[2]<-"v[I]"
colnames(prediction[[1]])[3]<-"v[R]"

colnames(prediction[[1]])[4]<-"v[D]"
colnames(prediction[[2]])[2]<-"v[I]"
colnames(prediction[[2]])[3]<-"v[R]"
colnames(prediction[[2]])[4]<-"v[D]"

myIAKS <- wide2long(prediction[[1]])
myIAKS$time<-myIAKS$time+myinitime
myIRKI <-wide2long(prediction[[2]])
myIRKI$time<-myIRKI$time+myinitime
myAKS <- subset(myIAKS, name%in%c("v[I]", "v[R]", "v[D]"))
myRKI <- subset(myIRKI, name%in%c("v[I]", "v[R]", "v[D]"))

myAKS$name <- factor(myAKS$name, levels=c("v[I]", "v[R]", "v[D]"))
myRKI$name <- factor(myRKI$name, levels=c("v[I]", "v[R]", "v[D]"))
data$name <- factor(data$name, levels=c("v[I]", "v[R]", "v[D]"))

cols <- c("ODE trajectory with AKS estimated R" = dMod:::dMod_colors[1],
          "Data" = dMod:::dMod_colors[4],
          "ODE trajectory with Incidence-based s=4 days R" = dMod:::dMod_colors[2])

plot<-ggplot() + 
  geom_line(data=rbind(cbind(myAKS, method="ODE trajectory with AKS estimated R")), aes(x=as.Date(time, "2020-01-22"), y=value, color=method))+
  geom_line(data=rbind(cbind(myRKI, method="ODE trajectory with Incidence-based s=4 days R")), aes(x=as.Date(time, "2020-01-22"), y=value, color=method))+
  geom_point(data=rbind(cbind(data, method="Data")), aes(x=as.Date(time, "2020-01-22"), y=value, color=method),size = 1, shape=20)+
  facet_wrap(~name, scales="free", labeller = label_parsed)  + scale_y_log10()+xlab("Time ")+ theme_dMod(base_size = 10)+scale_x_date(date_breaks = "4 months", date_labels = "%b %y")+
  scale_color_manual(values = cols)+scale_y_log10()+ theme(legend.title = element_blank(),legend.position = "bottom",plot.title = element_text(hjust = 0.5))
ggsave("Ringschluss.pdf",plot,height = 5 , width = 10)
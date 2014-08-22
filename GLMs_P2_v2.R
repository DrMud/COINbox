
---
  #  title: "GLM_A&C"
  #author: "Rafael S. de Souza, ..."
  #date: "27 de junho de 2014"
  #output: pdf_document
  ---
  
  
  
  
  # Required libraries
  
require(MASS)
require(grDevices)
require(plyr)
require(msme)
require(ggplot2)
#require(corrplot)
#require(corrgram)
#require(logistf)
require(boot)
require(ggthemes)
library(scales)
require(MCMCpack)
require(arm)
library(pROC)
library(grid)
require(CosmoGLM)


#Reading the data 

#Set working directory to the data folder (replace by your own directory)

data_path<-"/Users/killedar/Documents/Research/projects/CosmoStatsIAA/GLM/PopIIIstars/"
#data_path<-"/Users/rafael/Dropbox/artigos/Meusartigos/IAA-WGC/GLMs/Simulation/data/"



Biffi_data<-read.table(file=paste(data_path,"Biffi2014.csv",sep=""),
                       header=TRUE)

Biffi_data_original<-Biffi_data
# Problem 1: sSFR, fstar, fgas to form an explanatory variable. 
# Z is the response variable


#Transforming  into a binary problem. SF formation activity or not 

Biffi_data$sSFR  <- Biffi_data$SFR/Biffi_data$Mdm
Biffi_data$fstar <- Biffi_data$Mstar/Biffi_data$Mdm
Biffi_data$fgas  <- Biffi_data$Mgas/Biffi_data$Mdm
Biffi_data$Q2    <- Biffi_data$sSFR*Biffi_data$fstar
Biffi_data$StarToGas    <- Biffi_data$fstar/Biffi_data$fgas
Biffi_data$lfg   <- log(Biffi_data$fgas,10)
#Biffi_data$Q3    <- Biffi_data$fstar*Biffi_data$lfg
                        
Biffi_data$Z[which(Biffi_data$Z<1e-4)]<-0
#Biffi_data$Z[which(Biffi_data$Z>=1e-4 & Biffi_data$Z!=1)]<-0
Biffi_data$Z[which(Biffi_data$Z>=1e-4)]<-1
Biffi_data$Z <-as.numeric(Biffi_data$Z)
Biffi_data$Zf<-as.factor(Biffi_data$Z)

#Transforming  variable into numeric (required by GLM packages) 

Biffi_data$sSFR <-as.numeric(Biffi_data$sSFR)
Biffi_data$fstar<-as.numeric(Biffi_data$fstar)
Biffi_data$fgas <-as.numeric(Biffi_data$fgas)
Biffi_data$Q2   <-as.numeric(Biffi_data$Q2)
Biffi_data$StarToGas   <-as.numeric(Biffi_data$StarToGas)
Biffi_data$Mgas <-as.numeric(Biffi_data$Mgas)
Biffi_data$SFR  <-as.numeric(Biffi_data$SFR)
Biffi_data$lfg  <-as.numeric(Biffi_data$lfg)

# Define ranges for later use

fstar_range <- seq(0,1.2e-13,3e-15)
fgas_range  <- seq(0.01,0.14,0.003)
StG_range <- seq(0,1.8e-12,3e-14)

#------------------------------------------------------------------------

# Scatterplot of Metallicity vs 2 of the 3 underlying explanatory variables

scat <- ggplot(data=Biffi_data,aes(x=fstar*1e14,y=fgas))
scat <- scat + xlab("10^14*fstar") + ylab("fgas")
scat <- scat + geom_point(size=3,aes(color=StarToGas*1e13,shape=Zf))
scat <- scat + theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  theme(plot.title = element_text(hjust=0.5),axis.title.y=element_text(vjust=0.75),
        axis.title.x=element_text(vjust=-0.25),text = element_text(size=20),
        legend.key.size=unit(2.0,"cm"), legend.key.height=unit(0.5,"cm"),
        legend.text=element_text(size=16),legend.title=element_text(size=15))+
  #scale_colour_gradient_tableau(name="")
  scale_colour_gradient_tableau()
ggsave("Scatter_P2_Z1.pdf", width=7, height=6)



scat <- ggplot(data=Biffi_data,aes(x=fstar*1e14,y=fgas))
scat <- scat + xlab("10^14*fstar") + ylab("fgas")
scat <- scat + geom_point(size=3,aes(color=StarToGas,shape=Zf))
scat <- scat + theme_economist_white(gray_bg = F, base_size = 11, base_family = "sans")+
  theme(plot.title = element_text(hjust=0.5),axis.title.y=element_text(vjust=0.75),
        axis.title.x=element_text(vjust=-0.25),text = element_text(size=20),
        legend.key.size=unit(2.0,"cm"), legend.key.height=unit(0.5,"cm"),
        legend.text=element_text(size=16),legend.title=element_text(size=15))+
  scale_colour_gradient_tableau()
  scale_y_continuous(expand=c(0.1,0),trans = 'log10',
                   breaks=trans_breaks("log10",function(x) 10^x),
                   labels=trans_format("log10",math_format(10^.x)))
#scale_x_continuous(trans = 'log10',
#                   breaks=trans_breaks("log10",function(x) 10^x),
#                   labels=trans_format("log10",math_format(10^.x)))
ggsave("Scatter_P2_Z2.pdf", width=7, height=6)



scat <- ggplot(data=Biffi_data,aes(x=fstar,y=fgas))
scat <- scat + xlab("fstar") + ylab("fgas")
scat <- scat + geom_point(size=3,aes(color=Q2,shape=Zf))
scat <- scat + theme_economist_white(gray_bg = F, base_size = 10, base_family = "sans")+
  theme(plot.title = element_text(hjust=0.5),axis.title.y=element_text(vjust=0.75),
        axis.title.x=element_text(vjust=-0.25),text = element_text(size=20),
        legend.key.size=unit(2.0,"cm"), legend.key.height=unit(0.5,"cm"),
        legend.text=element_text(size=16),legend.title=element_text(size=15))+
  scale_colour_gradient_tableau()+
ggsave("Scatter_P2_Z3.pdf", width=8, height=7)


# Order of importance: fstar, fgas, sSFR, Q2, Q1


#------------------------------------------------------------------------

#Frequentist Cloglog link
#Fglm.clog <- glm(Z~scale(sSFR)+scale(fstar)+scale(fgas),family=binomial(link="cloglog"),
#           data = Biffi_data)

#plotProb(Fglm.clog)+ylab("Predicted Probabilities for host Pop II/I stars ")+xlab("")+
#  ggtitle("Frequentist Logistic Regression")+coord_cartesian(ylim = c(0,1.05))


#Bayesian Cloglog link
#Bglm.clog <- bayesglm(Z~scale(sSFR)+scale(fstar)+scale(fgas),family=binomial(link="cloglog"),scaled=TRUE,
#           data = Biffi_data,print.unnormalized.log.posterior=T)


#------------------------------------------------------#
#Logit link

#Frequentist
#Fglm.logit <- glm(Z~fstar+StarToGas+Q2+fgas,family=binomial(link="logit"),data = Biffi_data)
Fglm.logit <- glm(Z~StarToGas,family=binomial(link="logit"),data = Biffi_data)
summary(Fglm.logit)
stepres <- step(Fglm.logit)
summary(stepres)
Z_pred_Flog   <- predict(Fglm.logit, list(StarToGas=StG_range), type="response")


#GLM for binomial response using msme from Andrew Robinson & Joseph Hilbe

#Bayesian (Cauchy prior is default?)
#Bglm.logit <- bayesglm(Z~fstar,family=binomial(link="logit"),data = Biffi_data)
Bglm.logit <- bayesglm(Z~StarToGas,family=binomial(link="logit"),data = Biffi_data)
#Bglm.logit <- bayesglm(Z~fstar+sSFR+Q2+fgas+SFR+Mstar+Mgas,family=binomial(link="logit"),data = Biffi_data)
summary(Bglm.logit)
stepres <- step(Bglm.logit)
summary(stepres)
#stargazer(type="latex",style="mnras",summary=TRUE,out="coeff.tex")
Z_pred_Blog <- predict(Bglm.logit, list(StarToGas=StG_range), type="response")


pdf("PredPr_P2_Blog.pdf",width=8,height=7)
plot_Prob(Bglm.logit)+ggtitle("Bayesian Logistic Regression")+ylab("Predicted Probabilities for Z")+xlab("")+
  coord_cartesian(ylim = c(0,1.05))
with(Bglm.logit, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
dev.off()

#pp <- ggplot(data=Biffi_data, aes(x=fstar,y=Z)) + xlab("fstar") + ylab("Z")
#pp <- pp + geom_point(size=3,aes(shape=Zf))
#pp <- pp + geom_line(df2,aes(x=fstar_range,y=Z_pred_Blog),colour="red",size=1.5)
#pp
#ggsave("PredPr_P2_fstar2.pdf",width=6,height=5)

pdf("PredPr_P2_Blog_StG.pdf",width=5,height=4)
plot(Z~StarToGas, data=Biffi_data, pch=c(16,17)[as.factor(Z)], xlab="Stellar/Gas", ylab="Z", cex=1.3)
lines(StG_range,Z_pred_Blog,col="red",lwd=2)
dev.off()



#------------------------------------------------------#
# PROBIT link

#Frequentist
Fglm.probit <- glm(Z~StarToGas,family=binomial(link="probit"),data = Biffi_data)
summary(Fglm.probit)
stepres <- step(Fglm.probit)
summary(stepres)
Z_pred_Fprob   <- predict(Fglm.probit, list(StarToGas=StG_range), type="response")

#Bayesian (Cauchy prior is default?)
#Bglm.probit <- bayesglm(Z~fstar+fgas+Q2+StarToGas+lfg,family=binomial(link="probit"),data = Biffi_data)
Bglm.probit <- bayesglm(Z~StarToGas,family=binomial(link="probit"),data = Biffi_data)
summary(Bglm.probit)
stepres <- step(Bglm.probit)
summary(stepres)
Z_pred_Bprob   <- predict(Bglm.probit, list(StarToGas=StG_range), type="response")

pdf("PredPr_P2_Bprob_StG.pdf",width=5,height=4)
plot(Z~StarToGas, data=Biffi_data, pch=c(16,17)[as.factor(Z)], xlab="Star/Gas", ylab="Z", cex=1.3)
lines(StG_range,Z_pred_Bprob,col="red",lwd=2)
dev.off()

pdf("PredPr_P2_StG.pdf",width=4,height=4)
plot(Z~StarToGas, data=Biffi_data, pch=c(16,17)[as.factor(Z)], xlab="Star/Gas", ylab="Z", cex=1.0)
lines(StG_range,Z_pred_Bprob,col="red",lwd=1.8)
lines(StG_range,Z_pred_Blog, col="blue",lwd=1.8)
lines(StG_range,Z_pred_Fprob,col="green",lwd=0.7)
lines(StG_range,Z_pred_Flog,col="orange",lwd=0.7)
dev.off()


#------------------------------------------------------#

# ROC Curve

require(caret)

folds <- createFolds(Biffi_data$Zf, k=10)
AUC<-c()
for(i in 1:10){
  training <- Biffi_data[-folds[[i]], ]
  testing <- Biffi_data[folds[[i]], ]
  glm.roc <- bayesglm(Zf~fstar*fgas,family=binomial(link="logit"),data = training)
  ROCF<- data.frame(True=training$Zf,predicted=predict(glm.roc, newdata=training,type = "response"))
  F1 <-roc(ROCF$True,ROCF$predicted)
  ROCF2<- data.frame(True=testing$Zf,predicted=predict(glm.roc, newdata=testing,type = "response"))
  F2 <-roc(ROCF2$True,ROCF2$predicted)
  AUC<-append(AUC,F2$auc) 
}
AUC
mean(AUC)


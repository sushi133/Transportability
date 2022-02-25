#packages
library(dplyr)
library(plyr)
library(survival)
library(timereg)
library(aftgee)
library(table1)
library(addhazard)

# R code for reading in Dan Sargent's colon cancer dataset
colon <- read.csv("/Users/sushi5824907/Desktop/RA 2021/Debashis/colon040909.csv",header=F,skip=2)
names(colon) <- c("study","tx","clintx","age","stage","gender","daterand","recurstat","recurtime","dfsstat","dfstime","dstat","dtime")
#subset stage III 
colon0<-colon[which(colon$stage==3),]
#number of subjects in each treatment group
colon0 %>% 
  group_by(tx)%>%
  dplyr::summarise(n = n())
#convert to character
colon0$tx <- as.factor(colon0$tx)
colon0$study <- as.factor(colon0$study)
colon0$gender <- as.factor(colon0$gender)
#convert time to year
colon0$dyear <- colon0$dtime/365
#5 years of follow-up is considered in the analysis
colon0$dyear5<-ifelse(colon0$dyear<=5,colon0$dyear,5)
colon0<-colon0[colon0$dyear5 > 0,]
colon0$dstat<-ifelse(colon0$dyear5<colon0$dyear & colon0$dstat==1,0,colon0$dstat)
#subject ID
colon0$ID <- seq_along(colon0[,1])

#covariates summary table by study
table1(~ age+gender | study, data=colon0)
#KM plots
par(mfrow=c(2,5),mar = c(2.5, 2.5, 2.5, 2.5))
for (i in 1:10) {
  plot(survfit(Surv(dyear5,dstat)~tx,data=colon0,subset=(colon0$study == levels(colon0$study)[i])),lty=1:2,col=c('black','red'),cex.axis=1.5)
  tmp <- survdiff(Surv(dyear5,dstat)~tx,data=colon0,subset=(colon0$study == levels(colon0$study)[i]))
  title(main=paste(levels(colon0$study)[i],"\nLog-rank \nstatistic = ",round(tmp$chisq,2)),line=-20,cex.main = 2)
}

#one-trial-at-a-time transport

#Cox Proportional Hazards Model
tables <- NULL
for (k in 1:10) {
  pred <- numeric(10)
  true <- numeric(10)
  study <- numeric(10)
  fit1<-coxph(Surv(dyear5,dstat)~tx+age+gender,data=colon0[colon0$study == levels(colon0$study)[k],])
  for (i in 1:10) {
    c2<-colon0[colon0$study == levels(colon0$study)[i],]
    c2$hr<-1-survfit(fit1,newdata=c2)$surv[nrow(survfit(fit1,newdata=c2)$surv),]
    c3<-aggregate(c2$hr, list(c2$tx), FUN=mean)
    pred[i]<-c3[2,2]-c3[1,2]
    fit2<-coxph(Surv(dyear5,dstat)~tx+age+gender,data=colon0[colon0$study == levels(colon0$study)[i],])
    c2$hr<-1-survfit(fit2,newdata=c2)$surv[nrow(survfit(fit2,newdata=c2)$surv),]
    c3<-aggregate(c2$hr, list(c2$tx), FUN=mean)
    true[i]<-c3[2,2]-c3[1,2]
  }
  tables <- rbind(tables, cbind(pred[-k],true[-k],k,study[-k]))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("HD_pred", "HD_true",'model','target')
tables$PE<-(tables$HD_pred-tables$HD_true)^2
#study level PE
tables %>%
  group_by(model) %>%
  dplyr::summarise(mean = mean(PE))
#overall PE
tables %>%
  dplyr::summarise(mean = mean(PE))

#Additive Hazards Regression Model
tables <- NULL
for (k in 1:10) {
  pred <- numeric(10)
  true <- numeric(10)
  study <- numeric(10)
  colon0$dyear5_ <- colon0$dyear5 + runif(dim(colon0)[1],0,1)*1e-8
  fit1<-ah(Surv(dyear5_,dstat)~tx+age+gender,data=colon0[colon0$study == levels(colon0$study)[k],],robust=F,ties=F)
  for (i in 1:10) {
    c2<-colon0[colon0$study == levels(colon0$study)[i],]
    c2$hd<-fit1$coef[1]*I(c2$tx==2)+fit1$coef[2]*c2$age+fit1$coef[3]*I(c2$gender==1)
    c3<-aggregate(c2$hd, list(c2$tx), FUN=mean)
    pred[i]<-c3[2,2]-c3[1,2]
    fit2<-ah(Surv(dyear5_,dstat)~tx+age+gender,data=colon0[colon0$study == levels(colon0$study)[i],],robust=F,ties=F)
    true[i]<-fit2$coef[1]
  }
  tables <- rbind(tables, cbind(pred[-k],true[-k],k,study[-k]))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("HD_pred", "HD_true",'model','target')
tables$PE<-(tables$HD_pred-tables$HD_true)^2
#study level PE
tables %>%
  group_by(model) %>%
  dplyr::summarise(mean = mean(PE))
#overall PE
tables %>%
  dplyr::summarise(mean = mean(PE))

#AFT Model-lognormal
tables <- NULL
for (k in 1:10) {
  pred <- numeric(10)
  true <- numeric(10)
  study <- numeric(10)
  fit1<-survreg(Surv(dyear5,dstat)~tx+age+gender,dist="lognormal",data=colon0,subset=(colon0$study == levels(colon0$study)[k]))
  for (i in 1:10) {
    c2<-colon0[colon0$study == levels(colon0$study)[i],]
    pre<-predict(fit1,newdata=c2)
    c3<-aggregate(pre, list(c2$tx), FUN=mean)
    pred[i]<-c3[2,2]-c3[1,2]

    fit2<-survreg(Surv(dyear5,dstat)~tx+age+gender,dist="lognormal",data=colon0,subset=(colon0$study == levels(colon0$study)[i]))
    pre<-predict(fit2)
    c3<-aggregate(pre, list(c2$tx), FUN=mean)
    true[i]<-c3[2,2]-c3[1,2]
    study[i]<-i
  }
  tables <- rbind(tables, cbind(pred[-k],true[-k],k,study[-k]))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("TTD_pred", "TTD_true",'model','target')
tables$PE<-(tables$TTD_pred-tables$TTD_true)^2
#study level PE
tables %>%
  group_by(model) %>%
  dplyr::summarise(mean = mean(PE))
#overall PE
tables %>%
  dplyr::summarise(mean = mean(PE))

#AFT Model-weibull
tables <- NULL
for (k in 1:10) {
  pred <- numeric(10)
  true <- numeric(10)
  study <- numeric(10)
  fit1<-survreg(Surv(dyear5,dstat)~tx+age+gender,data=colon0,subset=(colon0$study == levels(colon0$study)[k]))
  for (i in 1:10) {
    c2<-colon0[colon0$study == levels(colon0$study)[i],]
    pre<-predict(fit1,newdata=c2)
    c3<-aggregate(pre, list(c2$tx), FUN=mean)
    pred[i]<-c3[2,2]-c3[1,2]
    
    fit2<-survreg(Surv(dyear5,dstat)~tx+age+gender,data=colon0,subset=(colon0$study == levels(colon0$study)[i]))
    pre<-predict(fit2)
    c3<-aggregate(pre, list(c2$tx), FUN=mean)
    true[i]<-c3[2,2]-c3[1,2]
    study[i]<-i
  }
  tables <- rbind(tables, cbind(pred[-k],true[-k],k,study[-k]))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("TTD_pred", "TTD_true",'model','target')
tables$PE<-(tables$TTD_pred-tables$TTD_true)^2
#study level PE
tables %>%
  group_by(model) %>%
  dplyr::summarise(mean = mean(PE))
#overall PE
tables %>%
  dplyr::summarise(mean = mean(PE))

#AFT Model-semiparametric
tables <- NULL
for (k in 1:10) {
  pred <- numeric(10)
  true <- numeric(10)
  study <- numeric(10)
  fit1<-aftgee(Surv(dyear5,dstat)~tx+age+gender,id = ID,corstr = "ind",B=0,data=colon0,subset=(colon0$study == levels(colon0$study)[k]))
  for (i in 1:10) {
    c2<-colon0[colon0$study == levels(colon0$study)[i],]
    xb<-coef(fit1)[1]+coef(fit1)[2]*I(c2$tx==2)+coef(fit1)[3]*c2$age+coef(fit1)[4]*I(c2$gender==1)
    c2$ettf<-exp(xb)
    c3<-aggregate(c2$ettf, list(c2$tx), FUN=mean)
    pred[i]<-c3[2,2]-c3[1,2]
    
    fit2<-aftgee(Surv(dyear5,dstat)~tx+age+gender,id = ID,corstr = "ind",B=0,data=colon0,subset=(colon0$study == levels(colon0$study)[i]))
    xb<-coef(fit2)[1]+coef(fit2)[2]*I(c2$tx==2)+coef(fit2)[3]*c2$age+coef(fit2)[4]*I(c2$gender==1)
    c2$ettf<-exp(xb)
    c3<-aggregate(c2$ettf, list(c2$tx), FUN=mean)
    true[i]<-c3[2,2]-c3[1,2]
    study[i]<-i
  }
  tables <- rbind(tables, cbind(pred[-k],true[-k],k,study[-k]))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("TTD_pred", "TTD_true",'model','target')
tables$PE<-(tables$TTD_pred-tables$TTD_true)^2
#study level PE
tables %>%
  group_by(model) %>%
  dplyr::summarise(mean = mean(PE))
#overall PE
tables %>%
  dplyr::summarise(mean = mean(PE))

#left-one-trial-out transport

#Cox Proportional Hazards Model
tables <- NULL
for (k in 1:10) {
    fit1<-coxph(Surv(dyear5,dstat)~tx+age+gender,data=colon0[colon0$study != levels(colon0$study)[k],])
    c2<-colon0[colon0$study == levels(colon0$study)[k],]
    c2$hr<-1-survfit(fit1,newdata=c2)$surv[nrow(survfit(fit1,newdata=c2)$surv),]
    c3<-aggregate(c2$hr, list(c2$tx), FUN=mean)
    pred<-c3[2,2]-c3[1,2]
    fit2<-coxph(Surv(dyear5,dstat)~tx+age+gender,data=colon0[colon0$study == levels(colon0$study)[k],])
    c2$hr<-1-survfit(fit2,newdata=c2)$surv[nrow(survfit(fit2,newdata=c2)$surv),]
    c3<-aggregate(c2$hr, list(c2$tx), FUN=mean)
    true<-c3[2,2]-c3[1,2]
    tables <- rbind(tables, cbind(pred,true,-k,k))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("HD_pred", "HD_true",'model','target')
tables$PE<-(tables$HD_pred-tables$HD_true)^2
#study level PE
tables %>%
  group_by(model) %>%
  dplyr::summarise(mean = mean(PE))%>%
  arrange(desc(model))
#overall PE
tables %>%
  dplyr::summarise(mean = mean(PE))

#Additive Hazards Regression Model
tables <- NULL
for (k in 1:10) {
  fit1<-ah(Surv(dyear5_,dstat)~tx+age+gender,data=colon0[colon0$study != levels(colon0$study)[k],],robust=F,ties=F)
  c2<-colon0[colon0$study == levels(colon0$study)[k],]
  c2$hd<-fit1$coef[1]*I(c2$tx==2)+fit1$coef[2]*c2$age+fit1$coef[3]*I(c2$gender==1)
  c3<-aggregate(c2$hd, list(c2$tx), FUN=mean)
  pred<-c3[2,2]-c3[1,2]
  fit2<-ah(Surv(dyear5_,dstat)~tx+age+gender,data=colon0[colon0$study == levels(colon0$study)[k],],robust=F,ties=F)
  true<-fit2$coef[1]
  tables <- rbind(tables, cbind(pred,true,-k,k))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("HD_pred", "HD_true",'model','target')
tables$PE<-(tables$HD_pred-tables$HD_true)^2
#study level PE
tables %>%
  group_by(model) %>%
  dplyr::summarise(mean = mean(PE)) %>%
  arrange(desc(model))
#overall PE
tables %>%
  dplyr::summarise(mean = mean(PE))


#AFT Model-lognormal
tables <- NULL
for (k in 1:10) {
  fit1<-survreg(Surv(dyear5,dstat)~tx+age+gender,dist="lognormal",data=colon0,subset=(colon0$study != levels(colon0$study)[k]))
  c2<-colon0[colon0$study == levels(colon0$study)[k],]
  pre<-predict(fit1,newdata=c2)
  c3<-aggregate(pre, list(c2$tx), FUN=mean)
  pred<-c3[2,2]-c3[1,2]
  
  fit2<-survreg(Surv(dyear5,dstat)~tx+age+gender,dist="lognormal",data=colon0,subset=(colon0$study == levels(colon0$study)[k]))
  pre<-predict(fit2)
  c3<-aggregate(pre, list(c2$tx), FUN=mean)
  true<-c3[2,2]-c3[1,2]
    
  tables <- rbind(tables, cbind(pred,true,-k,k))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("TTD_pred", "TTD_true",'model','target')
tables$PE<-(tables$TTD_pred-tables$TTD_true)^2
#study level PE
tables %>%
  group_by(model) %>%
  dplyr::summarise(mean = mean(PE)) %>%
  arrange(desc(model))
#overall PE
tables %>%
  dplyr::summarise(mean = mean(PE))

#AFT Model-weibull
tables <- NULL
for (k in 1:10) {
  fit1<-survreg(Surv(dyear5,dstat)~tx+age+gender,data=colon0,subset=(colon0$study != levels(colon0$study)[k]))
  c2<-colon0[colon0$study == levels(colon0$study)[k],]
  pre<-predict(fit1,newdata=c2)
  c3<-aggregate(pre, list(c2$tx), FUN=mean)
  pred<-c3[2,2]-c3[1,2]
  
  fit2<-survreg(Surv(dyear5,dstat)~tx+age+gender,data=colon0,subset=(colon0$study == levels(colon0$study)[k]))
  pre<-predict(fit2)
  c3<-aggregate(pre, list(c2$tx), FUN=mean)
  true<-c3[2,2]-c3[1,2]
  
  tables <- rbind(tables, cbind(pred,true,-k,k))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("TTD_pred", "TTD_true",'model','target')
tables$PE<-(tables$TTD_pred-tables$TTD_true)^2
#study level PE
tables %>%
  group_by(model) %>%
  dplyr::summarise(mean = mean(PE)) %>%
  arrange(desc(model))
#overall PE
tables %>%
  dplyr::summarise(mean = mean(PE))

#AFT Model-semiparametric
tables <- NULL
for (k in 1:10) {
  fit1<-aftgee(Surv(dyear5,dstat)~tx+age+gender,id = ID,corstr = "ind",B=0,data=colon0,subset=(colon0$study != levels(colon0$study)[k]))
  c2<-colon0[colon0$study == levels(colon0$study)[k],]
  xb<-coef(fit1)[1]+coef(fit1)[2]*I(c2$tx==2)+coef(fit1)[3]*c2$age+coef(fit1)[4]*I(c2$gender==1)
  c2$ettf<-exp(xb)
  c3<-aggregate(c2$ettf, list(c2$tx), FUN=mean)
  pred<-c3[2,2]-c3[1,2]
    
  fit2<-aftgee(Surv(dyear5,dstat)~tx+age+gender,id = ID,corstr = "ind",B=0,data=colon0,subset=(colon0$study == levels(colon0$study)[k]))
  xb<-coef(fit2)[1]+coef(fit2)[2]*I(c2$tx==2)+coef(fit2)[3]*c2$age+coef(fit2)[4]*I(c2$gender==1)
  c2$ettf<-exp(xb)
  c3<-aggregate(c2$ettf, list(c2$tx), FUN=mean)
  true<-c3[2,2]-c3[1,2]
  
  tables <- rbind(tables, cbind(pred,true,-k,k))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("TTD_pred", "TTD_true",'model','target')
tables$PE<-(tables$TTD_pred-tables$TTD_true)^2
#study level PE
tables %>%
  group_by(model) %>%
  dplyr::summarise(mean = mean(PE)) %>%
  arrange(desc(model))
#overall PE
tables %>%
  dplyr::summarise(mean = mean(PE))

#meta-analysis

#Cox Proportional Hazards Model
tables <- NULL
for (k in 1:10) {
  fit1<-coxph(Surv(dyear5,dstat)~tx+age+gender,data=colon0[colon0$study != levels(colon0$study)[k],])
  c2<-colon0[colon0$study != levels(colon0$study)[k],]
  c2$hr<-1-survfit(fit1,newdata=c2)$surv[nrow(survfit(fit1,newdata=c2)$surv),]
  c3<-aggregate(c2$hr, list(c2$tx), FUN=mean)
  pred<-c3[2,2]-c3[1,2]
  fit2<-coxph(Surv(dyear5,dstat)~tx+age+gender,data=colon0[colon0$study == levels(colon0$study)[k],])
  c2<-colon0[colon0$study == levels(colon0$study)[k],]
  c2$hr<-1-survfit(fit2,newdata=c2)$surv[nrow(survfit(fit2,newdata=c2)$surv),]
  c3<-aggregate(c2$hr, list(c2$tx), FUN=mean)
  true<-c3[2,2]-c3[1,2]
  tables <- rbind(tables, cbind(pred,true,-k,k))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("HD_pred", "HD_true",'model','target')
tables$PE<-(tables$HD_pred-tables$HD_true)^2
#study level PE
tables %>%
  group_by(model) %>%
  dplyr::summarise(mean = mean(PE)) %>%
  arrange(desc(model))
#overall PE
tables %>%
  dplyr::summarise(mean = mean(PE))


#Additive Hazards Regression Model
tables <- NULL
for (k in 1:10) {
  fit1<-ah(Surv(dyear5_,dstat)~tx+age+gender,data=colon0[colon0$study != levels(colon0$study)[k],],robust=F,ties=F)
  pred<-fit1$coef[1]
  fit2<-ah(Surv(dyear5_,dstat)~tx+age+gender,data=colon0[colon0$study == levels(colon0$study)[k],],robust=F,ties=F)
  true<-fit2$coef[1]
  tables <- rbind(tables, cbind(pred,true,-k,k))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("HD_pred", "HD_true",'model','target')
tables$PE<-(tables$HD_pred-tables$HD_true)^2
#study level PE
tables %>%
  group_by(model) %>%
  dplyr::summarise(mean = mean(PE)) %>%
  arrange(desc(model))
#overall PE
tables %>%
  dplyr::summarise(mean = mean(PE))

#AFT Model-lognormal
tables <- NULL
for (k in 1:10) {
  fit1<-survreg(Surv(dyear5,dstat)~tx+age+gender+study,dist="lognormal",data=colon0,subset=(colon0$study != levels(colon0$study)[k]))
  c2<-colon0[colon0$study != levels(colon0$study)[k],]
  pre<-predict(fit1)
  c3<-aggregate(pre, list(c2$tx), FUN=mean)
  pred<-c3[2,2]-c3[1,2]
  
  fit2<-survreg(Surv(dyear5,dstat)~tx+age+gender,dist="lognormal",data=colon0,subset=(colon0$study == levels(colon0$study)[k]))
  c2<-colon0[colon0$study == levels(colon0$study)[k],]
  pre<-predict(fit2)
  c3<-aggregate(pre, list(c2$tx), FUN=mean)
  true<-c3[2,2]-c3[1,2]
  
  tables <- rbind(tables, cbind(pred,true,-k,k))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("TTD_pred", "TTD_true",'model','target')
tables$PE<-(tables$TTD_pred-tables$TTD_true)^2
#study level PE
tables %>%
  group_by(model) %>%
  dplyr::summarise(mean = mean(PE)) %>%
  arrange(desc(model))
#overall PE
tables %>%
  dplyr::summarise(mean = mean(PE))

#AFT Model-weibull
tables <- NULL
for (k in 1:10) {
  fit1<-survreg(Surv(dyear5,dstat)~tx+age+gender+study,data=colon0,subset=(colon0$study != levels(colon0$study)[k]))
  c2<-colon0[colon0$study != levels(colon0$study)[k],]
  pre<-predict(fit1)
  c3<-aggregate(pre, list(c2$tx), FUN=mean)
  pred<-c3[2,2]-c3[1,2]
  
  fit2<-survreg(Surv(dyear5,dstat)~tx+age+gender,data=colon0,subset=(colon0$study == levels(colon0$study)[k]))
  c2<-colon0[colon0$study == levels(colon0$study)[k],]
  pre<-predict(fit2)
  c3<-aggregate(pre, list(c2$tx), FUN=mean)
  true<-c3[2,2]-c3[1,2]
  
  tables <- rbind(tables, cbind(pred,true,-k,k))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("TTD_pred", "TTD_true",'model','target')
tables$PE<-(tables$TTD_pred-tables$TTD_true)^2
#study level PE
tables %>%
  group_by(model) %>%
  dplyr::summarise(mean = mean(PE)) %>%
  arrange(desc(model))
#overall PE
tables %>%
  dplyr::summarise(mean = mean(PE))

#AFT Model-semiparametric
tables <- NULL
for (k in 1:10) {
  fit1<-aftgee(Surv(dyear5,dstat)~tx+age+gender,id = study,corstr = "ind",B=0,data=colon0,subset=(colon0$study != levels(colon0$study)[k]))
  c2<-colon0[colon0$study != levels(colon0$study)[k],]
  pre<-unlist(predict(fit1))
  c3<-aggregate(pre, list(c2$tx), FUN=mean)
  pred<-c3[2,2]-c3[1,2]
  
  fit2<-aftgee(Surv(dyear5,dstat)~tx+age+gender,id = ID,corstr = "ind",B=0,data=colon0,subset=(colon0$study == levels(colon0$study)[k]))
  c2<-colon0[colon0$study == levels(colon0$study)[k],]
  pre<-unlist(predict(fit2))
  c3<-aggregate(pre, list(c2$tx), FUN=mean)
  true<-c3[2,2]-c3[1,2]
  
  tables <- rbind(tables, cbind(pred,true,-k,k))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("TTD_pred", "TTD_true",'model','target')
tables$PE<-(tables$TTD_pred-tables$TTD_true)^2
#study level PE
tables %>%
  group_by(model) %>%
  dplyr::summarise(mean = mean(PE)) %>%
  arrange(desc(model))
#overall PE
tables %>%
  dplyr::summarise(mean = mean(PE))
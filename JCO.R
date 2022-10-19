#packages
library(dplyr)
library(plyr)
library(survival)
library(timereg)
library(aftgee)
library(table1)
library(ciTools)
library(scales)
library(readxl)
library(lubridate)
set.seed(123)

#function for survival of AFT model
survreg_calc_probs <- function(df, fit, q, comparison){
  form <- formula(fit)
  m <- model.frame(form, df)
  mat <- model.matrix(form, m)
  
  dist <- survival::survreg.distributions[[fit$dist]][["dist"]]
  fn_list <- survival::survreg.distributions[[dist]][["density"]]
  cdf <- function(x) fn_list(x)[,1]
  
  pred <- predict(fit, df, type = "linear")
  scale <- fit$scale
  zeta <- (log(q) - pred) / scale
  
  if (comparison == "<" || comparison == "<=")
    F <- cdf(zeta)
  else if (comparison == ">" || comparison == ">=")
    F <- 1 - cdf(zeta)
  else
    stop("invalid comparison")
  
  return(list(
    F = F,
    mat = mat,
    zeta = zeta
  ))
}

# R code for reading in Dan Sargent's colon cancer dataset
dat1 <- read_excel("/Users/sushi5824907/Desktop/RA 2021/Debashis/C06_C07_04.29.2009.xls")
dat2 <- read_excel("/Users/sushi5824907/Desktop/RA 2021/Debashis/NCCTG_NSABP_09APR09.xls")
dat12<-rbind(dat1,dat2)
names(dat12) <- c("study","tx","clintx","age","stage","gender","daterand","recurstat","recurtime","dfsstat","dfstime","dstat","dtime")
#exclude the row for labeling
colon<-dat12[!is.na(dat12$study),]

#stage II and III 
colon0<-colon[which(colon$stage==2 | colon$stage==3),]
#exclude subject with missing rand date
colon0<-colon0[!is.na(colon0$daterand),]

#convert to factor
colon0$study <- as.factor(colon0$study)
colon0$tx <- as.factor(colon0$tx)
colon0$gender <- as.factor(colon0$gender)
colon0$stage <- as.factor(colon0$stage)
#convert to numeric
colon0$dstat <- as.numeric(colon0$dstat)

#8 years of follow-up is considered in the analysis
colon0$dyear <- colon0$dtime/365
colon0$dyear8<-ifelse(colon0$dyear<=8,colon0$dyear,8)
colon0$dstat8<-ifelse(colon0$dyear8<colon0$dyear & colon0$dstat==1,0,colon0$dstat)

#5 years of median follow-up
colon0$daterand_num<-as.numeric(ymd(colon0$daterand))
#median rand date for each trial
rand<-colon0 %>% 
         group_by(study) %>% 
         dplyr::summarise(median=median(daterand_num))
#median 5 year FU
rand$tday<-rand$median+5*365
colon0 <- left_join(colon0, rand, by = 'study')
#different FU for each subject
colon0$tyear<-(colon0$tday-colon0$daterand_num)/365

colon0$dyear5<-ifelse(colon0$dyear8<=colon0$tyear,colon0$dyear8,colon0$tyear)
#exclude patients with a randomization date but without any follow-up information
colon0<-colon0[colon0$dyear5 > 0,]
colon0$dstat5<-ifelse(colon0$dyear5<colon0$dyear8 & colon0$dstat8==1,0,colon0$dstat8)

#number of subjects in each treatment group
colon0 %>% 
  group_by(tx)%>%
  dplyr::summarise(n = n())

#covariates summary table by study
sum_tab <- colon0 %>% 
  group_by(study) %>% 
  dplyr::summarize(n = n(), 
                   stage = paste0(sum(stage == 3),' (',round((sum(stage == 3)/n),3)*100,'%)'),
                   age = paste0(round(mean(age),1),' (',sd = round(sd(age),1),')'),
                   sex = paste0(sum(gender == 1),' (',round((sum(gender == 1)/n),3)*100,'%)'))

#KM plots
par(mfrow=c(3,4),mar = c(2.9,2.9,2.9,2.9))
for (i in 1:12) {
  surv<-survfit(Surv(dyear5,dstat5)~tx,data=colon0,subset=(colon0$study == levels(colon0$study)[i]))
  plot(surv,lty=1:1,col=c('black','red'),cex.axis=1,ymin=0.4,xlim=c(0,8))
  tmp <- survdiff(Surv(dyear5,dstat5)~tx,data=colon0,subset=(colon0$study == levels(colon0$study)[i]))
  title(main=paste(levels(colon0$study)[i]),line=-10,cex.main = 1.5)
  title(xlab='Time (year)',ylab='Survival Probability',line=1.8,cex.lab = 1)
}

#vars for modeling
colon0$tx_age<-I(colon0$tx==2)*colon0$age
colon0$tx_gender<-I(colon0$tx==2)*I(colon0$gender==1)
colon0$tx_stage<-I(colon0$tx==2)*I(colon0$stage==3)
colon0 <- subset(colon0, select = c(dyear5,dstat5,tx,age,gender,stage,study,tx_age,tx_gender,tx_stage))

#one-trial-at-a-time analysis

#Cox Proportional Hazards Model
tables <- NULL
for (k in 1:12) {
  pred <- numeric(12)
  true <- numeric(12)
  target <- numeric(12)
  fit1<-coxph(Surv(dyear5,dstat5)~tx*age+tx*gender+tx*stage,data=colon0[colon0$study == levels(colon0$study)[k],])
  for (i in 1:12) {
    target[i]<-i
    c2<-colon0[colon0$study == levels(colon0$study)[i],]
    c2$hr<-as.numeric(unlist(summary(survfit(fit1,newdata=c2),times=5)$surv))
    c3<-aggregate(c2$hr, list(c2$tx), FUN=mean)
    pred[i]<-c3[2,2]-c3[1,2]
    fit2<-coxph(Surv(dyear5,dstat5)~tx*age+tx*gender+tx*stage,data=colon0[colon0$study == levels(colon0$study)[i],])
    c2$hr<-as.numeric(unlist(summary(survfit(fit2,newdata=c2),times=5)$surv))
    c3<-aggregate(c2$hr, list(c2$tx), FUN=mean)
    true[i]<-c3[2,2]-c3[1,2]
  }
  tables <- rbind(tables, cbind(pred[-k],true[-k],k,target[-k]))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("pred", "true",'study','target')
tables$pe<-(tables$pred-tables$true)^2
#study level PE
m1_s1<-tables %>%
  group_by(study) %>%
  dplyr::summarise(mean = mean(pe))
#overall PE
t_m1_s1<-tables %>%
  dplyr::summarise(mean = mean(pe))

#Additive Hazards Regression Model
tables <- NULL
for (k in 1:12) {
  pred <- numeric(12)
  true <- numeric(12)
  target <- numeric(12)
  fit1<-aalen(Surv(dyear5,dstat5==1)~const(tx)+const(age)+const(gender)+const(stage)+const(tx_age)+const(tx_gender)+const(tx_stage),
              data=colon0[colon0$study == levels(colon0$study)[k],],resample.iid=1)
  for (i in 1:12) {
    target[i]<-i
    c2<-colon0[colon0$study == levels(colon0$study)[i],]
    c2$hd<-predict(fit1,c2,times=5)$S0
    c3<-aggregate(c2$hd, list(c2$tx), FUN=mean)
    pred[i]<-c3[2,2]-c3[1,2]
    fit2<-aalen(Surv(dyear5,dstat5==1)~const(tx)+const(age)+const(gender)+const(stage)+const(tx_age)+const(tx_gender)+const(tx_stage),
                data=colon0[colon0$study == levels(colon0$study)[i],],resample.iid=1)
    c2$hd<-predict(fit2,c2,times=5)$S0
    c3<-aggregate(c2$hd, list(c2$tx), FUN=mean)
    true[i]<-c3[2,2]-c3[1,2]
  }
  tables <- rbind(tables, cbind(pred[-k],true[-k],k,target[-k]))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("pred","true",'study','target')
tables$pe<-(tables$pred-tables$true)^2
#study level PE
m1_s2<-tables %>%
  group_by(study) %>%
  dplyr::summarise(mean = mean(pe))
#overall PE
t_m1_s2<-tables %>%
  dplyr::summarise(mean = mean(pe))

#AFT Model-lognormal
tables <- NULL
for (k in 1:12) {
  pred <- numeric(12)
  true <- numeric(12)
  target <- numeric(12)
  fit1<-survreg(Surv(dyear5,dstat5)~tx*age+tx*gender+tx*stage,dist="lognormal",data=colon0,subset=(colon0$study == levels(colon0$study)[k]))
  for (i in 1:12) {
    target[i]<-i
    c2<-colon0[colon0$study == levels(colon0$study)[i],]
    pre<-add_probs(c2,fit1,q=5,comparison = ">")$prob_greater_than5
    c3<-aggregate(pre, list(c2$tx), FUN=mean)
    pred[i]<-c3[2,2]-c3[1,2]
    
    fit2<-survreg(Surv(dyear5,dstat5)~tx*age+tx*gender+tx*stage,dist="lognormal",data=colon0,subset=(colon0$study == levels(colon0$study)[i]))
    pre<-add_probs(c2,fit2,q=5,comparison = ">")$prob_greater_than5
    c3<-aggregate(pre, list(c2$tx), FUN=mean)
    true[i]<-c3[2,2]-c3[1,2]
  }
  tables <- rbind(tables, cbind(pred[-k],true[-k],k,target[-k]))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("pred", "true",'study','target')
tables$pe<-(tables$pred-tables$true)^2
#study level PE
m1_s3<-tables %>%
  group_by(study) %>%
  dplyr::summarise(mean = mean(pe))
#overall PE
t_m1_s3<-tables %>%
  dplyr::summarise(mean = mean(pe))


#left-one-trial-out analysis

#Cox Proportional Hazards Model
tables <- NULL
for (k in 1:12) {
  
  #weight observation by trial size
  p_colon0<-colon0[colon0$study != levels(colon0$study)[k],]
  study_n<-as.data.frame(table(p_colon0$study))
  colnames(study_n)[1] <- "study"
  p_colon0 <- left_join(p_colon0, study_n, by = 'study')
  p_colon0$wt<-(nrow(p_colon0)/length(unique(p_colon0$study)))/p_colon0$Freq
  
  fit1<-coxph(Surv(dyear5,dstat5)~tx*age+tx*gender+tx*stage,weights=wt,data=p_colon0)
  c2<-colon0[colon0$study == levels(colon0$study)[k],]
  c2$hr<-as.numeric(unlist(summary(survfit(fit1,newdata=c2),times=5)$surv))
  c3<-aggregate(c2$hr, list(c2$tx), FUN=mean)
  pred<-c3[2,2]-c3[1,2]
  fit2<-coxph(Surv(dyear5,dstat5)~tx*age+tx*gender+tx*stage,data=colon0[colon0$study == levels(colon0$study)[k],])
  c2$hr<-as.numeric(unlist(summary(survfit(fit2,newdata=c2),times=5)$surv))
  c3<-aggregate(c2$hr, list(c2$tx), FUN=mean)
  true<-c3[2,2]-c3[1,2]
  tables <- rbind(tables, cbind(pred,true,-k,k))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("pred", "true",'study','target')
tables$pe<-(tables$pred-tables$true)^2
#study level PE
m2_s1<-tables %>%
  group_by(study) %>%
  dplyr::summarise(mean = mean(pe))%>%
  arrange(desc(study))
#overall PE
t_m2_s1<-tables %>%
  dplyr::summarise(mean = mean(pe))

#Additive Hazards Regression Model
tables <- NULL
for (k in 1:12) {
  
  #weight observation by trial size
  p_colon0<-colon0[colon0$study != levels(colon0$study)[k],]
  study_n<-as.data.frame(table(p_colon0$study))
  colnames(study_n)[1] <- "study"
  p_colon0 <- left_join(p_colon0, study_n, by = 'study')
  p_colon0$wt<-(nrow(p_colon0)/length(unique(p_colon0$study)))/p_colon0$Freq
  
  fit1<-aalen(Surv(dyear5,dstat5==1)~const(tx)+const(age)+const(gender)+const(stage)+const(tx_age)+const(tx_gender)+const(tx_stage),
              weights=p_colon0$wt,data=p_colon0,resample.iid=1)
  c2<-colon0[colon0$study == levels(colon0$study)[k],]
  c2$hd<-predict(fit1,c2,times=5)$S0
  c3<-aggregate(c2$hd, list(c2$tx), FUN=mean)
  pred<-c3[2,2]-c3[1,2]
  fit2<-aalen(Surv(dyear5,dstat5==1)~const(tx)+const(age)+const(gender)+const(stage)+const(tx_age)+const(tx_gender)+const(tx_stage),
              data=colon0[colon0$study == levels(colon0$study)[k],],resample.iid=1)
  c2$hd<-predict(fit2,c2,times=5)$S0
  c3<-aggregate(c2$hd, list(c2$tx), FUN=mean)
  true<-c3[2,2]-c3[1,2]
  tables <- rbind(tables, cbind(pred,true,-k,k))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("pred", "true",'study','target')
tables$pe<-(tables$pred-tables$true)^2
#study level PE
m2_s2<-tables %>%
  group_by(study) %>%
  dplyr::summarise(mean = mean(pe)) %>%
  arrange(desc(study))
#overall PE
t_m2_s2<-tables %>%
  dplyr::summarise(mean = mean(pe))


#AFT Model-lognormal
tables <- NULL
for (k in 1:12) {
  
  #weight observation by trial size
  p_colon0<-colon0[colon0$study != levels(colon0$study)[k],]
  study_n<-as.data.frame(table(p_colon0$study))
  colnames(study_n)[1] <- "study"
  p_colon0 <- left_join(p_colon0, study_n, by = 'study')
  p_colon0$wt<-(nrow(p_colon0)/length(unique(p_colon0$study)))/p_colon0$Freq
  
  fit1<-survreg(Surv(dyear5,dstat5)~tx*age+tx*gender+tx*stage,dist="lognormal",weights=wt,data=p_colon0)
  c2<-colon0[colon0$study == levels(colon0$study)[k],]
  pre<-survreg_calc_probs(c2,fit1,q=5,comparison = ">")$F
  c3<-aggregate(pre, list(c2$tx), FUN=mean)
  pred<-c3[2,2]-c3[1,2]
  
  fit2<-survreg(Surv(dyear5,dstat5)~tx*age+tx*gender+tx*stage,dist="lognormal",data=colon0,subset=(colon0$study == levels(colon0$study)[k]))
  pre<-survreg_calc_probs(c2,fit2,q=5,comparison = ">")$F
  c3<-aggregate(pre, list(c2$tx), FUN=mean)
  true<-c3[2,2]-c3[1,2]
  
  tables <- rbind(tables, cbind(pred,true,-k,k))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("pred", "true",'study','target')
tables$pe<-(tables$pred-tables$true)^2
#study level PE
m2_s3<-tables %>%
  group_by(study) %>%
  dplyr::summarise(mean = mean(pe)) %>%
  arrange(desc(study))
#overall PE
t_m2_s3<-tables %>%
  dplyr::summarise(mean = mean(pe))


#meta-analysis

#Cox Proportional Hazards Model
tables <- NULL
for (k in 1:12) {
  
  #weight observation by trial size
  p_colon0<-colon0[colon0$study != levels(colon0$study)[k],]
  study_n<-as.data.frame(table(p_colon0$study))
  colnames(study_n)[1] <- "study"
  p_colon0 <- left_join(p_colon0, study_n, by = 'study')
  p_colon0$wt<-(nrow(p_colon0)/length(unique(p_colon0$study)))/p_colon0$Freq
  
  fit1<-coxph(Surv(dyear5,dstat5)~tx*age+tx*gender+tx*stage,weights=wt,data=p_colon0)
  c2<-p_colon0
  c2$hr<-as.numeric(unlist(summary(survfit(fit1,newdata=c2),times=5)$surv))
  c3<-aggregate(c2$hr, list(c2$tx), FUN=mean)
  pred<-c3[2,2]-c3[1,2]
  fit2<-coxph(Surv(dyear5,dstat5)~tx*age+tx*gender+tx*stage,data=colon0[colon0$study == levels(colon0$study)[k],])
  c2<-colon0[colon0$study == levels(colon0$study)[k],]
  c2$hr<-as.numeric(unlist(summary(survfit(fit2,newdata=c2),times=5)$surv))
  c3<-aggregate(c2$hr, list(c2$tx), FUN=mean)
  true<-c3[2,2]-c3[1,2]
  tables <- rbind(tables, cbind(pred,true,-k,k))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("pred","true",'study','target')
tables$pe<-(tables$pred-tables$true)^2
#study level PE
m3_s1<-tables %>%
  group_by(study) %>%
  dplyr::summarise(mean = mean(pe)) %>%
  arrange(desc(study))
#overall PE
t_m3_s1<-tables %>%
  dplyr::summarise(mean = mean(pe))


#Additive Hazards Regression Model
tables <- NULL
for (k in 1:12) {
  
  #weight observation by trial size
  p_colon0<-colon0[colon0$study != levels(colon0$study)[k],]
  study_n<-as.data.frame(table(p_colon0$study))
  colnames(study_n)[1] <- "study"
  p_colon0 <- left_join(p_colon0, study_n, by = 'study')
  p_colon0$wt<-(nrow(p_colon0)/length(unique(p_colon0$study)))/p_colon0$Freq
  
  fit1<-aalen(Surv(dyear5,dstat5==1)~const(tx)+const(age)+const(gender)+const(stage)+const(tx_age)+const(tx_gender)+const(tx_stage),
              weights=p_colon0$wt,data=p_colon0,resample.iid=1)
  c2<-p_colon0
  c2$hd<-predict(fit1,c2,times=5)$S0
  c3<-aggregate(c2$hd, list(c2$tx), FUN=mean)
  pred<-c3[2,2]-c3[1,2]
  fit2<-aalen(Surv(dyear5,dstat5==1)~const(tx)+const(age)+const(gender)+const(stage)+const(tx_age)+const(tx_gender)+const(tx_stage),
              data=colon0[colon0$study == levels(colon0$study)[k],],resample.iid=1)
  c2<-colon0[colon0$study == levels(colon0$study)[k],]
  c2$hd<-predict(fit2,c2,times=5)$S0
  c3<-aggregate(c2$hd, list(c2$tx), FUN=mean)
  true<-c3[2,2]-c3[1,2]
  tables <- rbind(tables, cbind(pred,true,-k,k))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("pred", "true",'study','target')
tables$pe<-(tables$pred-tables$true)^2
#study level PE
m3_s2<-tables %>%
  group_by(study) %>%
  dplyr::summarise(mean = mean(pe)) %>%
  arrange(desc(study))
#overall PE
t_m3_s2<-tables %>%
  dplyr::summarise(mean = mean(pe))

#AFT Model-lognormal
tables <- NULL
for (k in 1:12) {
  
  #weight observation by trial size
  p_colon0<-colon0[colon0$study != levels(colon0$study)[k],]
  study_n<-as.data.frame(table(p_colon0$study))
  colnames(study_n)[1] <- "study"
  p_colon0 <- left_join(p_colon0, study_n, by = 'study')
  p_colon0$wt<-(nrow(p_colon0)/length(unique(p_colon0$study)))/p_colon0$Freq
  
  fit1<-survreg(Surv(dyear5,dstat5)~tx*age+tx*gender+tx*stage,dist="lognormal",weights=wt,data=p_colon0)
  c2<-p_colon0
  pre<-survreg_calc_probs(c2,fit1,q=5,comparison = ">")$F
  c3<-aggregate(pre, list(c2$tx), FUN=mean)
  pred<-c3[2,2]-c3[1,2]
  
  fit2<-survreg(Surv(dyear5,dstat5)~tx*age+tx*gender+tx*stage,dist="lognormal",data=colon0,subset=(colon0$study == levels(colon0$study)[k]))
  c2<-colon0[colon0$study == levels(colon0$study)[k],]
  pre<-survreg_calc_probs(c2,fit2,q=5,comparison = ">")$F
  c3<-aggregate(pre, list(c2$tx), FUN=mean)
  true<-c3[2,2]-c3[1,2]
  
  tables <- rbind(tables, cbind(pred,true,-k,k))
}
tables<-as.data.frame(tables)
colnames(tables) <- c("pred","true",'study','target')
tables$pe<-(tables$pred-tables$true)^2
#study level PE
m3_s3<-tables %>%
  group_by(study) %>%
  dplyr::summarise(mean = mean(pe)) %>%
  arrange(desc(study))
#overall PE
t_m3_s3<-tables %>%
  dplyr::summarise(mean = mean(pe))

#plot
m1_s1$method<-'One-trial-at-a-time'
m1_s1$survival<-'Proportional Hazards'
m1_s2$method<-'One-trial-at-a-time'
m1_s2$survival<-'Additive Hazards'
m1_s3$method<-'One-trial-at-a-time'
m1_s3$survival<-'Log-normal'
m2_s1$method<-'Leave-one-trial-out'
m2_s1$survival<-'Proportional Hazards'
m2_s2$method<-'Leave-one-trial-out'
m2_s2$survival<-'Additive Hazards'
m2_s3$method<-'Leave-one-trial-out'
m2_s3$survival<-'Log-normal'
m3_s1$method<-'Meta-analysis'
m3_s1$survival<-'Proportional Hazards'
m3_s2$method<-'Meta-analysis'
m3_s2$survival<-'Additive Hazards'
m3_s3$method<-'Meta-analysis'
m3_s3$survival<-'Log-normal'

fdat<-rbind(m1_s1,m1_s2,m1_s3,m2_s1,m2_s2,m2_s3,m3_s1,m3_s2,m3_s3)
fdat$study<-abs(fdat$study)

fdat$method <- factor(fdat$method, levels=c("One-trial-at-a-time","Leave-one-trial-out","Meta-analysis"))
fdat$survival <- factor(fdat$survival, levels=c("Proportional Hazards","Additive Hazards","Log-normal"))

ggplot(fdat, aes(x=study, y=mean, color=method, shape=method)) + 
  geom_point(size=2.5) +
  #geom_line(linetype = "dashed",size=1) +
  labs(x = '', y = 'Prediction Errors') +
  scale_x_discrete(limits=c("C01","C02","C03","C04","C05","C06","C07",
                            "INT-0035","NCCTG-78-48-52","NCCTG-87-46-51","NCCTG-89-46-51","NCCTG-91-46-53")) +
  scale_y_continuous(breaks=seq(0,0.012,0.002),limits= c(0,0.012),expand = c(0.001,0.001),labels = comma) +
  #scale_y_continuous(labels = comma) +
  scale_color_manual(values = c("#0099f8","red3","green4")) +
  scale_shape_manual(values = c(16,3,4)) +
  facet_wrap(vars(survival), nrow = 1, scales = "free") +
  theme_bw() +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
        legend.position = 'bottom',
        legend.text = element_text(size = 15, margin = margin(r = 0.4, unit = 'cm')),
        legend.key.size = unit(2, 'lines'),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90),
        strip.text = element_text(size = 15),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank())


#table1

tab1<-cbind(rbind(t_m1_s1,t_m1_s2,t_m1_s3),rbind(t_m2_s1,t_m2_s2,t_m2_s3),rbind(t_m3_s1,t_m3_s2,t_m3_s3))
tab1

#table2

#proportional hazards
df_s1 = data.frame()
for (i in 1:12) {
  df<-rbind(m1_s1[i,],m2_s1[i,],m3_s1[i,])
  df_s1 <- rbind(df_s1,df)
}

#additive hazards
df_s2 = data.frame()
for (i in 1:12) {
  df<-rbind(m1_s2[i,],m2_s2[i,],m3_s2[i,])
  df_s2 <- rbind(df_s2,df)
}

#log-normal
df_s3 = data.frame()
for (i in 1:12) {
  df<-rbind(m1_s3[i,],m2_s3[i,],m3_s3[i,])
  df_s3 <- rbind(df_s3,df)
}

tab2<-cbind(subset(df_s1,select='mean'),subset(df_s2,select='mean'),subset(df_s3,select='mean'))
tab2


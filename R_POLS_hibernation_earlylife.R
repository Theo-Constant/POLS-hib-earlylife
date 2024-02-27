
############################################################
#XXXXXXXXXXXXX
#workspace prep
############################################################

#set directory
setwd(" ")

###load packages

library(psych)
library(car)
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(lmtest)
library(FactoMineR)
library(lme4)  
library(nlme) 
library(factoextra)
library(vcd)
library(missMDA)
library(MuMIn)
library(survival)
library(Hmisc)
library(glmm)
library(glmmTMB)
library(Hmisc)
library(Rtools)
library(cowplot)
library(MASS)
library(lmPerm)
library(GGally)
library(rptR)
library(fitdistrplus)


#Biological data 

Data<-read.csv2("data.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)
Data_beha<-read.csv2("data_behaviors.csv",header=TRUE,sep=";",dec=",", stringsAsFactors = FALSE)

######################################
#PCA hibernation
ACPhib <- Data[,c(10,16,17,18,23)]
ACPhib
R<-cor(ACPhib)
cortest.bartlett(R,n=34)
KMO(R)
#kmo for the date of first torpor is too low
#we make a new PCA without the date of the first torpor
ACPhib2 <- Data[,c(10,17,18,23)]
R2<-cor(ACPhib2)
cortest.bartlett(R2,n=34)
KMO(R2)
#kmo for temperature in torpor  is too low
#we make a new PCA without the temperature in torpor
ACPhib3 <- Data[,c(10,17,18)]
R3<-cor(ACPhib3)
cortest.bartlett(R3,n=34)
KMO(R3)
res.PCA.hib<-PCA(ACPhib3, graph = TRUE)

eig.val <- get_eigenvalue(res.PCA.hib)
eig.val

PC1hib<- res.PCA.hib$ind$coord[,1]
Data$PC1hib <- PC1hib
Data$PC1hib


#Extract variability in telomere erosion that is not explained by telomere length at the beginning of the period
DatadT <- Data[,c(2,3,34,37)]
DatadT<-na.omit(DatadT)
DatadT
test1_1<-lme(Delta_telomere_experiment~ Telomere_pre_hibernation*Sex , random=~1|Mother,method="ML", na.action = "na.fail", data=DatadT)
summary(model.avg(dredge(test1_1, fixed=~(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test1_1<-lme(Delta_telomere_experiment~ Telomere_pre_hibernation , random=~1|Mother,method="ML", na.action = "na.fail", data=DatadT)
summary(test1_1)
Res_delta_telomere<-residuals(test1_1)
#For simplicity, the Res_delta_telomere data was directly added to the excel file "data_POLS.csv" under the name "Res_delta_telomere".


#Extract variability in food intake during hibernation that is not explained by body mass variation or body mass at the beginning of the period
test2_1<-lme(Food_intake_hib~ PC1hib+Variation_mass_hibernation+Body_mass_post_hib, random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
vif(test2_1)
#body mass variation during hibernation and body mass after hibernation show a too high level of vif
#we take the body mass variation during hibernation
test2_2<-lme(Food_intake_hib~PC1hib*Sex+Variation_mass_hibernation*Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(model.avg(dredge(test2_2, fixed=~(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test2_3<-lme(Food_intake_hib~ PC1hib+Variation_mass_hibernation, random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
hist(resid(test2_3))
qqPlot(resid(test2_3)) 
bptest(test2_3)
summary(test2_3)
res_food_intake_hib<-residuals(test2_3)
Data$res_food_intake_hib <- res_food_intake_hib


###########################
#Creating a PCA for Female reproduction

DataF<-subset(Data,Data$Sex=="F")
ACPreproF <- DataF[,c(29,30,31,33)]
ACPreproF
R4<-cor(ACPreproF, use="na.or.complete")
R4
cortest.bartlett(R4,n=14)
KMO(R4)
#kmo for total offspring number is too low
#we make a new PCA without the total offspring number

ACPreproF1 <- DataF[,c(30,31,33)]
ACPreproF1
R5<-cor(ACPreproF1, use="na.or.complete")
cortest.bartlett(R5,n=16)
KMO(R5)
res.PCA.reproF<-PCA(ACPreproF1, scale.unit = TRUE, graph = TRUE)
eig.val <- get_eigenvalue(res.PCA.reproF)
eig.val
#For simplicity, the PC1repro data was directly added to the excel file "data_POLS.csv" under the name "PC1repro".

#################################
#behaviors

##########
#Time_latency

plot(fitdist(Data$Time_latency,"pois"))
plot(fitdist(Data$Time_latency,"norm"))
plot(fitdist(Data$Time_latency,"nbinom"))
#Time_latency followed a negative binomial distribution
#To take into account this type of distribution, we used 
#the method of Nakagawa et al., (2017) based on the function glmer.nb 
#from the package lme4 (version 1.1-29.).


r1<-glmer.nb(Time_latency ~Session+Sex+(1|Individual)+(1|Mother), data=Data_beha)
summary(r1)
#no sex effect so we make another model without the sex factor
r2<-glmer.nb(Time_latency ~Session+(1|Individual)+(1|Mother), data=Data_beha)
summary(r2)

thetaN2 <- getME(r2, "glmer.nb.theta")
lambda2 <- as.numeric(exp(fixef(r2) + 0.5 * (as.numeric(VarCorr(r2)$Individu))))
VarOlN2 <- log(1 + (1/lambda2) + (1/thetaN2)) # log-normal approximation
VarOlN2
c(VarOlN2 = VarOlN2)
ICCrawPop2 <- as.numeric(VarCorr(r2)$Individu)/(sum(as.numeric(VarCorr(r2))) +
                                                 VarOlN2)
c(ICCrawPop2 = ICCrawPop2)

#low repeatability of 0.26 
#latency is not considered as personality traits

#########
#Number_grooming 

plot(fitdist(Data$Number_grooming ,"pois"))
plot(fitdist(Data$Number_grooming ,"norm"))
plot(fitdist(Data$Number_grooming ,"nbinom"))
#Number_grooming followed a negative binomial distribution
#To take into account this type of distribution, we used 
#the method of Nakagawa et al., (2017) based on the function glmer.nb 
#from the package lme4 (version 1.1-29.).

r3<-glmer.nb(Number_grooming ~Session+Sex+(1|Individual)+(1|Mother), data=Data_beha)
summary(r3)
#no sex ans session effect so we make another model without the sex factor 
r4<-glmer.nb(Number_grooming ~(1|Individual)+(1|Mother), data=Data_beha)
summary(r4)

thetaN4 <- getME(r4, "glmer.nb.theta")
lambda4 <- as.numeric(exp(fixef(r4) + 0.5 * (as.numeric(VarCorr(r4)$Individu))))
VarOlN4 <- log(1 + (1/lambda4) + (1/thetaN4)) # log-normal approximation
VarOlN4
c(VarOlN4 = VarOlN4)
ICCrawPop4 <- as.numeric(VarCorr(r4)$Individu)/(sum(as.numeric(VarCorr(r4))) +
                                                 VarOlN4)
c(ICCrawPop4 = ICCrawPop4)

#no repeatability
#Number_grooming is not considered as personality traits

######
#Number_transition

plot(fitdist(Data$Number_transition,"pois"))
plot(fitdist(Data$Number_transition,"norm"))
plot(fitdist(Data$Number_transition,"nbinom"))
#Number_transition followed a negative binomial distribution

r5<-glmer.nb(Number_transition ~Session+Sex+(1|Individual), data=Data_beha)
summary(r5)
#no sex and session effect so we make another model without these factors 
r6<-glmer.nb(Number_transition ~(1|Individual), data=Data_beha)
summary(r6)

thetaN6 <- getME(r6, "glmer.nb.theta")
lambda6 <- as.numeric(exp(fixef(r6) + 0.5 * (as.numeric(VarCorr(r6)$Individu))))
VarOlN6 <- log(1 + (1/lambda6) + (1/thetaN6)) # log-normal approximation
VarOlN6
c(VarOlN6 = VarOlN6)
ICCrawPop6 <- as.numeric(VarCorr(r6)$Individu)/(sum(as.numeric(VarCorr(r6))) +
                                                 VarOlN6)
c(ICCrawPop6 = ICCrawPop6)

#moderate repeatability of 0.37 
#Number_transition is considered as personality traits

########
#Number_rearing

plot(fitdist(Data$Number_rearing,"pois"))
plot(fitdist(Data$Number_rearing,"norm"))
plot(fitdist(Data$Number_rearing,"nbinom"))
#Number_rearing followed a negative binomial distribution

r7<-glmer.nb(Number_rearing ~Session+Sex+(1|Individual)+(1|Mother), data=Data_beha)
summary(r7)
#no sex and session effect so we make another model without these factors 
r8<-glmer.nb(Number_rearing ~(1|Individual)+(1|Mother), data=Data_beha)
summary(r8)

thetaN8 <- getME(r8, "glmer.nb.theta")
lambda8 <- as.numeric(exp(fixef(r8) + 0.5 * (as.numeric(VarCorr(r8)$Individu))))
VarOlN8 <- log(1 + (1/lambda8) + (1/thetaN8)) # log-normal approximation
VarOlN8
c(VarOlN8 = VarOlN8)
ICCrawPop8 <- as.numeric(VarCorr(r8)$Individu)/(sum(as.numeric(VarCorr(r8))) +
                                                 VarOlN8)
c(ICCrawPop8 = ICCrawPop8)

#no repeatability
#Number_rearing is not considered as personality traits

##########################


#We make a correlation matrix to select the variables that will be tested subsequently with more complex models.
#we test males and females separately due to sex differences

matrixF <- DataF[,c(6,38,39,40,45,50,51,52,53)] 
matrixF
rcorr(as.matrix(matrixF))

DataM<-subset(Data,Data$Sex=="M")
matrixM <- DataM[,c(6,38,39,40,41,42,43,44,45,50,51,52,53)] 
matrixM
rcorr(as.matrix(matrixM))


###############################
#test for covariation between selected factors 

test3_1<-glmmTMB(Number_transition~res_food_intake_hib*Sex+(1|Mother),family=nbinom2,na.action = "na.fail", data=Data)
#we want to see all the best models that contain food intake and the random variable mother
summary(model.avg(dredge(test3_1, fixed=~cond(res_food_intake_hib)+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test3_2<-glmmTMB(Number_transition~res_food_intake_hib*Sex+(1|Mother),family=nbinom2,na.action = "na.fail", data=Data)
summary(test3_2)
hist(resid(test3_2))
qqPlot(resid(test3_2)) 
bptest(test3_2)

test4_1<-glmmTMB(Number_transition~PC1hib*Sex+(1|Mother),family=nbinom2,na.action = "na.fail", data=Data)
summary(model.avg(dredge(test4_1, fixed=~cond(PC1hib)+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test4_2<-glmmTMB(Number_transition~PC1hib+Sex+(1|Mother),family=nbinom2,na.action = "na.fail", data=Data)
summary(test4_2)
hist(resid(test4_2))
qqPlot(resid(test4_2)) 
bptest(test4_2)


test5_1<-glmmTMB(Number_transition~Growth_rate*Sex+(1|Mother),family=nbinom2,na.action = "na.fail", data=Data)
summary(model.avg(dredge(test5_1, fixed=~cond(Growth_rate)+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test5_2<-glmmTMB(Number_transition~Growth_rate+Sex+(1|Mother),family=nbinom2,na.action = "na.fail", data=Data)
summary(test5_2)
hist(resid(test5_2))
qqPlot(resid(test5_2)) 
bptest(test5_2)


DatadT <- Data[,c(2,3,45,50)]
DatadT<-na.omit(DatadT)
DatadT
test6_1<-glmmTMB(Number_transition~Res_delta_telomere*Sex+(1|Mother),family=nbinom2,na.action = "na.fail", data=DatadT)
summary(model.avg(dredge(test6_1, fixed=~cond(Res_delta_telomere)+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test6_2<-glmmTMB(Number_transition~Res_delta_telomere+(1|Mother),family=nbinom2,na.action = "na.fail", data=DatadT)
summary(test6_2)
hist(resid(test6_2))
qqPlot(resid(test6_2)) 
bptest(test6_2)

test7_1<-glmmTMB(Number_transition~PC1repro+(1|Mother),family=nbinom2,na.action = "na.fail", data=DataF)
summary(test7_1)
hist(resid(test7_1))
qqPlot(resid(test7_1)) 
bptest(test7_1)

test8_1<-lme(PC1hib ~Growth_rate*Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(model.avg(dredge(test8_1, fixed=~Growth_rate+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test8_2<-lme(PC1hib ~  Growth_rate      , random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(test8_2)
hist(resid(test8_2))
qqPlot(resid(test8_2)) 
bptest(test8_2)

test9_1<-lme(PC1hib ~  res_food_intake_hib*Sex      , random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(model.avg(dredge(test9_1, fixed=~res_food_intake_hib+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test9_2<-lme(PC1hib ~  res_food_intake_hib  , random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(test9_2)
hist(resid(test9_2))
qqPlot(resid(test9_2)) 
bptest(test9_2)


DatadT <- Data[,c(2,3,50,52)]
DatadT<-na.omit(DatadT)
DatadT
test10_1<-lme(PC1hib~Res_delta_telomere*Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=DatadT)
summary(model.avg(dredge(test10_1, fixed=~Res_delta_telomere+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test10_2<-lme(PC1hib ~  Res_delta_telomere+Sex     , random=~1|Mother,method="ML", na.action = "na.fail", data=DatadT)
summary(test10_2)
hist(resid(test10_2))
qqPlot(resid(test10_2)) 
bptest(test10_2)

test11_1<-lme(PC1hib~PC1repro, random=~1|Mother,method="ML", na.action = "na.fail", data=DataF)
summary(test11_1)
hist(resid(test11_1))
qqPlot(resid(test11_1)) 
bptest(test11_1)

DatadT <- DataF[,c(2,3,50,51)]
DatadT<-na.omit(DatadT)
DatadT
test12_1<-lme(PC1repro~Res_delta_telomere, random=~1|Mother,method="ML", na.action = "na.fail", data=DatadT)
summary(test12_1)
hist(resid(test12_1))
qqPlot(resid(test12_1)) 
bptest(test12_1)

test13_1<-lme(PC1repro~res_food_intake_hib, random=~1|Mother,method="ML", na.action = "na.fail", data=DataF)
summary(test13_1)
hist(resid(test13_1))
qqPlot(resid(test13_1)) 
bptest(test13_1)

test14_1<-lme(PC1repro~Growth_rate , random=~1|Mother,method="ML", na.action = "na.fail", data=DataF)
summary(test14_1)
hist(resid(test14_1))
qqPlot(resid(test14_1)) 
bptest(test14_1)

test15_1<-lme(Growth_rate~res_food_intake_hib*Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(model.avg(dredge(test15_1, fixed=~res_food_intake_hib+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test15_2<-lme(Growth_rate~res_food_intake_hib+Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(test15_2)
hist(resid(test15_2))
qqPlot(resid(test15_2)) 
bptest(test15_2)

DatadT <- Data[,c(2,3,6,50)]
DatadT<-na.omit(DatadT)
DatadT
test16_1<-lme(Growth_rate  ~Res_delta_telomere*Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=DatadT)
summary(model.avg(dredge(test16_1, fixed=~Res_delta_telomere+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test16_2<-lme(Growth_rate~Res_delta_telomere+Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=DatadT)
summary(test16_2)
hist(resid(test16_2))
qqPlot(resid(test16_2)) 
bptest(test16_2)

DatadT <- Data[,c(2,3,50,53)]
DatadT<-na.omit(DatadT)
DatadT
test17_1<-lme(res_food_intake_hib     ~Res_delta_telomere*Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=DatadT)
summary(model.avg(dredge(test17_1, fixed=~Res_delta_telomere+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test17_2<-lme(res_food_intake_hib~Res_delta_telomere, random=~1|Mother,method="ML", na.action = "na.fail", data=DatadT)
summary(test17_2)
hist(resid(test17_2))
qqPlot(resid(test17_2)) 
bptest(test17_2)



#creation of graphics for the article 

split.screen(c(2,2))

plot1<-ggplot(Data, aes(x=res_food_intake_hib, y=Number_transition, color=as.character(Sex)))+
  geom_point(size=3)+
  theme_classic()+
geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, linetype="dotted",
              aes(y=Number_transition, x=res_food_intake_hib,color=as.character(Sex)))+

  scale_color_manual(values = c("F" = "#4E84C4", "M" = "#293352"), name="Sex",labels=c("Female","Male"))+
  labs(y = "Exploratory behavior",
       x = "Relative food intake")
 
plot2<-ggplot(Data, aes(x=Res_delta_telomere, y=Number_transition, color=as.character(Sex)))+
  geom_point(size=3)+
  theme_classic()+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color="Black")+
  scale_color_manual(values = c("F" = "#4E84C4", "M" = "#293352"), name="Sex",labels=c("Female","Male"))+
  labs(y = "Exploratory behavior",
     x = "Relative telomere erosion")

plot3<-ggplot(Data, aes(x=PC1repro, y=Number_transition))+
  geom_point(size=3, color="#4E84C4")+
  theme_classic()+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color="Black")+
  labs(y = "Exploratory behavior",
       x = "PC1 reproduction")

plot4<-ggplot(Data, aes(x=Growth_rate, y=Number_transition, color=as.character(Sex)))+
  geom_point(size=3)+
  theme_classic()+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color="Black")+
  scale_color_manual(values = c("F" = "#4E84C4", "M" = "#293352"), name="Sex",labels=c("Female","Male"))+
  labs(y = "Exploratory behavior",
       x = "Growth rate")

plot_grid(plot1, plot2, plot3, plot4, labels = "AUTO")
 
#######################################"
#test the effect of early life conditions on trait covariation

#We used a Principal Component Analysis (PCA) to derive new non-correlated variables for trait covariation in males
#we take into account only those factors that show a significant influence in males
DataM_pace <- DataM[,c(6,45,50,53)]
DataM_pace
DataM_pace<-na.omit(DataM_pace)
R9<-cor(DataM_pace, use="na.or.complete")
cortest.bartlett(R9,n=17)
KMO(R9)
#KMO Growth rate is too low
#we make a new PCA without Growth rate
DataM_pace2 <- DataM[,c(45,50,53)]
R10<-cor(DataM_pace2, use="na.or.complete")
cortest.bartlett(R10,n=17)
KMO(R10)
res.PCA.PACEM<-PCA(DataM_pace2, graph = TRUE)
eig.val <- get_eigenvalue(res.PCA.PACEM)
eig.val

#We tested the effects of early life on trait covariation (first dimension of the PCA) in males
PC1polsM<- res.PCA.PACEM$ind$coord[,1]
DataM$PC1covariationM <- PC1polsM
DataM$PC1covariationM
test18_1<-lme(PC1covariationM ~  Littersize+Birth      , random=~1|Mother,method="ML", na.action = "na.fail", data=DataM)
summary(model.avg(dredge(test18_1, fixed=~(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test18_2<-lme(PC1covariationM ~  Littersize , random=~1|Mother,method="ML", na.action = "na.fail", data=DataM)
summary(test18_2)
hist(resid(test18_2))
qqPlot(resid(test18_2)) 
bptest(test18_2)

#We used a Principal Component Analysis (PCA) to derive new non-correlated variables for trait covariation in females
#we take into account only those factors that show a significant influence in females

DataF_pace <- DataF[,c(6,45,50,51,53)]
DataF_pace
DataF_pace<-na.omit(DataF_pace)
R11<-cor(DataF_pace, use="na.or.complete")
cortest.bartlett(R11,n=15)
KMO(R11)

#KMO res_food_intake_hib is too low
#we make a new PCA without res_food_intake_hib 
DataF_pace1 <- DataF[,c(6,45,50,51)]
DataF_pace1
R12<-cor(DataF_pace1, use="na.or.complete")
cortest.bartlett(R12,n=15)
KMO(R12)
#KMO Res_delta_telomere is too low
#we make a new PCA without Res_delta_telomere
DataF_pace2 <- DataF[,c(6,45,51)]
DataF_pace2
R13<-cor(DataF_pace2, use="na.or.complete")
cortest.bartlett(R13,n=16)
KMO(R13)
res.PCA.PACEF<-PCA(DataF_pace2, graph = TRUE)
eig.val <- get_eigenvalue(res.PCA.PACEF)
eig.val


#We tested the effects of early life on trait covariation (first dimension of the PCA) in males

PC1polsF<- res.PCA.PACEF$ind$coord[,1]
DataF$PC1covariationF <- PC1polsF
DataF$PC1covariationF
test19_1<-lme(PC1covariationF ~  Littersize+Birth      , random=~1|Mother,method="ML", na.action = "na.fail", data=DataF)
summary(model.avg(dredge(test19_1, fixed=~(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test19_2<-lme(PC1covariationF ~  Littersize , random=~1|Mother,method="ML", na.action = "na.fail", data=DataF)
summary(test19_2)
hist(resid(test19_2))
qqPlot(resid(test19_2)) 
bptest(test19_2)

#creation of graphics for the article 

ggplot(DataF, aes(x=Littersize, y=PC1covariationF))+
  geom_point(size=4, color="#4E84C4")+
  theme_classic()+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color="Black")+
  labs(y = "Female traits covariation",
       x = "Littersize")



# we suspected that hibernation pattern was influenced by the location of the hamster cages due to a temperature gradient in the room
#We test the influence of hamster location on hibernation patterns of hamsters in the room.
#Minimum body temperature in torpor and date of first torpor was the only hibernation pattern 
#that was significantly influenced by the location of the hamster cages. 
#We created a corrected value of the minimum body temperature in torpor and date of the first torpor 
#which corresponds to the residuals of the linear regression between these parameters 
#and the cage location in the room. 
test20_1<-lme(Temperature_torpor~Hamster_location, random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(test20_1)
hist(resid(test20_1))
qqPlot(resid(test20_1))
bptest(test20_1)
res_temperature_torpor<-residuals(test20_1)
Data$res_temperature_torpor<-res_temperature_torpor

test22_1<-lme(Torpor_bout_duration~Hamster_location+Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(test22_1)
hist(resid(test22_1))
qqPlot(resid(test22_1))
bptest(test22_1)

test23_1<-lme(Date_first_torpor~Hamster_location, random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(test23_1)
hist(resid(test23_1))
qqPlot(resid(test23_1))
bptest(test23_1)
res_date_first_torpor<-residuals(test23_1)
Data$res_date_first_torpor<-res_date_first_torpor

test24_1<-lme(Date_last_torpor~Hamster_location+Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(test24_1)
hist(resid(test24_1))
qqPlot(resid(test24_1))
bptest(test24_1)


test25_1<-lme(Time_torpor~Hamster_location+Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(test25_1)
hist(resid(test25_1))
qqPlot(resid(test25_1))
bptest(test25_1)

#Hamster locatio showed signiificant effect on Temperature_torpor and Date_first_torpor  
#We used these corrected variables to create the hibernation PCA and test an influence on the covariation between hibernation PCA and other traits
ACPhib4 <- Data[,c(10,17,18,54,55)]
ACPhib4
R13<-cor(ACPhib4)
cortest.bartlett(R13,n=34)
KMO(R13)

#KMO date of corrected date of first torpor is too low (as with the original variable)
#we make a new PCA without corrected date of first torpor
ACPhib5 <- Data[,c(10,17,18,53)]
ACPhib5
R14<-cor(ACPhib5)
cortest.bartlett(R14,n=34)
KMO(R14)
res.PCA.hib2<-PCA(ACPhib5, graph = TRUE)

eig.val <- get_eigenvalue(res.PCA.hib2)
eig.val

PC1hib2<- res.PCA.hib2$ind$coord[,1]
Data$PC1hib2 <- PC1hib2
Data$PC1hib2

test26_1<-glmmTMB(Number_transition~PC1hib2*Sex+(1|Mother),family=nbinom2,na.action = "na.fail", data=Data)
summary(model.avg(dredge(test26_1, fixed=~cond(PC1hib2)+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test26_2<-glmmTMB(Number_transition~PC1hib2+Sex+(1|Mother),family=nbinom2,na.action = "na.fail", data=Data)
summary(test26_2)
hist(resid(test26_2))
qqPlot(resid(test26_2)) 
bptest(test26_2)

test27_1<-lme(PC1hib2 ~Growth_rate*Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(model.avg(dredge(test27_1, fixed=~Growth_rate+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test27_2<-lme(PC1hib2 ~  Growth_rate      , random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(test27_2)
hist(resid(test27_2))
qqPlot(resid(test27_2)) 
bptest(test27_2)

test28_1<-lme(PC1hib2 ~  res_food_intake_hib*Sex      , random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(model.avg(dredge(test28_1, fixed=~res_food_intake_hib+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test28_2<-lme(PC1hib2 ~  res_food_intake_hib     , random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(test28_2)
hist(resid(test28_2))
qqPlot(resid(test28_2)) 
bptest(test28_2)


DatadT <- Data[,c(2,3,50,56)]
DatadT<-na.omit(DatadT)
DatadT
test29_1<-lme(PC1hib2~Res_delta_telomere*Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=DatadT)
summary(model.avg(dredge(test29_1, fixed=~Res_delta_telomere+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test29_2<-lme(PC1hib2 ~  Res_delta_telomere+Sex     , random=~1|Mother,method="ML", na.action = "na.fail", data=DatadT)
summary(test29_2)
hist(resid(test29_2))
qqPlot(resid(test29_2)) 
bptest(test29_2)

test30_1<-lme(PC1hib~PC1repro, random=~1|Mother,method="ML", na.action = "na.fail", data=DataF)
summary(test30_1)
hist(resid(test30_1))
qqPlot(resid(test30_1)) 
bptest(test30_1)

#The results obtained are similar to those obtained with the uncorrected variables of torpor temperature 
#and date of first torpor.


#We tested whether taking daily torpor into account could affect the results obtained.
#We performed a hibernation PCA using the multiday torpor data only.

ACPhib6 <- Data[,c(11,14,15,21,22)]
ACPhib6
R15<-cor(ACPhib6)
cortest.bartlett(R15,n=34)
KMO(R15)

#KMO date of corrected date of first torpor is too low (as with the original variable)
#we make a new PCA without corrected date of first torpor
ACPhib7 <- Data[,c(11,15,21,22)]
ACPhib7
R16<-cor(ACPhib7)
cortest.bartlett(R16,n=34)
KMO(R16)
res.PCA.hib3<-PCA(ACPhib7, graph = TRUE)
PC1hib3<- res.PCA.hib3$ind$coord[,1]
Data$PC1hib3 <- PC1hib3
Data$PC1hib3


test31_1<-glmmTMB(Number_transition~PC1hib3*Sex+(1|Mother),family=nbinom2,na.action = "na.fail", data=Data)
summary(model.avg(dredge(test31_1, fixed=~cond(PC1hib3)+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test31_2<-glmmTMB(Number_transition~PC1hib3+Sex+(1|Mother),family=nbinom2,na.action = "na.fail", data=Data)
summary(test31_2)
hist(resid(test31_2))
qqPlot(resid(test31_2)) 
bptest(test31_2)

test32_1<-lme(PC1hib3 ~Growth_rate*Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(model.avg(dredge(test32_1, fixed=~Growth_rate+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test32_2<-lme(PC1hib3 ~  Growth_rate+Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(test32_2)
hist(resid(test32_2))
qqPlot(resid(test32_2)) 
bptest(test32_2)

test33_1<-lme(PC1hib3 ~  res_food_intake_hib*Sex      , random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(model.avg(dredge(test33_1, fixed=~res_food_intake_hib+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test33_2<-lme(PC1hib3 ~  res_food_intake_hib     , random=~1|Mother,method="ML", na.action = "na.fail", data=Data)
summary(test33_2)
hist(resid(test33_2))
qqPlot(resid(test33_2)) 
bptest(test33_2)


DatadT <- Data[,c(2,3,50,57)]
DatadT<-na.omit(DatadT)
DatadT
test34_1<-lme(PC1hib3~Res_delta_telomere*Sex, random=~1|Mother,method="ML", na.action = "na.fail", data=DatadT)
summary(model.avg(dredge(test34_1, fixed=~Res_delta_telomere+(1|Mother),rank="AICc",m.lim = c(NA,4)),delta<5))
test34_2<-lme(PC1hib3 ~  Res_delta_telomere+Sex     , random=~1|Mother,method="ML", na.action = "na.fail", data=DatadT)
summary(test34_2)
hist(resid(test34_2))
qqPlot(resid(test34_2)) 
bptest(test34_2)

DataF<-subset(Data,Data$Sex=="F")
test35_1<-lme(PC1hib3~PC1repro, random=~1|Mother,method="ML", na.action = "na.fail", data=DataF)
summary(test35_1)
hist(resid(test35_1))
qqPlot(resid(test35_1)) 
bptest(test35_1)

#The results obtained are similar to those obtained with the variable including daily torpor
### instead probably should model this as multinomial of retaining or losing. Subtract from the total

data<-read.csv("MultiTerritoryRepeatedMeasures.csv")
rm(list=ls())
load("/Volumes/hofmannlab/Shared/Manuscripts/N. multifasciatus social networks/Figures/code to produce figures/Figures5+6/shells_data.RDATA")head(data)
data<-data[,c("TreatComCont","Replicate","M1","M2","M3","M4","M5",
              "F1","F2","F3","F4","F5",
              "J1","J2","J3","J4","J5","filter_.")]

#install.packages('reshape2')
library(reshape2)
library(dplyr)

data.long<-melt(data,id.vars=c("TreatComCont","filter_.","Replicate"))
data.long$sex<-as.factor(unlist(strsplit(as.character(
  data.long$variable),split=NULL))[seq(from=1,to=5399,by=2)])

data.long$TreatComCont<-as.factor(data.long$TreatComCont)
data.long$TreatComCont<-plyr::revalue(data.long$TreatComCont,c("1"="R2","2"="R4","3"="A2","4"="A4","5"="AR2","6"="AR4","7"="Control"))

data.long$treatment_type<-plyr::revalue(data.long$TreatComCont,c("R2"="R","R4"="R","A2"="A","A4"="A","AR2"="AR","AR4"="AR","Control"="Control"))
data.long$amount<-as.numeric(as.character(plyr::revalue(data.long$TreatComCont,c("R2"="2","R4"="4","A2"="2","A4"="4","AR2"="2","AR4"="4","Control"="0"))))
data.long$observation<-as.numeric(unlist(strsplit(as.character(
  data.long$variable),split=NULL))[seq(from=2,to=5400,by=2)])
library(lme4)
library(ordinal)
test<-
  filter(data.long, TreatComCont != "Control") %>%
  filter(observation == 1) %>%
  mutate(value1 = value) %>%
  dplyr::select(TreatComCont,Replicate,sex,value1) %>%
  left_join(filter(data.long, TreatComCont != "Control"),.) %>%
  mutate(delta.ind = value - value1)

males<-filter(test,sex=="F")
males$delta.ind<-factor(males$delta.ind,ordered=T)
males$Replicate<-factor(males$Replicate,ordered=F)

males.mod<-clmm(delta.ind~ observation+TreatComCont + (1|Replicate),data=males, link = "logit")
males$delta.ind
levels(as.factor(data.long$value))

head(data.long)
males<-filter(data.long,sex=="M")
males$bottom<-as.factor(rownames(males))

summary(lmer(value~observation*TreatComCont+(1|Replicate)
                ,data=males))

#separate control and exp
data(wine)
data.long.exp<-filter(data.long,TreatComCont!="Control")
data.long.control<-filter(data.long,TreatComCont=="Control")



data.group.uncensored<-group_by(data.long.exp,TreatComCont,treatment_type,line_type,observation,sex,variable)
data.group.censored<-group_by(filter(data.long.exp,filter_.!=0),TreatComCont,treatment_type,line_type,observation,sex,variable)
data.group.uncensored.control<-group_by(data.long.control,observation,sex)
data.group.censored.control<-group_by(filter(data.long.control,filter_.!=0),observation,sex)
data.group.uncensored.control.IQR<-summarise(data.group.uncensored.control,mean=mean(value,na.rm=TRUE),
                                             sem= sd(value,na.rm=TRUE)/sqrt(length(value)))
data.group.censored.control.IQR<-summarise(data.group.censored.control,mean=mean(value,na.rm=TRUE),
                                             sem= sd(value,na.rm=TRUE)/sqrt(length(value)))

data.group.uncensored.control.IQR<-select(mutate(data.group.uncensored.control.IQR,lower=mean-sem,upper=mean+sem),observation,sex,lower,upper)
data.group.censored.control.IQR<-select(mutate(data.group.censored.control.IQR,lower=mean-sem,upper=mean+sem),observation,sex,lower,upper)


data.summ.uncensored<-dplyr::summarise(data.group.uncensored,mean=mean(value,na.rm=TRUE),
                             sem= sd(value,na.rm=TRUE)/sqrt(length(value)),
                             count=length(value[!is.na(value)]))
plot.data.uncensored<-left_join(data.summ.uncensored,data.group.uncensored.control.IQR,by=c("observation","sex"))

data.summ.censored<-dplyr::summarise(data.group.censored,mean=mean(value,na.rm=TRUE),
                                       sem= sd(value,na.rm=TRUE)/sqrt(length(value)),
                                       count=length(value[!is.na(value)]))
plot.data.censored<-left_join(data.summ.censored,data.group.censored.control.IQR,by=c("observation","sex"))

library(ggplot2)
plot.data.uncensored$sex<-factor(plot.data.uncensored$sex,levels=c("M","F","J"))
dodge <- position_dodge(width=0.4)
ggplot(plot.data.uncensored,aes(x=observation,y=mean,group=TreatComCont,ymax=mean+sem,ymin=mean-sem))+
  geom_point(position=dodge,aes(shape=line_type))+geom_line(aes(lty=line_type),position=dodge)+theme_bw()+
  geom_errorbar(position=dodge)+facet_grid(sex~treatment_type,scale="free_y")+
  geom_line(aes(y=upper),col="black")+geom_line(aes(y=lower),col="black")
ggsave("uncensored_data.pdf")

plot.data.censored$sex<-factor(plot.data.censored$sex,levels=c("M","F","J"))
ggplot(plot.data.censored,aes(x=observation,y=mean,group=TreatComCont,ymax=mean+sem,ymin=mean-sem))+
  geom_point(position=dodge,aes(shape=line_type))+geom_line(aes(lty=line_type),position=dodge)+theme_bw()+
  geom_errorbar(position=dodge)+facet_grid(sex~treatment_type,scale="free_y")+
  geom_line(aes(y=upper),col="black")+geom_line(aes(y=lower),col="black")
ggsave("censored_data.pdf")

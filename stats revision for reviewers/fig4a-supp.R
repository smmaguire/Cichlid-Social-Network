library(dplyr)
data<-read.csv('MultiTerritoryRepeatedMeasures.csv')
data$treatment_type<-plyr::revalue(data$TreatComCont,c("R2"="R","R4"="R","A2"="A","A4"="A","AR2"="AR","AR4"="AR","C"="Control"))
data$amount<-as.factor(as.character(plyr::revalue(data$TreatComCont,c("R2"="2","R4"="4","A2"="2","A4"="4","AR2"="2","AR4"="4","C"="0"))))
data$territory_bi<-as.factor(as.character(plyr::revalue(data$Territory,
                                                        c("Kept"="0","Evicted"="1"))))
data<-filter(data,Misplaced=="Tracked")
#data$TreatComCont<-factor(data$TreatComCont,levels = levels(data$TreatComCont)[c(5,1:4,6:7)])
data$Treat<-factor(data$Treat,levels = levels(data$Treat)[c(5:7,8,9,1:4)])

data$treatment_type<-factor(data$treatment_type,levels = levels(data$treatment_type)[c(3,1:2,4)])
data$territory_bi<-as.numeric(as.character((data$territory_bi)))

library(ggplot2)
#ggplot(data,aes(x=Treat,y=territory_bi,group=1)) +
# geom_point(shape=1,position=position_jitter(width=.05,height=.05))+
#stat_smooth(method='glm',family='binomial',se=F)

#levels(data$Treat) <- list(c1="Control I - No touch", 
#                          c2="Control II - Replace w broken shells", 
#                         c3="Control III - Add broken shells",
#                        r2="Replace 2", r4="Replace 4",
#                       a2="Add 2", a4="Add 4",
#                      ar2="Add/Replace 2",
#                     ar4="Add/Replace 4")


library(arm)
fit <- bayesglm(territory_bi~TreatComCont, data=data, family="binomial",
                drop.baseline=F)
summary(fit)
car::Anova(fit,type=3)

#library(multcomp)

#summary(glht(fit,linfct = mcp(Treat = c ("a2 - c2 = 0",
#                                         "a4 - c2 = 0",
#                                        "ar2 - c2:c3 = 0",
#                                       "ar4 - c2+c3 = 0",
#                                      "r2 - c3 = 0",
#                                     "r4 - c3 = 0",
#                                    "a4 - a2 = 0",
#                                   "ar4 - ar2 = 0",
#                                  "r4 - r2 = 0"))),
#test=adjusted(type="fdr"))
library(lsmeans)
lsmean.matrix<-lsmeans(fit,"TreatComCont",cont=list(my.own = list(a2.v.c = c(1,0,0,0,-1,0,0),
                                                           a4.v.c = c(0,1,0,0,-1,0,0),
                                                           ar2.v.c = c(0,0,1,0,-1,0,0),
                                                           ar4.v.c = c(0,0,0,1,-1,0,0),
                                                           r2.v.c = c(0,0,0,0,-1,1,0),
                                                           r4.v.c = c(0,0,0,0,-1,0,1),
                                                           a4.v.a2 = c(1,-1,0,0,0,0,0),
                                                           ar4.v.ar2 = c(0,0,1,-1,0,0,0),
                                                           r4.v.r2 = c(0,0,0,0,0,1,-1)
)),adjust="fdr")
summary(lsmean.matrix)


library(effects)
pdf(file="probabilities.pdf")
x<-summary(allEffects(fit))
plot(allEffects(fit))
table<-as.data.frame(effect("Treat",fit))
library(ggthemes)
ggplot(table,aes(x=Treat,y=fit))+geom_point(size=3)+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=.5)+
  ylim(0,1)+ylab("Probability of Territory Loss")+
  xlab("Treatment")+theme_few()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()
summary(mod)


#################### part B ########################
library(reshape2)
library(dplyr)

data<-read.csv("MultiTerritoryRepeatedMeasures.csv")
data$Treat<-factor(data$Treat,levels = levels(data$Treat)[c(5:7,8,9,1:4)])
#data<-filter(data,filter_.=="Selected")
data<-filter(data,Territory=="Kept",Misplaced=="Tracked")
data<-mutate(data,mdiff=M5-M1,fdiff=F5-F1,jdiff=J5-J1)
jmod<-(bayesglm(jdiff~Treat,data=data))
fmod<-(bayesglm(fdiff~Treat,data=data,drop.baseline=F))
mmod<-(bayesglm(mdiff~Treat,data=data,drop.baseline=F))

emale<-as.data.frame(effect("Treat",mmod))
emale$sex<-"M"
efemale<-as.data.frame(effect("Treat",fmod))
efemale$sex<-"F"
#ejuv<-as.data.frame(effect("Treat",jmod))
#ejuv$sex<-"J"
effects<-rbind(emale,efemale)

effects$sex<-factor(effects$sex,levels=c("M","F"))
effects$Treat<-factor(effects$Treat,levels=levels(effects$Treat)[c(5:7,8,9,1:4)])

pd <- position_dodge(0.7)
ggplot(effects,aes(x=Treat,y=fit,col=sex,shape=sex))+geom_point(position=pd)+
  geom_errorbar(aes(ymin=lower,ymax=upper),position=pd) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Average Change in the Number of Individuals")+scale_color_grey(start=0,end=.6)


male.lsmean<-lsmeans(fmod,"Treat",cont=list(my.own = list(a2.v.c2 = c(0,1,0,0,0,-1,0,0,0),
                                                          a4.v.c2 = c(0,1,0,0,0,0,-1,0,0),
                                                          ar2.v.cc = c(0,1,1,0,0,0,0,-1,0),
                                                          ar4.v.cc = c(0,1,1,0,0,0,0,0,-1),
                                                          r2.v.c3 = c(0,0,1,-1,0,0,0,0,0),
                                                          r4.v.c3 = c(0,0,1,0,-1,0,0,0,0),
                                                          a4.v.a2 = c(0,0,0,0,0,1,-1,0,0),
                                                          ar4.v.ar2 = c(0,0,0,0,0,0,0,1,-1),
                                                          r4.v.r2 = c(0,0,0,1,-1,0,0,0,0)
)),adjust="fdr")

car::Anova(mmod,type=3)
car::Anova(fmod,type=3)

male.lsmean<-lsmeans(mtest,"Treat",cont=list(my.own = list(a2.v.c2 = c(0,1,0,0,0,-1,0,0,0),
                                                           a4.v.c2 = c(0,1,0,0,0,0,-1,0,0),
                                                           ar2.v.cc = c(0,1,1,0,0,0,0,-1,0),
                                                           ar4.v.cc = c(0,1,1,0,0,0,0,0,-1),
                                                           r2.v.c3 = c(0,0,1,-1,0,0,0,0,0),
                                                           r4.v.c3 = c(0,0,1,0,-1,0,0,0,0),
                                                           a4.v.a2 = c(0,0,0,0,0,1,-1,0,0),
                                                           ar4.v.ar2 = c(0,0,0,0,0,0,0,1,-1),
                                                           r4.v.r2 = c(0,0,0,1,-1,0,0,0,0)
)),adjust="fdr")

summary(male.lsmean)
data %>%
  group_by(Treat) %>%
  summarise(n())





pdf("new_figure5.pdf",height=7,width=10)  
table$Treat<-factor(table$Treat,levels=c("c1","c2","c3","r2","r4","a2","a4",'ar2','ar4'))
ggplot(table,aes(x=Treat,y=fit))+geom_point(size=3)+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=.5)+
  ylim(0,1)+ylab("Probability of Territory Loss")+
  xlab("Treatment")+theme_few()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pd <- position_dodge(0.4)
ggplot(effects,aes(x=Treat,y=fit,col=sex,shape=sex))+geom_point(size=3,position=pd)+
  geom_errorbar(aes(ymin=lower,ymax=upper),position=pd,width=.5) + theme_few() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="none") +
  ylab("Change in the Number of Individuals (Day 5 - Day 1)")+
  scale_color_grey(start=0,end=.6)+geom_hline(yintercept=0)
dev.off()

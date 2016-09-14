#individual level analysis
rm(list=ls())
setwd("/Volumes/hofmannlab/Shared/All projects/N. multifasciatus Social Networks/individual level")
source('//files.ccbb.utexas.edu/hofmannlab/Shared/All projects/N. multifasciatus Social Networks/community level/refactored_functions_simulation.R')
source('//files.ccbb.utexas.edu/hofmannlab/Shared/All projects/N. multifasciatus Social Networks/raw data/iGraph data pullout_updated.R')
attributes<-getWTcom(networkList,matrixList)
template<-read.csv("//files.ccbb.utexas.edu/hofmannlab/Shared/All projects/N. multifasciatus Social Networks/individual level/template.csv")
for(gname in names(networkList)){
  m<-matrixList[[gname]]
  g<-networkList[[gname]]
  attr<-attributes[[gname]]
  degM<-degree_from_sex(m,"M",attr)
  degF<-degree_from_sex(m,"F",attr)
  cc<-transitivity(g,type="weighted")
  ccdf<-data.frame(unique=V(g)$name,cc=cc)
  pr<-data.frame(unique=names(page.rank(g)$vector),pageRank=page.rank(g)$vector)
  attr$day<-as.numeric(strsplit(gname,"t")[[1]][2])
  attr<-dplyr::left_join(attr,degM,by="unique")
  attr<-dplyr::left_join(attr,degF,by="unique")
  attr<-dplyr::left_join(attr,ccdf,by="unique")
  attr<-dplyr::left_join(attr,pr,by="unique")
  template<-rbind(template,attr)
}
ind.data<-template[-1,]
#manipulated group is A2, B3 and C5, D1, E5, F2.
ind.data<-ind.data[-which(ind.data$sex=="A"),]
ind.data<-mutate(ind.data,manip= ifelse(unique=="A2M1"|unique=="B3M1"|unique=="C5M1"|unique=="D1M1"|unique=="E5M1","manipulated","non-manipulated"))
ind.data<-filter(ind.data,sex=="M",day<=4)

pdf('ind_data.pdf',useDingbats = F)
ggplot(ind.data,aes(x=factor(day),y=outdegree_from_F,fill=manip))+geom_boxplot()+theme_bw()+ylab("Outdegree towards female")+
  scale_fill_grey(start=.5)+xlab("Day")+ylim(c(0,35))
ggplot(ind.data,aes(x=factor(day),y=indegree_from_F,fill=manip))+geom_boxplot()+theme_bw()+ylab("Indegree from female")+
  scale_fill_grey(start=.5)+xlab("Day")+ylim(c(0,35))
ggplot(ind.data,aes(x=factor(day),y=outdegree_from_M,fill=manip))+geom_boxplot()+theme_bw()+ylab("Outdegree towards males")+
  scale_fill_grey(start=.5)+xlab("Day")+ylim(c(0,35))
ggplot(ind.data,aes(x=factor(day),y=indegree_from_M,fill=manip))+geom_boxplot()+theme_bw()+ylab("indegree from male")+
  scale_fill_grey(start=.5)+xlab("Day")+ylim(c(0,35))
dev.off()

library(lme4)
require(lsmeans)
require(multcomp)
require(multcompView)
lsmean.matrix<-lsmeans(lmerm1,"Group.code",by="time")
tukey.result<-pairs(lsmean.matrix)
cld.result<-cld(lsmean.matrix)
ind.data$day<-as.factor(ind.data$day)
ind.data$manip<-as.factor(ind.data$manip)


out2F<-glmer(outdegree_from_F~day*manip+(1|unique)+(1|community),family=poisson,data=ind.data)
summary(out2F)
car::Anova(out2F,type=3)
lsmean.matrix<-lsmeans(out2F,"day",by="manip")
pairs(lsmean.matrix)
contrast(lsmean.matrix, "trt.vs.ctrl", ref = c(1,2,3))

inFromF<-glmer(indegree_from_F~day*manip+(1|unique)+(1|community),family=poisson,data=ind.data)

car::Anova(inFromF,type=3)
summary(inFromF)
lsmean.matrix<-lsmeans(inFromF,"day",by="manip")
pairs(lsmean.matrix)
contrast(lsmean.matrix, "trt.vs.ctrl", ref = c(1,2,3))

out2M<-glmer(outdegree_from_M~day*manip+(1|unique)+(1|community),family=poisson,data=ind.data)
summary(out2M)
lsmean.matrix<-lsmeans(out2M,"day",by="manip")
pairs(lsmean.matrix)
cld(lsmean.matrix)
contrast(lsmean.matrix, "trt.vs.ctrl", ref = c(1,2,3))


inFromM<-glmer(indegree_from_M~day*manip+(1|unique)+(1|community),family=poisson,data=ind.data)
summary(inFromM)
car::Anova(inFromM,type=3)
lsmean.matrix<-lsmeans(inFromM,"day",by="manip")
pairs(lsmean.matrix)
cld(lsmean.matrix)
contrast(lsmean.matrix, "trt.vs.ctrl", ref = c(1,2,3))



day.1.data<-filter(ind.data,day==1)
day.1.data$sub_group<-as.character(day.1.data$community)
in.group<-vector()
for(i in 1:nrow(day.1.data)){
  if(day.1.data$sub_group[i]=="A"&day.1.data$wtcom[i]==day.1.data$wtcom[day.1.data$unique=="A2M1"]) in.group<-c(in.group,"in.group") else 
    if(day.1.data$sub_group[i]=="A"&day.1.data$wtcom[i]!=day.1.data$wtcom[day.1.data$unique=="A2M1"]) in.group<-c(in.group,"out.group") else 
      if(day.1.data$sub_group[i]=="B"&day.1.data$wtcom[i]==day.1.data$wtcom[day.1.data$unique=="B3M1"]) in.group<-c(in.group,"in.group") else 
        if(day.1.data$sub_group[i]=="B"&day.1.data$wtcom[i]!=day.1.data$wtcom[day.1.data$unique=="B3M1"]) in.group<-c(in.group,"out.group") else 
          if(day.1.data$sub_group[i]=="C"&day.1.data$wtcom[i]==day.1.data$wtcom[day.1.data$unique=="C5M1"]) in.group<-c(in.group,"in.group") else 
            if(day.1.data$sub_group[i]=="C"&day.1.data$wtcom[i]!=day.1.data$wtcom[day.1.data$unique=="C5M1"]) in.group<-c(in.group,"out.group") else 
              if(day.1.data$sub_group[i]=="D"&day.1.data$wtcom[i]==day.1.data$wtcom[day.1.data$unique=="D1M1"]) in.group<-c(in.group,"in.group") else 
                if(day.1.data$sub_group[i]=="D"&day.1.data$wtcom[i]!=day.1.data$wtcom[day.1.data$unique=="D1M1"]) in.group<-c(in.group,"out.group") else 
                  if(day.1.data$sub_group[i]=="E"&day.1.data$wtcom[i]==day.1.data$wtcom[day.1.data$unique=="E5M1"]) in.group<-c(in.group,"in.group") else 
                    if(day.1.data$sub_group[i]=="E"&day.1.data$wtcom[i]!=day.1.data$wtcom[day.1.data$unique=="E5M1"]) in.group<-c(in.group,"out.group")
}
tojoin<-data.frame(unique=day.1.data$unique,in.group)
ind.data2<-left_join(ind.data,tojoin,by="unique")

pdf('in_or_out_group.pdf')
ggplot(ind.data2,aes(x=factor(day),y=indegree_from_M,col=manip))+geom_boxplot()+facet_wrap(sex~in.group)+ylab('indegree from males')
ggplot(ind.data2,aes(x=factor(day),y=outdegree_from_M,col=manip))+geom_boxplot()+facet_wrap(sex~in.group)+ylab('outdegree towards males')
ggplot(ind.data2,aes(x=factor(day),y=indegree_from_F,col=manip))+geom_boxplot()+facet_wrap(sex~in.group)+ylab('indegree from females')
ggplot(ind.data2,aes(x=factor(day),y=outdegree_from_F,col=manip))+geom_boxplot()+facet_wrap(sex~in.group)+ylab('outdegree towards females')
dev.off()

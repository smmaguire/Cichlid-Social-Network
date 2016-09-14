#Re-do for reviewers
library(dplyr)
setwd("/Volumes/hofmannlab/Shared/All projects/N. multifasciatus Social Networks/individual level")
source('/Volumes/hofmannlab/Shared/All projects/N. multifasciatus Social Networks/community level/refactored_functions_simulation.R')
source('/Volumes/hofmannlab/Shared/All projects/N. multifasciatus Social Networks/raw data/iGraph data pullout_updated_changed_for_reviewers.R')
#names:
#matrixList_col1
#matrixList_col2and3

template<-read.csv("/Volumes/hofmannlab/Shared/All projects/N. multifasciatus Social Networks/individual level/template.csv")
template<-template[,-12]
networkList<-networkList_col1
matrixList<-matrixList_col1
for(gname in names(networkList)){
  m<-matrixList[[gname]]
  g<-networkList[[gname]]
  attr<-attributes[[gname]]
  degM<-degree_from_sex(m,"M",attr)
  degF<-degree_from_sex(m,"F",attr)
  cc<-transitivity(g,type="weighted")
  ccdf<-data.frame(unique=V(g)$name,cc=cc)
  #pr<-data.frame(unique=names(page.rank(g)$vector),pageRank=page.rank(g)$vector)
  attr$day<-as.numeric(strsplit(gname,"t")[[1]][2])
  attr<-dplyr::left_join(attr,degM,by="unique")
  attr<-dplyr::left_join(attr,degF,by="unique")
  attr<-dplyr::left_join(attr,ccdf,by="unique")
  #attr<-dplyr::left_join(attr,pr,by="unique")
  template<-rbind(template,attr)
}
ind.data<-template[-1,]
#manipulated group is A2, B3 and C5, D1, E5, F2.
ind.data<-ind.data[-which(ind.data$sex=="A"),]
ind.data<-mutate(ind.data,manip= ifelse(unique=="A2M1"|unique=="B3M1"|unique=="C5M1"|unique=="D1M1"|unique=="E5M1","manipulated","non-manipulated"))
ind.data<-filter(ind.data,sex=="M",day<=4)

ind.data$manip<-factor(ind.data$manip,levels=c("non-manipulated","manipulated"))

p1<-ggplot(ind.data,aes(x=factor(day),y=outdegree_from_F,fill=manip))+geom_boxplot()+theme_bw()+ylab("Outdegree towards female")+
  scale_fill_grey(start=.5)+xlab("Day")+ylim(c(0,35))+theme(legend.position="none")
p2<-ggplot(ind.data,aes(x=factor(day),y=indegree_from_F,fill=manip))+geom_boxplot()+theme_bw()+ylab("Indegree from female")+
  scale_fill_grey(start=.5)+xlab("Day")+ylim(c(0,35))+theme(legend.position="none")
p3<-ggplot(ind.data,aes(x=factor(day),y=outdegree_from_M,fill=manip))+geom_boxplot()+theme_bw()+ylab("Outdegree towards males")+
  scale_fill_grey(start=.5)+xlab("Day")+ylim(c(0,35))+theme(legend.position="none")
p4<-ggplot(ind.data,aes(x=factor(day),y=indegree_from_M,fill=manip))+geom_boxplot()+theme_bw()+ylab("indegree from male")+
  scale_fill_grey(start=.5)+xlab("Day")+ylim(c(0,35))+theme(legend.position="none")

setwd("/Volumes/hofmannlab/Shared/Manuscripts/N. multifasciatus social networks/Supp figures")
pdf("supp_fig_3_ind_data_multiplot_col1.pdf",height=8,width=8)
multiplot(p3,p4,p1,p2,cols=2)
dev.off()

library(lme4)
require(lsmeans)
require(multcomp)
require(multcompView)
library(MASS)
#library(glmmADMB)
#lsmean.matrix<-lsmeans(lmerm1,"Group.code",by="time")
#tukey.result<-pairs(lsmean.matrix)
#cld.result<-cld(lsmean.matrix)
ind.data$day<-factor(ind.data$day,ordered=F)
ind.data$manip<-as.factor(ind.data$manip)
ind.data$bottom<-as.factor(as.numeric(rownames(ind.data)))

out2F<-glmer(outdegree_from_F~day*manip+
               (1|community/unique)+(1|bottom),
             family=poisson,data=ind.data)
summary(out2F)
car::Anova(out2F,type=3)

lsmean.matrix<-lsmeans(out2F,"day",by="manip")
pairs(lsmean.matrix)
contrast(lsmean.matrix, "trt.vs.ctrl", ref = c(1,2,3))

inFromF<-glmer(indegree_from_F~day*manip+(1|community/unique)+(1|bottom),family=poisson,data=ind.data)
summary(inFromF)
car::Anova(inFromF,type=3)
summary(inFromF)
lsmean.matrix<-lsmeans(inFromF,"day",by="manip")
pairs(lsmean.matrix)
contrast(lsmean.matrix, "trt.vs.ctrl", ref = c(1,2,3))

out2M<-glmer(outdegree_from_M~day*manip+(1|community/unique)+(1|bottom),family=poisson,data=ind.data)
summary(out2M)
car::Anova(out2M,type=3)
lsmean.matrix<-lsmeans(out2M,"day",by="manip")
pairs(lsmean.matrix)
cld(lsmean.matrix)
contrast(lsmean.matrix, "trt.vs.ctrl", ref = c(1,2,3))


inFromM<-glmer(indegree_from_M~day*manip+(1|community/unique)+(1|bottom),family=poisson,data=ind.data)

summary(inFromM)
car::Anova(inFromM,type=3)
lsmean.matrix<-lsmeans(inFromM,"day",by="manip")
pairs(lsmean.matrix)
cld(lsmean.matrix)
contrast(lsmean.matrix, "trt.vs.ctrl", ref = c(1,2,3))








####### col 2 and 3 #############
template<-read.csv("/Volumes/hofmannlab/Shared/All projects/N. multifasciatus Social Networks/individual level/template.csv")
template<-template[,-12]
networkList<-networkList_col2and3
matrixList<-matrixList_col2and3
for(gname in names(networkList)){
  m<-matrixList[[gname]]
  g<-networkList[[gname]]
  attr<-attributes[[gname]]
  degM<-degree_from_sex(m,"M",attr)
  degF<-degree_from_sex(m,"F",attr)
  cc<-transitivity(g,type="weighted")
  ccdf<-data.frame(unique=V(g)$name,cc=cc)
  #pr<-data.frame(unique=names(page.rank(g)$vector),pageRank=page.rank(g)$vector)
  attr$day<-as.numeric(strsplit(gname,"t")[[1]][2])
  attr<-dplyr::left_join(attr,degM,by="unique")
  attr<-dplyr::left_join(attr,degF,by="unique")
  attr<-dplyr::left_join(attr,ccdf,by="unique")
  #attr<-dplyr::left_join(attr,pr,by="unique")
  template<-rbind(template,attr)
}
ind.data<-template[-1,]
#manipulated group is A2, B3 and C5, D1, E5, F2.
ind.data<-ind.data[-which(ind.data$sex=="A"),]
ind.data<-mutate(ind.data,manip= ifelse(unique=="A2M1"|unique=="B3M1"|unique=="C5M1"|unique=="D1M1"|unique=="E5M1","manipulated","non-manipulated"))
ind.data<-filter(ind.data,sex=="M",day<=4)

ind.data$manip<-factor(ind.data$manip,levels=c("non-manipulated","manipulated"))

p1<-ggplot(ind.data,aes(x=factor(day),y=outdegree_from_F,fill=manip))+geom_boxplot()+theme_bw()+ylab("Outdegree towards female")+
  scale_fill_grey(start=.5)+xlab("Day")+ylim(c(0,35))+theme(legend.position="none")
p2<-ggplot(ind.data,aes(x=factor(day),y=indegree_from_F,fill=manip))+geom_boxplot()+theme_bw()+ylab("Indegree from female")+
  scale_fill_grey(start=.5)+xlab("Day")+ylim(c(0,35))+theme(legend.position="none")
p3<-ggplot(ind.data,aes(x=factor(day),y=outdegree_from_M,fill=manip))+geom_boxplot()+theme_bw()+ylab("Outdegree towards males")+
  scale_fill_grey(start=.5)+xlab("Day")+ylim(c(0,35))+theme(legend.position="none")
p4<-ggplot(ind.data,aes(x=factor(day),y=indegree_from_M,fill=manip))+geom_boxplot()+theme_bw()+ylab("indegree from male")+
  scale_fill_grey(start=.5)+xlab("Day")+ylim(c(0,35))+theme(legend.position="none")

setwd("/Volumes/hofmannlab/Shared/Manuscripts/N. multifasciatus social networks/Supp figures")
pdf("supp_fig4_ind_data_multiplot_col2and3.pdf",height=8,width=8)
multiplot(p3,p4,p1,p2,cols=2)
dev.off()

ind.data$day<-factor(ind.data$day,ordered=F)
ind.data$manip<-as.factor(ind.data$manip)
ind.data$bottom<-as.factor(as.numeric(rownames(ind.data)))

out2F<-glmer(outdegree_from_F~day*manip+
               (1|community/unique)+(1|bottom),
             family=poisson,data=ind.data)
summary(out2F)
car::Anova(out2F,type=3)

lsmean.matrix<-lsmeans(out2F,"day",by="manip")
pairs(lsmean.matrix)
contrast(lsmean.matrix, "trt.vs.ctrl", ref = c(1,2,3))

inFromF<-glmer(indegree_from_F~day*manip+(1|community/unique)+(1|bottom),family=poisson,data=ind.data)
summary(inFromF)
car::Anova(inFromF,type=3)
summary(inFromF)
lsmean.matrix<-lsmeans(inFromF,"day",by="manip")
pairs(lsmean.matrix)
contrast(lsmean.matrix, "trt.vs.ctrl", ref = c(1,2,3))

out2M<-glmer(outdegree_from_M~day*manip+(1|community/unique)+(1|bottom),family=poisson,data=ind.data)
summary(out2M)
car::Anova(out2M,type=3)
lsmean.matrix<-lsmeans(out2M,"day",by="manip")
pairs(lsmean.matrix)
cld(lsmean.matrix)
contrast(lsmean.matrix, "trt.vs.ctrl", ref = c(1,2,3))


inFromM<-glmer(indegree_from_M~day*manip+(1|community/unique)+(1|bottom),family=poisson,data=ind.data)

summary(inFromM)
car::Anova(inFromM,type=3)
lsmean.matrix<-lsmeans(inFromM,"day",by="manip")
pairs(lsmean.matrix)
cld(lsmean.matrix)
contrast(lsmean.matrix, "trt.vs.ctrl", ref = c(1,2,3))


###################### QAP #######################

lapply(c("sna","doMC","dplyr"),require,character.only=TRUE)
registerDoMC()
filters<-list(group = "as.numeric(attributes$sub_group)==sg",
              sex = "as.character(attributes$sex)==sx",
              group_sex = "as.numeric(attributes$sub_group)==sg & as.character(attributes$sex)==sx",
              group_adjgroup = "as.numeric(attributes$sub_group)<=sg+1|as.numeric(attributes$sub_group)>=sg-1",
              group_adjgroup_sex = "(as.numeric(attributes$sub_group)<=sg+1|as.numeric(attributes$sub_group)>=sg-1) & as.character(attributes$sex)==sx",
              random="random",
              walktrap="as.character(attributes$wtcom)==wtcom",
              walktrap_sex="as.character(attributes$wtcom)==wtcom & as.character(attributes$sex)==sx"
)
networkList<-networkList_col1
matrixList<-matrixList_col1
groupA<-returnArray(matrixList[1:7])
groupB<-returnArray(matrixList[8:14])
groupC<-returnArray(matrixList[15:21])
groupD<-returnArray(matrixList[22:28])
groupE<-returnArray(matrixList[29:35])
wtcom_A<-getWTcom(networkList[1:7],matrixList[1:7])
wtcom_B<-getWTcom(networkList[8:14],matrixList[8:14])
wtcom_C<-getWTcom(networkList[15:21],matrixList[15:21])
wtcom_D<-getWTcom(networkList[22:28],matrixList[22:28])
wtcom_E<-getWTcom(networkList[29:35],matrixList[29:35])

groupA.walktrap.col1<-runLoop_master_walktrap(groupA,"At",filters$walktrap,wtcom_A)
groupB.walktrap.col1<-runLoop_master_walktrap(groupB,"Bt",filters$walktrap,wtcom_B)
groupC.walktrap.col1<-runLoop_master_walktrap(groupC,"Ct",filters$walktrap,wtcom_C)
groupD.walktrap.col1<-runLoop_master_walktrap(groupD,"Dt",filters$walktrap,wtcom_D)
groupE.walktrap.col1<-runLoop_master_walktrap(groupE,"Et",filters$walktrap,wtcom_E)

networkList<-networkList_col2and3
matrixList<-matrixList_col2and3
groupA<-returnArray(matrixList[1:7])
groupB<-returnArray(matrixList[8:14])
groupC<-returnArray(matrixList[15:21])
groupD<-returnArray(matrixList[22:28])
groupE<-returnArray(matrixList[29:35])
wtcom_A<-getWTcom(networkList[1:7],matrixList[1:7])
wtcom_B<-getWTcom(networkList[8:14],matrixList[8:14])
wtcom_C<-getWTcom(networkList[15:21],matrixList[15:21])
wtcom_D<-getWTcom(networkList[22:28],matrixList[22:28])
wtcom_E<-getWTcom(networkList[29:35],matrixList[29:35])

groupA.walktrap.col2and3<-runLoop_master_walktrap(groupA,"At",filters$walktrap,wtcom_A)
groupB.walktrap.col2and3<-runLoop_master_walktrap(groupB,"Bt",filters$walktrap,wtcom_B)
groupC.walktrap.col2and3<-runLoop_master_walktrap(groupC,"Ct",filters$walktrap,wtcom_C)
groupD.walktrap.col2and3<-runLoop_master_walktrap(groupD,"Dt",filters$walktrap,wtcom_D)
groupE.walktrap.col2and3<-runLoop_master_walktrap(groupE,"Et",filters$walktrap,wtcom_E)

plot.perm.fnc.violin<-function(perms,title){
  require(ggplot2)
  plotDataGroup<-data.frame("firstQ"=vector(),"thirdQ"=vector(),"correlation"=vector(),
                            "group"=vector(),"pval"=vector())
  count<-1
  dist<-vector()
  group.dist<-vector()
  comp.dist<-vector()
  for (perm in perms) {
    firstQ<-vector()
    thirdQ<-vector()
    testVal<-vector()
    name.count<-1
    for(comp in perm[[1]]){
      firstQ<-c(firstQ,summary(comp$dist)[2])
      thirdQ<-c(thirdQ,summary(comp$dist)[5])
      testVal<-c(testVal,comp$testval)
      dist<-c(dist,comp$dist)
      group.dist<-c(group.dist,rep(LETTERS[count],length(comp$dist)))
      comp.dist<-c(comp.dist,rep(names(perm[[1]])[name.count],length(comp$dist)))
      name.count<-name.count+1
    }
    plotDataTemp<-(cbind(firstQ,thirdQ,correlation=testVal,group=LETTERS[count],pval=perm[[2]],deparse.level=1))
    rownames(plotDataTemp)<-NULL
    plotDataTemp<-as.data.frame(plotDataTemp)
    plotDataGroup<-rbind(plotDataGroup,plotDataTemp)
    count<-count+1
  }
  plotDataGroup$firstQ<-as.numeric(as.character(plotDataGroup$firstQ))
  plotDataGroup$thirdQ<-as.numeric(as.character(plotDataGroup$thirdQ))
  plotDataGroup$correlation<-as.numeric(as.character(plotDataGroup$correlation))
  
  plotDataGroup$timeComp<-as.factor(c('1 2','2 3','3 4'))
  dist.data<-data.frame(dist,timeComp=comp.dist,group=group.dist)
  violin<-ggplot(dist.data,aes(x=timeComp,y=dist,col=group))
  vplot<-violin+geom_violin()+geom_point(data=plotDataGroup,aes(y=correlation,x=timeComp))+facet_wrap(~group)+ggtitle(title)
  return(vplot)
}

p1<-plot.perm.fnc.violin(list(groupA.walktrap.col1,
                          groupB.walktrap.col1,
                          groupC.walktrap.col1,
                          groupD.walktrap.col1,
                          groupE.walktrap.col1),
                          "COL1 behaviors")
p2<-plot.perm.fnc.violin(list(groupA.walktrap.col2and3,
                          groupB.walktrap.col2and3,
                          groupC.walktrap.col2and3,
                          groupD.walktrap.col2and3,
                          groupE.walktrap.col2and3),
                     "COL 2 and 3 behaviors")
pdf("supp_fig_2.pdf",height=10,width=8)
multiplot(p1,p2)
dev.off()

2 group averages/ group
6 groups/ comm
5 comm
4 days
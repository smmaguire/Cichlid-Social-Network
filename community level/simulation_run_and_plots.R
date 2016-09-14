#should check that groups (adj/groups are handled correctly)
load('permutation_data_with_walktrap.RDATA')
source('/Volumes/hofmannlab/Shared/All projects/N. multifasciatus Social Networks/community level/refactored_functions_simulation.R')
source('/Volumes/hofmannlab/Shared/All projects/N. multifasciatus Social Networks/raw data/iGraph data pullout_updated.R')
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

groupA<-returnArray(matrixList[1:7])
groupB<-returnArray(matrixList[8:14])
groupC<-returnArray(matrixList[15:21])
groupD<-returnArray(matrixList[22:28])
groupE<-returnArray(matrixList[29:35])
groupA.group<-runLoop_master(groupA,attr_A,filt=filters$group)
groupB.group<-runLoop_master(groupB,attr_B,filt=filters$group)
groupC.group<-runLoop_master(groupC,attr_C,filt=filters$group)
groupD.group<-runLoop_master(groupD,attr_D,filt=filters$group)
groupE.group<-runLoop_master(groupE,attr_E,filt=filters$group)


groupA.sex<-runLoop_master(groupA,attr_A,filt=filters$sex)
groupB.sex<-runLoop_master(groupB,attr_B,filt=filters$sex)
groupC.sex<-runLoop_master(groupC,attr_C,filt=filters$sex)
groupD.sex<-runLoop_master(groupD,attr_D,filt=filters$sex)
groupE.sex<-runLoop_master(groupE,attr_E,filt=filters$sex)

groupA.group.sex<-runLoop_master(groupA,attr_A,filt=filters$group_sex)
groupB.group.sex<-runLoop_master(groupB,attr_B,filt=filters$group_sex)
groupC.group.sex<-runLoop_master(groupC,attr_C,filt=filters$group_sex)
groupD.group.sex<-runLoop_master(groupD,attr_D,filt=filters$group_sex)
groupE.group.sex<-runLoop_master(groupE,attr_E,filt=filters$group_sex)

groupA.adj.group<-runLoop_master(groupA,attr_A,filt=filters$group_adjgroup)
groupB.adj.group<-runLoop_master(groupB,attr_B,filt=filters$group_adjgroup)
groupC.adj.group<-runLoop_master(groupC,attr_C,filt=filters$group_adjgroup)
groupD.adj.group<-runLoop_master(groupD,attr_D,filt=filters$group_adjgroup)
groupE.adj.group<-runLoop_master(groupE,attr_E,filt=filters$group_adjgroup)

groupA.adj.group.sex<-runLoop_master(groupA,attr_A,filt=filters$group_adjgroup_sex)
groupB.adj.group.sex<-runLoop_master(groupB,attr_B,filt=filters$group_adjgroup_sex)
groupC.adj.group.sex<-runLoop_master(groupC,attr_C,filt=filters$group_adjgroup_sex)
groupD.adj.group.sex<-runLoop_master(groupD,attr_D,filt=filters$group_adjgroup_sex)
groupE.adj.group.sex<-runLoop_master(groupE,attr_E,filt=filters$group_adjgroup_sex)

groupA.random<-runLoop_master(groupA,attr_A,filt=filters$random)
groupB.random<-runLoop_master(groupB,attr_B,filt=filters$random)
groupC.random<-runLoop_master(groupC,attr_C,filt=filters$random)
groupD.random<-runLoop_master(groupD,attr_D,filt=filters$random)
groupE.random<-runLoop_master(groupE,attr_E,filt=filters$random)
#
groupA<-returnArray(matrixList[1:7])
groupB<-returnArray(matrixList[8:14])
groupC<-returnArray(matrixList[15:21])
groupD<-returnArray(matrixList[22:28])
groupE<-returnArray(matrixList[29:35])
#
wtcom_A<-getWTcom(networkList[1:7],matrixList[1:7])
wtcom_B<-getWTcom(networkList[8:14],matrixList[8:14])
wtcom_C<-getWTcom(networkList[15:21],matrixList[15:21])
wtcom_D<-getWTcom(networkList[22:28],matrixList[22:28])
wtcom_E<-getWTcom(networkList[29:35],matrixList[29:35])

groupA.walktrap<-runLoop_master_walktrap(groupA,"At",filters$walktrap,wtcom_A)
groupB.walktrap<-runLoop_master_walktrap(groupB,"Bt",filters$walktrap,wtcom_B)
groupC.walktrap<-runLoop_master_walktrap(groupC,"Ct",filters$walktrap,wtcom_C)
groupD.walktrap<-runLoop_master_walktrap(groupD,"Dt",filters$walktrap,wtcom_D)
groupE.walktrap<-runLoop_master_walktrap(groupE,"Et",filters$walktrap,wtcom_E)


groupA.walktrap.sex<-runLoop_master_walktrap(groupA,"At",filters$walktrap_sex,wtcom_A)
groupB.walktrap.sex<-runLoop_master_walktrap(groupB,"Bt",filters$walktrap_sex,wtcom_B)
groupC.walktrap.sex<-runLoop_master_walktrap(groupC,"Ct",filters$walktrap_sex,wtcom_C)
groupD.walktrap.sex<-runLoop_master_walktrap(groupD,"Dt",filters$walktrap_sex,wtcom_D)
groupE.walktrap.sex<-runLoop_master_walktrap(groupE,"Et",filters$walktrap_sex,wtcom_E)
# Label propagation instead
wtcom_A<-getLPcom(networkList[1:7],matrixList[1:7])
wtcom_B<-getLPcom(networkList[8:14],matrixList[8:14])
wtcom_C<-getLPcom(networkList[15:21],matrixList[15:21])
wtcom_D<-getLPcom(networkList[22:28],matrixList[22:28])
wtcom_E<-getLPcom(networkList[29:35],matrixList[29:35])

groupA.lp<-runLoop_master_walktrap(groupA,"At",filters$walktrap,wtcom_A)
groupB.lp<-runLoop_master_walktrap(groupB,"Bt",filters$walktrap,wtcom_B)
groupC.lp<-runLoop_master_walktrap(groupC,"Ct",filters$walktrap,wtcom_C)
groupD.lp<-runLoop_master_walktrap(groupD,"Dt",filters$walktrap,wtcom_D)
groupE.lp<-runLoop_master_walktrap(groupE,"Et",filters$walktrap,wtcom_E)


groupA.lp.sex<-runLoop_master_walktrap(groupA,"At",filters$walktrap_sex,wtcom_A)
groupB.lp.sex<-runLoop_master_walktrap(groupB,"Bt",filters$walktrap_sex,wtcom_B)
groupC.lp.sex<-runLoop_master_walktrap(groupC,"Ct",filters$walktrap_sex,wtcom_C)
groupD.lp.sex<-runLoop_master_walktrap(groupD,"Dt",filters$walktrap_sex,wtcom_D)
groupE.lp.sex<-runLoop_master_walktrap(groupE,"Et",filters$walktrap_sex,wtcom_E)
#plots
plot.perm.fnc<-function(perms,title){
  require(ggplot2)
  plotDataGroup<-data.frame("firstQ"=vector(),"thirdQ"=vector(),"correlation"=vector(),
                            "group"=vector(),"pval"=vector())
  count<-1
  for (perm in perms) {
    firstQ<-vector()
    thirdQ<-vector()
    testVal<-vector()
    for(comp in perm[[1]]){
      firstQ<-c(firstQ,summary(comp$dist)[2])
      thirdQ<-c(thirdQ,summary(comp$dist)[5])
      testVal<-c(testVal,comp$testval)
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
  
  plotDataGroup$timeComp<-as.factor(rep(c('1-2','2-3','3-4','4-5','5-6','6-7'),length(perms)))
  polydata<-as.data.frame(cbind(times=c(rep(c('1-2','2-3','3-4','4-5','5-6','6-7'),length(perms)),rep(rev(c('1-2','2-3','3-4','4-5','5-6','6-7')),length(perms))),points=c(plotDataGroup$firstQ,rev(plotDataGroup$thirdQ)),group=c(as.character(plotDataGroup$group),rev(as.character(plotDataGroup$group)))))
  polydata$points<-as.numeric(as.character(polydata$points))
  plotDataGroup$pval<-as.numeric(as.character(plotDataGroup$pval))
  plotDataGroup$pvalPlot<-ifelse(plotDataGroup$pval>.01,.05,NA)
  groupPlot<-ggplot(plotDataGroup,aes(x=timeComp,y=correlation,col=group,group=group))
  groupOnly<-groupPlot+geom_point()+geom_line()+geom_polygon(data=polydata,aes(x=times,y=points,group=group,fill=group),alpha=.5)+ggtitle(title)+xlab("Days Compared")+geom_point(aes(y=pvalPlot+correlation),col='black',shape=8)
  return(groupOnly)
}


plot.perm.fnc(list(groupA.group,groupB.group,groupC.group,groupD.group,groupE.group),"group only")
plot.perm.fnc(list(groupA.group.sex,groupB.group.sex,groupC.group.sex,groupD.group.sex,groupE.group.sex),"group and sex")
plot.perm.fnc(list(groupA.adj.group,groupB.adj.group,groupC.adj.group,groupD.adj.group,groupE.adj.group),"group and adj group")
plot.perm.fnc(list(groupA.adj.group.sex,groupB.adj.group.sex,groupC.adj.group.sex,groupD.adj.group.sex,groupE.adj.group.sex),"group,adj group and sex")
plot.perm.fnc(list(groupA.sex,groupB.sex,groupC.sex,groupD.sex,groupE.sex),"sex only")
plot.perm.fnc(list(groupA.random,groupB.random,groupC.random,groupD.random,groupE.random),"Random")

save.image("permutation_data_with_walktrap.RDATA")
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
  
  plotDataGroup$timeComp<-as.factor(rep(c('1 2','2 3','3 4','4 5','5 6','6 7'),length(perms)))
  dist.data<-data.frame(dist,timeComp=comp.dist,group=group.dist)
  violin<-ggplot(dist.data,aes(x=timeComp,y=dist,col=group))
  vplot<-violin+geom_violin()+geom_point(data=plotDataGroup,aes(y=correlation,x=timeComp))+facet_wrap(~group)+ggtitle(title)
  return(vplot)
}

pdf("violin_plots_permutation.pdf")
plot.perm.fnc.violin(list(groupA.group,groupB.group,groupC.group,groupD.group,groupE.group),"group only")
plot.perm.fnc.violin(list(groupA.group.sex,groupB.group.sex,groupC.group.sex,groupD.group.sex,groupE.group.sex),"group and sex")
plot.perm.fnc.violin(list(groupA.adj.group,groupB.adj.group,groupC.adj.group,groupD.adj.group,groupE.adj.group),"group and adj group")
plot.perm.fnc.violin(list(groupA.adj.group.sex,groupB.adj.group.sex,groupC.adj.group.sex,groupD.adj.group.sex,groupE.adj.group.sex),"group,adj group and sex")
plot.perm.fnc.violin(list(groupA.sex,groupB.sex,groupC.sex,groupD.sex,groupE.sex),"sex only")
plot.perm.fnc.violin(list(groupA.random,groupB.random,groupC.random,groupD.random,groupE.random),"Random")
plot.perm.fnc.violin(list(groupA.walktrap,groupB.walktrap,groupC.walktrap,groupD.walktrap,groupE.walktrap),"Walktrap")
plot.perm.fnc.violin(list(groupA.walktrap.sex,groupB.walktrap.sex,groupC.walktrap.sex,groupD.walktrap.sex,groupE.walktrap.sex),"Walktrap+sex")
plot.perm.fnc.violin(list(groupA.lp,groupB.lp,groupC.lp,groupD.lp,groupE.lp),"label propogation")
plot.perm.fnc.violin(list(groupA.lp.sex,groupB.lp.sex,groupC.lp.sex,groupD.lp.sex,groupE.lp.sex),"label propogation & sex")
dev.off()

#walktrap stats for paper
plotDataGroup$adj.pvals<-p.adjust(as.numeric(as.character(plotDataGroup$pval)),method="hochberg")
library(dplyr)

walktrap.first4days<-filter(plotDataGroup,timeComp %in% c("1 2","2 3","3 4"))
mmWt<-lme4::lmer(correlation~timeComp+(1|group),data=walktrap.first4days)
car::Anova(mmWt,type=3)
library(multcomp)
contr <- rbind("pre - post" = c(1, 1, -1),
               "pre - pre" = c(1,-1,0))

summary(glht(mmWt, linfct = mcp(timeComp = contr)))

dist.f4<-filter(dist.data,timeComp %in% c("1 2","2 3","3 4"))
violin<-ggplot(dist.f4,aes(x=timeComp,y=dist))
violin+geom_violin(fill="gray90")+geom_point(data=walktrap.first4days,aes(y=correlation,x=timeComp))+facet_grid(~group)+theme_bw()+ylab("Matrix Correlation Coefficient")+xlab("Days Correlated")


groups<-split(wtcom_E[[1]]$unique,wtcom_E[[1]]$sub_group)
colors<-ifelse(wtcom_E[[1]]$sex=="M", "#0072B2", "#D55E00")

tkplot(networkList[["Et1"]])
coords<-tkplot.getcoords(1)


group.colors <- c("#00000075", "#E69F0075", "#56B4E975", "#009E7375", "#F0E44275", "#CC79A775") #colorblind safe
group.border <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7") #colorblind safe

pdf("networkplots.pdf")
plot.igraph(networkList[["Et1"]],mark.groups=groups,layout=coords,main="Day 1",edge.width=get.edge.attribute(networkList[["Et1"]],"weight")/2,
            edge.arrow.size=.5,edge.arrow.width=.7,mark.col=group.colors,mark.border=group.border,vertex.color=colors,vertex.label="",
            vertex.frame.color=NA,edge.color='black')
plot.igraph(networkList[["Et2"]],mark.groups=groups,layout=coords,main="Day 2",edge.width=get.edge.attribute(networkList[["Et2"]],"weight")/2,
            edge.arrow.size=.5,edge.arrow.width=.7,mark.col=group.colors,mark.border=group.border,vertex.color=colors,vertex.label="",
            vertex.frame.color=NA,edge.color='black')
plot.igraph(networkList[["Et3"]],mark.groups=groups,layout=coords,main="Day 3",edge.width=get.edge.attribute(networkList[["Et3"]],"weight")/2,
            edge.arrow.size=.5,edge.arrow.width=.7,mark.col=group.colors,mark.border=group.border,vertex.color=colors,vertex.label="",
            vertex.frame.color=NA,edge.color='black')
plot.igraph(networkList[["Et4"]],mark.groups=groups,layout=coords,main="Day 4",
            edge.arrow.size=.5,edge.arrow.width=.7,mark.col=group.colors,mark.border=group.border,vertex.color=colors,vertex.label="",
            vertex.frame.color=NA,edge.color='black')
dev.off()

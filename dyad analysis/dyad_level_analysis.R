#g1<-graph.adjacency(group[day,,],mode='directed',weighted=TRUE,diag=FALSE)
#manipulated group is A2, B3 and C5, D1, E5, F2.

#dyads
library(igraph)
library(ggplot2)
source('//files.ccbb.utexas.edu/hofmannlab/Shared/All projects/N. multifasciatus Social Networks/community level/refactored_functions_simulation.R')
source('//files.ccbb.utexas.edu/hofmannlab/Shared/All projects/N. multifasciatus Social Networks/raw data/iGraph data pullout_updated.R')
groupA<-returnArray(matrixList[1:7])
groupB<-returnArray(matrixList[8:14])
groupC<-returnArray(matrixList[15:21])
groupD<-returnArray(matrixList[22:28])
groupE<-returnArray(matrixList[29:35])
#first dim is the time.
makeDyads<-function(group){
  require(dplyr)
  src<-vector();tar<-vector();srctar<-vector();weight<-vector();time<-vector()
  for(day in 1:dim(group)[1]){
    graph<-group[day,,]
    #go through each day, get source, target, source+target unique and the weight on the interaction
    for(animal1 in rownames(graph)){
      for(animal2 in colnames(graph)){
        #if(animal1==animal2) next
        src<-c(src,animal1)
        tar<-c(tar,animal2)
        srctar<-c(srctar,paste0(animal1,animal2))
        weight<-c(weight,graph[animal1,animal2])
        time<-c(time,day)
      }
    }
  }
data.frame(Source=src,Target=tar,Day=time,Dyad=srctar,Weight=weight)
}

dyadsA<-makeDyads(groupA)
dyadsB<-makeDyads(groupB)
dyadsC<-makeDyads(groupC)
dyadsD<-makeDyads(groupD)
dyadsE<-makeDyads(groupE)


#add other info
addAttr<-function(dyadT,attr){
for(type in c("Source","Target")){
attrT<-attr
names(attrT)<-paste(type,names(attrT),sep="_")
names(attrT)[1]<-type
dyadT<-left_join(dyadT,y=attrT,by=type)
}
temp<-mutate(dyadT,dyad_type=paste0(Source_sex,Target_sex),Within_group=ifelse(test=Source_sub_group==Target_sub_group,"In_group","Out_group"))
filter(temp,as.numeric(Source_sub_group)<=as.numeric(Target_sub_group)+1,
       as.numeric(Source_sub_group)>=as.numeric(Target_sub_group)-1,
       Source!=Target)
}
length(dyadsA[,1])
dyadsA<-addAttr(dyadsA,attr_A)
dyadsB<-addAttr(dyadsB,attr_B)
attr_C$sex[which(attr_C$sex=="A")]<-"M"
dyadsC<-addAttr(dyadsC,attr_C)
dyadsD<-addAttr(dyadsD,attr_D)
dyadsE<-addAttr(dyadsE,attr_E)
#gdyadA<-group_by(dyadsA,dyad_type,Source_sub_group,Day)

pdf('dyads.pdf',height=10,width=11)
ggplot(data=dyadsA,aes(x=factor(Day),y=Weight,col=dyad_type))+geom_boxplot()+facet_grid(dyad_type~Source_sub_group+Within_group)
ggplot(data=dyadsB,aes(x=factor(Day),y=Weight,col=dyad_type))+geom_boxplot()+facet_grid(dyad_type~Source_sub_group+Within_group)
ggplot(data=dyadsC,aes(x=factor(Day),y=Weight,col=dyad_type))+geom_boxplot()+facet_grid(dyad_type~Source_sub_group+Within_group)
dev.off()

add.manipulated<-function(Source_sub_group,group_number){
  if(as.numeric(Source_sub_group)==group_number) return ('manipulated') else 
    if(as.numeric(Source_sub_group)==group_number+1|as.numeric(Source_sub_group)==group_number-1) return("adj_manipulated") else 
      return("unmanipulated")
}
add.manipulated.df<-function(dyad_data,group_number){
dyad_group<-group_by(dyad_data,Source_sub_group)
temp<-mutate(dyad_group,manipulation=add.manipulated(Source_sub_group,group_number))
return(temp)
}

add.manipulated2<-function(Source_sub_group,group_number){
  if(as.numeric(Source_sub_group)==group_number) return ('manipulated') else return("unmanipulated")
}
add.manipulated.df2<-function(dyad_data,group_number){
  dyad_group<-group_by(dyad_data,Source_sub_group)
  dyad_group2<-group_by(dyad_data,Target_sub_group)
  temp<-mutate(dyad_group,source.man=add.manipulated2(Source_sub_group,group_number))
  temp2<-mutate(dyad_group2,target.man=add.manipulated2(Target_sub_group,group_number))
  temp3<-left_join(temp,temp2)
  return(temp3)
}

#need to make the dyads unique and not just pasted together I.E. A1M1A1M2 should be the same as A1M2A1M1

dyad.unique<-function(dyad.data){
  udy<-vector()
  for(row in 1:nrow(dyad.data)){
    temp1<-paste(dyad.data$Source[row],dyad.data$Target[row],sep="-")
    temp2<-paste(dyad.data$Target[row],dyad.data$Source[row],sep="-")
    if(temp1 %in% udy) udy<-c(udy,temp1) else if(temp2 %in% udy) udy<-c(udy,temp2) else udy<-c(udy,temp1)  
  }
  dyad.data$dyad.unique<-udy
  return(dyad.data)
}

dyadsA<-dyad.unique(dyadsA)
dyadsB<-dyad.unique(dyadsB)
dyadsC<-dyad.unique(dyadsC)
dyadsD<-dyad.unique(dyadsD)
dyadsE<-dyad.unique(dyadsE)

dyadsA.man<-add.manipulated.df2(dyadsA,2)
dyadsB.man<-add.manipulated.df2(dyadsB,3)
dyadsC.man<-add.manipulated.df2(dyadsC,5)
dyadsD.man<-add.manipulated.df2(dyadsD,1)
dyadsE.man<-add.manipulated.df2(dyadsE,5)

fullset<-as.data.frame(rbind(dyadsA.man,dyadsB.man,dyadsC.man,dyadsD.man,dyadsE.man))
head(fullset)
fullset$man_type<-as.factor(paste(fullset$source.man,fullset$target.man,sep="-"))
#females- unmanipulated on day 3----> and their summed interactions with their own male
fem.d3<-filter(fullset,Day==3,dyad_type=="MF"|dyad_type=="FM",man_type=="unmanipulated-unmanipulated",Within_group=="In_group")
fem.group<-group_by(fem.d3,dyad.unique)
fem.summ<-summarise(fem.group,weight.summ=sum(Weight))
fem.summ$female.involved<-unlist(strsplit(fem.summ$dyad.unique,split="-"))[seq(from=2,to=194,by=2)]
fem.summ$female.involved[fem.summ$female.involved=="C7AM1"]<-c("C7F1","C7F2","C7F3","C7F4")
fem.summ<-group_by(fem.summ,female.involved)
fem.summ<-summarise(fem.summ,weight.sumF=sum(weight.summ))
# Those same females and their summed interactions with the manipulated male
day4.fem.male<-filter(fullset,Source %in% fem.summ$female.involved | Target %in% fem.summ$female.involved,Day==4,dyad_type=="MF"|dyad_type=="FM",man_type=="unmanipulated-manipulated"| man_type== "manipulated-unmanipulated")
day4.fem.male.group<-group_by(day4.fem.male,dyad.unique)
day4.fem.male.summ<-summarise(day4.fem.male.group,weight.summ=sum(Weight))
day4.fem.male.summ$dyad1<-unlist(strsplit(day4.fem.male.summ$dyad.unique,split="-"))[seq(from=2,to=60,by=2)]
day4.fem.male.summ$dyad2<-unlist(strsplit(day4.fem.male.summ$dyad.unique,split="-"))[seq(from=1,to=59,by=2)]
write.csv(day4.fem.male.summ,"day4_male-fem.csv")
#fixed the csv
day4.man.data<-read.csv("day4_male-fem.csv")
names(day4.man.data)[2]<-"Interactions with manipulated Male"
names(fem.summ)[2]<-"Interactions with Own males"
joined<-left_join(day4.man.data,fem.summ,by="female.involved")
ggplot(joined,aes(x=`Interactions with Own males`,y=`Interactions with manipulated Male`))+xlab("Interactions within group males on Day3")+ylab("Interactions with the manipulated male on Day 4")+geom_point()

pdf('dyads_grouped_by_manipulation.pdf',height=10,width=11)
ggplot(data=fullset,aes(x=factor(Day),y=Weight,col=dyad_type))+geom_boxplot()+facet_grid(dyad_type~manipulation+Within_group)+ggtitle('all communities')
dev.off()

data2<-filter(fullset,manipulation=="manipulated"&dyad_type=="FM"&Day<=4&Within_group=="In_group")
x<-car::Anova(lm(Weight~factor(Day),data=data2),type=3)
TukeyHSD(aov(Weight~factor(Day),data=data2))

wtcom_A<-getWTcom(networkList[1:7],matrixList[1:7])
wtcom_B<-getWTcom(networkList[8:14],matrixList[8:14])
wtcom_C<-getWTcom(networkList[15:21],matrixList[15:21])
wtcom_D<-getWTcom(networkList[22:28],matrixList[22:28])
wtcom_E<-getWTcom(networkList[29:35],matrixList[29:35])

dyadsA<-makeDyads(groupA)
dyadsB<-makeDyads(groupB)
dyadsC<-makeDyads(groupC)
dyadsD<-makeDyads(groupD)
dyadsE<-makeDyads(groupE)

dyadsA<-addAttr(dyadsA,wtcom_A[[1]])
dyadsB<-addAttr(dyadsB,wtcom_B[[1]])
wtcom_C[[1]]$sex[which(wtcom_C[[1]]$sex=="A")]<-"M"
dyadsC<-addAttr(dyadsC,wtcom_C[[1]])
dyadsD<-addAttr(dyadsD,wtcom_D[[1]])
dyadsE<-addAttr(dyadsE,wtcom_E[[1]])
get.manipulated.wtcom<-function(man,attrs){
  return(unique(dplyr::filter(attrs,Day==1 & Source_sex == "M" & as.numeric(as.character(Source_sub_group)) == man)$Source_wtcom))
}
add.manipulated.wt<-function(Source_wtcom,group_number){
  if(as.numeric(Source_wtcom)==group_number) 
    return ('manipulated') 
  else 
    return("unmanipulated")
}

add.manipulated.df.wt<-function(dyad_data,group_number){
  dyad_group<-group_by(dyad_data,Source_sub_group)
  temp<-mutate(dyad_group,manipulation=add.manipulated.wt(Source_wtcom,group_number))
  return(temp)
}

dyadsA.man.wt<-add.manipulated.df.wt(dyadsA,get.manipulated.wtcom(2,dyadsA))
dyadsB.man.wt<-add.manipulated.df.wt(dyadsB,get.manipulated.wtcom(3,dyadsB))
dyadsC.man.wt<-add.manipulated.df.wt(dyadsC,get.manipulated.wtcom(5,dyadsC))
dyadsD.man.wt<-add.manipulated.df.wt(dyadsD,get.manipulated.wtcom(1,dyadsD))
dyadsE.man.wt<-add.manipulated.df.wt(dyadsE,get.manipulated.wtcom(5,dyadsE))

fullset<-rbind(dyadsA.man.wt,dyadsB.man.wt,dyadsC.man.wt,dyadsD.man.wt,dyadsE.man.wt)
ggplot(data=fullset,aes(x=factor(Day),y=Weight,col=dyad_type))+geom_boxplot()+facet_grid(dyad_type~manipulation+Within_group)+ggtitle('all communities')

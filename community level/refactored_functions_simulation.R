restricted_sample_master<-function(matrix,attributes,filt){
  require(dplyr)
  ids<-rownames(matrix)
  swapV<-rep(NA,length(ids))
  alreadySwaped<-vector()
  rand_ids<-sample(ids) # randomize order of ids
  for(ind in rand_ids){
    if(any(alreadySwaped==ind))(next) #skip ahead if that individual has already swaped
    indatt<-dplyr::filter(attributes,unique==ind)
    sg<-as.numeric(indatt$sub_group)
    sx<-as.character(indatt$sex)
    comm<-as.character(indatt$community)
    wtcom<-as.character(indatt$wtcom)
    #create a subset of the data that meets all the requierments (only within sex, only +1 or -1 of subgroup)
    # filter that subdata to take out animals that already swaped
    #ADD HERE TO FILTER DIFFERENT
    if(filt!="random") subdata<-dplyr::filter(attributes,eval(parse(text=filt))) else subdata<-attributes
    subdata<-subdata[!(Reduce('|', lapply(subdata, '%in%', alreadySwaped))),]
    # make a swap
    swaper<-sample(subdata$unique,1)
    swaper_number<-which(ids==swaper)
    ind_number<-which(ids==ind)
    swapV[ind_number]<-swaper_number
    swapV[swaper_number]<-ind_number
    # once the swap happens both the original indvidual and the swaper are out of the swaping pool
    alreadySwaped<-c(alreadySwaped,as.character(swaper),ind)
  }
  return(swapV)
}

restricted_rmperm_master<-function(m,attributes,filt){
  final<- array(dim=c(dim(m)[1],dim(m)[3],dim(m)[3]))
  new_order<-apply(X=m,MARGIN=1,FUN=restricted_sample_master,attributes=attributes,filt=filt) #change here to impliment differet sampling function
  for(i in 1:ncol(new_order)){
    final[i,,]<- array(m[i,new_order[,i],new_order[,i]])
  }
  final
}

new_hacked_QAP_master<-function (dat, reps = 1000,all_attr,g1,g2,filt) 
{
require(foreach)
  out <- list()
  dat<-dat[c(g1,g2),,]
  out$testval <- gcor(dat,g1=1,g2=2)
  out$dist<-foreach(i=1:reps,.combine='c') %dopar% gcor(restricted_rmperm_master(dat,all_attr,filt=filt),g1=1,g2=2)
  out$pgreq <- mean(as.numeric(out$dist >= out$testval))
  out$pleeq <- mean(as.numeric(out$dist <= out$testval))
  class(out) <- c("qaptest", "qap")
  out
}

runLoop_master<-function(groupA,all_attr,filt){
  p.value<-vector()
  #combos<-t(combn(1:7,2)) # all comparisons
  combos<-cbind(c(1,2,3,4,5,6),c(2,3,4,5,6,7)) #just sequential 
  pb<-txtProgressBar(min = 0, max = length(combos[,1]), style = 3)
  qapGroupA<-list()
  for(j in 1:length(combos[,1])){
    qap1<-new_hacked_QAP_master(groupA,reps=1000,all_attr,g1=combos[j,1],g2=combos[j,2],filt=filt)
    qapGroupA[[paste(combos[j,1],combos[j,2])]]<-qap1
    p.value<-c(p.value,qap1$pgreq)
    setTxtProgressBar(pb, j)
  }
  return(list(qapGroupA,p.value))
}


runLoop_master_walktrap<-function(groupA,groupCode,filt,wtcom_list){
  p.value<-vector()
  #combos<-t(combn(1:7,2)) # all comparisons
  combos<-cbind(c(1,2,3),c(2,3,4))
  pb<-txtProgressBar(min = 0, max = length(combos[,1]), style = 3)
  qapGroupA<-list()
  for(j in 1:length(combos[,1])){
    nameCom<-paste0(groupCode,combos[j,1])
    all_attr<-wtcom_list[nameCom][[1]]
    qap1<-new_hacked_QAP_master(groupA,reps=1000,all_attr,g1=combos[j,1],g2=combos[j,2],filt=filt)
    qapGroupA[[paste(combos[j,1],combos[j,2])]]<-qap1
    p.value<-c(p.value,qap1$pgreq)
    setTxtProgressBar(pb, j)
  }
  return(list(qapGroupA,p.value))
}

#example-community A
#wtcom_A<-getWTcom(networkList[1:7],matrixList[1:7])
#runLoop_master_walktrap(groupA,"At",filters$walktrap,wtcom_A)
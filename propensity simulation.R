source("Project propensity.R")

#set scenarios
scenarios<-c("A","B","C","D","E","F","G")

#set target parameter
par<-"ATE"

#set performance metrics
perfmetrics<-names(funsim( covfun(1000,"A"),"logit",par))

#number of replications
R<-1000

####### Simulation cycle 
for(size in c(500)){
  
  
  
  ## ps model results
  
  #  logit
  lgresults <- matrix(0,length(perfmetrics),1)
  timestart<-Sys.time()
  for(i in 1:length(scenarios)){
    pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
    partialresults<-matrix(0,length(perfmetrics),1)
    for(j in 1:R)
    {
      partialresults <-partialresults+try(as.matrix(simplify2array(funsim(
        covfun(size,substr(scenarios[i],1,1)),"logit",par))))
      setTxtProgressBar(pb, i)  
    }
    lgresults<-cbind(lgresults,partialresults/R)
  }
  lgresults<-lgresults[,-1];colnames(lgresults)<-paste(scenarios)
  
  
  timeend<-Sys.time();exectime<-timeend-timestart
  exectime 
  write.table(lgresults, file = file.path("results",paste("logit","R",R,"size",size,".txt",sep="")))
  save.image(file=file.path("results",paste("logit","R",R," size",size," par",par,".Rdata",sep="")))
  
  
  #bag
  bagresults <- matrix(0,length(perfmetrics),1)
  timestart<-Sys.time()
  for(i in 1:length(scenarios)){
    pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
    partialresults<-matrix(0,length(perfmetrics),1)
    for(j in 1:R)
    {
      partialresults <-partialresults+try(as.matrix(simplify2array(funsim(
        covfun(size,substr(scenarios[i],1,1)),"bag",par))))
      setTxtProgressBar(pb, i)  
    }
    bagresults<-cbind(bagresults,partialresults/R)
  }
  bagresults<-bagresults[,-1];colnames(bagresults)<-paste(scenarios)
  
  
  timeend<-Sys.time();exectime<-timeend-timestart
  exectime 
  write.table(bagresults, file = file.path("results",paste("bag","R",R,"size",size,".txt",sep="")))
  save.image(file=file.path("results",paste("bag","R",R," size",size," par",par,".Rdata",sep="")))
  
  
  # tree
  treeresults <- matrix(0,length(perfmetrics),1)
  timestart<-Sys.time()
  for(i in 1:length(scenarios)){
    pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3)
    #partialresults<-matrix(0,length(perfmetrics),1)
    partialresults<-list()
    for(j in 1:R){
      set.seed(j)
      prova<-covfun(size,substr(scenarios[i],1,1))
      
      partialresults[[j]]<- tryCatch(funsim(prova,"tree",par),
                                     error=function(e) paste("tree error at ",j,sep=""))
      setTxtProgressBar(pb, i)  
    }
    
    partialresults <- partialresults[lapply(partialresults,function(x) class(x))=="list"]
    partialresults <- sapply(partialresults,simplify2array)
    #partialresults <- sapply(partialresults,simplify2array)
    treeresults <- cbind(treeresults,apply(as.matrix(partialresults),1,mean))
  }
  treeresults<-treeresults[,-1];colnames(treeresults)<-paste(scenarios)
  
  write.table(treeresults, file = file.path("results",paste("tree","R",R,"size",size,".txt",sep="")))
  
  save.image(file=file.path("results",paste("tree","R",R," size",size," par",par,".Rdata",sep="")))
  
  
  
  # random forest
  rfresults <- matrix(0,length(perfmetrics),1)
  timestart<-Sys.time()
  for(i in 1:length(scenarios)){
    pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3)
    #partialresults<-matrix(0,length(perfmetrics),1)
    partialresults<-list()
    for(j in 1:R){
      set.seed(j)
      prova<-covfun(size,substr(scenarios[i],1,1))
      
      partialresults[[j]]<- tryCatch(funsim(prova,"randomforest",par),
                                     error=function(e) paste("rf error at ",j,sep=""))
      setTxtProgressBar(pb, i)  
    }
    
    partialresults <- partialresults[lapply(partialresults,function(x) class(x))=="list"]
    partialresults <- sapply(partialresults,simplify2array)
    #partialresults <- sapply(partialresults,simplify2array)
    rfresults <- cbind(rfresults,apply(as.matrix(partialresults),1,mean))
  }
  rfresults<-rfresults[,-1];colnames(rfresults)<-paste(scenarios)
  
  write.table(rfresults, file = file.path("results",paste("rf","R",R,"size",size,".txt",sep="")))
  
  save.image(file=file.path("results",paste("rf","R",R," size",size," par",par,".Rdata",sep="")))
  
  #boosted trees
  twresults <- matrix(0,length(perfmetrics),1)
  for(i in 1:length(scenarios)){
    pb <- txtProgressBar(min = 0, max = length(scenarios), style = 3) 
    partialresults<-matrix(0,length(perfmetrics),1)
    for(j in 1:R)
    {
      partialresults <-partialresults+try(as.matrix(simplify2array(
        funsim(
          covfun(size,substr(scenarios[i],1,1)),"gbmtwang",par))))
      setTxtProgressBar(pb, i)  
    }
    twresults<-cbind(twresults,partialresults/R)
  }
  twresults<-twresults[,-1];colnames(twresults)<-paste(scenarios)
  
  
  write.table(twresults, file = file.path("results",paste("twang","R",R,"size",size,".txt",sep="")))
  
  #### save workspace after gbm twang
  
  save.image(file=file.path("results",paste("R",R," size",size," par",par,".Rdata",sep="")))
  
  }   
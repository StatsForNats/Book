### FUNCTION TO CALCULATE THE WAIC OF A STAN MODEL

waic_sn<-function(x){
  log_lik <- extract (x, "log_lik")$log_lik
  lppd<-sum(apply(log_lik,2,function(x) log(mean(exp(x)))))
  differ<-sweep(log_lik,2,colMeans(log_lik),FUN="-")
  pwaic <- sum(colMeans(differ^2)*nrow(log_lik)/(nrow(log_lik)-1))
  waic<- -2*(lppd-pwaic)
  return(waic)
}

### FUNCTION TO CREATE A WAIC COMPARATIVE TABLE

WAICtab<-function(...){
  ModelList<-list(...)
  mnames <- match.call()
  mnames <- as.character(mnames)[2:(length(ModelList) + 1)]
  waicvals<-c()
  for(i in 1:length(ModelList)){
    waicval<-waic_sn(ModelList[[i]])
    waicvals<-c(waicvals,waicval)
  }
  ord<-order(waicvals)
  ModCompare<-array(NA,dim = c(length(ModelList),3),dimnames=list(c(),c("WAIC","dWAIC","weight")))
  rownames(ModCompare)<-mnames[ord]
  ModCompare[,1]<-round(waicvals[ord],2)
  ModCompare[,2]<-round(waicvals[ord]-min(waicvals),2)
  TotWeigh<-sum(sapply(as.numeric(ModCompare[,2]),function(x)exp(-0.5*x)))
  for(i in 2:length(ModelList)){
    ModCompare[i,3]<-(exp(-0.5*ModCompare[i,2]))/TotWeigh
  }
  ModCompare[1,3]<-1-sum(ModCompare[2:length(ModelList),3])
  ModCompare[,3]<-round(ModCompare[,3],2)
  return(ModCompare)
}
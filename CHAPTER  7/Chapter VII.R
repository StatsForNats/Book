#####################################################################
########        STATISTICAL MODELING FOR NATURALISTS       #########
########                    CHAPTER VII                   #########
########    QUINTANA-ASCENCIO, LOPEZ BORGHESI, MENGES    #########
#################################################################

##### SETTING REQUIREMENTS #####

rm(list=ls())  # this line clears the environment - helps prevent confusion

# Call and set up rstan and other necessary libraries.
library(rstan)
rstan_options(auto_write = TRUE) # it saves the compiled instance in our directory
options(mc.cores = parallel::detectCores()) # it runs the chains in parallel

# Set working directory to source file 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


### LOADING AND PREPARING DATA ###

orig_data <- read.csv("Hypericum_cumulicola.csv", header=T)
dd  <- subset(orig_data, orig_data$rep>0 & time<1997) # Subset of the data to be used
dd$lnrep <- log(dd$rep) # Creating log of number of fruits


### DATA VISUALIZATION ###
b<-seq(0,6100,50) # Set breaks for histograms

par(mfrow=c(1,2),mar=c(4,4,2,1))

hist(dd$rep[dd$surv_fin==0],breaks=b,xlim=c(0,1200),ylim=c(0,0.01),main="Survival = 0",
     xlab="Number of fruits",col="gray50",freq=F)

hist(dd$rep[dd$surv_fin==1],breaks=b,xlim=c(0,1200),ylim=c(0,0.01),main="Survival = 1",
     xlab="Number of fruits",col="gray50",freq=F)

### MODEL BUILDING ###

# Gathering the data for Stan

stan_data <- list(surv_fin=dd$surv_fin, rep= dd$rep,lnrep= dd$lnrep, N=nrow(dd) )


## MODEL 7.1 - logistic model in natural scale ##

# Setting the parameters monitored 
params <- c("a","b","dev","log_lik")

# Defining the MCMC settings 
ni <- 6000 	# number of iterations
nt <- 2	# interval for data collection
nb <- 2000	# number of iterations discarded
nc <- 3	# number of chains

# Defining initial values
inits <- lapply(1:nc, function(i) {
  list(a = rnorm(1,0,5), b = runif(1,0,5))})

# Calling for Stan from R. Model6.1
Model_7_1  <- stan("Chapter_VII_natural_scale.stan",
                   data = stan_data, init = inits, pars = params,
                   chains = nc, iter = ni, warmup = nb, thin = nt,
                   seed = 1, open_progress = FALSE)

# Summarizing and plotting the posteriors 
print(Model_7_1, digits = 4)

pairs(Model_7_1,pars=c("a","b"))
traceplot(Model_7_1,pars=c("a","b"))

## MODEL 7.2 - power model ##

# Data, monitored parameters, MCMC settnigs and initial values all remain the same as before

Model_7_2  <- stan("Chapter_VII_power_model.stan",
                   data = stan_data, init = inits, pars = params,
                   chains = nc, iter = ni, warmup = nb, thin = nt,
                   seed = 1, open_progress = FALSE)

# Summarizing and plotting the posteriors 
print(Model_7_2, digits = 4,pars=c("a","b"))

pairs(Model_7_2,pars=c("a","b"))
traceplot(Model_7_2,pars=c("a","b"))

###################################################################

### COMPARING MODELS ###

## CREATE FUNCTIONS ##
# First, we will generate a function to calculate the WAIC of a model (waic_sn)
# and another function to make a table comparing the WAIC of each model (WAICtab) 
# these functions can also be found in the "extra material" folder for future reference

waic_sn<-function(x){
  log_lik <- extract (x, "log_lik")$log_lik
  lppd<-sum(apply(log_lik,2,function(x) log(mean(exp(x)))))
  differ<-sweep(log_lik,2,colMeans(log_lik),FUN="-")
  pwaic <- sum(colMeans(differ^2)*nrow(log_lik)/(nrow(log_lik)-1))
  waic<- -2*(lppd-pwaic)
  return(waic)
}

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

## USE WAICtab TO COMPARE BOTH MODELS ##

WAICtab(Model_7_1,Model_7_2)

#####################################

### PLOTING MODELS 7.1 AND 7.2 ###

## CALCULATING MODEL PREDICTIONS ##

x_pred<-seq(min(dd$rep),max(dd$rep),1) # Series of x values for prediction
lnx_pred<-log(x_pred)  # log of x

# Model 1:
post1 <- extract(Model_7_1) #Gets all values of parameters (6000 iterations)

Mu1<- (exp(rep(post1$a,length(x_pred))+outer(post1$b,x_pred,FUN = '*')))/
  (1+exp(rep(post1$a,length(x_pred))+outer(post1$b,x_pred,FUN = '*'))) # We calculate P(survival) using the link function

mu.mean1 <- apply(Mu1,2,mean)  #Calculate mean prediction
mu.PI1U <- mu.mean1+(qt(0.975,nrow(Mu1)-1)*apply(Mu1,2,sd)) # Calculate upper CI values
mu.PI1L <- mu.mean1+(qt(0.025,nrow(Mu1)-1)*apply(Mu1,2,sd)) # Calculate lower CI values

# Model 2:
post2 <- extract(Model_7_2) 

Mu2<- (exp(rep(post2$a,length(lnx_pred))+outer(post2$b,lnx_pred,FUN = '*')))/
  (1+exp(rep(post2$a,length(lnx_pred))+outer(post2$b,lnx_pred,FUN = '*')))

mu.mean2 <- apply(Mu2,2,mean)
mu.PI2U <- mu.mean2+(qt(0.975,nrow(Mu2)-1)*apply(Mu2,2,sd))
mu.PI2L <- mu.mean2+(qt(0.025,nrow(Mu2)-1)*apply(Mu2,2,sd))

## PLOT ##

par(mfrow=c(1,1),mar=c(4,4,2,1))
plot(dd$rep, dd$surv_fin,type="n",ylab="P(survival)",xlab="Number of fruits",
     xlim=c(1,1500),ylim=c(0.2,0.8))

b <-c(seq(0,150,10),seq(200,500,50),seq(1000,6500,200))
z <- cut(dd$rep, b)
prebyden <- tapply(dd$surv_fin,z,sum)
tab <- table(z)
probs <-prebyden/tab
probs <- as.vector(probs)
for (i in 1:length(b-1)){
  points(b[i],probs[i],pch=1,cex=(prebyden[i]+1)^(0.4),col="black")}

lines(x_pred,mu.mean1,col=rgb(0.4,0,0,1),lwd=3)
polygon(c(x_pred,rev(x_pred)),c(mu.PI1U,rev(mu.PI1L)),border=NA,col=rgb(0.8,0,0,0.3))

lines(x_pred,mu.mean2,col=rgb(0,0,0.4,1),lwd=3)
polygon(c(x_pred,rev(x_pred)),c(mu.PI2U,rev(mu.PI2L)),border=NA,col=rgb(0,0,0.7,0.3))

legend(1000,0.8,c("Model 7.1","Model 7.2"),bty="n",lty=1,col=c(rgb(0.4,0,0,1),rgb(0,0,0.4,1)),
       lwd=2,x.intersp=0.6,y.intersp=0.8,seg.len=0.6)

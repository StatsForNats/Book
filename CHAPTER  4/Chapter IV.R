#####################################################################
########        STATISTICAL MODELING FOR NATURALISTS       #########
########                      CHAPTER IV                  #########
########    QUINTANA-ASCENCIO, LOPEZ BORGHESI, MENGES    #########
#################################################################

##### SETTING REQUIREMENTS #####

rm(list=ls())  # this line clears the environment - helps prevent confusion

# Call rstan library to allow communication with 'stan' and set up
library(rstan)
rstan_options(auto_write = TRUE) # it saves the compiled instance in our directory
options(mc.cores = parallel::detectCores()) # it runs the chains in parallel

# Set working directory to source file 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# DATA RETRIEVAL AND VISUALIZATIOn
Hc_data <- read.csv("Hypericum_cumulicola.csv", header=T)

hgts<-na.omit(Hc_data$init_height[Hc_data$time==2014 & Hc_data$stage != "sg"]) # vector of heights
# NOTE: we use na.omit to avoid accidentally selecting NA values in the subsample. Stan does work with NAs

hist(hgts,25,main="Histogram of Hypericum cumulicola height (cm)")
segments(mean(hgts),0,mean(hgts),110, col="black",lwd=3)


##### MODEL BUILDING  #####

n<-10  # determines sample size
x <- sample(hgts,n) # create subsample of n = 10

stan_data <- list(height=x, N=length(x) ) # put data in format useful for stan

# 1. MODEL WITH DIFFUSE PRIORS

## Define parameters monitored
params <- c("a","sigma")

## MCMC settings
ni <- 5000 # total iterations per chain
nt <- 1    # interval of data collection
nb <- 1000 # burned iterations - first 1000 are discarded
nc <- 3 # number of mcmc chains

## Initial values
inits <- lapply(1:nc, function(i) {
  list(a = rnorm(1,0,30), sigma = runif(1,0,15))})

## Call Stan from R
chap3dif  <- stan("Stan averages difuse.stan",
                  data = stan_data, init = inits, pars = params,
                  chains = nc, iter = ni, warmup = nb, thin = nt,
                  seed = 1,
                  open_progress = FALSE)

### Visual assessment of results
pairs(chap3dif,pars=c("a","sigma")) # Shows distributions of each parameter and their correlation

traceplot(chap3dif) # Shows how parameter estimation change from iteration to iteration (ignores burned iterations)

## Obtain summary of model
print(chap3dif, digits = 3) # provides mean values and estimates of uncertainty for each parameter
# (a and sigma) and loglikelihood (lp_)

posteriors<-extract(chap3dif) # extract posteriors


# 2. MODEL WITH INFORME PRIORS

# NOTE: The definition of parameters, MCMC settings, and initial values remain the same

## Call Stan from R
chap3inf  <- stan("Stan averages informed.stan",
                  data = stan_data, init = inits, pars = params,
                  chains = nc, iter = ni, warmup = nb, thin = nt,
                  seed = 1,
                  open_progress = FALSE)


### Visual assessment of results
pairs(chap3inf,pars=c("a","sigma"))
traceplot(chap3inf)

## Obtain summary of model
print(chap3inf, digits = 3)

posteriors2<-extract(chap3inf) 

# Checking values
mean(x) - qt(0.975, n-1)*sqrt(var(x)/n)
mean(x) + qt(0.975, n-1)*sqrt(var(x)/n)

##### COMPARING PRIORS AND POSTERIORS #####

pop_mean <- tapply(Hc_data$init_height[Hc_data$stage != "sg" & !is.na(Hc_data$init_height) & Hc_data$time < 1997 ],
                   Hc_data$site[Hc_data$stage != "sg" & !is.na(Hc_data$init_height) & Hc_data$time < 1997 ],mean)

pop.mean.mean <- mean(pop_mean)   # informed prior mean
pop.mean.sd <- sd(pop_mean)       # informed prior standard deviation

overmean<-mean(na.omit(Hc_data$init_height[Hc_data$stage != "sg"],))
freqmean<-mean(x)

avAdiff<-summary(chap3dif)$summary[1,1]   # posterior mean (diffuse priors model)
sdAdiff<-summary(chap3dif)$summary[1,3]   # posterior standard deviation (diffuse priors model)
avAinf<-summary(chap3inf)$summary[1,1]    # posterior mean (informed priors model)
sdAinf<-summary(chap3inf)$summary[1,3]    # posterior standard deviation (informed priors model)


PostAdiff<-density(rnorm(1000000,avAdiff,sdAdiff))
PostAinf<-density(rnorm(1000000,avAinf,sdAinf))
PriorAinf<-density(rnorm(1000000,pop.mean.mean,pop.mean.sd))

par(mfrow=c(1,1))
plot(NULL,xlim=c(15,50),ylim=c(0,0.22),cex.main=1.8,cex.axis=1.6,cex.lab=1.8,ylab="Frequency",
     mgp=c(2.4,0.7,0),xlab="Height (cm)",main="Distributions of the mean")

lines(PriorAinf,lwd=3,col="red")
lines(PostAdiff,lwd=3,col="blue")
lines(PostAinf,lwd=3,col="darkgreen")

segments(pop.mean.mean,0,pop.mean.mean,0.18,lwd=3,col="gray65")
segments(overmean,-0,overmean,0.18,lwd=3)


legend(36,0.21,c("Prior (informed)","Posterior (diffused)","Posterior (informed)",
                 "Overall mean", "Population mean"),bty="n",lwd=2,
       col=c("red","blue","darkgreen","black","gray65"),
       cex=1.3,x.intersp=0.6,y.intersp=0.6)


############ LOOPS FOR MULTIPLE SAMPLES   #################
# IN LINES 137 TO 279 WE REPEAT THE ABOVE PROCESS 10 TIMES FOR RANDOM SAMPLES OF N = 10 AND
# FOR 10 RANDOM SAMPLES OF N = 100

### SAMPLES OF 10
size <- 10

rawdata<-array(0,c(size,10)) #creates an place to save all the sampling events
colnames(rawdata)<-c("Event1","Event2","Event3","Event4","Event5","Event6","Event7","Event8","Event9","Event10")
rownames(rawdata)<-c(1:size)
rawdata

for (i in 1:10){
  x <- sample(na.omit(Hc_data$init_height[Hc_data$time==2014 & Hc_data$stage != "sg"]),size)
  rawdata[,i]<-x
}

BDResults10 <- array(0,c(10,8))
colnames(BDResults10) <- c("Mean", "StdDev", "2.5%", "97.5%", "Sigma", "Stdev", "2.5%", "97.5%")
rownames(BDResults10)<-c(1:10)
BDResults10

BIResults10 <- array(0,c(10,8))
colnames(BIResults10) <- c("Mean", "StdDev", "2.5%", "97.5%", "Sigma", "Stdev", "2.5%", "97.5%")
rownames(BIResults10)<-c(1:10)
BIResults10

for (i in 1:10){
  ## Parameters monitored
  params <- c("a","sigma")
  
  ## MCMC settings
  ni <- 5000
  nt <- 1
  nb <- 1000
  nc <- 3
  
  ## Initial values
  inits <- lapply(1:nc, function(i) {
    list(a = rnorm(1,0,30), sigma = runif(1,0,15))})
  stan_data <- list(height=rawdata[,i], N=length(rawdata[,i]))
  ## Call Stan from R - diffuse
  loop10diff  <- stan("Stan averages difuse.stan",
                      data = stan_data, init = inits, pars = params,
                      chains = nc, iter = ni, warmup = nb, thin = nt,
                      seed = 1,
                      open_progress = FALSE)
  ## Call Stan from R - informed
  loop10inf  <- stan("Stan averages informed.stan",
                     data = stan_data, init = inits, pars = params,
                     chains = nc, iter = ni, warmup = nb, thin = nt,
                     seed = 1,
                     open_progress = FALSE)
  ### Saving values from each loop
  BDResults10[i,1]<-summary(loop10diff)$summary[1,1]
  BDResults10[i,2]<-summary(loop10diff)$summary[1,3]
  BDResults10[i,3]<-summary(loop10diff)$summary[1,4]
  BDResults10[i,4]<-summary(loop10diff)$summary[1,8]
  BDResults10[i,5]<-summary(loop10diff)$summary[2,1]
  BDResults10[i,6]<-summary(loop10diff)$summary[2,3]
  BDResults10[i,7]<-summary(loop10diff)$summary[2,4]
  BDResults10[i,8]<-summary(loop10diff)$summary[2,8]
  
  BIResults10[i,1]<-summary(loop10inf)$summary[1,1]
  BIResults10[i,2]<-summary(loop10inf)$summary[1,3]
  BIResults10[i,3]<-summary(loop10inf)$summary[1,4]
  BIResults10[i,4]<-summary(loop10inf)$summary[1,8]
  BIResults10[i,5]<-summary(loop10inf)$summary[2,1]
  BIResults10[i,6]<-summary(loop10inf)$summary[2,3]
  BIResults10[i,7]<-summary(loop10inf)$summary[2,4]
  BIResults10[i,8]<-summary(loop10inf)$summary[2,8]
}

write.csv(BDResults10,'BDResults10.csv') # SAVE RESULTS OF MODELS WITH SAMPLES N = 10 AND DIFFUSED PRIORS
write.csv(BIResults10,'BIResults10.csv') # SAVE RESULTS OF MODELS WITH SAMPLES N = 10 AND INFORMED PRIORS


### SAMPLES OF 100
size <- 100

rawdata<-array(0,c(size,10)) #creates an place to save all the sampling events
colnames(rawdata)<-c("Event1","Event2","Event3","Event4","Event5","Event6","Event7","Event8","Event9","Event10")
rownames(rawdata)<-c(1:size)
rawdata

for (i in 1:10){
  x <- sample(na.omit(Hc_data$init_height[Hc_data$time==2014 & Hc_data$stage != "sg"]),size)
  rawdata[,i]<-x
}

BDResults100 <- array(0,c(10,8))
colnames(BDResults100) <- c("Mean", "StdDev", "2.5%", "97.5%", "Sigma", "Stdev", "2.5%", "97.5%")
rownames(BDResults100)<-c(1:10)
BDResults100

BIResults100 <- array(0,c(10,8))
colnames(BIResults100) <- c("Mean", "StdDev", "2.5%", "97.5%", "Sigma", "Stdev", "2.5%", "97.5%")
rownames(BIResults100)<-c(1:10)
BIResults100

for (i in 1:10){
  ## Parameters monitored
  params <- c("a","sigma")
  
  ## MCMC settings
  ni <- 5000
  nt <- 1
  nb <- 1000
  nc <- 3
  
  ## Initial values
  inits <- lapply(1:nc, function(i) {
    list(a = rnorm(1,0,30), sigma = runif(1,0,15))})
  stan_data <- list(height=rawdata[,i], N=length(rawdata[,i]))
  ## Call Stan from R - diffuse
  loop100diff  <- stan("Stan averages difuse.stan",
                       data = stan_data, init = inits, pars = params,
                       chains = nc, iter = ni, warmup = nb, thin = nt,
                       seed = 1,
                       open_progress = FALSE)
  ## Call Stan from R - informed
  loop100inf  <- stan("Stan averages informed.stan",
                      data = stan_data, init = inits, pars = params,
                      chains = nc, iter = ni, warmup = nb, thin = nt,
                      seed = 1,
                      open_progress = FALSE)
  ### Saving values from each loop
  BDResults100[i,1]<-summary(loop100diff)$summary[1,1]
  BDResults100[i,2]<-summary(loop100diff)$summary[1,3]
  BDResults100[i,3]<-summary(loop100diff)$summary[1,4]
  BDResults100[i,4]<-summary(loop100diff)$summary[1,8]
  BDResults100[i,5]<-summary(loop100diff)$summary[2,1]
  BDResults100[i,6]<-summary(loop100diff)$summary[2,3]
  BDResults100[i,7]<-summary(loop100diff)$summary[2,4]
  BDResults100[i,8]<-summary(loop100diff)$summary[2,8]
  
  BIResults100[i,1]<-summary(loop100inf)$summary[1,1]
  BIResults100[i,2]<-summary(loop100inf)$summary[1,3]
  BIResults100[i,3]<-summary(loop100inf)$summary[1,4]
  BIResults100[i,4]<-summary(loop100inf)$summary[1,8]
  BIResults100[i,5]<-summary(loop100inf)$summary[2,1]
  BIResults100[i,6]<-summary(loop100inf)$summary[2,3]
  BIResults100[i,7]<-summary(loop100inf)$summary[2,4]
  BIResults100[i,8]<-summary(loop100inf)$summary[2,8]
}

write.csv(BDResults100,'BDResults100.csv') # SAVE RESULTS OF MODELS WITH SAMPLES N = 100 AND DIFFUSED PRIORS
write.csv(BIResults100,'BIResults100.csv') # SAVE RESULTS OF MODELS WITH SAMPLES N = 100 AND INFORMED PRIORS

#########################################

### LOAD DATA CREATED IN THE LOOPS ABOVE

BDdata10<-read.csv("BDResults10.csv",header = T)
BIdata10<-read.csv("BIResults10.csv",header = T)
BDdata100<-read.csv("BDResults100.csv",header = T)
BIdata100<-read.csv("BIResults100.csv",header = T)

#### PLOT MULTIPLE RUNS 

### MEANS: PLOT POSTERIOR OF PARAMETER "A" OF MODEL FOR ALL SAMPLES AND BOTH PRIORS

par(fig=c(0,50,50,95)/100,mar=c(5,4,2,1))
plot(NULL,xlim=c(15,50),ylim=c(0,0.40),cex.main=1.9,cex.axis=1.7,cex.lab=1.9,ylab="Probability",
     mgp=c(2.4,0.7,0),xlab="Height (cm)",main="Diffuse priors (N = 10)",yaxs="i")
for(i in 1:10){
  diff10<-density(rnorm(1000000,BDdata10$Mean[i],BDdata10$StdDev[i]))
  lines(diff10$x,diff10$y,col="red")
}
segments(pop.mean.mean,0,pop.mean.mean,0.36,lwd=4,lty=1)
segments(overmean,0,overmean,0.36,lwd=4,lty=2)


par(fig=c(50,100,50,95)/100,mar=c(5,4,2,1),new=T)
plot(NULL,xlim=c(15,50),ylim=c(0,0.40),cex.main=1.9,cex.axis=1.7,cex.lab=1.9,ylab="Probability",
     mgp=c(2.4,0.7,0),xlab="Height (cm)",main="Diffuse priors (N = 100)",yaxs="i")

for(i in 1:10){
  diff10<-density(rnorm(1000000,BDdata100$Mean[i],BDdata100$StdDev[i]))
  lines(diff10$x,diff10$y,col="blue")
}
segments(pop.mean.mean,0,pop.mean.mean,0.36,lwd=4,lty=1)
segments(overmean,0,overmean,0.36,lwd=4,lty=2)

par(fig=c(0,50,5,50)/100,mar=c(5,4,2,1),new=T)
plot(NULL,xlim=c(15,50),ylim=c(0,0.40),cex.main=1.9,cex.axis=1.7,cex.lab=1.9,ylab="Probability",
     mgp=c(2.4,0.7,0),xlab="Height (cm)",main="Informed priors (N = 10)",yaxs="i")

for(i in 1:10){
  diff10<-density(rnorm(1000000,BIdata10$Mean[i],BIdata10$StdDev[i]))
  lines(diff10$x,diff10$y,col="red")
}
segments(pop.mean.mean,0,pop.mean.mean,0.36,lwd=4,lty=1)
segments(overmean,0,overmean,0.36,lwd=4,lty=2)

par(fig=c(50,100,5,50)/100,mar=c(5,4,2,1),new=T)
plot(NULL,xlim=c(15,50),ylim=c(0,0.40),cex.main=1.9,cex.axis=1.7,cex.lab=1.9,ylab="Probability",
     mgp=c(2.4,0.7,0),xlab="Height (cm)",main="Informed priors (N = 100)",yaxs="i")

for(i in 1:10){
  diff10<-density(rnorm(1000000,BIdata100$Mean[i],BIdata100$StdDev[i]))
  lines(diff10$x,diff10$y,col="blue")
}
segments(pop.mean.mean,0,pop.mean.mean,0.36,lwd=4,lty=1)
segments(overmean,0,overmean,0.36,lwd=4,lty=2)

############################

### VARIANCE: PLOT POSTERIOR OF PARAMETER "SIGMA" OF MODEL FOR ALL SAMPLES AND BOTH PRIORS

oversd<-sqrt(var(na.omit(Hc_data$init_height[Hc_data$stage != "sg"],)))

par(fig=c(0,50,50,95)/100,mar=c(5,4,2,1))
plot(NULL,xlim=c(0,20),ylim=c(0,0.9),cex.main=1.9,cex.axis=1.7,cex.lab=1.9,ylab="Probability",
     mgp=c(2.4,0.7,0),xlab="sd of variance in height (cm)",main="Diffuse priors (N = 10)",yaxs="i")
for(i in 1:10){
  diff10<-density(rnorm(1000000,BDdata10$Sigma[i],BDdata10$Stdev[i]))
  lines(diff10$x,diff10$y,col="red")
}
segments(sdpop_var,0,sdpop_var,0.85,lwd=4,lty=1)
segments(oversd,0,oversd,0.85,lwd=4,lty=2)


par(fig=c(50,100,50,95)/100,mar=c(5,4,2,1),new=T)
plot(NULL,xlim=c(0,20),ylim=c(0,0.9),cex.main=1.9,cex.axis=1.7,cex.lab=1.9,ylab="Probability",
     mgp=c(2.4,0.7,0),xlab="sd of variance in height (cm)",main="Diffuse priors (N = 100)",yaxs="i")

for(i in 1:10){
  diff10<-density(rnorm(1000000,BDdata100$Sigma[i],BDdata100$Stdev[i]))
  lines(diff10$x,diff10$y,col="blue")
}
segments(sdpop_var,0,sdpop_var,0.85,lwd=4,lty=1)
segments(oversd,0,oversd,0.85,lwd=4,lty=2)


par(fig=c(0,50,5,50)/100,mar=c(5,4,2,1),new=T)
plot(NULL,xlim=c(0,20),ylim=c(0,0.9),cex.main=1.9,cex.axis=1.7,cex.lab=1.9,ylab="Probability",
     mgp=c(2.4,0.7,0),xlab="sd of variance in height (cm)",main="Informed priors (N = 10)",yaxs="i")

for(i in 1:10){
  diff10<-density(rnorm(1000000,BIdata10$Sigma[i],BIdata10$Stdev[i]))
  lines(diff10$x,diff10$y,col="red")
}
segments(sdpop_var,0,sdpop_var,0.85,lwd=4,lty=1)
segments(oversd,0,oversd,0.85,lwd=4,lty=2)


par(fig=c(50,100,5,50)/100,mar=c(5,4,2,1),new=T)
plot(NULL,xlim=c(0,20),ylim=c(0,0.9),cex.main=1.9,cex.axis=1.7,cex.lab=1.9,ylab="Probability",
     mgp=c(2.4,0.7,0),xlab="sd of variance in height (cm)",main="Informed priors (N = 100)",yaxs="i")

for(i in 1:10){
  diff10<-density(rnorm(1000000,BIdata100$Sigma[i],BIdata100$Stdev[i]))
  lines(diff10$x,diff10$y,col="blue")
}
segments(sdpop_var,0,sdpop_var,0.85,lwd=4,lty=1)
segments(oversd,0,oversd,0.85,lwd=4,lty=2)


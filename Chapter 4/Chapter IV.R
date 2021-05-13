#################################################
########             CHAPTER IV        #########
###############################################

##### SETTING REQUIREMENTS #####

rm(list=ls())  # this line clears the environment - helps prevent confusion

# Call rstan library to allow communication with 'stan' and set up
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Set working directory to source file 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data and obtain random sub-sample
Hc_data <- read.csv("Hypericum_cumulicola.csv", header=T)

n<-10  # determines sample size

hgts<-na.omit(Hc_data$init_height[Hc_data$time==2014 & Hc_data$stage != "sg"]) # vector of heights
# NOTE: we use na.omit to avoid accidentally selecting NA values in the 
# subsample. Stan does work with NAs

x <- sample(hgts,n) # create subsample of n = 10

stan_data <- list(height=x, N=length(x) ) # put data in format useful for stan


##### MODEL BUILDING  #####

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

pop.mean.mean <- mean(pop_mean)
overmean<-mean(na.omit(Hc_data$init_height[Hc_data$stage != "sg"],))
freqmean<-mean(x)

avAdiff<-summary(chap3dif)$summary[1,1]
sdAdiff<-summary(chap3dif)$summary[1,3]
avAinf<-summary(chap3inf)$summary[1,1]
sdAinf<-summary(chap3inf)$summary[1,3]


PostAdiff<-density(rnorm(1000000,avAdiff,sdAdiff))
PostAinf<-density(rnorm(1000000,avAinf,sdAinf))
PriorAinf<-density(rnorm(1000000,pop.mean.mean,pop.mean.sd))

par(mfrow=c(1,1))
plot(NULL,xlim=c(15,50),ylim=c(-0.04,0.20),cex.main=1.8,cex.axis=1.6,cex.lab=1.8,ylab="Frequency",
           mgp=c(2.4,0.7,0),xlab="Height (cm)",main="Distributions of the mean")

lines(PriorAinf,lwd=6,col="gray65")
lines(PostAdiff,lwd=6,lty=3,col="gray25")
lines(PostAinf,lwd=6,lty=2,col="gray15")
abline(h=-0.005)

segments(pop.mean.mean,-0.01,pop.mean.mean,-0.044,lwd=4,col="gray65")
segments(overmean,-0.01,overmean,-0.044,lwd=4)
segments(freqmean,-0.01,freqmean,-0.044,lwd=4,lty=3)


legend(32,0.21,c("Prior (informed)","Posterior (diffused)","Posterior (informed)"),
       bty="n",lty=c(1,3,2),lwd=3,col=c("gray65","gray25","gray15"),
       cex=1.6,x.intersp=0.6,y.intersp=0.6)

#### DIFFUSE

PostA<-as.data.frame(extract(chap3dif,permuted=FALSE,inc_warmup=TRUE)[,,1])
PostSigma<-as.data.frame(extract(chap3dif,permuted=FALSE,inc_warmup=TRUE)[,,2])


### TRACEPLOT WITH CHAINS SEPARATED
## MEAN

par(mfrow=c(2,3),mar=c(4,5,2,2))
plot(1:1000,PostA$`chain:1`[1:1000],type="l",ylim=c(-10,60),xlim=c(0,5000),
     ylab="Mean (a)",xlab="",main="chain 1",col="gray60",cex.main=2.8,cex.axis=2.2,
     cex.lab=2.8,mgp=c(3,1,0))
lines(1001:5000,PostA$`chain:1`[1001:5000],col="gray25")
abline(v=1000,lty=3,lwd=3)
plot(1:1000,PostA$`chain:2`[1:1000],type="l",ylim=c(-10,60),xlim=c(0,5000),
     ylab="",xlab="",main="chain 2",col="gray60",cex.main=2.8,cex.axis=2.2,
     cex.lab=2.8,mgp=c(3,1,0))
lines(1001:5000,PostA$`chain:2`[1001:5000],col="gray25")
abline(v=1000,lty=3,lwd=3)
plot(1:1000,PostA$`chain:3`[1:1000],type="l",ylim=c(-10,60),xlim=c(0,5000),
     ylab="",xlab="",main="chain 3",col="gray60",cex.main=2.8,cex.axis=2.2,
     cex.lab=2.8,mgp=c(3,1,0))
lines(1001:5000,PostA$`chain:3`[1001:5000],col="gray30")
abline(v=1000,lty=3,lwd=3)
## SIGMA

plot(1:1000,PostSigma$`chain:1`[1:1000],type="l",ylim=c(0,50),xlim=c(0,5000),
     ylab="Variance (sigma)",xlab="",main="chain 1",col="gray60",cex.main=2.8,cex.axis=2.2,
     cex.lab=2.8,mgp=c(3,1,0))
lines(1001:5000,PostSigma$`chain:1`[1001:5000],col="gray25")
abline(v=1000,lty=3,lwd=3)
plot(1:1000,PostSigma$`chain:2`[1:1000],type="l",ylim=c(0,50),xlim=c(0,5000),
     ylab="",xlab="",main="chain 2",col="gray60",cex.main=2.8,cex.axis=2.2,
     cex.lab=2.8,mgp=c(3,1,0))
lines(1001:5000,PostSigma$`chain:2`[1001:5000],col="gray25")
abline(v=1000,lty=3,lwd=3)
plot(1:1000,PostSigma$`chain:3`[1:1000],type="l",ylim=c(0,50),xlim=c(0,5000),
     ylab="",xlab="",main="chain 3",col="gray60",cex.main=2.8,cex.axis=2.2,
     cex.lab=2.8,mgp=c(3,1,0))
lines(1001:5000,PostSigma$`chain:3`[1001:5000],col="gray25")
abline(v=1000,lty=3,lwd=3)

#### TRACEPLOT
par(mfrow=c(1,2))
plot(PostA$`chain:1`,type="l",ylim=c(-5,55),xlim=c(0,5000),
     ylab="Mean",xlab="",main="Mean (a)",col="gray20")
lines(PostA$`chain:2`,col="gray40")
lines(PostA$`chain:3`,col="gray60")
abline(v=1000,lty=3)

plot(PostSigma$`chain:1`,type="l",ylim=c(5,50),xlim=c(0,5000),
     ylab="Mean",xlab="",main="Sigma",col="gray20")
lines(PostSigma$`chain:2`,col="gray40")
lines(PostSigma$`chain:3`,col="gray60")
abline(v=1000,lty=3)

### PAIR PLOTS

allA<-c(PostA[1001:5000,1],PostA[1001:5000,2],PostA[1001:5000,3])
allSigma<-c(PostSigma[1001:5000,1],PostSigma[1001:5000,2],PostSigma[1001:5000,3])
      plot(allA,allSigma,cex=0.6,col=alpha("gray30",0.25),pch=19)
plot(density(allA))
polygon(density(allA)$x,density(allA)$y,col="red")


############ LOOPS FOR MULTIPLE SAMPLES   #################

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

write.csv(BDResults10,'BDResults10.csv')
write.csv(BIResults10,'BIResults10.csv')


### SAMPLES OF 10
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

write.csv(BDResults100,'BDResults100.csv')
write.csv(BIResults100,'BIResults100.csv')

#########################################

BDdata10<-read.csv("BDResults10.csv",header = T)
BIdata10<-read.csv("BIResults10.csv",header = T)
BDdata100<-read.csv("BDResults100.csv",header = T)
BIdata100<-read.csv("BIResults100.csv",header = T)

#### PLOT MULTIPLE RUNS 

### MEANS

par(fig=c(0,50,50,95)/100,mar=c(5,4,2,1))
plot(NULL,xlim=c(15,50),ylim=c(0,0.40),cex.main=1.9,cex.axis=1.7,cex.lab=1.9,ylab="Probability",
     mgp=c(2.4,0.7,0),xlab="Height (cm)",main="Diffuse priors (N = 10)",yaxs="i")
for(i in 1:10){
diff10<-density(rnorm(1000000,BDdata10[i,1],BDdata10[i,2]))
polygon(diff10$x,diff10$y,col=alpha("gray65",0.5),border=alpha("gray65",0.5))
lines(diff10$x,diff10$y,col="gray50")
}
segments(pop.mean.mean,0,pop.mean.mean,0.36,lwd=4,lty=1)
segments(overmean,0,overmean,0.36,lwd=4,lty=2)


par(fig=c(50,100,50,95)/100,mar=c(5,4,2,1),new=T)
plot(NULL,xlim=c(15,50),ylim=c(0,0.40),cex.main=1.9,cex.axis=1.7,cex.lab=1.9,ylab="Probability",
     mgp=c(2.4,0.7,0),xlab="Height (cm)",main="Diffuse priors (N = 100)",yaxs="i")

for(i in 1:10){
  diff10<-density(rnorm(1000000,BDdata100[i,1],BDdata100[i,2]))
  polygon(diff10$x,diff10$y,col=alpha("gray65",0.5),border=alpha("gray65",0.5))
  lines(diff10$x,diff10$y,col="gray50")
}
segments(pop.mean.mean,0,pop.mean.mean,0.36,lwd=4,lty=1)
segments(overmean,0,overmean,0.36,lwd=4,lty=2)

par(fig=c(0,50,5,50)/100,mar=c(5,4,2,1),new=T)
plot(NULL,xlim=c(15,50),ylim=c(0,0.40),cex.main=1.9,cex.axis=1.7,cex.lab=1.9,ylab="Probability",
     mgp=c(2.4,0.7,0),xlab="Height (cm)",main="Informed priors (N = 10)",yaxs="i")

for(i in 1:10){
  diff10<-density(rnorm(1000000,BIdata10[i,1],BIdata10[i,2]))
  polygon(diff10$x,diff10$y,col=alpha("gray65",0.5),border=alpha("gray65",0.5))
  lines(diff10$x,diff10$y,col="gray50")
}
segments(pop.mean.mean,0,pop.mean.mean,0.36,lwd=4,lty=1)
segments(overmean,0,overmean,0.36,lwd=4,lty=2)

par(fig=c(50,100,5,50)/100,mar=c(5,4,2,1),new=T)
plot(NULL,xlim=c(15,50),ylim=c(0,0.40),cex.main=1.9,cex.axis=1.7,cex.lab=1.9,ylab="Probability",
     mgp=c(2.4,0.7,0),xlab="Height (cm)",main="Informed priors (N = 100)",yaxs="i")

for(i in 1:10){
  diff10<-density(rnorm(1000000,BIdata100[i,1],BIdata100[i,2]))
  polygon(diff10$x,diff10$y,col=alpha("gray65",0.5),border=alpha("gray65",0.5))
  lines(diff10$x,diff10$y,col="gray50")
}
segments(pop.mean.mean,0,pop.mean.mean,0.36,lwd=4,lty=1)
segments(overmean,0,overmean,0.36,lwd=4,lty=2)

############################

oversd<-sqrt(var(na.omit(Hc_data$init_height[Hc_data$stage != "sg"],)))

### VARIANCE

par(fig=c(0,50,50,95)/100,mar=c(5,4,2,1))
plot(NULL,xlim=c(0,20),ylim=c(0,0.9),cex.main=1.9,cex.axis=1.7,cex.lab=1.9,ylab="Probability",
     mgp=c(2.4,0.7,0),xlab="sd of variance in height (cm)",main="Diffuse priors (N = 10)",yaxs="i")
for(i in 1:10){
  diff10<-density(rnorm(1000000,BDdata10[i,5],BDdata10[i,6]))
  polygon(diff10$x,diff10$y,col=alpha("gray65",0.5),border=alpha("gray65",0.5))
  lines(diff10$x,diff10$y,col="gray50")
}
segments(sdpop_var,0,sdpop_var,0.85,lwd=4,lty=1)
segments(oversd,0,oversd,0.85,lwd=4,lty=2)


par(fig=c(50,100,50,95)/100,mar=c(5,4,2,1),new=T)
plot(NULL,xlim=c(0,20),ylim=c(0,0.9),cex.main=1.9,cex.axis=1.7,cex.lab=1.9,ylab="Probability",
     mgp=c(2.4,0.7,0),xlab="sd of variance in height (cm)",main="Diffuse priors (N = 100)",yaxs="i")

for(i in 1:10){
  diff10<-density(rnorm(1000000,BDdata100[i,5],BDdata100[i,6]))
  polygon(diff10$x,diff10$y,col=alpha("gray65",0.5),border=alpha("gray65",0.5))
  lines(diff10$x,diff10$y,col="gray50")
}
segments(sdpop_var,0,sdpop_var,0.85,lwd=4,lty=1)
segments(oversd,0,oversd,0.85,lwd=4,lty=2)


par(fig=c(0,50,5,50)/100,mar=c(5,4,2,1),new=T)
plot(NULL,xlim=c(0,20),ylim=c(0,0.9),cex.main=1.9,cex.axis=1.7,cex.lab=1.9,ylab="Probability",
     mgp=c(2.4,0.7,0),xlab="sd of variance in height (cm)",main="Informed priors (N = 10)",yaxs="i")

for(i in 1:10){
  diff10<-density(rnorm(1000000,BIdata10[i,5],BIdata10[i,6]))
  polygon(diff10$x,diff10$y,col=alpha("gray65",0.5),border=alpha("gray65",0.5))
  lines(diff10$x,diff10$y,col="gray50")
}
segments(sdpop_var,0,sdpop_var,0.85,lwd=4,lty=1)
segments(oversd,0,oversd,0.85,lwd=4,lty=2)


par(fig=c(50,100,5,50)/100,mar=c(5,4,2,1),new=T)
plot(NULL,xlim=c(0,20),ylim=c(0,0.9),cex.main=1.9,cex.axis=1.7,cex.lab=1.9,ylab="Probability",
     mgp=c(2.4,0.7,0),xlab="sd of variance in height (cm)",main="Informed priors (N = 100)",yaxs="i")

for(i in 1:10){
  diff10<-density(rnorm(1000000,BIdata100[i,5],BIdata100[i,6]))
  polygon(diff10$x,diff10$y,col=alpha("gray65",0.5),border=alpha("gray65",0.5))
  lines(diff10$x,diff10$y,col="gray50")
}
segments(sdpop_var,0,sdpop_var,0.85,lwd=4,lty=1)
segments(oversd,0,oversd,0.85,lwd=4,lty=2)




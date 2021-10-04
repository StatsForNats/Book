#####################################################################
########        STATISTICAL MODELING FOR NATURALISTS       #########
########                    CHAPTER VI                    #########
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

### LOAD DATA ###
# Subsamples were created elsewhere and are available in the website

reg_dataFull<-read.csv("SampFull.csv")  # Sample of 625
reg_data100<-read.csv("Samp100.csv")    # Sample of 100
reg_data40<-read.csv("Samp40.csv")      # sample of 40

### VISUALIZE FULL SAMPLE DATA ###

par(fig=c(10,80,5,76)/100,mar=c(5.1,4.1,0,0))
plot(reg_dataFull$lf,reg_dataFull$lg,pch=19,ylab="Height (cm)",xlab="Number of tillers",
     col="gray50",cex=0.5,cex.lab=1.2,xlim=c(0,14))

par(fig=c(10,80,76,95)/100,mar=c(0,4.1,0.7,0),new=TRUE)
barplot(hist(reg_dataFull$lf,breaks=14,plot=FALSE)$density,
        horiz = F,yaxt='n',col="gray75",xlim=c(0,14))

par(fig=c(80,99,5,76)/100,mar=c(5.1,0,0,0.7),new=TRUE)
barplot(hist(reg_dataFull$lg,breaks=20,plot=FALSE)$density,
        space=0,horiz = T,xaxt='n',col="gray75")


### PREPARING DATA FOR SAMPLE OF 100 FOR ANALYISIS ###
# For exemplification purposes, we will show how both models are run using the data
# with a sample = 100 (we live it up to the reader to try the other two subsamples)

max <- 17  # establish limit number of tillers (this removes outlier plant)
reg_data100$x<- reg_data100$lf-mean(reg_data100$lf)  # mean-centered the predictor
reg_data100$x2 <- reg_data100$x^2  # square of predictor
reg_data100$y <- reg_data100$lg
reg_data100$x0 <- na.omit(reg_data100$x)
ord <- order(reg_data100$x0)
reg_data100 <-subset(reg_data100,lf<max)


### MODEL BUILDING ###
# using full sample data (n=625)

# Gathering the data for Stan
stan_data <- list( y = reg_data100$y, x=reg_data100$x, x2=reg_data100$x2,
                   N = length(reg_data100$y)) 

## MODEL 6.1 - with diffuse priors ##

## Define parameters monitored
params <- c("a","b1","b2","sigma")

## MCMC settings
ni <- 6000  # total iterations per chain
nt <- 2     # interval of data collection
nb <- 2000  # burned iterations
nc <- 3     # number of mcmc chains

## Initial values
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,100),b1= rnorm(1,0,10),b2= rnorm(1,0,10),sigma=runif(1,1,10)) })

## Call Stan to run model_VI_1.stan
model1 <- stan("model_VI_1.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(model1, digits=3)

# Visual assessment of output
pairs(model1,pars=c("a","b1","b2","sigma"))
traceplot(model1)


## MODEL 6.2 - with informed priors ##

# Data, monitored parameters, MCMC settnigs and initial values all remain the same as before

## Call Stan to run model_VI_2.stan
model2 <- stan("model_VI_2.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(model2, digits=3)

## Visual assessment of output
pairs(model2,pars=c("a","b1","b2","sigma"))
traceplot(model2)


### YOU CAN RUN BOTH THE DIFFUSE AND INFORMED MODELS WITH THE OTHER SAMPLE SIZES
# (n=625 and n=40) BY REPEATING LINES 43 TO 101 BUT REPLACING reg_data100 WITH THE
# APPROPIATE DATA SET

#######################################

### PLOTING MODELS 6.1 AND 6.2 FOR N = 100 ###

## CALCULATING MODEL PREDICTIONS ##

tils <- seq(-6, 6, 0.1) # Creating a series of number of tiller values to predict responses
tils2<-tils^2 # Squared number of tillers for prediction

# Model 1:
post1 <- extract(model1) #Gets all values of parameters (6000 iterations)

Mu1<- rep(post1$a,length(tils))+outer(post1$b1,tils,FUN = '*')+ 
  outer(post1$b2,tils2,FUN = '*')#Calculates predicted values for tillers tils under each iteration

mu.mean1 <- apply(Mu1,2,mean)  #Calculate mean prediction
mu.PI1U <- mu.mean1+(qt(0.975,nrow(Mu1)-1)*apply(Mu1,2,sd)) # Calculate upper CI values
mu.PI1L <- mu.mean1+(qt(0.025,nrow(Mu1)-1)*apply(Mu1,2,sd)) # Calculate lower CI values

# Model 2:
post2 <- extract(model2) 

Mu2<- rep(post2$a,length(tils))+outer(post2$b1,tils,FUN = '*')+ 
  outer(post2$b2,tils2,FUN = '*')

mu.mean2 <- apply(Mu2,2,mean)
mu.PI2U <- mu.mean2+(qt(0.975,nrow(Mu2)-1)*apply(Mu2,2,sd))
mu.PI2L <- mu.mean2+(qt(0.025,nrow(Mu2)-1)*apply(Mu2,2,sd))


## PLOTS ##

## Model predictions

par(mfrow=c(1,2),mar=c(4,4,2,1))

plot(reg_data100$x,reg_data100$y,col="gray50",main="Diffuse Priors (n=100)",
     ylab="height (cm)",xlab="scaled number of tillers")
lines(tils,mu.mean1,col=rgb(0.4,0,0,1),lwd=2,lty=1)
polygon(c(tils,rev(tils)),c(mu.PI1U,rev(mu.PI1L)),border=NA,col=rgb(0.8,0,0,0.3))

plot(reg_data100$x,reg_data100$y,col="gray50",main="Informed Priors (n=100)",
     ylab="height (cm)",xlab="scaled number of tillers")
lines(tils,mu.mean2,col=rgb(0,0,0.4,1),lwd=2,lty=1)
polygon(c(tils,rev(tils)),c(mu.PI2U,rev(mu.PI2L)),border=NA,col=rgb(0,0,0.8,0.3))


## Prior vs posterior distributions for each parameter under both models

par(mfcol=c(1,3),mar=c(3,4,3,1))

# parameter 'a'
s_p <- 10000
prior_m_no <-rnorm(s_p, 0,100)
prior_m_in <-rnorm(s_p,29.32,11.59)
post_m_no <- post1$a
post_m_in <- post2$a

plot(density(prior_m_no), ylim=c(0,0.5),xlim=c(5,50),col="red",
     ylab="density",main ="Parameter A", xlab="",lty=1,lwd=3)
lines(density(prior_m_in),col="blue",lwd=3,lty=1)
lines(density(post_m_no),col="red",lwd=3,lty=3)
lines(density(post_m_in),col="blue",lwd=3,lty=3)

legend(6,0.48,c("Model 6.1 Prior","Model 6.2 Prior","Model 6.1 Posterior",
                   "Model 6.2 Posterior"), bty="n",lty=c(1,1,2,2),cex=1,
       col=c("red","blue","red","blue"),x.intersp=0.6,y.intersp=0.6,seg.len=1)

# parameter 'b1'
s_p <- 3000
prior_m_no <-rnorm(s_p, 0,sqrt(100))
prior_m_in <-rnorm(s_p,3.37,1.85)
post_m_no <- post1$b1
post_m_in <- post2$b1

plot(density(prior_m_no), ylim=c(0,1.5),xlim=c(-1,7),col="red",
     ylab="density",main ="Parameter B1", xlab="",lty=1,lwd=3)
lines(density(prior_m_in),col="blue",lwd=3,lty=1)
lines(density(post_m_no),col="red",lwd=3,lty=3)
lines(density(post_m_in),col="blue",lwd=3,lty=3)

# parameter 'b2'
s_p <- 3000
prior_m_no <-rnorm(s_p, 0,sqrt(10))
prior_m_in <-rnorm(s_p,-0.16,0.13)
post_m_no <- post1$b2
post_m_in <- post2$b2

plot(density(prior_m_no), ylim=c(0,6),xlim=c(-1,0.4),col="red",
     ylab="density",main ="Parameter B2", xlab="",lty=1,lwd=3)
lines(density(prior_m_in),col="blue",lwd=3,lty=1)
lines(density(post_m_no),col="red",lwd=3,lty=3)
lines(density(post_m_in),col="blue",lwd=3,lty=3)

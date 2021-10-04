#####################################################################
########        STATISTICAL MODELING FOR NATURALISTS       #########
########                      CHAPTER X                   #########
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


### LOADING AND PREPARING DATA ###
# remember the data in this chapter is simulated - see ChapterX_SimulatingData.R to see the process

cg_seed <- read.csv("cg_seed.csv", header=T)


### DATA VISUALIZATION ###

par(mfrow=c(2,2),mar=c(4,4,2,2),mgp=c(2.4,0.7,0))
hist(cg_seed$seed[cg_seed$treat=="auto"], main="auto",breaks = 0:35,
     ylim=c(0,400),col="gray40",xlab="")
hist(cg_seed$seed[cg_seed$treat=="self"], main="self",breaks = 0:35,
     ylim=c(0,400),col="gray40",xlab="")
hist(cg_seed$seed[cg_seed$treat=="cross"], main="cross",breaks = 0:35,
     ylim=c(0,400),col="gray40",xlab="")
hist(cg_seed$seed[cg_seed$treat=="control"], main="control",breaks = 0:35,
     ylim=c(0,400),col="gray40",xlab="")

### MODEL BUILDING ###

# Create dummy variables
cg_seed$self <-0
cg_seed$self[cg_seed$treat=="self" ] <- 1
cg_seed$cross <-0
cg_seed$cross[cg_seed$treat=="cross" ] <- 1
cg_seed$control <-0
cg_seed$control[cg_seed$treat=="control" ] <- 1


# Gathering the data for Stan
stan_data <- list( seed = cg_seed$seed, self=cg_seed$self, cross=cg_seed$cross,
                   control=cg_seed$control, plant = cg_seed$plant, 
                   N_plant=length(unique(cg_seed$plant)),
                   N = length(cg_seed$self))

## MODEL 10_1 - Poisson model without zero inflation ##

# Setting the parameters monitored 
params <- c("a","b1","b2","b3","sigmap","log_lik")

# Defining the MCMC settings 
ni <- 6000 	# number of iterations
nt <- 2	# interval for data collection
nb <- 3000	# number of iterations discarded
nc <- 3	# number of chains

# Defining initial values
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,100),b1= rnorm(1,0,10),b2= rnorm(1,0,10),b3= rnorm(1,0,10),sigmap=runif(1,1,10)) })

# Calling for Stan from R. Model_10_1
Model_10_1 <- stan("model_X_1.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

# Summary of model and visual assessment
print(Model_10_1, digits=3)

pairs(Model_10_1,pars=c("a","b1","b2","b3","sigmap"))
traceplot(Model_10_1,pars=c("a","b1","b2","b3","sigmap"))

## MODEL 10_2 - Zero-inflated Poisson model ##

# Setting the parameters monitored
params <- c("a","b1","b2","b3","p1","sigmap1",
            "q","z1","z2","z3","p2","sigmap2")

# Defining initial values
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,1),b1= rnorm(1,0,10),b2= rnorm(1,0,10),b3= rnorm(1,0,10),sigmap1=runif(1,1,10),
        q= rnorm(1,0,1),z1= rnorm(1,0,10),z2= rnorm(1,0,10),z3= rnorm(1,0,10),sigmap2=runif(1,1,10)) })

# MCMC settings remain unchanged

# Calling for Stan from R. Model_10_2
Model_10_2 <- stan("model_X_2.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

# Summary of model and visual assessment
print(Model_10_2, digits=3)

pairs(Model_10_2,pars=c("a","b1","b2","b3","sigmap1"))
traceplot(Model_10_2,pars=c("a","b1","b2","b3","sigmap1"))

pairs(Model_10_2,pars=c("q","z1","z2","z3","sigmap2"))
traceplot(Model_10_2,pars=c("q","z1","z2","z3","sigmap2"))

###################

### PLOT PREDICTED POSTERIORS ###

# Extract posteriors for each model

post1<-extract(Model_10_1)
post2<-extract(Model_10_2)

# Load table of observed means

obs_means<-read.csv("ObservedMeans.csv",header=T)

## FIGURE 10.2 - Predicted posterior for average seed counts (and observed 
# values underneath) for each of the 2 models

par(mfrow=c(1,2),mar=c(4,4,2,1),mgp=c(2.4,0.7,0))

plot(density(exp(post1$a)),xlim=c(0,17),ylim=c(-1,5),xlab ="Number of seeds",
     main="Poisson",col="darkgreen",lwd=2)
abline(h=-0.05,col="gray")
lines(density(exp(post1$a+post1$b1)),col="red",lwd=2)
lines(density(exp(post1$a+post1$b2)), col="blue",lwd=2)
lines(density(exp(post1$a+post1$b3)),col="black",lwd=2,lty=2)

segments(obs_means[2,1],-0.1,obs_means[2,1],-2,lwd=2,col="darkgreen")
segments(obs_means[2,2],-0.1,obs_means[2,2],-2,lwd=2,col="red")
segments(obs_means[2,3],-0.1,obs_means[2,3],-2,lwd=2,col="blue")
segments(obs_means[2,4],-0.1,obs_means[2,4],-2,lwd=2,lty=2,col="black")

legend(4,4.5,c("auto", "self", "cross", "control"),lty=c(1,1,1,2),lwd=2,bty="n",
       col=c("darkgreen","red","blue","black"),x.intersp=0.6,y.intersp = 0.6)

par(mar=c(4,3,2,2))
plot(density(exp(post2$q)),xlim=c(0,17),ylim=c(-0.3,1.5),xlab ="Number of seeds",
     main="Zero Inflated",col="darkgreen",lwd=2,ylab="")
abline(h=-0.015,col="gray")
lines(density(exp(post2$q+post2$z1)),col="red",lwd=2)
lines(density(exp(post2$q+post2$z2)), col="blue",lwd=2)
lines(density(exp(post2$q+post2$z3)),lwd=2,col="black",lty=2)

segments(obs_means[2,1],-0.03,obs_means[2,1],-2,lwd=2,col="darkgreen")
segments(obs_means[2,2],-0.03,obs_means[2,2],-2,lwd=2,col="red")
segments(obs_means[2,3],-0.03,obs_means[2,3],-2,lwd=2,col="blue")
segments(obs_means[2,4],-0.03,obs_means[2,4],-2,lwd=2,col="black",lty=2)



#### FIGURE 9.3 - Predicted posterior for probability of fruit set for Model_10_2 - zero-inflated model

par(mfrow=c(1,1))

plot(density(1-(1/(1 + (1/exp(post2$a))))),xlim=c(0,1),ylim=c(-10,50),
     xlab ="Probability of fruit mooring", main="",col="darkgreen",lwd=2,ylab="")
abline(h=-0.5)
lines( density(1-(1/(1 + (1/exp(post2$a+post2$b1))))),col="red",lwd=2)
lines( density(1-(1/(1 + (1/exp(post2$a+post2$b2))))), col="blue",lwd=2)
lines( density(1-(1/(1 + (1/exp(post2$a+post2$b3))))),lwd=2,col="black",lty=2)

segments(obs_means[1,1],-1,obs_means[1,1],-20,lwd=2,col="darkgreen")
segments(obs_means[1,2],-1,obs_means[1,2],-20,lwd=2,col="red")
segments(obs_means[1,3],-1,obs_means[1,3],-20,lwd=2,col="blue")
segments(obs_means[1,4],-1,obs_means[1,4],-20,lwd=2,col="black",lty=2)

legend(0.7,50,c("auto", "self", "cross", "control"),lty=c(1,1,1,2),lwd=2,bty="n",
       col=c("darkgreen","red","blue","black"),x.intersp=0.6,y.intersp = 0.6)


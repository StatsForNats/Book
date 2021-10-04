#####################################################################
########        STATISTICAL MODELING FOR NATURALISTS       #########
########                      CHAPTER V                   #########
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


## read the Hypericum data from file
orig_data <- read.table("hypericum_data_94_07.txt", header=T)
## change variables to include only reproductive individuals

dt <- subset(orig_data, !is.na(ht_init) & rp_init > 0 & year<1997 )
yr <- unique(dt$year)
dt$height <- dt$ht_init
dt$rep_structures <- dt$rp_init
dt$id <- dt$tag
sdfruits <- sd(dt$rep_structures)


##### MODEL BUILDING  #####

## MODEL 1  ##

stan_data <- list(height = dt$height, rep_structures=dt$rep_structures, 
                  N = length(dt$height))

## Define parameters monitored
params <- c("a","b","sigma")

## MCMC settings
ni <- 5000  # total iterations per chain
nt <- 2     # interval of data collection
nb <- 1000  # burned iterations - first 1000 are discarded
nc <- 3     # number of mcmc chains

## Initial values
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,500),b= rnorm(1,0,100),sigma=runif(1,1,400)) })

## Call Stan to run model_V_1.stan
model1 <- stan("model_V_1.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(model1, digits=3)

# Visual assessment of output
pairs(model1,pars=c("a","b","sigma"))
traceplot(model1)


## MODEL 2 ##

# Monitored parameters and MCMC settings remain the same as above, so there is no need to repeat them

# Both height and number of reproductive structures will be log-transformed for this model

dt$lnrep <- log(dt$rep_structures)
dt$lnht <-  log(dt$height)
sdlnfruits <- sd(dt$lnrep)

stan_data <- list(lnht = dt$lnht, lnrep=dt$lnrep , 
                  N = length(dt$lnht))
## Initial values
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,1),b= rnorm(1,0,1),sigma=runif(1,1,10)) })

## Call Stan to run model_V_2.stan
model2 <- stan("model_V_2.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(model2, digits=3)

# Visual assessment of output
pairs(model2)
traceplot(model2)

## MODEL 3 ##

stan_data <- list(height = dt$height, rep_structures=dt$rep_structures, 
                  N = length(dt$height))

params <- c("a","b","sigma")

ni <- 20000  # total iterations per chain
nt <- 2      # interval of data collection
nb <- 10000  # burned iterations
nc <- 3      # number of mcmc chains

## Initial values
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,10),b= rnorm(1,0,10),sigma=runif(1,1,10)) })

## Call Stan to run model_V_3.stan
model3 <- stan("model_V_3.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(model3, digits=4)

# Visual assessment of output
pairs(model3)
traceplot(model3)

## MODEL 4 ##

# Monitored parameters, MCMC settings and initial values remain the same as above, so there is no need to repeat them

## Call Stan to run model_V_3.stan
model4 <- stan("model_V_4.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(model4, digits=4)

# Visual assessment of output
pairs(model4)
traceplot(model4)


## MODEL 5 ##

## Define parameters monitored

params <- c("a","b","scale")

inits <- lapply(1:nc,function(i) {
  list (a= runif(1,0,10),b= runif(1,0,10),scale=runif(1,0.1,10)) })

model5 <- stan("model_V_5.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(model5, digits=3)

# Visual assessment of output
pairs(model5)
traceplot(model5)

#################################################

### PLOTTING AND COMPARING MODELS ###

## CALCULATING MODEL PREDICTIONS ##

ht <- seq(0.1, max(dt$ht_init), 0.1) # Creating a series of height values to predict responses
lnht<-log(ht) # Log of height

# Model 1:
post1 <- extract(model1) #Gets all values of parameters (6000 iterations)

Mu1<- rep(post1$a,length(ht))+outer(post1$b,ht,FUN = '*') #Calculates predicted values for heights ht under each iteration

mu.mean1 <- apply(Mu1,2,mean)  #Calculate mean prediction
mu.PI1U <- mu.mean1+(qt(0.975,nrow(Mu1)-1)*apply(Mu1,2,sd)) # Calculate upper CI values
mu.PI1L <- mu.mean1+(qt(0.025,nrow(Mu1)-1)*apply(Mu1,2,sd)) # Calculate lower CI values

# Model 2:
post2 <- extract(model2)

Mu2<- exp(rep(post2$a,length(lnht))+outer(post2$b,lnht, FUN = '*'))

mu.mean2 <- apply(Mu2,2,mean)
mu.PI2U <-mu.mean2+(qt(0.975,nrow(Mu2)-1)*apply(Mu2,2,sd))
mu.PI2L <- mu.mean2+(qt(0.025,nrow(Mu2)-1)*apply(Mu2,2,sd))

# Model 3:
post3 <- extract(model3)

Mu3<- rep(post3$a,length(ht))*t(outer(ht,post3$b,FUN = '^'))

mu.mean3 <- apply(Mu3,2,mean)
mu.PI3U <-mu.mean3+(qt(0.975,nrow(Mu3)-1)*apply(Mu3,2,sd))
mu.PI3L <- mu.mean3+(qt(0.025,nrow(Mu3)-1)*apply(Mu3,2,sd))

# Model 4:
post4 <- extract(model4)

Mu4<- exp(rep(post4$a,length(ht))*t(outer(ht,post4$b,FUN = '^'))) 

mu.mean4 <- apply(Mu4,2,mean)
mu.PI4U <-mu.mean4+(qt(0.975,nrow(Mu4)-1)*apply(Mu4,2,sd))
mu.PI4L <- mu.mean4+(qt(0.025,nrow(Mu4)-1)*apply(Mu4,2,sd))

# Model 5:
post5 <- extract(model5)

Mu5<- exp(rep(post5$a,length(ht))+outer(post5$b,ht,FUN = '*'))

mu.mean5 <- apply(Mu5,2,mean)
mu.PI5U <-mu.mean5+(qt(0.975,nrow(Mu5)-1)*apply(Mu5,2,sd))
mu.PI5L <- mu.mean5+(qt(0.025,nrow(Mu5)-1)*apply(Mu5,2,sd))

## PLOT ##

par(fig=c(0,100,5,95)/100,mar=c(4,5,2,1))
plot(dt$height,dt$rep_structures, main="Comparing models fit",xlab="Height (cm)",
     ylab="Number of fruits",pch=1,cex=1.2, xlim=c(0,80), ylim= c(0,2000),
     mgp=c(2.8,1,0),cex.main=1.5,cex.lab=1.2,cex.axis=1.2,col="gray60")

lines(ht,mu.mean1,col=rgb(0.4,0,0,1),lwd=2,lty=1)
polygon(c(ht,rev(ht)),c(mu.PI1U,rev(mu.PI1L)),border=NA,col=rgb(0.8,0,0,0.3))

lines(ht,mu.mean2,col=rgb(0,0,0.4,1),lwd=2,lty=2)
polygon(c(ht,rev(ht)),c(mu.PI2U,rev(mu.PI2L)),border=NA,col=rgb(0,0,0.7,0.3))

lines(ht,mu.mean3,col=rgb(0,0.4,0,1),lwd=2,lty=1)
polygon(c(ht,rev(ht)),c(mu.PI3U,rev(mu.PI3L)),border=NA,col=rgb(0,0.8,0,0.3))

lines(ht,mu.mean4,col=rgb(0,0.2,0,1),lwd=2,lty=2)
polygon(c(ht,rev(ht)),c(mu.PI4U,rev(mu.PI4L)),border=NA,col=rgb(0,0.7,0,0.3))

lines(ht,mu.mean5,col=rgb(0.1,0.1,0.1,1),lwd=2,lty=4)
polygon(c(ht,rev(ht)),c(mu.PI5U,rev(mu.PI5L)),border=NA,col=rgb(0.1,0.1,0.1,0.3))

legend(0,1900, c("Model 5.1","Model 5.2","Model 5.3","Model 5.4","Model 5.5"),
       lty=c(1,2,1,2,4), lwd=2, col=c(rgb(0.4,0,0,1),rgb(0,0,0.4,1),rgb(0,0.4,0,1),
       rgb(0,0.2,0,1),rgb(0.1,0.1,0.1,1)), bty="n",cex=1.4,x.intersp=0.6,
       y.intersp=0.5,seg.len=1)

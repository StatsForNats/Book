#####################################################################
########        STATISTICAL MODELING FOR NATURALISTS       #########
########                      CHAPTER IX                 #########
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
cd <- read.csv("seeds.csv", header=T)
cd <- subset(cd, seeds >0)

names(cd)
pop <- unique(cd$site)
yr <- unique(cd$year)
pop <- pop[order(pop)]
m_pop <- rep(0, length(pop))
m_yr <- rep(0, length(yr))

### VISUAL ASSESSMENT ###

m1 <- glm.nb(cd$seeds~1)
cf1 <- coef(m1)
m2 <- glm(cd$seeds~1,family = "poisson")
cf2 <- coef(m2)

b = seq(-1,60,1)

par(mfrow=c(1,2))
H1 <- hist(cd$seeds,breaks =b, prob=T,xlab="Seeds",ylab="density",
       main="Poisson",ylim=c(0,0.15),xlim=c(0,25),col="gray70",border="gray70")
H2 <- hist(rpois(200000,exp(cf2[1])),breaks =b, plot=FALSE)
lines(H2$mids,H2$density, type="h",lwd=5)
text(12,0.14,paste("mean= ",round(exp(cf2[1]),2)))

H1 <- hist(cd$seeds,breaks =b, prob=T,xlab="Seeds",ylab="density",main="Negative binomial",
           ylim=c(0,0.15),xlim=c(0,25),col="gray70",border="gray70")
H3 <- hist(rnegbin(200000,exp(cf1[1]),m1$theta),breaks =b, plot=FALSE)
lines(H3$mids,H3$density, type="h",lwd=4)

text(12,0.14,paste("mean= ",round(exp(cf1[1]),2)))
text(12,0.13,paste("k= ", round(m1$theta,2)))


### MODEL BUILDING ###

# First we must set population and year to be used as random effects

cd$population <- 1
## models with different families
for(i in 1:nrow(cd)){ 
  cd$population[i] <- which(pop==cd$site[i])
}

cd$yrs <- 1
## models with different families
for(i in 1:nrow(cd)){ 
  cd$yrs[i] <- which(yr==cd$year[i])
}


## MODEL 9.1 - Poisson with population and year as random effects ##

# Gathering the data for Stan

stan_data <- list(seeds=cd$seeds, N_population= length(pop),
                  N_yrs=length(yr), N=nrow(cd),population= cd$population,
                  yrs = cd$yrs)

# Setting the parameters monitored 
params <- c("b0","p","y","sigmap","sigmay","xsqp","lambda")

# Defining the MCMC settings 
ni <- 6000 	# number of iterations
nt <- 2	# interval for data collection
nb <- 2000	# number of iterations discarded
nc <- 3	# number of chains

# Defining initial values
inits <- lapply(1:nc, function(i) {
  list(b0 = rnorm(1,0,5), sigmap = runif(1,0,5), sigmay = runif(1,0,5))})

# Calling for Stan from R. Model_9_1 (Poisson.stan)
Model_9_1  <- stan("Poisson.stan",
                   data = stan_data, init = inits, pars = params,
                   chains = nc, iter = ni, warmup = nb, thin = nt,
                   seed = 1, open_progress = FALSE)

# 
print(Model_9_1 ,digits=3,depth=2)
pairs(Model_9_1,pars=c("b0","sigmap","sigmay"))
traceplot(Model_9_1,pars=c("b0","sigmap","sigmay"))
traceplot(Model_9_1,pars=c("p[1]","p[2]","p[3]","p[4]","p[5]","p[6]"))
traceplot(Model_9_1,pars=c("p[7]","p[8]","p[9]","p[10]","p[11]","p[12]"))
traceplot(Model_9_1,pars=c("y[1]","y[2]","y[3]","y[4]"))

## CACULATING DISPERSAL STATISTIC - POISSON

post1 <- extract(Model_9_1) # Extract posteriors

post1$xsqp # parameter xsqp allows us to calculate dispersion
N<-length(cd$seeds)
dispersal_p<-rep(NA,6000)
dispersal_p<-post1$xsqp/(N-3)
mean(dispersal_p) # Value is NOT close to 1, indicating the assumtpion that mean = variance is not met



## MODEL 9.2 - Negaive Binomial with population and year as random effects ##

# Data, MCMC settings, and initial values remain the same as before

# Define monitored parameters
params <- c("b0","p","y","sigmap","sigmay","scale","xsq","pbar")

# Calling for Stan from R. Model_9_2 (NB.stan)
Model_9_2  <- stan("NB.stan",
                   data = stan_data, init = inits, pars = params,
                   chains = nc, iter = ni, warmup = nb, thin = nt,
                   seed = 1, open_progress = FALSE)

print(Model_9_2 ,digits=3,depth=2,pars=c("b0","sigmap","sigmay"))

traceplot(Model_9_2,pars=c("b0","sigmap","sigmay"))
traceplot(Model_9_2,pars=c("p[1]","p[2]","p[3]","p[4]","p[5]","p[6]"))
traceplot(Model_9_2,pars=c("p[7]","p[8]","p[9]","p[10]","p[11]","p[12]"))
traceplot(Model_9_2,pars=c("y[1]","y[2]","y[3]","y[4]"))

## CACULATING DISPERSAL STATISTIC - NEGATIVE BINOMIAL

post2<-extract(Model_9_2)

post2$xsq
N <-length(cd$seeds)
dispersal_b<-rep(NA,6000)
dispersal_b<-post2$xsq/(N-4)
mean(dispersal_b)


### PLOTS ###

# Fixed effects - Poisson Vs Negative binomial

par(mfrow=c(1,2),mgp=c(1.7,0.8,0),mar=c(4,4,2,2))

plot(density(post1$b0),col="gray60",xlab="",ylab="Fixed effect",
     main="Poisson",xlim=c(2,3.2))
polygon(density(post1$b0)$x,density(post1$b0)$y,col="gray50",border = F)
abline(v=mean(post1$b0),lwd=3,lty=2,col="blue")

plot(density(post2$b0),col="gray60",xlab="",ylab="",
     main="Negative binomial",xlim=c(2,3.2))
polygon(density(post2$b0)$x,density(post2$b0)$y,col="gray50",border = F)
abline(v=mean(post2$b0),lwd=3,lty=2,col="blue")


# Random effects - Poisson Vs Negative binomial

par(mfrow=c(1,2),mgp=c(3,0.7,0),mar=c(5,5,2,2))

boxplot(cbind(post1$y,rep(100,nrow(post1$p)),post1$p),horizontal=T,
        main="Poisson", ylim=c(2-mean(post1$b0),3.2-mean(post1$b0)),
        xlab="",ylab="Random effects",outline=FALSE,col="gray50",axes=F)
box()
abline(v=0,lwd=3,lty=2,col="blue")
axis(2, labels = c("Y[1]","Y[2]","Y[3]","Y[4]","","P[1]","P[2]","P[3]",
                   "P[4]","P[5]","P[6]","P[7]","P[8]","P[9]","P[10]",
                   "P[11]","P[12]"), at=1:17,las = 1,tick=t)
axis(1,labels=c("-0.6 ","-0.4 ","-0.2 ","x = 2.64","+0.2 ","+0.4 ","+0.6 "),
     at=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6),las=2)

boxplot(cbind(post2$y,rep(100,nrow(post2$p)),post2$p),horizontal=T,
        main="Negative binomial", ylim=c(2-mean(post1$b0),3.2-mean(post1$b0)),
        xlab="",ylab="Random effects",outline=FALSE,col="gray50",axes=F)
box()
abline(v=0,lwd=3,lty=2,col="blue")
axis(2, labels = c("Y[1]","Y[2]","Y[3]","Y[4]","","P[1]","P[2]","P[3]",
                   "P[4]","P[5]","P[6]","P[7]","P[8]","P[9]","P[10]",
                   "P[11]","P[12]"), at=1:17,las = 1,tick=t)
axis(1,labels=c("-0.6 ","-0.4 ","-0.2 ","x = 2.64","+0.2 ","+0.4 ","+0.6 "),
     at=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6),las=2)


## Posterior distribution of the variance of the deviations 
# due to population and year.

par( mfrow=c(1,2),mar=c(4,4,2,2),mgp=c(2,0.7,0))

plot(density(post1$sigmap),lty=2, main="Poisson",xlab="Standard deviation",
    xlim=c(0,0.6),ylim=c(0,16),lwd=3,col="red")
lines(density(post1$sigmay),lwd=3,col="blue")

plot(density(post2$sigmap),lty=2, main="Negative binomial",xlab="Standard deviation",ylab="",
     xlim=c(0,0.6),ylim=c(0,16),lwd=3,col="red")
lines(density(post2$sigmay),lwd=3,col="blue")

legend(0.13,15,c("Population (sigmap)","Year (sigmay)"),lty=c(2,1),
       lwd=2,bty="n",col=c("red","blue"))


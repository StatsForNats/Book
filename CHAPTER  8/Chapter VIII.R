#####################################################################
########        STATISTICAL MODELING FOR NATURALISTS       #########
########                    CHAPTER VIII                  #########
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

orig_data <- read.table("hypericum_data_94_07.txt", header=T)
dt <- subset(orig_data, !is.na(ht_init) & rp_init > 0 & year<1997) ## change variables to include only reproductive individuals
yr <- unique(dt$year)
dt$lgh <- log(dt$ht_init)
dt$lfr <- log(dt$rp_init)
site <- unique(dt$bald)


### MODEL BUILDING AND PLOTTING ###


## MODEL 8.1 - No mixed effects - pooled data  ##

# Gathering the data for Stan
stan_data <- list( lfr = dt$lfr, lgh=dt$lgh,
                   N = length(dt$lfr))

# Setting the parameters monitored 
params <- c("a","b","sigma","dev","log_lik")

# Defining the MCMC settings 
ni <- 6000 	# number of iterations
nt <- 2	# interval for data collection
nb <- 2000	# number of iterations discarded
nc <- 3	# number of chains

# Defining initial values
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,100),b= rnorm(1,0,10),sigma=runif(1,1,10)) })

# Calling for Stan from R. Model8.1
model_8_1 <- stan("model_VIII_1.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

# Summarizing and plotting the posteriors 
print(model_8_1, digits = 4)

pairs(model_8_1,pars=c("a","b","sigma"))
traceplot(model_8_1,pars=c("a","b","sigma"))

# Extracting posteriors and calculating mean and CI for prediction
lnhgt <- seq(1, log(max(dt$ht_init)), 0.01) # a predictor vector for plotting predictions

post1 <- extract(model_8_1)
Mu1<- rep(post1$a,length(lnhgt))+outer(post1$b,lnhgt, FUN = '*')

mu.mean1 <- apply(Mu1,2,mean)
mu.PI1U <-mu.mean1+(qt(0.975,nrow(Mu1)-1)*apply(Mu1,2,sd))
mu.PI1L <- mu.mean1+(qt(0.025,nrow(Mu1)-1)*apply(Mu1,2,sd))

# plot

par(fig=c(0,100,5,95)/100,mar=c(4,5,2,1))
plot(dt$lgh,dt$lfr, main="",xlab="log(height(cm))",
     ylab="log (number of fruits)",pch=16,cex=0.8, xlim=c(2,4), 
     ylim= c(0,8), col="gray60")

lines(lnhgt,mu.mean1,lwd=1.5)
polygon(c(lnhgt,rev(lnhgt)), c(mu.PI1U,rev(mu.PI1L)),
        border=NA,col=rgb(0,0,0.4,0.3))


###############
## MODEL(S) 8.2 - no mixed effects, but one model per population ##

# Here well do only a few populations (1, 57, and 91) for exemplification purposes

sel_pops <- c(1,57,91) # Selects populations for modeling and plotting
sel_popsi<-as.numeric(match(sel_pops,site)) # finds the position in vector "site"

pops <- array(0,c(4,length(sel_pops))) # creates an array to save the mean of each parameter for each population
colnames(pops) <- sel_pops
rownames(pops) <- c("Beta1","Beta2","Sigma","N")

ylabel<-c("log(number of fruits)","","")

# Monitored parameters, MCMC settings and initial values will remain the same
# in all models (and the same as Model8.1). Only the data would change.

lnhgt <- seq(1, log(max(dt$ht_init)), 0.01) # a predictor vector for plotting predictions

par(mfrow=c (1,3),mar=c(5,5,2,2))
for (i in 1:3){
  dtp <- subset(dt,bald==site[sel_popsi[i]]) #subset data
  stan_data <- list( lfr = dtp$lfr, lgh=dtp$lgh,
                     N = length(dtp$lfr))  
  model_8_2 <- stan("model_VIII_1.stan", 
             data = stan_data, init =inits, pars = params, chains = nc,
             iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)
  
  post2 <- extract(model_8_2)
  
  #Save parameter values:
  pops[1,i]<-round(mean(post2$a),2)
  pops[2,i]<-round(mean(post2$b),2)
  pops[3,i]<-round(mean(post2$sigma),2)
  pops[4,i]<-length(dtp$lfr)
  
  #calculating mean and CI for prediction

  Mu2<- rep(post2$a,length(lnhgt))+outer(post2$b,lnhgt, FUN = '*')
  mu.mean2 <- apply(Mu2,2,mean)
  mu.PI2U <-mu.mean2+(qt(0.975,nrow(Mu2)-1)*apply(Mu2,2,sd))
  mu.PI2L <- mu.mean2+(qt(0.025,nrow(Mu2)-1)*apply(Mu2,2,sd))

  # plot
  
  plot(dtp$lgh,dtp$lfr, main="",xlab="log(height(cm))",
       ylab="log (number of fruits)",pch=16,cex=0.8, xlim=c(2,4), 
       ylim= c(0,8), col="gray60",cex.lab=1.5,cex.axis=1.4)
  text(2.5,8,paste("population"," = ", sel_pops[i]),cex=1.5,pos=4)
  lines(lnhgt,mu.mean2)
  polygon(c(lnhgt,rev(lnhgt)), c(mu.PI2U,rev(mu.PI2L)),
          border=NA,col=rgb(0,0,0.4,0.3))
  
}

pops


###############
## MODEL 8.3 - Mixed effects, random intercepts but a common slope  ##

# First assign an id to each population (from 1 to 14)
pjlev<-1:length(site)
dt$pj<-as.numeric(pjlev[match(dt$bald,site)])

# Gathering the data for Stan
stan_data <- list( lfr = dt$lfr, lgh=dt$lgh, pj = dt$pj,
                   N = length(dt$lfr), N_pj = length(unique(dt$pj)))

# Defining initial values
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,100),b= rnorm(1,0,10),sigma=runif(1,1,10)) })

# Setting the parameters monitored
params <- c("a","b","p","sigma","dev","log_lik")

#Note: MCMC settings remain unchanged

# Calling for Stan from R. Model8.3
model_8_3 <- stan("model_VIII_3.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

# Visual assessment of model results
pairs(model_8_3,pars=c("a","b","sigma"))
traceplot(model_8_3,pars=c("a","b","sigma"))

## Extracting posteriors 
post3<-extract(model_8_3)

# Calculating mean for overall prediction (fixed effects only) - CI could be calculated as before
lnhgt <- seq(1, log(max(dt$ht_init)), 0.01) # a predictor vector for plotting predictions

Mu3<- rep(post3$a,length(lnhgt))+outer(post3$b,lnhgt, FUN = '*')
mu.mean3 <- apply(Mu3,2,mean)

# Calculating mean and CI for predictions of each population (fixed + random effects),
# To compare with previous model, we will do populations 1, 57 and 91

sel_pops <- c(1,57,91)
sel_popsi<-as.numeric(match(sel_pops,site))

par(mfrow=c (1,3),mar=c(5,5,2,2))
for (i in 1:3){
  dtp <- subset(dt,bald==site[sel_popsi[i]])
  #calculating mean and CI for each population prediction
  Mu3P<- rep(post3$a+post3$p[,sel_popsi[i]],length(lnhgt))+
          outer(post3$b,lnhgt, FUN = '*')
  mu.mean3P <- apply(Mu3P,2,mean)
  mu.PI3PU <-mu.mean3P+(qt(0.975,nrow(Mu3P)-1)*apply(Mu3P,2,sd))
  mu.PI3PL <- mu.mean3P+(qt(0.025,nrow(Mu3P)-1)*apply(Mu3P,2,sd))
  
  # plot
  
  plot(dtp$lgh,dtp$lfr, main="",xlab="log(height(cm))",
       ylab="log (number of fruits)",pch=16,cex=0.8, xlim=c(2.5,4), 
       ylim= c(0,7), col="gray60",cex.lab=1.5,cex.axis=1.4)
  text(2.5,7,paste("population"," = ", sel_pops[i]),cex=1.5,pos=4)
  lines(lnhgt,mu.mean3,lty=2,lwd=2)
  
  polygon(c(lnhgt,rev(lnhgt)), c(mu.PI3PU,rev(mu.PI3PL)),
          border=NA,col=rgb(0,0,0.5,0.4))
  lines(lnhgt,mu.mean3P,col=rgb(0,0,0.5,1),lwd=2)
  
}

################
## MODEL 8.4 - Mixed effects, random intercepts and slopes  ##

# Assign an id to each population (from 1 to 14) -if not done in previous step
pjlev<-1:length(site)
dt$pj<-as.numeric(pjlev[match(dt$bald,site)])

# Gathering the data for Stan
stan_data <- list( lfr = dt$lfr, lgh=dt$lgh, pj = dt$pj,
                   N = length(dt$lfr), N_pj = length(unique(dt$pj)))

# Defining initial values
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,100),b= rnorm(1,0,10),sigma=runif(1,1,10)) })

# Setting the parameters monitored
params <- c("a","b","p","bp","sigma","dev","log_lik")

# Calling for Stan from R. Model8.3
model_8_4 <- stan("model_VIII_4.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)


# Visual assessment of model results
pairs(model_8_4,pars=c("a","b","sigma"))
traceplot(model_8_4,pars=c("a","b","sigma"))

## Extracting posteriors 
post4<-extract(model_8_4)

# Calculating mean for overall prediction (fixed effects only) - CI could be calculated as before
lnhgt <- seq(1, log(max(dt$ht_init)), 0.01) # a predictor vector for plotting predictions

Mu4<- rep(post4$a,length(lnhgt))+outer(post4$b,lnhgt, FUN = '*')
mu.mean4 <- apply(Mu4,2,mean)

# Calculating mean and CI for predictions of each population (fixed + random effects),
# To compare with previous models, we will do populations 1, 57 and 91

sel_pops <- c(1,57,91)
sel_popsi<-as.numeric(match(sel_pops,site))

par(mfrow=c (1,3),mar=c(5,5,2,2))
for (i in 1:3){
  dtp <- subset(dt,bald==site[sel_popsi[i]])
  #calculating mean and CI for each population prediction
  Mu4P<- rep(post4$a+post4$p[,sel_popsi[i]],length(lnhgt))+
    outer(post4$b+post4$bp[,sel_popsi[i]],lnhgt, FUN = '*')
  mu.mean4P <- apply(Mu4P,2,mean)
  mu.PI4PU <-mu.mean4P+(qt(0.975,nrow(Mu4P)-1)*apply(Mu4P,2,sd))
  mu.PI4PL <- mu.mean4P+(qt(0.025,nrow(Mu4P)-1)*apply(Mu4P,2,sd))
  
  # plot
  
  plot(dtp$lgh,dtp$lfr, main="",xlab="log(height(cm))",
       ylab="log (number of fruits)",pch=16,cex=0.8, xlim=c(2,4), 
       ylim= c(0,8), col="gray60",cex.lab=1.5,cex.axis=1.4)
  text(2.5,8,paste("population"," = ", sel_pops[i]),cex=1.5,pos=4)
  lines(lnhgt,mu.mean4,lty=2,lwd=2)
  
  polygon(c(lnhgt,rev(lnhgt)), c(mu.PI4PU,rev(mu.PI4PL)),
          border=NA,col=rgb(0,0,0.5,0.4))
  lines(lnhgt,mu.mean4P,col=rgb(0,0,0.5,1),lwd=2)
  
}

### MODEL COMPARISON ###

#In order to do the model comparison, firts you need to run (load) both 
# the waic_sn and WAICtab functions from the 'Functions.R' script located
# in the Extra folder.

WAICtab(model_8_1,model_8_3,model_8_4)





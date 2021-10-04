#####################################################################
########        STATISTICAL MODELING FOR NATURALISTS       #########
########                    CHAPTER XI                    #########
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
orig_data <- read.csv("Hypericum_cumulicola.csv", header=T)

## subset for plants with the reproductive structures and complete data
dt <- subset(orig_data, !is.na(init_height) & !is.na(stems) & surv_fin!=2)

yr <- unique(dt$time) #studied years
site <- unique(dt$site) # studied sites

dt$lgh <- log(dt$init_height) #log of height
dt$lgh2 <- dt$lgh^2
dt$lghc <- dt$lgh- mean(dt$lgh) # centered for interpretation

dt$stems[dt$stems>8] <- 8

## Time-since-fire 
dt$TSF <- dt$time -dt$burn_yr
tsf <- unique(dt$TSF)
tsf <- tsf[order(tsf)]

### DATA VISUALIZATION ###

par(mfrow=c(1,2))
# Survival by height
b<-seq(0,95,5)
hist(dt$init_height[dt$surv_fin==1],breaks=b,xlab="Height (cm)",
     main="Surviving plants",ylim=c(0,2000))
hist(dt$init_height[dt$surv_fin==0],breaks=b,xlab="Height (cm)",
     main="Non surviving plants",ylim=c(0,2000))

# Survival by number of stems
b<-seq(0,8,1)
hist(dt$stems[dt$surv_fin==1],breaks=b,xlab="Number of stems",
     main="Surviving plants",ylim=c(0,4000))
hist(dt$stems[dt$surv_fin==0],breaks=b,xlab="Number of stems",
     main="Non surviving plants",ylim=c(0,4000))

# Survival by height
b<-seq(0,50,5)
hist(dt$TSF[dt$surv_fin==1],breaks=b,xlab="TSF (years)",
     main="Surviving plants",ylim=c(0,3000))
hist(dt$TSF[dt$surv_fin==0],breaks=b,xlab="TSF (years)",
     main="Non surviving plants",ylim=c(0,3000))


### MODEL BUILDING ###

# NOTE: These models can take a long time to run (some taking up to over an hour)
# If you prefer to load models provided in the website (extension .RData),
# please go to line 292

## Model 11.1 ##

# Gathering the data for Stan
stan_data <- list( surv_fin = dt$surv_fin , lgh=dt$lgh , TSF=dt$TSF ,
                   stems=dt$stems, lgh2=dt$lgh2,
                   N = length(dt$surv_fin))

# Defining the MCMC settings 
ni <- 6000
nt <- 2
nb <- 3000
nc <- 3

# Setting the parameters monitored 
params <- c("a","b","c","d","bc","bd","cd","dev","log_lik")

# Defining initial values
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,100),b= rnorm(1,0,1),c= rnorm(1,0,1),d= rnorm(1,0,1),
        bd= rnorm(1,0,1),bc= rnorm(1,0,1),cd= rnorm(1,0,1)) })

# Calling for Stan from R. Model_11_1
Model_11_1 <- stan("model_XI_1.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

# Summary of model and visual assessment
print(Model_11_1, digits=3)

traceplot(Model_11_1,pars=c("a","b","c"))
traceplot(Model_11_1,pars=c("d","bc","bd"))

#####      #####      ######      ######

## Model 11.2 ##

stan_data <- list( surv_fin = dt$surv_fin , lgh=dt$lgh , TSF=dt$TSF , stems=dt$stems,
                   N = length(dt$surv_fin))

params <- c("a","b","c","d","dev","log_lik")

inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,100),b= rnorm(1,0,1),c= rnorm(1,0,1)) })

Model_11_2 <- stan("model_XI_2.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(Model_11_2, digits=3)

traceplot(Model_11_2,pars=c("a","b","c","d"))

#####      #####      ######      ######

## Model 11.3 ##

stan_data <- list( surv_fin = dt$surv_fin , lgh=dt$lgh , TSF=dt$TSF ,
                   stems=dt$stems, lgh2=dt$lgh2,
                   N = length(dt$surv_fin))

params <- c("a","b","bb","c","d","bc","bd","cd","bbc","bbd","dev","log_lik")

inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,10),b= rnorm(1,0,1),c= rnorm(1,0,1),cd=rnorm(1,0,1),d= rnorm(1,0,1),
        bd= rnorm(1,0,1),bbd= rnorm(1,0,1),
        bc= rnorm(1,0,1),bb= rnorm(1,0,1),bbc= rnorm(1,0,1)) })

Model_11_3 <- stan("model_XI_3.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(Model_11_3, digits=3,pars=c("a","b","c","d","bc","bd","cd","bb","bbc","bbd"))

traceplot(Model_11_3,pars=c("a","b","c","d","bb"))
traceplot(Model_11_3,pars=c("bc","bd","cd","bbc","bbd"))

#####      #####      ######      ######

## Model 11.4 ## 

stan_data <- list( surv_fin = dt$surv_fin , lgh=dt$lgh , TSF=dt$TSF ,
                   stems=dt$stems, lgh2=dt$lgh2,
                   N = length(dt$surv_fin))

params <- c("a","b","bb","c","d","dev","log_lik")
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,100),b= rnorm(1,0,1),bb= rnorm(1,0,1),c= rnorm(1,0,1),d= rnorm(1,0,1)) })

Model_11_4 <- stan("model_XI_4.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(Model_11_4, digits=3,pars=c("a","b","c","d","bb"))

traceplot(Model_11_4,pars=c("a","b","c","d","bb"))

#####      #####      ######      ######

## Model 11.5 ## 

stan_data <- list( surv_fin = dt$surv_fin , lgh=dt$lgh , TSF=dt$TSF ,
                   stems=dt$stems, lgh2=dt$lgh2,
                   N = length(dt$surv_fin))

params <- c("a","b","bb","c","dev","log_lik")
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,100),b= rnorm(1,0,1),bb= rnorm(1,0,1),c= rnorm(1,0,1)) })

Model_11_5 <- stan("model_XI_5.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(Model_11_5, digits=3,pars=c("a","b","c","bb"))

traceplot(Model_11_5,pars=c("a","b","c","bb"))

#####      #####      ######      ######

## Model 11.6  Stan ##

stan_data <- list( surv_fin = dt$surv_fin , lgh=dt$lgh , TSF=dt$TSF ,
                   stems=dt$stems, lgh2=dt$lgh2,
                   N = length(dt$surv_fin))

params <- c("a","b","bb","d","dev","log_lik")
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,100),b= rnorm(1,0,1),bb= rnorm(1,0,1),d= rnorm(1,0,1)) })

Model_11_6 <- stan("model_XI_6.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(Model_11_6, digits=3,pars=c("a","b","d","bb"))

traceplot(Model_11_6,pars=c("a","b","d","bb"))

#####      #####      ######      ######

## Model 11.7 ## 

stan_data <- list( surv_fin = dt$surv_fin , lgh=dt$lgh , TSF=dt$TSF ,
                   stems=dt$stems, lgh2=dt$lgh2,
                   N = length(dt$surv_fin))

params <- c("a","b","bb","c","bc","bbc","dev","log_lik")
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,100),b= rnorm(1,0,1),c= rnorm(1,0,1),
        bc= rnorm(1,0,1),bb= rnorm(1,0,1),bbc= rnorm(1,0,1)) })

Model_11_7 <- stan("model_XI_7.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(Model_11_7, digits=3,pars=c("a","b","c","bb","bc","bbc"))

traceplot(Model_11_7,pars=c("a","b","c"))
traceplot(Model_11_7,pars=c("bb","bc","bbc"))

## Model 11.8 ## 


stan_data <- list( surv_fin = dt$surv_fin , lgh=dt$lgh , TSF=dt$TSF ,
                   stems=dt$stems, lgh2=dt$lgh2,
                   N = length(dt$surv_fin))

params <- c("a","b","bb","d","bd","bbd","dev","log_lik")
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,100),b= rnorm(1,0,1),bb= rnorm(1,0,1),d= rnorm(1,0,1),
        bd= rnorm(1,0,1),bbd= rnorm(1,0,1)) })

Model_11_8 <- stan("model_XI_8.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(Model_11_8, digits=3,pars=c("a","b","d","bb","bd","bbd"))

traceplot(Model_11_8,pars=c("a","b","d"))
traceplot(Model_11_8,pars=c("bb","bd","bbd"))

## Model 11.9 ##

stan_data <- list( surv_fin = dt$surv_fin , lgh=dt$lgh , TSF=dt$TSF ,
                   stems=dt$stems, lgh2=dt$lgh2,
                   N = length(dt$surv_fin))

params <- c("a","b","bb","dev","log_lik")
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,100),b= rnorm(1,0,1),
        bb= rnorm(1,0,1)) })

Model_11_9 <- stan("model_XI_9.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(Model_11_9, digits=3,pars=c("a","b","bb"))

traceplot(Model_11_9,pars=c("a","b","bb"))

#####      #####      ######      ######

## Model 11.10 ##


stan_data <- list( surv_fin = dt$surv_fin , lgh=dt$lgh , TSF=dt$TSF ,
                   stems=dt$stems, lgh2=dt$lgh2,
                   N = length(dt$surv_fin))

params <- c("a","b","dev","log_lik")
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,100),b= rnorm(1,0,1)) })

Model_11_10 <- stan("model_XI_10.stan", 
                data = stan_data, init =inits, pars = params, chains = nc,
                iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(Model_11_10, digits=3,pars=c("a","b"))

traceplot(Model_11_10,pars=c("a","b"))


### MODEL SELECTION ###

## Comparing WAIC scores

#In order to perfrom waic comparison, first you need to run (load) both 
# the waic_sn and WAICtab functions from the 'Functions.R' script located
# in the Extra folder.

WAICtab(Model_11_1,Model_11_2,Model_11_3,Model_11_4,
        Model_11_5,Model_11_6,Model_11_7,Model_11_8,
        Model_11_9,Model_11_10)


###########

### PLOTS ###

## Extract posteriors from models winning model and from 2 simplest models
#(linear - model 11.10, and quadratic - model 11.9) ##

post3<-extract(Model_11_3)
post9<-extract(Model_11_9)
post10<-extract(Model_11_10)

## FIGURE 11.1 : Linear vs quadratic models (simplest models) ##

# Calculate means and CI #
x_pred<-seq(1,max(dt$init_height),0.1)
lnx_pred<-log(x_pred)  # log of x
lnx2_pred <- lnx_pred^2

# Model 10:
Mu10<- (exp(rep(post10$a,length(lnx_pred))+outer(post10$b,lnx_pred,FUN = '*')))/
  (1+exp(rep(post10$a,length(lnx_pred))+outer(post10$b,lnx_pred, FUN = '*'))) # We calculate P(survival) using the link function

mu.mean10 <- apply(Mu10,2,mean)  #Calculate mean prediction
mu.PI10U <- mu.mean10+(qt(0.975,nrow(Mu10)-1)*apply(Mu10,2,sd)) # Calculate upper CI values
mu.PI10L <- mu.mean10+(qt(0.025,nrow(Mu10)-1)*apply(Mu10,2,sd))

# Model 9:
Mu9<- (exp(rep(post9$a,length(lnx_pred))+outer(post9$b,lnx_pred,FUN = '*')+
              outer(post9$bb,lnx2_pred,FUN = '*')))/
  (1+exp(rep(post9$a,length(lnx_pred))+outer(post9$b,lnx_pred, FUN = '*')+
           outer(post9$bb,lnx2_pred,FUN = '*'))) # We calculate P(survival) using the link function

mu.mean9 <- apply(Mu9,2,mean)  #Calculate mean prediction
mu.PI9U <- mu.mean9+(qt(0.975,nrow(Mu9)-1)*apply(Mu9,2,sd)) # Calculate upper CI values
mu.PI9L <- mu.mean9+(qt(0.025,nrow(Mu9)-1)*apply(Mu9,2,sd))

## Plot model predictios vs observed data

par(mfrow=c(1,1),mar=c(4,4,2,2),mgp=c(2.4,0.7,0))
plot(exp(lnx_pred),seq(0,1,length=length(lnx_pred)),type="n",ylab= "P(survival)",
     xlab="Height (cm)",log="x",ylim=c(0.4,0.9))

# Observed values
b <- c(-2,seq(0,max(dt$lgh),0.3)) # values for binning survival data
z <- cut(dt$lgh, b) # binned survival data
prebyden <- tapply(dt$surv_fin,z,sum)
tab <- table(z)
probs <-prebyden/tab
probs <- as.vector(probs)

for (k in 1:(length(prebyden)-1)){
  points(exp(b[k]),probs[k],pch=1,cex=log(prebyden[k]),
         col="black")}

# Predicted values
lines(x_pred,mu.mean10,col="red")
polygon(c(x_pred,rev(x_pred)),c(mu.PI10U,rev(mu.PI10L)),
        border=NA,col=rgb(0.8,0,0,0.3))

lines(x_pred,mu.mean9,col="blue")
polygon(c(x_pred,rev(x_pred)),c(mu.PI9U,rev(mu.PI9L)),
        border=NA,col=rgb(0,0,0.8,0.3))

legend(23,0.9,c("Model 9", "Model 10"),col=c("blue","red"),bty="n",cex=0.8,lty=1)
#####

## FIGURE 11.2 : Model 11.3 coefficients

post3b<-as.data.frame(do.call(cbind,post3))

par(mfrow=c(1,1),mar=c(4,6,2,2),mgp=c(1,0.5,0))

boxplot(post3b[,c(7,9,8,6,5,3,10,4,2,1)],horizontal=T,lty=1,outline=FALSE,axes=F)
box()
axis(2, labels = c("TSF:stems","height2:stems","height2:TSF","height:stems","height:TSF",
                   "height2","stems","TSF","height","intercept"), at=1:10,las = 1,tick=T)
axis(1,labels=c("-0.2","0.0","0.2","0.4","0.6","0.8","1.0","1.2"),
     at=c(-0.2,0,0.2,0.4,0.6,0.8,1,1.2),las=2)

#########
## FIGURE 11.3 : Some Model 11.3 predictions (TST = 2 , Stems =1,2,3)

# Values for model prediction:
x_pred<-seq(1,max(dt$init_height),0.1) # X
lnx_pred<-log(x_pred)  # log of x
lnx2_pred <- lnx_pred^2 # Square of log of x
tsf<-2 # Change this if you want to test a different TSF
st<-c(1,2,4)

# Visual elements
color_lines <- c("red","blue","black")
dottypes<-c(0,1,2)
MLab<-paste("TST =",tsf)
b <- c(-2,seq(0,max(dt$lgh),0.3)) # values for binning survival data

par(mfrow=c(1,1),mar=c(4,4,2,2),mgp=c(2.4,0.7,0)) 
plot(exp(lnx_pred),seq(0,1,length=length(lnx_pred)),type="n",ylim=c(0,1),
ylab="P(survival",xlab="Height (cm)",main=MLab,log="x")

# Loop to add observed data
for (i in 1:length(st)){
  z <- cut(dt$lgh[dt$stems==st[i] & dt$TSF>(tsf-3) & dt$TSF<(tsf+3)],b)
  prebyden <- tapply(dt$surv_fin[dt$stems==st[i] & dt$TSF>(tsf-3) & dt$TSF<(tsf+3)],z,sum)
  tab <- table(z)
  probs <-prebyden/tab
  probs <- as.vector(probs)

  for (k in 1:(length(prebyden)-1)){
    points(exp(b[k]),probs[k],pch=dottypes[i],cex=log(prebyden[k])*0.8,
           col=color_lines[i])}
}

# Loop to add predictions:
for(i in 1:length(st)){
  # First we should calculate the values of each part
  parA<-rep(post3$a,length(lnx_pred))
  parB<-outer(post3$b,lnx_pred,FUN = '*')
  parBB<- outer(post3$bb,lnx2_pred,FUN = '*')
  parC<-rep(post3$c*tsf,length(lnx_pred))
  parBC<-outer(post3$bc,lnx_pred,FUN = '*')*tsf
  parBBC<-outer(post3$bbc,lnx2_pred,FUN = '*')*tsf
  parD<-rep(post3$d*st[i],length(lnx_pred))
  parBD<-outer(post3$bd,lnx_pred,FUN = '*')*st[i]
  parBBD<-outer(post3$bbd,lnx2_pred,FUN = '*')*st[i]
  parCD<-rep(post3$cd*tsf*st[i],length(lnx_pred))
  
  Mu<- (exp(parA+parB+parBB+parC+parBC+parBBC+parD+parBD+
        parBBD+parCD))/(1+exp(parA+parB+parBB+parC+parBC+
        parBBC+parD+parBD+parBBD+parCD))
  
  mu.mean <- apply(Mu,2,mean)

  lines(x_pred,mu.mean,col=color_lines[i],lwd=3)
  }

legend(1,0.4,c("1 stem","2 stems","4 stems"),lty=1,lwd=2,bty="n",pch=c(0,1,2),
col=c("red","blue","black"),x.intersp=0.6,y.intersp = 0.6,seg.len=0.8)

  




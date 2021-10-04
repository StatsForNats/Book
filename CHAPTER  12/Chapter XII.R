#####################################################################
########        STATISTICAL MODELING FOR NATURALISTS       #########
########                     CHAPTER XII                  #########
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

orig_data <- read.table("hypericum_data_94_07.txt", header=T)
dist_data <- read.table("Bald_mat_dis.txt", header=T)
coord_data <- read.table("Bald_coordinates.txt", header=T)
coord_data <- as.data.frame(coord_data)

# Seting site information
site <- unique(orig_data$bald)
ran_sites <-sample(site,14)
ran_sites <- ran_sites[order(ran_sites)]

# Creating distance matrix
dis_pop_ran <- as.matrix(dist_data/1000,14,14)
colnames(dis_pop_ran) <- ran_sites
rownames(dis_pop_ran) <- ran_sites

# Subseting data and setting variables
dt <- subset(orig_data, !is.na(ht_init) & !is.na(st_init) & rp_init > 0 & 
               year==1995 & bald %in% ran_sites )
yr <- unique(dt$year)

dt$lgh <- log(dt$ht_init)
dt$lfr <- log(dt$rp_init)
dt$tsf <- dt$year-dt$fire_year

dt$lgh_s <- scale(dt$lgh)
dt$tsf_s <- scale(dt$tsf)

###############

### MODEL BUILDING ###

## Model_12_1 : Model with no random effects ##

# Gathering data for Stan
stan_data <- list( lfr = log(dt$rp_init), lgh_s=dt$lgh_s[,1], tsf_s=dt$tsf_s[,1],
                   N = length(log(dt$rp_init)))

# Setting the parameters monitored 
params <- c("a","b","c", "cc","sigma","dev","log_lik")

# Defining the MCMC settings 
ni <- 10000 # number of iterations
nt <- 2 # interval for data collection
nb <- 3000 # number of iterations discarded
nc <- 3 # number of chains

# Defining initial values
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,1),b= rnorm(1,0,1),c= rnorm(1,0,1),cc= rnorm(1,0,1),sigma=runif(1,1,5)) })

# Calling for Stan from R. Model_12_1
Model_12_1 <- stan("model_XII_1.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

# Summary of model and visual assessment
print(Model_12_1, digits=3,pars=c("a","b","c", "cc","sigma"))

traceplot(Model_12_1,pars=c("a","b","c", "cc","sigma"))

########
## Model_12_2 : Random intercepts due to population ##

# Assign populations indices
dt$pj <- 1
for (i in 2: length(ran_sites)){
  dt$pj[dt$bald==ran_sites[i]] <-i}

stan_data <- list( lfr = log(dt$rp_init), lgh_s=dt$lgh_s[,1], tsf_s=dt$tsf_s[,1],pj=dt$pj,
                   N = length(log(dt$rp_init)),N_pj=length(unique(dt$pj)))

params <- c("a","b","c", "cc","sigma","sigmap","dev","log_lik")

inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,1),b= rnorm(1,0,1),c= rnorm(1,0,1),cc= rnorm(1,0,1),
        sigma=runif(1,1,5),sigmap=runif(1,1,5)  ) })

Model_12_2 <- stan("model_XII_2.stan", 
               data = stan_data, init =inits, pars = params, chains = nc,
               iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(Model_12_2, digits=3,pars=c("a","b","c", "cc","sigma","sigmap"))

traceplot(Model_12_2,pars=c("a","b","c", "cc","sigma","sigmap"))


########
## Model_12_3 : correlated random effects among populations
# and with interaction between time-since fire and height. ##

stan_data <- list( lfr = log(dt$rp_init), lgh_s=dt$lgh_s[,1], tsf_s=dt$tsf_s[,1],pj=dt$pj,
                   Dmat = dis_pop_ran,
                   N = length(log(dt$rp_init)),N_pj=length(unique(dt$pj)))

params <- c("a","b","c", "cc","p","etasq","rhosq","sigmap","sigma","dev","log_lik")

inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,1),b= rnorm(1,0,1),c= rnorm(1,0,1),cc= rnorm(1,0,1),
        sigma=runif(1,1,5),sigmap=runif(1,1,5),etasq=runif(1,1,5),rhosq=runif(1,1,5)    ) })

Model_12_3<- stan("model_XII_3.stan", 
              data = stan_data, init =inits, pars = params, chains = nc,
              iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(Model_12_3, digits=3,pars=c("a","b","c","cc","p","etasq",
                              "rhosq","sigmap","sigma"))

traceplot(Model_12_3,pars=c("a","b","c", "cc","etasq",
                            "rhosq","sigmap","sigma"))

#In order to perfrom waic comparison, first you need to run (load) both 
# the waic_sn and WAICtab functions from the 'Functions.R' script located
# in the Extra folder.

WAICtab(Model_12_1,Model_12_2,Model_12_3)

########
## Model_12_4 : correlated random effects among populations
# but without interaction between time-since fire and height.

params <- c("a","b","c","etasq","rhosq","sigmap","sigma","dev","log_lik")

Model_12_4<- stan("model_XII_4.stan", 
              data = stan_data, init =inits, pars = params, chains = nc,
              iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(Model_12_4, digits=3,pars=c("a","b","c","etasq",
                                  "rhosq","sigmap","sigma"))

traceplot(Model_12_4,pars=c("a","b","c","etasq","rhosq",
                            "sigmap","sigma"))


########
## Model_12_5 : correlated random effects among populations
# but without interaction between time-since fire and height.

params <- c("a","b","etasq","rhosq","sigmap","sigma","dev","log_lik")

Model_12_5<- stan("model_XII_5.stan", 
              data = stan_data, init =inits, pars = params, chains = nc,
              iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(Model_12_5, digits=3,pars=c("a","b","etasq","rhosq","sigmap","sigma"))

traceplot(Model_12_5,pars=c("a","b","etasq","rhosq","sigmap","sigma"))


### MODEL COMPARISON ###

## Comparing WAIC scores

WAICtab(Model_12_1,Model_12_2,Model_12_3,Model_12_4,
        Model_12_5)


### PLOTTING POSTERIORS ###

## Figure 12.3 : Posterior distribution of spatial covariance between populations
par(mfrow=c (1,1)) 
x <- seq(0.1,14,0.1)
post3 <- extract(Model_12_3)
plot (x,median(post3$etasq)*exp(-median(post3$rhosq)*x^2),xlab="Distance (km)", type="n",
      ylab="Covariance",lwd=3,ylim=c(0,2.5),xlim=c(0,14),bty="n")
for (i in 1:100){
  lines(x,post3$etasq[i]*exp(-post3$rhosq[i]*x^2),
        col="grey60",lwd=1)
}
lines(x,median(post3$etasq)*exp(-median(post3$rhosq)*x^2),col="red",lwd=2)


## Figure 12.5 - Posterior distribution of coefficients

# Extract posterior from models
post1<-extract(Model_12_1)
post2<-extract(Model_12_2)
post3<-extract(Model_12_3)

par(mfrow=c(2,2),mar=c(3,4,2,2))
plot(NA,xlim=c(3.8,4.6),ylim=c(0,16),main="Intercept, a",xlab="",ylab="density")
lines(density(post1$a),col="red")
lines(density(post2$a),col="blue")
lines(density(post3$a),col="black")

plot(NA,xlim=c(0.95,1.2),ylim=c(0,20),main="Height, b",xlab="",ylab="density")
lines(density(post1$b),col="red")
lines(density(post2$b),col="blue")
lines(density(post3$b),col="black")

plot(NA,xlim=c(-0.5,0.5),ylim=c(0,18),main="TSF, c",xlab="",ylab="density")
lines(density(post1$c),col="red")
lines(density(post2$c),col="blue")
lines(density(post3$c),col="black")

plot(NA,xlim=c(-0.2,0.1),ylim=c(0,18),main="TSF, c",xlab="",ylab="density")
lines(density(post1$cc),col="red")
lines(density(post2$cc),col="blue")
lines(density(post3$cc),col="black")

legend(0,17,c("Model 12.1","Model 12.2","Model 12.3"),lty=1,
       col=c("red","blue","black"),bty="n",cex=0.7,y.intersp = 0.4)

## Figure 12.6 - Predicted number of fruits per bald

post3<-extract(Model_12_3)
# Set visual parameters
coltsf <- c("red","springgreen","royalblue")

tsfc <- c(3,2,3,1,1,1,3,3,2,2,2,2,3,2)
nsize = c(25, 42)
size <- log(nsize)-mean(dt$lgh)

par(mfrow=c(1,2),mar=c(5,5,2,1))
for(i in 1:2){
  plot(coord_data[,2],coord_data[,3],type="n", ylab= "Coordinates (N)",xlab="Coordinates (E)",
       xlim=c(-0.15,1.2), ylim=c(-0.25,9), main = paste( nsize[i], "cm in Height"))
  
  for (k in 1:14){
    TSF <- unique(dt$tsf_s[dt$bald==ran_sites[k]])
    pred_S<-median(post3$a+(post3$b*size[i])+post3$p[,k]+
                     (post3$c*TSF)+(post3$cc*size*TSF))
    points(coord_data[k,2],coord_data[k,3],
           cex=exp(pred_S)/30,pch=16,col=coltsf[tsfc[k]])
    
    if (i == 1) {legend(0.6,8, legend=c("recently burned", "Intermediate","long-unburned"), 
                    pch= c(16,16,16), col=coltsf,bty="n",cex=0.8,
                    x.intersp = 0.3,y.intersp = 0.3)  
      
      points(0.9,5.5,cex=(100/30),pch=16)
      text(0.9,5,"100 Fruits")
    }
  }}


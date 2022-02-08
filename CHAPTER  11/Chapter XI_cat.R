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
dt$lgh2 <- dt$lgh^2  # square of log of height (for quadratic response)
dt$lghc <- dt$lgh- mean(dt$lgh) # centered for interpretation

dt$stems[dt$stems>8] <- 8

## Time-since-fire 
dt$TSF <- dt$time -dt$burn_yr
tsf <- unique(dt$TSF)
tsf <- tsf[order(tsf)]


####################################################################
## SETTING DUMMY VARIABLES ##
# The reference represents 6 or less years since last fire
# TSF2 represents intermediate burns (6 - 10 years) and 
# TSF3 long-unburned (over 20 years since last fire)

dt$TSF2 <- 0   
dt$TSF2[dt$TSF>6 & dt$TSF<=10] <- 1
dt$TSF3 <- 0
dt$TSF3[dt$TSF>20] <- 1

####################################################################

## MODEL 11.3 ##

#WARNING: This model can take several hours to run, even in powerful
# computers. For a simpler model that demonstrates the use of dummy
# variables, skip to the next model (Model 11.5)

# Defining the MCMC settings 
ni <- 6000
nt <- 2
nb <- 3000
nc <- 3

# Gathering the data for Stan
stan_data <- list( surv_fin = dt$surv_fin , lgh=dt$lgh , TSF2=dt$TSF2 ,
                   TSF3=dt$TSF3 ,stems=dt$stems, lgh2=dt$lgh2,
                   N = length(dt$surv_fin))

# Setting the parameters monitored 
params <- c("a","b","bb","c2","c3","d","bc2","bc3","bd","c2d","c3d",
            "bbc2","bbc3","bbd","dev","log_lik")

# Defining initial values
inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,10),b= rnorm(1,0,1),c2= rnorm(1,0,1),c3= rnorm(1,0,1),
        c2d=rnorm(1,0,1),c3d=rnorm(1,0,1),d= rnorm(1,0,1),
        bd= rnorm(1,0,1),bbd= rnorm(1,0,1),
        bc2= rnorm(1,0,1),bc3= rnorm(1,0,1),bb= rnorm(1,0,1),bbc2= rnorm(1,0,1),
        bbc3= rnorm(1,0,1)) })

# Calling for Stan from R. Model_11_3_cat
Model_11_3_cat <- stan("model_XI_3_cat.stan", 
                   data = stan_data, init =inits, pars = params, chains = nc,
                   iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

# Summary of model and visual assessment
print(Model_11_3_cat, digits=3,pars=c("a","b","c2","c3","d","bc2","bc3","bd","c2d",
                                  "c3d","bb","bbc2","bbc3","bbd"))

traceplot(Model_11_3_cat,pars=c("a","b","c2","c3"))
traceplot(Model_11_3_cat,pars=c("bd","bbd","d","bb"))
traceplot(Model_11_3_cat,pars=c("bc2","bc3","c2d","c3d"))
traceplot(Model_11_3_cat,pars=c("bbc2","bbc3"))

####################################################################

## MODEL 11.5 ##

ni <- 6000
nt <- 2
nb <- 3000
nc <- 3

stan_data <- list( surv_fin = dt$surv_fin , lgh=dt$lgh , TSF2=dt$TSF2 ,
                   TSF3=dt$TSF3 , lgh2=dt$lgh2, N = length(dt$surv_fin))

params <- c("a","b","bb","c2","c3","dev","log_lik")

inits <- lapply(1:nc,function(i) {
  list (a= rnorm(1,0,10),b= rnorm(1,0,1),
        c2= rnorm(1,0,1),c3= rnorm(1,0,1)) })

Model_11_5_cat <- stan("model_XI_5_cat.stan", 
                       data = stan_data, init =inits, pars = params, chains = nc,
                       iter = ni, warmup = nb, thin = nt, seed =1, open_progress = FALSE)

print(Model_11_3_cat, digits=3,pars=c("a","b","c2","c3","bb"))

traceplot(Model_11_3_cat,pars=c("a","b","c2","c3","bb"))

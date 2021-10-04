#####################################################################
########        STATISTICAL MODELING FOR NATURALISTS       #########
########           CHAPTER X - DATA SIMULATION            #########
########    QUINTANA-ASCENCIO, LOPEZ BORGHESI, MENGES    #########
#################################################################

##### SETTING REQUIREMENTS #####

rm(list=ls())  # this line clears the environment - helps prevent confusion

# Set working directory to source file 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### LOADING DATASETS AND SIMULATING DATA ###
# remember the data for this chapter is simulated from parameters found in 
# Evans et al. (2003). This is done to ensure the sources of zeroes are known.

# Read the original field data from file (to calculate observed parameters)
fruit_data <- read.csv("HCfruits.csv", header=T)
seed_data <- read.csv("HCseeds.csv", header=T)

# Prepare arrays for saving observed mean probability of fruit set and mean number
# of seeds per fruit
T_P <- array(0, c(2,4))
colnames(T_P) <- c("auto", "self", "cross", "control")
rownames(T_P) <- c("fruit_set","seeds")

# Generate the extra zeros
zero <- array(0,c(20,4)) # Array to store number of extra zeros per plant 
colnames(zero) <- c("auto", "self", "cross", "control")

s <- 0
for (j in 1:4){
  T_P[1,j] <- round(mean(fruit_data$nofull[fruit_data$treat==j+s]/ #calculates mean probability of fruit set
                           fruit_data$nostruct[fruit_data$treat==j+s]),3)
  T_P[2,j] <- round(mean(seed_data$NOFULL[seed_data$treat==j+s]/ #calculates mean number of seeds per fruit
                           seed_data$X_FREQ_[seed_data$treat==j+s]),1)
  zero[,j] <- 20-20*(rbinom(20,20,T_P[1,j])/20)  #calculates extra zeros per plant per treatment
  s <-1
}

# Generate seed counts and group them with extra zeroes per plant per treatment
yes_o  <- yes_c <- yes_s <- yes_a <- NULL
for(j in 1:20){
  yes_a <-  c(yes_a,rpois(20-zero[j,1],T_P[2,1]),rep(0,zero[j,1]))
  yes_s <-  c(yes_s,rpois(20-zero[j,2],T_P[2,2]),rep(0,zero[j,2]))
  yes_c <-  c(yes_c,rpois(20-zero[j,3],T_P[2,3]),rep(0,zero[j,3]))
  yes_o <-  c(yes_o,rpois(20-zero[j,4],T_P[2,4]),rep(0,zero[j,4]))
}

#Put everything together in simulated dataset
id <- rep(1:80, each = 20, len = 1600)
treat <- rep(c("auto","self","cross","control"),each=20*20,len=20*20*4)

cg_seed <- cbind(treat,id,c(yes_a,yes_s,yes_c,yes_o))
colnames(cg_seed) <- c("treat","plant","seed")
cg_seed <- as.data.frame(cg_seed)

# Save simulated data to csv

write.csv(cg_seed,"cg_seed.csv",row.names = F)
write.csv(T_P,"ObservedMeans.csv",row.names = F)


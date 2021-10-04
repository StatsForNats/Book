#####################################################################
########        STATISTICAL MODELING FOR NATURALISTS       #########
########          CHAPTER VI - BUILDING PRIORS            #########
########    QUINTANA-ASCENCIO, LOPEZ BORGHESI, MENGES    #########
#################################################################

##### SETTING REQUIREMENTS #####

rm(list=ls())  # this line clears the environment - helps prevent confusion

# Set working directory to source file 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### LOAD DATA

reg_data<-read.table("wg_priors.txt",header=T) # file with data from other studies used for priors
reg_data1 <- read.table("wg_data.txt", header=T) # our data
reg_data1<-reg_data1[!duplicated(reg_data1$tag),]


### VISUAL ASSESSMENT - COMPARING DATA FROM PRIOR STUDIES TO OUR OWN

studyname<-c("DHt122","DHt155","Dh187","DHt222","Our Study",
             "DHt264","DHt332","DHt412","KK") #Name of each study (for plot tags)

ylabels<-c("Height","","","Height","","","Height","","")
xlabels<-c("","","","","","","Number of tillers",
           "Number of tillers","Number of tillers")

co<-c(1,3,5,7,9,10,11,13,15) # List of columns in reg_data containing number of tiller information
# for each prior study

par(mfrow=c(3,3),mar=c(5,5,1,1))

for(i in 1:9) {
  
  if(i==5){ #This plots our study's data
    plot(reg_data1$lf, reg_data1$lg,ylim=c(0,70),xlim=c(0,25),
         pch=16, cex=1.3 ,main="",xlab=xlabels[i],ylab=ylabels[i])
    
  }else{  #This plots data of each study used for priors
    x1<-reg_data[,co[i]] # These lines retrieve number of tiller data for each study
    x1<-x1[!is.na(x1)]
    x<-x1[x1<=max(reg_data1$lf)]
    
    y1<-reg_data[,co[i]+1] # These lines retrieve height data for each study
    y1<-y1[!is.na(y1)]
    y<-y1[x1<=max(reg_data1$lf)]
    
    plot(x,y,main="",xlab=xlabels[i],ylab=ylabels[i],cex.lab=2,
         ylim=c(0,70),xlim=c(0,25),pch=16, cex=1.3,col="blue")
    
  }
  text(25,65,studyname[i],pos=2,cex=2)
}

### SET AN ARRAY TO SAVE RESULTS FOR EACH PRIOR DATA SET

names(reg_data) # Checking name of columns in compiled dataframe

PriorParams<-array(0,c(8,3)) #an array to save the model parameters calculated for each data set

colnames(PriorParams)<-c("b0","b1","b2")

names<-c("DHt122","DHt155","Dh187","DHt222","DHt264","DHt332","DHt412","KK")

rownames(PriorParams)<-names

co<-c(1,3,5,7,9,11,13,15) # List of columns containing number of tiller information
# for each prior study

for(i in 1:8) {
  
  x1<-reg_data[,co[i]] # Retrieves number of tiller data for each prior study
  x1<-x1[!is.na(x1)]
  x<-x1[x1<=max(reg_data1$lf)] # Keep only plants with tillers equal or less to the maximum number in our study
  
  
  y1<-reg_data[,co[i]+1] # Retrieves plant height data for each prior study
  y1<-y1[!is.na(y1)]
  y<-y1[x1<=max(reg_data1$lf)]  # Keep only plants with tillers equal or less to the maximum number in our study
  
  xs<- scale(x,scale=FALSE) # mean-center the number of tillers (just like in our models)
  x2<-x^2 # create the square of the predictor (as we suspect a quadratic relationship)
  
  model<-lm(y~xs+x2)  #use lm to calculate parameters of a quadratic model for each prior data set
  
  PriorParams[i,]<- as.numeric(model$coeff)	# Save said coefficients to our result array
  
}

### CALCULATING VALUES FOR INFORMED PRIOR!!!
PriorParamsMeans<-apply(PriorParams,2,mean) # Calculate prior mean for each parameter
PriorParamsSD<-apply(PriorParams,2,sd) # Calculate prior standard deviation for each parameter

##########################################################################################

#####             Statistical Modeling for Naturalists                #####

### CHAPTER 8: Random effects of populations: revisiting reproductive outputs ###


This folder contains the following material:

    - Chaper VIII.R:
            The code to run the examples in this chapter. This script uses "hypericum_data_94_07.txt" data and calls the 3 stan models: 'model_VIII_1.stan', 'model_VIII_3.stan' and 'model_VIII_4.stan'

    - model_VIII_1.stan
            Stan model to perform a model with no random effects - on pooled data (Model 8.1) and on each population (Models 8.2).
 
    - model_VIII_3.stan
            Stan model to perform a mixed effects model with random variation around the intercept by population.

    - model_VIII_4.stan
            Stan model to perform a mixed effects model with random variation around both the intercept and the slope by population.

    - hypericum_data_94_07.txt

All 3 models are power models (log-log) - linear model with normal likelihoods on log transformed data.



Instructions:

To run this example, save all files in the same folder, open the R script 'Chaper VIII.R' and run all lines.


NOTE: Before being able to run the model comparison in lines 284, you would need to open the file "Funcions.R" located in the 'Extra' folder, and run it without closing your RStudio session- This will create the 2 functions needed for the model comparison (waic_sn, and WAICtab).

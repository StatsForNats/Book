#####             Statistical Modeling for Naturalists                #####

### CHAPTER 11: Model selection for mixed models: Effects of size, reproduction, and time-since-fire on survival ###


WARNING: Some of the models presented in this chapter are computantionally intensive (some taking over an hour in relatively powerful computers).


This folder contains the following material:

    - Chaper XI.R:
            The code to run the examples in this chapter. This script uses "Hypericum_cumulicola.csv" and calls all the stan files (10 files total).

    - model_XI_1.stan
            Stan code incorporating height, number of stems, and time-since-fire variables , and 2-way interactions between each pair of variables.

    - model_XI_2.stan
            Stan code incorporating height, number of stems, and time-since-fire variables but no interactions

    - model_XI_3.stan
            Stan code incorporating height, number of stems, and time-since-fire variables , with a quadratic effect due to height, and 2-way interactions between each pair of variables.

    - model_XI_4.stan
            Stan code incorporating height, number of stems, and time-since-fire variables , with a quadratic effect due to height, but no interactions.

    - model_XI_5.stan
            Stan code incorporating height and time-since-fire variables, with a quadratic effect due to height, and no interactions.

    - model_XI_6.stan
            Stan code incorporating height and number of stems variables, with a quadratic effect due to height, and no interactions.

    - model_XI_7.stan
            Stan code incorporating height and time-since-fire variables, with a quadratic effect due to height, and interaction between both variables.

    - model_XI_8.stan
            Stan code incorporating height and number of stems variables, with a quadratic effect due to height, and interaction between both variables.

    - model_XI_9.stan
            Stan code incorporating just height with a quadratic effect.

    - model_XI_10.stan
            Stan code incorporating just height without a quadratic effect.


    - Hypericum.cumulicola.csv


All models are linear models with a binomal likelihood and a logit link function - i.e. the logarithm of the odds of survival is a linear function of different combinations of predictors.


Instructions:

Download "Chapter XI.R", "Hypericum_cumulicola.csv" and all 10 .stan files in a single chapter. Run all lines.


NOTE: Before being able to run the model comparison in lines 318-320, you would need to open the file "Funcions.R" located in the 'Extra' folder, and run it without closing your RStudio session- This will create the 2 functions needed for the model comparison (waic_sn, and WAICtab).

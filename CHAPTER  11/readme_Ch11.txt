#####             Statistical Modeling for Naturalists                #####

### CHAPTER 11: Model selection for mixed models: Effects of size, reproduction, and time-since-fire on survival ###


WARNING: Some of the models presented in this chapter are computantionally intensive (some taking several hours in relatively powerful computers).

This folder contains 2 sets of scripts, one used to run all competing models using time-since-fire (TSF) as a quantitative variable (as presented in the book) and another to show how the models could be run with TSF as a categorical variable using dummy variables. The folder also contains the data needed to complete this exercise (Hypericum.cumulicola.csv).


Files for models with TSF as quantitative variable:

    - Chaper XI.R:
            The code to run the examples in this chapter. This script uses "Hypericum_cumulicola.csv" and calls the correspinding stan files (10 files total).

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


Files for models with TSF as categorical variable:

    - Chaper XI_cat.R:
            The code to run the examples in this chapter. This script uses "Hypericum_cumulicola.csv" and calls .

    - model_XI_3_cat.stan
            Stan code for the winning model from the chapter (Model 11.3) but using TSF as categorical - 3 categories. (This model is computationally intensive, taking several hours to run)

    - model_XI_5_cat.stan
            Stan code for simplest model (Model 11.5) using TSF as categorical (This model is provided as simpler alternative to show how dummy variables work).


All models are linear models with a binomal likelihood and a logit link function - i.e. the logarithm of the odds of survival is a linear function of different combinations of predictors.



Instructions:

- For TSF as quantitative:
Download "Chapter XI.R", "Hypericum_cumulicola.csv" and all 10 .stan files in a single folder. Run all lines.

NOTE: Before being able to run the model comparison in lines 318-320, you would need to open the file "Funcions.R" located in the 'Extra' folder, and run it without closing your RStudio session- This will create the 2 functions needed for the model comparison (waic_sn, and WAICtab).

- For TSF as categorical variable:
Download "Chapter XI_cat.R", "Hypericum_cumulicola.csv" and all 2 .stan files in a single folder. Model 11.3 is computationally intensive, so skip if necessary.


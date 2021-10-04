#####             Statistical Modeling for Naturalists                #####

### CHAPTER 5: The importance of assumptions: Evaluating reproductive output ###


This folder contains the following material:

    - Chaper V.R:
            The code to run the examples in this chapter. This script uses 'hypericum_data_94_07.txt'  data and calls the 5 stan models: 'model_V_1.stan','model_V_2.stan','model_V_3.stan','model_V_4.stan', and 'model_V_5.stan'

    - model_V_1.stan
            Stan model to perform linear model with normal likelihood

    - model_V_2.stan
            Stan model to perform a power model (log-log) - linear model with normal likelihood on log transformed data.

    - model_V_3.stan
            Stan model to perform a non-linear model with normal likelihood (and diffuse priors)

    - model_V_4.stan
            Stan model to perform a non-linear model with log-normal likelihood (and diffuse priors)

    - model_V_5.stan
            Stan model to perform a linear model with negative binomial likelihood and a logarithmic link function.

    - hypericum_data_94_07.txt


Note: all 5 models were performed with diffuse priors


Instructions:

To run this example, save all files in the same folder, open the R script 'Chaper V.R' and run all lines



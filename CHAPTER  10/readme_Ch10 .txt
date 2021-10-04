#####             Statistical Modeling for Naturalists                #####

### CHAPTER 10: Count data: Seeds per fruit ###


This folder contains the following material:

    - Chaper X.R:
            The code to run the examples in this chapter. This script uses "cg_seed.csv" and "ObservedMeans.csv" data and calls the 2 stan models: 'model_X_1.stan', and 'model_X_2.stan'

    - model_X_1.stan
            Stan model to perform a Poisson model with no zero inflation
 
    - model_X_2.stan
            Stan model to perform a mixed zero-inflated linear model with Binomial likelihoods for the probability of producing fruit, and Poisson likelihood for the number of seeds produced.
    
    - cg_seed.csv
            Simulated data of seeds per fruit 

    - ObservedMeans.csv
            Observed means for each treatment

    - Chapter X_SimulatingData.R
            The code to generate simulated data using published data for parameterization. It uses "HCfruits.csv" and "HCseeds.csv"

    - HC_fruits.csv
            Data on fruit set from Evans et. al 2003

    - wg_priors.txt
            Data on seed production from Evans et. al 2003



Instructions:

To run this example, save files 'Chaper X.R', 'model_X_1.stan', 'model_X_2.stan', 'cg_seed.csv', and 'ObservedMeans.csv' in the same folder, open the R script 'Chaper VI.R' and run all lines.

To see how simulated data was generated, save files 'Chapter X_SimulatingData.R', 'HC_fruits.csv', and 'HCseeds.csv' in the same foldar, open the R script and run all lines.


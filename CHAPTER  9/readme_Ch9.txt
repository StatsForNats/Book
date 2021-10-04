#####             Statistical Modeling for Naturalists                #####

### CHAPTER 9: Count data: Seeds per fruit ###


This folder contains the following material:

    - Chaper IX.R:
            The code to run the examples in this chapter. This script uses "seeds.csv" data and calls the 2 stan models: 'Poisson.stan',and 'NBstan'

    - Poisson.stan
            Stan model to perform a model with Poisson likelihood and logarithmic link function
 
    - NB.stan
            Stan model to perform a model with Negative Binomial likelihood and logarithmic link function

    
    - seeds.csv

Both models have random effects by population and by year.


Instructions:

To run this example, save all files in the same folder, open the R script 'Chaper IX.R' and run all lines.



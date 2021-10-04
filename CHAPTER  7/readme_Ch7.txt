#####             Statistical Modeling for Naturalists                #####

### CHAPTER 7: Analyzing binary responses: Trade-offs between reproduction and survival ###


This folder contains the following material:

    - Chaper VII.R:
            The code to run the examples in this chapter. This script uses "Hypericum_cumulicola.csv" data and calls the 2 stan models: 'Chapter_VII_natural_scale.stan', and 'Chapter_VII_power_model.stan'

    - Chapter_VII_natural_scale.stan
            Stan model to perform a linear model with binomial likelihood and a logit link function - the logarithm of the odds of survival changes as a linear function of the number of fruits

    - Chapter_VII_power_model.stan
            Stan model to perform a power model with binomial likelihood and logit link function - the logarithm of the odds of survival changes as power function of the number of fruits by log-transforming number of fruits.

    - Hypericum_cumulicola.csv



Instructions:

To run this example, save all files in the same folder, open the R script 'Chaper VII.R' and run all lines.





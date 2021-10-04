#####             Statistical Modeling for Naturalists                #####

### CHAPTER 6: The use of priors: Plant height and number of tillers in wire grass ###


This folder contains the following material:

    - Chaper VI.R:
            The code to run the examples in this chapter. This script uses 'Samp40.csv','Samp100.csv' and 'SampFull.csv'  data and calls the 2 stan models: 'model_VI_1.stan', and 'model_VI_2.stan'

    - model_VI_1.stan
            Stan model to perform linear model with normal likelihood

    - model_VI_2.stan
            Stan model to perform a power model (log-log) - linear model with normal likelihood on log transformed data.

    - Samp40.csv
            Contains a subsample of data with n=40

    - Samp100.csv
            Contains a subsample of data with n=100

    - SampFull.csv
            Contains a full sample of data with n=625

    - ChapterVI_CauculatingPriors.R
            The code to calculate prior parameter distributions from previous studies found in the literature. It uses 'wg_data.txt' and 'wg_priors.txt'

    - wg_data.txt
            Data from our study on wiregrass

    - wg_priors.txt
            Data from previous studies on wiregrass



Instructions:

To run this example, save files 'Chaper VI.R', 'model_VI_1.stan', 'model_VI_2.stan', 'Samp40.csv', 'Samp100.csv', 'SampFull.csv' in the same folder, open the R script 'Chaper VI.R' and run all lines.

To see how priors were calculated, save files 'ChapterVI_CauculatingPriors.R', 'wg_data.txt', and 'wg_priors.txt'. open the R script and run all lines.




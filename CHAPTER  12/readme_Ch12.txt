#####             Statistical Modeling for Naturalists                #####

### CHAPTER 12: Spatial and temporal variation: Effects of fire and habitat structure on reproduction ###



This folder contains the following material:

    - Chaper XII.R:
            The code to run the examples in this chapter. This script uses "hypericum_data_94_07.txt" as well as "Bald_Coordinates.txt" and "Bald_mat_dis.txt". It calls 'model_XII_1.stan', 'model_XII_2.stan', 'model_XII_3.stan', 'model_XII_4.stan', and 'model_XII_5.stan' 

    - model_XII_1.stan
            Model with no random effects where the normal likelihood of the log of fruits is a function of the log of height, time since fire and the interaction between both.

    - model_XII_2.stan
            Model with random intercetps due to populations where the normal likelihood of the log of fruits is a function of the log of height, time since fire and the interaction between both.

    - model_XII_3.stan
            Model with correlated random effects among populations where the normal likelihood of the log of fruits is a function of the log of height, time since fire and the interaction between both.

    - model_XII_4.stan
            Model with correlated random effects among populations where the normal likelihood of the log of fruits is a function of the log of height, time since fire but no interactions.

    - model_XII_5.stan
            Model with correlated random effects among populations where the normal likelihood of the log of fruits is a function of the log of height (no effect of time-since-fire).

    - hypericum_data_94_07.txt
            Main data set

    - Bald_coordinates.txt
            Coordinates of all 14 populations

    - Bald_mat_dis.txt
            A matrix with pair-wise distances between all 14 populations


Instructions:

To run this example, download Chapter XII.R, model_XII_1.stan, model_XII_2.stan, model_XII_3.stan, model_XII_4.stan, model_XII_5.stan, hypericum_data_94_07.txt, Bald_coordinates.txt, and Bald_mat_dis.txt into the same folder. Open the R script and run every line.



NOTE: Before being able to run the model comparisons in line 137 and in line 175, you would need to open the file "Funcions.R" located in the 'Extra' folder, and run it without closing your RStudio session- This will create the 2 functions needed for the model comparison (waic_sn, and WAICtab).

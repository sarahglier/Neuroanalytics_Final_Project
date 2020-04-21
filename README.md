# Neuroanalytics_Final_Project
Simulated Datasets for UNC Neuroanalytics Course with accompanying script


These are simulated cortisol (cort_wide.csv) and participant demographic information (id_data.csv) datasets.  The cortisol wide contains cortisol concentrations across time and the id_data contains sex and age information.  All these datasets were created in R using a randomize function and setting the min and max boundaries.


You'll find the accompanying R script titled Final_Project in which I have coded an analysis for this simulated dataset. Specifically, it contains code to call and load in .csv files, extensive data munging and tidying, examining the distributions of raw data, multiple imputation, piecewise growth curve modeling with landmark registration, and comparisons against standard multilevel models and repeated measures ANOVA.  

The script is designed for anyone to follow along and to try to make growth curve modeling more accessible to the larger research community, even for someone with basic R scripting knowledge.  I also tried to add in helpful notes if someone wants to use their own dataset and follow along with the script.  


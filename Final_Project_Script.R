#Author: Sarah Glier, MD/PhD Candidate
#Institution: UNC Chapel Hill
#Date: April 21, 2020
#Neuroanalytics Final Project

#https://github.com/sarahglier/Neuroanalytics_Final_Project
        #this is a repository where I've stored the script, README files, and the simulated datasets

#The overall fixed effect formula and landmark registration technique was derived from the 2014 Lopez-Duran paper: https://www.ncbi.nlm.nih.gov/pubmed/24754834
      #I will refer to it in notes throughout and provide the link to the paper above. 

#Below are all of the packages used in this script to analyze data, perform multiple imputation, plot, run statistics, etc.
#just use install.packages(" insert package name here ") if you do not already have any of these packages installed 

library(jomo)
library(VIM)
library(mice)
library(psych)
library(ggplot2)
library(tidyverse)
library(lme4)
library(nlme)
library(multcomp)
library(readr)
library(afex)

#load in my simulated cortisol dataset
urlfile="https://raw.githubusercontent.com/sarahglier/Neuroanalytics_Final_Project/master/cort_wide.csv"
cort_wide<-read_csv(url(urlfile))


#see the variable names - we have subj ID in first column and then 6 cortisol samples in the successive columns
head(cort_wide)

###Creating variable of data without the SubjID column before impute data
cort2 <- cort_wide[,-1]

#Several subjects have missing saliva samples (not enough saliva, contamination, etc.) 
#So I'll look at the patterns of missing data 
md.pattern(cort2)
##more helpful visual representation of pattern above
#Provides a histogram of the missing data.  Note that there are no samples with greater than 5% of data missing!
aggr_plot <- aggr(cort2, col=c('navyblue','red'), numbers=TRUE, 
                  sortVars=TRUE, labels=names(cort2), 
                  cex.axis=.7, gap=3, 
                  ylab=c("Histogram of missing data","Pattern")) 



##Imputing the missing data using mice package. Excellent documentation and specifics on the statistics in link below
#https://www.r-bloggers.com/imputing-missing-data-with-r-mice-package/ 
tempData <- mice(cort2,m=5,meth='pmm',seed=500)
  #m=5 refers to the number of imputed datasets. 5 is the default value.
  #meth='pmm' refers to the imputation method. Here we use predictive mean matching (pmm) as imputation method; type methods(mice) for a list of available imputation methods if you don't want pmm.
  #normally the default seed is set to seed=NA. BUT I changed this to 500 here (randomly chose it) because I want someone to be able to closely replicate what I get in this step
summary(tempData)  

#check the imputation method used for each variable
#should see nothing under cort_1 and cort_3 as they had no missing observations
tempData$meth

######Inspecting the distribution of original and imputed data#######
xyplot(tempData,cort_0 ~ cort_1 + cort_2 + cort_3 + cort_4 + cort_5,pch=18,cex=1)
    #use a scatterplot and plot cort_0 sample against all the other variables
    # magenta points (imputed) 
    # blue ones (observed)

stripplot(tempData, pch = 20, cex = 1.2)
    #Another useful visual take on the distributions can be obtained using the stripplot() function 
    #that shows the distributions of the variables as individual points
##########END plotting of imputed vs observed data points

#Now we can get back the completed dataset 
completedData <- mice::complete(tempData,1) 
    #I had to use mice:: to call its complete function because tidyverse (tidyr) has the same function and since I loaded it first, if I just put complete - it will not call the right function for me

###Export Completed Dataset with Subj ID back
cort_wide_2 <- cbind(id = cort_wide$`id`, completedData) 
    #I am calling the id column from the original loaded dataset and binding it to our complete imputed dataset
    #all within a new dataframe called cort_wide_2
    #now we have a dataset that looks exactly like the original except with imputed values

md.pattern(cort_wide_2)#Just to reassure ourselves that there are no missing values

###NOTE!!! Technically you do NOT need to impute ANY data when performing growth curve / multilevel modeling
  #I did this so I could also compare against other methods like the area under the curve analyses or a RM-ANOVA
      #both of which cannot handle missing values (they completely throw out an entire subject's data in the analysis if there's one missing sample)
  #I also used MICE as this is a coding / Neuroanalytics class and I thought it would be a useful exercise as part of this final project!


#The data is in "wide" format, where all the repeated measures have their own column.
    #we need to get this data into "long" format for our analysis where one column signifies time (e.g. cort_0,1,2)
    #and the other column indicates cortisol concentration values
    #this means that the id column will have REPEATED study ids 
    #we will use the reshape function to go from wide > long
cort_long <- reshape(data=cort_wide_2, #calling our wide dataframe 
                         timevar=c("time"), #creating a col called time which will contain the repeated meas names AKA column names (e.g. cort_0,1,2) 
                         idvar="id", #specifying the id col 
                         varying=c("cort_0","cort_1","cort_2","cort_3",
                                   "cort_4","cort_5"), #telling it where our repeated measure values are (cort concentrations)
                         direction="long", sep="_") #telling it we want our data to reshape to long

head(cort_long,12)
    #this 12 added will allow me to see the top 12 rows
    #see that we now have 3 columns: id, time, and cort
    #but it's not very organized
    #we can't see all of subject1's data

cort_long <- cort_long[order(cort_long$id,cort_long$time), ]
    #I'll order my dataset by id and time to organize it better
    #by writing the code this way it tells it to FIRST order by ID and THEN order by time
    #This way all values from one subject in order can be seen

head(cort_long)
    #now we see subject 1 and all cort conc @ times 0 - 5


#I am just cleaning up my workspace and reducing memory usage now that we are done imputing and data munging
  rm(aggr_plot, completedData, cort_wide, tempData, cort2)
  #rm() just removes the prior matrices and dataframes from memory

  
describe(cort_wide_2)
    #This is now going to let me look at the basic descriptives of my dataset (means, stand dev., etc).
    #Note here that I am calling the WIDE dataset 
    #I like to look at this because we can see a general trend from the mean column across time
    #There seems to be a curvilinear pattern which we expect the cortisol stress response to look like!

#Create a boxplot of this data to look at distribution of cortisol concentrations within each timepoint
ggplot(data=cort_long, aes(x=factor(time), y=cort)) + #x axis time, y axis cort concentration
  geom_boxplot(notch = TRUE) + #telling it we want a boxplot and telling notch=TRUE will allow us to see the MEDIAN in our boxplot where it "notches"
  stat_summary(fun.y="mean", geom="point", shape=23, size=3, fill="white") +  #adding summary statistics to my boxplot - the shape=23 will put a diamond in the boxplot where the mean is located and I also specified what color I wanted it to be - this can be changed
  labs(x = "Time", y = "Cortisol Concentration (ug/dL)") #setting axis titles
        #can see overlap of means and medians
        #there doesn't appear to be any crazy outlying concentrations anywhere
        #can also visualize this curvilinear trajectory across time that we noticed in the means

#Another pretty visualization of how the cortisol changes across the stressor / time
    #we can plot density curves for each sample
ggplot(data=cort_long, aes(x=cort)) + 
  geom_density(aes(group=factor(time), colour=factor(time), fill=factor(time)), alpha=0.3) #this is telling it that I want a separate denisty curve for each time point and I want them to be different colors

      
#One reason why RM_ANOVAs or standard multiple regression is sub-optimal to multilevel modeling is because it doesn't handle multicollinearity very well
    #example: a single subject's cort concentrations are likely to be very correlated with one another, on top of that there are likely to be sample time points (eg. cort 2 & 3) that are highly correlated
      #so we have colinearity within and between subjects
      #below, let's just take a look at the correlations between time points
pairs.panels(cort_wide_2[,c("cort_0","cort_1","cort_2","cort_3","cort_4","cort_5"),])
      #this function comes from the psych package


#NOW, I want to plot the trajectories of individuals' stress responses or a visual of intraindividual change in cortisol trajectories across time (stress)
    ggplot(data = cort_long, aes(x = time, y = cort, group = id)) + #call long dataset, x axis time, y axis cort conc, and each individual will have their own line
      geom_point(color="black") + #plot each point
      geom_line(color="black") +  #connect each point by a line
      xlab("Time") + #label x axis
      ylab("Cortisol Concentration (ug/dL)") + ylim(0,30) + #label y axis and set limits
      scale_x_continuous(breaks=seq(0,8,by=1)) 
      #we can see that there's a lot of individual variation in where someone starts, what pattern of response they might follow, and where they end

 

###LANDMARK REGISTRATION (LR)
    #in the original Lopez-Duran 2014 paper, they had a small dataset and manually plotted and scanned for peak cortisol concentrations in excel to identify the "peak landmark"
    #Then they made a new column called "newPeakTime" which was the landmark registration step - identifying each individual's peak
    #This is incredibly time inefficient if you have a large dataset - below we will automate finding the peak for each individual across time and storing that time in a new column
            
    LR.list <- colnames(cort_wide_2[2:7])[max.col(cort_wide_2[2:7],ties.method="first")]
          #LR.list <- colnames(cort_wide_2[2:7])[apply(cort_wide_2[2:7],1,which.max)]  # this also works, but it will print BOTH cols if there are 2 maxima AKA a "peaked plateau" - I only want the first peak not cort 2 and 3, etc.
    LR.identity <- data.frame(Reduce(rbind, LR.list)) #making this list into a dataframe
          names(LR.identity)[1]<-"LR.character" 
    #now recoding the "cort 2/3/4" time points to their actual times and storing in new column (LR.time) #NOTE - simulated data will not have any peaks outside of cort 2 or 3
          LR.identity$LR.time[LR.identity$LR.character=="cort_0"] <- 0
          LR.identity$LR.time[LR.identity$LR.character=="cort_1"] <- 1
          LR.identity$LR.time[LR.identity$LR.character=="cort_2"] <- 2
          LR.identity$LR.time[LR.identity$LR.character=="cort_3"] <- 3
          LR.identity$LR.time[LR.identity$LR.character=="cort_4"] <- 4
          LR.identity$LR.time[LR.identity$LR.character=="cort_5"] <- 5
    #now I want to bind id with this LR dataset
          LR <- cbind(id = cort_wide_2$'id', LR.identity)
    #now I want to merge the long dataset and this dataset by id - the LR will be repeated long ways to match subject id
          long_LR <- merge(cort_long, LR, by="id", all=T) 
  #Creation of the adjusted time variable is as follows and will be called time_to_peak; time_to_peak = (LR.time - Time)*-1
          long_LR$time_to_peak <- ((long_LR$LR.time - long_LR$time)*-1)
              #this has now created a time variable that is shifted so each individual's peak lands on the same time point
  #Now we will create 2 more time variables: (1) time before individual peak OR TBP (2) time after individual peak or TAP
      #this can be done using a conditional formulas:
        #(1) TBP 
          #If time_to_peak <0 then TBP = time_to_peak 
          #Else TBP = 0
          long_LR$TBP <- ifelse(long_LR$time_to_peak < 0, long_LR$time_to_peak, ifelse(long_LR$time_to_peak >= 0, 0, NA))
          
        #(2) TAP
          #If time_to_peak > 0 then TAP = time_to_peak
          #Else TAP = 0
          long_LR$TAP <- ifelse(long_LR$time_to_peak > 0, long_LR$time_to_peak, ifelse(long_LR$time_to_peak <= 0, 0, NA))

#Now that we've finished coding our time variables, we can add in our demographic / ID variables to our dataset
    #This simulated example data contains random ages (between 8-18) and sex (M=0 v F=1)
          urlfile="https://raw.githubusercontent.com/sarahglier/Neuroanalytics_Final_Project/master/id_data.csv"
          id_data<-read_csv(url(urlfile))
    #sex & age will be added as fixed effects or time invarying covariates to the analysis
    #first we have to merge it with our long dataset
          cort_long_model <- merge(long_LR, id_data, by="id", all=T)

########
#MODEL TIME - all fixed effect models specified below are derived from the Lopez-Duran 2014 paper
          str(cort_long_model) #checking the structure of our variables - problematic b/c it thinks sex is an integer and not a factor
                  cort_long_model$sex <- factor(cort_long_model$sex) #we've now changed sex to a factor not int. Can run above line again if want to confirm this
      
        #A note on the general form of the equation:
        # Y (dependent var / outcome) ~ X (predictor / independent variable) + X1*X2 (these are interaction terms you can specify between your predictors) + X3 (some covariate you want to control for)
        # the ~ represents the = or Y "regressed" on X
        # the 1|id notes the random effect - this can be expanded... for example say you have randomized participants to protocols/days of the week or some study design you want to account for.
                        #the way to write in that random effect would be: random = ~ 1 + study_design_random|id
                                #it can be hard to conceptualize what constitutes a random or fixed effect.  Below are 2 resources that have helped me.
                                #Drs. Curran and Bauer @ UNC: https://curranbauer.org/wp-content/uploads/2016/SRA/Curran-Bauer-SRA-4APRIL2016-2-up.pdf
                                #https://stats.stackexchange.com/questions/4700/what-is-the-difference-between-fixed-effect-random-effect-and-mixed-effect-mode
                                #http://statweb.stanford.edu/~jtaylo/courses/stats203/notes/fixed+random.pdf
      model.sex <- lme(cort ~ 1 + sex + TBP + TAP + sex*TBP + sex*TAP, 
                        random = ~ 1|id,
                        data = cort_long_model) #note if you have NAs, just add in na.action = "na.exclude" to this formula
      summary(model.sex)  
          #The way to interpret the fixed effects is as follows:
              #Bo = intercept (peak) ; the 1 in the above lme model
              #B1 = impact of sex on the peak ; sex term in ^^model
              #B2 = is the activation slope ; TBP term above
              #B3 = is the regulation slope ; TAP term above
              #B4 = is the impact of sex on activation slope ; sex*TBP term
              #B5 = is the impact of sex on regulation slope ; sex*TAP term
              #random effect - I allow intercept & slope to vary by participant id as not all individuals will show the same starting place (intercept) or same rates of change (slope)
  
      #Now lets add in our covariate age!
      model.covar <- lme(cort ~ 1 + sex + TBP + TAP + sex*TBP + sex*TAP + age, 
                         random = ~ 1|id,
                         data = cort_long_model) #note if you have NAs, just add in na.action = "na.exclude" to this formula
      summary(model.covar)
      
##Let's visualize and plot predicted trajectories based on the above model
      #getting predicted scores for individuals
      cort_long_model$predict_sex <- predict(model.sex)
      
      #plotting predicted trajectories - intraindividual change trajetories plotted
      ggplot(data = cort_long_model, aes(x = time, y = predict_sex, group = id)) +
        #geom_point(color="black") + 
        geom_line(color="black")


      
      
#COMPARISONS - what if we were JUST to run a repeated measure ANOVA?
      ####TIME rescaling
      #let's 1st rescale the time variable - so instead of going from 0-5, it's rescaled from 0-1
      cort_long_model$time_sc <- (cort_long_model$time - 0)/5
          #Here I am calling our original time column, subtracting 0 then dividing by 5
          #I then store this new scaled time in a new column in our dataframe called time_sc
      head(cort_long)
          #can now see that time goes from 0-1 in 0.2 increments
      #now, lets make another quick time transformation - a timesq variable so we can approximate a quadratic fit (which we know exists from examining the raw data)
      cort_long_model$timesq <- cort_long_model$time_sc^2
      
  #RM ANOVA Analysis 
      require(afex)
      
      #RM_ANOVA below using the quadratic fit (time squared variable)
              #Note the RM-ANOVA techniques have huge flaws & this will not work if you are applying your own dataset that contains NA or missing values
                  #if you do want to compare with your own dataset - use na.omit(data) on your dataframe before running below if you have NA values 
      ez.glm("id", "cort", cort_long_model, within=c("timesq"), 
             between= c("sex", "age"))

  #Standard Multilevel Model Approach WITHOUT Landmark Registration
      cort_quad <- lme(cort ~ 1 + sex*timesq + age,  
                        random= ~1 + timesq|id,
                        data = cort_long_model)
      summary(cort_quad)    
      

      
      
      

      
      
      
      
      
      
              
      
      
      
      
      
    

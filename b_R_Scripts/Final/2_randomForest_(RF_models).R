#------------------------------------------ RANDOMFOREST LIPIDOMIC DATA ------------------------------------------#

# Rheumatoid arthritis is a inflammatory disease that characterized for not having a resolution phase in the 
# inflammatory procces. One of the most used treatments is modifying antirheumatic drugs (DMARDs), that is not 
# always efective. 

# We have a dataset of the lipid profiles of rheumatoid arthritis patients that responded or not to the treatment.

# This is the first approach to creat a model using "randomForest" (decision trees).

#---> LIBRARY LOAD:

library(randomForest)
library("mltools")
set.seed(415) # To get same results even with the random part.

#---> DATA MANIPULATION: 

# TRAINING SET:

# Data uses to create the model!

# Open the txt file with the profiles information. Make sure that the path is correct:

# The dataset consist in a tab-delimited file in .txt format with the follow specifications: 
# Columns: The different lipid mediators plus a column called "responses" that contains information about the 
# "Responder" and "Non_Responder". 
# Rows: The different samples (each patient data).

# See a_Toy_Data/2_randomForest_(RF_models)/2_randomForest_(RF_models)_toy_data.txt

lm_profiles <- read.table(
  file = "GitHub/2018_Machine_Learning_MTX_treatment_in_RA_patients/a_Toy_Data/2_randomForest_(RF_models)/2_randomForest_(RF_models)_toy_data.txt",
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t")

lm_profiles_scale <- as.data.frame(scale(lm_profiles[, -1], center = FALSE, scale = TRUE))

# VALIDATION SET:

# Data use to test the model (independent cohort)!

# The dataset consist in a tab-delimited file in .txt format with the follow specifications: 
# Columns: The different lipid mediators plus a column called "responses" that contains information about the 
# "Responder" and "Non_Responder". 
# Rows: The different samples (each patient data).




#---> DATA PREPARATION TRAINING DATA: 

# Add the classification variable to the data frame (Responder and non responder):

# Getting the explanatory (x) and response (y) variable. By explanatory, it means all the data that can explain why a
# patient response or not to the treatment (the lipid meadiator profiles) and the response variable is if the 
# patients response or not to the treatment. In random Forest you have to create a formula where the response 
# variable is explain in terms of the explanatory variable (responses ~ everything else).

lm_profiles_scale$responses <- lm_profiles$responses

# Make sure that column names do not represent a problem to randomForest making them a valid name to R.

names(lm_profiles_scale) <- make.names(names(lm_profiles_scale))

# Separates the profiles data by lipid mediators types. 

# By Substrates:

dha <- lm_profiles_scale[ ,c(1:24, 56)]
n_three_DPA <- lm_profiles_scale[ , c(25:34, 56)]
epa <- lm_profiles_scale[ , c(35:37, 56)]
aa <- lm_profiles_scale[ , c(38:55, 56)]

#---> MACHINE LEARNING (randomForest R):

# In Random Forests the idea is to decorrelate the several trees which are generated on the different bootstrapped 
# samples from training Data and then reduce the variance in the trees by averaging them.

# Averaging the trees also improve the perfomance of decision trees on Test Set and eventually avoid overfitting.

# The idea is to build lots of trees in such a way to make the correlation between the trees smaller.

# BEST MTRY:
# mtry is the number of variables available for splitting at each tree node. Random Forest creates several trees, 
# each one using different variables to create the best version of it. With mtry we can define how many variables 
# the data is split to create the different trees.
# More: https://stats.stackexchange.com/questions/102867/random-forest-mtry-question

# In this case we defined the lipid mediators to create a loop to define which mtry is the best one for our model. 

oob_error <- double(ncol(lm_profiles_scale) - 1) #Define number of variable. -1 is because the last column is responses.

# Loop to select the best mtry. 

for (mtry in 1:(ncol(lm_profiles_scale) - 1)) {
  
  #importance = TRUE creates the graph of the important variables. 
  rf_lm_profiles_scales <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = mtry, 
                                        importance = TRUE, ntree = 10000)
  
  oob_error[mtry] <- 100 - ((rf_lm_profiles_scales$err.rate[10000])*100)
  
}

# Define the best mtry according to the best prediction value. 

final_mtry <- which.max(oob_error)

# Run the model again with the right mtry value. 

rf_lm_profiles_final <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = final_mtry, 
                                     importance = TRUE, ntree = 10000)

# Save relevan Plots:

pdf(file = "GitHub/2018_Machine_Learning_MTX_treatment_in_RA_patients/c_Expected_Output/2_randomForest_(RF_models)/2_RF_lmprofiles_scale.pdf",
    width = 14, height = 10, onefile = TRUE)
varImpPlot(rf_lm_profiles_final, sort = TRUE, main = "lm_profiles") # Importance of the variable for the model. 
plot(rf_lm_profiles_final, main = "lm_profiles") # Decreasing of the error base on the number of tres. 
legend("topright", inset=.05, title="Curves:",
       c("Non responder", "OOB error", "Responder"), fill=c("red", "black", "green"),
       text.font = 3, cex = 1)
dev.off()

#---> VALIDATION TEST:

# "prediction" takes the created models and the validation dataset to try to predict which samples belongs to the
# responder and non-responder. 



# % OF ACCURACY:

# An error estimate is made for the cases which were not used while building the tree. That is called an 
# OOB (Out-of-bag) error estimate which is mentioned as a percentage.

# In this case, then, if a model has 1% error will means that we can predict with 99% accuracy.

no_oob_error_table <- data.frame(groups = "scalated lm profiles",
                                 percentage_accuracy = 100 - ((rf_lm_profiles_final$err.rate[10000])*100),
                                 stringsAsFactors = FALSE)

# Save the models as an R object:
# To avoid running all the scripts to obtain the models, it can be saved as R objects. If you want to use it, 
# you can call it using "readRDS". 

saveRDS(rf_lm_profiles_final, 
        file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/4_random_forest/models_RF/RF_lmprofiles_scale.R", 
        ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)

#---> RF PER GROUP: 

# Run the same analysis but automatically for the rest of the subgroups. 

# Create a list with the names of all the subgroups:
groups <- list(dha, n_three_DPA, aa, resolvins_d, epa, protectins, pctr, maresins, mctr, rvt, 
               rvd, pd, mar, lx, ltb, lt, pg, tx, dha_n_dpa_epa)

# Create a vector with the names associated to all the elements in the list: 
names(groups) <- c("DHA", "n-3 DPA", "AA", "Resolvins_d", "EPA", "Protectins", 
                   "PCTRS", "Maresins", "MCTR", "RVT", "D-Series Resolvins", "PD", "Mar", "LX", 
                   "LTB", "Leukotrienes", "Prostaglandins", "TX", "DHA_N-3DPA_EPA")

# The loop goes through all the elements in gropus and creates a model for each of them. The model can not work with
# only one row, so the "if" makes sure that only the subgroups with more than one row are analyzed. 

for (lm in 1:length(groups)) {
  
  if (ncol(groups[[lm]]) > 2) {
    
    # Identifying the best mtry:
    
    oob_error_all <- double(ncol(groups[[lm]]) - 1)
    
    for (mtry_all in 1:(ncol(groups[[lm]]) - 1)) {
      
      random_forest <- randomForest(responses ~ ., data = groups[[lm]], mtry = mtry_all, 
                                            importance = TRUE, ntree = 10000)
      oob_error_all[mtry_all] <- 100 - ((random_forest$err.rate[10000])*100)
      
    }
    
    # Creating the final model: 
    
    final_mtry_all <- which.max(oob_error_all)
    random_forest_final <- randomForest(responses ~ ., data = groups[[lm]], mtry = final_mtry_all, 
                                         importance = TRUE, ntree = 10000)
    
    # Model Graphs:
    
    pdf(file = paste("C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/4_random_forest/graphs_RF/RF_",
                     names(groups)[[lm]], ".pdf", sep = ""), 
        width = 14, height = 10, onefile = TRUE)
    varImpPlot(random_forest_final, sort = TRUE, main = names(groups)[[lm]]) # Importance of the variable for the model. 
    plot(random_forest_final, main = names(groups)[[lm]]) # Decreasing of the error base on the number of tres. 
    legend("topright", inset=.05, title="Curves:",
           c("Non responder", "OOB error", "Responder"), fill=c("red", "black", "green"),
           text.font = 3, cex = 1)
    dev.off()
    
    # Accuracy table: 
    
    accuracy_table <- data.frame(groups = names(groups)[[lm]],
                                 percentage_accuracy = 100 - ((random_forest_final$err.rate[10000])*100),
                                 stringsAsFactors = FALSE)
    
    no_oob_error_table <- rbind(no_oob_error_table, accuracy_table)
    
    # Save final model:
    
    saveRDS(random_forest_final, 
            file = paste("C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/4_random_forest/models_RF/RF_",
                         names(groups)[[lm]], ".R", sep = ""),
            ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
  
  }
  
  else { next }
  
}
    
# Write table with accuracy values: 

write.table(no_oob_error_table, 
            file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/4_random_forest/no_oob_error_table.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)  

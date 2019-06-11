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

#---> INPUT AND OUTPUT:

# In this section please specify where are the input files and where you want to save the output files.
# In the input and output variable you can see what is the path expected from the user to write.

input <- "C:/Users/hhy270/Documents/GitHub/2018_Machine_Learning_MTX_treatment_in_RA_patients/a_Toy_Data/2_randomForest_(RF_models)/"
output <- "C:/Users/hhy270/Documents/GitHub/2018_Machine_Learning_MTX_treatment_in_RA_patients/c_Expected_Output/2_randomForest_(RF_models)/"

# !!!! IMPORTANT: For this script to work the training dataset has to be called: 2_randomForest_(RF_models)_toy_data.txt
# !!!! IMPORTANT: For this script to work the validation dataset has to be called: 2_randomForest_(RF_models)_data_validation.txt

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
  file = paste(input, "2_randomForest_(RF_models)_toy_data.txt", sep = ""),
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t")

lm_profiles_scale <- as.data.frame(scale(lm_profiles[, -1], center = FALSE, scale = TRUE))

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

# VALIDATION SET:

# Data use to test the model (independent cohort)!

# The dataset consist in a tab-delimited file in .txt format with the follow specifications: 
# Columns: The different lipid mediators plus a column called "responses" that contains information about the 
# "Responder" and "Non_Responder". 
# Rows: The different samples (each patient data).

val_lm_profiles <- read.table(
  file = paste(input, "2_randomForest_(RF_models)_data_validation.txt", sep = ""), 
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t")

val_lm_profiles_scale <- as.data.frame(scale(val_lm_profiles[, -1], center = FALSE, scale = TRUE))
names(val_lm_profiles_scale) <- make.names(names(val_lm_profiles_scale))

# By Substrates:

val_dha <- val_lm_profiles_scale[ ,c(1:24)]
val_n_three_DPA <- val_lm_profiles_scale[ , c(25:34)]
val_epa <- val_lm_profiles_scale[ , c(35:37)]
val_aa <- val_lm_profiles_scale[ , c(38:55)]

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
  
  # NOTE: 
  # importance = TRUE creates the plot of the important variables, that can gave us an idea, based on the
  # decrease of the accuracy of the models, what lipid mediators are contributing to make a better model. 
  
  rf_lm_profiles_scales <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = mtry, 
                                        importance = TRUE, ntree = 10000)
  
  oob_error[mtry] <- 100 - ((rf_lm_profiles_scales$err.rate[10000])*100)
  
}

# Define the best mtry according to the best prediction value. 

final_mtry <- which.max(oob_error)

# Run the model again with the right mtry value. 

rf_lm_profiles_final <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = final_mtry, 
                                     importance = TRUE, ntree = 10000)

# Save relevant Plots:

pdf(file = paste(output, "2_RF_lmprofiles_scale.pdf", sep = ""), 
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

pred_rf_lm_profiles <- as.data.frame(predict(rf_lm_profiles_final, val_lm_profiles_scale))
names(pred_rf_lm_profiles) <- c("prediction")

# In order to further evaluate the predictivenss of this approach we next calculated  Matthews correlation 
# coefficient (MCC), which represents the accuracy of the model at predicting outcome. Very helpful when you
# have imbalance data. 

mcc_value = mcc(preds = pred_rf_lm_profiles$prediction, actuals = val_lm_profiles$responses) # From the mltools package.

# Creates a table with all the models, the %CC and the MCC.

# An error estimate is made for the cases which were not used while building the tree. That is called an 
# OOB (Out-of-bag) error estimate which is mentioned as a percentage.

# In this case, then, if a model has 1% error will means that we can predict with 99% accuracy.

no_oob_error_table <- data.frame(groups = "scalated lm profiles",
                                 percentage_accuracy = 100 - ((rf_lm_profiles_final$err.rate[10000])*100),
                                 MCC = mcc_value,
                                 stringsAsFactors = FALSE)

# Save the models as an R object:
# To avoid running all the scripts to obtain the models, it can be saved as R objects. If you want to use it, 
# you can call it using "readRDS". 

saveRDS(rf_lm_profiles_final, 
        file = paste(output, "2_RF_lmprofiles_scale.R", sep = ""),  
        ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)

#---> RF PER GROUP: 

# Create a list with the names of all the subgroups:

# Modifying the original lipid mediator file, including other data (such as clinical scores) and dividing the table
# in different sections is possible to create automatically all the models you want. YOU NEED TO MAKE SURE WHATSOEVER
# that the "groups" list and the names of the "groups" list is updated. 

groups <- list(dha, n_three_DPA, epa, aa) # Update in case you want to create other models. 
val_groups <- list(val_dha, val_n_three_DPA, val_epa, val_aa) # Update in case you want to create other models.

# Create a vector with the names associated to all the elements in the list: 

names(groups) <- c("DHA", "n-3 DPA", "EPA", "AA") # Update in case you want to create other models. 

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
    
    pdf(file = paste(output, "2_RF_",
                     names(groups)[[lm]], ".pdf", sep = ""), 
        width = 14, height = 10, onefile = TRUE)
    varImpPlot(random_forest_final, sort = TRUE, main = names(groups)[[lm]]) # Importance of the variable for the model. 
    plot(random_forest_final, main = names(groups)[[lm]]) # Decreasing of the error base on the number of tres. 
    legend("topright", inset=.05, title="Curves:",
           c("Non responder", "OOB error", "Responder"), fill=c("red", "black", "green"),
           text.font = 3, cex = 1)
    dev.off()
    
    # Validation: 
    
    pred_val <- as.data.frame(predict(random_forest_final, val_groups[[lm]]))
    names(pred_val) <- c("prediction")
    
    mcc_val = mcc(preds = pred_val$prediction, actuals = val_lm_profiles$responses)
    
    # Accuracy table: 
    
    accuracy_table <- data.frame(groups = names(groups)[[lm]],
                                 percentage_accuracy = 100 - ((random_forest_final$err.rate[10000])*100),
                                 MCC = mcc_val,
                                 stringsAsFactors = FALSE)
    
    no_oob_error_table <- rbind(no_oob_error_table, accuracy_table)
    
    # Save final model:
    
    saveRDS(random_forest_final, 
            file = paste(output, "2_RF_",
                         names(groups)[[lm]], ".R", sep = ""),
            ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
  
  }
  
  else { next }
  
}
    
# OUTPUT: 

# Besides the functional models, the OUTPUT of this script is a table with Accuracy score of the model and the
# validation of the model based on the MCC value. 

# The resulting plots also gave us an idea of the performance of the models, that includes the "importance" analysis,
# and the improvement of the model based on the number of tress used to create it. 

write.table(no_oob_error_table, 
            file = paste(output, "2_accuracy_table_RF.txt", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)  

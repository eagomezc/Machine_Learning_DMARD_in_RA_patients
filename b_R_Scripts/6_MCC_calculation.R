#------------------------------- randomForest VALIDATION AND ROC CURVES  --------~----------------------------------#

# Rheumatoid arthritis is a inflammatory disease that characterized for not having a resolution phase in the inflammatory
# procces. One of the most used treatments is Methotrexate (MTX), that is not always efective. 

# We have a dataset of the lipid profiles of rheumatoid arthritis patients that responded or not to the treatment.

# This is the first approach to compare all the accuracy results from all the models using Matthews Correlation 
# Coefficient (MCC).

#---> LIBRARY LOAD:

library("mltools")

#---> DATA MANIPULATION: 

# Create data frames where the validation results will be saved: 

mcc_tables_val <- data.frame(groups = character(),
                             mcc_value = numeric())

mcc_tables_tra <- data.frame(groups = character(),
                             mcc_value = numeric())

# Creates a vector with all the validations sets:

groups <- c("n-3 DPA", "AA", "Resolvins_d", "Protectins", "PCTRS", "Maresins", "MCTR", "RVT", "D-Series Resolvins", 
            "PD", "LX", "LTB", "Leukotrienes", "Prostaglandins", "lmprofile", "dha", "epa", "dha_dpa_epa")

#---> MODEL VALIDATION:

# With Validation set: 

# Creates a loop that open all the validation results for each model and then calculates the MCC (Mathews Coefficient
# Correlation).

for (lm in 1:length(groups)) {
  
# Open the corresponding prediction file (CLASSYFIRE):  

predictions <- read.table(
  file = paste("C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/6_MCC_calculations/validation sets/pred_",
               groups[lm], ".txt", sep = ""),
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t")

# Change the name of the column from the prediction file:

names(predictions) <- c("prediction")

# Creates a data frame that contains the real classification of the validation set samples:

response_val_set <- data.frame(row.names = row.names(predictions),
                               responses = c(rep("responder", 37), rep("non-responder", 8)))

# Creates a data frame with the MCC value of the specific model. 

mcc_result <- data.frame(groups = groups[lm],
                         mcc_value = mcc(preds = predictions$prediction,        # This bit is the code that calculates
                                         actuals = response_val_set$responses)) # the MCC value. 

# Add the result to the final table: 

mcc_tables_val <- rbind(mcc_tables_val, mcc_result)

# Open the corresponding prediction file (RANDOM FOREST):  

predictions <- read.table(
  file = paste("C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/6_MCC_calculations/validation sets/rf_pred_",
               groups[lm], ".txt", sep = ""),
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t")

# Change the name of the column from the prediction file:

names(predictions) <- c("prediction")

# Creates a data frame that contains the real classification of the validation set samples:

response_val_set <- data.frame(row.names = row.names(predictions),
                               responses = c(rep("responder", 37), rep("non-responder", 8)))

# Creates a data frame with the MCC value of the specific model. 

mcc_result <- data.frame(groups = paste("rf_", groups[lm], sep = ""),
                         mcc_value = mcc(preds = predictions$prediction,        # This bit is the code that calculates
                                         actuals = response_val_set$responses)) # the MCC value. 

# Add the result to the final table: 

mcc_tables_val <- rbind(mcc_tables_val, mcc_result)

}

# With Training set: 

# Do the same that the previous loop except that it uses the training set. It's only to confirm that the model is
# working properly. 

for (lm in 1:length(groups)) {
  
  predictions_tra <- read.table(
    file = paste("C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/6_MCC_calculations/training_sets/pred_",
                 groups[lm], "_training.txt", sep = ""),
    header = TRUE,
    row.names = 1, # Specify that the first column is row names. 
    sep = "\t")
  
  names(predictions_tra) <- c("prediction")
  
  response_tra_set <- data.frame(row.names = row.names(predictions_tra),
                                 responses = c(rep("responder", 30), rep("non-responder", 22)))
  
  
  mcc_result_tra <- data.frame(groups = groups[lm],
                           mcc_value = mcc(preds = predictions_tra$prediction, 
                                           actuals = response_tra_set$responses))
  
  mcc_tables_tra <- rbind(mcc_tables_tra, mcc_result_tra)
  
  predictions_tra <- read.table(
    file = paste("C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/6_MCC_calculations/training_sets/rf_pred_",
                 groups[lm], "_training.txt", sep = ""),
    header = TRUE,
    row.names = 1, # Specify that the first column is row names. 
    sep = "\t")
  
  names(predictions_tra) <- c("prediction")
  
  response_tra_set <- data.frame(row.names = row.names(predictions_tra),
                                 responses = c(rep("responder", 30), rep("non-responder", 22)))
  
  mcc_result_tra <- data.frame(groups = paste("rf_", groups[lm], sep = ""),
                           mcc_value = mcc(preds = predictions_tra$prediction,        
                                           actuals = response_tra_set$responses)) 
  
  mcc_tables_tra <- rbind(mcc_tables_tra, mcc_result_tra)
  
}

write.table(mcc_tables_val, 
           file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/6_MCC_calculations/all_val_mcc.txt",
           sep = "\t",
           quote = FALSE,
           row.names = TRUE)

write.table(mcc_tables_tra, 
            file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/6_MCC_calculations/all_tra_mcc.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)
#------------------------------------------ CLASSYFIRE LIPIDOMIC DATA --------------------------------------------#

# Rheumatoid arthritis is a inflammatory disease that characterized for not having a resolution phase in the 
# inflammatory procces. One of the most used treatments is modifying antirheumatic drugs (DMARDs), that is not 
# always efective. 

# We have a dataset of the lipid profiles of rheumatoid arthritis patients that respond or not to the MTX treatment.
# This is the first approach to create a model using machine learning that predicts if a patient will react or not 
# to the treatment with DMARDs.

#---> LIBRARY LOAD:

library(classyfire)
library(mltools)
library(pROC)

#---> INPUT AND OUTPUT:

# In this section please specify where are the input files and where you want to save the output files.
# In the input and output variable you can see what is the path expected from the user to write.

input <- "C:/Users/hhy270/Documents/GitHub/2018_Machine_Learning_MTX_treatment_in_RA_patients/a_Data/1_classyfire_(SVM_models)/"
output <- "C:/Users/hhy270/Documents/GitHub/2018_Machine_Learning_MTX_treatment_in_RA_patients/c_Expected_Output/1_classyfire_(SVM_models)/"

# !!!! IMPORTANT: For this script to work the training dataset has to be called: 1_classyfire_(SVM_models)_data.txt
# !!!! IMPORTANT: For this script to work the test dataset has to be called: 2_classyfire_(SVM_models)_data_validation.txt

#---> DATA MANIPULATION: 

# TRAINING SET:

# Data uses to create the model!

# Open the txt file with the profiles information. Make sure that the path is correct:

# The dataset consist in a tab-delimited file in .txt format with the follow specifications: 
# Columns: The different samples (each patient data)
# Row number 1: Class row that contains the word "Responder" if the patient respond to treatment and "Non_Responder"
# otherwhise.
# Following rows: All the lipid mediators used to create the model.

# See a_Data/1_classyfire_(SVM_models)/1_classyfire_(SVM_models)_data.txt

lm_profiles <- read.table(
 file = paste(input, "1_classyfire_(SVM_models)_data.txt", sep = ""), 
 header = TRUE,
 row.names = 1,
 sep = "\t")

# Separates the profiles data by lipid mediators types. 

# by Substrates:

dha <- lm_profiles[2:24, ]
n_three_DPA <- lm_profiles[25:34, ]
epa <- lm_profiles[35:37, ]
aa <- lm_profiles[38:55, ]

# VALIDATION SET:

# Data use to test the model (independent cohort)!

# Open the txt file with the profiles information. Make sure that the path is correct:

# The dataset consist in a tab-delimited file in .txt format with the follow specifications: 
# Columns: The different samples (each patient data)
# Row number 1: Class row that contains the word "Responder" if the patient respond to treatment and "Non_Responder"
# otherwhise.
# Following rows: All the lipid mediators used to create the model (NOTE: They have to be in the same order as
# the training dataset).

# See a_Data/1_classyfire_(SVM_models)/2_classyfire_(SVM_models)_data_validation.txt

val_lm_profiles <- read.table(
  file = paste(input, "2_classyfire_(SVM_models)_data_validation.txt", sep = ""), 
  header = TRUE,
  row.names = 1,
  sep = "\t")

# by Substrates:

val_dha <- val_lm_profiles[2:24, ]
val_n_three_DPA <- val_lm_profiles[25:34, ]
val_epa <- val_lm_profiles[35:37, ]
val_aa <- val_lm_profiles[38:55, ]

#---> DATA PREPARATION TRAINING DATA:

# Getting the explanatory (x) and response (y) variable. By explanatory, it means all the data that can explain 
# why a patient respond or not to the treatment (the lipid meadiator profiles), and the response variable is if 
# the patients respond or not to the treatment. 

# Transpose columns and rows:
lm_profiles_transpose <- t(lm_profiles)

# Explanatory and Response variable:
response <- data.frame(row.names = row.names(lm_profiles_transpose), # Create the Response variable.
                       responses = lm_profiles_transpose[, 1]) 

response_matrix <- as.matrix(response)

explanatorys <-  lm_profiles_transpose[, -1] # Delete the first column that contains the Response variable.

# Because the response variable was in the data.frame all the elements were saved as factors, and classyfire 
# requieres a numeric matrix. here we make it:

explanatory <- matrix(as.numeric(unlist(explanatorys)),nrow=nrow(explanatorys))
rownames(explanatory) <- rownames(explanatorys)
colnames(explanatory) <- colnames(explanatorys)

# If you are working with machine learning, the best method of scalation is standarization. 
# Scale data prevents that the modelfrom being based on variables with high normal values.

explanatory_scale <- scale(explanatory, center = FALSE, scale = TRUE)
explanatory_scale[is.na(explanatory_scale)] <- 0

#---> DATA PREPARATION TEST DATA:

# Similar process to the training dataset.

# Transpose columns and rows:
val_lm_profiles_transpose <- t(val_lm_profiles)

# Explanatory and Response variable:
val_response <- data.frame(row.names = row.names(val_lm_profiles_transpose), # Create the Response variable.
                       responses = val_lm_profiles_transpose[, 1]) 

val_explanatorys <-  val_lm_profiles_transpose[, -1] # Delete the first row that contains the Response variable.

# Because the response variable was in the data.frame all the elements were saved as factors, and classyfire 
# requieres a numeric matrix. here we make it:

val_explanatory <- matrix(as.numeric(unlist(val_explanatorys)),nrow=nrow(val_explanatorys))
rownames(val_explanatory) <- rownames(val_explanatorys)
colnames(val_explanatory) <- colnames(val_explanatorys)

# Since  the model is created using scalation and transpose data, the test datasets has to be scaled and 
# transpose as well. 

validation_scale <- scale(val_explanatory, center = FALSE, scale = TRUE) 
validation_scale[is.na(validation_scale)] <- 0

#---> MACHINE LEARNING (Classyfire R): 

# Classyfire uses Support Vector Machine to create the Machine Learning Model. SVM classifies, makes a reggresion,
# and creates a novelty detection for the creation of the model. 

# The idea is to create several models and see which one fits the best. The models will be based on the whole
# lipid profiles and the different groups based on substrates. 

# "cfBuild" to create the SVM:

support_lmprofiles_scale <- cfBuild(explanatory_scale, response_matrix, 
                                    bootNum = 70,ensNum = 70, cpus = 4) 

# Expositional figures to get info of the model: 

ggClassPred(support_lmprofiles_scale, displayAll = TRUE, fillBrewer = TRUE, showText = TRUE) 
ggEnsTrend(support_lmprofiles_scale, ylims = c(50,100))

# Model accuracy: 

getAvgAcc(support_lmprofiles_scale)$Test # Get the %CC (Overall percentage of correctly classified test objects)
getConfMatr(support_lmprofiles_scale) # Get a table of the consensus classification of the best model. 

#---> MODEL VALIDATION: 

# "cfPredict" takes the created models and the test dataset to try to predict which samples belongs to the
# responder and non-responder. It creates a data frame with the identifications and % of accuracy. 

prediction_validation <- cfPredict(support_lmprofiles_scale, validation_scale)
names(prediction_validation) <- c("prediction", "Coef Score")

# For ROC curves we need the likehood of each sample to belong to one group of another. We add two new columns for 
# that purpose:

prediction_validation$Responder[prediction_validation$prediction == "Responder"] <- 
  prediction_validation$`Coef Score`[prediction_validation$prediction == "Responder"]/100 

prediction_validation$Responder[prediction_validation$prediction == "Non_Responder"] <- 
  1 - (prediction_validation$`Coef Score`[prediction_validation$prediction == "Non_Responder"]/100)

prediction_validation$Non_Responder <- 1 - prediction_validation$Responder

# In order to further evaluate the predictivenss of this approach we next calculated  Matthews correlation 
# coefficient (MCC), which represents the accuracy of the model at predicting outcome. Very helpful when you
# have imbalance data. 

mcc_value = mcc(preds = prediction_validation$prediction, actuals = val_response$responses) # From the mltools package.

# ROC curve calculation:

roc_value = roc(val_response$responses, prediction_validation$Non_Responder)

# Creates a table with all the models, the %CC and the MCC.

accuracy_table <- data.frame(groups = "scalated lm profiles",
                             percentage_accuracy = getAvgAcc(support_lmprofiles_scale)$Test,
                             MCC = mcc_value,
                             AUC = roc_value$auc,
                             stringsAsFactors = FALSE)

# Save the models as an R object:
# To avoid running all the scripts to obtain the model, it can be saved as R objects. If you want to use it, 
# you can call it using "readRDS". 

# Make sure that the you specify the path were you want to save your model: 

saveRDS(support_lmprofiles_scale, 
        file = paste(output, "1_SVM_lmprofiles_scale.R",sep = ""),
        ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)

#---> SVM PER GROUP: 

# Run the same analysis but automatically for the rest of the groups by Substrates. 

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
  
  if (nrow(groups[[lm]]) > 1) {
    
  # Training data set:
  
  transpose <- matrix(as.numeric(t(groups[[lm]])),nrow=nrow(t(groups[[lm]]))) # Transpose and create matrix.
  rownames(transpose) <- rownames(t(groups[[lm]]))
  colnames(transpose) <- colnames(t(groups[[lm]]))
  
  scale <- scale(transpose, center = FALSE, scale = TRUE) # All the model will use scale data
  scale[is.na(scale)] <- 0
  
  # test dataset: 
  
  val_transpose <- matrix(as.numeric(t(val_groups[[lm]])),nrow=nrow(t(val_groups[[lm]]))) # Transpose and create matrix.
  rownames(val_transpose) <- rownames(t(val_groups[[lm]]))
  colnames(val_transpose) <- colnames(t(val_groups[[lm]]))
  
  val_scale <- scale(val_transpose, center = FALSE, scale = TRUE) # All the model will use scale data
  val_scale[is.na(val_scale)] <- 0
  
  # Support Vector Machine Model:
  
  support_vm <- cfBuild(scale, response_matrix, bootNum = 65, ensNum = 65, cpus = 4) # Creates model
  
  # Validation of the models:
  
  pred_val <- cfPredict(support_vm, val_scale)
  names(pred_val) <- c("prediction", "Coef Score")
  
  pred_val$Responder[pred_val$prediction == "Responder"] <- 
    pred_val$`Coef Score`[pred_val$prediction == "Responder"]/100 
  
  pred_val$Responder[pred_val$prediction == "Non_Responder"] <- 
    1 - (pred_val$`Coef Score`[pred_val$prediction == "Non_Responder"]/100)
  
  pred_val$Non_Responder <- 1 - pred_val$Responder
  
  # MCC:
  
  mcc_val = mcc(preds = pred_val$prediction, actuals = val_response$responses)
  
  # ROC curves: 
  
  roc_val = roc(val_response$responses, pred_val$Non_Responder)
  
  average_right <- data.frame(groups = names(groups)[[lm]], 
                              percentage_accuracy = getAvgAcc(support_vm)$Test, 
                              MCC = mcc_val,
                              AUC = roc_val$auc,
                              stringsAsFactors = FALSE) 
  
  accuracy_table <- rbind(accuracy_table, average_right) # Append %CC, AUC and MCC to the accuracy table. 
  
  # Save model as a R object. 
  
  saveRDS(support_vm, 
          file = paste(output, "1_SVM_",
                       names(groups)[[lm]], ".R", sep = ""),
          ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
  
  }
  
  else { next }
  
}

# OUTPUT: 

# Besides the functional models, the OUTPUT of this script is a table with Accuracy score of the model and the
# validation of the model based on the MCC value. 

accuracy_table[, c(2)] <- round(accuracy_table[, c(2)], digits = 0) 
accuracy_table[, c(3, 4)] <- round(accuracy_table[, c(3, 4)], digits = 2) 

write.table(accuracy_table, 
            file = paste(output, "1_accuracy_table.txt", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE) 


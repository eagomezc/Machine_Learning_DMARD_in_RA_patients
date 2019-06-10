#------------------------------------------ CLASSYFIRE LIPIDOMIC DATA --------------------------------------------#

# Rheumatoid arthritis is a inflammatory disease that characterized for not having a resolution phase in the 
# inflammatory procces. One of the most used treatments is modifying antirheumatic drugs (DMARDs), that is not 
# always efective. 

# We have a dataset of the lipid profiles of rheumatoid arthritis patients that respond or not to the MTX treatment.
# This is the first approach to create a model using machine learning that predicts if a patient will react or not 
# to the treatment with DMARDs.

#---> LIBRARY LOAD:

library(classyfire)

#---> DATA MANIPULATION: 

# TRAINING SET:

# Open the txt file with the profiles information. Make sure that the path is correct:

# The dataset consist in a tab-delimited file in .txt format with the follow specifications: 
# Columns: The different samples (each patient data)
# Row number 1: Class row that contains the word "Responder" if the patient respond to treatment and "Non_Responder"
# otherwhise.
# Following rows: All the lipid mediators used to create the model.

# See a_Toy_Data/1_classyfire_(SVM models)_toy_data.txt

lm_profile <- read.table(
 file = "GitHub/2018_Machine_Learning_MTX_treatment_in_RA_patients/a_Toy_Data/1_classyfire_(SVM models)_toy_data.txt",
 header = TRUE,
 row.names = 1,
 sep = "\t")

# Separates the profiles data by lipid mediators types. 

# by Substrates:

dha <- lm_profiles[2:25, ]
n_three_DPA <- lm_profiles[26:35, ]
epa <- lm_profiles[36:38, ]
aa <- lm_profiles[39:56, ]

#---> DATA PREPARATION:

# Getting the explanatory (x) and response (y) variable. By explanatory, it means all the data that can explain 
# why a patient respond or not to the treatment (the lipid meadiator profiles), and the response variable is if 
# the patients respond or not to the treatment. 

# Transpose columns and rows:
lm_profiles_transpose <- t(lm_profiles)

# Explanatory and Response variable:
response <- data.frame(row.names = row.names(lm_profiles_transpose), # Create the Response variable.
                       responses = lm_profiles_transpose[, 1]) 
response_matrix <- as.matrix(response)

explanatorys <-  lm_profiles_transpose[, -1] # Delete the first row that contains the Response variable.

# Because the response variable was in the data.frame all the elements were saved as factors, and classyfire 
# requieres a numeric matrix. here we make it:

explanatory <- matrix(as.numeric(unlist(explanatorys)),nrow=nrow(explanatorys))
rownames(explanatory) <- rownames(explanatorys)
colnames(explanatory) <- colnames(explanatorys)

# If you are working with machine learning, the best method of scalation is standarization. 
# Scale data prevents that the modelfrom being based on variables with high normal values.
explanatory_scale <- scale(explanatory, center = FALSE, scale = TRUE)

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

# Accuracy table: 

getAvgAcc(support_lmprofiles_scale)$Test # Get the %CC (Overall percentage of correctly classified test objects)
getConfMatr(support_lmprofiles_scale) # Get a table of the consensus classification of the best model. 

# Creates a table with all the models and the %CC.

accuracy_table <- data.frame(groups = "scalated lm profiles",
                             percentage_accuracy = getAvgAcc(support_lmprofiles_scale)$Test,
                             stringsAsFactors = FALSE)

# Save the models as an R object:
# To avoid running all the scripts to obtain the model, it can be saved as R objects. If you want to use it, 
# you can call it using "readRDS". 

# Make sure that the you specify the path were you want to save your model: 

saveRDS(support_lmprofiles_scale, 
        file = "GitHub/2018_Machine_Learning_MTX_treatment_in_RA_patients/c_Expected_Output/1_SVM_lmprofiles_scale.R", 
        ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)

#---> SVM PER GROUP: 

# Run the same analysis but automatically for the rest of the groups by Substrates. 

# Create a list with the names of all the subgroups:

# Modifying the original lipid mediator file, including other data (such as clinical scores) and dividing the table
# in different sections is possible to create automatically all the models you want. YOU NEED TO MAKE SURE WHATSOEVER
# that the "groups" list and the names of the "groups" list is updated. 

groups <- list(dha, n_three_DPA, epa, aa) # Update in case you want to create other models. 

# Create a vector with the names associated to all the elements in the list: 

names(groups) <- c("DHA", "n-3 DPA", "EPA", "AA") # Update in case you want to create other models. 

# The loop goes through all the elements in gropus and creates a model for each of them. The model can not work with
# only one row, so the "if" makes sure that only the subgroups with more than one row are analyzed. 

for (lm in 1:length(groups)) {
  
  if (nrow(groups[[lm]]) > 1) {
  
  transpose <- matrix(as.numeric(t(groups[[lm]])),nrow=nrow(t(groups[[lm]]))) # Transpose and create matrix.
  rownames(transpose) <- rownames(t(groups[[lm]]))
  colnames(transpose) <- colnames(t(groups[[lm]]))
  
  
  scale <- scale(transpose, center = FALSE, scale = TRUE) # All the model will use scale data
  
  support_vm <- cfBuild(scale, response_matrix, bootNum = 65, ensNum = 65, cpus = 4) # Creates model
  
  average_right <- data.frame(groups = names(groups)[[lm]], 
                              percentage_accuracy = getAvgAcc(support_vm)$Test) # Save %CC value.
  
  accuracy_table <- rbind(accuracy_table, average_right) # Append %CC to the accuracy table. 
  
  # Save model as a R object. 
  
  saveRDS(support_vm, 
          file = paste("GitHub/2018_Machine_Learning_MTX_treatment_in_RA_patients/c_Expected_Output/1_SVM_",
                       names(groups)[[lm]], ".R", sep = ""),
          ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
  
  }
  
  else { next }
  
}

# Write table with accuracy values: 

write.table(accuracy_table, 
            file = "GitHub/2018_Machine_Learning_MTX_treatment_in_RA_patients/c_Expected_Output/1_accuracy_table.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE) 


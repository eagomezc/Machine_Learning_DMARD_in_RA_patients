#------------------------------------------ RANDOMFOREST LIPIDOMIC DATA ------------------------------------------#

# Rheumatoid arthritis is a inflammatory disease that characterized for not having a resolution phase in the 
# inflammatory procces. One of the most used treatments is modifying antirheumatic drugs (DMARDs), that is not 
# always efective. 

# We have a dataset of the lipid profiles of rheumatoid arthritis patients that responded or not to the treatment.

# This is the first approach to creat a model using "randomForest" (decision trees).

#---> LIBRARY LOAD:

library(randomForest)
library(ggplot2)
library(caret)
library(pROC)

set.seed(415) # To get same results even with the random part.

#---> INPUT AND OUTPUT:

# In this section please specify where are the input files and where you want to save the output files.
# In the input and output variable you can see what is the path expected from the user to write.

input <- "../a_Data/2_randomForest_(RF_models)/"
output <- "../c_Expected_Output/2_randomForest_(RF_models)/"

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
  file = paste(input, "2_randomForest_(RF_models)_data.txt", sep = ""),
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t")

lm_profiles_scale <- as.data.frame(scale(lm_profiles[, -1], center = FALSE, scale = TRUE))

# If all the values from one column are the same, the scalation will give you NA. For those cases, to avoid errors,
# replace the NA for zeros. 

lm_profiles_scale[is.na(lm_profiles_scale)] <- 0

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

dha <- lm_profiles_scale[ ,c(1:23, 55)]
n_three_DPA <- lm_profiles_scale[ , c(24:33, 55)]
epa <- lm_profiles_scale[ , c(34:36, 55)]
aa <- lm_profiles_scale[ , c(37:54, 55)]

# Best LM:

val_top_four_model <- lm_profiles_scale[ ,c(32, 33, 4, 40)]
val_top_six_model <- lm_profiles_scale[ ,c(32, 33, 4, 40, 18, 43)]

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
val_lm_profiles_scale[is.na(val_lm_profiles_scale)] <- 0
names(val_lm_profiles_scale) <- make.names(names(val_lm_profiles_scale))

# By Substrates:

val_dha <- val_lm_profiles_scale[ ,c(1:23)]
val_n_three_DPA <- val_lm_profiles_scale[ , c(24:33)]
val_epa <- val_lm_profiles_scale[ , c(34:36)]
val_aa <- val_lm_profiles_scale[ , c(37:54)]

# Best LM:

top_four_model <- val_lm_profiles_scale[ ,c(32, 33, 4, 40)]
top_four_model$responses <- val_lm_profiles$responses
top_six_model <- val_lm_profiles_scale[ ,c(32, 33, 4, 40, 18, 43)]
top_six_model$responses <- val_lm_profiles$responses

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

# Get the confusion matrix of the model, sensitivity and specificity: 

confusion_lm_profiles <- as.data.frame(rf_lm_profiles_final$confusion)
confusion_lm_profiles[is.na(confusion_lm_profiles)] <- 0

# Calculates sensitivity, specificity and AUC.

sensitivity_lm_profiles <- confusion_lm_profiles[2, 2]/(confusion_lm_profiles[2, 2] + confusion_lm_profiles[2, 1])
specificity_lm_profiles <- confusion_lm_profiles[1, 1]/(confusion_lm_profiles[1, 1] + confusion_lm_profiles[1, 2])

# Save relevant Plots:

pdf(file = paste(output, "2_RF_lmprofiles_scale.pdf", sep = ""), 
    width = 14, height = 10, onefile = TRUE)
varImpPlot(rf_lm_profiles_final, sort = TRUE, main = "lm_profiles") # Importance of the variable for the model. 
plot(rf_lm_profiles_final, main = "lm_profiles") # Decreasing of the error base on the number of tres. 
legend("topright", inset=.05, title="Curves:",
       c("Non responder", "OOB error", "Responder"), fill=c("red", "black", "green"),
       text.font = 3, cex = 1)
dev.off()

#---> MODEL EVALUATION:

# "prediction" takes the created models and the validation dataset to try to predict which samples belongs to the
# responder and non-responder. 

pred_rf_lm_profiles <- as.data.frame(predict(rf_lm_profiles_final, val_lm_profiles_scale))
names(pred_rf_lm_profiles) <- c("prediction")

# In order to further evaluate the predictivenss of this approach we next calculated  ROC curves (AUC), 
# which represents the accuracy of the model at predicting outcome.

pred_rf_lm_profiles_p <- as.data.frame(predict(rf_lm_profiles_final, val_lm_profiles_scale, type = "prob"))

roc_val = roc(val_lm_profiles$responses, pred_rf_lm_profiles_p$Non_Responder) 

# Confusion Matrix for the validation:

confusion_val <- confusionMatrix(factor(pred_rf_lm_profiles$prediction, levels = c("Responder", "Non_Responder")),
                                 factor(val_lm_profiles$responses, levels = c("Responder", "Non_Responder")))

table_parameters <- as.data.frame(confusion_val$byClass)

# Creates a table with all the models, the %CC and the AUC scores.

# An error estimate is made for the cases which were not used while building the tree. That is called an 
# OOB (Out-of-bag) error estimate which is mentioned as a percentage.

# In this case, then, if a model has 1% error will means that we can predict with 99% accuracy.

no_oob_error_table <- data.frame(groups = "scalated lm profiles",
                                 percentage_accuracy = 100 - ((rf_lm_profiles_final$err.rate[10000])*100),
                                 sensitivity_validation = sensitivity_lm_profiles,
                                 specificity_validation = specificity_lm_profiles,
                                 TP_per = (1 - confusion_lm_profiles[2, 3])*100,
                                 FP_per = confusion_lm_profiles[1, 3]*100,
                                 TN_per = (1 -confusion_lm_profiles[1, 3])*100,
                                 FN_per = confusion_lm_profiles[2, 3]*100,
                                 AUC = roc_val$auc,
                                 sensitivity_evaluation = table_parameters[1,1],
                                 specificity_evaluation = table_parameters[2,1], 
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

groups <- list(dha, n_three_DPA, epa, aa, 
               top_four_model, top_six_model) # Update in case you want to create other models. 
val_groups <- list(val_dha, val_n_three_DPA, val_epa, val_aa, val_top_four_model, val_top_six_model) # Update in case you want to create other models.
actual_group <- list(val_lm_profiles, val_lm_profiles, val_lm_profiles, val_lm_profiles, 
                     lm_profiles, lm_profiles)
# Create a vector with the names associated to all the elements in the list: 

names(groups) <- c("DHA", "n-3 DPA", "EPA", "AA", "RvD4 10S-17S-diHDPA 15R-LXA4 MaR1n-3 DPA", 
                   "RvD4 10S-17S-diHDPA 15R-LXA4 5S12S-diHETE 4-14-diHDHA MaR1n-3 DPA") 

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
    
    # Get the confusion matrix of the model, sensitivity and specificity: 
    
    confusion <- as.data.frame(random_forest_final$confusion)
    confusion[is.na(confusion)] <- 0
    
    # Calculates sensitivity and specificity.
    
    sensitivity <- confusion[2, 2]/(confusion[2, 2] + confusion[2, 1])
    specificity <- confusion[1, 1]/(confusion[1, 1] + confusion[1, 2])
    
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
    
    # Evaluation: 
   
    pred_val <- as.data.frame(predict(random_forest_final, val_groups[[lm]]))
    names(pred_val) <- c("prediction")
    
    actuals <- actual_group[[lm]]
    
    # ROC = 
    
    pred_val_p <- as.data.frame(predict(random_forest_final, val_groups[[lm]], type = "prob"))
    
    roc_val_all = roc(actuals$responses, pred_val_p$Non_Responder) 
    
    roc_curve <- paste("roc_", names(groups)[[lm]], sep = "")
    
    assign(roc_curve, roc_val_all)
    
    # Confusion Matrix: 
    
    confusion_val_all <- confusionMatrix(factor(pred_val$prediction, levels = c("Responder", "Non_Responder")),
                                     factor(actuals$responses, levels = c("Responder", "Non_Responder")))
    
    table_parameters_all <- as.data.frame(confusion_val_all$byClass)
    
    # Accuracy table: 
    
    accuracy_table <- data.frame(groups = names(groups)[[lm]],
                                 percentage_accuracy = 100 - ((random_forest_final$err.rate[10000])*100),
                                 sensitivity_validation = sensitivity,
                                 specificity_validation = specificity,
                                 TP_per = (1 - confusion[2, 3])*100,
                                 FP_per = confusion[1, 3]*100,
                                 TN_per = (1 -confusion[1, 3])*100,
                                 FN_per = confusion[2, 3]*100,
                                 AUC = roc_val_all$auc,
                                 sensitivity_evaluation = table_parameters_all[1,1],
                                 specificity_evaluation = table_parameters_all[2,1],
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
    
# ---> OUTPUT: 

# SCORE TABLE:

# Besides the functional models, the OUTPUT of this script is a table with Accuracy score of the model and the
# validation of the model based on the MCC value. 

# The resulting plots also gave us an idea of the performance of the models, that includes the "importance" analysis,
# and the improvement of the model based on the number of tress used to create it.

no_oob_error_table[, c(-1, -3, -4, -9, -10, -11)] <- round(no_oob_error_table[, c(-1, -3, -4, -9, -10, -11)], digits = 0)
no_oob_error_table[, c(3, 4, 9, 10, 11)] <- round(no_oob_error_table[, c(3, 4, 9, 10, 11)], digits = 2)

write.table(no_oob_error_table, 
            file = paste(output, "2_accuracy_table_RF_final.txt", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)  

# IMPORTANCE FIGURE: 

# Another OUTPUT is a Importance Plot showing the lipid mediators from the lipid mediator model that contributes
# the most in the prediction step:

# Create an importance list using the importance function from randomForest.

importance_list <- as.data.frame(importance(rf_lm_profiles_final, type = 1))
importance_list$lipid_mediators <- rownames(importance_list)

# Modify lipid mediators names to make them ready for publishing: 

# Everything that starts with a number will have removed the X and it will be in '' to avoid been wrongfully called
# in the figure: 

importance_list$lipid_mediators[grep("^X", 
             importance_list$lipid_mediators)] <- paste(gsub("^X", "'", importance_list$lipid_mediators[grep("^X", 
                                                              importance_list$lipid_mediators)]), "'", sep = "")

# Special lipid mediator cases that cann't be change automatically: 

importance_list$lipid_mediators[importance_list$lipid_mediators == "RvD4"] <- "'RVD4'"
importance_list$lipid_mediators[importance_list$lipid_mediators == "RvT4"] <- "'RVT4'"
importance_list$lipid_mediators[importance_list$lipid_mediators == "PGF2a"] <- "PGF[2][alpha]"
importance_list$lipid_mediators[importance_list$lipid_mediators == "TXB2"] <- "TXB[2]"
importance_list$lipid_mediators[importance_list$lipid_mediators == "Maresin1"] <- "MaR1"
importance_list$lipid_mediators[importance_list$lipid_mediators == "Maresin2"] <- "MaR2"
importance_list$lipid_mediators[importance_list$lipid_mediators == "'15R.LXA4'"] <- "'15R-LXA'[4]"

# Everything that finish with a 4 will be finish with a subscript 4 (in parse language []): 

importance_list$lipid_mediators <- gsub("4$", "[4]", importance_list$lipid_mediators)

# Everything that finish with a B4 and starts with a number will be replace with the subscript 4: 

importance_list$lipid_mediators <- gsub("B4'", "B'[4]", importance_list$lipid_mediators)

# Replace "." with "-" when requieres: 

importance_list$lipid_mediators <- gsub("\\.", "-", importance_list$lipid_mediators)

# Transform n-3 DPA as a subscript: 

importance_list$lipid_mediators <- gsub("n.3.DPA", "[n-3~DPA]", importance_list$lipid_mediators)

# Add the information about the LM substrates from each lipid mediator to be able to add color after:

importance_list$fatty_acids <- c(1:54)
importance_list$fatty_acids[1:23] <- "DHA"      
importance_list$fatty_acids[24:33] <- "n-3 DPA" 
importance_list$fatty_acids[34:36] <- "EPA" 
importance_list$fatty_acids[37:54] <- "AA" 

# Order the data in increasing order based on the Mean Decrease Accuracy and create a factor column
# (this with the purpose of get the LM in decreasing order in the figure):

importance_list <- importance_list[order(importance_list$MeanDecreaseAccuracy), ]

importance_list$lipid_mediators <- factor(importance_list$lipid_mediators, 
                                          levels = importance_list$lipid_mediators[
                                            order(importance_list$MeanDecreaseAccuracy)])

# Create a vector with the parse names of the lipid mediators son they can be publishing printed in 
# the figure: 

lipids <- as.character(importance_list$lipid_mediators)

# As in heatmaps. Add the neccesary information to add colors to the x and y axis.

lm_classes <- as.factor(importance_list$fatty_acids)
names(lm_classes) <- importance_list$lipid_mediators

dha_index<-which((lm_classes=="DHA")== TRUE)
n_three_index<-which((lm_classes=="n-3 DPA")== TRUE)
epa_index<-which((lm_classes=="EPA")== TRUE)
aa_index<-which((lm_classes=="AA")== TRUE)

lm_colors <- NULL
lm_colors[dha_index] <- "blue"
lm_colors[n_three_index] <- "brown"
lm_colors[epa_index] <- "darkgoldenrod1"
lm_colors[aa_index] <- "darkslategray"

# Save the plot as an object so we can save the pdf and png files: 

importance_plot <- ggplot(data = importance_list, mapping = aes(x = lipid_mediators, y = MeanDecreaseAccuracy, color = fatty_acids)) +
  geom_point(size = 8) +
  scale_y_continuous(name = "Mean Decrease Accuracy") +
  labs(x = "Lipid Mediators") +
  scale_x_discrete(labels = parse(text = lipids)) + 
  scale_color_manual(name = "Lipid Mediators Metabolomes:", values = c("darkslategray", "blue", "darkgoldenrod1", 
                                                                       "brown")) +
  coord_flip()  +
  theme(axis.title = element_text(size = 40),
        axis.text.x  = element_text(size = 40, hjust = 0.5),
        axis.text.y  = element_text(size = 40, hjust = 1, colour = lm_colors),
        legend.title = element_text(size = 40),
        legend.text  = element_text(size = 35), 
        legend.position = "top",
        aspect.ratio = 2/1,
        legend.key = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black")) 

# Export the plot:

pdf(file = paste(output, "importance_figure_final.pdf", sep = ""),
    width = 25, height = 35, onefile = TRUE)

importance_plot

dev.off()


# ACCURACY FIGURE: 

# Barplot with LM models information: 

# Creates a factor column with the correct names of the LM models: 

no_oob_error_table$models <- factor(c("All LM", "DHA", "n-3 DPA", "EPA", "AA", "RvD4, 10S,17S-diHDPA, \n 15R-LXA4, MaR1n-3 DPA", 
                                      "RvD4, 10S,17S-diHDPA, 15R-LXA4, \n 5S,12S-diHETE, 4,14-diHDHA, MaR1n-3 DPA"),
                                    levels = c("All LM", "DHA", "n-3 DPA", "EPA", "AA", "RvD4, 10S,17S-diHDPA, \n 15R-LXA4, MaR1n-3 DPA", 
                                               "RvD4, 10S,17S-diHDPA, 15R-LXA4, \n 5S,12S-diHETE, 4,14-diHDHA, MaR1n-3 DPA"))

# Creates the figure: 

accuracy_plot <- ggplot(data = no_oob_error_table, aes(x = models, y = percentage_accuracy)) +
  geom_bar(stat = "identity",  position = position_dodge(w = 0.5), 
           fill = c(rep("goldenrod1", 5), rep("hotpink", 2)), colour = "black", size = 2) + # Two colors and border
  geom_text(position = position_dodge(w = 0.5), vjust = -0.7,
            aes(label = paste(percentage_accuracy, "%", sep = "")), size = 15) +
  scale_y_continuous(name = "% Accuracy Score", breaks = seq(from = 0, to = 100, by = 20), 
                     expand = c(0, 1)) + # Expand to delete the space between zero and the x-axis. 
  coord_cartesian(ylim = c(1, 100)) +
  theme(axis.title = element_text(size = 40),
        axis.title.x = element_blank(),
        axis.text.x  =  element_text(size = 40, hjust = 1, angle = 45, colour = "black"), # Put color to the labels
        axis.text.y  = element_text(size = 40, hjust = 1, colour = "black"), # Put color to the labels
        axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
        axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(0.4, "cm"))  # Size of every axis sep.  

# Export the plot:

pdf(file = paste(output, "accuracy_figure_final.pdf", sep = ""),
    width = 20, height = 20, onefile = TRUE)

accuracy_plot

dev.off()

# CONFUSIO PLOT:

confusio_table <- data.frame(InitClass = c(rep("Responder", times = 1), rep("Non Responder", times = 2), 
                                                 rep("Responder", times = 1)),
                                   PredClass = c(rep("Responder", times = 2), rep("Non Responder", times = 2)),
                                   Percentage = c(no_oob_error_table[no_oob_error_table$groups == "n-3 DPA", ]$TP_per, 
                                                  no_oob_error_table[no_oob_error_table$groups == "n-3 DPA", ]$FP_per,
                                                  no_oob_error_table[no_oob_error_table$groups == "n-3 DPA", ]$TN_per,
                                                  no_oob_error_table[no_oob_error_table$groups == "n-3 DPA", ]$FN_per))

# Make the plot:

confusio_plot <- ggplot(data = confusio_table, mapping = aes(x = InitClass, y = Percentage, fill = PredClass)) +
  geom_bar(stat = "identity", colour = "black", size = 2) +
  geom_text(data = confusio_table, aes(x = InitClass, y = Percentage, label = paste(Percentage,"%",sep="")), 
            size = 12, position = "stack", colour = c("chartreuse3", "white", "black", "blue3"), 
            alpha = c(0, 1, 1, 0), vjust = 1.2) +
  scale_fill_manual(values = alpha(c("blue3", "chartreuse3"), 0.7), labels = c("Non Resp", "Resp")) +
  scale_y_continuous(name = "Class Predictions (%)", breaks = seq(from = 0, to = 100, by = 20),
                     expand = c(0, 1)) +
  scale_x_discrete(labels = c("Non Resp", "Resp")) + 
  labs(x = "Classes") +
  coord_cartesian(ylim = c(1, 100)) +
  theme(plot.title = element_text(size = 40, hjust = 0.5),
        axis.title = element_text(size = 40),
        axis.text.x  = element_text(size = 38, hjust = 0.5, colour = "black"),
        axis.text.y  = element_text(size = 40, hjust = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.key.size = unit(1.3, "cm"), 
        legend.text  = element_text(size = 40),
        axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
        axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(0.4, "cm"))

pdf(file = paste(output, "confusio_figure_n-3DPA.pdf", sep = ""),
    width = 11, height = 12, onefile = TRUE)

confusio_plot

dev.off()

# ROC Figures: 

roc_best_lm <- ggroc(list(`roc_RvD4 10S-17S-diHDPA 15R-LXA4 MaR1n-3 DPA`,
                          `roc_RvD4 10S-17S-diHDPA 15R-LXA4 5S12S-diHETE 4-14-diHDHA MaR1n-3 DPA`), legacy.axes = TRUE, size = 3) + 
  scale_color_manual(labels = c(paste("RvD4, 10S,17S-diHDPA, \n 15R-LXA4, MaR1n-3 DPA (AUC=", 
                                      round(`roc_RvD4 10S-17S-diHDPA 15R-LXA4 MaR1n-3 DPA`$auc, 2), ")", sep = ""),
                                paste("RvD4, 10S,17S-diHDPA, 15R-LXA4, 5S,12S-diHETE, \n 4,14-diHDHA, MaR1n-3 DPA (AUC=", 
                                      round(`roc_RvD4 10S-17S-diHDPA 15R-LXA4 5S12S-diHETE 4-14-diHDHA MaR1n-3 DPA`$auc, 2), 
                                      ")", sep = "")), 
                     values = c("chocolate1", "blue4")) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed", size = 3) +
  theme(axis.title.y = element_text(size = 75, colour = "black", margin = unit(c(0, 1, 0, 1), "cm")),
        axis.title.x = element_text(size = 75, colour = "black", margin = unit(c(1, 0, 1, 0), "cm")),
        axis.text.x  = element_text(size = 60, colour = "black"), 
        axis.text.y  = element_text(size = 60, hjust = 1, colour = "black"), # Put color to the labels
        axis.line = element_line(colour = 'black', size = 3), # Color and thickness of axis
        axis.ticks = element_line(colour = "black", size = 3), # Color and thickness of every axis sep. 
        panel.background = element_rect(fill = "white"),
        axis.ticks.length = unit(0.4, "cm"),
        legend.title = element_blank(),
        legend.text  = element_text(size = 40), 
        legend.position = c(0.63, 0.15),
        legend.background = element_rect(colour = "black", size = 3),
        legend.key = element_rect(fill = "white"),
        legend.key.width = unit(2, "cm"))


pdf(file = paste(output, "roc_best_lm_final.pdf", sep = ""),
    width = 25, height = 20, onefile = TRUE)

roc_best_lm

dev.off()


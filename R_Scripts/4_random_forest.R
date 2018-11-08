#------------------------------------- randomForest MODEL LIPIDOMIC DATA  ---------------------------------------------#

# Rheumatoid arthritis is a inflammatory disease that characterized for not having a resolution phase in the inflammatory
# procces. One of the most used treatments is Methotrexate (MTX), that is not always efective. 

# We have a dataset of the lipid profiles of rheumatoid arthritis patients that responded or not to the treatment.

# This is the first approach to creat a model using "randomForest" (decision trees).

# From Seth methodology: 
# ntree = 100000
# 10 test were performed 
# Out of the bag error (OOB) taken as the predective accuracy.
# Number of ensembles = 50 (we are using the same as before)

#---> LIBRARY LOAD:

library(randomForest)
set.seed(415) # To get same results even with the random part.

#---> DATA MANIPULATION: 

# TRAINING SET:

# Open the txt file with the profiles information. 

lm_profiles <- read.table(
  file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/4_random_forest/rheumatoid_arthritis_LM_profile.txt",
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t")

lm_profiles <- lm_profiles[-c(16), ] # All the values of Marasein 1 are equal to zero, so it does not add info. 
lm_profiles <- lm_profiles[ ,-c(41, 42)] # The samples PL041 and PL042 were eliminated because there was a problem with
                                        # their extraction.

# Normalization of the data (scalation). The standarization has to be made with the transpose data frame (IMP).

lm_profiles_scale <- scale(t(lm_profiles), center = FALSE, scale = TRUE)

# Add the classification variable to the data frame (Responder and non responder):

# Getting the explanatory (x) and response (y) variable. By explanatory, it means all the data that can explain why a
# patient response or not to the treatment (the lipid meadiator profiles) and the response variable is if the patients
# response or not to the treatment. In random Forest you have to create a formula where the response variable is
# explain in terms of the explanatory variable (responses ~ everything else).

response <- data.frame(row.names = colnames(lm_profiles),
                         responses = c(rep("responder", 30), rep("non-responder", 22)))

lm_profiles_scale <- cbind(lm_profiles_scale, response)
names(lm_profiles_scale) <- make.names(names(lm_profiles_scale))

# Separates the profiles data by lipid mediators types. 

# Fatty acids origin:

dha <- lm_profiles_scale[ ,c(1:23, 55)]
n_three_DPA <- lm_profiles_scale[ , c(24:33, 55)]
aa <- lm_profiles_scale[ , c(37:54, 55)]
epa <- lm_profiles_scale[ , c(34:36, 55)]

# Lipid Mediators groups:

resolvins_d <- lm_profiles_scale[, c(1:8, 55)]
protectins <- lm_profiles_scale[ , c(9:12, 55)]
pctr <- lm_profiles_scale[ , c(13:15, 55)]
maresins <- lm_profiles_scale[ , c(16:20, 55)]
mctr <- lm_profiles_scale[ , c(21:23, 55)]
rvt <- lm_profiles_scale[, c(24:27, 55)]
rvd <- lm_profiles_scale[ , c(28:30, 55)]
pd <- lm_profiles_scale[ , c(31:32, 55)]
mar <- lm_profiles_scale[ , c(33, 55)]
lx <- lm_profiles_scale[ , c(37:41, 55)]
ltb <- lm_profiles_scale[ , c(42:47, 55)]
lt <- lm_profiles_scale[ , c(48:50, 55)]
pg <- lm_profiles_scale[ , c(51:53, 55)]
tx <- lm_profiles_scale[, c(54, 55)]

# Combination: 

# Seth did this one, so I did it as well to validate it. 

dha_n_dpa_epa <- lm_profiles_scale[ , c(1:33, 34:36, 55)]

#---> MACHINE LEARNING (randomForest R):

# In Random Forests the idea is to decorrelate the several trees which are generated on the different bootstrapped 
# samples from training Data.And then we simply reduce the Variance in the Trees by averaging them.

# Averaging the Trees helps us to reduce the variance and also improve the Perfomance of Decision Trees on Test 
# Set and eventually avoid Overfitting.

# The idea is to build lots of Trees in such a way to make the Correlation between the Trees smaller.

# BEST MTRY:
# mtry is number of variables available for splitting at each tree node. Random Forest creates several trees, each
# one using different variables to create the best version of it. With mtry we can define how many variables the 
# data is split into to create the different trees.
# More: https://stats.stackexchange.com/questions/102867/random-forest-mtry-question

# In this case we defined the 54 variables to create a loop to define which mtry is the best one for our model. 

oob_error <- double(ncol(lm_profiles_scale) - 1) #Define number of variable. -1 is because the last column is classes

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

# Save graphs:

pdf(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/4_random_forest/graphs_RF/RF_lmprofiles_scale.pdf",
    width = 14, height = 10, onefile = TRUE)
varImpPlot(rf_lm_profiles_final, sort = TRUE, main = "lm_profiles") # Importance of the variable for the model. 
plot(rf_lm_profiles_final, main = "lm_profiles") # Decreasing of the error base on the number of tres. 
legend("topright", inset=.05, title="Curves:",
       c("Non responder", "OOB error", "Responder"), fill=c("red", "black", "green"),
       text.font = 3, cex = 1)
dev.off()

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

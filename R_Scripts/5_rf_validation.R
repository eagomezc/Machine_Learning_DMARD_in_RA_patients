#-------------------------------- randomForest VALIDATION AND ROC CURVES  ---------------------------------------------#

# Rheumatoid arthritis is a inflammatory disease that characterized for not having a resolution phase in the inflammatory
# procces. One of the most used treatments is Methotrexate (MTX), that is not always efective. 

# We have a dataset of the lipid profiles of rheumatoid arthritis patients that responded or not to the treatment.

# This is the first approach to test the models created using random Forest and a validation dataset. 

#---> LIBRARY LOAD:

library(randomForest)
set.seed(415) # To get same results even with the random part.

# VALIDATION SET: 

# From the previous analysis, we saw that the combination DHA-N3DPA-EPA (84.61), DHA (82.69) and EPA (E-Series Resolvin)
# (78.84) were the best models. The whole LM profiles (78.84) model will be test as well. 

# Open the txt file with the profiles information. 

val_lm_profiles <- read.table(
  file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/5_rf_validation/rheumatoid_arthritis_LM_profile.txt",
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t")

# Separates the profiles data by lipid mediators types. 

# Fatty acids origin:

val_dha <- val_lm_profiles[1:23, ]
val_epa <- val_lm_profiles[34:36, ]

# Combination:

val_dha_n_dpa_epa <- val_lm_profiles[c(1:33, 34:36), ]

# MODELS:

# Open the already saved models using the function readRDS. 

rf_dha_dpa_epa_scale <- readRDS(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/5_rf_validation/RF_DHA_N-3DPA_EPA.R", 
                                 refhook = NULL)
rf_dha_scale <- readRDS(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/5_rf_validation/RF_DHA.R", 
                         refhook = NULL)
rf_EPA_scale <- readRDS(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/5_rf_validation/RF_EPA.R", 
                         refhook = NULL)
rf_lmprof_scale <- readRDS(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/5_rf_validation/RF_lmprofiles_scale.R", 
                            refhook = NULL)

# DATA PREPARATION:

# Since all the models were created using scalation and transpose data, the validation datasets has to be scaled and 
# transpose as well. 

val_dha_scale <- as.data.frame(scale(t(val_dha), center = FALSE, scale = TRUE))
val_dha_scale[, 13] <- seq(0, by = 54)
val_dha_scale[, 21] <- seq(0, by = 54)
names(val_dha_scale) <- make.names(names(val_dha_scale))

val_epa_scale <- as.data.frame(scale(t(val_epa), center = FALSE, scale = TRUE)) 
names(val_epa_scale) <- make.names(names(val_epa_scale))

# When there is a row with the same number, the scalation does not work. In this case all the values comes back to 0,
# the original value. 

val_dha_dpa_epa_scale <- as.data.frame(scale(t(val_dha_n_dpa_epa), center = FALSE, scale = TRUE))
val_dha_dpa_epa_scale[, 13] <- seq(0, by = 54)
val_dha_dpa_epa_scale[, 21] <- seq(0, by = 54)
names(val_dha_dpa_epa_scale) <- make.names(names(val_dha_dpa_epa_scale))

val_lmprofile_scale <- as.data.frame(scale(t(val_lm_profiles), center = FALSE, scale = TRUE))
val_lmprofile_scale[, 13] <- seq(0, by = 54)
val_lmprofile_scale[, 21] <- seq(0, by = 54)
names(val_lmprofile_scale) <- make.names(names(val_lmprofile_scale))

#---> VALIDATION TEST:

# "prediction" takes the created models and the validation dataset to try to predict which samples belongs to the
# responder and non-responder. 

pred_dha_n_dpa_epa <- as.data.frame(predict(rf_dha_dpa_epa_scale, val_dha_dpa_epa_scale))
pred_dha <- as.data.frame(predict(rf_dha_scale, val_dha_scale))
pred_epa <- as.data.frame(predict(rf_EPA_scale, val_epa_scale))
pred_lmprofile <- as.data.frame(predict(rf_lmprof_scale, val_lmprofile_scale))

#---> SAVE TABLES:

write.table(pred_dha_n_dpa_epa, 
            file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/5_rf_validation/rf_pred_dha_n_dpa_epa_training.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE) 

write.table(pred_dha, 
            file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/5_rf_validation/rf_pred_dha_training.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE) 

write.table(pred_epa, 
            file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/5_rf_validation/rf_pred_epa_training.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE) 

write.table(pred_lmprofile, 
            file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/5_rf_validation/rf_pred_lmprofile_training.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE) 
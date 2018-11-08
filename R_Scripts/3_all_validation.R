#---------------------------------------------------- VALIDATION  TEST ---------------------------------------------#

# Rheumatoid arthritis is a inflammatory disease that characterized for not having a resolution phase in the inflammatory
# procces. One of the most used treatments is Methotrexate (MTX), that is not always efective. 

# We have a dataset of the lipid profiles of rheumatoid arthritis patients that responded or not to the treatment.

# This is the first approach to test the models created by classyfire using a validation dataset.

#---> LIBRARY LOAD:

library(classyfire)

#---> DATA MANIPULATION: 

# VALIDATION SET: 

# Open the txt file with the profiles information. 

val_lm_profiles <- read.table(
  file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/3_validation_and_permutation_tests/data_validation.txt",
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t")

# Separates the profiles data by lipid mediators types. 

# Fatty acids origin:

val_dha <- val_lm_profiles[1:23, ]
val_n_three_DPA <- val_lm_profiles[24:33, ]
val_aa <- val_lm_profiles[37:54, ]
val_epa <- val_lm_profiles[34:36, ]

# Lipid Mediators groups:

val_maresins <- val_lm_profiles[16:20, ]
val_rvt <- val_lm_profiles[24:27, ]
val_resolvins_d <- val_lm_profiles[1:8, ]
val_protectins <- val_lm_profiles[9:12, ]
val_pctr <- val_lm_profiles[13:15, ]
val_mctr <- val_lm_profiles[21:23, ]
val_rvd <- val_lm_profiles[28:30, ]
val_pd <- val_lm_profiles[31:32, ]
val_lx <- val_lm_profiles[37:41, ]
val_ltb <- val_lm_profiles[42:47, ]
val_lt <- val_lm_profiles[48:50, ]
val_pg <- val_lm_profiles[51:53, ]

# Combination:

val_dha_n_dpa_epa <- val_lm_profiles[c(1:33, 34:36), ]

#---> LOOP FOR ALL THE MODELS:
# Make sure that the names of the lipid mediators are the same as the files in the folder. Add or delete the 
# lm that you are not interested. For examples, here, we dont have the lm profile since the analysis is already
# made. 

groups <- list(val_dha, val_n_three_DPA, val_epa, val_resolvins_d, val_protectins, val_pctr, 
               val_mctr, val_rvd, val_pd, val_lx, val_ltb, val_lt, val_pg)

names(groups) <- c("DHA", "n-3 DPA", "EPA", "Resolvins_d", "Protectins", "PCTRS", 
                   "MCTR", "D-Series Resolvins", "PD", "LX", "LTB", "Leukotrienes", "Prostaglandins")

for (lm in 1:length(groups)) {

# MODEL:

# Open the already saved model using the function readRDS. 

model <- readRDS(paste(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/3_validation_and_permutation_tests/SVM_",
                       names(groups)[[lm]], ".R", sep = ""), refhook = NULL)

# DATA PREPARATION:

# Since all the models were created using scalation and transpose data, the validation datasets has to be scaled and 
# transpose as well. 

validation_scale <- scale(t(groups[[lm]]), center = FALSE, scale = TRUE) 

# When there is a row with the same number, the scalation does not work. In this case all the values comes back to 0,
# the original value. 

validation_scale[is.na(validation_scale)] <- 0

# VALIDATION TEST:

# "cfPredict" takes the created models and the validation dataset to try to predict which samples belongs to the
# responder and non-responder. It creates data frames with the identifications and % of accuracy. 

prediction_validation <- cfPredict(model, validation_scale)

# SAVE TABLES: 

write.table(prediction_validation, 
            file = paste("C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/3_validation_and_permutation_tests/pred_",
                         names(groups)[[lm]], ".txt", sep = ""),             
            sep = "\t",
            quote = FALSE,
            row.names = TRUE) 

}
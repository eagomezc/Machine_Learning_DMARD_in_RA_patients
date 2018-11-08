#-------------------------------- randomForest VALIDATION   -------------------------------------------------------#

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

# Create a list with the names of all the subgroups:
groups <- list(val_n_three_DPA, val_aa, val_resolvins_d, val_protectins, val_pctr, val_maresins, val_mctr, val_rvt, 
               val_rvd, val_pd, val_lx, val_ltb, val_lt, val_pg)

# Create a vector with the names associated to all the elements in the list: 
names(groups) <- c("n-3 DPA", "AA", "Resolvins_d", "Protectins", "PCTRS", "Maresins", "MCTR", "RVT", 
                   "D-Series Resolvins", "PD", "LX", "LTB", "Leukotrienes", "Prostaglandins")

for (lm in 1:length(groups)) {

# MODELS:

# Open the already saved models using the function readRDS. 
  
model <- readRDS(paste(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/5_rf_validation/RF_",
                         names(groups)[[lm]], ".R", sep = ""), refhook = NULL)

# DATA PREPARATION:

# Since all the models were created using scalation and transpose data, the validation datasets has to be scaled and 
# transpose as well. 

validation_scale <- as.data.frame(scale(t(groups[[lm]]), center = FALSE, scale = TRUE)) 
validation_scale[is.na(validation_scale)] <- 0
names(validation_scale) <- make.names(names(validation_scale))


# VALIDATION TEST:

# "prediction" takes the created models and the validation dataset to try to predict which samples belongs to the
# responder and non-responder. 

prediction_validation <- as.data.frame(predict(model, validation_scale))

# SAVE TABLES:

write.table(prediction_validation, 
            file = paste("C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/5_rf_validation/rf_pred_",
                         names(groups)[[lm]], ".txt", sep = ""),  
            sep = "\t",
            quote = FALSE,
            row.names = TRUE) 

}
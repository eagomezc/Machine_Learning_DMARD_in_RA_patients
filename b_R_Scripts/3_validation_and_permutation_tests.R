#------------------------------------- VALIDATION AND PERMUTATION TEST ---------------------------------------------#

# Rheumatoid arthritis is a inflammatory disease that characterized for not having a resolution phase in the inflammatory
# procces. One of the most used treatments is Methotrexate (MTX), that is not always efective. 

# We have a dataset of the lipid profiles of rheumatoid arthritis patients that responded or not to the treatment.

# This is the first approach to test the models created by classyfire using a validation dataset.

#---> LIBRARY LOAD:

library(classyfire)

#---> DATA MANIPULATION: 

# VALIDATION SET: 

# From the previous analysis, we saw that E-Series Resolvin (71.58), Maresins (70.68) and RVT (68.24) were the 
# best models.We are to test this models using a validation sets. The whole LM profiles (63.28) and the
# DPA-N3DPA-EPA (65.43) models will be test as well. 

# Open the txt file with the profiles information. 

val_lm_profiles <- read.table(
  file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/3_validation_and_permutation_tests/data_validation.txt",
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t")

# Separates the profiles data by lipid mediators types. 

# Fatty acids origin:

val_epa <- val_lm_profiles[34:36, ]

# Lipid Mediators groups:

val_maresins <- val_lm_profiles[16:20, ]
val_rvt <- val_lm_profiles[24:27, ]

# Combination:

val_dha_n_dpa_epa <- val_lm_profiles[c(1:33, 34:36), ]

# MODELS:

# Open the already saved models using the function readRDS. 

spv_EPA_scale <- readRDS(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/3_validation_and_permutation_tests/SVM_EPA.R", 
                         refhook = NULL)
spv_mar_scale <- readRDS(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/3_validation_and_permutation_tests/SVM_Maresins.R", 
                         refhook = NULL)
spv_RVT_scale <- readRDS(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/3_validation_and_permutation_tests/SVM_RVT.R", 
                         refhook = NULL)
spv_lmprof_scale <- readRDS(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/3_validation_and_permutation_tests/support_lmprofiles_scale.R", 
                            refhook = NULL)
spv_dha_dpa_epa_scale <- readRDS(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/3_validation_and_permutation_tests/SVM_DHA_N-3DPA_EPA.R", 
                                 refhook = NULL)


# DATA PREPARATION:

# Since all the models were created using scalation and transpose data, the validation datasets has to be scaled and 
# transpose as well. 

val_epa_scale <- scale(t(val_epa), center = FALSE, scale = TRUE) 

val_maresins_scale <- scale(t(val_maresins), center = FALSE, scale = TRUE)

val_rvt_scale <- scale(t(val_rvt), center = FALSE, scale = TRUE)

# When there is a row with the same number, the scalation does not work. In this case all the values comes back to 0,
# the original value. 

val_lmprofile_scale <- scale(t(val_lm_profiles), center = FALSE, scale = TRUE)
val_lmprofile_scale[, 13] <- seq(0, by = 54)
val_lmprofile_scale[, 21] <- seq(0, by = 54)

val_dha_dpa_epa_scale <- scale(t(val_dha_n_dpa_epa), center = FALSE, scale = TRUE)
val_dha_dpa_epa_scale[, 13] <- seq(0, by = 54)
val_dha_dpa_epa_scale[, 21] <- seq(0, by = 54)

#---> VALIDATION TEST:

# "cfPredict" takes the created models and the validation dataset to try to predict which samples belongs to the
# responder and non-responder. It creates data frames with the identifications and % of accuracy. 

pred_epa <- cfPredict(spv_EPA_scale, val_epa_scale)
pred_maresins <- cfPredict(spv_mar_scale, val_maresins_scale)
pred_rvt <- cfPredict(spv_RVT_scale, val_rvt_scale)
pred_lmprofile <- cfPredict(spv_lmprof_scale, val_lmprofile_scale)
pred_dha_dpa_epa <- cfPredict(spv_dha_dpa_epa_scale, val_dha_dpa_epa_scale)

#---> PERMUTATION TEST: 

# MATRICES:

# Open the matrices used to create the models: 

mx_scale_EPA <- as.matrix(read.table(
  file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/3_validation_and_permutation_tests/mx_scale_EPA.txt",
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t"))
mx_scale_Maresins <- as.matrix(read.table(
  file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/3_validation_and_permutation_tests/mx_scale_Maresins.txt",
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t"))
mx_scale_RVT <- as.matrix(read.table(
  file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/3_validation_and_permutation_tests/mx_scale_RVT.txt",
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t"))
mx_scale_lmprof <- as.matrix(read.table(
  file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/3_validation_and_permutation_tests/mx_scale_lmprof.txt",
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t"))
mx_scale_DHA_DPA_EPA <- as.matrix(read.table(
  file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/3_validation_and_permutation_tests/mx_scale_DHA_N-3DPA_EPA.txt",
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t"))

# RESPONSE DATA: 

response_matrix <- as.matrix(read.table(
  file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/3_validation_and_permutation_tests/response_data.txt",
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t"))

# TEST: 

# The model building and testing process is reapeated hundreds of times in an attempt to map objects to randomly 
# permuted classes - a model performance that does not differ substantially form that achieved for the random 
# permutations cannot be considered significant.

# The performances of the models generated for random permutations are compared with the SVM model using a t-test
# to see if they differ. 

permu_test_epa <- cfPermute(mx_scale_EPA, response_matrix, bootNum = 70, ensNum = 70, permNum=70, cpus = 4)
permu_test_mar <- cfPermute(mx_scale_Maresins, response_matrix, bootNum = 70, ensNum = 70, permNum=70, cpus = 4)
permu_test_rtv <- cfPermute(mx_scale_RVT, response_matrix, bootNum = 70, ensNum = 70, permNum=70, cpus = 4)
permu_test_lmprof <- cfPermute(mx_scale_lmprof, response_matrix, bootNum = 70, ensNum = 70, permNum=70, cpus = 4)
permu_test_three <- cfPermute(mx_scale_DHA_DPA_EPA, response_matrix, bootNum = 70, ensNum = 70, permNum=70, cpus = 4)

#Support Graphs of the permutation models:

pdf(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/3_validation_and_permutation_tests/permutation_graphs.pdf",
    width = 14, height = 10, onefile = TRUE)
ggPermHist(permu_test_epa, density = TRUE, percentiles = TRUE, mean = TRUE, median = TRUE)
ggPermHist(permu_test_mar, density = TRUE, percentiles = TRUE, mean = TRUE, median = TRUE)
ggPermHist(permu_test_rtv, density = TRUE, percentiles = TRUE, mean = TRUE, median = TRUE)
ggPermHist(permu_test_lmprof, density = TRUE, percentiles = TRUE, mean = TRUE, median = TRUE)
ggPermHist(permu_test_three, density = TRUE, percentiles = TRUE, mean = TRUE, median = TRUE)
dev.off()

# t-test to compare the average results of Overall percentage of correctly classified test objects of our model and
# the permutation models

t_test_epa <- t.test(spv_EPA_scale$testAcc,permu_test_epa$avgAcc)
t_test_mar <- t.test(spv_mar_scale$testAcc,permu_test_mar$avgAcc)
t_test_rtv <- t.test(spv_RVT_scale$testAcc,permu_test_rtv$avgAcc)
t_test_lmprof <- t.test(spv_lmprof_scale$testAcc,permu_test_lmprof$avgAcc)
t_test_three <- t.test(spv_dha_dpa_epa_scale$testAcc,permu_test_three$avgAcc)

t_tests_results <- data.frame(models = c("EPA", "Maresins", "RTV", "LM_profile", "DHA, N-3DPA, EPA"),
                              t = c(t_test_epa$statistic, t_test_mar$statistic, t_test_rtv$statistic,
                                    t_test_lmprof$statistic, t_test_three$statistic),
                              df = c(t_test_epa$parameter, t_test_mar$parameter, t_test_rtv$parameter,
                                     t_test_lmprof$parameter, t_test_three$parameter), 
                              p_value = c(t_test_epa$p.value, t_test_mar$p.value, t_test_rtv$p.value,
                                          t_test_lmprof$p.value, t_test_three$p.value),
                              averages = c(t_test_epa$estimate, t_test_mar$estimate, t_test_rtv$estimate,
                                           t_test_lmprof$estimate, t_test_three$estimate))

#---> SAVE TABLES: 

write.table(t_tests_results, 
            file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/3_validation_and_permutation_tests/t_test_results.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE) 

write.table(pred_epa, 
            file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/3_validation_and_permutation_tests/pred_epa.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE) 

write.table(pred_maresins, 
            file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/3_validation_and_permutation_tests/pred_maresins.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE) 

write.table(pred_lmprofile, 
            file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/3_validation_and_permutation_tests/pred_lmprofile.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE) 

write.table(pred_dha_dpa_epa, 
            file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/3_validation_and_permutation_tests/pred_dha_dpa_epa.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE) 

write.table(pred_rvt, 
            file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/3_validation_and_permutation_tests/pred_rvt.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE) 

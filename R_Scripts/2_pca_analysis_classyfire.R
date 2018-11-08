#------------------------------------------ PCA AND CLASSYFIRE LIPIDOMICS DATA ------------------------------------------#

# Rheumatoid arthritis is a inflammatory disease that characterized for not having a resolution phase in the inflammatory
# procces. One of the most used treatments is Methotrexate (MTX), that is not always efective. 

# We have a dataset of the lipid profiles of rheumatoid arthritis patients that responded or not to the treatment.

# This is the first approach to create a model using machine learning that predicts if a patient will react or not to the
# treatment with MTX. 

#---> LIBRARY LOAD:

library(classyfire)

#---> DATA MANIPULATION: 

# TRAINING SET:

# Open the txt file with the profiles information. 

lm_profiles <- read.table(
  file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/2_pca_and_classyfire/rheumatoid_arthritis_LM_profile.txt",
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t")

lm_profiles <- lm_profiles[-c(16), ] # All the values of Marasein 1 are equal to zero, so it does not add info. 
lm_profiles <- lm_profiles[ ,-c(41, 42)] # The samples PL041 and PL042 were eliminated because there was a problem with
# their extraction.

# Separates the profiles data by lipid mediators types. 

# Fatty acids origin:

dha <- lm_profiles[1:23, ]
n_three_DPA <- lm_profiles[24:33, ]
aa <- lm_profiles[37:54, ]
epa <- lm_profiles[34:36, ]

# Lipid Mediators groups:

resolvins_d <- lm_profiles[1:8, ]
protectins <- lm_profiles[9:12, ]
pctr <- lm_profiles[13:15, ]
maresins <- lm_profiles[16:20, ]
mctr <- lm_profiles[21:23, ]
rvt <- lm_profiles[24:27, ]
rvd <- lm_profiles[28:30, ]
pd <- lm_profiles[31:32, ]
mar <- lm_profiles[33, ]
lx <- lm_profiles[37:41, ]
ltb <- lm_profiles[42:47, ]
lt <- lm_profiles[48:50, ]
pg <- lm_profiles[51:53, ]
tx <- lm_profiles[54, ]

# Combination: 

# Seth did this one, so I did it as well to validate it. 

dha_n_dpa_epa <- lm_profiles[c(1:33, 34:36), ]

#---> DATA PREPARATION:

# Getting the explanatory (x) and response (y) variable. By explanatory, it means all the data that can explain why a
# patient response or not to the treatment (the lipid meadiator profiles) and the response variable is if the patients
# response or not to the treatment. 

# Transpose columns and rows:
lm_profiles_transpose <- t(lm_profiles)


# Explanatory and Response variable:
explanatory <- as.matrix(lm_profiles_transpose)
response <- data.frame(row.names = row.names(lm_profiles_transpose),
                       responses = c(rep("responder", 30), rep("non-responder", 22)))
response_matrix <- as.matrix(response)

# Define a specific color for each of the populations:
classes <- as.factor(response$responses)
names(classes) <- row.names(lm_profiles_transpose)

responder_index<-which((classes=="responder")== TRUE)
non_responder_index<-which((classes=="non-responder")== TRUE)

colors <- NULL
colors[responder_index] <- "green"
colors[non_responder_index] <- "red"

#---> PRINCIPAL COMPONENT ANALYSIS:

# Data Pre-processing:
# Scaling without centering and using the SD as scaling instrument. Since the range of values of raw data varies widely,
# in some machine learning algorithms, objtive functions will not work properly without normalization. For example, 
# if one of the features has a broad range of values, the distance will be governed by this particular feature. 

# If you are working with machine learning, the best method of scalation is standarization. 
explanatory_scale <- scale(explanatory, center = FALSE, scale = TRUE)
write.table(explanatory_scale,
             file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/2_pca_and_classyfire/matrix_SVM/mx_scale_lmprof.txt",
             sep = "\t",
             quote = FALSE)

# PCA:
# prcomp function to run the PCA analysis. 

pca_lm_profile <- prcomp(explanatory_scale)
pca_sum <- summary(pca_lm_profile)

# Explained Variance of each PC:
# "prcom" clusters the features according to Principal Components. Each component then represents part of the data.
# The next part of the scripts shows how much each PC represents the whole data. 

explained_variance <- pca_sum$importance[2,]*100

pdf(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/2_pca_and_classyfire/PC_lm_profiles_graphs.pdf",
    width = 25, height = 12, onefile = TRUE)
barplot(explained_variance, las = 2, col = "steelblue", ylim = c(0,20), cex.axis = 1.5, cex.names = 1.2)
title(xlab = "Principal Components", line = 4, cex.lab = 2)
title(ylab = "Explained Variance (%)", line = 3, cex.lab = 2)


# Cumulative variance of PC:
# See how much Principal components explains more of the data. 

cumulative_variance <- pca_sum$importance[3,]*100
plot(cumulative_variance, type = "o", xlab = "Principal Components", ylab = "Cumulative Variance",
     col = "black", pch = 21, bg = "blue", cex = 0.8, cex.lab = 2, cex.axis = 1.5)
abline(h = 95, col= "red") # Line that shows when the 95% of the data is explained. 
dev.off()

# PC comparison:
# Comparing the values between the PC to see if there is a distribution.

# Based on the Cumulative Variance table, get the number of PC that represents 95% of the data. 

pc_hg95 <- cumulative_variance[cumulative_variance >= 95] # Get the list with the PCs that represents more 95%
min_pc <- names(pc_hg95[1]) # Get the first PC from the new list.
min_pc <- strsplit(min_pc, "PC") # Get the name of PC and separates the number from it.
min_pc <- as.numeric(min_pc[[1]][2]) # Convert the number name as a number. 

concentration_pc <- pca_lm_profile$x # Get the concentrations for the machine learning model.

pdf(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/2_pca_and_classyfire/PC_comparisons.pdf",
    width = 25, height = 12, onefile = TRUE)
pairs(concentration_pc[, 1:(min_pc/2)], pch = 21, bg = colors, cex=0.7, cex.lab=0.7, cex.axis=0.7)
pairs(concentration_pc[, (min_pc/2):min_pc], pch = 21, bg = colors, cex=0.7, cex.lab=0.7, cex.axis=0.7)
dev.off()

#---> MACHINE LEARNING (Classyfire R): 

# Classyfire uses Support Vector Machine to create the Machine Learning Model. SVM classifies, makes a reggresion,
# and creates a novelty detection for the creation of the model. 

# The idea is to create several models and see which one fits the best. The models will be based on the whole
# lipid profiles, the different groups independently and the values from the PCA. 

# "cfBuild" to create the SVM:

support_vm_pc <- cfBuild(concentration_pc[,1:min_pc], response_matrix, bootNum = 70,ensNum = 70, cpus = 4) #PC
support_lmprofiles_scale <- cfBuild(explanatory_scale, response_matrix, bootNum = 70,ensNum = 70, cpus = 4) #scale lm

# Expositional graphs to get info of the model: 

# PCA Model

pdf(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/2_pca_and_classyfire/graphs_SVM/SVM_pc_lm_profiles.pdf",
    width = 14, height = 10, onefile = TRUE)
ggClassPred(support_vm_pc, displayAll = TRUE, fillBrewer = TRUE, showText = TRUE) # Graphs distribution of prediction.
ggEnsTrend(support_vm_pc, ylims = c(50,100)) # Graphs showing how repitions were need to stabilize a prediction %.
dev.off()

# LM profile model

pdf(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/2_pca_and_classyfire/graphs_SVM/SVM_lmprofiles_scale.pdf",
    width = 14, height = 10, onefile = TRUE)
ggClassPred(support_lmprofiles_scale, displayAll = TRUE, fillBrewer = TRUE, showText = TRUE) 
ggEnsTrend(support_lmprofiles_scale, ylims = c(50,100))
dev.off()

# Accuracy table: 

getAvgAcc(support_vm_pc)$Test # Performance average of the model. % of accuracy. 
getConfMatr(support_vm_pc) # Distribution of the responder and non-responder and their predictions 

getAvgAcc(support_lmprofiles_scale)$Test # Get the %CC (Overall percentage of correctly classified test objects)
getConfMatr(support_lmprofiles_scale) # Get a table of the consensus classification of the best model. 

# Creates a table with all the models and the %CC.

accuracy_table <- data.frame(groups = c("PCA", "scalated lm profiles"),
                             percentage_accuracy = c(getAvgAcc(support_vm_pc)$Test, 
                                                     getAvgAcc(support_lmprofiles_scale)$Test),
                             stringsAsFactors = FALSE)


# Save the models as an R object:
# To avoid running all the scripts to obtain the models, it can be saved as R objects. If you want to use it, 
# you can call it using "readRDS". 

saveRDS(support_vm_pc, 
        file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/2_pca_and_classyfire/models_SVM/SVM_pc.R", 
        ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)

saveRDS(support_lmprofiles_scale, 
        file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/2_pca_and_classyfire/models_SVM/SVM_lmprofiles_scale.R", 
        ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)

#---> SVM PER GROUP: 

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
  
  if (nrow(groups[[lm]]) > 1) {
  
  transpose <- as.matrix(t(groups[[lm]])) # Transpose and create matrix.
  scale <- scale(transpose, center = FALSE, scale = TRUE) # All the model will use scale data
  write.table(scale,
               paste(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/2_pca_and_classyfire/matrix_SVM/mx_scale_",
                     names(groups)[[lm]], ".txt", sep = ""),
              sep = "\t",
              quote = FALSE)
  support_vm <- cfBuild(scale, response_matrix, bootNum = 65, ensNum = 65, cpus = 4) # Creates model
  average_right <- data.frame(groups = names(groups)[[lm]], 
                              percentage_accuracy = getAvgAcc(support_vm)$Test) # Save %CC value.
  
  accuracy_table <- rbind(accuracy_table, average_right) # Append %CC to the accuracy table. 
  
  # Model graphs:
  
  pdf(file = paste("C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/2_pca_and_classyfire/graphs_SVM/SVM_",
                   names(groups)[[lm]], ".pdf", sep = ""), 
      width = 14, height = 10, onefile = TRUE)
  print(ggClassPred(support_vm, displayAll = TRUE, fillBrewer = TRUE, showText = TRUE))
  print(ggEnsTrend(support_vm, ylims = c(50,100)))
  dev.off()
  
  # Save model as a R object. 
  
  saveRDS(support_vm, 
          file = paste("C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/2_pca_and_classyfire/models_SVM/SVM_",
                       names(groups)[[lm]], ".R", sep = ""),
          ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
  
  }
  
  else { next }
  
}

# Write table with accuracy values: 

write.table(accuracy_table, 
            file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/2_pca_and_classyfire/accuracy_table.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE) 

# Write the response data for future uses:

write.table(response, 
            file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/2_pca_and_classyfire/response_data.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE) 

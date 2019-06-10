#---------------------------------- EXPLORATORY ANALYSIS OF LIPIDOMICS DATA --------------------------------------------#

# Rheumatoid arthritis is a inflammatory disease that characterized for not having a resolution phase in the inflammatory
# procces. One of the most used treatments is Methotrexate (MTX), that is not always efective. 

# We have a dataset of the lipid profiles of rheumatoid arthritis patients that responded or not to the treatment.

# This is the first approach to analysis this data. 

#---> LIBRARY LOAD:

library(dendextend)
library(gplots)
library(limma)

#---> DATA LOAD: 
# Open the txt file with the profiles information. 

lm_profiles <- read.table(
  file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/input/1_exploratory_analysis/rheumatoid_arthritis_LM_profile.txt",
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t")

lm_profiles <- lm_profiles[-c(16), ] # All the values of Marasein 2 are equal to zero, so it does not add info. 
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

summary(lm_profiles)

#---> DATA PREPARATION:

# Get a list with the names of all the lipid mediators from the table. 

lipid_mediators <- rownames(lm_profiles)

# Create a new data frame with the samples names and the information about if they are responder or  not. 

clasification <- data.frame(samples_names = names(lm_profiles),
                            response = c(1:52))
clasification$response[1:30] <- "responder"      #In this case the samples are in order so the association
clasification$response[31:52] <- "non_responder" # where made mannually. 

# Define the response to the treatment as a factor element and associate it with the samples names. 

classes <- as.factor(clasification$response)
names(classes) <- names(lm_profiles)

# Define a specific color for each of the populations.

responder_index<-which((classes=="responder")== TRUE)
non_responder_index<-which((classes=="non_responder")== TRUE)

colors <- NULL
colors[responder_index] <- "green"
colors[non_responder_index] <- "red"

# Some visualization plots:

# Plot with all the concentrations from the lipid mediators: 

pdf(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/1_exploratory_analysis/lm_profile_boxplot.pdf",
    width = 25, height = 12, onefile = TRUE)
boxplot(lm_profiles, main = "Lm profiles", col = colors, las = 2, xaxt="n")
title(xlab = "Samples", line = 4, cex.lab = 1.5)
title(ylab = "Lipid mediators (pg/ml)", line = 3, cex.lab = 1.5)
axis(side = 1, at = 1:30, col.axis = "green", labels = names(lm_profiles)[1:30], las = 2, cex.axis = 1.2)
axis(side = 1, at = 31:52, col.axis = "red", labels = names(lm_profiles)[31:52], las = 2, cex.axis = 1.2)
legend("topleft", inset=.05, title="Response to treatment:",
       c("Responder", "Non responder"), fill=c("green", "red"),
       text.font = 3, cex = 1.5)
dev.off()

# Plot zoom in of the bloxplot:

#boxplot(lm_profiles, col = colors, las = 2, xaxt="n", ylim = c(0, 100))
#title(xlab = "Samples", line = 4, cex.lab = 1.5)
#title(ylab = "Lipid mediators (pg/ml)", line = 3, cex.lab = 1.5)
#axis(side = 1, at = 1:30, col.axis = "green", labels = names(lm_profiles)[1:30], las = 2, cex.axis = 1.2)
#axis(side = 1, at = 31:52, col.axis = "red", labels = names(lm_profiles)[31:52], las = 2, cex.axis = 1.2)
#legend("topleft", inset=.05, title="Response to treatment:",
#       c("Responder", "Non responder"), fill=c("green", "red"),
#       text.font = 3, cex = 1.5)

#---> CLUSTER ANALYSIS: 

# Hierarchical clustering analysis using a distance and a linkage algorithms: 

# Combination with the best distribution:
# Euclidean distance calculates the distance between point A to B. 
# Ward.D methodology uses the euclidean distances to create the cluster and the linkages. 

hcew <- hclust(dist(t(lm_profiles),method = "euclidean"),method = "ward.D")

# Other combinations used: 

# hcew2 <- hclust(dist(t(lm_profiles),method = "euclidean"),method = "ward.D2")
# hces <- hclust(dist(t(lm_profiles),method = "euclidean"),method = "single")
# hcec <- hclust(dist(t(lm_profiles),method = "euclidean"),method = "complete")
# hcms <- hclust(dist(t(lm_profiles),method = "maximum"),method = "single")
# hcmc <- hclust(dist(t(lm_profiles),method = "maximum"),method = "complete")
# hcmhs <- hclust(dist(t(lm_profiles),method = "manhattan"),method = "single")
# hcmhc <- hclust(dist(t(lm_profiles),method = "manhattan"),method = "complete")

# Dendrogram: 

dend <- as.dendrogram(hcew)

# To put colors to the samples names according to their association:
color_codes <- c(non_responder = "red", responder = "green")
labels_colors(dend) <- color_codes[classes][order.dendrogram(dend)]

pdf(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/1_exploratory_analysis/lm_profile_dendogram.pdf",
    width = 25, height = 12, onefile = TRUE)
par(cex = 1.3)
plot(dend, main="Lm profiles - Distances = Euclidean / Linkages = Ward.D", ylab = "Height")
legend("topright", inset=.05, title="Response to treatment:",
       c("Responder", "Non responder"), fill=c("green", "red"),
       text.font = 3, cex = 1)
dev.off()

# Heatmap: 

# Create a data frame with the lipid mediators groups: 

lm_groups <- data.frame(l_mediators = row.names(lm_profiles),
                            fatty_acids = c(1:54))
lm_groups$fatty_acids[1:23] <- "DHA"      
lm_groups$fatty_acids[24:33] <- "n-3 DPA" 
lm_groups$fatty_acids[34:36] <- "EPA" 
lm_groups$fatty_acids[37:54] <- "AA" 

# Define the fatty acids as a factor element and associate it with the lipid mediators. 

lm_classes <- as.factor(lm_groups$fatty_acids)
names(lm_classes) <- row.names(lm_profiles)

# Heatmap requires a matrix. So, since we are using dataframes, we need to transform them to matrix. 

lm_matrix <- data.matrix(lm_profiles)

# Define a specific color for each of the fatty acids groups:

dha_index<-which((lm_classes=="DHA")== TRUE)
n_three_index<-which((lm_classes=="n-3 DPA")== TRUE)
epa_index<-which((lm_classes=="EPA")== TRUE)
aa_index<-which((lm_classes=="AA")== TRUE)

lm_colors <- NULL
lm_colors[dha_index] <- "blue"
lm_colors[n_three_index] <- "brown"
lm_colors[epa_index] <- "darkgoldenrod1"
lm_colors[aa_index] <- "darkslategray"

# Key color for the heatmap:

breaks <- c(seq(from = 0, to = 100, by = 1), 500, 1000) # To change the format of the Color Key. 
                                                                    # Higher [] clearer the color. 

pdf(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/1_exploratory_analysis/lm_profile_heatmap.pdf",
    width = 25, height = 12, onefile = TRUE)
par(xpd = TRUE)
heatmap.2(lm_matrix,
          hclustfun = function(lm_matrix) hclust(dist(t(lm_matrix),method = "euclidean"),method = "ward.D2"),
          breaks = breaks, # Use the breaks specify before.
          col = redblue(102),
          key.par = list(cex = 1.5),
          dendrogram = c("both"),
          colCol = colors, # Colors from columns
          colRow = lm_colors, # Colors from rows
          density.info = "none",
          trace = "none",
          cexCol = 1.5,
          cexRow = 1.1,
          margins = c(9, 9))
legend(x = 0.75, y =  0.99, inset = .05, title = "Fatty acids:",  bty = "n",
       c("DHA", "n-3 DPA", "EPA", "AA"), fill=c("blue", "brown", "darkgoldenrod1", "darkslategray"),
       text.font = 1, cex = 1.3)
legend(x = 0.85, y =  0.99, inset=.05, title="Response to treatment:",  bty = "n",
       c("Responder", "Non responder"), fill=c("green", "red"),
       text.font = 1, cex = 1.3)
rect(0.75, 0.84, 0.96, 0.99)
dev.off()

#---> DIFFERENTIAL CONCENTRATION:

# Paired t-tests between the population were made to compare the concentration expression. Those Lipid mediators
# associated with low p-values (p-value < 0.05) are defined as lipid mediators with concentrations statistically 
# different between the two groups: responder vs no responder. 

#Create the design using all the responses:

design <- model.matrix(~0+classes) # Create a design where responder are 0 and non-responder are 1. 
fit <- lmFit(lm_profiles, design)  # Organized the values according to the design.
contrasts <- makeContrasts(comp = classesresponder-classesnon_responder,
                           levels = design)  # Make the pairwise comparison using a t-test 
fit2 <- contrasts.fit(fit,contrasts)         
fit3 <- eBayes(fit2)

#Top ten genes that show the most significant difference in expression:

toptable <- topTable(fit3, sort.by = "p",number = 52, adjust.method="BH", coef = "comp")
toptable <- rbind(c("#", "#", "lm_profiles", "#", "#", "#"), toptable)
write.table(toptable, 
            file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/1_exploratory_analysis/lm_diff_conc.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

#---> GRAPHS PER GROUP: 

# Run the same analysis but automatically for the rest of the subgroups. 

# Create a list with the names of all the subgroups:
groups <- list(dha, n_three_DPA, epa, aa,  resolvins_d, protectins, pctr, maresins, mctr, rvt, 
               rvd, pd, mar, lx, ltb, lt, pg, tx)

# Create a vector with the names associated to all the elements in the list: 
names(groups) <- c("DHA", "n-3 DPA", "EPA", "AA", "Resolvins_d", "Protectins", 
                   "PCTRS", "Maresins", "MCTR", "RVT", "D-Series Resolvins", "PD", "Mar", "LX", 
                   "LTB", "Leukotrienes", "Prostaglandins", "TX")

# Because is better to have a single pdf with all the graphs, I open the pdf before the loop: 

pdf(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/1_exploratory_analysis/per_groups_boxplot.pdf",
    width = 25, height = 12, onefile = TRUE)

# The loop goes through all the data frames in the list to create the different plots. It has to be used the [[]]
# since R put the list inside of another list. 

for (lm in 1:length(groups)) {
  
  boxplot(groups[[lm]], main = names(groups)[[lm]], col = colors, las = 2, xaxt="n")
  title(xlab = "Samples", line = 4, cex.lab = 1.5)
  title(ylab = "Lipid mediators (pg/ml)", line = 3, cex.lab = 1.5)
  axis(side = 1, at = 1:30, col.axis = "green", labels = names(lm_profiles)[1:30], las = 2, cex.axis = 1.2)
  axis(side = 1, at = 31:52, col.axis = "red", labels = names(lm_profiles)[31:52], las = 2, cex.axis = 1.2)
  legend("topleft", inset=.05, title="Response to treatment:",
         c("Responder", "Non responder"), fill=c("green", "red"),
         text.font = 3, cex = 1.5)
  
}

dev.off()

# Same thing for the dendrograms.

pdf(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/1_exploratory_analysis/per_groups_dend.pdf",
    width = 25, height = 12, onefile = TRUE)

for (lm in 1:length(groups)) {
  
  hcew <- hclust(dist(t(groups[[lm]]),method = "euclidean"),method = "ward.D")
  dend <- as.dendrogram(hcew)
  color_codes <- c(non_responder = "red", responder = "green")
  labels_colors(dend) <- color_codes[classes][order.dendrogram(dend)]
  
  par(cex = 1.3)
  plot(dend, main = names(groups)[[lm]], ylab = "Height")
  legend("topright", inset=.05, title="Response to treatment:",
         c("Responder", "Non responder"), fill=c("green", "red"),
         text.font = 3, cex = 1)
  
}

dev.off()

# Same thing for the heatmaps:

pdf(file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/1_exploratory_analysis/per_groups_heatmap.pdf",
    width = 25, height = 12, onefile = TRUE)

for (lm in 1:length(groups)) {
  
  # Heatmaps needs more than one row, so it is necessary not to room the program in the data frames with only one
  # Lipid Mediator. 
  
  if (nrow(groups[[lm]]) > 1) {
    
    lm_matrix <- data.matrix(groups[[lm]])
    breaks <- seq(from = 0, to = max(groups[[lm]]), by = max(groups[[lm]])/102)
    par(xpd = TRUE)
    heatmap.2(lm_matrix,
              main = names(groups)[[lm]],
              hclustfun = function(lm_matrix) hclust(dist(t(lm_matrix),method = "euclidean"),method = "ward.D2"),
              breaks = breaks, # Use the breaks specify before.
              col = redblue(102),
              key.par = list(cex = 1.5),
              dendrogram = c("both"),
              colCol = colors, # Colors from columns
              density.info = "none",
              trace = "none",
              cexCol = 1.5,
              cexRow = 1.1,
              margins = c(9, 9))
    legend(x = 0.85, y =  0.99, inset=.05, title="Response to treatment:",  bty = "n",
          c("Responder", "Non responder"), fill=c("green", "red"),
          text.font = 1, cex = 1.3)
  }
  else { next }
}

dev.off()

# Same for the tables with the differential concentrations: 

for (lm in 1:length(groups)) {
  
  if (nrow(groups[[lm]]) > 1) {
  
    design <- model.matrix(~0+classes)
    fit <- lmFit(groups[[lm]], design)
    contrasts <- makeContrasts(comp = classesresponder-classesnon_responder,
                               levels = design)
    fit2 <- contrasts.fit(fit,contrasts)
    fit3 <- eBayes(fit2)
    toptable <- topTable(fit3, sort.by = "p",number = 52, adjust.method="BH", coef = "comp")
    toptable <- rbind(c("#", "#", names(groups)[[lm]], "#", "#", "#"), toptable)
    write.table(toptable, 
                file = "C:/Users/hhy270/Dropbox/Rheumatoid arthritis/output/1_exploratory_analysis/lm_diff_conc.txt",
                sep = "\t",
                quote = FALSE,
                row.names = TRUE,
                append = TRUE) # Append adds new information to previously made files. 
    
  }
  else { next }
}


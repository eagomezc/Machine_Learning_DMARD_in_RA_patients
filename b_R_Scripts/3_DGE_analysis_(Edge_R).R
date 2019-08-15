#---------------------------------- DIFFERENTIAL GENE EXPRESSION ANALYSIS  ----------------------------------------#

# Rheumatoid arthritis is a inflammatory disease that characterized for not having a resolution phase in the 
# inflammatory procces. One of the most used treatments is modifying antirheumatic drugs (DMARDs), that is not 
# always efective. 

# A dataset of Gene Expression level (READ COUNTS) of meaningful enzymes in the lipid mediator pathways of 
# Responder and Non Responder Patients is the input of this script along with a classification table of patients. 

# This is the first approach to analysis of this data. 

#---> LIBRARY LOAD:

library(edgeR)
library(ggplot2)

#---> INPUT AND OUTPUT:

# In this section please specify where are the input files and where you want to save the output files.
# In the input and output variable you can see what is the path expected from the user to write.

input <- "C:/Users/hhy270/Documents/GitHub/2018_Machine_Learning_MTX_treatment_in_RA_patients/a_Toy_Data/3_DGE_analysis_(Edge_R)/"
output <- "C:/Users/hhy270/Documents/GitHub/2018_Machine_Learning_MTX_treatment_in_RA_patients/c_Expected_Output/3_DGE_analysis_(Edge_R)/"

# !!!! IMPORTANT: For this script to work the training dataset has to be called: 3_DGE_analysis_(Edge_R)_read_counts.txt
# !!!! IMPORTANT: For this script to work the validation dataset has to be called: 3_DGE_analysis_(Edge_R)_class_reads.txt

#---> DATA LOAD: 

# Open the txt file with the gene expression information (read counts). 

# The dataset consist in a tab-delimited file in .txt format with the follow specifications: 
# Columns: The different samples (each patient data).
# Rows: The interested enzyme read counts for each patient. 

# See a_Toy_Data/3_DGE_analysis_(Edge_R)/3_DGE_analysis_(Edge_R)_read_counts.txt

counts_blood <- read.table(
  file = paste(input, "3_DGE_analysis_(Edge_R)_read_counts.txt", sep = ""),
  header = TRUE,
  row.names = 1, # Specify that the first column is row names. 
  sep = "\t")

# Open the txt file with the classification table (Responder and Non Responder).

# The table consist in two columns: column "Samples" with the different sample IDs and the column "Response"
# with information for each patient about whether they are responders or not. 

# See a_Toy_Data/3_DGE_analysis_(Edge_R)/3_DGE_analysis_(Edge_R)_class_reads.txt

classification_blood <- read.table(
  file = paste(input, "3_DGE_analysis_(Edge_R)_class_reads.txt", sep = ""), 
  header = TRUE,
  sep = "\t")

#---> DATA PREPARATION:

# Define the response to the treatment as a factor elements and associate them with the samples names.
# This will be helpful at the moment to design the model who is going to identify statistical differents
# between the two groups. 

classes_blood <- as.factor(classification_blood$Response)
names(classes_blood) <- names(counts_blood)

#---> GENE EXPRESSION ANALYSIS (EDGE R):

# Edge R does differential expression analysis of RNA-seq expression profiles with biological replication. 
# Implements a range of statistical methodology based on the negative binomial distributions, including empirical 
# Bayes estimation, exact tests, generalized linear models and quasi-likelihood tests. 

# It takes raw read counts to performance the analysis so pre-normalization steps are not requiered. 

DGE_file_blood <- DGEList(counts = counts_blood, group = classes_blood) # Creates the DGE file. 

# Gene to be expressed at a reasonable level in a sample if two counts per each million mapped reads in that sample.
# Gene should be expressed in at least one condition. 

keep <- rowSums(cpm(DGE_file_blood)>2) >= 10
DGE_file_blood <- DGE_file_blood[keep, , keep.lib.sizes = FALSE]

DGE_file_blood <- calcNormFactors(DGE_file_blood) # TMM Normalization 

# Since it is not a paired comparison and there is not a batch effect, the comparison design is only based on the 
# groups. 

design_edge_blood <- model.matrix(~classes_blood) 

# Estimate the dispersion. robust = TRUE -> Robustified against potential outlier genes.

DGE_file_blood_d <- estimateDisp(DGE_file_blood, design_edge_blood,  robust = TRUE)  

# The quasi-likelihood method was used to calculated differences between the groups.

fit_blood_r <- glmQLFit(DGE_file_blood_d, design_edge_blood) 
lrt_blood <- glmQLFTest(fit_blood_r, coef = 2) # Coef = 2 specify comparison Non Responder vs Responder. 

# The results of the DGE analysis can be seen as a data frame tha cointains for each gene the next information:
# LogFC (Log(FoldChange)).
# LogCPM (Log(CountPerMillion)).
# F value.
# p value and adjust p value (FDR), in this case using the BH correction. 

edge_fc_table_blood <- as.data.frame(topTags(lrt_blood, n = Inf, adjust.method = "BH", sort.by = "PValue"))

# The Output of this section is a tab-delimited table with the information of DGE results: 

write.table(edge_fc_table_blood, 
            file = paste(output, "3_DGE_results.txt", sep = ""),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

# ---> VIOLIN PLOTS:

# This section creates Violin plots for the Alox related genes (ALOX12, ALOX15, ALOX15B and ALOX5).

# Get only the interested genes. 

counts_blood <- counts_blood[c(2:5), ]

# Calculates the Log(cpm) for each gene:

immune_reads <- cpm(counts_blood)

immune_reads_log <- cpm(counts_blood, log =TRUE)
immune_reads_log <- as.data.frame(t(immune_reads_log)) 

immune_reads <- as.data.frame(t(immune_reads))
immune_reads$groups <- classification_blood$Response

# Create a data frame that will contains the requieres information for the violin plots:

violin_table <- data.frame(groups = rep(classification_blood$Response, 4))
genes <- NULL
counts <- NULL

# Small loop that get the information from the counts_genes data frame to the violin_table data frame:

for (i in 1:ncol(immune_reads_log)) {
  genes[i] <- list(rep(colnames(immune_reads_log)[i], ncol(counts_blood)))
  counts[i] <- list(immune_reads_log[, i])
}

violin_table$genes <- unlist(genes)
violin_table$counts <- unlist(counts)

# Specify the factors of the comparison to put them in the violin plot: 

violin_table$genes <- factor(violin_table$genes, levels = colnames(immune_reads_log))
violin_table$groups <- factor(violin_table$groups, levels = c("responder", "non_responder"))

# VIOLIN PLOTS:
# Creates and saves as a pdf the violin plots:

pdf(file = paste(output, "3_Violin_plot_ALOX.pdf", sep = ""), 
    width = 25, height = 12, onefile = TRUE)
ggplot(violin_table, aes(y = counts, x = groups)) +
  geom_violin(aes(fill = groups), trim = FALSE) +
  scale_fill_manual(values = c("lightgreen", "firebrick2")) +
  facet_wrap("genes", scales = "free_y", ncol = 2) +
  scale_y_continuous(name = "Log(CPM)") +
  labs(x = "Response to treatment") +
  theme(plot.title = element_text(size = 30, hjust = 0.5),
        axis.title = element_text(size = 30),
        axis.text.x  = element_text(size = 20, hjust = 0.5),
        axis.text.y  = element_text(size = 30, hjust = 1),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        legend.text  = element_text(size = 20),
        legend.position = "none",
        panel.grid.major = element_blank())
dev.off()


#setwd("/Users/tharun kota/OneDrive/Desktop/Dr.vadim/NTP/immuno/immuno/")

# ICC stim classifying genes
library(readxl)
stim_genes <- read_excel('stim_genes.xlsx')

# CMScaller requires  geneexpressionset  as input
library(Biobase)
library(CMScaller)


load_and_assign <- function(file, new_name) {
  temp_env <- new.env()           # Create a temporary environment
  load(file, envir = temp_env)    # Load data into the temporary environment
  obj_name <- ls(temp_env)        # Get the object name in the file (assuming one object per file)
  assign(new_name, get(obj_name, envir = temp_env), envir = .GlobalEnv)  # Assign to a new name in the global environment
} #function to load data
save_gene_expression <- function(data) {
  # Ensure input data has the required components
  if (!all(c("genes", "samples", "mRNA") %in% names(data))) {
    stop("Input data must contain 'genes', 'samples', and 'mRNA' components.")
  }
  
  # Extract genes, samples, and mRNA data
  genes <- data$genes
  samples <- data$samples
  mRNA <- data$mRNA
  
  # Create a data frame with genes in the first column and mRNA values
  gene_expression_df <- data.frame(Gene = genes, mRNA)
  
  # Set column names, including sample names
  colnames(gene_expression_df)[2:ncol(gene_expression_df)] <- samples
  
  # Return the gene expression data frame
  return(gene_expression_df)
} #  function to extract and save gene expression data
createExpressionSet <- function(data_samples) {
  # Convert the data frame to a matrix, excluding the first column
  exprss_matrix <- as.matrix(data_samples[,-1])
  
  # Set the row names as the gene column for clarity
  rownames(exprss_matrix) <- data_samples$Gene
  
  # Create an ExpressionSet object from the matrix
  minimalset <- ExpressionSet(assayData = exprss_matrix)
  
  # Return the ExpressionSet
  return(minimalset)
} # function to convert data_samples into an ExpressionSet
createExpressionSet2 <- function(data_samples) {
  # Check and fix any issues in the Gene column
  data_samples$Gene[is.na(data_samples$Gene) | data_samples$Gene == ""] <- paste0("gene_", seq_len(sum(is.na(data_samples$Gene) | data_samples$Gene == "")))
  data_samples$Gene <- make.names(data_samples$Gene, unique = TRUE)
  
  # Convert the data frame to a matrix, excluding the first column
  exprss_matrix <- as.matrix(data_samples[,-1])
  
  # Set the row names as the Gene column
  rownames(exprss_matrix) <- data_samples$Gene
  
  # Create an ExpressionSet object
  minimalset <- ExpressionSet(assayData = exprss_matrix)
  
  return(minimalset)
}# if createExpressionset doesn't works
createExpressionSet3 <- function(data_samples) {
  # Check and fix any issues in the Gene column
  data_samples$Gene[is.na(data_samples$Gene) | data_samples$Gene == ""] <- paste0("gene_", seq_len(sum(is.na(data_samples$Gene) | data_samples$Gene == "")))
  data_samples$Gene <- make.names(data_samples$Gene, unique = TRUE)
  
  # Convert the data frame to a matrix, excluding the first column
  exprss_matrix <- as.matrix(data_samples[,-1])
  
  # Set the row names as the Gene column
  rownames(exprss_matrix) <- data_samples$Gene
  
  # Check for NA values in the expression matrix and remove rows with any NAs
  if (any(is.na(exprss_matrix))) {
    exprss_matrix <- exprss_matrix[!rowSums(is.na(exprss_matrix)), ]
  }
  
  # Create an ExpressionSet object
  minimalset <- ExpressionSet(assayData = exprss_matrix)
  
  return(minimalset)
}
filter_by_fdr <- function(data_list, p_val_threshold = 0.05, top_n = 10) {
  # Apply the selection criteria to each dataframe in the list
  lapply(data_list, function(df) {
    # Filter for adjusted p-value < threshold, then select top up- and downregulated genes
    upregulated <- df %>%
      filter(adj.P.Val < p_val_threshold) %>%
      arrange(desc(logFC)) %>%
      slice_head(n = top_n)
    
    downregulated <- df %>%
      filter(adj.P.Val < p_val_threshold) %>%
      arrange(logFC) %>%
      slice_head(n = top_n)
    
    # Combine up- and downregulated genes into one dataframe
    bind_rows(upregulated, downregulated)
  })
}
save_dataframes_to_excel <- function(data_list, file_name = "data_list_output.xlsx", file_path = ".") {
  # Create a new workbook
  wb <- createWorkbook()
  
  # Loop through each dataframe in the list
  for (name in names(data_list)) {
    df <- data_list[[name]]
    
    # Add the dataframe as a new sheet, including row names
    addWorksheet(wb, sheetName = name)
    writeData(wb, sheet = name, x = df, rowNames = TRUE)
  }
  
  # Construct the full file path
  full_file_path <- file.path(file_path, file_name)
  
  # Save workbook to the specified file
  saveWorkbook(wb, file = full_file_path, overwrite = TRUE)
  cat("Data saved to", full_file_path, "\n")
}
plot_stim_cluster_counts <- function(counts, main = "Count of Each STIM Cluster",las) {
  # Create the bar plot and save the bar midpoints
  bar_midpoints <- barplot(counts, 
                           main = main,
                           col = "steelblue",
                           las = las,width = 1,space = 1,
                           ylim = c(0, max(counts)*1.5)) # Adjust y-axis limit to reduce bar height
    
  
  # Add counts above each bar
  text(x = bar_midpoints, 
       y = counts, 
       label = counts, 
       pos = 3,  # Position text above the bars
       cex = 1,  # Text size
       col = "black")
}
plot_stim_cluster_counts(counts = table(stim_genes$`5 STIM clusters`),main = "STIM genes",las = 1)



# Gide_data analysis
load_and_assign('dat.Gide.RData','Gide_data')
Gide_data_samples <- save_gene_expression(Gide_data)
Gide_minimal_set <- createExpressionSet(Gide_data_samples)
Gide_template <- subset(stim_genes, `Gene ID` %in% Gide_data_samples$Gene)
colnames(Gide_template)<- c('probe','class')
Gide_template <- as.data.frame(Gide_template)
# CMS caller
par(mfrow=c(1,1))
Gide_res <- CMScaller(Gide_minimal_set, RNAseq=TRUE, doPlot=TRUE,templates = Gide_template,seed = 1000)
title(main= "Gide_heatmap",outer = FALSE) # heatmap
### limma differential gene expression analysis and visualization
Gide_deg <- subDEG(emat=Gide_minimal_set, class= Gide_res$prediction, doVoom=TRUE,doPairwise = FALSE)
subVolcano(Gide_deg,geneID = 'rownames')
Gide_deg_FDR <- filter_by_fdr(Gide_deg)
save_dataframes_to_excel(Gide_deg_FDR,file_name = "GIDE_deg.xlsx",file_path = '/Users/tharun kota/OneDrive/Desktop/Dr.vadim/NTP/immuno/immuno/Results/Gide_data_results/')
write.xlsx(Gide_res, rowNames = TRUE,colNames = TRUE,file = '/Users/tharun kota/OneDrive/Desktop/Dr.vadim/NTP/immuno/immuno/Results/Gide_data_results/Gide_predictions.xlsx' )
plot_stim_cluster_counts(counts = table(Gide_template$class),main = 'STIM genes in Gide',las = 1)
# Differential gene expression  subDEG wit,h pairwise comparisons
#deg_liu_pr <- subDEG(emat = minimalset_liu, class = res_liu$prediction, doVoom = TRUE, doPairwise = TRUE)
#subVolcano(deg_liu_pr,geneID = 'rownames')



# liu_data analysis
load_and_assign('dat.liu.RData','liu_data')
liu_data_samples <- save_gene_expression(liu_data)
liu_minimal_set <- createExpressionSet(liu_data_samples)
liu_template <- subset(stim_genes, `Gene ID` %in% liu_data_samples$Gene)
colnames(liu_template)<- c('probe','class')
liu_template <- as.data.frame(liu_template)
# CMS caller
par(mfrow=c(1,1))
liu_res <- CMScaller(liu_minimal_set, RNAseq=TRUE, doPlot=TRUE,templates = liu_template,seed = 1000)
title(main= "liu_heatmap",outer = FALSE) # heatmap
### limma differential gene expression analysis and visualization
liu_deg <- subDEG(emat=liu_minimal_set, class= liu_res$prediction, doVoom=TRUE,doPairwise = FALSE)
subVolcano(liu_deg,geneID = 'rownames')
liu_deg_FDR <- filter_by_fdr(liu_deg)
save_dataframes_to_excel(liu_deg_FDR,file_name = 'Liu_Deg.xlsx',file_path = '/Users/tharun kota/OneDrive/Desktop/Dr.vadim/NTP/immuno/immuno/Results/Liu_data_results/')
write.xlsx(liu_res, rowNames = TRUE,colNames = TRUE,file = '/Users/tharun kota/OneDrive/Desktop/Dr.vadim/NTP/immuno/immuno/Results/liu_data_results/Liu_predictions.xlsx' )
plot_stim_cluster_counts(counts = table(liu_template$class),main = 'STIM genes in liu',las = 1)


# riaz_data analysis
load_and_assign('dat.Riaz.RData','Riaz_data')
Riaz_data_samples <- save_gene_expression(Riaz_data)
Riaz_minimal_set <- createExpressionSet3(Riaz_data_samples)
Riaz_template <- subset(stim_genes, `Gene ID` %in% Riaz_data_samples$Gene)
colnames(Riaz_template)<- c('probe','class')
Riaz_template <- as.data.frame(Riaz_template)
# CMS caller
par(mfrow=c(1,1))
Riaz_res <- CMScaller(Riaz_minimal_set, RNAseq=TRUE, doPlot=TRUE,templates = Riaz_template,seed = 1000)
title(main= "Riaz_heatmap",outer = FALSE) # heatmap
### limma differential gene expression analysis and visualization
Riaz_deg <- subDEG(emat=Riaz_minimal_set, class= Riaz_res$prediction, doVoom=TRUE,doPairwise = FALSE)
subVolcano(Riaz_deg,geneID = 'rownames')
Riaz_deg_FDR <- filter_by_fdr(Riaz_deg)
save_dataframes_to_excel(Riaz_deg_FDR, file_name = 'Riaz_deg.xlsx',file_path = '/Users/tharun kota/OneDrive/Desktop/Dr.vadim/NTP/immuno/immuno/Results/Riaz_data_results/') 
write.xlsx(Riaz_res, rowNames = TRUE,colNames = TRUE,file = '/Users/tharun kota/OneDrive/Desktop/Dr.vadim/NTP/immuno/immuno/Results/Riaz_data_results//Riaz_predictions.xlsx' )
plot_stim_cluster_counts(counts = table(Riaz_template$class),main = 'STIM genes in Riaz',las = 1)
















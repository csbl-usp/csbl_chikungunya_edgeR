# Cleaning all workspace 
rm(list = ls())

# Importing libraries
# Define your work directory 
workDir <- "~/<change here>"
if(!require("edgeR")){ # Version 3.24.3
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("edgeR", version = "3.8")
}

datafile <- paste0(workDir, "/fccounts.matrix.infected")
if(file.exists(datafile)) {
  data = read.table(datafile, header=T, row.names=1, com='')
} else {
  stop(paste0("File: ", workDir, "/fccounts.matrix.infected don't exist!"))
}

samplesinfofile <- paste0(workDir, "/samplesinfo.tsv")
if(file.exists(samplesinfofile)) {
  samplesinfo <- read.table(samplesinfofile, header = T, sep = "\t")
} else {
  stop(paste0("File: ", samplesinfofile, " don't exist!"))
}

samplesinfo <- samplesinfo[ c( which( samplesinfo$Class == "CTRL" ),
                               which( samplesinfo$Class == "INF" ) ), ]

samplesinfo <- samplesinfo[ samplesinfo$Sample %in% colnames( data ), ]
data <- data[ , c( as.character( samplesinfo$Sample ) ) ]
samplesinfo <- samplesinfo[ samplesinfo$Sample %in% colnames( data ), ]

# Parameters to determine DEGs
adjPcut <- 0.05
logFCcut <- 1
disp <- 0.1

DEG_analysis <- function(data, samplesinfo, disp, method = c("edgeR"), verbose = F) {
  
  if(method == "edgeR") {
    
    if(verbose) message("Filtering data expression")
    colnames(samplesinfo)
    rnaseqMatrix <- data[ , c( as.character( samplesinfo$Sample ) ) ]
    rnaseqMatrix = round(rnaseqMatrix)
    rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
    repControl <- dim(samplesinfo[ which( samplesinfo$Class == "CTRL" ), ])[1]
    repInfected <- dim(samplesinfo[ which( samplesinfo$Class == "INF" ), ])[1]
    
    if(verbose) message("Creating group factor")
    conditions = factor(c(rep("CTRL", repControl), rep("INF", repInfected)))
    
    if(verbose) message("Creating DGEList object and normalize expression")
    exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
    exp_study = calcNormFactors(exp_study)
    
    if(verbose) message("Analyzing DEG with dispersion")
    et = exactTest(exp_study, pair=c("INF", "CTRL"), dispersion=0.1)
    tTags = topTags(et,n=NULL)
    result_table = tTags$table
    result_table = data.frame(sampleA="INF", sampleB="CTRL", result_table)
    result_table$logFC = -1 * result_table$logFC
    
    return(result_table)
  }
}

table <- DEG_analysis(data = data, samplesinfo = samplesinfo, disp = 0.1, method = "edgeR", verbose = F)

symbols_unique <- read.table(file = paste0(workDir, "/symbols_ensembl_id.tsv"), header = T, sep = "\t")
table$Symbol <- rownames(table)
table <- merge(table, symbols_unique, by = "Symbol", all.x = T)

table_degs      <- table[table$FDR < adjPcut,]
table_degs_up   <- table_degs[table_degs$logFC > logFCcut,]
table_degs_down <- table_degs[table_degs$logFC < -logFCcut,]

write.table(x = data.frame(Symbol = table$geneSymbol.1, table[, c(4:7)]), file = paste0(workDir, "/edgeR_all_FDR_", adjPcut, "_log2FC_", logFCcut, "_infected_VS_control_all.tsv"), quote = F, sep = "\t", row.names = F)
write.table(x = data.frame(Symbol = table_degs_up$geneSymbol.1, table_degs_up[, c(4:7)]), file = paste0(workDir, "/edgeR_UP_FDR_", adjPcut, "_log2FC_", logFCcut, "_infected_VS_control_all.tsv"), quote = F, sep = "\t", row.names = F)
write.table(x = data.frame(Symbol = table_degs_down$geneSymbol.1, table_degs_down[, c(4:7)]), file = paste0(workDir, "/edgeR_DOWN_FDR_", adjPcut, "_log2FC_", logFCcut, "_infected_VS_control_all.tsv"), quote = F, sep = "\t", row.names = F)

message(paste0("Criteria DEG - AdjP: ", adjPcut, " and Log2FC: ", logFCcut,
               "\nNumber of Up-regulated genes: ", dim(table_degs_up)[1],
               "\nNumber of Down-regulated genes: ", dim(table_degs_down)[1]))

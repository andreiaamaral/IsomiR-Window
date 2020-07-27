#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("This R script needs at least 3 arguments to work", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[4] = "out.txt"
}

source("http://bioconductor.org/biocLite.R")

#Load the packages
library(DESeq2)
library("pheatmap")
library("gplots")
library(grid)
library(reshape)

# Define the variables for the arguments that will be used
ID <- args[1]
deseq_file <- args[2]
path_results <- args[4]
paired_unpaired <- args[5]
pvalue <- args[6]

cts <- as.matrix(read.csv(deseq_file,sep="\t",row.names="IsomiR_ID"))
cts[is.na(cts)] = 0

experimental_design_file <- args[3]

coldata <- read.table(experimental_design_file, row.names=1)


  if (paired_unpaired == "Unpaired_Sample") {
    coldata <- rename(coldata, c("V2"="fileName", "V3"="condition"))
    coldata <- coldata[, c("fileName", "condition")]

    dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design =~ condition)
    
    
  } else if (paired_unpaired == "Paired_Sample"){
    
    coldata <- rename(coldata, c("V2"="fileName", "V3"="condition", "V4"="subject"))
    coldata <- coldata[, c("fileName", "condition", "subject")]
    
    #print(coldata)
    
    dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design =~ subject + condition)
    
  }
   

  # Perform a minimal pre-filtering to remove rows that have only 0 or 1 read
  dds <- dds[ rowSums(counts(dds)) > 1, ]

  dds <- DESeq(dds)

  res<-results(dds)

  # Order our results table by the smallest adjusted p value:
  res<-res[order(res$padj),]

  # Define the name of the complete table
  # This table is allways produced
  name_file_table_DE_complete <- paste('CompleteTableDE_',ID,'.txt', sep='')

  # Define the file path
  path_to_table_DE_complete <- paste(path_results, name_file_table_DE_complete, sep='/')

# Write the output DE table

 dds <- estimateSizeFactors(dds)
   norm<-data.frame(counts(dds, normalized=TRUE))
   res_data<-as.data.frame(res)

    res_data$IsomiRs<-row.names(res_data)
    norm$IsomiRs<-row.names(norm)
    table_merge_res<-merge(res_data, norm, by="IsomiRs")

    write.table(table_merge_res, path_to_table_DE_complete, quote=F, row.names=F, col.names=T, sep="\t")
  # Select a subset with the p-value less than 0.1 or 0.05
  if(pvalue == "0.1") {
    
    #differ<-subset(res, res$padj < 0.1)
    differ<-subset(res, res$padj < 0.1 | res$padj != 'NA')
    if (nrow(differ)>=15) {
      top15<-differ[1:15,]
      test<-dds[rownames(top15)]
    }
    else {
      test<-matrix(, nrow=0, ncol=0)
    }
    

  } else if (pvalue == "0.05"){    
    #differ<-subset(res, res$padj < 0.05)
    differ<-subset(res, res$padj < 0.05 | res$padj != 'NA')
    if (nrow(differ)>=15) {
      top15<-differ[1:15,]
      test<-dds[rownames(top15)]
    }
    else {
      test<-matrix(, nrow=0, ncol=0)
    }
    
  }
  

  if (nrow(test)==0) {
    
    # Produce a file justifying the absence of a heatmap
    name_file_no_heatmap <- paste('NoHeatmap_',ID,'.txt', sep='')
    
    # Define the file path
    path_to_file_no_heatmap <- paste(path_results, name_file_no_heatmap, sep='/')
    
    message <- paste("There are no differentially expressed miRNAs/isomiRs for the value", 
                     " of p-value selected. Only the first 100 entries of the complete DE", 
                     " table are shown. No heatmap was produced.",sep='')
    
    # Write the output DE table
    writeLines(message, path_to_file_no_heatmap)
    
    
  } else {
    
    
    # Define the name of the table
    # This table is only produced when there are  
    # differently expressed miRNAs
    name_file_table_DE_test <- paste('TableDE_',ID,'.txt', sep='')
    
    # Define the file path
    path_to_table_DE_test <- paste(path_results, name_file_table_DE_test, sep='/')
    

    # add the normalized counts
    dds <- estimateSizeFactors(dds)
    norm<-data.frame(counts(dds, normalized=TRUE))
    differ_data<-as.data.frame(differ)

    differ_data$IsomiRs<-row.names(differ_data)
    norm$IsomiRs<-row.names(norm)
    table_merge<-merge(differ_data, norm, by="IsomiRs")

    write.table(table_merge, path_to_table_DE_test, quote=F, row.names=F, col.names=T, sep="\t")

    
    ###### MA-plot ####
    
    # Define the names of the plots
    name_file_MAplot_png <- paste('MAplot_',ID,'.png', sep='')
    name_file_MAplot_pdf <- paste('MAplot_',ID,'.pdf', sep='')
    name_file_Heatmap_png <- paste('heatmap_',ID,'.png', sep='')
    
    
    # Define the path of plot
    path_to_MAplot_png <- paste(path_results, name_file_MAplot_png, sep='/')
    path_to_MAplot_pdf <- paste(path_results, name_file_MAplot_pdf, sep='/')
    path_to_Heatmap_png <- paste(path_results, name_file_Heatmap_png, sep='/')
    
    # Create the MA plots in png and in pdf formats
    png(path_to_MAplot_png, width=5*300, height=5*300, res=300, pointsize=8)
    plotMA(res, main="DESeq2", ylim=c(-2,2))
    dev.off()
    pdf(path_to_MAplot_pdf)
    plotMA(res, main="DESeq2", ylim=c(-2,2))
    dev.off()
    
    
    # Heatmap

    vsd <- varianceStabilizingTransformation(test, blind=FALSE)   
    df <- as.data.frame(colData(dds)[,c("condition")])
    
    row.names(df)<-row.names(colData(dds))
    colnames(df)<-("Conditions")
    
    # Write the heatmap to a file
    png(path_to_Heatmap_png, width=1050, height=1000, res=100)
    pheatmap(assay(test), cluster_rows=TRUE,show_colnames = TRUE, show_rownames=TRUE, cluster_cols=TRUE, fontsize=9, fontsize_row = 8, fontsize_col = 8, annotation_col=df)
    dev.off()
  } 

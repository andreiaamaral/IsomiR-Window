![](https://github.com/andreiaamaral/IsomiR-Window/blob/master/IsomirPNGcopy_logo.jpg)


 # Isomir Window is a system for analyzing small-RNA-seq data in an integrative and user-friendly manner

MiRNAs are small-non-coding-RNA molecules with ~22 nt in length that control gene expression. IsomiRs are miRNA variants that vary in length and/or sequence when compared to their canonical forms and can display additions or deletions of one or more nucleotides at the 5' or 3' ends or both, internal editings or 3' end tailings. Consequently, isomiRs might regulate the expression of different targets in comparison to their respective canonical form, affecting the pathways of many known biological processes. 
The IsomiR Window performs the integrated analysis of small-RNA-seq datasets from animals and plants with a user-friendly graphical interface, enabling the discovery of all types of miRNA variants through a Web application. 

The IsomiR Window web-interface allows:

1- Detection and quantification of animal or plant sncRNAs in two different experimental conditions, allowing for the processing of multiple datasets in simultaneous;
2-Identification of all possible types of isomiRs and prediction of their targets and consequent functional impact;
3-Generation of user-interactive plots and charts that can be customized and exported;


 ## Instalation and usage

The IsomiR Window is deployed in a Virtual Machine (download it at http://isomir.fc.ul.pt). For full details about how to install and use the IsomiR Window VM, please download and read the IsomiRWindow_install.pdf

### Prerequisites

There are **NO PREREQUISITES** for **IsomiR Window VM**, all software is already installed within the VM as well as all the IsomiR Window source code. Just follow the  installation manual: IsomiR Window_install.pdf

But in case you want to run the IsomiR Window pipeline in standalone mode without installing the VM you will need the following third software:

BEDTools 2.26.0

Bowtie 1.2.2

DESeq2 1.24.0

GATK 3.8-1-0

HTSeq 0.6.1p1

miRanda v3.3a

miRDeep2 2.0.1.2 (it is advised to download the miRDeep2 2.0.1.2 that we make available here, it is corrected for some bugs that we have found)

miRDP2 1.1.2

SAMtools 1.5

TargetFinder

TargetScan 7.0

topGO 2.26.0



## Using the IsomiR Window standalone pipeline 

The annotation module has the follwoing main features:
-maps reads to reference genome;
-quantifies ncRNAs, miRNAs and IsomiRs
-tests for differential expression
-predicts novel miRNAs
**For step by step instructions go to the IsomiR Window wiki **


# IsomiR-Window
Isomir Window A system for analyzing small-RNA-seq data in an integrative and user-friendly manner

MiRNAs are small-non-coding-RNA molecules with ~22 nt in length that control gene expression. IsomiRs are miRNA variants that vary in length and/or sequence when compared to their canonical forms and can display additions or deletions of one or more nucleotides at the 5' or 3' ends or both, internal editings or 3' end tailings. Consequently, isomiRs might regulate the expression of different targets in comparison to their respective canonical form, affecting the pathways of many known biological processes. 
The IsomiR Window performs the integrated analysis of small-RNA-seq datasets from animals and plants with a user-friendly graphical interface, enabling the discovery of all types of miRNA variants through a Web application. 

The IsomiR Window web-interface allows:

1- Detection and quantification of animal or plant sncRNAs in two different experimental conditions, allowing for the processing of multiple datasets in simultaneous;
2-Identification of all possible types of isomiRs and prediction of their targets and consequent functional impact;
3-Generation of user-interactive plots and charts that can be customized and exported;


 ## Instalation and usage

The IsomiR Window is deployed in a Virtual Machine (download it at http://isomir.fc.ul.pt). For full details about how to install and use the IsomiR Window VM, please read the IsomiRWindow_install.pdf

### Prerequisites

There are NO PREREQUISITES for IsomiR Window VM, all software is already installed within the VM as well as all the IsomiR Window source code. Just follow the  installation manual: IsomiR Window_install.pdf

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



## Using the IsomiR Window standalone pipeline for the annotation module

The annotation module has the follwoing main features:
-maps reads to reference genome;
-quantifies ncRNAs, miRNAs and IsomiRs
-tests for differential expression
-predicts novel miRNAs

### Step 1 - create the directory tree required for output and log processing files

/results_IsomiR_window

then inside this one create subdirectories

/bowtie_stats --> files with bowtie alignment statistics

/Log_processing_files --> files with alignment log

/Log_done_files --> log for rRNA filtering

/filtered_SAM --> in this directory  results of the alignment file are produced. *filtered.sam: alignment file filtered for rRNA reads

/individual_SAM_results --> original alignment file before filtering for rRNA

/Htseq_outputs/individual_ncRNAs --> counts per sncRNA feature. miRNAs are not counted at this stage

/Htseq_outputs/counts_ncRNAs --> counts per snRNA category.

### Step 2 - create the directory tree required for input files 

/fastq_file --> directory with fastq files (files should be trimmed for adapters and advised to be already filtered for read quality)

/Knowledge_bases --> directory with genome assemblies (fasta and indexed genomes), RNAcentral.gff3 and miRBase.gff (to receive these knowledge bases please send email to isomiR_window1.0@gmail.com write species code in subject. Then you will receive a link to transfer all the required files.

Species code 
            
            Bos_taurus
            Canis_familiaris
            Capra_hircus
            Danio_rerio
            Drosophila_melanogaster
            Equus_caballus
            Gallus_gallus
            Homo_sapiens
            Mus_musculus
            Ovis_aries
            Pan_troglodytes
            Rattus_norvegicus
            Sus_scrofa
            
            Arabidopsis_thaliana
            Brachypodium_distachyon
            Oryza_sativa
            Solanum_lycopersicum
            Sorghum_bicolor
            Vitis_vinifera
            Zea_mays


### Step 3 - download and run find_ncRNAs.pl
 It requires Species.pl and species_hash.txt
 NM= number of mismatches for alignment
 NMH= number of allowed multiple hits for aligment
 Bowtie runs in --best strata-- setting
 
 running command
 
    perl <directoryforscript>/find_ncRNAs.pl  <pathto>/Knowledge_bases./ <pathto>/fastq_file <pathto>/results_IsomiR_window/bowtie_stats, <pathto>/Log_processing_files <pathto>/Log_done_files <pathto>/filtered_SAM <pathto>/individual_SAM_results <pathto>/Htseq_outputs/individual_ncRNAs <pathto>/Htseq_outputs/counts_ncRNAs <pathto>/fastq_file/<name_fastq_file> <NM> <NMH> <species code>
 
 
### Step 4 -create the directory tree required for result outputs of find_isomiRs.pl
/find_isomirs_results
/find_isomirs_results/isomirs
/find_isomirs_results/isomirs/ambiguous
/find_isomirs_results/isomirs/not_ambiguous
/find_isomirs_results/SAMs_VCFs

### Step 5 -download and run find_isomiRs.pl
 It requires Species.pl and species_hash.txt
 
 running command
 
    perl <directoryforscript>/find_isomiRs.pl <pathto>/Knowledge_bases <pathto>/filtered_SAM <pathto>/individual_SAM_results <pathto>/not_ambiguous <pathto>/ambiguous <pathto>/filtered_SAM/<nameof_filtered_SAMfile> <species code>
    
    
### Step 6 - Download scripts: Deseq.pl and Rdeseq.R. These scripts prepare matrix of isomiR counts and of ncRNA counts accross all replicates in two experimental conditions and perform testing for differential expression using Deseq. 

this script has been created to generate read count matrices for ncRNAs and for IsomiRs of two experimental conditions so it requires the existance of the following directories:

C1/find_isomirs_results/isomirs/not_ambiguous
C2/find_isomirs_results/isomirs/not_ambiguous
/C1/Htseq_outputs/individual_ncRNAs
/C2/Htseq_outputs/individual_ncRNAs
/input_DESeq
/log_DESeq



running command

     perl <directoryforscript>/deseq.pl <jobID> <pathto>/C1/find_isomirs_results/isomirs/not_ambiguous <pathto>/C2/find_isomirs_results/isomirs/not_ambiguous <pathto>/C1/Htseq_outputs/individual_ncRNAs <pathto>/C2/Htseq_outputs/individual_ncRNAs /input_DESeq /log_DESeq <paired or unpaired>

after this step the required inputs for differential expression testing for IsomiRs and for ncRNAs (written at/input_DESeq) are prepared
create diretory to send Deseq analysis output.
/output_DESeq
 
                                                                                                                                                   
then
running command

    Rscript Rdeseq.R <jobID> /input_DESeq /output_DESeq <choose design: paired or unpaired> <choose p-value: 0.05 or 0.10>

The outputs are:
a) file CompleteTableDE_'<jobID>'.txt --> a Tab delimited matrix with the results of the Differential Expression Test order by p-value (smallest to biggest) and including the normalized counts of each replicate in each condition.
b) file TableDE_'<jobID>'.txt' --> a Tab delimited matrix with the results of the Differential Expression Test, for IsomiRs with the p-value <= value set for the analysis including the normalized counts of each replicate in each condition.
                                                                                                                                                    c) file MAplot_'<jobID>'.png --> MA plot
d) file MAplot_'<jobID>'.pdf --> MA plot
e) file heatmap_'<jobID>'.png --> heatmap displaying hierarquical clustering of samples including the top15 differentially expressed IsomiRs.
    
### Step 7 - Download charts.pl and run it to summary stats of previous analysis: read length distribution, ncRNA counts, isomiRs Biogenesis, frequency of IsomiR types and miRNA normalized count matrix 

this script allows to create summary tables that can be used to produce bar charts of read length distribution, piecharts of the proportion of ncRNAs and of miRNA biogenesis (3' or 5' arm) and it further creates a matrix of normalized counts at miRNA level aggregating the complexity of isomIRs.

requires species.pl

requires creating a directory for results:
/charts_C1
/charts_C2
/summary_tables

command

     perl charts.pl <jobID> <pathto>/C1/find_isomirs_results/isomirs/not_ambiguous <pathto>/C2/find_isomirs_results/isomirs/not_ambiguous <pathto>/C1/find_isomirs_results/isomirs/ambiguous <pathto>/C2/find_isomirs_results/isomirs/ambiguous 
	<pathto>/C1/Htseq_outputs/individual_ncRNAs <pathto>/C2/Htseq_outputs/individual_ncRNAs <pathto>/charts_C1 <pathto>/charts_C2 <pathto>/summary_tables <speciescode>
 
 The outputs are:
 
 Total_ncRNAs_piechart_C#.txt
 Types_isomiRs_C#.txt
 Counts_3p5p_C#.txt
 Counts_editings_C#.txt
 Counts_tailings_C#.txt
 Counts_families_C#.txt
 Final_TableDE_<jobID>.txt --> summary statistics of isomiR counts for each condition (mean and stdev)
 miRNA_normalized_counts_<jobID>.txt --> normalized miRNA counts for each replicate in each condition

### Step 8 -Download script mirdeep.pl and merged_arf_and_collapsed_id.pl to perform miRNA prediction and incorporation of novel miRNAs in the species miRNA annotation
The IsomiR Window pipeline provides a script that can be used before step 5. This script enables de identification of novel miRNAs based in your data and allows if you desire you can add these novel miRNAs to the annotation of miRNAs of the corresponding species. In case you execute this script before find_isomiRs.pl, you will be able to identify the isomiRs that have been generated for these novel miRNAs as well. The script executes automatically miRdeep2 in caso of animal genomes or miRD2 in case of plant genomes. These predictions sometimes predict overlapping miRA percursors that overlap >=90%, the script creates an output that removes this redundancy, displaying in this case only the predictions with better score and better sequencing coverage. Nevertheless the full set of results is kept.

You will see that in the command flag options <pathto>/filtered_SAM is repeated twice, this is due to the fact that the IsomiR Window web interface is designed for the comparison of two experimental conditions. In case you have that design  you will use <pathto>/C1/filtered_SAM and <pathto>/C2/filtered_SAM, otherwise you just repeat the same path twice. 
	
Requires: merged_arf_and_collapsed_id.pl

Before you run the script you must create two new folders:

/find_isomirs_results/miRprediction_input --> here it will be written the input files required for miRDeep2 (miRNA prediction for animal genomes) or miRDP2 (miRNA prediction for plant genomes) and that are created from the SAM file input 

/find_isomirs_results/miRprediction_results--> here it will be written the results of the novel miRNA prediction 

running command


    perl <directoryforscript>/mirdeep.pl <jobID> <pathto>/Knowledge_bases <pathto>/filtered_SAM <pathto>/filtered_SAM /find_isomirs_results/miRprediction_input /find_isomirs_results/miRprediction_results yes <species code> <'C1', 'C2' ou 'both'> <minimum read stack (10, 50 ou 100)> <use annotation in prediction ("yes", "no")> <add novel miRNAs to annotation ("yes", "no")>
    
## Using the IsomiR Window standalone pipeline for the functional module

### Step 1 -Download functional.pl to perform miRNA target and gene set enrichmnet analysis on the set of identified targets 
This script performs the following automated analysis: a) target prediction of miRNAs or isomiRs by invoking Miranda and/or TargetScan in the case of animal genomes or by invoking Targetfinder in the case of plant genomes.

requires:

Species.pl
R_topGO.R
targetfinder.pl
targetscan_70.pl

   `perl functional.pl `

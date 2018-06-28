# Genome-wide repressive capacity of promoter DNA methylation is revealed through epigenomic manipulation

This repository contains all of the code to produce results and figures in our 
manuscript "Genome-wide repressive capacity of promoter DNA methylation is revealed through epigenomic manipulation." Our paper is a re-analysis of the data presented by Ford et al., 
"Frequent lack of repressive capacity of promoter DNA methylation identified 
through genome-wide epigenomic manipulation" (
available on [BioRxiv](https://www.biorxiv.org/content/early/2017/09/20/170506)).

We come to a contradictory conclusion using the data made available by the 
authors. Specifically, we observe that forced promoter methylation overwhelmingly
results in suppression of gene expression, as well as reduced complementary 
chromatin marks of active transcription.

## File contents

- `Rmd`: Main analysis files (`.Rmd` format) as well as the R script for formating the 
main manuscript figures 

    - `reproduceFigure5.Rmd`: Using processed data in supplementary tables of Ford et al. (2017), reproduce figure 5
    
    - `mCG-RNAseq-analysis.Rmd`: Using methylation counts and RNA-Seq counts from GEO, investigate the relationship between promoter DNA methylation and gene expression
        
    - `mCG-ChIPbs-analysis.Rmd`: Using raw ChIP-bisulfite sequencing reads, investigate the relationship between promoter DNA methylation and active chromatin marks (H3K4me3 and RNA PolII)
    
    - `figures.R`: Compile final figures for manuscript using outputs of the three `Rmd` files

- `R`: Helper scripts with functions for files in `Rmd` 

    - `read.lister`: read methylation count data from GEO into R 
    
    - `meanDiff.dss.R`: calculate mean methylation differences for DSS regions
    
    - `removeGeom.R`: remove elements from ggplot2 objects 

- `img`: Relevant figures from Ford et al. (2017) for comparison

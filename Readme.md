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
main manuscript figures are in the `Rmd` directory. 

- `R`: Helper scripts with functions for reading in raw data, calculating mean methylation
differences over regions, and removing elements from ggplot2 objects are found in the
`R` directory.

- `img`: Relevant figures from Ford et al. (2017) are included in the `img` directory.

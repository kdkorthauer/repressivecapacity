#!/bin/bash
#SBATCH -J knit-lister              # A single job name for the array
#SBATCH -n 4                        # Number of cores
#SBATCH -N 1                        # All cores on one machine
#SBATCH -p serial_requeue                   # Partition - shared
#SBATCH --mem 200G                  # Memory request - 125G
#SBATCH -t 2-0:00                  # Maximum execution time (D-HH:MM)- 0-18:00
#SBATCH -o ../SLURM/render-%j.out      # Standard output
#SBATCH -e ../SLURM/render-%j.err      # Standard error
 
export RSTUDIO_PANDOC="/n/sw/fasrcsw/apps/Core/rstudio/0.98.1103-fasrc01/bin/pandoc/"
ulimit -u 2000 
module load sratoolkit
module load macs2/2.1.1.20160309-fasrc02
module load fastqc
module load cutadapt
module load jdk/1.8.0_45-fasrc01
module load samtools
module load bismark

# change filename to Rmd to be knitted. 
# make sure ncores in Rmd matches -n batch param above

#R -e "rmarkdown::render('methylationCountAnalysis.Rmd', clean = TRUE)"
R -e "rmarkdown::render('reproduceFigure3.Rmd', clean = TRUE)"
#R -e "rmarkdown::render('reproduceFigure5.Rmd', clean = TRUE)"

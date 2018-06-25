
# set up workspace
library(ggplot2)
library(cowplot)
library(R.utils)

setwd("../plots/")

# rscript directory
rdir <- "../R/"
sourceDirectory(rdir)

### FIGURE 1
# load sub figures
A1 <- readRDS(file="A1.rds") + ggtitle("Gene - UMR pairs") +
             ylab("log2 mRNA Fold Change (+dox/no-dox)") +
             xlab(expression(paste(Delta, "mCG in UMR (+dox - no-dox)")))
A2 <- readRDS(file="A2.rds") + ggtitle("Gene - DMR pairs")  +
             ylab("log2 mRNA Fold Change (+dox/no-dox)") +
             xlab("dmrseq test statistic (+dox - no-dox)")
B1 <- readRDS(file="B1.rds") +
             ggtitle(expression(paste("Gene - UMR pairs (", Delta, "mCG > 0.3)"))) +
             xlab("log2 mRNA Fold Change (+dox/no-dox)")
B2 <- readRDS(file="B2.rds") +
             ggtitle("Gene (FDR < 0.05) - DMR pairs (FDR < 0.01)") +
             xlab("log2 mRNA Fold Change (+dox/no-dox)") 

B1bw <- readRDS(file="B1bw.rds") +   ylab("Density") +
  ggtitle(expression(paste("Gene - UMR pairs (", Delta, "mCG > 0.3)")))
B2bw <- readRDS(file="B2bw.rds") +   ylab("Density") +
  ggtitle("Gene (FDR < 0.05)  - DMR pairs (FDR < 0.01)") 

leg_A <- get_legend(A1)
leg_B <- get_legend(B1)

# remove smoothing line 
A2 <- remove_geom(A2, "GeomSmooth")

# plot 4 panel with black and white histograms
plot_grid(A1 + theme(legend.position="none"), B1 + theme(legend.position="none"),
          A2 + theme(legend.position="none"), B2+ theme(legend.position="none"),
          ncol = 2, labels = LETTERS[1:4])
ggsave("fig1.pdf", width=8.6, height=7)


### FIGURE 2
# load subfigures
A <- readRDS(file="mCG_vs_H3K4me3_scatter_DMR.rds") + xlim(0,0.8) +
      ylab("log2 Fold Change (+dox/no-dox)") +
      ggtitle("H3K4me3 ChIP-BS") +
      xlab(expression(paste(Delta, "mCG in DMR (+dox - no-dox)")))
B <- readRDS(file="stat_vs_H3K4me3_scatter_DMR.rds") + xlim(0,40) +
      ylab("log2 Fold Change (+dox/no-dox)") +
      ggtitle("H3K4me3 ChIP-BS")
C <- readRDS(file="mCG_vs_RNApolII_scatter_DMR.rds") + xlim(0,0.8) +
      ylab("log2 Fold Change (+dox/no-dox)") +
      ggtitle("RNA PolII ChIP-BS") +
      xlab(expression(paste(Delta, "mCG in DMR (+dox - no-dox)")))
D <- readRDS(file="stat_vs_RNApolII_scatter_DMR.rds")  + xlim(0,40) +
       ylab("log2 Fold Change (+dox/no-dox)") +
       ggtitle("RNA PolII ChIP-BS")

plot_grid(A, B, C, D,
          ncol = 2, labels = LETTERS[1:4])
ggsave("fig2.pdf", width=8.6, height=7)

### FIGURE 3
# load subfigures
A <- readRDS(file="mCG_H3K4me3_boxplot_byDMR.rds")
B <- readRDS(file="mCG_RNApolII_boxplot_byDMR.rds")

leg_A <- get_legend(A + labs(color = "Treatment") + 
                      scale_color_manual(values = c("darkblue", "darkred"),
                                         labels = c("no-dox", "+dox")))

p3 <- plot_grid(A + theme(legend.position="none"), B + theme(legend.position="none"),
          ncol = 1, labels = LETTERS[1:2])
plot_grid(p3, leg_A, ncol = 2, rel_widths = c(1, 0.25))
ggsave("fig3.pdf", width=7.5, height=7)

###
### SUPPLEMENTARY
### 


### FIGURE S1

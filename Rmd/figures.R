
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

### FIGURE 4
f4 <- readRDS(file="../plots/barplot_CGI.rds")
f4 + xlab("CG density of methylated promoter")
ggsave("fig4.pdf", width=5, height=3)


###
### SUPPLEMENTARY
### 


### FIGURE S1
# load subfigures
A <- readRDS(file="expdist_raw.rds") + ggtitle("")
B <- readRDS(file="expdist_global.rds") + ggtitle("")
C <- readRDS(file="expdist_control.rds") + ggtitle("")

leg_A <- get_legend(A + theme(legend.position="bottom"))

pS1 <- plot_grid(A + theme(legend.position="none"), 
                 B + theme(legend.position="none"),
                 C + theme(legend.position="none"),
                ncol = 1, labels = LETTERS[1:3])
plot_grid(pS1, leg_A, ncol = 1, rel_heights = c(1, 0.05))
ggsave("figS1.pdf", width=4.5, height=9)

### FIGURE S2 - scatterplot meth expr faceted by CGI status
A <- readRDS(file="A2_cgi.rds") + ylab("log2 Fold Change (+dox/no-dox)") 
B <- readRDS(file="B2_cgi.rds") + xlab("log2 mRNA Fold Change (+dox/no-dox)") 
C <- readRDS(file="A2_noncgi.rds") + ylab("log2 Fold Change (+dox/no-dox)") 
D <- readRDS(file="B2_noncgi.rds") + xlab("log2 mRNA Fold Change (+dox/no-dox)") 

leg_A <- get_legend(A + theme(legend.position="bottom"))

pS2 <- plot_grid(A + theme(legend.position="none") + ggtitle("") + expand_limits(x=36), 
                 B + theme(legend.position="none") + ggtitle(""),
                 C + theme(legend.position="none") + ggtitle("") + expand_limits(x=36),
                 D + theme(legend.position="none") + ggtitle(""),
                 ncol = 2, labels = c("A", "", "B", ""))
pS2 + draw_label("non-CG island promoters", size = 13, x=0.5, y=0.475) + 
      draw_label("CG island promoters", size = 13, x=0.5, y=0.975)
ggsave("figS2.pdf", width=8.6, height=7)

### FIGURE S3 - scatterplot h3k4 and RNA pol II faceted by CGI status
# load subfigures
A <- readRDS(file="stat_vs_H3K4me3_scatter_DMR_CGI.rds") + xlim(0,40) +
  ylab("log2 Fold Change (+dox/no-dox)") +
  ggtitle("H3K4me3 ChIP-BS, CG Islands only")
B <- readRDS(file="stat_vs_H3K4me3_scatter_DMR_nonCGI.rds")  + xlim(0,40) +
  ylab("log2 Fold Change (+dox/no-dox)") +
  ggtitle("H3K4me3 ChIP-BS, non CG Islands only")
C <- readRDS(file="stat_vs_RNApolII_scatter_DMR_CGI.rds") + xlim(0,40) +
  ylab("log2 Fold Change (+dox/no-dox)") +
  ggtitle("RNA PolII ChIP-BS, CG Islands only")
D <- readRDS(file="stat_vs_RNApolII_scatter_DMR_nonCGI.rds")  + xlim(0,40) +
  ylab("log2 Fold Change (+dox/no-dox)") +
  ggtitle("RNA PolII ChIP-BS, non CG Islands only")

plot_grid(A, B, C, D,
          ncol = 2, labels = LETTERS[1:4])
ggsave("figS3.pdf", width=8.6, height=7)


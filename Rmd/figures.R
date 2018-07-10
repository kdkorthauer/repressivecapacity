
# set up workspace
library(ggplot2)
library(cowplot)
library(R.utils)
library(dplyr)

setwd("../plots/")

# results directory
resdir <- "../RESULTS"

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
  ggtitle("Gene - DMR pairs (FDR < 0.01)") +
  xlab("log2 mRNA Fold Change (+dox/no-dox)") 


leg_A <- get_legend(A1)
leg_B <- get_legend(B1)

# remove smoothing line 
#A2 <- remove_geom(A2, "GeomSmooth")

# plot 4 panel with black and white histograms
plot_grid(A1 + theme(legend.position="none"), B1 + theme(legend.position="none"),
          A2 + theme(legend.position="none"), B2+ theme(legend.position="none"),
          ncol = 2, labels = LETTERS[1:4])
ggsave("fig1.pdf", width=8.6, height=7)


### FIGURE 2
A <- readRDS("../plots/odds.plot.list.rds") 

A[[2]] + ggtitle("") +
  coord_cartesian(ylim = c(1, 60)) 
ggsave("fig2.pdf", width = 6, height = 4)


## FIGURE 3 
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
ggsave("fig3.pdf", width=8.6, height=7)

### FIGURE 4
# load subfigures
A <- readRDS(file="mCG_H3K4me3_boxplot_byDMR.rds")
B <- readRDS(file="mCG_RNApolII_boxplot_byDMR.rds")

leg_A <- get_legend(A + labs(color = "Treatment") + 
                      scale_color_manual(values = c("darkblue", "darkred"),
                                         labels = c("no-dox", "+dox")))

p3 <- plot_grid(A + theme(legend.position="none"), B + theme(legend.position="none"),
                ncol = 1, labels = LETTERS[1:2])
plot_grid(p3, leg_A, ncol = 2, rel_widths = c(1, 0.25))
ggsave("fig4.pdf", width=7.5, height=7)

### FIGURE 5
f4 <- readRDS(file="../plots/barplot_CGI.rds")
f4 + xlab("CG density of methylated promoter")
ggsave("fig5.pdf", width=5, height=3)


###
### SUPPLEMENTARY
### FIGURE S1 - figure 2 by cutoff for expressed 
df <- readRDS(file.path(resdir, "odds.table.rds"))
df <- df %>% 
  mutate(count.on = paste0("Min count expressed genes: ", count.on))
df$count.on <- factor(df$count.on)
df$count.on <- factor(df$count.on, levels = levels(df$count.on)[c(1,6,2:5)])

ggplot(df,
       aes(x = fdr.de, y = odds, group = fdr.dmr, color = fdr.dmr, 
           fill = fdr.dmr)) +
  scale_y_continuous(trans= "log2" ) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_ribbon(aes(ymin = ci.lower, ymax = ci.upper), alpha = 0.3, linetype=0) +
  geom_point() + 
  geom_line() +
  theme_bw() +
  xlab("DESeq2 FDR") +
  ylab("Odds of gene repression") +
  labs(color = "dmrseq FDR", fill= "dmrseq FDR")  +
  facet_wrap(~ count.on, ncol = 2) +
  coord_cartesian(ylim = c(1, 60)) 
ggsave("../plots/figS1.pdf", width = 7, height = 7)


### FIGURE S2 - figure 2 by number of genes
df <- readRDS(file.path(resdir, "odds.table.rds"))
df <- df %>% 
  mutate(count.on = paste0("Min count expressed genes: ", count.on))
df$count.on <- factor(df$count.on)
df$count.on <- factor(df$count.on, levels = levels(df$count.on)[c(1,6,2:5)])

ggplot(df,
       aes(x = fdr.de, y = n, group = fdr.dmr, color = fdr.dmr, 
           fill = fdr.dmr)) +
  geom_point() + 
  geom_line() +
  theme_bw() +
  xlab("DESeq2 FDR") +
  ylab("Number of genes") +
  labs(color = "dmrseq FDR", fill= "dmrseq FDR") +
  facet_wrap(~ count.on, ncol = 2)
ggsave("../plots/figS2.pdf", width = 7, height = 7)

### FIGURE S3 - expression distributions by normalization method
# load subfigures
A <- readRDS(file="expdist_raw.rds") + ggtitle("")
B <- readRDS(file="expdist_global.rds") + ggtitle("")
C <- readRDS(file="expdist_control.rds") + ggtitle("")

leg_A <- get_legend(A + theme(legend.position="bottom"))

pS3 <- plot_grid(A + theme(legend.position="none"), 
                 B + theme(legend.position="none"),
                 C + theme(legend.position="none"),
                 ncol = 1, labels = LETTERS[1:3])
plot_grid(pS3, leg_A, ncol = 1, rel_heights = c(1, 0.05))
ggsave("figS3.pdf", width=4.5, height=9)

### FIGURE S4 - dist of fold changes of extra DE genes
A <- readRDS(file.path("../plots/extraDEgenes.rds"))
A
ggsave("figS4.pdf", height=3, width =7)

### FIGURE S5 - scatterplot meth expr faceted by CGI status
A <- readRDS(file="A2_cgi.rds") + ylab("log2 Fold Change (+dox/no-dox)") 
B <- readRDS(file="B2_cgi.rds") + xlab("log2 mRNA Fold Change (+dox/no-dox)") 
C <- readRDS(file="A2_noncgi.rds") + ylab("log2 Fold Change (+dox/no-dox)") 
D <- readRDS(file="B2_noncgi.rds") + xlab("log2 mRNA Fold Change (+dox/no-dox)") 

leg_A <- get_legend(A + theme(legend.position="bottom"))

pS4 <- plot_grid(A + theme(legend.position="none") + ggtitle("") + expand_limits(x=36), 
                 B + theme(legend.position="none") + ggtitle(""),
                 C + theme(legend.position="none") + ggtitle("") + expand_limits(x=36),
                 D + theme(legend.position="none") + ggtitle(""),
                 ncol = 2, labels = c("A", "", "B", ""))
pS4 + draw_label("non-CG island promoters", size = 13, x=0.5, y=0.475) + 
  draw_label("CG island promoters", size = 13, x=0.5, y=0.975)
ggsave("figS5.pdf", width=8.6, height=7)

### FIGURE S6 - figure 2 by CGI status
df <- readRDS(file.path(resdir, "odds.table.cgi.rds"))
df <- df %>% 
  mutate(ci.upper = pmin(ci.upper, 60),
         cgi = ifelse(cgi, "CG Island", "non-CG Insland")) 
ggplot(df %>% filter(count.on == 50),
       aes(x = fdr.de, y = odds, group = fdr.dmr, color = fdr.dmr, 
           fill = fdr.dmr)) +
  scale_y_continuous(trans= "log2") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_ribbon(aes(ymin = ci.lower, ymax = ci.upper), alpha = 0.3, linetype=0) +
  geom_point() + 
  geom_line() +
  theme_bw() +
  xlab("DESeq2 FDR") +
  ylab("Odds of gene repression") +
  labs(color = "dmrseq FDR", fill= "dmrseq FDR") +
  facet_wrap(~ cgi)  +
  coord_cartesian(ylim = c(1, 60)) 
ggsave("../plots/figS6.pdf", width = 8.6, height = 4)

### FIGURE S7 - scatterplot h3k4 and RNA pol II faceted by CGI status
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
ggsave("figS7.pdf", width=8.6, height=7)


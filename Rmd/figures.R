
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

# transf fxn
neglog10 <- function(x) -log10(x)
neglog10_inv <- function(x) 10^(-x)

nl_breaks <- function() {
  function(x) c(0.005, 0.01, 0.05, 0.10)
}

nl_format <- function(x) {
  as.numeric(x)
}

scale_type.nl <- function(x) "neglog10"

scale_nl <- function(aesthetics, ...) {
  continuous_scale(aesthetics, "neglog10", identity,
                   guide = "none", trans = neglog10(), ...)
}

scale_x_nl <- function(...) {
  scale_nl(aesthetics = c("x", "xmin", "xmax", "xend"), ...)
}

scale_y_nl <- function(...) {
  scale_nl(aesthetics = c("y", "ymin", "ymax", "yend"), ...)
}


nl_trans <- function() {
  scales::trans_new("neglog10", transform = neglog10, inverse = neglog10_inv,
            breaks = nl_breaks(), format = nl_format)
}

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
df <- readRDS(file.path(resdir, "odds.table.rds")) %>% 
  filter(count.on == 50) %>% 
  filter(fdr.de %in% c(0.01, 0.05, 0.10, 0.5, 1) )
df$fdr.dmr <- as.numeric(as.character(df$fdr.dmr))
df$fdr.de[df$fdr.de == 1] <- "None (all genes)"
df$fdr.de <- factor(df$fdr.de)

ggplot(df,
       aes(x = fdr.dmr, y = prop, group = fdr.de, color = fdr.de, 
           fill = fdr.de)) +
  #scale_y_continuous(trans= "log2") +
  #geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_ribbon(aes(ymin = p.ci.lower, ymax = p.ci.upper), alpha = 0.3, linetype=0) +
  geom_point() + 
  geom_line() +
  viridis::scale_color_viridis(discrete = TRUE, direction = -1) +
  viridis::scale_fill_viridis(discrete = TRUE, direction = -1) +
  scale_x_continuous(trans=nl_trans()) +
  theme_bw() +
  xlab("dmrseq FDR level") +
  ylab("Probability of being significantly repressed") +
  labs(color = "DESeq2\nFDR level", fill= "DESeq2\nFDR level")

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
df <- readRDS(file.path(resdir, "odds.table.rds")) %>% 
  filter(fdr.de == 1 )
df$count.on <- as.factor(df$count.on)
df$fdr.dmr <- as.numeric(as.character(df$fdr.dmr))
ggplot(df,
       aes(x = fdr.dmr, y = prop, group = count.on, color = count.on, 
           fill = count.on)) +
  #scale_y_continuous(trans= "log2") +
  #geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_ribbon(aes(ymin = p.ci.lower, ymax = p.ci.upper), alpha = 0.3, linetype=0) +
  geom_point() + 
  geom_line() +
  viridis::scale_color_viridis(discrete = TRUE, direction = -1) +
  viridis::scale_fill_viridis(discrete = TRUE, direction = -1) +
  scale_x_continuous(trans=nl_trans()) +
  theme_bw() +
  xlab("dmrseq FDR level") +
  ylab("Probability of being repressed (log2 fold change < 0)") +
  labs(color = "On genes\nthreshold", fill= "On genes\nthreshold") 
ggsave("../plots/figS1.pdf", width = 5.5, height = 4.25)


### FIGURE S2 - example regions
# created and saved in `mCG-RNAseq-analysis.Rmd`
ggsave("../plots/figS2.pdf", width = 10.5, height = 3.5)

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
df <- readRDS(file.path(resdir, "odds.table.cgi.rds")) %>% 
  filter(count.on == 50) %>% 
  filter(fdr.de %in% c(0.01, 0.05, 0.10, 0.5, 1) )
df$fdr.dmr <- as.numeric(as.character(df$fdr.dmr))
df$fdr.de[df$fdr.de == 1] <- "None (all genes)"
df$fdr.de <- factor(df$fdr.de)
df <- df %>% 
  mutate(cgi = ifelse(cgi, "CG Island", "non-CG Insland")) 

ggplot(df,
       aes(x = fdr.dmr, y = prop, group = fdr.de, color = fdr.de, 
           fill = fdr.de)) +
  #scale_y_continuous(trans= "log2") +
  #geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_ribbon(aes(ymin = p.ci.lower, ymax = p.ci.upper), alpha = 0.3, linetype=0) +
  geom_point() + 
  geom_line() +
  viridis::scale_color_viridis(discrete = TRUE, direction = -1) +
  viridis::scale_fill_viridis(discrete = TRUE, direction = -1) +
  scale_x_continuous(trans=nl_trans()) +
  theme_bw() +
  xlab("dmrseq FDR level") +
  ylab("Probability of being significantly repressed") +
  labs(color = "DESeq2\nFDR level", fill= "DESeq2\nFDR level") +
  facet_grid(~ cgi)

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


### FIGURE S8 - figure 2 by number of genes
df <- readRDS(file.path(resdir, "odds.table.rds")) %>%
  filter(fdr.de %in% c(0.01, 0.05, 0.10, 0.5, 1) )
df$fdr.dmr <- as.numeric(as.character(df$fdr.dmr))
df$fdr.de[df$fdr.de == 1] <- "None (all genes)"
df$fdr.de <- factor(df$fdr.de)

ggplot(df,
       aes(x = fdr.dmr, y = o.n, group = fdr.de, color = fdr.de, 
           fill = fdr.de)) +
  #scale_y_continuous(trans= "log2") +
  #geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_ribbon(aes(ymin = p.ci.lower, ymax = p.ci.upper), alpha = 0.3, linetype=0) +
  geom_point() + 
  geom_line() +
  viridis::scale_color_viridis(discrete = TRUE, direction = -1) +
  viridis::scale_fill_viridis(discrete = TRUE, direction = -1) +
  scale_x_continuous(trans=nl_trans()) +
  theme_bw() +
  xlab("dmrseq FDR level") +
  ylab("Number of genes") +
  labs(color = "DESeq2\nFDR level", fill= "DESeq2\nFDR level") +
  facet_wrap(~ count.on, ncol = 3)

ggsave("../plots/figS8.pdf", width = 10.5, height = 3.5)

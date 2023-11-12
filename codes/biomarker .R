
# read the saved phyloseq object
pooled <- readRDS(file = "data/R_objects/pooledN_soil_phyloseq.rds")

library("coin")
library("stringr")
library("ggtree") # BiocManager::install("ggtree")
library("ggtreeExtra") #install.packages("ggExtra")
library('MicrobiotaProcess')

set.seed(1024)
deres <- diff_analysis(obj = pooledN_soil, 
                       classgroup = "Status",
                       mlfun = "lda",
                       ratio=0.70,
                       filtermod = "pvalue",
                       firstcomfun = "kruskal.test",
                       padjust = "fdr",
                       firstalpha = 0.05,
                       strictmod = TRUE,
                       secondcomfun = "wilcox_test",
                       clmin = 10,
                       subclmin = 5,
                       subclwilc = TRUE,
                       ldascore = 4,
                       secondalpha = 0.05,
                       lda = 4,
                       normalization = 1e+08,
                       bootnums=999,
                       ci=0.95,
                       type = "species")
deres
str(deres)

diff_ana_results <- deres@result

str(diff_ana_results)


diff_ana_results[c('order', 'Species')] <- str_split_fixed(diff_ana_results$f, '__', 2)

diff_ana_results_S <- subset(diff_ana_results, order == 's')

ggplot(diff_ana_results_S,
       aes(x=LDAmean, y=reorder(Species, LDAmean)))+
  geom_segment(aes(yend=Species), xend=0, colour= "grey", linetype="dashed")+
  geom_errorbarh(aes(xmin=diff_ana_results_S$LDAlower, 
                      xmax=diff_ana_results_S$LDAupper,
                      height = .2))+
  geom_point(size=3, shape=21, aes(fill=Status))+
  scale_fill_manual(values=c("#990038", 
                             "#01AED9",
                             "#999999"))+
  theme_bw()+
  theme(panel.grid.major.y = element_blank())+
  theme(legend.position="top")+
  xlab(bquote(Log[10](LDA)))+
  ylab("Indicator species")


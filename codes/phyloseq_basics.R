###### packages ######

library("readxl")      
library("tibble")
library("vegan")
library("DESeq2") # BiocManager::install("DESeq2")
library("speedyseq") # remotes::install_github("mikemc/speedyseq") 


library("patchwork")


ape
cowplot
futile.logger
ggplot2
ggpubr
ggthemes
ggtree
ggtreeExtra 
MicrobiotaProcess
permute
phyloseq
vegan  
VennDiagram



##### create phyloseq object fom qiime2 outputs #####
library("qiime2R") # devtools::install_github("jbisanz/qiime2R")

pooled<-qza_to_phyloseq(
  features="data/from_QIIME2/table.qza",
  tree="data/from_QIIME2/rooted-tree.qza",
  taxonomy="data/from_QIIME2/Blast_taxonomy.qza",
  metadata = "data/oom_fir_metadata.tsv"
)
pooled

# Because P. ramorum was used as a PCR positive control
# ASVs associated with P. ramorum were removed from the dataset
library("phyloseq")

badTaxa = c("0789711810c5ac7ca64435418f6365fb",
            "12de6547fe64fd5c1f6322d6a05d5227",
            "161ff56a7342614f6ef4b8edf0a842db",
            "16add5903a825067d9ff71ca6fe38c5d",
            "17b6c8b48385f2a636be7740b8b511f4",
            "2c80d4b1b2029a973ba049a285ccfda1",
            "3295670d82066ee3a561e94073a7edc5",
            "3abb831e13ad2705665fcc47aff54887",
            "4ed05997f8e24927ed9e11bcb033d8bf",
            "4fb40bedcea0ceed7ccf0558244ccf89",
            "4fd790b2f8577717ed3ab4d05d96cc80",
            "5119b64c0df5a69d4816db6c02dea414",
            "65511ffe4fb30d148be5a536f2f94245",
            "676940d284db09340a2ca892707a3535",
            "6fe439d356df2d3999cf92353a1c06de",
            "74bfc4a34d427c08ed09913f55688db5",
            "75d11632dd44fef7a3aecaaa5bd46fec",
            "7b11f7761104a66e79c36c2ee5a53908",
            "928a575752e68f7a6f1718e562c064ad",
            "92c733f04e60297354ac09cff9c0a9f5",
            "9b6d7b7fe6df6a1ddcd31f51f3d32b6b",
            "a0338fd9a4d185a6bd8002606c610554",
            "a3893768727a77b236d2e94f3ff1e58d",
            "c253e538354f29f62863d8a4c0565ede",
            "c491f8a59eaadc5b3f91957811f04ae9",
            "cc2ee63e92bb2d8869626d4ddbf5e8ba",
            "e3a03b0bc3007108138f5422fc2ca183",
            "e7e247458bf3b54c57e19fcb005bef93",
            "ff9a1f5dc04b12ace376d1b33ff1cd2d"
            )
goodTaxa <- setdiff(taxa_names(pooled), badTaxa)
pooled <- prune_taxa(goodTaxa, pooled)

pooled

# verify files 
sample_names(pooled)
rank_names(pooled)
sample_variables(pooled)

#save the basic Phyloseq object
saveRDS(pooled, file = "data/R_objects/pooled_phyloseq.rds")

######### Summary stats ########

library(microbiome) # BiocManager::install("microbiome")
library(microbiomeutilities) #remotes::install_github("microsud/microbiomeutilities")

summarize_phyloseq(pooled)

Dep1<-plot_read_distribution(pooled, groups = "Status", 
                       plot.type = "histogram")+
                       theme_biome_utils()+
                       scale_x_continuous(trans='log10',limits=c(1, 1000000))+
                       scale_fill_manual(values=c("#111111"))+
                       geom_vline(xintercept = 400, colour = "black", linetype="dashed")+
                       theme(legend.position="none")+
  labs(x = "", y = "Count")

pooled <- prune_samples(sample_sums(pooled) >= 400, pooled)

summarize_phyloseq(pooled)

Dep2<-plot_read_distribution(pooled, groups = "Status", 
                      plot.type = "histogram")+
                      theme_biome_utils()+
                      scale_x_continuous(trans='log10', limits=c(1, 1000000))+
                      scale_fill_manual(values=c("#111111"))+ 
                      geom_vline(xintercept = 400, colour = "black", linetype="dashed")+
                      theme(legend.position="none")+
  labs(x = "Reads per samples", y = "Count")


library("cowplot")
depth<-plot_grid(Dep1+theme(legend.position="none"),
                  Dep2+theme(legend.position="none"), 
                  align="vh",
                  labels = c("A", "B"),
                  hjust = -1,
                  vjust= 2,
                  nrow = 2)

depth_final<-plot_grid(depth, ncol = 1, rel_heights = c(0.8, .05))
depth_final

library("ggpubr")
ggsave(file="figures/FigS2_depth_final.pdf", 
       width=8, height=5, units="in", dpi=300)


# Standardize number of reads in each sample using median sequencing depth ####
total = median(sample_sums(pooled))
standf = function(x, t=total) round(t * (x / sum(x)))
pooledN = transform_sample_counts(pooled, standf)

summarize_phyloseq(pooledN)

# bar_plot of relative abundance by sample type ####
library("ggtree") # BiocManager::install("ggtree")
library("ggtreeExtra") #install.packages("ggExtra")
library('MicrobiotaProcess') # BiocManager::install("MicrobiotaProcess")

library(ggthemes)

  classtaxa_g <- get_taxadf(obj=pooled, taxlevel=6)

  fclass_genus <- ggbartax(obj=classtaxa_g, 
                           facetNames="Type", 
                           plotgroup=TRUE, 
                           topn=15) +
    xlab(NULL) +
    ylab("Relative abundance (%)") +
    guides(fill= guide_legend(keywidth = 0.5, 
                              keyheight = 0.5, 
                              ncol=4))+
    theme(legend.text=element_text(size=8))+
    theme(axis.text.x=element_text(size=rel(1.2), angle=0,
                                   hjust = 0.5))+
    theme(axis.text.y=element_text(size=rel(1.2)))
  
  fclass_genus
  
  ggsave(file="figures/Fig2_rel_abund.pdf", 
         width=6, height=4, units="in", dpi=300)
  
# you may need to detach microbiotaprocess 
  detach("package:MicrobiotaProcess")
  
#### plot tree using the tax_glom function to merge ASVs with same taxon ####
  
# Subset the soil samples  
pooledN_soil <- subset_samples(pooledN, Type =="Soil")
pooledN_soil

saveRDS(pooledN_soil, file = "data/R_objects/pooledN_soil_phyloseq.rds")
  
### plot trees for Pythicaea ####
# subset the Pythicaea for soil samples
  
Soil_Pythiacae <- subset_taxa(pooledN_soil, Family =="Pythiaceae")
rank_names(Soil_Pythiacae)

A1<-plot_tree(tax_glom(Soil_Pythiacae, 
                       taxrank="Species"),
              method = "sampledodge",
              ladderize="left",
              nodelabf=nodeplotblank, 
              color="Status", 
              label.tips="Species", 
              text.size=3, 
              base.spacing=0.01,
              justify="jagged",
              shape="Status",
              size="abundance",
              plot.margin =0.9)+
  scale_size_continuous(range = c(0.0001, 4)) +
  scale_color_manual(values=c("#990038", 
                              "#01AED9",
                              "#999999"))+
  facet_wrap(~Type, scales="free_x")
A1

# Subset the Baiting samples  
pooledN_Baiting <- subset_samples(pooledN, Type =="Baiting")
pooledN_Baiting

# subset the Pythicaea for baiting samples
Bait_Pythiacae <- subset_taxa(pooledN_Baiting, Family =="Pythiaceae")
rank_names(Bait_Pythiacae)

B1<-plot_tree(tax_glom(Bait_Pythiacae, 
                   taxrank="Species"),
          method = "sampledodge",
          ladderize="left",
          nodelabf=nodeplotblank, 
          color="Status", 
          label.tips="Species", 
          text.size=3, 
          base.spacing=0.01,
          justify="jagged",
          shape="Status",
          size="abundance",
          plot.margin =0.9)+
  scale_size_continuous(range = c(0.0001, 4)) +
  scale_color_manual(values=c("#990038", 
                              "#01AED9",
                              "#999999"))+
  facet_wrap(~Type, scales="free_x")
B1

# Subset the Root samples  
pooledN_Root <- subset_samples(pooledN, Type =="Root")
pooledN_Root

# subset the Pythicaea for Root samples
Root_Pythiacae <- subset_taxa(pooledN_Root, Family =="Pythiaceae")
rank_names(Root_Pythiacae)

C1<-plot_tree(tax_glom(Root_Pythiacae, 
                       taxrank="Species"),
              method = "sampledodge",
              ladderize="left",
              nodelabf=nodeplotblank, 
              color="Status", 
              label.tips="Species", 
              text.size=3, 
              base.spacing=0.01,
              justify="jagged",
              shape="Status",
              size="abundance",
              plot.margin =0.9)+
  scale_size_continuous(range = c(0.0001, 4)) +
  scale_color_manual(values=c("#990038", 
                              "#01AED9",
                              "#999999"))+
  facet_wrap(~Type, scales="free_x")
C1

legend_b <- get_legend(
  C1+ guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom"))

library("cowplot")
Tree_1<-plot_grid(A1+theme(legend.position="none"),
          B1+theme(legend.position="none"), 
          C1+theme(legend.position="none"), 
          align="vh",
          labels = c("A", "B", "C"),
          hjust = -1,
          vjust= 2,
          nrow = 1)

plot_grid(Tree_1, legend_b, ncol = 1, rel_heights = c(0.8, .05))

ggsave(file="figures/fig4_tree_Pythiacae.pdf", 
       width=18, height=9, units="in", dpi=300)


#### Plot tree for everything but Pythiacae ####  

Soil_other <- subset_taxa(pooledN_soil, Family !="Pythiaceae")

A2<-plot_tree(tax_glom(Soil_other, 
                       taxrank="Species"),
              method = "sampledodge",
              ladderize="left",
              nodelabf=nodeplotblank, 
              color="Status", 
              label.tips="Species", 
              text.size=3, 
              base.spacing=0.01,
              justify="jagged",
              shape="Status",
              size="abundance",
              plot.margin =0.9)+
  scale_size_continuous(range = c(0.0001, 4)) +
  scale_color_manual(values=c("#990038", 
                                       "#01AED9",
                                       "#999999"))+
                                         facet_wrap(~Type, scales="free_x")
A2

Bait_other <- subset_taxa(pooledN_Baiting, Family !="Pythiaceae")

B2<-plot_tree(tax_glom(Bait_other, 
                       taxrank="Species"),
              method = "sampledodge",
              ladderize="left",
              nodelabf=nodeplotblank, 
              color="Status", 
              label.tips="Species", 
              text.size=3, 
              base.spacing=0.01,
              justify="jagged",
              shape="Status",
              size="abundance",
              plot.margin =0.9)+
  scale_size_continuous(range = c(0.0001, 4)) +
  scale_color_manual(values=c("#990038", 
                                       "#01AED9",
                                       "#999999"))+
                                         facet_wrap(~Type, scales="free_x")
B2


Root_other <- subset_taxa(pooledN_Root, Family !="Pythiaceae")

C2<-plot_tree(tax_glom(Root_other, 
                       taxrank="Species"),
              method = "sampledodge",
              ladderize="left",
              nodelabf=nodeplotblank, 
              color="Status", 
              label.tips="Species", 
              text.size=3, 
              base.spacing=0.01,
              justify="jagged",
              shape="Status",
              size="abundance",
              plot.margin =0.9)+
  scale_size_continuous(range = c(0.0001, 4)) +
  scale_color_manual(values=c("#990038", 
                                       "#01AED9",
                                       "#999999"))+
                                         facet_wrap(~Type, scales="free_x")
C2

legend_b <- get_legend(
  C2+ guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom"))

Tree_1<-plot_grid(A2+theme(legend.position="none"),
                  B2+theme(legend.position="none"), 
                  C2+theme(legend.position="none"), 
                  align="vh",
                  labels = c("A", "B", "C"),
                  hjust = -1,
                  vjust= 2,
                  nrow = 1)


plot_grid(Tree_1, legend_b, ncol = 1, rel_heights = c(0.8, .05))

ggsave(file="figures/fig5_tree_NO-Pythiacae.pdf", 
       width=18, height=9, units="in", dpi=300)


# Alpha diversity for soil samples -> microbiota process    #####

alphaobj <- get_alphaindex(pooledN_soil)
alphaobj2<-as.data.frame(alphaobj)
head(alphaobj2)

# ##### wilcoxon_test ######

pairwise.wilcox.test(alphaobj2$Observe, 
                     alphaobj2$Status, 
                     p.adjust.method = "hochberg")

pairwise.wilcox.test(alphaobj2$Pielou, 
                     alphaobj2$Status, 
                     p.adjust.method = "hochberg")

pairwise.wilcox.test(alphaobj2$Shannon, 
                     alphaobj2$Status, 
                     p.adjust.method = "hochberg")
#________________________________________________________

pairwise.wilcox.test(alphaobj2$Observe, 
                     alphaobj2$Specie, 
                     p.adjust.method = "hochberg")

pairwise.wilcox.test(alphaobj2$Pielou, 
                     alphaobj2$Specie, 
                     p.adjust.method = "hochberg")

pairwise.wilcox.test(alphaobj2$Shannon, 
                     alphaobj2$Specie, 
                     p.adjust.method = "hochberg")
#________________________________________________________

pairwise.wilcox.test(alphaobj2$Observe, 
                     alphaobj2$Region, 
                     p.adjust.method = "hochberg")

pairwise.wilcox.test(alphaobj2$Pielou, 
                     alphaobj2$Region, 
                     p.adjust.method = "hochberg")

pairwise.wilcox.test(alphaobj2$Shannon, 
                     alphaobj2$Region, 
                     p.adjust.method = "hochberg")

# ordination #####

# If the phyloseq tree need to be 're-rooted':

library("ape")
ps_tree = phy_tree(pooledN_soil)
sprintf("Is tree/ binary: %s", is.binary(ps_tree))
phy_tree(pooledN_soil) = multi2di(ps_tree)
sprintf("Is tree binary: %s", is.binary(phy_tree(pooledN_soil)))

library('MicrobiotaProcess') 
pcoares <- get_pcoa(obj=pooledN_soil, 
                    distmethod="wunifrac", 
                    method="hellinger")

# Visulizing the result by Status #####
pcaplot1 <- ggordpoint(obj=pcoares, biplot=TRUE, speciesannot=FALSE,
                       factorNames=c("Status", "Specie"), 
                       ellipse=TRUE, 
                       ellipse_pro=0.7,
                       ellipse_linewd = 0.7,
                       ellipse_lty = 2,
                       topn=5)+
  scale_color_manual(values=c("#990038", 
                              "#01AED9",
                              "#999999")) +
  scale_fill_manual(values=c("#990038", 
                             "#01AED9",
                             "#999999"))+
  theme(legend.position="none")+ 
  theme(plot.title = element_blank())

# pc = c(1, 3) to show the first and third principal components.
pcaplot2 <- ggordpoint(obj=pcoares, pc=c(1, 3), biplot=TRUE, speciesannot=FALSE,
                       factorNames=c("Status", "Specie"), 
                       ellipse=TRUE, 
                       ellipse_pro=0.7, 
                       ellipse_linewd = 0.7,
                       ellipse_lty = 2,
                       topn=5) +
  scale_color_manual(values=c("#990038", 
                              "#01AED9",
                              "#999999")) +
  scale_fill_manual(values=c("#990038", 
                             "#01AED9",
                             "#999999"))+ 
  theme(plot.title = element_blank())+
  ylim(-0.25, 0.30)

pcaplot1 | pcaplot2

ggsave(file="figures/Fig8_PCOA.pdf", 
       width=8, height=3, units="in", dpi=300)

###### Visulizing the result by Region ####### 

pcaplot3 <- ggordpoint(obj=pcoares, biplot=TRUE, speciesannot=FALSE,
                       factorNames=c("Region", "Specie"), 
                       ellipse=TRUE, 
                       ellipse_pro=0.7,
                       ellipse_linewd = 0.7,
                       ellipse_lty = 2,
                       topn=5)+
  scale_color_manual(values=c("#003f5c", 
                              "#e18745",
                              "#7baa68")) +
  scale_fill_manual(values=c("#003f5c", 
                             "#e18745",
                             "#7baa68"))+
  theme(legend.position="none")+ 
  theme(plot.title = element_blank())

# pc = c(1, 3) to show the first and third principal components.
pcaplot4 <- ggordpoint(obj=pcoares, pc=c(1, 3), biplot=TRUE, speciesannot=FALSE,
                       factorNames=c("Region", "Specie"), 
                       ellipse=TRUE, 
                       ellipse_pro=0.7,
                       ellipse_linewd = 0.7,
                       ellipse_lty = 2,
                       topn=5) +
  scale_color_manual(values=c("#003f5c", 
                              "#e18745",
                              "#7baa68")) +
  scale_fill_manual(values=c("#003f5c", 
                             "#e18745",
                             "#7baa68"))+ 
  theme(plot.title = element_blank())

pcaplot3 | pcaplot4

ggsave(file="figures/FigS4_PCOA.pdf", 
       width=8, height=3, units="in", dpi=300)

####Permutational Multivariate Analysis of Variance with vegan ADONIS #### 
library(vegan)

distme <- get_dist(pooledN_soil, distmethod ="wunifrac", method="hellinger")
sampleda <- data.frame(sample_data(pooledN_soil), check.names=FALSE)
sampleda <- sampleda[match(colnames(as.matrix(distme)),rownames(sampleda)),,drop=FALSE]
sampleda$Status <- factor(sampleda$Status)
sampleda$Specie <- factor(sampleda$Specie)
sampleda$Region <- factor(sampleda$Region)
sampleda$Year <- factor(sampleda$Year)

str(sampleda)

set.seed(1024)
adonis_pooledN_soil <- adonis2(distme ~ Status+Region+Year+Specie+Status*Region+Status*Year+Status*Specie+Status*Year*Region,  
                       data=sampleda,
                       permutation=999,
                       parallel = 16)
adonis_pooledN_soil

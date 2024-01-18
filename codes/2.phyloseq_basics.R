
library("phyloseq")
library("tidyverse")

pooledN<-read_rds("data/R_objects/pooledN_phyloseq.rds")

###############################################################################
# bar_plot of relative abundance by sample type ####
pooledN_dat <- psmelt(pooledN)

str(pooledN_dat)

predefined_Genus = c("Globisporangium",
                     "Pythium",
                     "Saprolegnia",
                     "Phytophthora",
                     "Phytopythium",
                     "Elongisporangium",
                     "Leptolegnia",
                     "Peronospora",
                     "Achlya",
                     "Pythiopsis",
                     "Aphanomyces")

#################################### 

library("tidyverse")
library("ggtext")
library("dplyr")

pooledN_dat <- pooledN_dat %>%
  mutate(Genus = case_when(
    Genus %in% 
      predefined_Genus ~ Genus,  # Keep valid species as is
    TRUE ~ "Other"  # Replace other species with "other"
  ))

otu_rel_abund <- pooledN_dat %>%
  group_by(Type) %>%
  mutate(rel_abund = Abundance / sum(Abundance)) %>%
  ungroup() %>%
  select(-Abundance) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
               names_to="level",
               values_to="taxon")


###################################

Type_abund<-otu_rel_abund %>%
  filter(level=="Genus") %>%
  group_by(Type, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Type, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  mutate(taxon = factor(taxon, 
                        levels=c("Globisporangium",
                                 "Pythium",
                                 "Saprolegnia",
                                 "Phytophthora",
                                 "Phytopythium",
                                 "Elongisporangium",
                                 "Leptolegnia",
                                 "Peronospora",
                                 "Achlya",
                                 "Pythiopsis",
                                 "Aphanomyces", 
                                 "Other")))

summary <- Type_abund %>% 
  group_by(Type, taxon) %>%
  summarize(mean_rel_abund) %>% 
  pivot_wider(names_from = Type, values_from = mean_rel_abund);summary


get_cols <- function (n){
  col <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
           "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
           "#ccebc5", "#ffed6f")
  
  col2 <- c("#1f78b4", "#ffff33", "#c2a5cf", "#ff7f00", "#810f7c",
            "#a6cee3", "#006d2c", "#4d4d4d", "#8c510a", "#d73027",
            "#78c679", "#7f0000", "#41b6c4", "#e7298a", "#54278f")
  
  col3 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
            "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
            "#ffff99", "#b15928")
  colorRampPalette(col2)(n)
}

# plot the stacked bar chart 
PooledN_stacked<-ggplot(data=Type_abund, 
                        aes(x=Type, 
                            y=mean_rel_abund, 
                            fill=taxon)) +
  theme_bw()+
  geom_col(colour = "black", width=0.8, linewidth=0.1) +
  theme(legend.title=element_blank())+
  labs(x=NULL,
       y="Relative Abundance (%)") +
  theme(legend.text = element_text(face="italic"))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=get_cols(12))+
  theme(legend.position="bottom")+
  guides(fill= guide_legend(keywidth = 0.6, 
                            keyheight = 0.7, 
                            ncol=3))+
  theme(axis.text.x=element_text(size=rel(1.2), hjust = 0.5))+
  theme(axis.text.y=element_text(size=rel(1.2)));PooledN_stacked

  ggsave(file="figures/Fig2_rel_abund.pdf", 
         width=4.7, height=5, units="in", dpi=300)

#############################################################################  
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

library("cowplot")
legend_b <- get_legend(
  C1+ guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom"))

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

A12<-plot_tree(tax_glom(Soil_other, 
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
A12

Bait_other <- subset_taxa(pooledN_Baiting, Family !="Pythiaceae")

B12<-plot_tree(tax_glom(Bait_other, 
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
B12


Root_other <- subset_taxa(pooledN_Root, Family !="Pythiaceae")

C12<-plot_tree(tax_glom(Root_other, 
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
C12

legend_c <- get_legend(
  C12+ guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom"))

Tree_1<-plot_grid(A12+theme(legend.position="none"),
                  B12+theme(legend.position="none"), 
                  C12+theme(legend.position="none"), 
                  align="vh",
                  labels = c("A", "B", "C"),
                  hjust = -1,
                  vjust= 2,
                  nrow = 1)


plot_grid(Tree_1, legend_c, ncol = 1, rel_heights = c(0.8, .05))

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

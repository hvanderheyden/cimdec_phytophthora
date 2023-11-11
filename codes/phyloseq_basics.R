###### packages ######

library("phyloseq")
library("readxl")      
library("tibble")
library("vegan")
library("DESeq2") # BiocManager::install("DESeq2")
library("speedyseq") # remotes::install_github("mikemc/speedyseq") 
library("ape")
library("ggstar")
library("forcats")
library("patchwork")
library("ggpubr")
library("plotROC")
library("viridis")
library("cowplot")

setwd("C:/Users/vanderheydenh/OneDrive - AGR-AGR/Projets/2023/Biovigilance/Pooled_paper/Manuscript")

# ##### color legend:
# Region 
# Beauce: "#003f5c", 
# Estrie: "#e18745", 
# Québec: "#7baa68"

# #### Status 
# Diseased_plantation: "#990038"
# Healthy_forest: "#01AED9"
# Healthy_plantation: "#999999"

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
# we removed P. ramorum from the dataset

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

library('microeco') # install.packages("microeco")
library("file2meco") # install.packages("file2meco", repos = BiocManager::repositories())

HV_2023_04_06 <- phyloseq2meco(pooled)

saveRDS(HV_2023_04_06, file = "HV_2023_04_06.rds")
saveRDS(pooled, file = "HV_2023_04_06_PS.rds")

# verify files 
sample_names(pooled)
rank_names(pooled)
sample_variables(pooled)

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

depth<-plot_grid(Dep1+theme(legend.position="none"),
                  Dep2+theme(legend.position="none"), 
                  align="vh",
                  labels = c("A", "B"),
                  hjust = -1,
                  vjust= 2,
                  nrow = 2)

depth_final<-plot_grid(depth, ncol = 1, rel_heights = c(0.8, .05))
depth_final


ggsave(file="depth_final.pdf", width=8, height=5, units="in", dpi=300)


# Standardize number of reads in each sample using median sequencing depth ####
total = median(sample_sums(pooled))
standf = function(x, t=total) round(t * (x / sum(x)))
pooledN = transform_sample_counts(pooled, standf)

summarize_phyloseq(pooledN)

pooled_soil <- subset_samples(pooled, Type =="Soil")
pooled_soil

pooled_Root <- subset_samples(pooled, Type =="Root")
pooled_Root

pooled_Baiting <- subset_samples(pooled, Type =="Baiting")
pooled_Baiting



pooledN_soil <- subset_samples(pooledN, Type =="Soil")
pooledN_soil

pooledN_Root <- subset_samples(pooledN, Type =="Root")
pooledN_Root

pooledN_Baiting <- subset_samples(pooledN, Type =="Baiting")
pooledN_Baiting

########## Need to subset the PS_objects by Status for the network analysis ####

pooledN_soil_dis <- subset_samples(pooledN_soil, Status =="Diseased_plantation")
pooledN_soil_dis

pooledN_soil_healthP <- subset_samples(pooledN_soil, Status =="Healthy_plantation")
pooledN_soil_healthP

pooledN_soil_healthF <- subset_samples(pooledN_soil, Status =="Healthy_forest")
pooledN_soil_healthF


# bar_plot ####

test<-plot_bar(pooled, "Status", "Abundance", "Family") + 
  geom_bar(aes(fill=Family), stat="identity", position="stack")+ 
  facet_grid(.~Type, scales= "free_y")+
  scale_fill_brewer("dark2")

data_test<-(test$data)
str(data_test)
head(data_test)


library(pivottabler)

pt <- PivotTable$new()
pt$addData(data_test)
pt$addColumnDataGroups("Type")
pt$addColumnDataGroups("Status")
pt$addRowDataGroups("Order")
pt$addRowDataGroups("Species")
pt$defineCalculation(calculationName="TotalTrains", 
                     summariseExpression="sum(Abundance)")
pt$renderPivot()

# not used
library("ggtree") # BiocManager::install("ggtree")
library("ggtreeExtra") #install.packages("ggExtra")
library('MicrobiotaProcess') # BiocManager::install("MicrobiotaProcess")

# you may need to detach microbiotaprocess 
detach("package:MicrobiotaProcess")

library(ggthemes)

  classtaxa_g <- get_taxadf(obj=pooledN, taxlevel=6)

  fclass_genus <- ggbartax(obj=classtaxa_g, 
                           facetNames="Type", 
                           plotgroup=TRUE, 
                           topn=11) +
    xlab(NULL) +
    ylab("Relative abundance (%)") +
    guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5, ncol=4))
  
  fclass_genus
  
#### plot tree using the tax_glom function to merge ASVs with same taxon ####

### plot trees for Pythicaea ####
  
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

Tree_1<-plot_grid(A1+theme(legend.position="none"),
          B1+theme(legend.position="none"), 
          C1+theme(legend.position="none"), 
          align="vh",
          labels = c("A", "B", "C"),
          hjust = -1,
          vjust= 2,
          nrow = 1)
        

plot_grid(Tree_1, legend_b, ncol = 1, rel_heights = c(0.8, .05))

ggsave(file="tree_final.pdf", width=11, height=8, units="in", dpi=300) 


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

################################### MicrobiotaProcess ###################

library("ggtree") # BiocManager::install("ggtree")
library("ggtreeExtra") #install.packages("ggExtra")
library('MicrobiotaProcess') # BiocManager::install("MicrobiotaProcess")
library("tidytree") # install.packages("tidytree")

# you need to you can detach a package
detach("package:MicrobiotaProcess")

# Alpha diversity -> microbiota process    #####

alphaobj <- get_alphaindex(pooledN_soil)
alphaobj2<-as.data.frame(alphaobj)
head(alphaobj2)

write.csv(alphaobj2, "alpha_obj.csv")

#############Pivot_table #####

library(pivottabler)

pt2 <- PivotTable$new()
pt2$addData(alphaobj2)
pt2$addRowDataGroups("Region")
pt2$defineCalculation(calculationName="Mean", 
                      summariseExpression="mean(Observe)")
pt2$defineCalculation(calculationName="Median", 
                     summariseExpression="median(Observe)")
pt2$defineCalculation(calculationName="min", 
                     summariseExpression="min(Observe)")
pt2$defineCalculation(calculationName="SE", 
                      summariseExpression="max(Observe)")
pt2$renderPivot()


pt1 <- PivotTable$new()
pt1$addData(alphaobj2)
pt1$addRowDataGroups("Status")
pt1$defineCalculation(calculationName="Mean", 
                      summariseExpression="mean(Observe)")
pt1$defineCalculation(calculationName="Median", 
                      summariseExpression="median(Observe)")
pt1$defineCalculation(calculationName="min", 
                      summariseExpression="min(Observe)")
pt1$defineCalculation(calculationName="max", 
                      summariseExpression="max(Observe)")
pt1$defineCalculation(calculationName="count", 
                      summariseExpression="n()")

pt1$renderPivot()


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

##### violon_plots (Status) ######

Observe_status <- ggbox(alphaobj, 
                 geom="violin", 
                 factorNames="Status",   
                 compare = TRUE,
                 testmethod = "wilcox.test",
                 signifmap = TRUE,
                 indexNames="Observe")+ 
  theme(aspect.ratio = 0.5)+
  theme( legend.position="none")+
  theme(strip.background = element_rect(colour=NA, fill="grey"))+
  scale_fill_manual(values=c("#990038", 
                             "#01AED9",
                             "#999999",
                             "#1B9E77", 
                             "#000033", 
                             "#FD9347"))

Shannon_status<- ggbox(alphaobj, 
              geom="violin", 
              factorNames="Status",   
              compare = TRUE, 
              testmethod = "wilcox.test",
              signifmap = TRUE,
              indexNames="Shannon")+ 
  theme(legend.position="Bottom")+
  theme(aspect.ratio = 0.5)+
  theme(strip.background = element_rect(colour=NA, fill="grey"))+
  scale_fill_manual(values=c("#990038", 
                             "#01AED9",
                             "#999999",
                             "#1B9E77", 
                             "#000033", 
                             "#FD9347"))
                             

Pielou_status<- ggbox(alphaobj, 
                       geom="violin", 
                       factorNames="Status",   
                       compare = TRUE, 
                       testmethod = "wilcox.test",
                       signifmap = TRUE,
                       indexNames="Pielou")+ 
  theme(legend.position="Bottom")+
  theme(aspect.ratio = 0.5)+
  theme(strip.background = element_rect(colour=NA, fill="grey"))+
  scale_fill_manual(values=c("#990038", 
                                      "#01AED9",
                                      "#999999",
                                      "#1B9E77", 
                                      "#000033", 
                                      "#FD9347"))




plot_grid(Observe_status+theme(legend.position="none"),
          Shannon_status+theme(legend.position="none"),
          Pielou_status+theme(legend.position="none"),
                  align="vh",
                  labels = c("A", "B", 'C'),
                  hjust = -1,
                  vjust= 3.5,
                  ncol=1)


##### violon_plots (Region) ######

Observe_Region <- ggbox(alphaobj, 
                        geom="violin", 
                        factorNames="Region",   
                        compare = TRUE,
                        testmethod = "wilcox.test",
                        signifmap = TRUE,
                        indexNames="Observe")+ 
  theme(aspect.ratio = 0.5)+
  theme( legend.position="none")+
  theme(strip.background = element_rect(colour=NA, fill="grey"))+
  scale_fill_manual(values=c("#990038", 
                                      "#01AED9",
                                      "#999999",
                                      "#1B9E77", 
                                      "#000033", 
                                      "#FD9347"))
                                      
Shannon_Region<- ggbox(alphaobj, 
                       geom="violin", 
                       factorNames="Region",   
                       compare = TRUE, 
                       testmethod = "wilcox.test",
                       signifmap = TRUE,
                       indexNames="Shannon")+ 
  theme(legend.position="Bottom")+
  theme(aspect.ratio = 0.5)+
  theme(strip.background = element_rect(colour=NA, fill="grey"))+
  scale_fill_manual(values=c("#990038", 
                                      "#01AED9",
                                      "#999999",
                                      "#1B9E77", 
                                      "#000033", 
                                      "#FD9347"))
                                      

Pielou_Region<- ggbox(alphaobj, 
                       geom="violin", 
                       factorNames="Region",   
                       compare = TRUE, 
                       testmethod = "wilcox.test",
                       signifmap = TRUE,
                       indexNames="Pielou")+ 
  theme(legend.position="Bottom")+
  theme(aspect.ratio = 0.5)+
  theme(strip.background = element_rect(colour=NA, fill="grey"))+
  scale_fill_manual(values=c("#990038", 
                                      "#01AED9",
                                      "#999999",
                                      "#1B9E77", 
                                      "#000033", 
                                      "#FD9347"))
                                      



plot_grid(Observe_Region+theme(legend.position="none"),
          Shannon_Region+theme(legend.position="none"),
          Pielou_Region+theme(legend.position="none"),
          align="vh",
          labels = c("A", "B", 'C'),
          hjust = -1,
          vjust= 3.5,
          ncol=1)


# ordination #####

#method  = “total”, “max”, “frequency”, “normalize”, “range”, “rank”, “rrank”, “standardize”, “pa”, 
# “chi.square”, “hellinger”, “log”, “clr”, “rclr”, “alr”

#distmethod = "unifrac",  "wunifrac", "manhattan", "euclidean", "canberra", "bray", "kulczynsk"

# If the phyloseq tree need to be 're-rooted':

ps_tree = phy_tree(pooledN_soil)
sprintf("Is tree binary: %s", is.binary(ps_tree))
phy_tree(pooledN_soil) = multi2di(ps_tree)
sprintf("Is tree binary: %s", is.binary(phy_tree(pooledN_soil)))

pcoares <- get_pcoa(obj=pooledN_soil, distmethod="unifrac", method="hellinger")

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

##### Adonis for diseased plantations vs forest  ####
diseased_forest<- subset_samples(pooledN_soil, Status !="Healthy_plantation")

distmeDF <- get_dist(diseased_forest, distmethod ="wunifrac", method="hellinger")
sampledaDF <- data.frame(sample_data(diseased_forest), check.names=FALSE)
sampledaDF <- sampledaDF[match(colnames(as.matrix(distmeDF)),rownames(sampledaDF)),,drop=FALSE]
sampledaDF$Status <- factor(sampledaDF$Status)

set.seed(1024)
adonis_DF <- adonis(distmeDF ~ Status, 
                               data=sampledaDF, 
                               permutation=9999)

pairwise_p<-numeric()
pairwise_p["DF"] <- adonis_DF[["aov.tab"]][["Pr(>F)"]][1]


##### Adonis for healthy plantations vs forest  #####
healthy_forest<- subset_samples(pooledN_soil, Status !="Diseased_plantation")

distmeHF <- get_dist(healthy_forest, distmethod ="wunifrac", method="hellinger")
sampledaHF <- data.frame(sample_data(healthy_forest), check.names=FALSE)
sampledaHF <- sampledaHF[match(colnames(as.matrix(distmeHF)),rownames(sampledaHF)),,drop=FALSE]
sampledaHF$Status <- factor(sampledaHF$Status)

set.seed(1024)
adonis_HF <- adonis(distmeHF ~ Status, 
                               data=sampledaHF, 
                               permutation=9999)
adonis_HF

pairwise_p["HF"] <- adonis_HF[["aov.tab"]][["Pr(>F)"]][1]

##### Adonis for diseased vs healthy plantations ####
diseased_healthy<- subset_samples(pooledN_soil, Status !="Healthy_forest")

distmeDH <- get_dist(diseased_healthy, distmethod ="wunifrac", method="hellinger")
sampledaDH <- data.frame(sample_data(diseased_healthy), check.names=FALSE)
sampledaDH <- sampledaDH[match(colnames(as.matrix(distmeDH)),rownames(sampledaDH)),,drop=FALSE]
sampledaDH$Status <- factor(sampledaDH$Status)

set.seed(1024)
adonis_DH <- adonis2(distmeDH ~ Status*Specie*Region, 
                    data=sampledaDH, 
                    permutation=9999)

adonis_DH

pairwise_p["D_H"] <- adonis_DH[["aov.tab"]][["Pr(>F)"]][1]

p.adjust(pairwise_p, method="hochberg")
















library('microeco') # install.packages("microeco")
library(file2meco) # install.packages("file2meco", repos = BiocManager::repositories())

mecoP_D <- phyloseq2meco(pooledN_soil)

t1 <- trans_abund$new(dataset = mecoP_D, taxrank = "Species", ntaxa = 60)
t1$plot_heatmap(facet = "Type", xtext_keep = FALSE, withmargin = FALSE)

mecoP_D_M <- mecoP_D$merge_samples(use_group = "Type")
# dataset1 is a new microtable object
# create trans_venn object
t1 <- trans_venn$new(mecoP_D_M, 
                     ratio = NULL)
venn<-t1$plot_venn()
venn

######### 

############# spiting diseased and healthy samples #######

pooledN_diseased <- subset_samples(pooledN, 
                                   Status_3 =="Diseased")
pooledN_diseased

pooledN_H <- subset_samples(pooledN, 
                            Status_3 =="Healthy")
pooledN_H

pooledN_diseased_phyto <- subset_taxa(pooledN_diseased, 
                                      Genus =="Phytophthora")

pooledN_H_phyto <- subset_taxa(pooledN_H, 
                               Genus =="Phytophthora")

pooledN_H_phyto <- subset_samples(pooledN_H_phyto, 
                               Type !="Root")



glom4 <- tax_glom(pooledN, 
                  taxrank = 'Species')

View(glom4@tax_table@.Data)

#Un diagramme de Venn a été utilisé pour répartir la distribution des Phytophthora spp. entre les différents types de sol.

percentages4 <- psmelt(glom4)
View(percentages4)
str(percentages4)
write.table(as.data.frame(percentages4), file = "percentage4.txt", sep = "\t")

library(VennDiagram)
vennlist <- get_vennlist(obj=pooledN_H_phyto, factorNames="Type")

vennp <- venn.diagram(vennlist,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c("#336600", "#663366"), #333366
                      cat.col=c("#336600", "#663366"), #333366
                      alpha = 0.45, 
                      resolution = 600,
                      fontfamily = "serif",
                      fontface = "bold",
                      cex = 0,
                      cat.cex = 1.3,
                      cat.default.pos = "outer",
                      cat.dist=0.1,
                      margin = 0.1, 
                      lwd = 3,
                      lty ='dotted',
                      imagetype = "svg")
grid::grid.draw(vennp)

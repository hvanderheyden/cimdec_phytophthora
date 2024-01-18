
##### create phyloseq object fom qiime2 outputs #####
library("qiime2R") # devtools::install_github("jbisanz/qiime2R")
library("phyloseq")

pooled<-qza_to_phyloseq(
  features="data/from_QIIME2/table.qza",
  tree="data/from_QIIME2/rooted-tree.qza",
  taxonomy="data/from_QIIME2/jan_2024_blast_taxonomy.qza",
  metadata = "data/oom_fir_metadata.tsv"
)
pooled

## Print a list of unique taxa 
get_taxa_unique(pooled, "Species")

# Clean the tax_table 
# remove the prefixes (k__, p__, etc.) 

library("tidyverse")
taxM <- data.frame(phyloseq::tax_table(pooled))
tax.cleanM <- data.frame(row.names = row.names(taxM),
                         Kingdom = str_replace(taxM[,1], "k__",""),
                         Phylum = str_replace(taxM[,2], "p__",""),
                         Class = str_replace(taxM[,3], "c__",""),
                         Order = str_replace(taxM[,4], "o__",""),
                         Family = str_replace(taxM[,5], "f__",""),
                         Genus = str_replace(taxM[,6], "g__",""),
                         Species = str_replace(taxM[,7], "s__",""),
                         stringsAsFactors = FALSE)
for (i in 1:7){ tax.cleanM[,i] <- as.character(tax.cleanM[,i])}

# remove the "_" and other corrections in species names
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "_", " "))

tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "x_", "x "))

tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "._", ". "))

# Actualize taxonomy for certain accessions  
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "CAL-2011b", "chlamydospora"))%>%
  mutate(Species=str_replace(Species, "sp. P1050", "sp. raspberry"))%>%
  mutate(Species=str_replace(Species, "aff. sansomeana", "sansomeana"))%>%
  mutate(Species=str_replace(Species, "sp. chlamydospora", "chlamydospora"))

# Concatenation of the species in clusters and 

####################################
###### Phytophthora species#########
####################################

# P. cactorum cluster (calde 1a)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "cactorum", "cactorum cluster"))%>%
  mutate(Species=str_replace(Species, "idaei", "cactorum cluster"))%>%
  mutate(Species=str_replace(Species, "alpina", "cactorum cluster"))%>%
  mutate(Species=str_replace(Species, "pseudotsugae", "cactorum cluster"))%>%
  mutate(Species=str_replace(Species, "aleatoria", "cactorum cluster"))

# P. citricola complex (clade 2c)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "citricola", "citricola complex"))%>%
  mutate(Species=str_replace(Species, "acerina", "citricola complex"))%>%
  mutate(Species=str_replace(Species, "plurivora", "citricola complex"))

# P. nemorosa cluster (clade 3a)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "nemorosa", "nemorosa cluster"))%>%
  mutate(Species=str_replace(Species, "ilicis", "nemorosa cluster"))%>%
  mutate(Species=str_replace(Species, "pseudosyringae", "nemorosa cluster"))%>%
  mutate(Species=str_replace(Species, "pluvialis", "nemorosa cluster"))

# P. cooljarloo cluster (clade 6a)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "aquae-cooljarloo", "cooljarloo cluster"))%>%
  mutate(Species=str_replace(Species, "rosacearum", "cooljarloo cluster"))%>%
  mutate(Species=str_replace(Species, "kwongonina", "cooljarloo cluster"))%>%
  mutate(Species=str_replace(Species, "pseudorosacearum", "cooljarloo cluster"))

# P. inundata cluster (clade 6a)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "inundata", "inundata cluster"))%>%
  mutate(Species=str_replace(Species, "humicola", "inundata cluster"))%>%
  mutate(Species=str_replace(Species, "condilina", "inundata cluster"))

# P. chlamydospora cluster (clade 6b)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "Phytophthora chlamydospora", "Phytophthora chlamydospora cluster"))%>%
  mutate(Species=str_replace(Species, "borealis", "chlamydospora cluster"))%>%
  mutate(Species=str_replace(Species, "mississippiae", "chlamydospora cluster"))%>%
  mutate(Species=str_replace(Species, "gonapodyides", "chlamydospora cluster"))

# P. megasperma cluster (clade 6b)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "Phytophthora megasperma", "Phytophthora megasperma cluster"))%>%
  mutate(Species=str_replace(Species, "Phytophthora crassamura", "Phytophthora megasperma cluster"))

# P. europaea cluster (clade 7a)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "europaea", "europaea cluster"))%>%
  mutate(Species=str_replace(Species, "abietivora", "europaea cluster"))%>%
  mutate(Species=str_replace(Species, "uliginosa", "europaea cluster"))%>%
  mutate(Species=str_replace(Species, "flexuosa", "europaea cluster"))%>%
  mutate(Species=str_replace(Species, "sp. cadmea", "europaea cluster"))

# P. pisi cluster (clade 7a)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "pisi", "pisi cluster"))%>%
  mutate(Species=str_replace(Species, "asiatica", "pisi cluster"))

# P. cryptogea complex (clade 8a)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "cryptogea", "cryptogea complex"))%>%
  mutate(Species=str_replace(Species, "pseudocryptogea", "cryptogea complex"))%>%
  mutate(Species=str_replace(Species, "erythroseptica", "cryptogea complex"))%>%
  mutate(Species=str_replace(Species, "kelmanii", "cryptogea complex"))

#####################################
###### Globisporangium species ######
#####################################

# G. attrantheridium complex (clade F)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "Globisporangium attrantheridium", "Globisporangium_attrantheridium complex"))%>%
  mutate(Species=str_replace(Species, "Pythium sp. balticum", "Globisporangium_attrantheridium complex"))

tax.cleanM<-tax.cleanM %>%
  mutate(Genus = if_else(Species == "Globisporangium_attrantheridium complex", "Globisporangium", Genus))

# G. irregulare complex (clade F)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "Globisporangium irregulare", "Globisporangium irregulare complex"))%>%
  mutate(Species=str_replace(Species, "Globisporangium cryptoirregulare", "Globisporangium irregulare complex"))%>%
  mutate(Species=str_replace(Species, "Globisporangium cylindrosporum", "Globisporangium irregulare complex"))

# G. Sylvaticum complex (clade F)
tax.cleanM<-tax.cleanM %>%
  mutate(Species=str_replace(Species, "Globisporangium sylvaticum", "Globisporangium sylvaticum complex"))%>%
  mutate(Species=str_replace(Species, "Globisporangium terrestre", "Globisporangium sylvaticum complex"))

tax.cleanM<-tax.cleanM %>%
  mutate(Genus = if_else(Species == "Globisporangium irregulare complex", "Globisporangium", Genus))

# G. heterothallicum complex (clade J)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "Globisporangium heterothallicum", "Globisporangium heterothallicum complex"))%>%
  mutate(Species=str_replace(Species, "Globisporangium glomeratum", "Globisporangium heterothallicum complex"))

# G. megalacanthum complex (clade J)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "Globisporangium megalacanthum", "Globisporangium megalacanthum complex"))%>%
  mutate(Species=str_replace(Species, "Globisporangium polymastum", "Globisporangium megalacanthum complex")) %>%
  mutate(Species=str_replace(Species, "Globisporangium buismaniae", "Globisporangium megalacanthum complex"))

#####################################
########## Pythium species ##########
#####################################

# P. arrhenomanes cluster (clade B1)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "arrhenomanes", "arrhenomanes cluster"))%>%
  mutate(Species=str_replace(Species, "aristoporum", "arrhenomanes cluster"))%>%
  mutate(Species=str_replace(Species, "phragmitis", "arrhenomanes cluster"))

# P. myriotylum complex (clade B1)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "myriotylum", "myriotylum complex"))%>%
  mutate(Species=str_replace(Species, "zingiberis", "myriotylum complex"))

# P. salpingophorum complex (clade B1)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "cf. salpingophorum/conidiophorum", "salpingophorum complex"))%>%
  mutate(Species=str_replace(Species, "conidiophorum", "salpingophorum complex"))

# P. dissotocum complex (clade B2)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "dissotocum", "salpingophorum complex"))%>%
  mutate(Species=str_replace(Species, "aff. dictyosporum", "salpingophorum complex"))%>%
  mutate(Species=str_replace(Species, "coloratum", "salpingophorum complex"))%>%
  mutate(Species=str_replace(Species, "diclinum", "salpingophorum complex"))%>%
  mutate(Species=str_replace(Species, "lutarium", "salpingophorum complex"))

# P. pachycaule complex (clade B2)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "pachycaule", "pachycaule complex"))%>%
  mutate(Species=str_replace(Species, "oopapillum", "pachycaule complex"))

# P. aquatile complex (clade B2)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "aquatile", "aquatile complex"))%>%
  mutate(Species=str_replace(Species, "sukuiense", "aquatile complex"))

# P. minus complex (clade E2)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "minus", "minus complex"))%>%
  mutate(Species=str_replace(Species, "aff. pleroticum", "minus complex"))


#####################################
###### Elongisporangium species #####
#####################################

# E. undulatum complex (clade H)
tax.cleanM<-tax.cleanM %>% 
  mutate(Species=str_replace(Species, "Elongisporangium undulatum", "Elongisporangium undulatum complex"))%>%
  mutate(Species=str_replace(Species, "Elongisporangium dimorphum", "Elongisporangium undulatum complex"))

# replace original tax_table with the clean tax_table 

phyloseq::tax_table(pooled) <- as.matrix(tax.cleanM) ; phyloseq::tax_table(pooled)

## Print a list of unique taxa 
get_taxa_unique(pooled, "Species")

# Because P. ramorum was used as a PCR positive control
# ASVs associated with P. ramorum were removed from the dataset


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

######################################################################################
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

#save the basic Phyloseq object
saveRDS(pooled, file = "data/R_objects/pooled_phyloseq.rds")

###############################################################################
# Standardize number of reads in each sample using median sequencing depth ####
total = median(sample_sums(pooled))
standf = function(x, t=total) round(t * (x / sum(x)))
pooledN = transform_sample_counts(pooled, standf)

summarize_phyloseq(pooledN)

#save the basic Phyloseq object
saveRDS(pooledN, file = "data/R_objects/pooledN_phyloseq.rds")



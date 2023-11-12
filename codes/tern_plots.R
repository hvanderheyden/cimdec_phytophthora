
library(microbiome)
library(microbiomeutilities)
library(ggtern) # install.packages("ggtern")
library("paletteer") # install.packages("paletteer")

pooledN_soilgg = subset_taxa(pooledN_soil, Species !="Unknown")

set.seed(1024)
tern_df_soil <- prep_ternary(pooledN_soilgg, group="Status", 
                        level= "Species", 
                        prev.thres=0.001)
head(tern_df_soil)

####
######
set.seed(1024)
ggtern(data = tern_df_soil, 
       aes(x=Diseased_plantation, 
           y=Healthy_plantation, 
           z=Healthy_forest)) +
  geom_point(aes(fill = Genus), 
             shape = 21,
             size = 2.5, 
             position= position_jitter_tern(x=0.01, 
                                            y=0.01, 
                                            z=0.01)) +
  theme_arrowcustomlength(0.05,
                          0.9) +   
  theme_ticksinside() +
  theme_showsecondary() +
  theme(legend.position = "bottom") + 
  theme(legend.title=element_blank()) +
  labs( x       = "D",
        xarrow  = "Diseased plantation",
        y       = "H",
        yarrow  = "Healthy plantation",
        z       = "F",
        zarrow  = "Healthy forest")+
  guides(color = guide_legend(override.aes = list(size = 10))) + 
  guides(fill=guide_legend(nrow=5))+
  facet_wrap(vars(Order), 
             ncol = 2, 
             nrow = 3) +
  theme(text = element_text(size = 10)) +
  theme(strip.text.x = element_text(size = 10)) +
  theme_showgrid_major() +
  theme_showgrid_minor() + 
  tern_limits(T=1.05, 
              L=1.05, 
              R=1.05) +
  scale_fill_paletteer_d("ggthemes::Classic_20")


ggsave(file="figures/FigS3_tern-plots.pdf", 
       width=5.5, height=8, units="in", dpi=300)

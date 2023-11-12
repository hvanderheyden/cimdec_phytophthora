![plot](https://github.com/hvanderheyden/cimdec_phytophthora/blob/main/figures/Fig1_Map.jpg?raw=true)
*Figure 1. Location of the sites sampled in this study*
# Oomycetes communities are influenced by land use and disease status in Christmas tree production in Southern Québec, Canada

#### *Hervé Van der Heyden, Marc-Olivier Duceppe, Guillaume Charron, Philippe Tanguay, and Guillaume J. Bilodeau*

#### Forests are threatened by many natural stressors intensified by climate change and anthropogenic activities, which tend to increase their susceptibility to pests and pathogens. Consequently, oomycetes-related forest decline or dieback cases are increasing in natural, urban, and agricultural landscapes. Christmas tree growers from Southern Québec, Canada, are experiencing root rot problems, with reported incidences up to 25%. In a previous study, seven Phytophthora spp. were associated with this root rot problem, but the overall diversity of oomycetes has not yet been investigated. Hence, in this study, we use a metabarcoding approach to provide an overview of the diversity, richness, and composition of the oomycetes community in fir plantations compared to surrounding natural forests.


![plot](https://github.com/hvanderheyden/cimdec_phytophthora/blob/main/figures/Graphical_abstract.png)
*Figure 2. Overview of the analytical pipeline used in this study*

## Description 
#### Here we provide data and R codes to perform the analysis presented in this study. 
#### 1. The pipeline used for taxonomic assignment was developed by DuceppeMO and can be found at https://github.com/duceppemo/QIIME2_ITS. 

#### The complete list of accessions used to train the QIIME2 classifier is provided in the data/accession_list.txt file, and the specific commands used are provided in the codes/QIIME2_pipeline_Phyto_paper.txt file. The obtained QIIME2 artifacts are provided in the data/from_QIIME2 folder.

#### 2. The R code used to import QIIME2 artifacts into phyloseq and perform analysis are provided in the codes/phyloseq_basics.R file. 
#### 3. The R code used to perform the LEFSe analysis is provided in the codes/biomarker.R.
#### 4. The R code used to perform the analysis of species cooccurence is provided in the codes/coocur_oom.R, whitch is used with the site_by_species.csv file.
#### 5. Finally, the R code used to generate the tern plots is provided in the codes/tern_plots.R file. 

## Citation 
### Upcoming 

## Funding 
#### This study was supported in part by a grant in the Cellule d'innovation des methodologies de diagnostic des ennemis de cultures from the Prime-Vert funding program of the Québec Ministry of Agriculture, Fisheries and Food.
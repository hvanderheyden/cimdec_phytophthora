library(sf) # the base package manipulating shapes
library(rgdal) # geo data abstraction library
library(geojsonio) # geo json input and output
library(spdplyr) # the `dplyr` counterpart for shapes
library(rmapshaper) # the package that allows geo shape transformation
library(magrittr) # data wrangling
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggspatial)
library(ggrepel)

setwd ("C:/Users/vanderheydenh/OneDrive - AGR-AGR/Projets/2023/Biovigilance/Pooled_paper")

#################### MAPPING ######### ############################
# transform the shpae files inot Geo JSON format to reduce the size of the file to plot 
# create a /data folder in the working directory (dsn = "data")
canada_raw = readOGR(dsn = "data", layer = "can", encoding = 'latin1') # 1
canada_raw_json <- geojson_json(canada_raw) # 2
canada_raw_sim <- ms_simplify(canada_raw_json) # 3
geojson_write(canada_raw_sim, file = "data/canada_cd_sim.geojson") # 4

#plotting the Canada MAP
#the geo JSON file create above should be in the /data folder within the wd 
canada_cd <- st_read("data/canada_cd_sim.geojson", quiet = TRUE) # 1
crs_string = "+proj=lcc +lat_1=49 +lat_2=77 +lon_0=-91.52 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs" # 2

# Define the maps' theme -- remove axes, ticks, borders, legends, etc.
theme_map <- function(base_size=9, base_family="") { # 3
  require(grid)
  theme_bw(base_size=base_size, base_family=base_family) %+replace%
    theme(axis.line=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid=element_blank(),
          panel.spacing=unit(0, "lines"),
          plot.background=element_blank(),
          legend.justification = c(0,0),
          legend.position = c(0,0)
    )
}

# Define the filling colors for each province; max allowed is 9 but good enough for the 13 provinces + territories
map_colors <- RColorBrewer::brewer.pal(9, "Pastel1") %>% rep(37) # 4

#Overlaying Custom Data
#using hardcoded coordinates, read in as a data frame with read.delim:
#2019####
conn_2019 = textConnection("Region	lat	long
Estrie	45.036552	-72.021969
Estrie	45.063416	-71.875822
Estrie	45.114804	-71.998998
Estrie	45.172889	-71.487605
Estrie	45.230306	-71.764462
Estrie	45.230506	-71.937717
Estrie	45.355868	-71.829743
Estrie	45.373306	-71.69727")

city_coords_2019 = read.delim(conn_2019, stringsAsFactors = F)

#Conforming the data frame to shape objects 
sf_cities_2019 = city_coords_2019 %>%
  select(long, lat) %>% # 1
  as.matrix() %>% # 2
  st_multipoint(dim = 'XY') %>% # 3
  st_sfc() %>% # 4
  st_set_crs(4269) # 5

#2020#####
conn_2020 = textConnection("Region	lat	long
Beauce	45.674133	-71.14798
Beauce	45.87638	-70.8572
Beauce	45.902924	-71.023241
Beauce	45.92416	-70.84666
Beauce	45.93305	-70.75379
Beauce	45.938861	-70.856917
Beauce	45.9544884	-70.8582454
Beauce	45.9544921	-70.8582454
Beauce	45.955	-70.82834
Beauce	45.96583	-70.79666
Beauce	45.991139	-71.485611
Beauce	46.01888	-70.87527
Beauce	46.02861	-70.87861
Beauce	46.056972	-71.124556
Beauce	46.05916	-70.80722
Beauce	46.068222	-70.929111
Beauce	46.171528	-70.848972")

city_coords_2020 = read.delim(conn_2020, stringsAsFactors = F)

#Conforming the data frame to shape objects 
sf_cities_2020 = city_coords_2020 %>%
  select(long, lat) %>% # 1
  as.matrix() %>% # 2
  st_multipoint(dim = 'XY') %>% # 3
  st_sfc() %>% # 4
  st_set_crs(4269) # 5

#2021####
conn_2021 = textConnection("Site	lat	long
Quebec	46.7132349	-71.2133818
Quebec	46.71446	-71.21265
Quebec	46.77147	-71.1087
Quebec	46.780222	-71.120361")

city_coords_2021 = read.delim(conn_2021, stringsAsFactors = F)

#Conforming the data frame to shape objects 
sf_cities_2021 = city_coords_2021 %>%
  select(long, lat) %>% # 1
  as.matrix() %>% # 2
  st_multipoint(dim = 'XY') %>% # 3
  st_sfc() %>% # 4
  st_set_crs(4269) # 5


villes <- data.frame(villes = c("Mt Katadin, Maine (USA)", "Sherbrooke", "Montréal", "Québec", "St-George de Beauce"), 
                     lat = c(45.941552027348216, 45.37674471, 45.40719381, 46.7802224369156, 46.129112569428784), 
                     lng = c(-68.91960321260669, -71.9284732729248,-73.93925767, -71.27733204, -70.67031424353718))

villes <- st_as_sf(villes, coords = c("lng", "lat"), remove = FALSE, 
                   crs = 4269, agr = "constant")
# Plot the maps####
ggplot() +
  geom_sf(aes(), color = "gray60", size = 0.3, data = canada_cd) +
  geom_sf(data = sf_cities_2019, aes(color= '#4A6C6F'), shape=20, size=5, alpha=0.9, show.legend ="point") +
  geom_sf(data = sf_cities_2020, aes(color = '#A84268'), shape=20, size=5, alpha=0.9, show.legend ="point") +
  geom_sf(data = sf_cities_2021, aes(color = '#FF5E5B'), shape=20, size=5, alpha=0.9, show.legend ="point")+
  theme(legend.position = c(0, 1))+
  geom_sf(data = villes)+
  geom_text_repel(data = villes, aes(x = lng, y = lat, label = villes), 
                  fontface = "bold", 
                  nudge_x = c(0.2, -0.3, -1, -1, 0.5), 
                  nudge_y = c(-0.45, 0.5, 0.25, 0.5, 0.4)) +
  coord_sf(xlim = c(-75, -65), ylim = c(45, 50), expand = TRUE) +
  scale_fill_manual(values = map_colors) +
  guides(fill = "none") +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                        style = north_arrow_nautical, 
                         height = unit(5, "cm"),
                         width = unit(5, "cm")) +
  theme_map() +
  scale_color_manual(values=c("#4A6C6F", "#FF5E5B", "#A84268"), 
                     labels=c("Estrie", "Beauce", "Québec"), name = "")+
  theme(legend.background = element_rect(fill='transparent'),
    panel.grid.major = element_line(color = "white"),
                legend.position = c(0.03, 0.81))+
  theme(text = element_text(size = 20))
                

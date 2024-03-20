############# Fisher et al. Figure 1 #############
#
# Figure 1. Map of the four collection sites in 
# Willapa Bay, WA, with the distribution of 
# carapace widths (in millimeters) at each site 
#
# M Fisher
#
#################################################


# Set up ------------------------------------------------------------------

# Stamen map tiles migrated to Stadia; follow this
#    to get the updated ggmap version needed 
#    to use 'get_stadiamaps' and 'register_stadiamaps'
#    First, Create a Stadia Maps account. I used an API key
# library(devtools)
# devtools::install_github("stadiamaps/ggmap") 


library(tidyverse)
library(here)
library(magrittr)
library(ggrepel)
library(cowplot)

# Data --------------------------------------------------------------------

# all crabs sequenced
meta.run <- read_csv(here('data','metadata','runs2-5_BF3_run_metadata.csv'))

# crabs with putative prey ID'd
dat <- read_csv(here('data','results','allRuns_BF3_filtered_FINAL_unique_taxa.csv'))

# site coordinates
site_coords <- read_csv(here('data/metadata/site_coords.csv')) %>%
  mutate(type=ifelse(type=="Sentinel","Slough",type))

# carapace widths
cw <- read_csv(here('data','metadata','Willapa Bay EGC Samples - Sample Data.csv'))

mapdat <- site_coords %>% dplyr::select(site_name,latitude,longitude,type) %>% distinct()

mapdat <- data.frame(site_name=c("Stackpole","Nahcotta","Long Beach","Oysterville"),
                     latitude=c(-124.03,-124.03,-124.017,-124.017),
                     longitude=c(46.598,46.492,46.434,46.554),
                     type=c("Slough","Clam Bed","Slough","Clam Bed"))

# DATA Site Map USGS National Map ----------------------------------------------

library(terrainr)
library(tmaptools)
library(sf)

# get northern part of bay
location2_of_interest <- tmaptools::geocode_OSM("Bay Center Washington")$coords

location2_of_interest <- data.frame(
  x = location2_of_interest[["x"]],
  y = location2_of_interest[["y"]]
)
location2_of_interest <- st_as_sf(
  location2_of_interest, 
  coords = c("x", "y"), 
  crs = 4326
)
output2_tiles <- get_tiles(set_bbox_side_length(location2_of_interest, 30000),
                           services = c("ortho"),
                           resolution = 70 # pixel side length in meters
)


# file.rename(from=output2_tiles$ortho,to=here('data','metadata','USGSNAIPPlus_ortho60m_2.tif'))

# get west part of bay
location3_of_interest <- tmaptools::geocode_OSM("Klipsan Beach Washington")$coords

location3_of_interest <- data.frame(
  x = location3_of_interest[["x"]],
  y = location3_of_interest[["y"]]
)
location3_of_interest <- st_as_sf(
  location3_of_interest, 
  coords = c("x", "y"), 
  crs = 4326
)
output3_tiles <- get_tiles(set_bbox_side_length(location3_of_interest, 30000),
                           services = c("ortho"),
                           resolution = 70 # pixel side length in meters
)
# file.rename(from=output3_tiles$ortho,to=here('data','metadata','USGSNAIPPlus_ortho60m_file1ebc3c826cfd.tif'))

# get east part of bay
location4_of_interest <- tmaptools::geocode_OSM("Naselle Washington")$coords

location4_of_interest <- data.frame(
  x = location4_of_interest[["x"]],
  y = location4_of_interest[["y"]]
)
location4_of_interest <- st_as_sf(
  location4_of_interest, 
  coords = c("x", "y"), 
  crs = 4326
)
output4_tiles <- get_tiles(set_bbox_side_length(location4_of_interest, 30000),
                           services = c("ortho"),
                           resolution = 70 # pixel side length in meters
)

location5_of_interest <- tmaptools::geocode_OSM("Raymond Washington")$coords

location5_of_interest <- data.frame(
  x = location5_of_interest[["x"]],
  y = location5_of_interest[["y"]]
)
location5_of_interest <- st_as_sf(
  location5_of_interest, 
  coords = c("x", "y"), 
  crs = 4326
)
output5_tiles <- get_tiles(set_bbox_side_length(location5_of_interest, 20000),
                           services = c("ortho"),
                           resolution = 70 # pixel side length in meters
)
# file.rename(from=output4_tiles$ortho,to=here('data','metadata','USGSNAIPPlus_ortho60m_file1ebc498433ea.tif'))

# WA state Map (1a inset) --------------------------------------------------

inset <- ggplot() +
  geom_polygon(data = map_data("county", region="Washington"), 
               aes(x = long, 
                   y = lat, 
                   group = group), fill="gray40",color="gray40") +
  geom_polygon(data = map_data("county", region="Oregon"), 
               aes(x = long, 
                   y = lat, 
                   group = group), fill="gray70",color="gray70") +
  geom_rect(aes(ymin=46.2, ymax=46.85,xmax=-123.7,xmin=-124.2), fill=NA,color="red",size=0.5) +
  coord_fixed(xlim=c(-124.7,-122),ylim=c(45.5,48.9)) +
  theme_void() + theme(panel.background=element_rect(fill="white",color="white"))

inset

# MAPPING Site Map USGS National Map (1a) ---------------------------------

library(terrainr)
library(sf)

# read back in map tiles
# output2_tiles<-terra::rast(here('data','metadata','USGSNAIPPlus_ortho_2.tif'))
# output3_tiles<-terra::rast(here('data','metadata','USGSNAIPPlus_ortho60m_file1ebc3c826cfd.tif'))
# output4_tiles<-terra::rast(here('data','metadata','USGSNAIPPlus_ortho60m_file1ebc498433ea.tif'))

mapdat_sf <- sf::st_as_sf(mapdat,coords = c("latitude","longitude"))
mapdat_sf <- sf::st_set_crs(mapdat_sf, 4326)

tiff(here('figs','Figure1a.tif'),res=300, height=2000,width=2700)
ggdraw(ggplot() +  
  geom_spatial_rgb(data = output2_tiles[["ortho"]],
                   aes(x = x, y = y, r = red, g = green, b = blue)) +
  geom_spatial_rgb(data = output3_tiles[["ortho"]],
                   aes(x = x, y = y, r = red, g = green, b = blue)) + 
  geom_spatial_rgb(data = output4_tiles[["ortho"]],
                   aes(x = x, y = y, r = red, g = green, b = blue))  +
    geom_spatial_rgb(data = output5_tiles[["ortho"]],
                     aes(x = x, y = y, r = red, g = green, b = blue))  +
  geom_rect(aes(ymin=46.375,ymax=46.56,xmin=-124.12,xmax=-124.0655), fill="#1D4054") +
  geom_sf(data=mapdat_sf, aes(color=type), size=6) + 
  geom_sf(data=mapdat_sf, pch=1,color="white", size=6) +
  ggrepel::geom_label_repel(
    data = head(mapdat_sf),
    aes(label = site_name, geometry = geometry),
    stat = "sf_coordinates", size=6,
    min.segment.length = 0,nudge_x=0.2,nudge_y=c(0.02,-0.02,-0.02,0.02),segment.colour = 'white'
  ) +
  ylab("Latitude") + xlab("Longitude") +
  scale_color_manual(values=c("#018571","#a6611a")) +
  coord_sf(xlim=c(-124.1, -123.75),ylim=c(46.3,46.79),crs = 4326) + 
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1,size=10),
                     axis.title=element_blank(),
                     axis.text.y=element_text(size=10),
                     legend.title=element_blank(),
                     legend.text=element_text(size=12), legend.position="top",
                     plot.margin=unit(c(t=0,r=-1,b=0,l=0), "cm"))
) + draw_plot(inset, x=0.6,y=0.1,width=0.1,height=0.13)
dev.off()




# Count crab at each step -------------------------------------------------
dissected <- site_coords %>% group_by(site_name) %>% summarise(n_dissected=sum(n_dissected)) %>%
  rename(Site=site_name)

sequenced <- meta.run %>%
  filter(!is.na(site_month)) %>%
  separate(site_month, into=c("Site","Month"), sep="-") %>%
  group_by(Site) %>%
  summarise(`sequenced`=length(unique(sample)))

success <- dat %>%
  filter(!is.na(site_month)) %>%
  separate(site_month, into=c("Site","Month"), sep="-") %>%
  group_by(Site) %>%
  summarise(`with prey`=length(unique(sample))) %>%
  mutate(`with prey`=ifelse(is.na(`with prey`),0,`with prey`))

ncrab <- full_join(dissected,sequenced,by="Site") %>%
  full_join(success,by="Site") %>%
  mutate(`seq with prey` = round(`with prey`/`sequenced`,2)) 

ncrab


# Carapace widths by step (1b) ---------------------------------------------

cw %<>%
  mutate(sequenced=ifelse(Sample_label %in% dat$sample, "Sequenced w/ Prey",
                          ifelse(Sample_label %in% meta.run$sample, "Sequenced","Trapped"))) %>%         # which crab were sequenced / with prey / just trapped?
  mutate(site=ifelse(grepl("Stackpole",Site_Name), "Stackpole",                                          # site infor for each crab
                           ifelse(grepl("Oysterville",Site_Name), "Oysterville",
                                  ifelse(grepl("Nahcotta",Site_Name),"Nahcotta",
                                         ifelse(grepl("Long Beach",Site_Name),"Long Beach", NA)))))
## calculate mean CW by site
cw %>% group_by(site,sequenced) %>% summarise(mean_cw=mean(CW_mm)) %>% pivot_wider(names_from=sequenced,values_from=mean_cw)




cw.prey <- cw %>% filter(sequenced=="Sequenced w/ Prey") %>% mutate(sequenced="with Prey Taxa")
cw.sequenced <- cw%>%filter(sequenced %in% c("Sequenced","Sequenced w/ Prey")) %>% mutate(sequenced="Sequenced")
plotdat <- cw %>% mutate(sequenced="Trapped") %>%
  bind_rows(cw.sequenced,cw.prey)

plotdat$site  <- factor(plotdat$site, levels=c("Stackpole","Oysterville","Nahcotta","Long Beach"))
plotdat$sequenced <- factor(plotdat$sequenced, levels=c("Trapped","Sequenced","with Prey Taxa"))

plot_CW <- ggplot(plotdat, aes(x=CW_mm, fill=sequenced)) +
  geom_histogram(position="dodge") +
  scale_fill_manual(values=c("grey60","black","firebrick3"),name="") +
  facet_grid(rows=vars(site), switch="y") + 
  labs(x="Carapace Width (mm)",y="Number of Green Crab") + 
  theme_bw() + theme(legend.position="top", strip.text=element_text(size=14),legend.text=element_text(size=12),
                     strip.background=element_rect(fill="white"),
                     axis.text.x=element_text(size=14), axis.title.x=element_text(size=14),
                     strip.placement="outside",
                     axis.title.y=element_blank(),
                     plot.margin=unit(c(0,0,0,0), "cm"))
plot_CW


# Final Figure ------------------------------------------------------------

tiff(here('figs','Figure1a.tif'),res=300, height=2000,width=2700)
ggplot() +  
  geom_spatial_rgb(data = output2_tiles[["ortho"]],
                   aes(x = x, y = y, r = red, g = green, b = blue)) +
  geom_spatial_rgb(data = output3_tiles[["ortho"]],
                   aes(x = x, y = y, r = red, g = green, b = blue)) + 
  geom_spatial_rgb(data = output4_tiles[["ortho"]],
                   aes(x = x, y = y, r = red, g = green, b = blue)) +
  ylab("Latitude") + xlab("Longitude") +
  coord_sf(xlim=c(-124.1, -123.8),ylim=c(46.3,46.8),crs = 4326) + 
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1,size=10),
                     axis.title=element_blank(),
                     axis.text.y=element_text(size=10),
                     legend.title=element_blank(),
                     legend.text=element_text(size=12), legend.position="top",
                     plot.margin=unit(c(t=0,r=-1,b=0,l=0), "cm"))
dev.off()

tiff(here('figs','Figure1a-2.tif'),res=300, height=2000,width=2700)
ggplot() +  
  geom_spatial_rgb(data = output2_tiles[["ortho"]],
                   aes(x = x, y = y, r = red, g = green, b = blue)) +
  geom_spatial_rgb(data = output3_tiles[["ortho"]],
                   aes(x = x, y = y, r = red, g = green, b = blue)) + 
  geom_spatial_rgb(data = output4_tiles[["ortho"]],
                   aes(x = x, y = y, r = red, g = green, b = blue)) +
  geom_point(data=mapdat, aes(x=longitude,y=latitude, fill=type), size=6) +
  geom_label_repel(data=mapdat, aes(x=longitude,y=latitude, label=site_name),
                   segment.colour = 'white',
                   nudge_x=0.22,nudge_y=0.01,
                   force=1, force_pull=0.2, segment.size=0.5, size=6) + 
  ylab("Latitude") + xlab("Longitude") +
  scale_color_manual(values=c("#018571","#a6611a")) +
  coord_sf(xlim=c(-124.1, -123.8),ylim=c(46.3,46.8),crs = 4326) + 
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1,size=10),
                     axis.title=element_blank(),
                     axis.text.y=element_text(size=10),
                     legend.title=element_blank(),
                     legend.text=element_text(size=12), legend.position="top",
                     plot.margin=unit(c(t=0,r=-1,b=0,l=0), "cm"))
dev.off()


png(here('figs','Figure1.png'),res=300, height=2000,width=2700)
plot_grid(sitemap,
           plot_CW,ncol=2, rel_widths=c(1,0.92))
dev.off()


tiff(here('figs','Figure1_rmaps.tif'),res=300, height=2000,width=2700)
plot_grid(sitemap,plot_CW,ncol=2, rel_widths=c(1,0.92))
dev.off()


# Figure S1 ---------------------------------------------------------------

cw.prey <- cw %>% filter(sequenced=="Sequenced w/ Prey")
cw.sequenced <- cw%>%filter(sequenced %in% c("Sequenced","Sequenced w/ Prey")) %>% mutate(sequenced="Sequenced")
plotdat <- cw %>% mutate(sequenced="Trapped") %>%
  bind_rows(cw.sequenced,cw.prey)

plotdat$site  <- factor(plotdat$site, levels=c("Stackpole","Oysterville","Nahcotta","Long Beach"))
plotdat$sequenced <- factor(plotdat$sequenced, levels=c("Trapped","Sequenced","Sequenced w/ Prey"))

figs1b <- ggplot(plotdat, aes(x=CW_mm, fill=sequenced)) +
  geom_histogram(position="dodge") +
  scale_fill_manual(values=c("grey60","black","firebrick3"),name="") +
  facet_grid(cols=vars(Sex)) + 
  labs(x="Carapace Width (mm)",y="Number of Green Crab") + theme_bw() + theme(legend.position="top")
figs1b

figs1a <- plotdat %>% group_by(sequenced,Sex) %>% summarise(ncrab=length(unique(Sample_label))) %>%
  group_by(sequenced) %>% mutate(totalCrab=sum(ncrab)) %>% ungroup() %>%
  mutate(pcrab=ncrab/totalCrab) %>%
  ggplot(aes(x=sequenced,y=pcrab,fill=Sex)) + 
  geom_col() +
  scale_fill_manual(values=c("darkolivegreen2","chartreuse4"), name="") +
  labs(x="",y="Proportion of Green Crab") +
  theme_bw() + theme(legend.position="top",axis.text.x=element_text(angle=45,hjust=1))
figs1a

png(here('figs','FigureS1_metadat.png'),res=300, height=1500,width=2100)
plot_grid(plotlist=list(figs1a,figs1b), rel_widths = c(0.4,1),labels=c("a","b"))
dev.off()




# Site Map  Stamen (1a) -----------------------------------------------------------

library(ggmap)

register_stadiamaps(skey, write = FALSE)
myLocation <- c(-124.5, 46, -123.5, 47.5)

## terrain
# mymap <- get_stadiamap(bbox=myLocation,maptype="stamen_terrain_background",crop=FALSE)
# saveRDS(mymap, file=here('data','stamen_map_terrain_wide_WAcoast.rds'))
mymap <- readRDS(file=here('data','stamen_map_terrain_wide_WAcoast.rds'))

## toner
# mymap <- get_stadiamap(bbox=myLocation,maptype="stamen_toner_background",crop=FALSE)
# saveRDS(mymap, file=here('data','stamen_map_toner_wide_WillapaBay.rds'))
mymap <- readRDS(here('data','stamen_map_toner_wide_WillapaBay.rds'))


attr_mymap <- attr(mymap, "bb")    # save attributes from original

# change background color in raster from black to dark blue
mymap[mymap == "#000000"] <- "#2f567f"

# correct class, attributes
class(mymap) <- c("ggmap", "raster")
attr(mymap, "bb") <- attr_mymap
ggmap(mymap)  # check map

sitemap <- ggmap(mymap) +
  geom_point(data=mapdat, aes(x=longitude,y=latitude, col=type), size=6) + 
  geom_label_repel(data=mapdat, aes(x=longitude,y=latitude, label=site_name),
                   segment.colour = 'black',
                   nudge_x=0.22,nudge_y=0.01,
                   force=1, force_pull=0.2, segment.size=0.5, size=6) +
  geom_text(aes(x=-123.98,y=46.64,label="italic('Willapa\n Bay')"), col="deepskyblue4", size=4, parse=TRUE) +
  geom_text(aes(x=-123.87,y=46.22,label="italic('Columbia River')"), col="deepskyblue4", size=4, parse=TRUE) +
  ylab("Latitude") + xlab("Longitude") +
  scale_color_manual(values=c("#018571","#a6611a")) +
  coord_fixed(xlim=c(-124.15,-123.7),ylim=c(46.15,46.9)) +
  theme(axis.title=element_blank(),
        axis.text=element_text(size=10),
        legend.title=element_blank(),
        legend.text=element_text(size=12), legend.position="top",
        plot.margin=unit(c(0,0,0,0), "cm"))
sitemap



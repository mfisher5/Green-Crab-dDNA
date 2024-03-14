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
library(ggplot2)
library(ggrepel)
library(ggmap)
library(janitor)
library(cowplot)




# Background map ----------------------------------------------------------

register_stadiamaps(skey, write = FALSE)

# myLocation <- c(-124.2, 46, -123.7, 47.4)
myLocation <- c(-124.5, 46, -123.5, 47.5)


## terrain
mymap <- get_stadiamap(bbox=myLocation,maptype="stamen_terrain_background",crop=FALSE)
saveRDS(mymap, file=here('data','stamen_map_terrain_wide_WAcoast.rds'))
mymap <- readRDS(file=here('data','stamen_map_terrain_wide_WAcoast.rds'))

## toner
mymap <- get_stadiamap(bbox=myLocation,maptype="stamen_toner_background",crop=FALSE)
saveRDS(mymap, file=here('data','stamen_map_toner_wide_WillapaBay.rds'))
mymap <- readRDS(here('data','stamen_map_toner_wide_WillapaBay.rds'))


attr_mymap <- attr(mymap, "bb")    # save attributes from original

# change background color in raster from black to dark blue
mymap[mymap == "#000000"] <- "#2f567f"

# correct class, attributes
class(mymap) <- c("ggmap", "raster")
attr(mymap, "bb") <- attr_mymap
ggmap(mymap)  # check map



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



# Site Map (1a) -----------------------------------------------------------

mapdat <- site_coords %>% dplyr::select(site_name,latitude,longitude,type) %>% distinct()

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

png(here('figs','Figure1.png'),res=300, height=2000,width=2700)
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






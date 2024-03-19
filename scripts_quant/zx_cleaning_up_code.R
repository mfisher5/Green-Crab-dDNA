
# reconcile changes to data list from moving code between z1,z2 for readability
# data - original
load(here('../','data','zoid',"qm_model_dataLIST_2023-8-29.rds"))


# make sure mocks have "Z marina" and samples have "Zosteraceae"
metabarcoding <- EGC_obs[[2]]
mock <- EGC_obs[[1]]

# make sure Mya is listed as Mya arenaria
mock %<>% mutate(species=ifelse(species=="Mya","Mya arenaria",species))
metabarcoding %<>% mutate(species=ifelse(species=="Mya","Mya arenaria",species))
unique(mock$species)
unique(metabarcoding$species)

# fix species --> genus from manual filter
manual.filter <- data.frame(species=c("Crangon sp. UF_53625","Elachista sp. 1fucicola BOLD-2018",
                                      "Scytosiphon sp. Mitopromiscuus","Arachnida sp. BOLD:ACM9770","Platynereis sp. CMC02",
                                      "Chordariaceae sp. 15 AP-2014","Ectocarpus sp. 1siliculosus",
                                      "Kurtiella bidentata","Neoporphyra haitanensis"),
                            taxon=c("Crangon","Elachista","Scytosiphon","Arachnida",
                                    "Platynereis","Chordariaceae","Ectocarpus",
                                    "Kurtiella","Neoporphyra"))
metabarcoding %<>% left_join(manual.filter) 
metabarcoding %<>% mutate(species=ifelse(is.na(taxon),species,taxon))
metabarcoding %<>% filter(species!="Appasus major")

unique(metabarcoding$species)[which(!(unique(metabarcoding$species) %in% sample.tax$taxon))]

# add in missing zeros
metabarcoding_reformat <- metabarcoding %>% 
  dplyr::select(-taxon) %>%
  mutate(Sample=as.character(Sample),station=as.character(station)) %>%
  unite(station, creek, time, tech, col = "Sample", remove = T, sep = "#") %>%
  pivot_wider(names_from = Sample, values_from = Nreads, values_fill = 0) %>% 
  pivot_longer(-c(species), names_to = "Sample", values_to = "Nreads") %>% 
  separate(Sample, into = c("station", "creek", "time", "tech"), sep = "#")



dim(mock)
mock.reformatted <- mock %>% 
  unite(Community, CommType, tech, col = "Sample", remove = T, sep = "#") %>% 
  dplyr::select(-totReads,-propReads) %>%
  pivot_wider(names_from = Sample, values_from = Nreads, values_fill = 0) %>%
  pivot_longer(-c(species, b_proportion, N_pcr_mock), names_to = "Sample", values_to = "Nreads") %>% 
  separate(Sample, into = c("Community", "CommType", "tech"), sep = "#") %>%
  left_join(mock%>%dplyr::select(Community,tech,totReads) %>% mutate(tech=as.character(tech)) %>% distinct())
mock.reformatted %<>% mutate(propReads=Nreads/totReads)
dim(mock.reformatted)

EGC_obs <- list(metabarcoding,mock)

saveRDS(EGC_obs,here('../','data','quant_modeling',"QM_dataLIST.rds"))


# Metabarcoding calibration model function
# Ryan P Kelly Aug 2022; Revised March 2023

################DEFINE SUB-FUNCTIONS

format_metabarcoding_data_forQM <- function(input_metabarcoding_data, input_mock_comm_data, N_pcr_mock,
                                      Level_1_treatment_envir, #e.g., unique sampling site
                                      Level_2_treatment_envir, #nested within level1 units, e.g., unique biological samples
                                      Level_3_treatment_envir,  #nested within level2 units, e.g., technical replicates of biological samples (NA if absent)
                                      Level_1_treatment_mock, #e.g., unique biol community
                                      Level_2_treatment_mock, #nested within level1 units, e.g., unique biological replicates
                                      Level_3_treatment_mock  #nested within level2 units, e.g., technical replicates of biological replicates (NA if absent)
){
  require(tidyverse)
  #require(rlang)
  
  # level1 <- rlang::enquo(Level_1_treatment) #e.g., unique sampling site
  # level2 <- rlang::enquo(Level_2_treatment) #nested within level1 units, e.g., unique biological samples
  # level3 <- rlang::enquo(Level_3_treatment) #nested within level2 units, e.g., technical replicates of biological samples
  
  Observation <- input_metabarcoding_data
  
  Mock <- input_mock_comm_data %>% 
    replace(is.na(.), 0) %>% 
    filter(b_proportion > 0) %>% 
    filter(Nreads > 0) %>% #omit things that are absent from the mocks
    filter(!is.na(Nreads)) %>% 
    rename(level1 = as.name(Level_1_treatment_mock),
           level2 = as.name(Level_2_treatment_mock))
  
  if (is.na(Level_3_treatment_mock)) {Mock$level3 <- 1} else {Mock <- Mock %>% rename(level3 = as.name(Level_3_treatment_mock))}  
  
  Mock <- Mock %>% 
    unite(c(level1, level2, level3), col = "Sample", sep = ".", remove = F)
  
  # only keep species present in both mocks and observations
  keepSpecies <- intersect(Mock$species, Observation$species)
  Observation <- Observation %>% 
    filter(species %in% keepSpecies) 
  Mock <- Mock %>% filter(species %in% keepSpecies)
  
  # index species to a common standard 
  sp_list <- data.frame(
    species = c(Mock$species, Observation$species) %>% unique(),
    species_idx = NA)
  sp_list$species_idx <- match(sp_list$species, unique(sp_list$species)) 
  
  #reindex and renormalize to deal with omitted species
  Mock <- Mock %>% 
    left_join(sp_list) %>% 
    mutate(level1 = match(level1, unique(level1)),
           speciesname = species,
           species = species_idx) %>% 
    group_by(Sample) %>% 
    mutate(b_proportion = b_proportion/sum(b_proportion)) %>% 
    ungroup()
  
  Observation <- Observation %>% 
    rename(level1 = as.name(Level_1_treatment_envir),
           level2 = as.name(Level_2_treatment_envir))
  
  if (is.na(Level_3_treatment_envir)) {Observation$level3 <- 1} else {Observation <- Observation %>% rename(level3 = as.name(Level_3_treatment_envir))}  
  
  
  Observation <- Observation %>% 
    left_join(sp_list) %>% 
    mutate(speciesname = species,
           species = species_idx) %>% 
    # mutate(station = case_when(stationname == "Dn" ~ 1,
    #                            stationname == "Up" ~ 2)) %>% 
    group_by(level1) %>% 
    mutate(level2 = match(level2, unique(level2))) %>% #reindex, if necess
    mutate(level3 = match(level3, unique(level3))) %>% #reindex, if necess
    unite(c(level1, level2, level3), col = "Sample", sep = ".", remove = F)
  
  station_list <- data.frame(
    station = Observation$level1 %>% unique(),
    station_idx = 1:length(unique(Observation$level1)))
  
  
  #list object containing named elements Observation, Mock, and Npcr
  
  return(
    metabarcoding_data <- list(
      Observation = Observation,
      Mock = Mock,
      N_pcr_mock = N_pcr_mock, 
      NSpecies = nrow(sp_list),
      station_list = station_list,
      sp_list = sp_list
    ))
}

# example
# library(here)
# metabarcoding_data <- format_metabarcoding_data(here("../example_data", "padden_sample.rds"),
#                           here("../example_data", "newmockdata.RDS"))

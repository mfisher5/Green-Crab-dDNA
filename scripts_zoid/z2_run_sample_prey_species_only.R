library(tidyverse)
library(magrittr)
library(here)
library(ggplot2)
library(ggrepel)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

library('devtools')
devtools::install_github("noaa-nwfsc/zoid")

source("calibrate_metabarcoding_RPK.R")


# A note about error messages ---------------------------------------------

# Max Tree Depth is a warning you can ignore.  It just means sampling wasn’t as efficient as it could have been, but the results are still totally valid.
# Divergent transitions are intermediate-serious warnings.  These mean that sampling isn’t working quite right, but if there are small numbers of divergent transitions, that’s not too bad. 
# Rhat warnings are serious and mean that the model isn’t converging; something is broken.
# Low sample-size warnings aren’t serious; you can either ignore them or else run the chains for longer.



# Set up Data -------------------------------------------------------------

EGC_obs <- readRDS(here('../','data','zoid',"zoid_dataLIST_2023-8-29.rds"))



n_pcr_samples <- read_csv(here('qm_model_SAMPLEpcr.csv'))

#reindex so tech reps always start with 1
# EGC_obs[[2]] <- EGC_obs[[2]] %>% 
#   group_by(Sample, species) %>% 
#   mutate(tech = match(tech, unique(tech)))

#check to make sure 10 species are matching (toss out Zosteraceae later)
# unique(intersect(EGC_obs[[1]]$species, EGC_obs[[2]]$species))
# EGC_obs[[1]]$species[which(!(EGC_obs[[1]]$species %in% EGC_obs[[2]]$species))]  # Cottus was not ID'd in diet data


EGC_obs[[2]] %>% 
  filter(station == "WACO21-097" & species == "Alitta")

unique(EGC_obs[[2]]$station)


# Check for Read Depth ----------------------------------------------------
#what are the read depths going into the quant model?
EGC_obs[[2]] %>%
  filter(species %in% EGC_obs[[1]]$species) %>%                # species in quant model
  group_by(station,tech) %>% summarise(total_reads=sum(Nreads)) %>%
  ggplot(aes(x=total_reads)) + geom_histogram() +theme_bw()

png(here('read-depth_qm_model_Zmar_2023-5-27.png'))
EGC_obs[[2]] %>%
  filter(species %in% EGC_obs[[1]]$species) %>%                # species in quant model
  group_by(station,tech) %>% summarise(total_reads=sum(Nreads)) %>%
  ggplot(aes(x=total_reads)) + geom_histogram() + 
  labs(x="Read Depth into Model",y="Number of Crab+Replicates") +theme_bw()
dev.off()

#which crab do we get rid of if we only focus on samples with > 100 reads?
EGC_obs[[2]] %>%
  filter(species %in% EGC_obs[[1]]$species) %>% 
  group_by(station,tech) %>% summarise(total_reads=sum(Nreads)) %>% filter(total_reads < 100)

low_read_crab <- EGC_obs[[2]] %>%
  filter(species %in% EGC_obs[[1]]$species) %>% 
  group_by(station,tech) %>% summarise(total_reads=sum(Nreads)) %>% filter(total_reads < 100)
# 1 WACO21-134 1              56
# 2 WACO21-134 2              36
# 4 WACO21-503 3              33
# 5 WACO21-512 1              69
# 6 WACO21-512 3               2
# 8 WACO21-645 1              57
# 9 WACO21-645 2              38

EGC_obs[[2]] <- EGC_obs[[2]] %>% anti_join(dplyr::select(low_read_crab,station,tech))

# separate sample data frame ----------------------------------------------
n_pcr_samples <- anti_join(n_pcr_samples %>% mutate(tech=as.character(tech)),
                           low_read_crab %>% 
                             dplyr::select(station,tech) %>% mutate(tech=as.character(tech)), by=c("sample"="station","tech"))
#format, sample subset had 45 PCR cycles
sub43.samples <- n_pcr_samples %>% filter(PCRcycles==43) %>% dplyr::select(sample) %>% distinct()
sub40.samples <- n_pcr_samples %>% filter(PCRcycles==40) %>% dplyr::select(sample) %>% distinct()
sub47.samples <- n_pcr_samples %>% filter(PCRcycles==47) %>% dplyr::select(sample) %>% distinct()


# Run for samples with Npcr=40 ---------------------------------------

#format 
my_dDNA_data <- format_metabarcoding_data(EGC_obs[[2]] %>% filter(station %in% sub40.samples$sample), N_pcr_mock=40,
                                          EGC_obs[[1]],
                                          Level_1_treatment_envir <- "station",
                                          Level_2_treatment_envir <- "tech",
                                          Level_3_treatment_envir <- NA,
                                          Level_1_treatment_mock <- "Community",
                                          Level_2_treatment_mock <- "tech",
                                          Level_3_treatment_mock <- NA)
mydata40 <- makeDesign(my_dDNA_data, N_pcr_cycles = 40)


#run models
QM_Bayes_40pcr <- QM_bayes("quant_metabar_multinom.stan", mydata40, NCHAINS = 3, ITER = 5000 )
saveRDS(QM_Bayes_40pcr,file="QM_Bayes_out_40pcr_100rmin_2023-07-21.rds")


# Warning messages:
#   1: There were 13500 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 12. See
# https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 2: Examine the pairs() plot to diagnose sampling problems


# Run for samples with Npcr=43 --------------------------------------------
my_dDNA_data43 <- format_metabarcoding_data(EGC_obs[[2]] %>% filter(station %in% sub43.samples$sample), N_pcr_mock=40,
                                            EGC_obs[[1]],
                                            Level_1_treatment_envir <- "station",
                                            Level_2_treatment_envir <- "tech",
                                            Level_3_treatment_envir <- NA,
                                            Level_1_treatment_mock <- "Community",
                                            Level_2_treatment_mock <- "tech",
                                            Level_3_treatment_mock <- NA)
mydata43 <- makeDesign(my_dDNA_data43, N_pcr_cycles = 43)



#run models
QM_Bayes_43pcr <- QM_bayes("quant_metabar_multinom.stan", mydata43, NCHAINS = 3, ITER = 10000)
saveRDS(QM_Bayes_43pcr,file="QM_Bayes_out_43pcr_100rmin_2023-07-21.rds")
# Warning messages:
# 1: There were 28500 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 12. See
# https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 2: Examine the pairs() plot to diagnose sampling problems


compare_Bayes_Reads(QM_Bayes_43pcr, EGC_obs[[2]])


# Run for samples with Npcr=47 --------------------------------------------
my_dDNA_data47 <- format_metabarcoding_data(EGC_obs[[2]] %>% filter(station %in% sub47.samples$sample), N_pcr_mock=40,
                                            EGC_obs[[1]],
                                            Level_1_treatment_envir <- "station",
                                            Level_2_treatment_envir <- "tech",
                                            Level_3_treatment_envir <- NA,
                                            Level_1_treatment_mock <- "Community",
                                            Level_2_treatment_mock <- "tech",
                                            Level_3_treatment_mock <- NA)
mydata47 <- makeDesign(my_dDNA_data47, N_pcr_cycles = 47)


#run models
QM_Bayes_47pcr <- QM_bayes("quant_metabar_multinom.stan", mydata47, NCHAINS = 3, ITER = 5000)
saveRDS(QM_Bayes_47pcr,file="QM_Bayes_out_47pcr_100rmin_2023-07-21.rds")


# Combine QM output -------------------------------------------------------

QM_Bayes_all <- bind_rows((QM_Bayes_40pcr$Bayes_estimates %>%
                             rownames_to_column("crab") %>%
                             pivot_longer(cols=2:(dim(QM_Bayes_40pcr$Bayes_estimates)[2]+1)) %>%
                             mutate(npcr=40)),
                          (QM_Bayes_43pcr$Bayes_estimates %>%
                             rownames_to_column("crab") %>%
                             pivot_longer(cols=2:(dim(QM_Bayes_43pcr$Bayes_estimates)[2]+1)) %>%
                             mutate(npcr=43)),
                          (QM_Bayes_47pcr$Bayes_estimates %>%
                             rownames_to_column("crab") %>%
                             pivot_longer(cols=2:(dim(QM_Bayes_47pcr$Bayes_estimates)[2]+1))) %>%
                            mutate(npcr=47))
bayes_clean <- filter(QM_Bayes_all, value > 0.001) %>%  mutate(name=ifelse(name=="Mya","Mya arenaria",name))

bayes_clean %>% ggplot(aes(x=crab,y=value,fill=name)) +
  geom_col() + 
  theme(axis.text.x=element_blank())


# Save Input Dataframes ---------------------------------------------------

#for fitting zoid!
all_dDNA_data_Observation <- my_dDNA_data$Observation %>% mutate(npcr=40) %>%
  bind_rows(my_dDNA_data43$Observation %>% mutate(npcr=43)) %>%
  bind_rows(my_dDNA_data47$Observation %>% mutate(npcr=47))

saveRDS(all_dDNA_data_Observation,here('../zoid_avg_crab/QM_Bayes_my_dDNA_data_observation.rds'))






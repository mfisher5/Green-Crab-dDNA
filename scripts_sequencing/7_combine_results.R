############ Data Wrangling: Results files ############
# 
# get blast outputs from different sequencing runs into one data frame.
# with metadata.
#
#######################################################


# Set up ------------------------------------------------------------------
library(tidyverse)
library(here)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(janitor)
library(magrittr)

# User inputs
indir   <- 'data/blast'
run.nums <- c(2,3,4,5)
marker  <- 'BF3'
metadata_dir <- 'data/metadata'



# Metadata ----------------------------------------------------------------

# metadata for run
for(r in run.nums){
  if(r==run.nums[1]){
    meta.run <- read.csv(here(metadata_dir,paste0("run",r,"_",marker,"_run_metadata.csv"))) %>% clean_names()}
  else{
    meta.run %<>% bind_rows(read.csv(here(metadata_dir,paste0("run",r,"_",marker,"_run_metadata.csv"))) %>%
                              clean_names())
    
  }
}

# metadata filtered (crab with prey info)
for(r in run.nums){
  if(r==run.nums[1]){
    meta <- read.csv(here(metadata_dir,paste0("run",r,"_",marker,"_filtered_metadata.csv"))) %>% clean_names()}
  else{
    meta %<>% bind_rows(read.csv(here(metadata_dir,paste0("run",r,"_",marker,"_filtered_metadata.csv")))
                        %>% clean_names())
    
  }
}

## not all filtered metadata has "source" column
meta %<>% mutate(source=ifelse(is.na(source) & grepl("WACO",sample), "Sample",
                               ifelse(is.na(source) & grepl("MC",sample), "Mock Community",source)))




# Blast data --------------------------------------------------------------
# read blast data for each run into a list
dat.list <- list()
for(i in seq(1,length(run.nums))){
  r <- run.nums[i]
  dat <- read_csv(here(indir, paste0("run",r,"_",marker,"_taxonomy_filtered.csv"))) %>%
    filter(!(taxon %in% c("Homo sapiens","Scomber japonicus","Sus scrofa"))) %>% mutate(MiSeqRun=r)
    dat.list[[i]] <- dat
}


lapply(dat.list,colnames)


## in common: [1] "Sample_name" "Locus"       "sample"      "tech"        "Hash"        "nReads"      "rank"       
# [9] "kingdom"     "phylum"      "class"       "order"       "family"      "genus"       "species"     "taxon" 

## runs 2-4 use "type", and run 5 uses "source"
dat.list[[4]] <- dat.list[[4]] %>% rename(type=source)

## runs 2 and 4 don't have mc_rm, because they didn't have mock communities
dat.list[[1]] <- dat.list[[1]] %>% mutate(mc_rm=0) %>% dplyr::select(-Locus.y) %>% rename("Locus"=Locus.x)
dat.list[[2]] <- dat.list[[2]] %>% dplyr::select(-Locus.y) %>% rename("Locus"=Locus.x)
dat.list[[3]] <- dat.list[[3]] %>% mutate(mc_rm=0)

dat.list.thin <- purrr::map(dat.list,dplyr::select,c("Sample_name","Locus","sample","tech", 
                                                     "Hash","nReads","rank","kingdom","phylum","class","order","family","genus","species","taxon","mc_rm","type","MiSeqRun"))


## ok, now make a data frame
dat <- bind_rows(dat.list.thin)


# Sample reruns -----------------------------------------------------------
## now deal with duplicated samples between runs 3 and 4, by re-naming the technical replicates in run 4 (to be 4-6)
dat.out <- dat %>% filter(MiSeqRun %in% c(2,3,5))

dat3 <- filter(dat,MiSeqRun==3)
dat4 <- filter(dat,MiSeqRun==4) %>%
  mutate(rerun=ifelse(sample %in% dat3$sample, 1,0)) %>%
  mutate(tech=ifelse(rerun==1,tech+3,tech))

View(filter(dat4,rerun==1)) # just checking
View(filter(dat4,rerun==0 & tech %in% c(4,5,6))) # just checking

dat4 %<>% dplyr::select(-rerun)

dat.out %<>% bind_rows(dat4)


## deal with the sample re-run between 2 and 3, 106....
tmp106 <- dat.out %>% filter(sample=="WACO21-106")
tmp106 %<>% mutate(tech=ifelse(MiSeqRun==3,match('g',letters),tech))

dat.out %<>% filter(sample!="WACO21-106")

dat.out %<>% bind_rows(tmp106)


# Add Metadata to Blast ---------------------------------------------------

dat.out %<>% left_join(meta %>% dplyr::select(sample,egc_id,site_month),by=c("sample"))
dim(dat.out %>% filter(type=="Sample" & is.na(egc_id)))
dim(dat.out %>% filter(type=="Sample" & is.na(site_month)))

dat.out %<>% distinct()

# write_csv(dat.out, here('results',paste0('runs',first(run.nums),'-',last(run.nums),'_',marker,'_CEG05blast_taxonomy_filtered.csv')))
write_csv(dat.out, here(indir,paste0('runs',first(run.nums),'-',last(run.nums),'_',marker,'_taxonomy_filtered.csv')))


############ Data Wrangling: Unique Taxa / Species IDs ############
# 
# split the final sequence data into species-level IDs, and unique
# taxa identifications inclusive of less specific results.
#
###################################################################


# get species-level prey for each crab
spdat <- dat.out %>% filter(rank=="species")
# write_csv(spdat, here('results','allRuns_BF3_CEG05Blast_filtered_species.csv'))
write_csv(spdat, here('data','results','allRuns_BF3_TESTING_filtered_species.csv'))


#  keep non-specific ids in data set only if that taxon is not represented for the given crab -- for all plant / algae IDs that only go to order thru genus, otherwise just through genus
nonspec_dat <- dat.out %>% 
  filter((phylum %in% c("Chlorophyta","Rhodophyta","Streptophyta") & rank=="order") | (class %in% c("Phaeophyceae") & rank=="order") | rank %in% c("genus","family")) %>%
  mutate(dup_sample=sample)

nonspec_dat %<>%
  group_by(sample) %>%
  nest()

remove_duplicate_taxonomy <- function(x){
  tmp_sample <- as.character(unique(x[,"dup_sample"]))
  
  if(tmp_sample %in% spdat$sample){
    green.filter <- x %>% filter(phylum %in% c("Chlorophyta","Rhodophyta","Streptophyta") | (class %in% c("Phaeophyceae","Chrysophyceae"))) %>%
      filter(!(taxon %in% (filter(spdat,sample==tmp_sample))$genus) & !(taxon %in% (filter(spdat,sample==tmp_sample))$family) & !(taxon %in% (filter(spdat,sample==tmp_sample))$order))
    x.filter <- x %>% filter(!(phylum %in% c("Chlorophyta","Rhodophyta","Streptophyta")) & (!(class %in% c("Phaeophyceae","Chrysophyceae")))) %>%
      filter(!(taxon %in% (filter(spdat,sample==tmp_sample))$genus) & !(taxon %in% (filter(spdat,sample==tmp_sample))$family)) %>% 
    bind_rows(green.filter)
    return(x.filter)
    
  } else{return(x)}
}
# check fx
# test <- (nonspec_dat$data[which(nonspec_dat$sample=="MC-B")])[[1]]
# tmp_sample <- as.character(unique(test[,"dup_sample"]))
# tmp_sample %in% spdat$sample
# test.filter <- test %>% mutate(non_green=!(phylum %in% c("Chlorophyta","Rhodophyta","Streptophyta")) & (!(class %in% c("Phaeophyceae","Chrysophyceae")))) %>%
#   mutate(keep_genus=ifelse(!(taxon %in% (filter(spdat,sample==tmp_sample))$genus),1,0)) %>%
#   mutate(keep_family=ifelse(!(taxon %in% (filter(spdat,sample==tmp_sample))$family),1,0))



nonspec_dat %<>% mutate(data.filter=purrr::map(data,remove_duplicate_taxonomy))

unqdat <- nonspec_dat %>% dplyr::select(-data) %>% unnest(c(sample,data.filter)) %>% dplyr::select(-dup_sample) %>%
  bind_rows(spdat)


# remove duplicates that have two + entries of higher level taxonomy!
unqdat %>% 
  filter(rank!="species") %>% 
  dplyr::select(sample,Hash,taxon,rank,species,genus,family,order,class) %>% distinct() %>%
  group_by(Hash, sample) %>% summarise(duplicated.ranks=length(unique(rank))) %>% filter(duplicated.ranks > 1)

# fix hash 6ab46a5ec84ed667fae6c421dc0b59c9cc86727f
tmpdat <- filter(unqdat, Hash=="6ab46a5ec84ed667fae6c421dc0b59c9cc86727f" & sample %in% c("WACO21-288","WACO21-607")) %>%
  filter(rank=="genus")

unqdat %<>% filter(  !( (Hash=="6ab46a5ec84ed667fae6c421dc0b59c9cc86727f") & (sample %in% c("WACO21-288","WACO21-607")) )  ) %>%
  bind_rows(tmpdat)


write_csv(unqdat, here('data','results','allRuns_BF3_TESTING_filtered_unique_taxa.csv'))

# summarize taxonomy ------------------------------------------------------


unq.summary <- unqdat %>%
  group_by(taxon,rank,phylum,class,family,genus,species,type) %>%
  mutate(type=ifelse(type=="Samples","Sample",type)) %>%    # fix inconsistent naming
  summarise(n_crab=length(unique(sample)),
            miseqruns=paste0(unique(MiSeqRun),collapse=","),
            site_months=paste0(unique(site_month),collapse=","))
write_csv(unq.summary, here('data','results','allRuns_BF3_TESTING_filtered_unique_taxa_summary_SUBMIT.csv'))











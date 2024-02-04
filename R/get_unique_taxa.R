get_unique_taxa <- function(taxa.df, level="all"){
  
  if(level=="all"){
  # get species-level prey for each crab
  spdat <- taxa.df %>% filter(rank=="species") %>% dplyr::select(taxon,genus,family,order,class) %>% distinct() %>%
    pivot_longer(cols=c(genus,family,order,class),names_to="tax.level",values_to="sp.taxonomy") %>%
    rename("sp.taxon"=taxon)
  
  # keep non-specific ids in data set only if that taxon is not represented across the entire dataset
  nonspec_dat <- taxa.df %>% filter(rank!="species") %>% dplyr::select(taxon,rank,species,genus,family,order,class) %>% distinct()
  
  nonspec_dat %<>% pivot_longer(cols=c(genus,family,order,class),names_to="tax.level",values_to="taxonomy") %>%
    filter(tax.level==rank) %>%
    left_join(spdat,by=c("tax.level","taxonomy"="sp.taxonomy")) %>% filter(!is.na(sp.taxon))
  
  dat <- taxa.df %>% anti_join(nonspec_dat,by=c("rank","taxon"="taxonomy"))
  
  # remove higher level taxonomy represented by genus
  gdat <- dat %>% filter(rank=="genus")
  dat <- dat %>%
    filter(!(taxon %in% gdat$family) & !(taxon %in% gdat$order))
  
  # remove higher level taxonomy represented by family
  fdat <- dat %>% filter(rank=="family")
  dat <- dat %>%
    filter(!(taxon %in% fdat$order))
  
  
  
  } else if(level=="site"){
    spdat <-taxa.df %>% filter(rank=="species") %>% 
      dplyr::select(site,taxon,genus,family,order,class) %>% distinct() %>%
      pivot_longer(cols=c(genus,family,order,class),names_to="tax.level",values_to="sp.taxonomy") %>%
      rename("sp.taxon"=taxon)
    
    nonspec_dat <-taxa.df %>% 
      filter(rank!="species") %>% 
      dplyr::select(site,taxon,rank,species,genus,family,order,class) %>% distinct() %>%
      pivot_longer(cols=c(genus,family,order,class),names_to="tax.level",values_to="taxonomy") %>%
      filter(tax.level==rank) %>%
      left_join(spdat,by=c("site","tax.level","taxonomy"="sp.taxonomy")) %>% filter(!is.na(sp.taxon))
    
   dat <-taxa.df %>%
      anti_join(nonspec_dat,by=c("site","rank","taxon"="taxonomy"))
   
   # remove higher level taxonomy represented by genus
   gdat <- dat %>% filter(rank=="genus") %>%
     dplyr::select(site,family,order) %>% distinct()
   keep_fodat <- dat %>% filter(rank %in% c("family","order")) %>%
     anti_join(gdat, by=c("site","family","order"))
   
   dat <- dat %>% filter(!(rank %in% c("family","order"))) %>%
     bind_rows(keep_fodat)
     
   # remove higher level taxonomy represented by family
   fdat <- dat %>% filter(rank=="family") %>%
     dplyr::select(site,order) %>% distinct()
   
   keep_odat <- dat %>% filter(rank %in% c("order")) %>%
     anti_join(fdat, by=c("site","order"))
   
   dat <- dat %>% filter(!(rank %in% c("order"))) %>%
     bind_rows(keep_odat)
   
   
   

   
  } else{stop("Please choose level = site or level = all")}
  
  return(dat)
  
}

get_unique_taxa <- function(taxa.df, level="all", return.removed=FALSE){
  
  if(level=="all"){
    
  # get species-level prey for all crab
  spdat <- taxa.df %>% filter(rank=="species") %>% dplyr::select(taxon,genus,family,order,class) %>% distinct() %>%
    pivot_longer(cols=c(genus,family,order,class),names_to="tax.level",values_to="sp.taxonomy") %>%
    rename("sp.taxon"=taxon)
  
  # keep non-specific ids in data set only if that taxon is not represented across the entire dataset
  nonspec_dat <- taxa.df %>% filter(rank!="species") %>% dplyr::select(taxon,rank,species,genus,family,order,class) %>% distinct()
  
  nonspec_dat %<>% pivot_longer(cols=c(genus,family,order,class),names_to="tax.level",values_to="taxonomy") %>%
    filter(tax.level==rank) %>%
    left_join(spdat,by=c("tax.level","taxonomy"="sp.taxonomy")) %>% filter(!is.na(sp.taxon))
  
  # remove the nonspecific IDs that matched to species IDs. The updated data set is 'dat'
  dat <- taxa.df %>% anti_join(nonspec_dat,by=c("rank","taxon"="taxonomy"))
  out <- taxa.df %>% anti_join(nonspec_dat,by=c("rank","taxon"="taxonomy"))
  if(return.removed==TRUE){
    outR <- taxa.df %>% left_join(nonspec_dat %>% 
                                    dplyr::select(taxonomy, rank) %>% mutate(matched="species"),by=c("rank","taxon"="taxonomy")) %>%
      filter(!is.na(matched))
    }
  
  # remove higher level taxonomy that is already represented by genus / family / order / class
  levels <- c("genus","family","order","class","phylum")
  for(t in seq(1:(length(levels)-1))){
    tdat <- dat %>% filter(rank==levels[t]) %>%  # get genus-level identifications
      dplyr::select(all_of(levels)) %>% distinct() %>%
      pivot_longer(cols=levels[(t+1):length(levels)], names_to="rank",values_to="taxon")
    
    out <- out %>% anti_join(tdat, by=c("rank","taxon"))
    
    if(return.removed==TRUE){outR <- outR %>% 
      bind_rows(left_join(taxa.df,tdat %>% dplyr::select(taxon,rank) %>% mutate(matched=levels[t]),
                          keep=FALSE,by=c("taxon","rank"))) %>%
        filter(!is.na(matched))}
    }
  
  } else if(level=="site"){
    # get species-level prey for all crab
    spdat <-taxa.df %>% filter(rank=="species") %>% 
      dplyr::select(site,taxon,genus,family,order,class) %>% distinct() %>%
      pivot_longer(cols=c(genus,family,order,class),names_to="tax.level",values_to="sp.taxonomy") %>%
      rename("sp.taxon"=taxon)
    
    # keep non-specific ids in data set only if that taxon is not represented by species IDs
    nonspec_dat <-taxa.df %>% 
      filter(rank!="species") %>% 
      dplyr::select(site,taxon,rank,species,genus,family,order,class) %>% distinct() %>%
      pivot_longer(cols=c(genus,family,order,class),names_to="tax.level",values_to="taxonomy") %>%
      filter(tax.level==rank) %>%
      left_join(spdat,by=c("site","tax.level","taxonomy"="sp.taxonomy")) %>% filter(!is.na(sp.taxon))
    
    
    # remove the nonspecific IDs that matched to species IDs. The updated data set is 'dat'
   dat <-taxa.df %>%
      anti_join(nonspec_dat,by=c("site","rank","taxon"="taxonomy"))
   if(return.removed==TRUE){
     outR <- taxa.df %>% left_join(nonspec_dat %>% 
                                     dplyr::select(taxonomy, rank) %>% mutate(matched="species"),
                                   by=c("rank","taxon"="taxonomy")) %>%
       filter(!is.na(matched))
   }
   
   # isolate non-specific IDs to compare to each other
   keep_higherdat <- dat %>% filter(rank %in% c("family","order","class","phylum"))
   
   # remove higher level taxonomy that is already represented by genus / family / order / class
   levels <- c("genus","family","order","class","phylum")
   for(t in seq(1:(length(levels)-1))){
     tdat <-  dat %>% filter(rank==levels[1]) %>%
       dplyr::select('site',levels[t:length(levels)]) %>% distinct() %>% # genus-level identifications by site
       pivot_longer(cols=levels[(t+1):length(levels)], names_to="rank",values_to="taxon")
     
     colnames(tdat)[2] <- paste0(levels[t],".matched")
     keep_higherdat <- keep_higherdat %>%  # get higher level IDs that are not represented by genus, by site
       left_join(tdat, by=c("site","rank","taxon"))
     
     if(t==1){cols_to_rmv <- c(paste0(levels[t],".matched"))
     } else{cols_to_rmv <- c(cols_to_rmv,paste0(levels[t],".matched"))}
     
   }
   
   ## what higher-level IDs will be removed? ##
   if(return.removed==TRUE){
     outR <- outR %>% bind_rows(
       keep_higherdat %>% 
         pivot_longer(ends_with("matched"),names_to="matched",values_to="tmp") %>%
         filter(!is.na(tmp)) %>% dplyr::select(-tmp)
       )
   }

   ## what higher-level IDs can be kept? ##
   keep_higherdat <- keep_higherdat %>% 
     filter_at(vars(ends_with("matched")), all_vars(is.na(.))) %>%
     dplyr::select(-all_of(cols_to_rmv))
   
   # remove all higher 
   out <- dat %>% filter(rank %in% c("genus","species")) %>%
     bind_rows(keep_higherdat)

   
  } else{stop("Please choose level: site, all")}
  
  
  if(return.removed==FALSE){
    return(out)
  }else{
    return(list(out,outR))
  }
  
}









# gdat <- dat %>% filter(rank=="genus")  # get genus-level identifications
# dat <- dat %>%
#   filter(!(taxon %in% gdat$family) & !(taxon %in% gdat$order) & 
#            !(taxon %in% gdat$class) & !(taxon %in% gdat$phylum))
# 
# # remove higher level taxonomy that is already represented by family IDs
# fdat <- dat %>% filter(rank=="family")
# dat <- dat %>%
#   filter(!(taxon %in% fdat$order) & 
#            !(taxon %in% fdat$class) & !(taxon %in% fdat$phylum))
# 
# # remove higher level taxonomy already represented by order IDs
# odat <- dat %>% filter(rank=="order")
# if(dim(odat)[1] > 0){
#   dat <- dat %>%
#     filter(!(taxon %in% odat$class) & !(taxon %in% odat$phylum))
# }
# 
# # remove phyla already represented by class
# cdat <- dat %>% filter(rank=="class")
# if(dim(cdat)[1] > 0){
#   dat <- dat %>%
#     filter(!(taxon %in% cdat$phylum))}




# gdat <- dat %>% filter(rank=="genus") %>%
#   dplyr::select(site,genus, family,order,class, phylum) %>% distinct() %>% # genus-level identifications by site
#   pivot_longer(cols=c(family,order,class,phylum), names_to="rank",values_to="taxon")
# 
# keep_higherdat <- keep_higherdat %>%  # get higher level IDs that are not represented by genus, by site
#   left_join(gdat %>% rename(genus.matched=genus), by=c("site","rank","taxon"))
# 
# # remove higher level taxonomy represented by family
# fdat <- dat %>% filter(rank=="family") %>%
#   dplyr::select(site,family,order,class, phylum) %>% distinct() %>%
#   pivot_longer(cols=c(order,class,phylum), names_to="rank",values_to="taxon")
# 
# keep_higherdat <- keep_higherdat %>%  # get higher level IDs that are not represented by genus, by site
#   left_join(fdat %>% rename(fam.matched=family), by=c("site","rank","taxon"))
# 
# # remove higher level taxonomy represented by order
# odat <- dat %>% filter(rank == "order") %>%
#   dplyr::select(site,order,class, phylum) %>% distinct() %>%
#   pivot_longer(cols=c(class,phylum), names_to="rank",values_to="taxon")
# 
# keep_higherdat <- keep_higherdat %>%  # get higher level IDs that are not represented by genus, by site
#   left_join(odat %>% rename(ord.matched=order), by=c("site","rank","taxon"))
# 
# # remove higher level taxonomy represented by class
# cdat <- dat %>% filter(rank == "class") %>%
#   dplyr::select(site,class, phylum) %>% distinct() %>%
#   pivot_longer(cols=c(phylum), names_to="rank",values_to="taxon")
# 
# keep_higherdat <- keep_higherdat %>%  # get higher level IDs that are not represented by genus, by site
#   left_join(cdat %>% rename(class.matched=class), by=c("site","rank","taxon"))

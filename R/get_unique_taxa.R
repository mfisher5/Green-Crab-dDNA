#' Get unique taxa
#'
#' `get_unique_taxa` removes higher-level taxa (i.e., non-specific identifications)
#' that may represent duplicates. For example, *Decapoda* would be removed if there is
#' an identification for *Cancer magister*, and *Crangon* would be removed if there is 
#' an identification for *Crangon franciscorum*.  
#'
#' 
#' @param taxa.df  a dataframe of of the taxa (rows) found each in sample (rows)
#' @param level    remove potentially duplicated non-specific IDs across the full data set, or within sites? 
#' @param return.removed return only a data frame without potentially duplicated non-specific IDs (FALSE), or return a list with that data frame plus a dataframe of what was removed (TRUE)

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



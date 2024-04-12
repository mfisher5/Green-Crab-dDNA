## Results

Final data files from DNA metabarcoding and mixture models, and files with summary information for main text / supplemental tables. 

### DNA metabarcoding 


**allRuns_BF3_filtered_FINAL_unique_taxa** : The Final Dataset. A file with all unique+ ASVs that were matched to taxonomy in BLAST, were filtered for predator / parasite / non-target DNA, were manually filtered, and were checked again for potentially duplicated higher-level taxon identifications. 

*allRuns_BF3_filtered_FINAL_unique_taxa_summary* : Same as above, grouped by taxon and summarized to show crabs per taxon detected plus associated MiSeq runs, site-months


*MockCommunities_BF_filtered_FINAL_input_taxa.csv* : The Final Dataset for the mock communities (used to calibrate metabarcoding data)



*+ We conservatively removed higher-level taxon IDs that potentially duplicated species-level IDs within each crab. For example, if a crab contained an identification of* Crangon *and* Crangon franciscorum, *the Crangon (genus) identification was removed. This was completed after each filtering step.*

<br>
<br>

### Mixture model for average diet

**zoid_fitted_vals_preyD100** : zoid output (proportion of calibrated species' DNA in an "average" green crab diet) when running fit_zoid with all crabs from all sites

**zoid_fitted_vals_CLAMBED_preyD100** : same as above, but only for green crabs trapped on clam beds

**zoid_fitted_vals_SLOUGH_preyD100** : same as above, but only for green crabs trapped at intertidal sloughs


<br>
<br>


#### Tables:

Table S4. Number of crabs per taxa, plus summary counts for higher taxonomic levels. Extension of Table 2 in main text.

Table S1. Summary of ASVs


<br>
<br>
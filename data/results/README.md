## Results

Final data files from DNA metabarcoding and mixture models, and files with summary information for main text / supplemental tables. 

### Sequencing


**allRuns_BF3_filtered_FINAL_unique_taxa** : The Final Dataset. A file with all unique+ ASVs that were matched to taxonomy in BLAST, were filtered for predator / parasite / non-target DNA, were manually filtered, and were checked again for potentially duplicated higher-level taxon identifications. 

*allRuns_BF3_filtered_FINAL_unique_taxa_summary* : Same as above, grouped by taxon and summarized to show crabs per taxon detected plus associated MiSeq runs, site-months

*+ We conservatively removed higher-level taxon IDs that potentially duplicated species-level IDs within each crab. For example, if a crab contained an identification of* Crangon *and* Crangon franciscorum, *the Crangon (genus) identification was removed. This was completed after each filtering step.*



### Mixture Model 

**zoid_fitted_vals_preyD100** : zoid output (proportion of calibrated species' DNA in an "average" green crab diet) when running fit_zoid with all crabs from all sites

**zoid_fitted_vals_CLAMBED_preyD100** : same as above, but only for green crabs trapped on clam beds

**zoid_fitted_vals_SLOUGH_preyD100** : same as above, but only for green crabs trapped at intertidal sloughs



#### Tables:

Table S4. Number of crabs per taxa, plus summary counts for higher taxonomic levels. Extension of Table 2 in main text.

Table S1. Summary of ASVs


<br>
<br>


*Versioning:* 

zoid output from 2023-08-29.

allRuns_BF3_filtered_unique_taxa: A file with all unique+ ASVs that were matched to taxonomy in BLAST and were filtered for predator / parasite / non-target DNA.

allRuns_BF3_filtered_unique_taxa_summary: Same as above, grouped by taxon and summarized to show crabs per taxon detected plus associated MiSeq runs, site-months

allRuns_BF3_filtered_unique_taxa_summary_blastPIdent: The mean, min, and max percent identity scores for each ASV - crab sample combination. Useful for applying the manual filter.
<br>

allRuns_BF3_filtered_unique_taxa_summary_MANUALfilter: A file with manual filtering information (including final taxon / rank, reason for change or removal) for each unique+ ASV from *allRuns_BF3_filtered_unique_taxa*.


allRuns_BF3_filtered_unique_taxa_manualFilter*: A file with all ASVs that were matched to taxonomy in BLAST, were filtered for predator / parasite / non-target DNA, and were manually filtered. Essentially *allRuns_BF3_filtered_unique_taxa* with the manual filter applied.



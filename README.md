# Green-crab-dDNA

### Revealing European green crab predation in a Washington State estuary with DNA metabarcoding (dDNA)
<br>
For [Fisher et al. in review](www.doi.org/10.1371/journal.pone.0302518

Authors: Mary Fisher, Emily Grason, Alex Stote, Ryan Kelly, Kate Litle, Sean McDonald

Abstract: Predation by invasive species can threaten local ecosystems and economies. The European green crab (Carcinus maenas), one of the most widespread marine invasive species, is an effective predator associated with clam and crab population declines outside of its native range. In the U.S. Pacific Northwest, green crab has recently increased in abundance and expanded its distribution, generating concern for estuarine ecosystems and associated aquaculture production. However, regionally-specific information on the trophic impacts of invasive green crab is very limited. We compared the stomach contents of green crabs collected on clam aquaculture beds versus intertidal sloughs in Willapa Bay, Washington, to provide the first in-depth description of European green crab diet at a particularly crucial time for regional management. We first identified putative prey items using DNA metabarcoding of stomach content samples. We compared diet composition across sites using prey presence/absence and an index of species-specific relative abundance. For eight prey species, we also calibrated metabarcoding data to quantitatively compare DNA abundance between prey taxa, and to describe an ‘average’ green crab diet at an intertidal slough versus a clam aquaculture bed. From the stomach contents of 61 green crabs, we identified 54 unique taxa belonging to nine phyla. The stomach contents of crabs collected from clam aquaculture beds were significantly different from the stomach contents of crabs collected at intertidal sloughs. Across all sites, arthropods were the most frequently detected prey, with the native hairy shore crab (Hemigrapsus oregonensis) the single most common prey item. Of the eight species calibrated with a quantitative model, two ecologically-important native species – the sand shrimp (Crangon franciscorum) and the Pacific staghorn sculpin (Leptocottus armatus) – had the highest average DNA abundance when detected in a stomach content sample. In addition to providing timely information on green crab diet, our research demonstrates the novel application of a recently developed model for more quantitative DNA metabarcoding. This represents another step in the ongoing evolution of DNA-based diet analysis towards producing the quantitative data necessary for modeling invasive species impacts.


*Bettina Thalinger, Georgina Cordone, Eily Allen, Erin D’Agnese, Megan Shaffer, and Maya Garber-Yonts shared invaluable lab and bioinformatic knowledge and assistance, and Eric Ward provided help troubleshooting the R package zoid.*  
___________________


### Workflow

The final, filtered file of ASVs used for diet composition analysis is: [data/results/allRuns_BF3_filtered_FINAL_unique_taxa](https://github.com/mfisher5/Green-crab-dDNA/blob/main/data/results/allRuns_BF3_filtered_FINAL_unique_taxa.csv)

![doc-worksflow-img](https://github.com/mfisher5/Green-crab-dDNA/blob/main/doc/analysis_workflow.png?raw=true)

The raw sequencing data, organized by MiSeq run number and demultiplexed, are available on 
[Zenodo](https://doi.org/10.5281/zenodo.10850508)



### Requirements

Scripts written in : R version 4.1.3 (2022-03-10), GNU bash version 4.4.23(2)-release (x86_64-pc-msys)

NCBI BLAST database matching with blastn was run on a high-performance compute (HPC) cluster (University of Washington's HYAK)


### Folders
- `R` custom R functions
- `bash` custom bash scripts 
- `data` all input, intermediate, and output data files
- `doc` lab protocols with template spreadsheets; analysis workflow diagram
- `figs` png or tif files with figures for main text and supplement
- `scripts_sequencing` scripts 1-10 for processing and analysis of DNA metabarcoding data
- `scripts_quant` scripts qm1-4 for running the Shelton et al. calibration model and the zoid R package mixture model




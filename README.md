# Green-crab-dDNA

### Revealing European green crab predation in a Washington State estuary with DNA metabarcoding (dDNA)
<br>
For [Fisher et al. (in review)](doi.org)

Authors: Mary Fisher, Emily Grason, Alex Stote, Ryan Kelly, Kate Litle, Sean McDonald

Abstract: Predation by invasive species can threaten local ecosystems and economies. The European green crab (Carcinus maenas), one of the most widespread marine invasive species, is an effective predator associated with clam and crab population declines outside of its native range. In the U.S. Pacific Northwest, green crab have recently experienced increases in abundance and expanding distributions, generating concern for estuarine ecosystems and associated aquaculture production. However, regionally-specific information on the trophic impacts of invasive green crab is highly limited. We compared the stomach contents of green crabs collected on aquaculture beds versus natural intertidal sloughs in Willapa Bay, Washington, to provide the first in-depth description of European green crab diet at a particularly crucial time for regional management. We first identified putative prey items using DNA metabarcoding of stomach content samples. DNA metabarcoding provided presence/absence data and was used to compare a given prey species’ relative abundance across site types. For eight key prey species, we also applied a quantitative model to compare DNA abundance between prey species, and to describe an ‘average’ green crab diet comprised of those species at an aquaculture bed and an intertidal slough. From the stomach contents of 62 green crabs, we identified 56 unique taxa belonging to nine phyla. The stomach contents of crabs collected from aquaculture beds were significantly different from the stomach contents of crabs collected at natural intertidal sloughs. Across all sites, arthropods were the most frequently detected prey, with the native hairy shore crab (Hemigrapsus oregonensis) the single most common prey item. Of the eight species included in the quantitative model, two ecologically-important native species – the sand shrimp (Crangon franciscorum) and the Pacific staghorn sculpin (Leptocottus armatus) – were the most abundant in crab stomach contents, when present. In addition to providing timely information on green crab diet, our research demonstrates the novel application of a recently developed model for quantitative DNA metabarcoding. This represents another step in the ongoing evolution of DNA-based diet analysis towards producing the quantitative data necessary for modeling invasive species impacts.


*Bettina Thallinger, Georgina Cordone, Eily Allen, Erin D’Agnese, Megan Schaffer, and Maya Garber-Yonts shared invaluable lab and bioinformatic knowledge and assistance, and Eric Ward provided help troubleshooting the R package zoid.*  
___________________


### Workflow

![doc-worksflow-img](https://github.com/mfisher5/Green-crab-dDNA/blob/main/doc/analysis_workflow.png?raw=true)

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




########## Pull data for non-native species ##########
# 
# Apr 27 2023 - M Fisher
#
#####################################################


# Set up  -----------------------------------------------------------------

library(tidyverse)
library(here)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(janitor)
library(magrittr)


# User inputs
dada2dir <- 'data/dada2'
blastdir <- 'data/blast'
indir   <- 'data/results'
run.nums <- c(2,3,4,5)
marker  <- 'BF3'
metadata_dir <- 'data/metadata'
blast_filename <- 'hash_key_clean_blast_2023-04-26.fasta'



# Read in data ------------------------------------------------------------

# raw blasted taxonomy data (for p ident and e values)
for(i in seq(1,length(run.nums))){
  r=run.nums[[i]]
  if(r==run.nums[1]){
    blastdat <- read_delim(here(blastdir, paste0("run",r,"_",blast_filename)), col_names=FALSE, col_types=rep("c",times=18)) %>% mutate(MiSeqRun=r)
  } else{
    blastdat %<>% bind_rows(read_delim(here(blastdir, paste0("run",r,"_",blast_filename)), col_names=FALSE, col_types=rep("c",times=18))
                            %>% mutate(MiSeqRun=r,X15=as.character(X15)))
  }
}
colnames(blastdat) <- c("qseqid", "sseqid", "sacc", "pident", "length",
                        "mismatch", "gapopen", "qcovus", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxids", "qlen", "sscinames", "sseq","MiSeqRun")

# filtered and processed blasted taxonomy data (for hashes)
dat <- read_csv(here(blastdir,paste0('runs',min(run.nums),'-',max(run.nums),'_',marker,'_taxonomy_filtered.csv')))


# Get sequence data for species of interest -------------------------------

s <- "Botrylloides violaceus"


asvdat <- dat %>% filter(taxon==s) %>% dplyr::select(Hash,MiSeqRun) %>% distinct()
dat.out <- asvdat %>% left_join(blastdat, by=c('Hash'='qseqid','MiSeqRun'))


write.csv(dat.out,file=here(paste0('tmpfile_',s,'.csv')))




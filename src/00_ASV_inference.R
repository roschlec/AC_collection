# ASV inference of 16S amplicon sequencing of epiphytic and endophytic bacterial 
# communities from Arabidopsis thaliana
# Author: Rudolf Schlechter

# Load libraries ----------------------------------------------------------
library(here)
library(tidyverse)
library(dada2)

# Files -------------------------------------------------------------------
path <- # use here() to add path of 16S read files
list.files(path)

# Forward and reverse fastq filenames have format: 
# SAMPLENAME_1.fastq.gz and SAMPLENAME_2.fastq.gz
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Quality Control ---------------------------------------------------------
# Exploratory
#plotQualityProfile(fnFs) ~250 bp
#plotQualityProfile(fnRs) ~200 bp

# Filter and Trim
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen = c(260, 210),
                     trimLeft = c(20, 20),
                     maxN = 0, 
                     maxEE = c(2, 4), 
                     truncQ = 2, 
                     rm.phix = TRUE,
                     compress = TRUE, 
                     multithread=TRUE)

# Error Rates -------------------------------------------------------------
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Sample Inference --------------------------------------------------------
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# Merge Paired Reads ------------------------------------------------------
merge_reads <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(merge_reads)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                    method = "consensus", 
                                    multithread = TRUE, 
                                    verbose = TRUE)

# Add ASV -----------------------------------------------------------------
asv_seqs <- colnames(seqtab.nochim)

# Assign ASV labels
asv_ids <- sprintf("ASV%04d", seq_along(asv_seqs))

# Set as new column names (for output tables)
colnames(seqtab.nochim) <- asv_ids

asv_seqs_df <- data.frame(ASV = asv_ids, Sequence = asv_seqs)

# Export ------------------------------------------------------------------
# ASV table
write_rds(seqtab.nochim, here("data_input", "seqtab.rds"))
write.csv(asv_seqs_df, here("data_input", "ASV_sequences.csv"), row.names=FALSE)

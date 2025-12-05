<p align="left">

<!-- Zenodo DOI badge -->
<a href="https://doi.org/INSERT_DOI_HERE">
  <img src="https://img.shields.io/badge/DOI-pending-blue.svg" alt="DOI">
</a>

<!-- License badge -->
<a href="https://github.com/USERNAME/REPO/blob/main/LICENSE">
  <img src="https://img.shields.io/badge/license-MIT-green.svg" alt="License: MIT">
</a>

<!-- Built with R -->
<img src="https://img.shields.io/badge/R-%3E%3D%204.3.0-blue.svg" alt="R >=4.3">

</p>

# Comparative genomics of epiphytic and endophytic bacterial culture collections from *Arabidopsis thaliana* in Ōtautahi (Christchurch), Aotearoa New Zealand
 
## Overview

This repository contains the R codes and analysis objects associated with the AC Collection, a curated phyllosphere bacterial culture collection isolated from *Arabidopsis thaliana* plants sampled in Christchurch, New Zealand.

The repository includes scripts for:
- Community ecology analyses
- Genome annotation parsing
- Phylogenetic reconstruction
- ANI clustering
- Genomica and functional annotation exploration (KEGG, CAZy, plasmids, phages)

All large datasets to reproduce the findings of our publication (colony images, isolate metadata, processed .rds objects) are stored in the associated Zenodo repository.
Genome assemblies and raw 16S rRNA gene amplicon sequencing data are available in the EMBL-ENI ENA BioProject accession PRJEB98369.

## Repository Structure
 ```
📂 data_input/                # Must be created locally and include the datasets from the Zenodo repository
📂 output/                    # Must be created locally to store output objects (plots)
📂 src/                       # Core analysis scripts and helpers
  📄 00_ASV_inference.R                # ASV inference and denoising
  📄 01_ASV_classification.R           # ASV taxonomy assignment and formatting
  📄 02_AC_taxonomy.R                  # Classification of isolates and GTDB/ANI integration
  📄 03_Colour_palette.R               # Centralised plotting styles and colour definitions
  📄 04_Community_analysis.R           # Community ecology and diversity analysis
  📄 05_Community_analysis_Figure1.R   # Figure 1 generation: community composition and diversity
  📄 06_ANI.R                          # ANI clustering and genome similarity workflows
  📄 07_AC_phylogeny_Figure2.R         # Figure 2: phylogenetic reconstruction and annotation
  📄 08_Taxa_distribution_Figure3.R    # Figure 3: taxonomic patterns and lineage distribution
  📄 09_Genome_structure_Figure4.R     # Figure 4: genome structural traits and statistics
  📄 10_Functional_analysis_Figure5.R  # Figure 5: functional trait analysis and statistics
  📄 helper_functions.R                # Utility functions shared across scripts
```

## Workflow
1. Clone this repository
```
git clone https://github.com/<username>/AC_collection.git
cd AC_collection
```

2. Create an input and output directory
```
mkdir -p data_input
mkdir -p output
```

3. Download the metadata from Zenodo
Visit the Zenodo repository: 10.5281/zenodo.17812155. Download the dataset archive and extract its contents into the ```data_input/``` directory.

4. Open RStudio and run the workflow

Execute scripts sequentially:
```
source("src/00_ASV_inference.R")
source("src/01_ASV_classification.R")
source("src/02_AC_taxonomy.R")
source("src/03_Colour_palette.R")
source("src/04_Community_analysis.R")
source("src/05_Community_analysis_Figure1.R")
source("src/06_ANI.R")
source("src/07_AC_phylogeny_Figure2.R")
source("src/08_Taxa_distribution_Figure3.R")
source("src/09_Genome_structure_Figure4.R")
source("src/10_Functional_analysis_Figure5.R")
```

Figures and processed outputs will be saved automatically into the ```output/``` directory.

## Reference

If you use this repository, please cite:

> Add Publication link

> Add Zenodo link

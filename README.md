

[![DOI](https://zenodo.org/badge/304080502.svg)](https://zenodo.org/badge/latestdoi/304080502)



# metabarcoding-primer-portfolio
Evaluation of 22 primer sets for metabarcoding marine and freshwater fishes, mostly

These are the analyses to accompany a manuscript testing the combined taxonomic recovery and resolution of a portfolio of primer pairs targeting short fragments (<300 bp) in multiple barcoding genes.


## Data

Data include the output from the initial bioinformatic workflow using cutadapt, dada2, and a blast search of a custom database for metazoa. These files are in the `/data/` directory. The primary data from which the analyses begin are:

1) Amplicon Sequence Variant (ASV) count tables - an R data object that contains the read counts per sample per ASV for each primer set.

2) ASV taxonomy outputs - a list of BLAST hits from NCBI for each FASTA sequence associated with each ASV, also includes the primer set that generated the ASV. 


## Analyses

Analyses are, for the most part, sequential, and numbered in the appropriate order to generate the dependencies for the next set of analyses.

Once the ASV tables have been decontaminated, the taxonomic hits from BLAST must be filtered before the data can be combined and concordance with reference pools assessed.

01-taxonomy-filter-BLAST-hits

02-read-counts-summary 

03-species-occupancy-detection-modeling 

04-filter-ASV-by-SODM 

05-Bray-Curtis-reference-DNA-pools

06-mock-feeds-SODM-Bray-Curtis-filters 

07-full-reference-DNA-pool-primer-evaluation

08-integrate-taxonomic-levels-within-loci

09-vouchered-false-positives

10-mock-feeds-analyses



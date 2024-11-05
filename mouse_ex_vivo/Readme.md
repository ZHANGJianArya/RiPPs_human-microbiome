
#### Metagenomic analysis
Raw reads from metagenomic sequencing were processed using Fastp v0.21.1 for quality control, with parameter of “--detect_adapter_for_pe -l 50 -5 3 -3 3”. High-quality metagenomic sequencing reads were further subjected to KneadData (https://github.com/biobakery/kneaddata) for detecting and removing reads belonging to the human genome, though searching against the mouse reference genome (GRCm39) from GENCODE87. MetaPhlAn version 4.1.1 was used to generate taxonomic profiles of metagenomes. 
Alpha diversity, as measured by Shannon diversity, and beta diversity, evaluated through Bray-Curtis dissimilarity, were computed for species-level taxonomic profiles utilizing the Vegan package in R. Differential abundance analysis between groups was conducted using MaAsLin2 (microbiome multivariable associations with linear models). Species with a false discovery rate (FDR)-adjusted p-value of less than 0.05 were deemed significantly different between the two groups.

#### Data availability
Raw read sequences of the shotgun metagenomic sequences were deposited at Sequence Read Archive (SRA, NCBI) under accession of PRJNA1166984.

# Code associated with the manuscript  
**"Microbial biosynthetic potential across global geothermal systems"**

This project explores the tectonic fingerprint on microbial secondary metabolite potential in geothermal environments along the South and Central American Volcanic Zone, Iceland and Italy. This repo contains the code behind the Docker Antismash pipeline used for automatising antismash looping through assembly FASTA files, along with the subsequent analysis and R code for figure generation. All the sequences analyzed in this study are available through ENA under the Umbrella Project CoEvolve PRJEB55081.

## This repo consists of:

1. **antiSMASH Pipeline folder**:  
   Contains the `raw_fasta_to_csv_antiSmash.ps1` script that recursively processes FASTA (.fa) files, runs antiSMASH via Docker for BGC detection, organises outputs into structured directories, and post-processes results by grouping `.gbk` files and running BiG-SCAPE for clustering and network analysis.

   - **SQL_for_Bigscape folder**:  
     Contains SQL scripts used to query the BiG-SCAPE SQLite database, generating sample-by-GCF presence/absence matrices and retrieving GCF product annotations. These outputs serve as the core input datasets for downstream analyses.

3. **Figures folder**:  
   Contains all the individual figures used in manuscript.

4. **Notebook folder**:  
   Contains the code to reproduce the analysis and figures in the manuscript, along with other exploratory plots. The folder also holds the core datasets used to build the phyloseq object for microbial community analysis: env_data_BGC (environmental metadata per sample), OTU_BGC (GCF OTU table), tax_BGC (taxonomy assignments for BiG-SCAPE GCFs), and the phyloseq object creation allowed for the analysis of BGC potential across the sampled regions.

5. **additional-analyses folder**:  
   Contains the subfolders for:  
   - **network_analysis**: with Cytoscape and network files for the four main BGC classes (TERPENE, PKS, RiPP, and NRPS), together with Python scripts for novelty analysis, MIBiG distance summarisation, conversion of network outputs for Gephi visualisation, and generation of cluster-level box plots
   - **CAT**: with files related to CAT-based taxonomic classification   
   - **phylogenetic_analysis**: with files and scripts for phylogenetic reconstruction  
   

## Please cite as:

Ana Clara Pelliciari Silva, Benoit de Pins, Francesco Montemagno, Flavia Migliaccio, Martina Cascone, Deborah Bastoni, Bernardo Barosa, Matteo Selci, Costantino Vetriani, Agostina Chiodi, Federico A. Vignale, Maria Garcia Alai, Alberto Vitale Brovarone, Gerdhard L. Jessen, Jenny M. Blamey, J. Maarten de Moor, Karen G. Lloyd, Peter Barry, Donato Giovannelli.  
**2025.** *Microbial biosynthetic potential across global geothermal systems.* Preprint

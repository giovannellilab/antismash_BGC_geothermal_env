# Code associated with the manuscript "Tectonic regime shapes microbial biosynthetic potential across global geothermal systems"

This project explores the tectonic fingerprint on microbial secondary metabolite potential in geothermal environments along the South and Central American Volcanic Zone, Iceland and Italy. This repo contains the code behind the Docker Antismash pipeline used for automatising antismash looping through assembly FASTA files, along with the subsequent analysis and R code for figure generation. All the sequences analyzed in this study are available through ENA under the Umbrella Project CoEvolve PRJEB55081.

## This repo consists of:

- antiSmash Pipeline folder : contains the 1) raw_fasta_to_csv_antiSmash.ps1 script that loops through FASTA files, uses a Docker container running antiSMASH to identify biosynthetic gene clusters (BGCs), and saves the results. 2) json_to_csv.py: A Python script that processes antiSMASH JSON outputs and extracts key information (node ID, product, and category) into CSV format for further analysis.
- Figures folder: contains all the individual figures used in manuscript.
- Notebook folder: contains the code to reproduce the analysis and figures in the manuscript, along with other exploratory plots. The folder also holds the core datasets used to build the phyloseq object for microbial community analysis: 1) env_data_BGC: environmental metadata per sample, 2) otu_REAL: OTU abundance table and 3) BGCs_tax_table: taxonomy assignments for antismash BGCs. The phyloseq object creation allowed for the analysis of BGC potential across the sampled regions.
- Resources: contains .json files with combined BGC data (e.g., cluster categories and sequence sizes) per site, used as inputs for the plots and statistical analysis in the notebook.

**Please cite as:**

Ana Clara Pelliciari Silva, Benoit de Pins, Agostina Chiodi, Gerdhard Jessen, J. Maarten de Moor, Karen Lloyd, Peter Barry, and Donato Giovannelli. 2025. Tectonic regime shapes microbial biosynthetic potential across global geothermal systems. _Preprint_

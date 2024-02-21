# Extended duration of gut microbiota depletion reverses protection from experimental autoimmune uveitis
This repository contains scripts/code used for analyzing the microbiome and nanostring data in this study.

- Authors: Vijay Nagarajan PhD
- Affiliation: Laboratory of Immunology, NEI/NIH
- Contact: nagarajanv@nih.gov
- Platform: The workflow BASH commands were developed to run on a Mac OS, but could be reused/reproduced in any unix environment with appropriate changes
- Publication: [Salvador R, Horai R, Zhang A, Jittayasothorn Y, Tang J, Gupta A, Nagarajan V, Caspi RR. Too Much of a Good Thing: Extended Duration of Gut Microbiota Depletion Reverses Protection From Experimental Autoimmune Uveitis. Invest Ophthalmol Vis Sci. 2023 Nov 1;64(14):43. doi: 10.1167/iovs.64.14.43. PMID: 38019490; PMCID: PMC10691388.](https://pubmed.ncbi.nlm.nih.gov/38019490/)
-------------------------------
1. [Tools/nsolver_extract.sh](Tools/nsolver_extract.sh) This bash script takes nanostring/nsolver DE file and extracts filtered DGE lists, filtered DGE names, ids, cleans ids, writes file for geo matrix and cytoscape input
2. [Tools/tabtogmt.R](Tools/tabtogmt.R) This R script reads in a tab delimited gene set file and generates the gmt format input file for gsea.
3. [Tools/gseaonly.sh](Tools/gseaonly.sh) This bash script reads in the gsea input, gmt file and generates gsea output for all samples in a folder
4. [Tools/enrichmentPlotsGsea.R](Tools/enrichmentPlotsGsea.R) This R script reads lists of GSEA generated enrichment TSV files and generates a plot with all the results in one panel.

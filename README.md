# Nanostring data analysis tools and methods
This repository contains scripts/code used for analyzing the nanostring data.
1. [Tools/nsolver_extract.sh](Tools/nsolver_extract.sh) This bash script takes nanostring/nsolver DE file and extracts filtered DGE lists, filtered DGE names, ids, cleans ids, writes file for geo matrix and cytoscape input
2. [Tools/enrichmentPlotsGsea.R](Tools/enrichmentPlotsGsea.R) This R script reads lists of GSEA generated enrichment TSV files and generates a plot with all the results in one panel.
3. [Tools/tabtogmt.R](Tools/tabtogmt.R) This R script reads in a tab delimited gene set file and generates the gmt format input file for gsea.

# Extended duration of gut microbiota depletion reverses protection from experimental autoimmune uveitis

This repository contains scripts/code used for analyzing the microbiome and nanostring data in this study.
1. [Tools/nsolver_extract.sh](Tools/nsolver_extract.sh) This bash script takes nanostring/nsolver DE file and extracts filtered DGE lists, filtered DGE names, ids, cleans ids, writes file for geo matrix and cytoscape input
2. [Tools/tabtogmt.R](Tools/tabtogmt.R) This R script reads in a tab delimited gene set file and generates the gmt format input file for gsea.
3. [Tools/gseaonly.sh](Tools/gseaonly.sh) This bash script reads in the gsea input, gmt file and generates gsea output for all samples in a folder
4. [Tools/enrichmentPlotsGsea.R](Tools/enrichmentPlotsGsea.R) This R script reads lists of GSEA generated enrichment TSV files and generates a plot with all the results in one panel.

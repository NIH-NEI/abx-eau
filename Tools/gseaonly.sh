# Author : Vijay Nagarajan PhD
# Institute : NEI/NIH
# This scripts reads in a tab delimited gene set file and generates the gmt format input file for gsea

# Sample list
samples="IEL_LT_d0.vs..PLN_NT_d0	IEL_ST_d0.vs..PLN_NT_d0	SILP_ST_d0.vs..PLN_NT_d0	SILP_ST_d7.vs..PLN_NT_d0"

# Input directory
datadirectory="/Users/gseainput/"

# Gmt directory
gmtdirectory="/Users/gmtfile/"

# Output directory
gseaoutput="/Users/results/GSEA_Results/"

# Run gsea for each sample
for sample in $samples;
  do
    echo ${sample}
    bash /Users/GSEA_4.1.0/gsea-cli.sh GSEAPreranked -gmx "${gmtdirectory}tcell_annotated_abx.gmt" -rnk "${datadirectory}${sample}_for_gsea.rnk" -out ../GSEA_Results -include_only_symbols true -rpt_label ${sample} -collapse No_Collapse -set_min 5
  done




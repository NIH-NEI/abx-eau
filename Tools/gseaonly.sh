# Author : Vijay Nagarajan PhD
# Institute : NEI/NIH
# This bash scripts reads in the gsea input, gmt file and generates gsea output for all samples in a folder

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




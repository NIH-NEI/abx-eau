#!/bin/bash
# Author : Vijay Nagarajan PhD
# Institute : NEI/NIH
# This bash script takes nanostring/nsolver DE file and extracts filtered DGE lists, filtered DGE names, ids, cleans ids, writes file for geo matrix and cytoscape input

samples="LiLP LP MLN SiLP SLN"

for sample in $samples;
do

# Set file/folder paths
datadirectory="/Users/AdvAnalysisReport_${sample}_PBS/results/DE"
resultspath="/Users/myresults"
TAB=$'\t'

# Change to data directory
cd "${datadirectory}"

# Covert csv to tsv
sed -E 's/("([^"]*)")?,/\2'"${TAB}"'/g' DE\ results\ -\ treatmentLAQ.csv | sed 's/"//g' > DE\ results\ -\ treatmentLAQ.txt

# Extract significantly differentially expressed genes and save files in results directory
awk -F"\t" '{ if(($9 <= 0.05) && ($2 <= -1)) { print }}' < DE\ results\ -\ treatmentLAQ.txt > "${resultspath}${sample}_sde_down.txt"
awk -F"\t" '{ if(($9 <= 0.05) && ($2 >= 1)) { print }}' < DE\ results\ -\ treatmentLAQ.txt > "${resultspath}${sample}_sde_up.txt"

# Extract geneids and genenames
awk -F"\t" '{ if(($9 <= 0.05) && ($2 <= -1)) { print }}' < DE\ results\ -\ treatmentLAQ.txt | cut -f13 | cut -d"." -f1 > "${resultspath}${sample}_downgeneids.txt"
awk -F"\t" '{ if(($9 <= 0.05) && ($2 <= -1)) { print }}' < DE\ results\ -\ treatmentLAQ.txt | awk -F"-mRNA" '{print $1}' > "${resultspath}${sample}_downgenenames.txt"
awk -F"\t" '{ if(($9 <= 0.05) && ($2 >= 1)) { print }}' < DE\ results\ -\ treatmentLAQ.txt | cut -f13 | cut -d"." -f1 > "${resultspath}${sample}_upgeneids.txt"
awk -F"\t" '{ if(($9 <= 0.05) && ($2 >= 1)) { print }}' < DE\ results\ -\ treatmentLAQ.txt | awk -F"-mRNA" '{print $1}' > "${resultspath}${sample}_upgenenames.txt"

# Combine up/down gene names and gene ids for each samples
cat "${resultspath}${sample}_downgeneids.txt" "${resultspath}${sample}_upgeneids.txt" > "${resultspath}${sample}_up_down_geneids.txt"
cat "${resultspath}${sample}_downgenenames.txt" "${resultspath}${sample}_upgenenames.txt" > "${resultspath}${sample}_up_down_genenames.txt"

cut -f2-14 "${resultspath}${sample}_sde_down.txt" > "${resultspath}tempdown"
cut -f2-14 "${resultspath}${sample}_sde_up.txt" > "${resultspath}tempup"

# Paste cleaned genename columns to data
paste "${resultspath}${sample}_downgenenames.txt" "${resultspath}tempdown" > "${resultspath}${sample}_sde_down.txt"
paste "${resultspath}${sample}_upgenenames.txt" "${resultspath}tempup" > "${resultspath}${sample}_sde_up.txt"

# Add sample name column to the results files
sed 's/^/'"${sample}"''"${TAB}"'/' "${resultspath}${sample}_sde_down.txt" > "${resultspath}testdown"
sed 's/^/'"${sample}"''"${TAB}"'/' "${resultspath}${sample}_sde_up.txt" > "${resultspath}testup"

# Paste all together
paste "${resultspath}${sample}_downgeneids.txt" "${resultspath}testdown" >> "${resultspath}LAQdegenescombinedfornetwork.txt"
paste "${resultspath}${sample}_upgeneids.txt" "${resultspath}testup" >> "${resultspath}LAQdegenescombinedfornetwork.txt"

# Cleanup
rm "${resultspath}tempdown"
rm "${resultspath}tempup"
rm "${resultspath}testdown"
rm "${resultspath}testup"

done

##### The below section of code generates log2normalized sample counts files, for up/down genes
# Set directory paths and make sure headers are correct/correspond to your samples in the datasets

for sample in $samples;
do

#sample="LP"
#samples="LiLP LP MLN SiLP SLN"
unset downgenenamefile
unset upgenenamefile

# Set file/folder paths
datadirectory="/Users/AdvAnalysisReport_${sample}_PBS/results/Normalization"
resultspath="/Users/myresults/"

MLN_header="Gene\tMLN_LAQ_1\tMLN_LAQ_2\tMLN_LAQ_3\tMLN_PBS_1\tMLN_PBS_2\tMLN_PBS_3"
LiLP_header="Gene\tLiLP_LAQ_1\tLiLP_LAQ_2\tLiLP_LAQ_3\tLiLP_PBS_1\tLiLP_PBS_2\tLiLP_PBS_3"
LP_header="Gene\tSiLP_LAQ_1\tSiLP_LAQ_2\tSiLP_LAQ_3\tSiLP_PBS_1\tSiLP_PBS_2\tSiLP_PBS_3\tLiLP_LAQ_1\tLiLP_LAQ_2\tLiLP_LAQ_3\tLiLP_PBS_1\tLiLP_PBS_2\tLiLP_PBS_3"
SiLP_header="Gene\tSiLP_LAQ_1\tSiLP_LAQ_2\tSiLP_LAQ_3\tSiLP_PBS_1\tSiLP_PBS_2\tSiLP_PBS_3"
SLN_header="Gene\tSLN_LAQ_1\tSLN_LAQ_2\tSLN_LAQ_3\tSLN_PBS_1\tSLN_PBS_2\tSLN_PBS_3"

TAB=$'\t'

# Change to data directory
cd "${datadirectory}"
#ls

# Covert csv to tsv
sed -E 's/("([^"]*)")?,/\2'"${TAB}"'/g' ALL_normalized_data.csv | sed 's/"//g' > ALL_normalized_data.txt
datamash --no-strict transpose < ALL_normalized_data.txt > "${resultspath}normtransposedtemp.txt"

# Block comment, if needed for troubleshooting
#: <<'END'

if [ $sample = LiLP ]
	then
	#echo -e $LiLP_header
	# Extract log2normalized data columns up/down heatmap
	cut -f2-8 "${resultspath}normtransposedtemp.txt" > "${resultspath}normtempcut"
	# Extract log2normalized data columns for combining ALL Samples ALL genes data
	cut -f2-8 "${resultspath}normtransposedtemp.txt" > "${resultspath}ALL_log2normalized_data-lilp.txt"
	downgenenamefile="${resultspath}${sample}_downgenenames.txt"
	upgenenamefile="${resultspath}${sample}_upgenenames.txt"

	# Extract up/down genes log2normalized data
	for updowngene in $(cat "$downgenenamefile" "$upgenenamefile"); 
		do
		grep -w "$updowngene" "${resultspath}normtempcut"
  		grep -w "$updowngene" "${resultspath}normtempcut" | cut -f2-7 >> "${resultspath}tempupdown"
		done
	# Add cleaned gene name column to extracted data
	cat "$downgenenamefile" "$upgenenamefile" | paste - "${resultspath}tempupdown" > "${resultspath}tempupdowngenenames"
	# Add column header to extracted data
	echo -e "$LiLP_header" | cat - "${resultspath}tempupdowngenenames"
	echo -e "$LiLP_header" | cat - "${resultspath}tempupdowngenenames" > "${resultspath}${sample}_log2normalized_for_heatmap.txt"

	# Sample LP
	elif [ $sample = LP ]
	then
	#echo -e $LP_header
	cut -f2-14 "${resultspath}normtransposedtemp.txt" > "${resultspath}normtempcut"	
	downgenenamefile="${resultspath}${sample}_downgenenames.txt"
	upgenenamefile="${resultspath}${sample}_upgenenames.txt"
	for updowngene in $(cat "$downgenenamefile" "$upgenenamefile"); 
		do
		grep -w "$updowngene" "${resultspath}normtempcut"
  		grep -w "$updowngene" "${resultspath}normtempcut" | cut -f2-13 >> "${resultspath}tempupdown"
		done
	cat "$downgenenamefile" "$upgenenamefile" | paste - "${resultspath}tempupdown" > "${resultspath}tempupdowngenenames"
	echo -e "$LP_header" | cat - "${resultspath}tempupdowngenenames"
	echo -e "$LP_header" | cat - "${resultspath}tempupdowngenenames" > "${resultspath}${sample}_log2normalized_for_heatmap.txt"
	
	# Sample MLN	
	elif [ $sample = MLN ]
	then
	#echo -e $MLN_header
	cut -f2-8 "${resultspath}normtransposedtemp.txt" > "${resultspath}normtempcut"
	cut -f3-8 "${resultspath}normtransposedtemp.txt" | paste "${resultspath}ALL_log2normalized_data-lilp.txt" - > "${resultspath}ALL_log2normalized_data-lilp-mln.txt"
	downgenenamefile="${resultspath}${sample}_downgenenames.txt"
	upgenenamefile="${resultspath}${sample}_upgenenames.txt"
	for updowngene in $(cat "$downgenenamefile" "$upgenenamefile"); 
		do
		grep -w "$updowngene" "${resultspath}normtempcut"
  		grep -w "$updowngene" "${resultspath}normtempcut" | cut -f2-7 >> "${resultspath}tempupdown"
		done
	cat "$downgenenamefile" "$upgenenamefile" | paste - "${resultspath}tempupdown" > "${resultspath}tempupdowngenenames"
	echo -e "$MLN_header" | cat - "${resultspath}tempupdowngenenames"
	echo -e "$MLN_header" | cat - "${resultspath}tempupdowngenenames" > "${resultspath}${sample}_log2normalized_for_heatmap.txt"

	# Sample SiLP
	elif [ $sample = SiLP ]
	then
	#echo -e $SiLP_header
	cut -f2-8 "${resultspath}normtransposedtemp.txt" > "${resultspath}normtempcut"
	cut -f3-8 "${resultspath}normtransposedtemp.txt" | paste "${resultspath}ALL_log2normalized_data-lilp-mln.txt" - > "${resultspath}ALL_log2normalized_data-lilp-mln-silp.txt"
	downgenenamefile="${resultspath}${sample}_downgenenames.txt"
	upgenenamefile="${resultspath}${sample}_upgenenames.txt"
	for updowngene in $(cat "$downgenenamefile" "$upgenenamefile"); 
		do
		grep -w "$updowngene" "${resultspath}normtempcut"
  		grep -w "$updowngene" "${resultspath}normtempcut" | cut -f2-7 >> "${resultspath}tempupdown"
		done
	cat "$downgenenamefile" "$upgenenamefile" | paste - "${resultspath}tempupdown" > "${resultspath}tempupdowngenenames"
	echo -e "$SiLP_header" | cat - "${resultspath}tempupdowngenenames"
	echo -e "$SiLP_header" | cat - "${resultspath}tempupdowngenenames" > "${resultspath}${sample}_log2normalized_for_heatmap.txt"
	
	# Sample SLN
	else 
	#echo -e $SLN_header
	cut -f2-8 "${resultspath}normtransposedtemp.txt" > "${resultspath}normtempcut"
	cut -f3-8 "${resultspath}normtransposedtemp.txt" | paste "${resultspath}ALL_log2normalized_data-lilp-mln-silp.txt" - > "${resultspath}ALL_log2normalized_data-lilp-mln-silp-sln.txt"
	downgenenamefile="${resultspath}${sample}_downgenenames.txt"
	upgenenamefile="${resultspath}${sample}_upgenenames.txt"
	for updowngene in $(cat "$downgenenamefile" "$upgenenamefile"); 
		do
		grep -w "$updowngene" "${resultspath}normtempcut"
  		grep -w "$updowngene" "${resultspath}normtempcut" | cut -f2-7 >> "${resultspath}tempupdown"
		done
	cat "$downgenenamefile" "$upgenenamefile" | paste - "${resultspath}tempupdown" > "${resultspath}tempupdowngenenames"
	echo -e "$SLN_header" | cat - "${resultspath}tempupdowngenenames"
	echo -e "$SLN_header" | cat - "${resultspath}tempupdowngenenames" > "${resultspath}${sample}_log2normalized_for_heatmap.txt"

fi

#END

# Cleanup
rm "${resultspath}tempupdowngenenames"	
rm "${resultspath}tempupdown"
rm "${resultspath}normtempcut"
rm "${resultspath}normtransposedtemp.txt"

done

rm "${resultspath}ALL_log2normalized_data-lilp.txt"
rm "${resultspath}ALL_log2normalized_data-lilp-mln.txt"
rm "${resultspath}ALL_log2normalized_data-lilp-mln-silp.txt"

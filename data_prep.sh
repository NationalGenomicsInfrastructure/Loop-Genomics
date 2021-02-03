#!/bin/bash

#------------Input Parameter-------------
file_name=/lupus/home/husain/Projects/LoopGenomics/ngi_projects/P14754/temp/counting/P14754_CGCTCATT_sample6_counting.tar.gz
output_file_name=/lupus/home/husain/Projects/LoopGenomics/ngi_projects/source_code/fusion_gene/data/P14754_CGCTCATT_sample6_counting

#------Extract blast files only----
#echo "Extarcting blast files...."
#mkdir blast_files_only
#tar -xzf $file_name -C blast_files_only --wildcards --no-anchored '*.blst'

#-----------------merge all blast files into one file---------
echo "Merge all blast files into a text file...."
find blast_files_only*/  -name "*.blst" -exec cat '{}' ';' > ${output_file_name}_all_blast.txt

#----------------save blast text files into RData format with their column names-------------
echo "Save text file into RData format with headers....."
Rscript blast_files_save_in_RData.R ${output_file_name}_all_blast

#----------Remove temporary folder blast_files_only-----------
echo "Remove temporary files....."
rm -r blast_files_only


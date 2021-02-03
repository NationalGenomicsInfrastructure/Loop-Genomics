# Loop-Genomics

Step1. Load and Prep data 
  1a. Use either Load_blast_files.R OR
  1b. data_prep.sh (this will internally call blast_files_save_in_RData.R )
Step2. Update blast files with genes and paired with shared molecules (pre_fusion.R)
Step3. Detect Fusion Genes (fusion.R)
Step4. For Visualization (fusion_plot.R)

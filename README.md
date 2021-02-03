# Loop-Genomics: Detect fusion genes

1. Step1: Load and Prep data
    1. Use either `Load_blast_files.R` OR
    2. `data_prep.sh` (this will internally call `blast_files_save_in_RData.R` )
2. Step2: Update blast files with genes and paired with shared molecules (`pre_fusion.R`)
3. Step3: Detect Fusion Genes (`fusion.R`)
4. Step4: For Visualization (`fusion_plot.R`)

# proteins_MDD
#(1) MR analysis
     # R_Computer_code: "Primary_MR_analysis.R"
     # Primary analysis
       # exposure: "CSF_protein_cis.xlsx"/ "Plasma_protein_cis.xlsx"
       # outcome: "MDD_SNPs.xlsx"
       # harmonise exposure and outcome data: "CSF_harmonise_data.xlsx"/ "Plasma_harmonise_data.xlsx"
     # External validation 1
       # outcome: "FinnGen_SNPs_1.xlsx"
     # External validation 2
       # exposure: "External_validation_2_proteins_cis.xlsx"
       # outcome: "MDD_SNPs_2.xlsx"/ "FinnGen_SNPs_2.xlsx"
#(2) Spearman’s rank correlation analysis
     # R_Computer_code: "Spearman.R"  
#(3) Bayesian colocalization analysis
     # R_Computer_code: "Bayesian_colocalization_analysis.R"  
     # GRch37->GRch38: "hg19ToHg38.over.chain"
#(4) SMR
     # software: "smr-1.3.1-win.zip"
     # Convert to BESD file: "BESD.R"
     # Preparation of the outcome document: "SMR_outcome.R"
     # SMR command: "SMR_command.txt"

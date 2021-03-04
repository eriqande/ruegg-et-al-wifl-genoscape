
# this first rule just ensures that all the RMarkdown files get run
rule all:
  input:
    "001-rad-seq-data-summaries-and-pca.html",
    "002-select-snps-for-assays-from-rad.html",
    "003-process-5-plates-of-birds-at-192-fluidigm-snps.html",
    "004-combine-RAD-and-fludigm-breeders-and-run-STRUCTURE.html",
    "005-choosing-96-SNPs-from-amongst-the-179.html",
    "006-make-the-genoscape.html",
    "007-process-fluidigm-plates-for-wintering-birds.html",
    "008-rubias-self-assign-and-wintering-birds.html"
    

rule stored_results_from_cluster:
  input:
  output:
    "stored_results/002/candidate_snps_tossed.csv",
    "stored_results/004/StructOuput_genos_slg_pipe.txt_dat001_k005_Rep001.txt_f",
    
    
rule input_data_files_etc:
  input:
  output:
    "data/rad_wifl_clean_175_105000.rds",
    "data/appendices/app1-SITable1_Breed_IndvAssigwgenos.csv",
    "data/rad_wifl_less_filtered_219_349014.rds",
    "data/wifl-snp-short-names.tsv",
    "data/wifl_snp_comparison_counts.txt",
    "data/fluidigm_panel_1_and_2_chips",
    "data/WIFL_fluidigm_assay_info.csv",
    "data/fluidigm_panel_3_chips",
    "data/fluidigm_additional_chips/WIFL.Plate16_PopID_CBdetailedResults.csv",
    "inputs/004/slg-pipe-locus-order.rds"


rule rad_seq_data_summaries_and_pca_001:
  input:
    "001-rad-seq-data-summaries-and-pca.Rmd",
    "data/rad_wifl_clean_175_105000.rds",
    "data/appendices/app1-SITable1_Breed_IndvAssigwgenos.csv",
  output:
    "001-rad-seq-data-summaries-and-pca.html",
    "outputs/001/rad-fst-table.csv",
    "outputs/001/rad-fst-sample-sizes-table.csv"
  shell:
    "Rscript --vanilla R/render-markdown.R {input[0]}"



rule select_assay_snps_from_rad_seq_002:
  input:
    "002-select-snps-for-assays-from-rad.Rmd",
    "data/rad_wifl_clean_175_105000.rds",
    "data/appendices/app1-SITable1_Breed_IndvAssigwgenos.csv",
    "data/rad_wifl_less_filtered_219_349014.rds",
    "stored_results/002/candidate_snps_tossed.csv",
    "data/wifl-snp-short-names.tsv",
    "data/wifl_snp_comparison_counts.txt"
  output:
    "002-select-snps-for-assays-from-rad.html",
    "outputs/002/top40_lfmm.csv",
    "outputs/002/candidate_snps_untossed.csv",
    "outputs/002/regions-for-candidates.txt",
    "outputs/002/relevant_variation.txt",
    "outputs/002/consensus_sequences_of_snps.txt",
    "outputs/002/wifl-assay-candidates-with-ranks.csv"
  shell:
    "Rscript --vanilla R/render-markdown.R {input[0]}"
    


rule process_5_plates_of_birds_at_192_fluidigm_snps_003:
  input:
    "003-process-5-plates-of-birds-at-192-fluidigm-snps.Rmd",
    "data/fluidigm_panel_1_and_2_chips"
  output:
    "003-process-5-plates-of-birds-at-192-fluidigm-snps.html",
    "outputs/003/564_birds_at_192_fluidigm_snps.rds",
    "outputs/003/final_fluidigm_breeders_179_loci_393_birds_long.rds"
  shell:
    "Rscript --vanilla R/render-markdown.R {input[0]}"


rule combine_RAD_and_fludigm_breeders_and_run_STRUCTURE_004:
  input:
    "004-combine-RAD-and-fludigm-breeders-and-run-STRUCTURE.Rmd",
    "data/rad_wifl_clean_175_105000.rds",
    "outputs/003/final_fluidigm_breeders_179_loci_393_birds_long.rds",
    "data/WIFL_fluidigm_assay_info.csv",
    "inputs/004/slg-pipe-locus-order.rds"

  output:
    "004-combine-RAD-and-fludigm-breeders-and-run-STRUCTURE.html",
    "outputs/004/rad_and_fluidigm_combined_568_breeders_at_179_snps.rds"
  shell:
    "Rscript --vanilla R/render-markdown.R {input[0]}"
  


rule choosing_96_SNPs_from_amongst_the_179_005:
  input:
    "005-choosing-96-SNPs-from-amongst-the-179.Rmd",
    "stored_results/004/StructOuput_genos_slg_pipe.txt_dat001_k005_Rep001.txt_f",
    "inputs/004/slg-pipe-locus-order.rds"
  output:
    "005-choosing-96-SNPs-from-amongst-the-179.html",
    "outputs/005/top96_loci.csv"
  shell:
    "Rscript --vanilla R/render-markdown.R {input[0]}"
  
rule make_the_genoscape_006:
  input:
    "006-make-the-genoscape.Rmd",
    "stored_results/004/Struct-q-values-k007_Rep004.txt",
    "data/appendices/app1-SITable1_Breed_IndvAssigwgenos.csv",
    "data/WIFL_RangeMap/WIFLrev.shp",
  output:
    "006-make-the-genoscape.html"
  shell:
    "Rscript --vanilla R/render-markdown.R {input[0]}"
  

rule process_fluidigm_plates_for_wintering_birds_007:
  input:
    "007-process-fluidigm-plates-for-wintering-birds.Rmd",
    "data/fluidigm_panel_3_chips",
    "data/fluidigm_additional_chips/WIFL.Plate16_PopID_CBdetailedResults.csv",
    "data/WIFL_fluidigm_assay_info.csv",
    "outputs/003/564_birds_at_192_fluidigm_snps.rds",
    "data/appendices/app1-winter-indivs-assignments-and-genos.csv"
  output:
    "007-process-fluidigm-plates-for-wintering-birds.html",
    "outputs/007/470_birds_at_96_panel3_fluidigm_snps.rds",
    "outputs/007/winter_rubias.rds"
  shell:
    "Rscript --vanilla R/render-markdown.R {input[0]}"

rule rubias_self_assign_and_wintering_birds_008:
  input:
    "008-rubias-self-assign-and-wintering-birds.Rmd",
    "data/appendices/app1-SITable1_Breed_IndvAssigwgenos.csv",
    "data/appendices/app1-winter-indivs-assignments-and-genos.csv",
    "outputs/007/winter_rubias.rds",
    "outputs/004/rad_and_fluidigm_combined_568_breeders_at_179_snps.rds"
  output:
    "008-rubias-self-assign-and-wintering-birds.html"
  shell:
    "Rscript --vanilla R/render-markdown.R {input[0]}"



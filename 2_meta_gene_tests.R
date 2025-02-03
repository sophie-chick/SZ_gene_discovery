# setup 
library(tidyverse) 
library(metap) 

# sample numbers 
n_case = 4650 
n_control = 5719 
n_male_case = 3329 
n_male_control = 2729 
n_female_case = 1321 
n_female_control = 2990 

# importing combined strata-level table and splitting into variant classes 
combined = read.delim("/scratch/c.c21070635/2_meta-analysis/output/3_meta-analysis/2_meta_rvas/1_combined_meta_gene_table.tsv") 
meta = list() 
meta[["ptv"]] = combined %>%
    filter(variant_class == "PTV") 
meta[["ptv_mpc3"]] = combined %>%
    filter(variant_class == "PTV+MPC3") 
meta[["ptv_mpc2"]] = combined %>%
    filter(variant_class == "PTV+MPC2") 
meta[["mpc2"]] = combined %>% 
    filter(variant_class == "MPC2") 

# nonpar gene list
nonpar_genes = read.delim("/scratch/c.c21070635/2_meta-analysis/data_processed/SZ_AD_nonpar_gene_list.tsv")$gene_id

# schema and mpc genes for filling in table 
schema_genes = read.delim("/scratch/c.c21070635/2_meta-analysis/data/gene_lists/SCHEMA_gene_results_paper_with_ids.tsv")
mpc2_genes = read.table("/scratch/c.c21070635/2_meta-analysis/data/gene_lists/MPC2_genes.tsv", header=TRUE) 
mpc3_genes = read.table("/scratch/c.c21070635/2_meta-analysis/data/gene_lists/MPC3_genes.tsv", header=TRUE) 

# AC list and "no" column lists ("no" = number of cases or controls not carrying variants)
# for filling in NAs after reshape
CC_list = c("EUR_exomes_case_CC", "EUR_exomes_control_CC", 
            "EUR_N_exomes_case_CC", "EUR_N_exomes_control_CC", 
            "AFR_genomes_case_CC", "AFR_genomes_control_CC", 
            "EAS_exomes_case_CC", "EAS_exomes_control_CC", 
            "AMR_exomes_case_CC", "AMR_exomes_control_CC", 
            "FIN_exomes_case_CC", "FIN_exomes_control_CC", 
            "ASJ_exomes_case_CC", "ASJ_exomes_control_CC", 
            "FIN_genomes_case_CC", "FIN_genomes_control_CC", 
            "EST_genomes_case_CC", "EST_genomes_control_CC", 
            "AFR_exomes_case_CC", "AFR_exomes_control_CC", 
            "SAS_exomes_case_CC", "SAS_exomes_control_CC",
            "in_house_case_CC", "in_house_control_CC", 
            "in_house_male_case_CC", "in_house_male_control_CC", 
            "in_house_female_case_CC", "in_house_female_control_CC") 
no_list = c("EUR_exomes_case_no", "EUR_exomes_control_no", 
            "EUR_N_exomes_case_no", "EUR_N_exomes_control_no", 
            "AFR_genomes_case_no", "AFR_genomes_control_no", 
            "EAS_exomes_case_no", "EAS_exomes_control_no", 
            "AMR_exomes_case_no", "AMR_exomes_control_no", 
            "FIN_exomes_case_no", "FIN_exomes_control_no", 
            "ASJ_exomes_case_no", "ASJ_exomes_control_no", 
            "FIN_genomes_case_no", "FIN_genomes_control_no", 
            "EST_genomes_case_no", "EST_genomes_control_no", 
            "AFR_exomes_case_no", "AFR_exomes_control_no", 
            "SAS_exomes_case_no", "SAS_exomes_control_no",
            "in_house_case_no", "in_house_control_no", 
            "in_house_male_case_no", "in_house_male_control_no", 
            "in_house_female_case_no", "in_house_female_control_no") 

# importing and reshaping
variant_list = list("ptv", "ptv_mpc3", "ptv_mpc2", "mpc2") 
for (variant in variant_list) {

  cat(paste((if (variant != "ptv") {"\n"} else {""}), "Converting ", variant, " from long to wide:\n", sep=""))
  cat(paste("Read in ", length(unique(meta[[variant]]$gene_id)), " genes\n", sep=""))

  # removing variant class column
  meta[[variant]]$variant_class = NULL 
  
  cat("Reshaping...\n")
  # preparing for long to wide conversion
  meta[[variant]]$group = gsub("-", "_", meta[[variant]]$group)
  meta[[variant]]$group = gsub(" ", "_", meta[[variant]]$group)
  meta[[variant]]$group = gsub("[(]", "", meta[[variant]]$group)
  meta[[variant]]$group = gsub("[)]", "", meta[[variant]]$group)
  
  # long to wide 
  meta[[variant]] = meta[[variant]] %>% pivot_wider( 
    id_cols = c(gene_id, gene_symbol), 
    names_from = group, 
    values_from = c(case_CC, control_CC, case_no, control_no), 
    names_sep = "_") 
  
  # changing prefixes to suffixes from reshape
  names(meta[[variant]]) = sapply(names(meta[[variant]]), function(x) 
    if (substr(x, 1, 7) == "case_CC") { 
      paste(substr(x, 9, nchar(x)), "case_CC", sep="_") } 
    else if (substr(x, 1, 10) == "control_CC") { 
      paste(substr(x, 12, nchar(x)), "control_CC", sep="_") } 
    else if (substr(x, 1, 7) == "case_no") { 
      paste(substr(x, 9, nchar(x)), "case_no", sep="_") } 
    else if (substr(x, 1, 10) == "control_no") { 
      paste(substr(x, 12, nchar(x)), "control_no", sep="_") } 
    else {x}) 
  
  # filling in genes after reshaping 
  if (variant == "ptv") { 
    meta[[variant]] = merge(schema_genes, meta[[variant]], by=c("gene_id","gene_symbol"), 
                            all.x=TRUE)}
  if (variant %in% c("ptv_mpc2","mpc2")) { 
    meta[[variant]] = merge(mpc2_genes, meta[[variant]], by=c("gene_id","gene_symbol"), 
                            all.x=TRUE)}
  if (variant == "ptv_mpc3") { 
    meta[[variant]] = merge(mpc3_genes, meta[[variant]], by=c("gene_id","gene_symbol"), 
                            all.x=TRUE)}
  cat(paste("Filling in ", length(unique(meta[[variant]]$gene_id)), " genes\n", sep=""))
  
  # filling in NAs that were introduced during reshaping/merge 
  cat("Replacing NAs and reordering columns...\n") 
  # can't filter for CC columns because might exclude in_house_male and female 
  # creating any columns that aren't there first 
  CC_cols = names(meta[[variant]])[str_detect(names(meta[[variant]]),"_CC")] 
  for (column in setdiff(CC_list,CC_cols)) {
    meta[[variant]][column] = 0
  }
  # then filling in NAs (using AC_cols which is cols before adding new ones)
  meta[[variant]][CC_cols][is.na(meta[[variant]][CC_cols])] = 0
  
  # strata numbers in "no" columns 
  # again creating columns that aren't there with NAs 
  no_cols = names(meta[[variant]])[str_detect(names(meta[[variant]]),"_no")] 
  for (column in setdiff(no_list,no_cols)) { 
    meta[[variant]][column] = NA 
  } 
  
  # filling in "no" columns 
  meta[[variant]]["EUR_exomes_case_no"][is.na(meta[[variant]]["EUR_exomes_case_no"])] = 8874  
  meta[[variant]]["EUR_N_exomes_case_no"][is.na(meta[[variant]]["EUR_N_exomes_case_no"])] = 7277 
  meta[[variant]]["AMR_exomes_case_no"][is.na(meta[[variant]]["AMR_exomes_case_no"])] = 1388 
  meta[[variant]]["FIN_exomes_case_no"][is.na(meta[[variant]]["FIN_exomes_case_no"])] = 944 
  meta[[variant]]["FIN_genomes_case_no"][is.na(meta[[variant]]["FIN_genomes_case_no"])] = 423 
  meta[[variant]]["EAS_exomes_case_no"][is.na(meta[[variant]]["EAS_exomes_case_no"])] = 1730 
  meta[[variant]]["AFR_exomes_case_no"][is.na(meta[[variant]]["AFR_exomes_case_no"])] = 127 
  meta[[variant]]["AFR_genomes_case_no"][is.na(meta[[variant]]["AFR_genomes_case_no"])] = 2245 
  meta[[variant]]["ASJ_exomes_case_no"][is.na(meta[[variant]]["ASJ_exomes_case_no"])] = 869 
  meta[[variant]]["EST_genomes_case_no"][is.na(meta[[variant]]["EST_genomes_case_no"])] = 261 
  meta[[variant]]["SAS_exomes_case_no"][is.na(meta[[variant]]["SAS_exomes_case_no"])] = 110 
  meta[[variant]]["in_house_case_no"][is.na(meta[[variant]]["in_house_case_no"])] = n_case 
  meta[[variant]]["in_house_male_case_no"][is.na(meta[[variant]]["in_house_male_case_no"])] = n_male_case 
  meta[[variant]]["in_house_female_case_no"][is.na(meta[[variant]]["in_house_female_case_no"])] = n_female_case  
  
  meta[[variant]]["EUR_exomes_control_no"][is.na(meta[[variant]]["EUR_exomes_control_no"])] = 42635 
  meta[[variant]]["EUR_N_exomes_control_no"][is.na(meta[[variant]]["EUR_N_exomes_control_no"])] = 11187 
  meta[[variant]]["AMR_exomes_control_no"][is.na(meta[[variant]]["AMR_exomes_control_no"])] = 15154 
  meta[[variant]]["FIN_exomes_control_no"][is.na(meta[[variant]]["FIN_exomes_control_no"])] = 11526 
  meta[[variant]]["FIN_genomes_control_no"][is.na(meta[[variant]]["FIN_genomes_control_no"])] = 655 
  meta[[variant]]["EAS_exomes_control_no"][is.na(meta[[variant]]["EAS_exomes_control_no"])] = 8413 
  meta[[variant]]["AFR_exomes_control_no"][is.na(meta[[variant]]["AFR_exomes_control_no"])] = 765 
  meta[[variant]]["AFR_genomes_control_no"][is.na(meta[[variant]]["AFR_genomes_control_no"])] = 1590 
  meta[[variant]]["ASJ_exomes_control_no"][is.na(meta[[variant]]["ASJ_exomes_control_no"])] = 2963 
  meta[[variant]]["EST_genomes_control_no"][is.na(meta[[variant]]["EST_genomes_control_no"])] = 2281 
  meta[[variant]]["SAS_exomes_control_no"][is.na(meta[[variant]]["SAS_exomes_control_no"])] = 153 
  meta[[variant]]["in_house_control_no"][is.na(meta[[variant]]["in_house_control_no"])] = n_control  
  meta[[variant]]["in_house_male_control_no"][is.na(meta[[variant]]["in_house_male_control_no"])] = n_male_control 
  meta[[variant]]["in_house_female_control_no"][is.na(meta[[variant]]["in_house_female_control_no"])] = n_female_control 
  
  # reorganising columns 
  meta[[variant]] = meta[[variant]] %>%
    select(gene_id, gene_symbol, 
           EUR_exomes_case_CC, EUR_exomes_control_CC, EUR_exomes_case_no, EUR_exomes_control_no,
           EUR_N_exomes_case_CC, EUR_N_exomes_control_CC, EUR_N_exomes_case_no, EUR_N_exomes_control_no,
           AFR_genomes_case_CC, AFR_genomes_control_CC, AFR_genomes_case_no, AFR_genomes_control_no,
           EAS_exomes_case_CC, EAS_exomes_control_CC, EAS_exomes_case_no, EAS_exomes_control_no,
           AMR_exomes_case_CC, AMR_exomes_control_CC, AMR_exomes_case_no, AMR_exomes_control_no,
           FIN_exomes_case_CC, FIN_exomes_control_CC, FIN_exomes_case_no, FIN_exomes_control_no,
           ASJ_exomes_case_CC, ASJ_exomes_control_CC, ASJ_exomes_case_no, ASJ_exomes_control_no,
           FIN_genomes_case_CC, FIN_genomes_control_CC, FIN_genomes_case_no, FIN_genomes_control_no,
           EST_genomes_case_CC, EST_genomes_control_CC, EST_genomes_case_no, EST_genomes_control_no,
           AFR_exomes_case_CC, AFR_exomes_control_CC, AFR_exomes_case_no, AFR_exomes_control_no,
           SAS_exomes_case_CC, SAS_exomes_control_CC, SAS_exomes_case_no, SAS_exomes_control_no,
           in_house_case_CC, in_house_control_CC, in_house_case_no, in_house_control_no,
           in_house_male_case_CC, in_house_male_control_CC, in_house_male_case_no, in_house_male_control_no,
           in_house_female_case_CC, in_house_female_control_CC, in_house_female_case_no, in_house_female_control_no)
  
}

# performing CMH tests

# strata lists for "apply" 
schema_strata_list = c("EUR_exomes_case_CC", "EUR_exomes_control_CC", "EUR_exomes_case_no", "EUR_exomes_control_no",
                       "EUR_N_exomes_case_CC", "EUR_N_exomes_control_CC", "EUR_N_exomes_case_no", "EUR_N_exomes_control_no",
                       "AFR_genomes_case_CC", "AFR_genomes_control_CC", "AFR_genomes_case_no", "AFR_genomes_control_no",
                       "EAS_exomes_case_CC", "EAS_exomes_control_CC", "EAS_exomes_case_no", "EAS_exomes_control_no",
                       "AMR_exomes_case_CC", "AMR_exomes_control_CC", "AMR_exomes_case_no", "AMR_exomes_control_no", 
                       "FIN_exomes_case_CC", "FIN_exomes_control_CC", "FIN_exomes_case_no", "FIN_exomes_control_no", 
                       "ASJ_exomes_case_CC", "ASJ_exomes_control_CC", "ASJ_exomes_case_no", "ASJ_exomes_control_no",
                       "FIN_genomes_case_CC", "FIN_genomes_control_CC", "FIN_genomes_case_no", "FIN_genomes_control_no", 
                       "EST_genomes_case_CC", "EST_genomes_control_CC", "EST_genomes_case_no", "EST_genomes_control_no", 
                       "AFR_exomes_case_CC", "AFR_exomes_control_CC", "AFR_exomes_case_no", "AFR_exomes_control_no", 
                       "SAS_exomes_case_CC", "SAS_exomes_control_CC", "SAS_exomes_case_no", "SAS_exomes_control_no")
in_house_autosomal_strata_list = c("in_house_case_CC", "in_house_control_CC", "in_house_case_no", "in_house_control_no")
in_house_nonpar_strata_list = c("in_house_male_case_CC", "in_house_male_control_CC", "in_house_male_case_no", "in_house_male_control_no",
                                "in_house_female_case_CC", "in_house_female_control_CC", "in_house_female_case_no", "in_house_female_control_no")
autosomal_strata_list = c(schema_strata_list, in_house_autosomal_strata_list) 
nonpar_strata_list = c(schema_strata_list, in_house_nonpar_strata_list) 

# functions to check for number of non-0 strata 
pass_2sa_schema = function(x) { # 11 strata
  pairs = list("1" = c(x[[1]], x[[2]]),
               "2" = c(x[[5]], x[[6]]),
               "3" = c(x[[9]], x[[10]]),
               "4" = c(x[[13]], x[[14]]),
               "5" = c(x[[17]], x[[18]]),
               "6" = c(x[[21]], x[[22]]),
               "7" = c(x[[25]], x[[26]]),
               "8" = c(x[[29]], x[[30]]),
               "9" = c(x[[33]], x[[34]]),
               "10" = c(x[[37]], x[[38]]),
               "11" = c(x[[41]], x[[42]])) 
  if (!all(pairs[[1]] == 0)) { # if eur is nonzero 
    bool = TRUE } 
  else if (sum(sapply(pairs[2:11], function(x) !all(x == 0))) >=2) { # if >=2 other strata
    bool = TRUE } 
  else { bool = FALSE }
  return(bool) 
}

pass_2sa_autosomal = function(x) { # 12 strata
  pairs = list("1" = c(x[[1]], x[[2]]),
               "2" = c(x[[5]], x[[6]]),
               "3" = c(x[[9]], x[[10]]),
               "4" = c(x[[13]], x[[14]]),
               "5" = c(x[[17]], x[[18]]),
               "6" = c(x[[21]], x[[22]]),
               "7" = c(x[[25]], x[[26]]),
               "8" = c(x[[29]], x[[30]]),
               "9" = c(x[[33]], x[[34]]),
               "10" = c(x[[37]], x[[38]]),
               "11" = c(x[[41]], x[[42]]),
               "12" = c(x[[45]], x[[46]]))
  if (!all(pairs[[1]] == 0)) { # if eur is nonzero 
    bool = TRUE } 
  else if (sum(sapply(pairs[2:12], function(x) !all(x == 0))) >=2) { # if >=2 other strata
    bool = TRUE } 
  else { bool = FALSE }
  return(bool) 
}

pass_2sa_nonpar = function(x) { # 13 strata
  pairs = list("1" = c(x[[1]], x[[2]]),
               "2" = c(x[[5]], x[[6]]),
               "3" = c(x[[9]], x[[10]]),
               "4" = c(x[[13]], x[[14]]),
               "5" = c(x[[17]], x[[18]]),
               "6" = c(x[[21]], x[[22]]),
               "7" = c(x[[25]], x[[26]]),
               "8" = c(x[[29]], x[[30]]),
               "9" = c(x[[33]], x[[34]]),
               "10" = c(x[[37]], x[[38]]),
               "11" = c(x[[41]], x[[42]]),
               "12" = c(x[[45]], x[[46]]),
               "13" = c(x[[49]], x[[50]])) 
  if (!all(pairs[[1]] == 0)) { # if eur is nonzero 
    bool = TRUE } 
  else if (sum(sapply(pairs[2:13], function(x) !all(x == 0))) >=2) { # if >=2 other strata
    bool = TRUE } 
  else { bool = FALSE }
  return(bool) 
} 

in_house_nonpar_count_non_0_strata = function(x) { # 2 strata
  return(sum(!sapply(list("1" = c(x[[1]], x[[2]]),
                          "2" = c(x[[5]], x[[6]])), 
                     function(x) all(x == 0)))) } 
in_house_nonpar_first_non_0_strata = function(x) { # 2 strata
  return((x[[1]] + x[[2]]) > 0) } 

for (variant in c("ptv","ptv_mpc3","ptv_mpc2","mpc2")) { 
  
  cat(paste(if (variant == "ptv") {""} else {"\n"}, "Testing ", variant, " variants:\n", sep = "")) 
  
  meta[[variant]] = meta[[variant]] %>% 
    mutate_if(is.integer,as.numeric) 
  
  cat(paste("CMH tests for ", nrow(meta[[variant]]), " genes\n",sep="")) 
  
  # SCHEMA discovery - OR and p-val
  meta[[variant]]$schema_cmh_or = apply(meta[[variant]][,schema_strata_list], 1, 
                                             function(x) 
                                               if (pass_2sa_schema(x)) { 
                                                 mantelhaen.test(array(x,c(2,2,11)), alternative = "two.sided", correct = TRUE)[[5]][[1]] 
                                               } else {
                                                 NA 
                                               }, 
                                             simplify=TRUE)
  meta[[variant]]$schema_cmh_or_lower = apply(meta[[variant]][,schema_strata_list], 1, 
                                                   function(x) 
                                                     if (pass_2sa_schema(x)) { 
                                                       mantelhaen.test(array(x,c(2,2,11)), alternative = "two.sided", correct = TRUE)[[4]][[1]] 
                                                     } else {
                                                       NA 
                                                     }, 
                                                   simplify=TRUE)
  meta[[variant]]$schema_cmh_or_upper = apply(meta[[variant]][,schema_strata_list], 1, 
                                                   function(x) 
                                                     if (pass_2sa_schema(x)) { 
                                                       mantelhaen.test(array(x,c(2,2,11)), alternative = "two.sided", correct = TRUE)[[4]][[2]] 
                                                     } else {
                                                       NA 
                                                     }, 
                                                   simplify=TRUE)
  meta[[variant]]$schema_cmh_p_val = apply(meta[[variant]][,schema_strata_list], 1, 
                                                function(x) 
                                                  if (pass_2sa_schema(x)) { 
                                                    mantelhaen.test(array(x,c(2,2,11)), alternative = "two.sided", correct = TRUE)[[3]][[1]] 
                                                  } else {
                                                    NA 
                                                  }, 
                                                simplify=TRUE) 
  meta[[variant]]$schema_pass_2sa = apply(meta[[variant]][,schema_strata_list], 1, 
                                               function(x) 
                                                 if (pass_2sa_schema(x)) { 
                                                   TRUE 
                                                 } else { 
                                                   FALSE
                                                 })
  
  # in-house discovery - autosomal genes  
  meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), "in_house_cmh_or"] = 
    apply(meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), in_house_autosomal_strata_list], 1, 
          function(x) 
            fisher.test(matrix(x, nrow=2), alternative = "two.sided")[[3]]) 
  meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), "in_house_cmh_lower_or"] = 
    apply(meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), in_house_autosomal_strata_list], 1, 
          function(x) 
            fisher.test(matrix(x, nrow=2), alternative = "two.sided")[[2]][[1]]) 
  meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), "in_house_cmh_upper_or"] = 
    apply(meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), in_house_autosomal_strata_list], 1, 
          function(x) 
            fisher.test(matrix(x, nrow=2), alternative = "two.sided")[[2]][[2]]) 
  meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), "in_house_cmh_p_val"] = 
    apply(meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), in_house_autosomal_strata_list], 1, 
          function(x) 
            fisher.test(matrix(x, nrow=2), alternative = "two.sided")[[1]]) 
  
  # in-house discovery - nonpar genes  
  meta[[variant]][(meta[[variant]]$gene_id %in% nonpar_genes), "in_house_cmh_or"] = 
    apply(meta[[variant]][(meta[[variant]]$gene_id %in% nonpar_genes), in_house_nonpar_strata_list], 1, 
          function(x) 
            if (in_house_nonpar_count_non_0_strata(x) == 2) {
              mantelhaen.test(array(x,c(2,2,2)), alternative = "two.sided", correct = TRUE)[[5]][[1]] 
            } else if (in_house_nonpar_count_non_0_strata(x) == 1) {
              if (in_house_nonpar_first_non_0_strata(x)) {
                fisher.test(matrix(x[1:4], nrow=2), alternative = "two.sided")[[3]]
              } else {
                fisher.test(matrix(x[5:8], nrow=2), alternative = "two.sided")[[3]]
              }
            } else { 
              fisher.test(matrix(x[1:4], nrow=2), alternative = "two.sided")[[3]]
            }
    )
  meta[[variant]][(meta[[variant]]$gene_id %in% nonpar_genes), "in_house_cmh_lower_or"] = 
    apply(meta[[variant]][(meta[[variant]]$gene_id %in% nonpar_genes), in_house_nonpar_strata_list], 1, 
          function(x) 
            if (in_house_nonpar_count_non_0_strata(x) == 2) {
              mantelhaen.test(array(x,c(2,2,2)), alternative = "two.sided", correct = TRUE)[[4]][[1]] 
            } else if (in_house_nonpar_count_non_0_strata(x) == 1) {
              if (in_house_nonpar_first_non_0_strata(x)) {
                fisher.test(matrix(x[1:4], nrow=2), alternative = "two.sided")[[2]][[1]]
              } else {
                fisher.test(matrix(x[5:8], nrow=2), alternative = "two.sided")[[2]][[1]]
              }
            } else { 
              fisher.test(matrix(x[1:4], nrow=2), alternative = "two.sided")[[2]][[1]]
            }
    )
  meta[[variant]][(meta[[variant]]$gene_id %in% nonpar_genes), "in_house_cmh_upper_or"] = 
    apply(meta[[variant]][(meta[[variant]]$gene_id %in% nonpar_genes), in_house_nonpar_strata_list], 1, 
          function(x) 
            if (in_house_nonpar_count_non_0_strata(x) == 2) {
              mantelhaen.test(array(x,c(2,2,2)), alternative = "two.sided", correct = TRUE)[[4]][[2]] 
            } else if (in_house_nonpar_count_non_0_strata(x) == 1) {
              if (in_house_nonpar_first_non_0_strata(x)) {
                fisher.test(matrix(x[1:4], nrow=2), alternative = "two.sided")[[2]][[2]]
              } else {
                fisher.test(matrix(x[5:8], nrow=2), alternative = "two.sided")[[2]][[2]]
              }
            } else { 
              fisher.test(matrix(x[1:4], nrow=2), alternative = "two.sided")[[2]][[2]]
            }
    )
  meta[[variant]][(meta[[variant]]$gene_id %in% nonpar_genes), "in_house_cmh_p_val"] = 
    apply(meta[[variant]][(meta[[variant]]$gene_id %in% nonpar_genes), in_house_nonpar_strata_list], 1, 
          function(x) 
            if (in_house_nonpar_count_non_0_strata(x) == 2) {
              mantelhaen.test(array(x,c(2,2,2)), alternative = "two.sided", correct = TRUE)[[3]][[1]] 
            } else if (in_house_nonpar_count_non_0_strata(x) == 1) {
              if (in_house_nonpar_first_non_0_strata(x)) {
                fisher.test(matrix(x[1:4], nrow=2), alternative = "two.sided")[[1]]
              } else {
                fisher.test(matrix(x[5:8], nrow=2), alternative = "two.sided")[[1]]
              }
            } else { 
              fisher.test(matrix(x[1:4], nrow=2), alternative = "two.sided")[[1]]
            }
    )
  
  # meta discovery - autosomal genes 
  # meta OR and p-val and two-strata fail 
  meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), "cmh_or"] = 
    apply(meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), autosomal_strata_list], 1, 
          function(x) 
            if (pass_2sa_autosomal(x)) { 
              mantelhaen.test(array(x,c(2,2,12)), alternative = "two.sided", correct = TRUE)[[5]][[1]] 
            } else {
              NA 
            }, 
          simplify=TRUE) 
  meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), "cmh_or_lower"] = 
    apply(meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), autosomal_strata_list], 1, 
          function(x) 
            if (pass_2sa_autosomal(x)) { 
              mantelhaen.test(array(x,c(2,2,12)), alternative = "two.sided", correct = TRUE)[[4]][[1]] 
            } else {
              NA 
            }, 
          simplify=TRUE) 
  meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), "cmh_or_upper"] = 
    apply(meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), autosomal_strata_list], 1, 
          function(x) 
            if (pass_2sa_autosomal(x)) { 
              mantelhaen.test(array(x,c(2,2,12)), alternative = "two.sided", correct = TRUE)[[4]][[2]] 
            } else {
              NA 
            }, 
          simplify=TRUE) 
  meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), "cmh_p_val"] = 
    apply(meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), autosomal_strata_list], 1, 
          function(x) 
            if (pass_2sa_autosomal(x)) { 
              mantelhaen.test(array(x,c(2,2,12)), alternative = "two.sided", correct = TRUE)[[3]][[1]] 
            } else { 
              NA 
            }, 
          simplify=TRUE) 
  meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), "pass_2sa"] = 
    apply(meta[[variant]][!(meta[[variant]]$gene_id %in% nonpar_genes), autosomal_strata_list], 1, 
          function(x) 
            if (pass_2sa_autosomal(x)) {
              TRUE
            } else {
              FALSE
            })
  
  # meta discovery - nonpar genes 
  meta[[variant]][meta[[variant]]$gene_id %in% nonpar_genes, "cmh_or"] = 
    apply(meta[[variant]][meta[[variant]]$gene_id %in% nonpar_genes, nonpar_strata_list], 1, 
          function(x) 
            if (pass_2sa_nonpar(x)) { 
              mantelhaen.test(array(x,c(2,2,13)), alternative = "two.sided", correct = TRUE)[[5]][[1]] 
            } else {
              NA 
            }, 
          simplify=TRUE) 
  meta[[variant]][meta[[variant]]$gene_id %in% nonpar_genes, "cmh_or_lower"] = 
    apply(meta[[variant]][meta[[variant]]$gene_id %in% nonpar_genes, nonpar_strata_list], 1, 
          function(x) 
            if (pass_2sa_nonpar(x)) { 
              mantelhaen.test(array(x,c(2,2,13)), alternative = "two.sided", correct = TRUE)[[4]][[1]] 
            } else {
              NA 
            }, 
          simplify=TRUE) 
  meta[[variant]][meta[[variant]]$gene_id %in% nonpar_genes, "cmh_or_upper"] = 
    apply(meta[[variant]][meta[[variant]]$gene_id %in% nonpar_genes, nonpar_strata_list], 1, 
          function(x) 
            if (pass_2sa_nonpar(x)) { 
              mantelhaen.test(array(x,c(2,2,13)), alternative = "two.sided", correct = TRUE)[[4]][[2]] 
            } else {
              NA 
            }, 
          simplify=TRUE) 
  meta[[variant]][meta[[variant]]$gene_id %in% nonpar_genes, "cmh_p_val"] = 
    apply(meta[[variant]][meta[[variant]]$gene_id %in% nonpar_genes, nonpar_strata_list], 1, 
          function(x) 
            if (pass_2sa_nonpar(x)) { 
              mantelhaen.test(array(x,c(2,2,13)), alternative = "two.sided", correct = TRUE)[[3]][[1]] 
            } else {
              NA
            }, 
          simplify=TRUE) 
  meta[[variant]][meta[[variant]]$gene_id %in% nonpar_genes, "pass_2sa"] = 
    apply(meta[[variant]][meta[[variant]]$gene_id %in% nonpar_genes, nonpar_strata_list], 1, 
          function(x) 
            if (pass_2sa_nonpar(x)) {
              TRUE
            } else {
              FALSE
            })
  
  # generating counts and reordering
  meta[[variant]] = meta[[variant]] %>% 
    mutate(schema_case_CC = rowSums(across(c(EUR_exomes_case_CC, EUR_N_exomes_case_CC, AFR_genomes_case_CC, 
                                             EAS_exomes_case_CC, AMR_exomes_case_CC, FIN_exomes_case_CC, 
                                             ASJ_exomes_case_CC, FIN_genomes_case_CC, EST_genomes_case_CC, 
                                             AFR_exomes_case_CC, SAS_exomes_case_CC))),
           schema_control_CC = rowSums(across(c(EUR_exomes_control_CC, EUR_N_exomes_control_CC, AFR_genomes_control_CC,
                                                EAS_exomes_control_CC, AMR_exomes_control_CC, FIN_exomes_control_CC,
                                                ASJ_exomes_control_CC, FIN_genomes_control_CC, EST_genomes_control_CC,
                                                AFR_exomes_control_CC, SAS_exomes_control_CC))), 
           in_house_case_CC = rowSums(across(c(in_house_case_CC, in_house_male_case_CC, in_house_female_case_CC))),
           in_house_control_CC = rowSums(across(c(in_house_control_CC, in_house_male_control_CC, in_house_female_control_CC)))) %>% 
    select(gene_id, gene_symbol, 
           schema_case_CC, schema_control_CC, 
           schema_cmh_or, schema_cmh_or_lower, schema_cmh_or_upper, schema_cmh_p_val, schema_pass_2sa, 
           in_house_case_CC, in_house_control_CC, 
           in_house_cmh_or, in_house_cmh_lower_or, in_house_cmh_upper_or, in_house_cmh_p_val, 
           cmh_or, cmh_or_lower, cmh_or_upper, cmh_p_val, pass_2sa) 
  
  cat(paste("Number of genes passing two-strata adjustment: ", 
            nrow(meta[[variant]] %>% filter(pass_2sa == TRUE)), "\n", sep="")) 
  
}

# de novo meta-analysis

# importing and splitting de novo counts 
de_novo = list() 
de_novo_counts = read.table("/scratch/c.c21070635/2_meta-analysis/data_processed/de_novo_gene_table.tsv", 
                            header=TRUE)
de_novo[["ptv"]] = de_novo_counts %>%
  select(gene_id, gene_symbol, PTV_count, PTV_p_val) %>%
  rename(de_novo_count = PTV_count, de_novo_p = PTV_p_val)
de_novo[["ptv_mpc3"]] = de_novo_counts %>%
  select(gene_id, gene_symbol, PTV_MPC3_count, PTV_MPC3_p_val) %>%
  rename(de_novo_count = PTV_MPC3_count, de_novo_p = PTV_MPC3_p_val)
de_novo[["ptv_mpc2"]] = de_novo_counts %>%
  select(gene_id, gene_symbol, PTV_MPC2_count, PTV_MPC2_p_val) %>%
  rename(de_novo_count = PTV_MPC2_count, de_novo_p = PTV_MPC2_p_val)
de_novo[["mpc2"]] = de_novo_counts %>%
  select(gene_id, gene_symbol, MPC2_count, MPC2_p_val) %>%
  rename(de_novo_count = MPC2_count, de_novo_p = MPC2_p_val)

for (variant in c("ptv","ptv_mpc3","ptv_mpc2","mpc2")) {
  
  cat(paste(if (variant == "ptv") {""} else {"\n"}, "Meta-analysing de novo for ", variant, ":\n", sep = "")) 
  
  cat("Merging...\n")
  meta[[variant]] = merge(meta[[variant]], de_novo[[variant]], 
                          by=c("gene_id","gene_symbol"), all.x=TRUE)
  
  # replacing de_novo cols with NAs if case-control p-val > 0.01 
  meta[[variant]] = meta[[variant]] %>%
    mutate(de_novo_count = if_else(cmh_p_val <= 0.01, de_novo_count, NA),
           de_novo_p = if_else(cmh_p_val <= 0.01, de_novo_p, NA))
  
  # number of cmh p-vals below 0.01
  cat(paste0("P-vals below 0.01: ", nrow(meta[[variant]] %>% filter(cmh_p_val <= 0.01)), "\n"))
  
  # making meta_p_val col 
  meta[[variant]]$meta_p_val = apply(meta[[variant]][,c("cmh_p_val","de_novo_p")], 1, 
                                     function(x) 
                                       if (is.na(x[1]) | all(is.na(c(x)))) {
                                         NA
                                       } else 
                                         if (is.na(x[2])) {
                                           x[1]
                                         } else
                                           if (x[1] <= 0.01) {
                                             # if case-control p < 0.01  
                                             min(x[1],sumlog(x)[[3]])
                                           } else {
                                             x[1]
                                           })
  # for schema + de novo 
  meta[[variant]]$schema_de_novo_p_val = apply(meta[[variant]][,c("schema_cmh_p_val","de_novo_p")], 1, 
                                               function(x) 
                                                 if (is.na(x[1]) | all(is.na(c(x)))) {
                                                   NA
                                                 } else 
                                                   if (is.na(x[2])) {
                                                     x[1]
                                                   } else
                                                     if (x[1] <= 0.01) {
                                                       # if the case-control p < 0.01  
                                                       min(x[1],sumlog(x)[[3]])
                                                     } else {
                                                       x[1]
                                                     })
  
}

# writing  
for (variant in c("ptv","ptv_mpc3","ptv_mpc2","mpc2")) { 
  write_tsv(meta[[variant]], paste("/scratch/c.c21070635/2_meta-analysis/output/3_meta-analysis/2_meta_rvas/2_meta_gene_results_", variant, ".tsv", sep=""))  
}

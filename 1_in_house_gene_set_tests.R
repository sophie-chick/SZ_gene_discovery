# setup 
library(tidyverse) 
library(data.table) 
library(logistf) 

# loading gene sets
constrained_genes = read.table('/scratch/c.c21070635/2_meta-analysis/data/gnomad_pLI_Ensemble.tsv', header=TRUE)$Ensembl_ID 
schema = read.delim('/scratch/c.c21070635/2_meta-analysis/data/SCHEMA_gene_results_paper_with_ids.tsv')
sig_genes = c((schema %>% filter(P.meta <= 2.14e-06 & !is.na(P.meta)))$gene_id, 
              c("ENSG00000023516","ENSG00000167978"))
sig_fdr = (schema %>% filter(Q.meta <= 0.05 & !is.na(Q.meta)))$gene_id 
fdr_genes = sig_fdr[!(sig_fdr %in% sig_genes)] 
constrained_only = constrained_genes[!(constrained_genes %in% sig_fdr)]

gene_id_lists = list("constrained" = constrained_genes,
                     "constrained_only" = constrained_only,
                     "sig" = sig_genes,
                     "fdr" = fdr_genes) 

# sample table for covariates 
sample_table = read.delim("/scratch/c.c21070635/2_meta-analysis/data/SZ_AD_hail0.2_final_sample_table_full.tsv")

# importing in-house gene sample tables for each variant class 
# testing for case enrichment at minor allele counts (MAC) of 1 (singletons not in gnomAD)
# and 5 (MAC <=5 in in-house sample and MAC <=5 in gnomAD)
# using Firth's penalised logistic regressions 

input = list(); counts = list(); output = list() 

for (mac in c("1","5")) {
  
  cat(paste("\nMAC", mac, "\n", sep="")) 
  # constructing list of lists
  input[[mac]] = list(); counts[[mac]] = list(); output[[mac]] = list() 
  
  for (variant in c("ptv","mpc3","mpc2_3","syn")) {
    
    cat(paste("Importing ", variant, "\n", sep="")) 
    
    mac = "5"; variant = "syn"
    input[[mac]][[variant]] = read.table(paste0("/scratch/c.c21070635/2_meta-analysis/data_processed/in_house_gene_sample_table_mac", mac, "/", variant, ".tsv.gz"), 
                                         header=T, sep="\t", check.names = FALSE) 
    # transposing so that samples are rows 
    input[[mac]][[variant]] = transpose(input[[mac]][[variant]], keep.names = "Sample", make.names = "gene_id") 
    View(input[[mac]][[variant]])
  } 
  
  # making table with all_syn column for use as covariate 
  input[[mac]][["syn_counts"]] = input[[mac]][["syn"]] %>% 
    mutate(all_syn = rowSums(across(where(is.numeric)))) %>% 
    select(Sample, all_syn) 
  
  for (variant in c("ptv","mpc3","mpc2_3","syn")) {
    
    # counting variants across gene sets 
    for (gene_set in c("constrained","constrained_only","sig","fdr")) {
      
      counts[[mac]][[variant]][[gene_set]] = input[[mac]][[variant]] %>%
        select(Sample, one_of(gene_id_lists[[gene_set]])) %>% # selects all the constrained genes that are in the columns 
        mutate(count = rowSums(across(where(is.numeric)))) %>% # produces sum column 
        select(Sample, count) # selects sum column with total for each sample  
    
      # excluding any tests with no variants (mpc >3 variants in fdr genes) 
      if (sum(counts[[mac]][[variant]][[gene_set]][["count"]]) > 0) {
        
        # merging each counts table with sample_table containing PCs, predicted_sex 
        counts[[mac]][[variant]][[gene_set]] = merge(
          counts[[mac]][[variant]][[gene_set]], sample_table, by="Sample", all=TRUE) 
        # merging each table with table containing all_syn 
        counts[[mac]][[variant]][[gene_set]] = merge(
          counts[[mac]][[variant]][[gene_set]], input[[mac]][["syn_counts"]], 
          by="Sample", all=TRUE) 
        
        # regressing phenotype on variant counts, covarying for PCs, all_syn and predicted_sex
        # except for synonymous variant tests where all_syn is not used as a covariate 
        # saving test output encompassing beta, SE, OR, lower, upper, P 
        if (variant != "syn") {
          result = summary(logistf(data = counts[[mac]][[variant]][[gene_set]], 
                                   pheno ~ count + UK_PC1 + UK_PC2 + UK_PC3 +
                                     UK_PC4 + UK_PC5 + UK_PC6 + UK_PC7 + UK_PC8 + 
                                     UK_PC9 + UK_PC10 + all_syn + predicted_sex)) 
        } else {
          result = summary(logistf(data = counts[[mac]][[variant]][[gene_set]], 
                                   pheno ~ count + UK_PC1 + UK_PC2 + UK_PC3 +
                                     UK_PC4 + UK_PC5 + UK_PC6 + UK_PC7 + UK_PC8 + 
                                     UK_PC9 + UK_PC10 + predicted_sex)) 
        }
        
        # extracting counts, beta, SE, OR, CI lower bound, CI upper bound, P 
        output[[mac]][[variant]][[gene_set]] = 
          data.frame("case_vars" = by(counts[[mac]][[variant]][[gene_set]][["count"]], 
                                      counts[[mac]][[variant]][[gene_set]][["pheno"]], sum)[["1"]], 
                     "control_vars" = by(counts[[mac]][[variant]][[gene_set]][["count"]], 
                                         counts[[mac]][[variant]][[gene_set]][["pheno"]], sum)[["0"]],
                     "beta" = result[["coefficients"]][["count"]], 
                     "SE" = sqrt(diag(vcov(result))[["count"]]), 
                     "OR" = exp(result[["coefficients"]][["count"]]), 
                     "lower" = exp(result[["ci.lower"]][["count"]]), 
                     "upper" = exp(result[["ci.upper"]][["count"]]), 
                     "P" = result[["prob"]][["count"]]) 
        
      } else {
        
        # filling in rows of table where there are no variants of that class in a gene set
        output[[mac]][[variant]][[gene_set]] = 
          data.frame("case_vars" = 0, "control_vars" = 0, "beta" = NA, "SE" = NA, 
                     "OR" = NA, "lower" = NA, "upper" = NA, "P" = NA) 
      
      }
      
      output[[mac]][[variant]][[gene_set]]$Gene_set = gene_set 
      
    }
    
    # variant level - rbinding results dataframes 
    output[[mac]][[variant]] = do.call(rbind, output[[mac]][[variant]])
    output[[mac]][[variant]]$Variant = variant 
    
  }
  
  # MAC level - rbinding results dataframes 
  output[[mac]] = do.call(rbind, output[[mac]])
  output[[mac]]$MAC = mac 
  
}

# rbinding results dataframes
output = do.call(rbind, output)
row.names(output) = NULL 

# adding variant rate columns 
output = output %>%
  mutate(case_var_rate = case_vars/4650, 
         control_var_rate = control_vars/5719) %>% 
  select(Gene_set, Variant, MAC, case_vars, case_var_rate, control_vars, control_var_rate, 
         beta, SE, OR, lower, upper, P) 

# writing 
write.table(output, 
            "/scratch/c.c21070635/2_meta-analysis/output/1_in_house_analysis/1_gene_set_table.tsv", sep="\t")

# setup 
library(tidyverse) 
library(data.table) 
#detach("package:purrr", unload = TRUE) # stop purrr::transpose masking data.table::transpose
#install.packages('/scratch/c.c21070635/2_meta-analysis/files/packages/mice_3.15.0.tar.gz', repos=NULL)
#install.packages('/scratch/c.c21070635/2_meta-analysis/files/packages/operator.tools_1.6.3.tar.gz', repos=NULL)
#install.packages('/scratch/c.c21070635/2_meta-analysis/files/packages/formula.tools_1.7.1.tar.gz', repos=NULL)
#install.packages('/scratch/c.c21070635/2_meta-analysis/files/packages/logistf_1.25.0.tar.gz', repos=NULL)
library(logistf) 
#install.packages('common') 
#library(common) # for superscripts 

# uses gene sample tables imported in loop 
  
# load gene sets
constrained_genes = read.table('/scratch/c.c21070635/2_meta-analysis/data/gnomad_pLI_Ensemble.tsv',
                               header=TRUE)$Ensembl_ID # 3,051 genes 
schema = read.delim('/scratch/c.c21070635/2_meta-analysis/data/SCHEMA_gene_results_paper_with_ids.tsv')
sig_genes = (schema %>% filter(P.meta <= 2.14e-06 & !is.na(P.meta)))$gene_id 
sig_genes = c(sig_genes,c("ENSG00000023516","ENSG00000167978")) #akap11 and srrm2 
sig_fdr = (schema %>% filter(Q.meta <= 0.05 & !is.na(Q.meta)))$gene_id 
fdr_genes = sig_fdr[!(sig_fdr %in% sig_genes)] 
constrained_only = constrained_genes[!(constrained_genes %in% sig_fdr)]

gene_id_lists = list("constrained" = constrained_genes,
                     "constrained_only" = constrained_only,
                     "sig" = sig_genes,
                     "fdr" = fdr_genes) 

# sample table for covariates 
sample_table = read.table('/scratch/c.c21070635/2_meta-analysis/data/SZ_AD_hail0.2_final_sample_table.tsv', 
                          header = T, sep = "\t", stringsAsFactors = F) %>% 
  filter(Final_sample_exclusions == 0) 
sample_table$pheno = ifelse(sample_table$SZ_AD_con == "SZ", 1, 0) 
# need to remove the additional 8 samples so it lines up 
clozuk_list = c("ClozUKII_Welsh_30101856_Ca","ClozUKII_London_30033686_Ca","ClozUKII_Welsh_30148457_Ca",
                "ClozUKII_London_30066697_Ca","ClozUKII_London_30035834_Ca","ClozUKII_London_30034061_Ca") 
contributing_list = c("ClozUKII_Scotland_30034228_Ca","F0273_03_S1_Ca")
combined_list = c(clozuk_list,contributing_list)
sample_table = sample_table %>% 
  filter(!(SAMPLE %in% combined_list)) %>%
  select(SAMPLE, pheno, predicted_sex, 
         UK_PC1, UK_PC2, UK_PC3, UK_PC4, UK_PC5, UK_PC6, UK_PC7, UK_PC8, UK_PC9, UK_PC10) %>%
  rename(Sample = SAMPLE)

# writing for PTVs first 

input = list() 
counts = list() 
output = list() 

for (mac in c("1","5")) {
  
  cat(paste("\nMAC", mac, "\n", sep="")) 
  input[[mac]] = list() 
  counts[[mac]] = list() 
  output[[mac]] = list() 
  
  for (variant in c("ptv","mpc3","mpc2_3","syn")) {
    
    cat(paste("Importing ", variant, "\n", sep="")) 
    input[[mac]][[variant]] = list() 
    counts[[mac]][[variant]] = list() 
    output[[mac]][[variant]] = list() 
    
    input[[mac]][[variant]] = read.table(paste("/scratch/c.c21070635/2_meta-analysis/data_processed/in_house_gene_sample_table_mac", mac, "/", variant, ".tsv.gz", sep=""), 
                                         header=T, sep="\t", check.names = FALSE) %>%
      filter(!is.na(gene_id)) 
    input[[mac]][[variant]] = transpose(input[[mac]][[variant]], keep.names = "Sample", make.names = "gene_id") 
    
  } 
  
  # making table with all_syn column 
  input[[mac]][["syn_counts"]] = input[[mac]][["syn"]] %>% 
    mutate(all_syn = rowSums(across(where(is.numeric)))) %>% 
    select(Sample, all_syn) 
  
  for (variant in c("ptv","mpc3","mpc2_3","syn")) {
    
    # counting across gene sets 
    for (gene_set in c("constrained","constrained_only","sig","fdr")) {
      
      counts[[mac]][[variant]][[gene_set]] = input[[mac]][[variant]] %>%
        select(Sample, one_of(gene_id_lists[[gene_set]])) %>% # selects all the constrained genes that are in the columns 
        mutate(count = rowSums(across(where(is.numeric)))) %>% # adds sum column 
        select(Sample, count) # selects sum column - this is for each sample still 
    
      # avoiding fdr mpc >3
      if (sum(counts[[mac]][[variant]][[gene_set]][["count"]]) > 0) {
        
        # merging on PCs, all_syn and predicted_sex 
        counts[[mac]][[variant]][[gene_set]] = merge(
          counts[[mac]][[variant]][[gene_set]], sample_table, by="Sample", all=TRUE) 
        counts[[mac]][[variant]][[gene_set]] = merge(
          counts[[mac]][[variant]][[gene_set]], input[[mac]][["syn_counts"]], 
          by="Sample", all=TRUE) 
        
        # getting beta, SE, OR, lower, upper, P 
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
        
        output[[mac]][[variant]][[gene_set]] = 
          data.frame("case_vars" = 0, "control_vars" = 0, "beta" = NA, "SE" = NA, 
                     "OR" = NA, "lower" = NA, "upper" = NA, "P" = NA) 
      
      }
      
      output[[mac]][[variant]][[gene_set]]$Gene_set = gene_set 
      
    }
    
    # variant level
    output[[mac]][[variant]] = do.call(rbind, output[[mac]][[variant]])
    output[[mac]][[variant]]$Variant = variant 
    
  }
  
  # MAC level 
  output[[mac]] = do.call(rbind, output[[mac]])
  output[[mac]]$MAC = mac 
  
}

output = do.call(rbind, output)
row.names(output) = NULL 

output = output %>%
  mutate(case_var_rate = case_vars/4650, 
         control_var_rate = control_vars/5719) %>% 
  select(Gene_set, Variant, MAC, case_vars, case_var_rate, control_vars, control_var_rate, 
         beta, SE, OR, lower, upper, P) 

write.table(output, 
            "/scratch/c.c21070635/2_meta-analysis/output/1_in_house_analysis/1_gene_set_table.tsv", sep="\t")

output = read.delim("/scratch/c.c21070635/2_meta-analysis/output/1_in_house_analysis/1_gene_set_table.tsv")
output = output %>%
  filter(MAC == 5) %>% 
  mutate(Variant = recode(Variant, "ptv"="PTV", "mpc3"="Missense (MPC >3)",
                          "mpc2_3"="Missense (MPC 2-3)","syn"="Synonymous"),
         Gene_set = recode(Gene_set, "constrained"="Constrained genes\n(n=3,051)",
                           "constrained_only"="Constrained only",
                           "sig"="Significant genes\n(n=12)",
                           "fdr"="FDR<5% genes\n(n=20)")) %>%
  filter(Gene_set %in% c(#"Constrained genes\n(n=3,051)",
                         "Significant genes\n(n=12)",
                         "FDR<5% genes\n(n=20)"
                         )) 
output = output[order(factor(output$Gene_set,levels=c(#"Constrained genes\n(n=3,051)",
                                                      "Significant genes\n(n=12)",
                                                      "FDR<5% genes\n(n=20)")),
                      factor(output$Variant,levels=c("PTV", "Missense (MPC >3)",
                                                     "Missense (MPC 2-3)","Synonymous"))),] 



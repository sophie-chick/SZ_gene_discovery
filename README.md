## SZ gene discovery 

### 1: Gene set enrichment analysis in the new sample 

This script performs gene set enrichment tests in the new sample for four classes of variant (PTV, MPC >3, MPC 2-3 and synonymous) at two minor allele count thresholds (MAC = 1 and MAC <= 5) across three gene sets (constrained genes, genes previously implicated at exome-wide significance, and genes previously implicated at FDR <5%). Gene-level tables retaining individual sample information are read in for each variant class and sample-level counts are generated across each gene set. Case-control status is then regressed on variant count using Firth’s logistic regression tests, controlling for the first ten principal components, sex, and total burden of synonymous variants. 

### 2: Single gene tests in the meta-analysis sample 

We perform gene-level association tests for four variant classes (PTV, PTV + MPC >3, PTV + MPC >2, and MPC >2) across the combined case-control sample, using variants filtered at MAC <=5 in the combined sample and MAC <= 5 in gnomAD. Strata-level tables aggregating variants in these classes within each gene are read in and CMH tests are applied. For autosomal and pseudoautosomal genes, 12-strata tests are performed; for non-pseudoautosomal genes, 13-strata tests are performed to accommodate male and female strata in the new sample. CMH tests only account for strata with non-zero case and control variant counts; to avoid false positives, we apply a "two-strata correction" wherein genes with less than two non-zero strata, excepting those with variants in the largest European stratum, are not tested. Case-control CMH p-values below p < 0.01 are meta-analysed with de novo p-values for that variant class using Fisher's combined tests. 


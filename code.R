setwd("/home/amelie/nas/PKU_Inter99/merged_kasper/model/stratified/sex_interactions/sex_scaled/")
library(BEDMatrix)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(parallel)

bed <- BEDMatrix("/home/amelie/nas/PKU_Inter99/merged_kasper/data/clumped_stratified/PKU_Inter99_clumped_h6.bed", simple_names=TRUE)
#[1] "BEDMatrix" attr(,"package") [1] "BEDMatrix"
fam <- fread("/home/amelie/nas/PKU_Inter99/merged_kasper/data/clumped_stratified/PKU_Inter99_clumped_h6.fam") %>%
  dplyr::select(FID=V1, IID=V2, PID=V3, MID=V4, SEX=V5, PHENO=V6)
# "data.table" "data.frame"
# add phenotypes to fam file
phenotypes <- fread("/home/amelie/nas/PKU_Inter99/merged_kasper/phenotypes.txt") %>% dplyr::rename(FID=V1, IID=V2, PHENO=V3)
fam <- left_join(fam, phenotypes, by=c("IID"="IID")) %>% dplyr::select(FID=FID.x, IID, PID, MID, SEX, PHENO=PHENO.y) 

covar <- fread("/home/amelie/nas/PKU_Inter99/merged_kasper/GWAS/covar.txt", header=TRUE)
X <- left_join(fam, covar, by=c("FID"="FID")) %>% dplyr::select(FID, IID.x, PID, MID, SEX, PHENO, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10)
X$SEX <- as.numeric(X$SEX)
covar <- X %>% dplyr::select(C1, C2, C3, C4, C5, C6, C7, C8, C9, C10) %>% as.matrix() %>% scale()
# "matrix" "array" 
pheno <- X$PHENO # "integer"
pheno <- pheno-1 # "numeric"

sex <- X$SEX-1 # "numeric"
snp <- "rs9442400" # "character"
model <- stan_model(file = "/nas/users/amelie/PKU_Inter99/merged_kasper/model/model_anders.stan")  
system.time(
  r <- mclapply(colnames(bed), function(snp) {
    
    g_SNP <- bed[,snp] # "integer"
      
    X <- cbind(sex, g_SNP) #  "matrix" "array" 
    #na_sum <- sum(is.na(X))
    X <- cbind(X, 0)
    
    X <- cbind(X, covar) 
    X <- cbind(X, pheno) %>% na.omit()
    
    #X[,"g_17Q"] <- X[,"g_17Q"] - mean(X[,"g_17Q"])
    X[,"g_SNP"] <- X[,"g_SNP"] - mean(X[,"g_SNP"])
    X[,"sex"] <- X[,"sex"] - mean(X[,"sex"])
    X[,3] <- X[,"g_SNP"] * X[,"sex"]
    
    
    data <- list(N = nrow(X), # number of samples (phenotypes)
                 pheno = X[,"pheno"], # phenotype
                 M = 3, # number of SNPs
                 P = 10, # number of covariates
                 X = X[,1:13]) # predictor matrix
    
    
    ##--## sampling from stan model
    #system.time(
    fit <- sampling(
      model,  # Stan program
      data = data,    # named list of data
      chains = 4,             # number of Markov chains
      warmup = 100,          # number of warmup iterations per chain
      iter = 1000,            # total number of iterations per chain
      cores = 4,              # number of cores (could use one per chain)
      refresh = 0             # no progress shown
    )#)

    fit_summary <- summary(fit)
    
    beta_sex <- round(as.numeric(fit_summary$summary[2,"mean"]), 2)
    sd_sex <- round(as.numeric(fit_summary$summary[2,"sd"]), 2)
    n_eff_sex <- round(as.numeric(fit_summary$summary[2,"n_eff"]), 2)
    rhat_sex <- round(as.numeric(fit_summary$summary[2,"Rhat"]), 2)
    
    beta_SNP <- round(as.numeric(fit_summary$summary[3,"mean"]), 2)
    sd_SNP<- round(as.numeric(fit_summary$summary[3,"sd"]), 2)
    n_eff_SNP <- round(as.numeric(fit_summary$summary[3,"n_eff"]), 2)
    rhat_SNP <- round(as.numeric(fit_summary$summary[3,"Rhat"]), 2)
    
    beta_int <- round(as.numeric(fit_summary$summary[4,"mean"]), 2)
    sd_int <- round(as.numeric(fit_summary$summary[4,"sd"]), 2)
    n_eff_int <- round(as.numeric(fit_summary$summary[4,"n_eff"]), 2)
    rhat_int <- round(as.numeric(fit_summary$summary[4,"Rhat"]), 2)
    
    print(c(snp, beta_sex, sd_sex, n_eff_sex, rhat_sex, beta_SNP, sd_SNP, n_eff_SNP, rhat_SNP, 
            beta_int, sd_int, n_eff_int, rhat_int))
    
  }, mc.cores=5)
)

result_df <- as.data.frame(matrix(unlist(r), ncol=13, byrow=TRUE))
colnames(result_df) <- c("snp", "beta_sex", "sd_sex", "n_eff_sex", "rhat_sex", 
                         "beta_SNP", "sd_SNP", "n_eff_SNP", "rhat_SNP", 
                         "beta_int", "sd_int", "n_eff_int", "rhat_int")
write.table(result_df, file = "result_sex_scaled_bayes_h6.txt", quote = FALSE, row.names = FALSE)

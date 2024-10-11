## A SCRIPT TO COMPARE DIFFERENT MPS QUANTIFICATION TECHNIQUES
## ISABEL SCHUURMANS
## SCRIPT STARTED ON 01-09-2024

# Load libraries
library(haven)
library(tidyverse)
library(ggpubr)
library(tibble)
library(glmnet)
library(caret)

################################################################################
# LOAD DATA
################################################################################
# here, we load in the sum stats for calculating the MPS, the DNAm data,
# all the DNAm specific covariates, the phenotype of interest,
# and a developmental outcome

# ------------------------------------------------------------------------------
# summary stats EWAS

## load data
wdir <- "L:/people/Leo/ADHD/"
setwd(wdir)

EA_sumstats <- read.csv("data/sumstats/mateduewas_summarystatistics_Model1_CordBlood_EurAncestries.csv")
load("data/sumstats/41398_2020_1058_MOESM7_ESM.Rdata")
ADHD_sumstats <- adhd_birth_results.data
rm(adhd_birth_results.data)
SDP_sumstats <- readxl::read_xlsx("data/sumstats/TableS1_SDP.xlsx") %>% as.data.frame()

# filter EA and ADHD sumstats to include suggestive CpG sites only
ADHD_sumstats$FDR <- p.adjust(ADHD_sumstats$p_random_HE, method = "fdr")

ADHD_sumstats <- ADHD_sumstats %>% dplyr::select(-"MAPINFO") %>%
  filter(FDR < 0.05)
names(ADHD_sumstats) <- c("CpG", "CHR", "BETA", "SE", "P", "N", "FDR")

EA_sumstats$FDR <- p.adjust(EA_sumstats$P, method = "fdr")
EA_sumstats <- EA_sumstats %>% dplyr::select(-"Details") %>%
  filter(FDR < 0.05)
names(EA_sumstats)[2] <- "BETA" 

names(SDP_sumstats)[c(1:6)] <- c("CpG", "BETA", "SE", "P", "FDR", "CHR")
SDP_sumstats <- SDP_sumstats %>% mutate_at(c("BETA", "SE", "P", "FDR", "CHR"), as.numeric)

gc()

# -------------------------------------------------------------------------------
# DNAm data (winsorized)

## load data
library(stringr)
wdir_DNAm <- "L:/data/genetic_data/MoBa_methylation/" 
vars <- list.files(wdir_DNAm)


# ------------------------------------------------------------------------------
# prepared data

## load data
data <- miceadds::load.Rdata2(filename = "data/data_DNAm_pheno.RData")

data_mps <- data %>% filter(!is.na(Sample_Name))
gc()

# -------------------------------------------------------------------------------
# covariates and phenotypes

## check data
dim(data_mps) # 7964 individuals, 75 variables
str(data_mps)

data2 <- data_mps %>% arrange(Sample_Name)

gc()



######
# load prediction models trained in GenR
load("data/fit.enet_ADHD.RData")
load("data/fit.enet_EA.RData")
load("data/fit.enet_SDP.RData")


##################################################

########################################################################################################
# CALCULATE MPS ### OPTION 1: THRESHOLDING
########################################################################################################

# THIS CODE WAS WRITTEN BY NICOLE CREASEY

#betas = beta matrix (cgs as rows)
#sumstats = csv. with two columns labeled "CpG" (cg names) and "Weights" (weight for a given CpG)
#filename = the name you want to give an output file with any CpGs that were not available for the GESTBIR
#note. the function will not run if there are missing values in the beta matrix (but entire cgs can be missing from the matrix, it will compute without them and output a list)

#function for computing GESTBIR, returns unscaled GESTBIR ("mps") and scaled GESTBIR ("zmps")

methGESTBIR <- function(betas, sumstats, filename){
  
  #extract the betas
  epgs_betas <- as.data.frame(betas[rownames(betas) %in% sumstats$CpG,])
  
  #identify missings CpGs and remove from sumstats
  missing <- setdiff(sumstats$CpG, row.names(epgs_betas))
  sumstats <- sumstats %>% filter(!CpG %in% (missing))
  myfilename <- filename
  write.table(missing,file=myfilename, row.names = FALSE, col.names = FALSE)
  
  #transpose the betas
  t_epgs_betas <- data.frame(t(epgs_betas))
  
  #prepare coefs
  coefs <- sumstats$Weights
  names(coefs) <- sumstats$CpG
  coefs <- sumstats$Weights
  names(coefs) <- sumstats$CpG
  
  #compute
  GESTBIR <- (t(epgs_betas)) %*% coefs
  
  #create dataframe
  GESTBIR <- cbind(rownames(GESTBIR), data.frame(GESTBIR, row.names=NULL))
  colnames(GESTBIR)[1] <- "Sample_Name" #change this to your matching variable name
  colnames(GESTBIR)[2] <- "mps"
  
  #z GESTBIR
  GESTBIR$zmps <- scale(GESTBIR$mps,center = TRUE,scale = TRUE)
  colnames(GESTBIR)[3] <- "zmps"
  
  #return
  GESTBIR
  
}

# Sumstats, only two columns
# ADHD  
adhd_summstats_fdr <- ADHD_sumstats[ADHD_sumstats$FDR < 0.05, c('CpG','BETA')]
names(adhd_summstats_fdr) <- c('CpG', 'Weights')

# SDP
sdp_summstats_fdr <- SDP_sumstats[SDP_sumstats$FDR < 0.05, c('CpG','BETA')]
names(sdp_summstats_fdr) <- c('CpG', 'Weights')

# EA
ea_summstats_fdr <- EA_sumstats[EA_sumstats$FDR < 0.05, c('CpG','BETA')]
names(ea_summstats_fdr) <- c('CpG', 'Weights')

# indicator for batches
batches <- names(table(data$DATASET))

# data sets 1 and 2 contain multiple batches


for(i in c(1:2)){
  
  df <- data2 %>% filter(DATASET == batches[i])
  
  j = 0
  while (j < 2) {
    j = j + 1
    temp <- readRDS(paste0(wdir_DNAm, word(vars[i]), "/QC/results/", word(vars[i]), "-0", j, "_13_bmiqed_beta_values.rds"))
    colnames(temp) <- toupper(colnames(temp))
    assign(paste0("DNAm_", i,"_", j), temp) 
    x <- temp
    rm(temp)
    
    gc()
    
    ## omit missings, and check dimension (should not change too much)
    dim(x)
    y <- na.omit(x)
    dim(y) # same dimensions
    rm(x)
    
    # transpose the data in order for it to be rows fors participants, and columns for cpgs
    winsdata <- t(y)
    rm(y)
    gc()
    
    
    # this code from Isa/Nicky does not work for me (0 individuals left)
    winsdata_adhd <- as.data.frame(winsdata[rownames(winsdata) %in% data$Sample_Name,]) #only keep the individuals you need
    gc()
    winsdata_adhd$Sample_Name <- rownames(winsdata_adhd)
    gc()
    
    df2 <- df %>% filter(Sample_Name %in% winsdata_adhd$Sample_Name)
    df2 <- df2 %>% arrange(Sample_Name)
    
    
    winsdata_adhd2 <- winsdata_adhd %>% arrange(Sample_Name)
    gc()
    
    winsdata_adhd <- winsdata_adhd2
    rm(winsdata_adhd2)
    gc()
    
    ## check if data and dnam are the same order
    all.equal(df2$Sample_Name, winsdata_adhd$Sample_Name)  # should be true
    ## 1097 mismatches, so we need to put them in the same order
    
    ## still empty ofcourse
    mps <- data.frame(df2$Sample_Name)
    
    # transpose beta matrix (needed, cpgs as rows)
    winsdata_adhd$Sample_Name <- NULL
    gc()
    
    winsdata_adhd_t <- t(winsdata_adhd)
    
    ## MPS calculation
    ## ADHD ##
    
    # -------------------------------------------------------------------------------
    ## OPTION 1 A: we use ONE preselected p-value threshold
    
    ### run function
    mps1a_fdr <- methGESTBIR(betas=winsdata_adhd_t, sumstats = adhd_summstats_fdr, filename = paste0("mymps_missingcgs_ADHD_", i, "_", j))
    mps$mps_adhd_fdr <- mps1a_fdr$mps
    
    ### run function
    mps1a_fdr2 <- methGESTBIR(betas=winsdata_adhd_t, sumstats = sdp_summstats_fdr, filename = paste0("mymps_missingcgs_SDP_", i, "_", j))
    mps$mps_sdp_fdr <- mps1a_fdr2$mps
    
    ### run function
    mps1a_fdr2 <- methGESTBIR(betas=winsdata_adhd_t, sumstats = ea_summstats_fdr, filename = paste0("mymps_missingcgs_EA_", i, "_", j))
    mps$mps_ea_fdr <- mps1a_fdr2$mps
    
    names(mps) <- c("Sample_Name", "MPS_ADHD_FDR", "MPS_SDP_FDR", "MPS_EA_FDR")
    assign(paste0("mps_", i,"_", j), mps) # different fit object for each batch
    
    
    ########################################################################################################
    # CALCULATE MPS ### OPTION 2: PENALIZED REGRESSION
    ########################################################################################################
    
    # THIS CODE WAS WRITTEN BY NICOLE CREASEY
    
    # Remove CpGs that were not used for prediction model
    
    gc()
    
    # fit <- fit_ADHD
    # 
    # include <- winsdata_adhd[, colnames(winsdata_adhd) %in% adhd_summstats_fdr$CpG]
    # 
    # # Remove white space to match CpG names
    # colnames(include) <- gsub(" ", "", colnames(include))
    # colnames(fit$ptype) <- gsub(" ", "", colnames(fit$ptype))
    
    # targetData <- include[, colnames(include) %in% colnames(fit$ptype)]
    
    # fit$ptype <- fit$ptype[colnames(fit$ptype) %in% colnames(targetData)]
    # fit$finalModel$beta@Dimnames[[1]] <- fit$finalModel$beta@Dimnames[[1]][fit$finalModel$beta@Dimnames[[1]] %in% colnames(targetData)]
    # fit$finalModel$beta@Dim[1] <- length(fit$finalModel$beta@Dimnames[[1]])
    # fit$finalModel$xNames <-  fit$finalModel$xNames[fit$finalModel$xNames %in% colnames(targetData)]
    # 
    # mps$mps_adhd_glm <- NA
    # assign(paste0("fit_", i,"_", j), fit) # different fit object for each batch
    # 
    # mps$mps_adhd_glm <- predict(fit, newdata = targetData, type = "raw")
    
  }
  
  
}

rm(DNAm_1_1, DNAm_1_2, DNAm_2_1, DNAm_2_2)
gc()

# data sets 6 and 7 contain one batch (no index j needed)
# NOTE: batch m006 is on position 4 of the 'batches' variable

for(i in c(6,7)){
  
  df <- data2 %>% filter(DATASET == batches[i-2]) # batches 6 & 7 on positions 4 & 5
  
  temp <- readRDS(paste0(wdir_DNAm, word(vars[i]), "/QC/results/", word(vars[i]), "_13_bmiqed_beta_values.rds"))
  colnames(temp) <- toupper(colnames(temp))
  assign(paste0("DNAm_", i), temp) 
  x <- temp
  rm(temp)
  
  gc()
  
  ## omit missings, and check dimension (should not change too much)
  dim(x)
  y <- na.omit(x)
  dim(y) # same dimensions
  rm(x)
  
  # transpose the data in order for it to be rows fors participants, and columns for cpgs
  winsdata <- t(y)
  rm(y)
  gc()
  
  
  # this code from Isa/Nicky does not work for me (0 individuals left)
  winsdata_adhd <- as.data.frame(winsdata[rownames(winsdata) %in% data$Sample_Name,]) #only keep the individuals you need
  gc()
  winsdata_adhd$Sample_Name <- rownames(winsdata_adhd)
  gc()
  
  df2 <- df %>% filter(Sample_Name %in% winsdata_adhd$Sample_Name)
  df2 <- df2 %>% arrange(Sample_Name)
  
  
  winsdata_adhd2 <- winsdata_adhd %>% arrange(Sample_Name)
  gc()
  
  winsdata_adhd <- winsdata_adhd2
  rm(winsdata_adhd2)
  gc()
  
  ## check if data and dnam are the same order
  all.equal(df2$Sample_Name, winsdata_adhd$Sample_Name)  # should be true
  ## 1097 mismatches, so we need to put them in the same order
  
  ## still empty ofcourse
  mps <- data.frame(df2$Sample_Name)
  
  # transpose beta matrix (needed, cpgs as rows)
  winsdata_adhd$Sample_Name <- NULL
  gc()
  
  winsdata_adhd_t <- t(winsdata_adhd)
  
  ## MPS calculation
  ## ADHD ##
  
  ### run function
  mps1a_fdr <- methGESTBIR(betas=winsdata_adhd_t, sumstats = adhd_summstats_fdr, filename = paste0("mymps_missingcgs_ADHD_", i))
  mps$mps_adhd_fdr <- mps1a_fdr$mps
  
  ### run function
  mps1a_fdr2 <- methGESTBIR(betas=winsdata_adhd_t, sumstats = sdp_summstats_fdr, filename = paste0("mymps_missingcgs_SDP_", i))
  mps$mps_sdp_fdr <- mps1a_fdr2$mps
  
  ### run function
  mps1a_fdr2 <- methGESTBIR(betas=winsdata_adhd_t, sumstats = ea_summstats_fdr, filename = paste0("mymps_missingcgs_EA_", i))
  mps$mps_ea_fdr <- mps1a_fdr2$mps
  
  names(mps) <- c("Sample_Name", "MPS_ADHD_FDR", "MPS_SDP_FDR", "MPS_EA_FDR")
  assign(paste0("mps_", i), mps) # different fit object for each batch
}

rm(DNAm_6, DNAm_7)

gc()

# data sets 4, 8 and 9, 10 contain parents, so there is a different data structure
j = 0
while (j < 5) {
  for(i in c(4)){
    
    df <- data2 %>% filter(DATASET == batches[3])
    
    j = j + 1
    temp <- readRDS(paste0(wdir_DNAm, word(vars[i]), "/QC/child/results/", word(vars[i]), "_child-0", j, "_13_bmiqed_beta_values.rds"))
    colnames(temp) <- toupper(colnames(temp))
    assign(paste0("DNAm_", i,"_", j), temp) 
    x <- temp
    rm(temp)
    
    gc()
    
    ## omit missings, and check dimension (should not change too much)
    dim(x)
    y <- na.omit(x)
    dim(y) # same dimensions
    rm(x)
    
    # transpose the data in order for it to be rows fors participants, and columns for cpgs
    winsdata <- t(y)
    rm(y)
    gc()
    
    
    # this code from Isa/Nicky does not work for me (0 individuals left)
    winsdata_adhd <- as.data.frame(winsdata[rownames(winsdata) %in% data$Sample_Name,]) #only keep the individuals you need
    gc()
    winsdata_adhd$Sample_Name <- rownames(winsdata_adhd)
    gc()
    
    df2 <- df %>% filter(Sample_Name %in% winsdata_adhd$Sample_Name)
    df2 <- df2 %>% arrange(Sample_Name)
    
    
    winsdata_adhd2 <- winsdata_adhd %>% arrange(Sample_Name)
    gc()
    
    winsdata_adhd <- winsdata_adhd2
    rm(winsdata_adhd2)
    gc()
    
    ## check if data and dnam are the same order
    all.equal(df2$Sample_Name, winsdata_adhd$Sample_Name)  # should be true
    ## 1097 mismatches, so we need to put them in the same order
    
    ## still empty ofcourse
    mps <- data.frame(df2$Sample_Name)
    
    # transpose beta matrix (needed, cpgs as rows)
    winsdata_adhd$Sample_Name <- NULL
    gc()
    
    winsdata_adhd_t <- t(winsdata_adhd)
    
    ## MPS calculation

    ### run function
    mps1a_fdr <- methGESTBIR(betas=winsdata_adhd_t, sumstats = adhd_summstats_fdr, filename = paste0("mymps_missingcgs_ADHD_", i, "_", j))
    mps$mps_adhd_fdr <- mps1a_fdr$mps
    
    ### run function
    mps1a_fdr2 <- methGESTBIR(betas=winsdata_adhd_t, sumstats = sdp_summstats_fdr, filename = paste0("mymps_missingcgs_SDP_", i, "_", j))
    mps$mps_sdp_fdr <- mps1a_fdr2$mps
    
    ### run function
    mps1a_fdr2 <- methGESTBIR(betas=winsdata_adhd_t, sumstats = ea_summstats_fdr, filename = paste0("mymps_missingcgs_EA_", i, "_", j))
    mps$mps_ea_fdr <- mps1a_fdr2$mps
    
    names(mps) <- c("Sample_Name", "MPS_ADHD_FDR", "MPS_SDP_FDR", "MPS_EA_FDR")
    assign(paste0("mps_", i,"_", j), mps) # different fit object for each batch
  }
}

rm(DNAm_4_1, DNAm_4_2, DNAm_4_3, DNAm_4_4, DNAm_4_5)

gc()

j = 0
while (j < 4) {
  for(i in c(8)){
    
    df <- data2 %>% filter(DATASET == batches[i-2])
    
    j = j + 1
    temp <- readRDS(paste0(wdir_DNAm, word(vars[i]), "/QC/child/results/", word(vars[i]), "_child-0", j, "_13_bmiqed_beta_values.rds"))
    colnames(temp) <- toupper(colnames(temp))
    assign(paste0("DNAm_", i,"_", j), temp) 
    x <- temp
    rm(temp)
    
    gc()
    
    ## omit missings, and check dimension (should not change too much)
    dim(x)
    y <- na.omit(x)
    dim(y) # same dimensions
    rm(x)
    
    # transpose the data in order for it to be rows fors participants, and columns for cpgs
    winsdata <- t(y)
    rm(y)
    gc()
    
    
    # this code from Isa/Nicky does not work for me (0 individuals left)
    winsdata_adhd <- as.data.frame(winsdata[rownames(winsdata) %in% data$Sample_Name,]) #only keep the individuals you need
    gc()
    winsdata_adhd$Sample_Name <- rownames(winsdata_adhd)
    gc()
    
    df2 <- df %>% filter(Sample_Name %in% winsdata_adhd$Sample_Name)
    df2 <- df2 %>% arrange(Sample_Name)
    
    
    winsdata_adhd2 <- winsdata_adhd %>% arrange(Sample_Name)
    gc()
    
    winsdata_adhd <- winsdata_adhd2
    rm(winsdata_adhd2)
    gc()
    
    ## check if data and dnam are the same order
    all.equal(df2$Sample_Name, winsdata_adhd$Sample_Name)  # should be true
    ## 1097 mismatches, so we need to put them in the same order
    
    ## still empty ofcourse
    mps <- data.frame(df2$Sample_Name)
    
    # transpose beta matrix (needed, cpgs as rows)
    winsdata_adhd$Sample_Name <- NULL
    gc()
    
    winsdata_adhd_t <- t(winsdata_adhd)
    
    ## MPS calculation
    
    ### run function
    mps1a_fdr <- methGESTBIR(betas=winsdata_adhd_t, sumstats = adhd_summstats_fdr, filename = paste0("mymps_missingcgs_ADHD_", i, "_", j))
    mps$mps_adhd_fdr <- mps1a_fdr$mps
    
    ### run function
    mps1a_fdr2 <- methGESTBIR(betas=winsdata_adhd_t, sumstats = sdp_summstats_fdr, filename = paste0("mymps_missingcgs_SDP_", i, "_", j))
    mps$mps_sdp_fdr <- mps1a_fdr2$mps

    ### run function
    mps1a_fdr2 <- methGESTBIR(betas=winsdata_adhd_t, sumstats = ea_summstats_fdr, filename = paste0("mymps_missingcgs_EA_", i, "_", j))
    mps$mps_ea_fdr <- mps1a_fdr2$mps
    
    names(mps) <- c("Sample_Name", "MPS_ADHD_FDR", "MPS_SDP_FDR", "MPS_EA_FDR")
    assign(paste0("mps_", i,"_", j), mps) # different fit object for each batch
    
    gc()
  }
}


rm(DNAm_8_1, DNAm_8_2, DNAm_8_3, DNAm_8_4)

gc()

for(i in c(9:10)){
  
  df <- data2 %>% filter(DATASET == batches[i-2])
  
  j = 0
  while (j < 3) {
    j = j + 1
    temp <- readRDS(paste0(wdir_DNAm, word(vars[i]), "/QC/child/results/", word(vars[i]), "_child-0", j, "_13_bmiqed_beta_values.rds"))
    colnames(temp) <- toupper(colnames(temp))
    assign(paste0("DNAm_", i,"_", j), temp) 
    x <- temp
    rm(temp)
    
    gc()
    
    ## omit missings, and check dimension (should not change too much)
    dim(x)
    y <- na.omit(x)
    dim(y) # same dimensions
    rm(x)
    
    # transpose the data in order for it to be rows fors participants, and columns for cpgs
    winsdata <- t(y)
    rm(y)
    gc()
    
    
    # this code from Isa/Nicky does not work for me (0 individuals left)
    winsdata_adhd <- as.data.frame(winsdata[rownames(winsdata) %in% data$Sample_Name,]) #only keep the individuals you need
    gc()
    winsdata_adhd$Sample_Name <- rownames(winsdata_adhd)
    gc()
    
    df2 <- df %>% filter(Sample_Name %in% winsdata_adhd$Sample_Name)
    df2 <- df2 %>% arrange(Sample_Name)
    
    
    winsdata_adhd2 <- winsdata_adhd %>% arrange(Sample_Name)
    gc()
    
    winsdata_adhd <- winsdata_adhd2
    rm(winsdata_adhd2)
    gc()
    
    ## check if data and dnam are the same order
    all.equal(df2$Sample_Name, winsdata_adhd$Sample_Name)  # should be true
    ## 1097 mismatches, so we need to put them in the same order
    
    ## still empty ofcourse
    mps <- data.frame(df2$Sample_Name)
    
    # transpose beta matrix (needed, cpgs as rows)
    winsdata_adhd$Sample_Name <- NULL
    gc()
    
    winsdata_adhd_t <- t(winsdata_adhd)
    
    ## MPS calculation
    
    ### run function
    mps1a_fdr <- methGESTBIR(betas=winsdata_adhd_t, sumstats = adhd_summstats_fdr, filename = paste0("mymps_missingcgs_ADHD_", i, "_", j))
    mps$mps_adhd_fdr <- mps1a_fdr$mps
    
    ### run function
    mps1a_fdr2 <- methGESTBIR(betas=winsdata_adhd_t, sumstats = sdp_summstats_fdr, filename = paste0("mymps_missingcgs_SDP_", i, "_", j))
    mps$mps_sdp_fdr <- mps1a_fdr2$mps
    
    ### run function
    mps1a_fdr2 <- methGESTBIR(betas=winsdata_adhd_t, sumstats = ea_summstats_fdr, filename = paste0("mymps_missingcgs_EA_", i, "_", j))
    mps$mps_ea_fdr <- mps1a_fdr2$mps
    
    names(mps) <- c("Sample_Name", "MPS_ADHD_FDR", "MPS_SDP_FDR", "MPS_EA_FDR")
    assign(paste0("mps_", i,"_", j), mps) # different fit object for each batch
  }
}

rm(DNAm_9_1, DNAm_9_2, DNAm_9_3,DNAm_10_1, DNAm_10_2, DNAm_10_3)

gc()


# merge MPS of each sample
mps_merged <- rbind(mps_1_1, mps_1_2, mps_2_1, mps_2_2, mps_4_1, mps_4_2, mps_4_3, mps_4_4, mps_4_5,
                    mps_6, mps_7, mps_8_1, mps_8_2, mps_8_3, mps_8_4, mps_9_1, mps_9_2, mps_9_3,
                    mps_10_1, mps_10_2, mps_10_3)


data <- full_join(data, mps_merged)
data_mps <- full_join(data_mps, mps_merged)

save(data_mps, file = "data/MPS_NPR_and_pheno.RData")
save(data, file = "data/MPS_NPR_all.RData")




################################################################################
################################################################################


## but potentially still duplicates, so let's check
## check is the maternal ID column unique to ensure that there are no siblings
"length(unique(data$MOTHER)) == nrow(data) #returns TRUE if there are no duplicates

## we have duplicate maternal IDs, we should keep only one sibling
## gives you a data frame with a list of IDs and the number of times they occurred
n_occur <- data.frame(table(data$MOTHER))
n_occur[n_occur$Freq > 1,]

## tells you which IDs occurred more than once
dups <- data[data$MOTHER %in% n_occur$Var1[n_occur$Freq > 1],]

## select the duplicate that occurs first
dups2 <- dups[!duplicated(dups$MOTHER),]

## remove the duplicate that occurs first
data <- anti_join(data, dups2)
rm(dups, dups2, n_occur)"


# -------------------------------------------------------------------------------


# test with subsample, 500 individuals, 500 cpgs + IDs
'test <- winsdata[1:500, c(1:500, ncol(winsdata))]'


## put the IDs in phenotype_file and in methylation_file into the same order
# data2 = data[data$Sample_ID %in% intersect(data$Sample_ID, rownames(winsdata_adhd)),]
# data3 = data2[match(rownames(winsdata_adhd), data$Sample_ID),]

sum(df2$Sample_Name == rownames(winsdata_adhd))/ nrow(winsdata_adhd) #check that the IDs are in the same order, should return 1




## Test associations between MPS and pheno (adjusting for batch and cell type)
summary(lm(sum_ADHD ~ MPS_FDR_ADHD + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK, data = data)) # negative association, non-significant?

# glm model not really tested in GenR but used for MoBa (test sample)
summary(lm(sum_ADHD ~ MPS_GLM_ADHD + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK, data = data))


## Educational attainment (EA) ##

mps <- data.frame(data$Sample_ID)

ea_summstats_fdr <- EA_sumstats[EA_sumstats$FDR < 0.05, c('CpG','BETA')]
names(ea_summstats_fdr) <- c('CpG', 'Weights')

mps1a_fdr <- methGESTBIR(betas=winsdata_adhd_t, sumstats = ea_summstats_fdr, filename="mymps_missingcgs_ea")
mps$mps1a_fdr <- mps1a_fdr$mps


## penalised regression

X_ea <- as.data.frame(winsdata_adhd)
X_ea <- X_ea[,names(X_ea) %in% ea_summstats_fdr$CpG]
dim(X_ea)

## set seed
set.seed(138) # makes random processes reproducible

# get phenotype on which we want to train our data
data$EDUCM <- as.numeric(data$EDUCM)
data_EA <- data[!is.na(data$EDUCM), ]
Y_EA <-  data_EA$EDUCM

# PARTITION 100 TRAIN 0 TESTING (Testing in MoBa)

## partition the data
in.train <- createDataPartition(
  y =  data_EA$EDUCM,
  ## the outcome data are needed
  p = 1,
  ## The percentage of data in the
  ## training set
  list = FALSE
)
in.train_ea <- as.numeric(in.train)

## data.subset
training_ea_100 <- data_EA[ in.train_ea,]
testing_ea_100  <- data_EA[-in.train_ea,]
nrow(training_ea_100)
nrow(testing_ea_100)

## fit.glmnet
set.seed(20)
fit.enet_EA <- caret::train(y = training_ea_100$EDUCM,
                         x = X_ea[in.train_ea,],
                         method = "glmnet",
                         trControl = trainControl(method = "cv",
                                                  n = 10))
mps$mps2_100 <- NA
mps$mps2_100 <- predict(fit.enet_EA, newdata = X_ea, type = "raw")
names(mps) <- c("Sample_ID", "MPS_FDR_EA", "MPS_GLM_EA")

#data <- left_join(data, mps)


## Test associations between MPS and pheno (adjusting for batch and cell type)
summary(lm(MPS_FDR_EA ~ EDUCM + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK, data = data)) # positive association, significant
summary(lm(MPS_GLM_EA ~ EDUCM + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK, data = data)) # positive association, non-significant



## Smoking during pregnancy (SDP) ##

mps <- data.frame(data$Sample_ID)

sdp_summstats_fdr <- SDP_sumstats[SDP_sumstats$FDR < 0.05, c('CpG','BETA')]
names(sdp_summstats_fdr) <- c('CpG', 'Weights')

### run function
mps1a_fdr <- methGESTBIR(betas=winsdata_adhd_t, sumstats = sdp_summstats_fdr, filename="mymps_missingcgs_sdp")
mps$mps1a_fdr <- mps1a_fdr$mps


## select probes

## dna m matrix, this time with cpgs as columns
X_sdp <- as.data.frame(winsdata_adhd)
X_sdp <- X_sdp[,names(X_sdp) %in% sdp_summstats_fdr$CpG]
dim(X_sdp)

## set seed
set.seed(138) # makes random processes reproducible

# get phenotype on which we want to train our data
data$SMOKE_ALL <- as.numeric(data$SMOKE_ALL)
data_SDP <- data[!is.na(data$SMOKE_ALL), ]
Y_SDP <-  data_SDP$SMOKE_ALL

# PARTITION 100 TRAIN 0 TESTING (Testing in MoBa)

## partition the data
in.train <- createDataPartition(
  y =  data_SDP$SMOKE_ALL,
  ## the outcome data are needed
  p = 1,
  ## The percentage of data in the
  ## training set
  list = FALSE
)
in.train_sdp <- as.numeric(in.train)

## data.subset
training_sdp_100 <- data_SDP[ in.train_sdp,]
testing_sdp_100  <- data_SDP[-in.train_sdp,]
nrow(training_sdp_100)
nrow(testing_sdp_100)

## fit.glmnet
set.seed(20)
fit.enet_sdp <- caret::train(y = training_sdp_100$SMOKE_ALL,
                         x = X_sdp[in.train_sdp,],
                         method = "glmnet",
                         trControl = trainControl(method = "cv",
                                                  n = 10))
mps$mps2_100 <- NA
mps$mps2_100 <- predict(fit.enet_sdp, newdata = X_sdp, type = "raw")
names(mps) <- c("Sample_ID", "MPS_FDR_SDP", "MPS_GLM_SDP")

data <- left_join(data, mps)


## Test associations between MPS and pheno (adjusting for batch and cell type)
summary(lm(MPS_FDR_SDP ~ SMOKE_ALL + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK, data = data)) # positive association, significant

# glm model not really tested in GenR but used for MoBa (test sample)
summary(lm(MPS_GLM_SDP ~ SMOKE_ALL + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK, data = data))


# save GLM prediction models to apply in MoBa
save(fit.enet_ADHD, file = "output/fit.enet_ADHD.RData")
save(fit.enet_EA, file = "output/fit.enet_EA.RData")
save(fit.enet_sdp, file = "output/fit.enet_SDP.RData")


## PLOTS ##
# using FDR MPS

data$MPS_FDR_ADHD <- scale(data$MPS_FDR_ADHD)
data$MPS_FDR_ADHD_r <- psych::reverse.code(keys = c(-1), 
                                        item = data$MPS_FDR_ADHD, 
                                        mini = min(data$MPS_FDR_ADHD),
                                        maxi = max(data$MPS_FDR_ADHD))

# check if reversing MPS worked
summary(lm(sum_ADHD ~ MPS_FDR_ADHD_r + Sample_Plate + Bcell + CD4T + CD8T + Gran + Mono + NK, data = data)) 
# positive association, p-value identical


plot_ADHD <- ggplot(aes(MPS_FDR_ADHD_r, sum_ADHD), data = data) +
  geom_point() + geom_smooth(method = "lm", color = "purple") +
  theme_minimal() + labs(x = "MPS for ADHD", y = "ADHD Sumscore at Age 8y") +
  scale_x_continuous(limits = c(-3.5, 3.5))



data_EA <- data %>% filter(EDUCM != 1 & !is.na(EDUCM))
data_EA$MPS_FDR_EA <- scale(data_EA$MPS_FDR_EA)

plot_EA <- ggplot(aes(as.factor(EDUCM), MPS_FDR_EA), data = data_EA) +
  geom_dotplot(aes(color = as.factor(EDUCM), fill = as.factor(EDUCM), alpha = 0),
               binwidth = 0.15, stackdir = "center",
               binaxis = 'y',
               stackratio = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 1.5, show.legend = F) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.1) +
  theme_minimal() + labs(y = "MPS for Educational Attainment", x = "Maternal Education") +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(-4, 4))



data_SDP <- data %>% filter(!is.na(SMOKE_ALL))
data_SDP$MPS_FDR_SDP <- scale(data_SDP$MPS_FDR_SDP)

plot_SDP <- ggplot(aes(as.factor(SMOKE_ALL), MPS_FDR_SDP), data = data_SDP) +
  geom_dotplot(aes(color = as.factor(SMOKE_ALL), fill = as.factor(SMOKE_ALL), alpha = 0),
               binwidth = 0.15, stackdir = "center",
               binaxis = 'y',
               stackratio = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 1.5, show.legend = F) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.1) +
  theme_minimal() + labs(y = "MPS for Smoking During Pregnancy", x = "Maternal Smoking During Pregnancy") +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(-4, 4))

# save plots
png(filename = "output/MPS_ADHD.png", width = 4000, height = 3000, type = "cairo", res = 500)
plot_ADHD
dev.off()

pdf("output/MPS_ADHD.pdf", width = 11)
plot_ADHD
dev.off()

png(filename = "output/MPS_EA.png", width = 4000, height = 3000, type = "cairo", res = 500)
plot_EA
dev.off()

pdf("output/MPS_EA.pdf", width = 11)
plot_EA
dev.off()

png(filename = "output/MPS_SDP.png", width = 4000, height = 3000, type = "cairo", res = 500)
plot_SDP
dev.off()

pdf("output/MPS_SDPD.pdf", width = 11)
plot_SDP
dev.off()

##########################################

# BELOW NOT RELEVANT AS WE'RE TESTING IN MOBA

##########################################

# -------------------------------------------------------------------------------
# PARTITION 25 TRAIN 75 TESTING

## partition the data
in.train <- createDataPartition(
    y =  data$GESTBIR,
    ## the outcome data are needed
    p = .25,
    ## The percentage of data in the
    ## training set
    list = FALSE
)
in.train_adhd <- as.numeric(in.train)

## data.subset
training_adhd_25te <- data[ in.train_adhd,]
testing_adhd_25te  <- data[-in.train_adhd,]
nrow(training_adhd_25te)
nrow(testing_adhd_25te)

## fit.glmnet
set.seed(20)
fit.enet <- caret::train(y = training_adhd_25te$GESTBIR,
                         x = X_adhd[in.train_adhd,],
                         method = "glmnet",
                         trControl = trainControl(method = "cv",
                                                  n = 10))
mps$mps2_25te <- NA
mps$mps2_25te[-in.train_adhd] <- predict(fit.enet, newdata = X_adhd[-in.train_adhd,], type = "raw")

# -------------------------------------------------------------------------------
# PARTITION 50 TRAIN 50 TESTING

## partition the data
in.train <- createDataPartition(
    y =  data$GESTBIR,
    ## the outcome data are needed
    p = .50,
    ## The percentage of data in the
    ## training set
    list = FALSE
)
in.train_adhd <- as.numeric(in.train)

## data.subset
training_adhd_50te <- data[ in.train_adhd,]
testing_adhd_50te  <- data[-in.train_adhd,]
nrow(training_adhd_50te)
nrow(testing_adhd_50te)

## fit.glmnet
set.seed(20)
fit.enet <- caret::train(y = training_adhd_50te$GESTBIR,
                         x = X_adhd[in.train_adhd,],
                         method = "glmnet",
                         trControl = trainControl(method = "cv",
                                                  n = 10))
mps$mps2_50te <- NA
mps$mps2_50te[-in.train_adhd] <- predict(fit.enet, newdata = X_adhd[-in.train_adhd,], type = "raw")

# -------------------------------------------------------------------------------
# PARTITION 75 TRAIN 25 TESTING

## partition the data
in.train <- createDataPartition(
    y =  data$GESTBIR,
    ## the outcome data are needed
    p = .75,
    ## The percentage of data in the
    ## training set
    list = FALSE
)
in.train_adhd <- as.numeric(in.train)

## data.subset
training_adhd_75te <- df_covs[ in.train_adhd,]
testing_adhd_75te  <- df_covs[-in.train_adhd,]
nrow(training_adhd_75te)
nrow(testing_adhd_75te)

## fit.glmnet
set.seed(20)
fit.enet <- caret::train(y = training_adhd_75te$GESTBIR,
                         x = X_adhd[in.train_adhd,],
                         method = "glmnet",
                         trControl = trainControl(method = "cv",
                                                  n = 10))
mps$mps2_75te <- NA
mps$mps2_75te[-in.train_adhd] <- predict(fit.enet, newdata = X_adhd[-in.train_adhd,], type = "raw")


########################################################################################################
# CALCULATE MPS ### OPTION 3: P+T CoMeBack
########################################################################################################

# scripts from Chen, 2023, epigenetics
# -------------------------------------------------------------------------------
# STEP ONE

# Load phenotypes file
#Need two columns, ID and ethnicity/race
pheno = as.data.frame(df_covs[,c("SampleID")])
pheno$ethnicity <- 1 # cause all are Ditcj
colnames(pheno) = c("ID", "race")
head(pheno)

# Load DNA methylation (beta) file
# rows need to be participants, cols are cpgs
str(winsdata_adhd)
#betas = slot(betas,"assayData")[["betas"]]
#betas = data.frame(t(betas))
winsdata_adhd$ID = rownames(winsdata_adhd)  

# Correct for ancestry first (regress out ancestry)
# Could skip this step and just use betas to calculate co-methylation regions
## is: not needed as only Dutch
#betas_final = merge(betas, pheno, by = "ID")
#res <- matrix(NA, ncol = ncol(betas) - 1, nrow = nrow(betas_final))
#rownames(res) <- betas_final$ID
#colnames(res) <- colnames(betas_final)[2:ncol(betas)]

#for(i in 1:(ncol(betas)-1)){
#  cpg = colnames(betas_final)[(i+1)]
#  dat2 = betas_final[,c(cpg, "race")]
#  colnames(dat2) = c("cpg", "race")
#  lm1 <- lm(cpg ~ race, data = dat2)
#  res[,i] = resid(lm1)
#}
#save(res, file="betas_test.RData")

# -------------------------------------------------------------------------------
# Run cmr function to generate Co-methylation regions

### Generate co-methylated regions ###
#install.packages("/home/i.schuurmans/MPS/comeback_1.0.1.tar.gz", repos = NULL, type="source")
#comeback_0.1.0.tar.gz can be found at https://bitbucket.org/flopflip/comeback.
library(comeback)

GenCoMeRegion = function(cmrs = NULL, beta, reference = NULL, Overlap = F){
    if (is.null(cmrs) & is.null(reference)) {
        print("No cmrs file or reference cmrs file. Please input at least one of cmrs file or reference cmrs file")
        exit()
    }
    if (is.null(beta)) {
        print("No DNA methylation (beta) file. Please input DNA methylation (beta) file")
        exit()
    }
    if (Overlap == T){
        if (is.null(cmrs) & !is.null(reference)) {
            print("No cmrs file. Please input cmrs file")
            exit()
        }else if (is.null(reference) & !is.null(cmrs)){
            print("No reference file. Please input reference cmrs file")
            exit()
        }else{
            chmr_refs_flt = unlist(reference,recursive = F) # this is the reference list of 24K CMRs
            cmrs_input_flt = unlist(cmrs,recursive = F) # the evan CMRs
            chmr_refs_prb = unlist(reference) # probes for reference list of 24K CMRs
            cmrs_input_prb = unlist(cmrs) # probes for evan CMRs
            which_ref_ovrlp = which(sapply(chmr_refs_flt, function(x){length(intersect(x,cmrs_input_prb))>0}))
            which_input_NOTovrlp = which(sapply(cmrs_input_flt, function(x){length(intersect(x,chmr_refs_prb))==0}))
            combined_cmr = list()
            combined_cmr = c(combined_cmr, chmr_refs_flt[which_ref_ovrlp], cmrs_input_flt[which_input_NOTovrlp])
            cmrsnew = combined_cmr
        }
    }else if (Overlap == F & !is.null(cmrs)){
        print("Use cmrs file to generate CoMeRegion")
        cmrsnew = unlist(cmrs,recursive = F)
    }else if (Overlap == F & !is.null(reference)){
        print("Use reference file to generate CoMeRegion")
        cmrsnew = unlist(reference,recursive = F)
    }
    
    CoMeRegion = NULL
    for (i in 1: length(cmrsnew)){
        #print(paste0("run this ", i, " times"))
        clustercpg = unlist(cmrsnew[i])
        data_cpg = data.frame(clustercpg, CoMethylCluster = i)
        CoMeRegion = rbind(CoMeRegion,data_cpg)
    }
    print("finished running 1")
    
    CpGnames = colnames(beta)
    print(length(CpGnames))
    
    CpGnames = CpGnames[!(CpGnames %in% CoMeRegion$clustercpg)]
    
    print(length(CpGnames))
    print("finished running 2")
    
    data_cpg = data.frame(clustercpg = CpGnames,
                          CoMethylCluster = (max(CoMeRegion$CoMethylCluster) +1):(length(CpGnames)  + max(CoMeRegion$CoMethylCluster)))
    CoMeRegion = rbind(CoMeRegion,data_cpg)
    print(dim(CoMeRegion))
    
    return(CoMeRegion)
}

GenMRS = function(beta, SS, Pthred, CoMeRegion, CoMeBack = T, weightSE = F){
    pvalueinfo = matrix(NA, length(Pthred),5)
    colnames(pvalueinfo) = c("Pvalue", "Number of CpG sites", "Numeber of CpG sites after matching with CoMeRegion", " Number of CpG sites after pruning", "Numeber of CpG sites after matching with DNA methylation data")
    if (weightSE == T){
        SS$BETA = SS$BETA/SS$SE
    }
    for (i in Pthred){
        pvalue = 5 * 10 ^ (-i)
        pvalueinfo[i-1,1] =  pvalue
        SS_sub = SS[SS$Pvalue < pvalue, ]
        pvalueinfo[i-1,2] =  nrow(SS_sub)
        
        if (CoMeBack){
            SS_final = merge(SS_sub, CoMeRegion, by.x = "Marker", by.y = "clustercpg")
            SS_final = data.frame(SS_final)
            pvalueinfo[i-1,3] =  nrow(SS_final)
            SS_final = SS_final %>%
                group_by(CoMethylCluster) %>%
                slice(which.min(Pvalue))
            pvalueinfo[i-1,4] =  nrow(SS_final)
        }else{
            SS_final = merge(SS_sub, CoMeRegion, by.x = "Marker", by.y = "clustercpg")
            SS_final = data.frame(SS_final)
        }
        betas_final = data.frame(t(beta[,SS_final$Marker, drop = F]))
        betas_final$Marker = rownames(betas_final)
        mat_final = merge(betas_final, SS_final, by = "Marker")
        pvalueinfo[i-1,5] =  nrow(mat_final)
        if (nrow(mat_final) > 0){
            MRS = data.frame(t(as.matrix(mat_final[,2:(ncol(betas_final))])) %*% (mat_final$BETA))
            colnames(MRS) = c(paste0("P",pvalue))
            MRS$ID = rownames(MRS)
            if (i == Pthred[1]){
                MRS_final = MRS
            }else{
                MRS_final = merge(MRS, MRS_final, by = "ID", all = T)
            }
        }
    }
    results = list(pvalueinfo = pvalueinfo, MRS = MRS_final)
    return(results)
}

#column are CpG sites, rows are for samples
#cmrs <- cmr(Mdata = winsdata_adhd, corlo = 0.3)
#cmrs is the output from cmr(), and saved all the co-methylation regions in a list
#save(cmrs, file="CMR.test.RData")

#read ref cmr
#cmrs_ref <- readRDS("ref.rds")

# To generate Co-methylation regions dataframe for later analysis,
# you can either input cmrs calculated from your own dataset,
#CoMeRegion = GenCoMeRegion(cmrs = cmrs, beta = winsdata_adhd, Overlap = F) # we are using this option
# or a cmrs calculated from a reference dataset,
#CoMeRegion = GenCoMeRegion(beta = res, reference = cmrs_ref, Overlap = F)
# or both and then overlap the cmrs
#CoMeRegion = GenCoMeRegion(cmrs = cmrs, beta = res, reference = cmrs_ref, Overlap = T)

#CoMeRegion is a matrix that assigned a unique number to each co-methylation region
#save(CoMeRegion, file = "CoMeRegion.rda")
load('CoMeRegion.rda')

# -------------------------------------------------------------------------------
# STEP TWO
# Calculate MRS

library(robustHD)

# Load real phenotype and methylation data
#DNAm = get(load("betas_test.RData"))
#DNAm = t(DNAm)
#DNAm[1:6,1:6] #check DNAm file

###Load smoking summary statistics from discovery dataset
#SS_newborn <- read.csv("SmokingSS_Adults_CHARGE.csv")
#The summary statistics file should have at least four columns: CpGs names, beta coefficient,
#standand errors and p-values
SS = adhd_summstats[,c('MarkerName', 'Effect', 'StdErr', 'P.value')]
#change the columns according to your dataset, 1st: CpGs, 2nd: coefficients, 3rd: SE and 4th: p-value
#Rename the column names as "Marker", "BETA" and "Pvalue"
colnames(SS) = c("Marker", "BETA", "SE", "Pvalue")
head(SS)
#      Marker    BETA    SE Pvalue
#1 cg05575921 -0.1 1  1.00e-2
#2 cg12803068  0.1 1  1.00e-2
#3 cg04180046  0.1 1  1.00e-2
#4 cg09935388 -0.1 1  1.00e-2
#5 cg25949550 -0.1 1  1.00e-2
#6 cg12876356 -0.1 1  1.00e-2

#Get the smallest p-value
minpvalue = min(SS$Pvalue[SS$Pvalue != 0])
minpvalue = sapply(strsplit(as.character(minpvalue), "-"), "[", 2)
###Load Co-methylation regions for newborns -> CoMeRegion
load("CoMeRegion.rda")
#Specify how p-value threshold, for example, if you want 5 * 10 ^ (-2), specify pthread to be 2
Pthred = 2:minpvalue
MRS = GenMRS(winsdata_adhd, SS, Pthred, CoMeRegion, CoMeBack = T, weightSE = F)
#if weightSE = T, weights = BETA/SE, where BETA is the effect size
#Basic information of MRS
write.csv(MRS$pvalueinfo, "MRS_pvalueinfo.csv", row.names = F)
write.csv(MRS$MRS, "MRS.csv", row.names = F)

mps$mps3 <- MRS$MRS
save(mps, file="mps_birthweight_20240111.RData")

########################################################################################################
# TEST THE DIFFERENT GESTBIRS
########################################################################################################

# load the MPSs
load("MPS/mps_birthweight_20240111.RData")

# fix mps and mps3 binding
mps3 <- as.data.frame(mps$mps3)
mps$mps3 <- NULL
names(mps3) <- paste0('mps3_', names(mps3))
mps3$mps3ID <- NULL
mps1 <- cbind(mps, mps3)
# correlations across MPSs
cor(mps1[,-1], use = 'pairwise.complete.obs')

## put the IDs in phenotype_file and in mps into the same order
df_covs = df_covs[df_covs$SampleID %in% intersect(df_covs$SampleID, mps$df_covs.SampleID),]
df_covs = df_covs[match(mps$df_covs.SampleID, df_covs$SampleID),]
sum(df_covs$SampleID == mps$df_covs.SampleID)/nrow(mps) #check that the IDs are in the same order, should return 1
all.equal(as.character(df_covs$SampleID),mps$df_covs.SampleID) # it is true now


# MPS1A

## the N of the test set
summary(mps$mps1a)
table(!is.na(mps$mps1a))
hist(mps$mps1a)

## assocation with phenotype
cor.test(mps$mps1a_00005, df_covs$GESTBIR)
model1 <- lm(scale(GESTBIR) ~ scale(mps$mps1a_00005) +
                 GENDER +
                 CD8T + CD4T + NK + Bcell + Mono + Gran + nRBC +
                 C1 + C2 + C3 + C4 + C5 +
                 Sample_Plate,
             data = df_covs)
summary(model1)

## association with outcome
cor.test(mps$mps1a_00005, df_covs$cbcl_total_problems)
model1 <- lm(scale(cbcl_total_problems) ~ scale(mps$mps1a_00005) + cbcl_age +
                 GENDER +
                 CD8T + CD4T + NK + Bcell + Mono + Gran + nRBC +
                 C1 + C2 + C3 + C4 + C5 +
                 Sample_Plate,
             data = df_covs)
summary(model1)

# MPS1B

# get training and testing set
set.seed(20)

## partition the data
in.train <- createDataPartition(
    y =  df_covs$GESTBIR,
    ## the outcome data are needed
    p = .50,
    ## The percentage of data in the
    ## training set
    list = FALSE
)
in.train_adhd <- as.numeric(in.train)

## data.subset
training_adhd_50te <- df_covs[ in.train_adhd,]
testing_adhd_50te  <- df_covs[-in.train_adhd,]
training_mps_50te <- mps[ in.train_adhd,]
testing_mps_50te  <- mps[-in.train_adhd,]

## first, we need to select which p-value threshold work best
pvalue_1b <- matrix(NA, 13, 3)

# writing a little loop to get correlations across all mpss
mps1b <- names(mps)[grep('1b',names(mps))]

for (i in 1:length(mps1b)){
    pvalue_1b[i,1] <- mps1b[i]
    pvalue_1b[i,2] <- summary(lm(scale((training_mps_50te[,mps1b[i]])) ~ scale(training_adhd_50te$GESTBIR)))$coef[2,1]
    pvalue_1b[i,3] <- summary(lm(scale((training_mps_50te[,mps1b[i]])) ~ scale(training_adhd_50te$GESTBIR)))$coef[2,4]}

pvalue_1b[(pvalue_1b[,3] == min(pvalue_1b[,3])),]
# mps1b_00000005 is the best one

summary(mps$mps1b_00000005)
table(!is.na(mps$mps1b_00000005))
hist(mps$mps1b_00000005)

## assocation with phenotype
cor.test(testing_mps_50te[, 'mps1b_00000005'], testing_adhd_50te$GESTBIR)
model1 <- lm(scale(GESTBIR) ~ scale(testing_mps_50te$mps1b_00000005) +
                 GENDER +
                 CD8T + CD4T + NK + Bcell + Mono + Gran + nRBC +
                 C1 + C2 + C3 + C4 + C5 +
                 Sample_Plate,
             data = testing_adhd_50te)
summary(model1)

## association with outcome
cor.test(testing_mps_50te$mps1b_00000005, testing_adhd_50te$cbcl_total_problems)
model1 <- lm(scale(cbcl_total_problems) ~ scale(testing_mps_50te$mps1b_00000005) + cbcl_age +
                 GENDER +
                 CD8T + CD4T + NK + Bcell + Mono + Gran + nRBC +
                 C1 + C2 + C3 + C4 + C5 +
                 Sample_Plate,
             data = testing_adhd_50te)
summary(model1)


# MPS2

##25/75

## the N of the test set
summary(mps$mps2_25te)
table(!is.na(mps$mps2_25te))
hist(mps$mps2_25te)

## assocation with phenotype
cor.test(mps$mps2_25te, df_covs$GESTBIR)
model1 <- lm(scale(GESTBIR) ~ scale(mps$mps2_25te) +
                 GENDER +
                 CD8T + CD4T + NK + Bcell + Mono + Gran + nRBC +
                 C1 + C2 + C3 + C4 + C5 +
                 Sample_Plate,
             data = df_covs)
summary(model1)

## association with outcome
cor.test(mps$mps2_25te, df_covs$cbcl_total_problems)
model1 <- lm(scale(cbcl_total_problems) ~ scale(mps$mps2_25te) + cbcl_age +
                 GENDER +
                 CD8T + CD4T + NK + Bcell + Mono + Gran + nRBC +
                 C1 + C2 + C3 + C4 + C5 +
                 Sample_Plate,
             data = df_covs)
summary(model1)

## 50/50

## the N of the test set
summary(mps$mps2_50te)
table(!is.na(mps$mps2_50te))
hist(mps$mps2_50te)

## assocation with phenotype
cor.test(mps$mps2_50te, df_covs$GESTBIR)
model1 <- lm(scale(GESTBIR) ~ scale(mps$mps2_50te) +
                 GENDER +
                 CD8T + CD4T + NK + Bcell + Mono + Gran + nRBC +
                 C1 + C2 + C3 + C4 + C5 +
                 Sample_Plate,
             data = df_covs)
summary(model1)

## association with outcome
cor.test(mps$mps2_50te, df_covs$cbcl_total_problems)
model1 <- lm(scale(cbcl_total_problems) ~ scale(mps$mps2_50te) + cbcl_age +
                 GENDER +
                 CD8T + CD4T + NK + Bcell + Mono + Gran + nRBC +
                 C1 + C2 + C3 + C4 + C5 +
                 Sample_Plate,
             data = df_covs)
summary(model1)

## 75/25

## the N of the test set
summary(mps$mps2_75te)
table(!is.na(mps$mps2_75te))
hist(mps$mps2_75te)

## assocation with phenotype
cor.test(mps$mps2_75te, df_covs$GESTBIR)
model1 <- lm(scale(GESTBIR) ~ scale(mps$mps2_75te) +
                 GENDER +
                 CD8T + CD4T + NK + Bcell + Mono + Gran + nRBC +
                 C1 + C2 + C3 + C4 + C5 +
                 Sample_Plate,
             data = df_covs)
summary(model1)

## association with outcome
cor.test(mps$mps2_75te, df_covs$cbcl_total_problems)
model1 <- lm(scale(cbcl_total_problems) ~ scale(mps$mps2_75te) + cbcl_age +
                 GENDER +
                 CD8T + CD4T + NK + Bcell + Mono + Gran + nRBC +
                 C1 + C2 + C3 + C4 + C5 +
                 Sample_Plate,
             data = df_covs)
summary(model1)

# MPS3

## data.subset
training_adhd_50te <- df_covs[ in.train_adhd,]
testing_adhd_50te  <- df_covs[-in.train_adhd,]
training_mps_50te <- mps1[ in.train_adhd,]
testing_mps_50te  <- mps1[-in.train_adhd,]

# writing a little loop to get correlations across all mpss
mps3 <- names(mps1)[grep('mps3_',names(mps1))]
mps3_names <- mps3[-1] # we dont need id in there)

## first, we need to select which p-value threshold work best
pvalue_3 <- matrix(NA, length(mps3_names), 3)


for (i in 1:length(mps3_names)){
    pvalue_3[i,1] <- mps3_names[i]
    pvalue_3[i,2] <- summary(lm(scale((training_mps_50te[,mps3_names[i]])) ~ scale(training_adhd_50te$GESTBIR)))$coef[2,1]
    pvalue_3[i,3] <- summary(lm(scale((training_mps_50te[,mps3_names[i]])) ~ scale(training_adhd_50te$GESTBIR)))$coef[2,4]}

pvalue_3[(pvalue_3[,3] == min(pvalue_3[,3])),]
# mps3_P5e-11 is the best one

summary(testing_mps_50te[, 'mps3_P5e-11'])
table(!is.na(testing_mps_50te[, 'mps3_P5e-11']))
hist(testing_mps_50te[, 'mps3_P5e-11'])

## assocation with phenotype
cor.test(testing_mps_50te[, 'mps3_P5e-11'], testing_adhd_50te$GESTBIR)
model1 <- lm(scale(GESTBIR) ~ scale(testing_mps_50te[, 'mps3_P5e-11']) +
                 GENDER +
                 CD8T + CD4T + NK + Bcell + Mono + Gran + nRBC +
                 C1 + C2 + C3 + C4 + C5 +
                 Sample_Plate,
             data = testing_adhd_50te)
summary(model1)

## association with outcome
cor.test(testing_mps_50te[, 'mps3_P5e-11'],
         testing_adhd_50te$cbcl_total_problems)
model1 <- lm(scale(cbcl_total_problems) ~ scale(testing_mps_50te[, 'mps3_P5e-11']) + cbcl_age +
                 GENDER +
                 CD8T + CD4T + NK + Bcell + Mono + Gran + nRBC +
                 C1 + C2 + C3 + C4 + C5 +
                 Sample_Plate,
             data = testing_adhd_50te)
summary(model1)
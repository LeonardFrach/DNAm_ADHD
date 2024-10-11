library(stringr)
wdir_DNAm <- "L:/data/genetic_data/MoBa_methylation/" 
vars <- list.files(wdir_DNAm)


# data sets 1 and 2 contain multiple batches


for(i in c(1:2)){
  
  j = 0
  while (j < 2) {
    j = j + 1
    temp <- read.csv(paste0(wdir_DNAm, word(vars[i]), "/QC/results/", word(vars[i]), "-0", j, "_9_estimated_cell_proportions.csv"))
    colnames(temp)[1] <- "Sample_Name"
    temp$Sample_Name <- toupper(temp$Sample_Name)
    assign(paste0("cell_", i,"_", j), temp) 
    rm(temp)
    
    gc()
  }
}

# data sets 6 and 7 contain one batch (no index j needed)
# NOTE: batch m006 is on position 4 of the 'batches' variable

for(i in c(6,7)){
  
  temp <- read.csv(paste0(wdir_DNAm, word(vars[i]), "/QC/results/", word(vars[i]), "_9_estimated_cell_proportions.csv"))
  colnames(temp)[1] <- "Sample_Name"
  temp$Sample_Name <- toupper(temp$Sample_Name)  
  assign(paste0("cell_", i), temp) 
  rm(temp)
  
  gc()
}
  
 
# data sets 4, 8 and 9, 10 contain parents, so there is a different data structure
j = 0
while (j < 5) {
  for(i in c(4)){
    
    j = j + 1
    temp <- read.csv(paste0(wdir_DNAm, word(vars[i]), "/QC/child/results/", word(vars[i]), "_child-0", j, "_9_estimated_cell_proportions.csv"))
    colnames(temp)[1] <- "Sample_Name"
    temp$Sample_Name <- toupper(temp$Sample_Name)    
    assign(paste0("cell_", i,"_", j), temp) 
    rm(temp)
    
    gc()
  }
}


j = 0
while (j < 4) {
  for(i in c(8)){
    
    j = j + 1
    temp <- read.csv(paste0(wdir_DNAm, word(vars[i]), "/QC/child/results/", word(vars[i]), "_child-0", j, "_9_estimated_cell_proportions.csv"))
    colnames(temp)[1] <- "Sample_Name"
    temp$Sample_Name <- toupper(temp$Sample_Name)    
    assign(paste0("cell_", i,"_", j), temp) 
    rm(temp)
    
    gc()
    
  }
}


for(i in c(9:10)){
  
    j = 0
  while (j < 3) {
    j = j + 1
    temp <- read.csv(paste0(wdir_DNAm, word(vars[i]), "/QC/child/results/", word(vars[i]), "_child-0", j, "_9_estimated_cell_proportions.csv"))
    colnames(temp)[1] <- "Sample_Name"
    temp$Sample_Name <- toupper(temp$Sample_Name)    
    assign(paste0("cell_", i,"_", j), temp) 
    rm(temp)
    
    gc()
    
  }
}


# merge cell type proportions of each sample
cells_merged <- rbind(cell_1_1, cell_1_2, cell_2_1, cell_2_2, cell_4_1, cell_4_2, cell_4_3, cell_4_4, cell_4_5,
                    cell_6, cell_7, cell_8_1, cell_8_2, cell_8_3, cell_8_4, cell_9_1, cell_9_2, cell_9_3,
                    cell_10_1, cell_10_2, cell_10_3)


# merge with data
setwd("N:/durable/people/Leo/DNAm_ADHD")
load("data/MPS_PGS_NPR_pheno.RData") 

alldata <- dplyr::left_join(alldata, cells_merged)

save(alldata, file = "data/MPS_PGS_NPR_pheno_cell.RData")

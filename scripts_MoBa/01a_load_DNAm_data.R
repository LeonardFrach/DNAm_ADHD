## read  and merge methylation data ##

library(stringr)
library(dplyr)
wdir <- "L:/data/genetic_data/MoBa_methylation/" 
vars <- list.files(wdir)


# data sets 1 and 2 contain multiple batches
gc()

for(i in c(1:2)){
  j = 0
  while (j < 2) {
    j = j + 1
    temp <- readRDS(paste0(wdir_DNAm, word(vars[i]), "/QC/results/", word(vars[i]), "-0", j, "_13_bmiqed_beta_values.rds"))
    temp <- temp[1,]
    temp$Sample_Name <- names(temp) 
    temp <- as.data.frame(temp)
    names(temp) <- toupper(names(temp))
    names(temp) <- gsub("X", "", names(temp))
    assign(paste0("DNAm_", i,"_", j), temp) 
    rm(temp)
  }
}


# data set 5 contains multiple batches
gc()

for(i in c(5)){
  j = 0
  while (j < 3) {
    j = j + 1
    temp <- readRDS(paste0(wdir, word(vars[i]), "/QC/results/", word(vars[i]), "-0", j, "_13_bmiqed_beta_values.rds"))
    temp <- temp[1,]
    temp$Sample_Name <- names(temp) 
    temp <- as.data.frame(temp)
    names(temp) <- toupper(names(temp))
    names(temp) <- gsub("X", "", names(temp))
    assign(paste0("DNAm_", i,"_", j), temp) 
    rm(temp)
  }
}

gc()

# data sets 6 and 7 contain one batch (no index j needed)
for(i in c(6,7)){
  temp <- readRDS(paste0(wdir, word(vars[i]), "/QC/results/", word(vars[i]), "_13_bmiqed_beta_values.rds"))
  temp <- temp[1,]
  temp$Sample_Name <- names(temp) 
  temp <- as.data.frame(temp)
  names(temp) <- toupper(names(temp))
  names(temp) <- gsub("X", "", names(temp))
  assign(paste0("DNAm_", i), temp) 
  rm(temp)
}

gc()

# data sets 4, 8 and 9,10 contain parents, so there is a different data structure
j = 0
while (j < 5) {
  for(i in c(4)){
    j = j + 1
    temp <- readRDS(paste0(wdir, word(vars[i]), "/QC/child/results/", word(vars[i]), "_child-0", j, "_13_bmiqed_beta_values.rds"))
    temp <- temp[1,]
    temp$Sample_Name <- names(temp) 
    temp <- as.data.frame(temp)
    names(temp) <- toupper(names(temp))
    names(temp) <- gsub("X", "", names(temp))
    assign(paste0("DNAm_", i,"_", j), temp) 
    rm(temp)
  }
}

gc()

j = 0
while (j < 9) {
  for(i in c(4)){
    j = j + 1
    temp <- readRDS(paste0(wdir, word(vars[i]), "/QC/parent/results/", word(vars[i]), "_parent-0", j, "_13_bmiqed_beta_values.rds"))
    temp <- temp[1,]
    temp$Sample_Name <- names(temp) 
    temp <- as.data.frame(temp)
    names(temp) <- toupper(names(temp))
    names(temp) <- gsub("X", "", names(temp))
    assign(paste0("DNAm_par_", i,"_", j), temp) 
    rm(temp)
  }
}

gc()

for(i in c(4)){
  j = 10
  temp <- readRDS(paste0(wdir, word(vars[i]), "/QC/parent/results/", word(vars[i]), "_parent-", j, "_13_bmiqed_beta_values.rds"))
  temp <- temp[1,]
  temp$Sample_Name <- names(temp) 
  temp <- as.data.frame(temp)
  names(temp) <- toupper(names(temp))
  names(temp) <- gsub("X", "", names(temp))
  assign(paste0("DNAm_par_", i,"_", j), temp) 
  rm(temp)
}


gc()

j = 0
while (j < 4) {
  for(i in c(8)){
    j = j + 1
    temp <- readRDS(paste0(wdir, word(vars[i]), "/QC/child/results/", word(vars[i]), "_child-0", j, "_13_bmiqed_beta_values.rds"))
    temp <- temp[1,]
    temp$Sample_Name <- names(temp) 
    temp <- as.data.frame(temp)
    names(temp) <- toupper(names(temp))
    names(temp) <- gsub("X", "", names(temp))
    assign(paste0("DNAm_", i,"_", j), temp) 
    rm(temp)
  }
}


gc()

j = 0
while (j < 7) {
  for(i in c(8)){
    j = j + 1
    temp <- readRDS(paste0(wdir, word(vars[i]), "/QC/parent/results/", word(vars[i]), "_parent-0", j, "_13_bmiqed_beta_values.rds"))
    temp <- temp[1,]
    temp$Sample_Name <- names(temp) 
    temp <- as.data.frame(temp)
    names(temp) <- toupper(names(temp))
    names(temp) <- gsub("X", "", names(temp))
    assign(paste0("DNAm_par_", i,"_", j), temp) 
    rm(temp)
  }
}


gc()

for(i in c(9:10)){
  j = 0
  while (j < 3) {
    j = j + 1
    temp <- readRDS(paste0(wdir, word(vars[i]), "/QC/child/results/", word(vars[i]), "_child-0", j, "_13_bmiqed_beta_values.rds"))
    temp <- temp[1,]
    temp$Sample_Name <- names(temp) 
    temp <- as.data.frame(temp)
    names(temp) <- toupper(names(temp))
    names(temp) <- gsub("X", "", names(temp))
    assign(paste0("DNAm_", i,"_", j), temp) 
    rm(temp)
  }
}


for(i in c(9)){
  temp <- readRDS(paste0(wdir, word(vars[i]), "/QC/parent/results/", word(vars[i]), "_parent", "_13_bmiqed_beta_values.rds"))
  temp <- temp[1,]
  temp$Sample_Name <- names(temp) 
  temp <- as.data.frame(temp)
  names(temp) <- toupper(names(temp))
  names(temp) <- gsub("X", "", names(temp))
  assign(paste0("DNAm_par_", i), temp) 
  rm(temp)
}


gc()

Sample_Names <- c(names(DNAm_1_1), names(DNAm_1_2),names(DNAm_2_1),names(DNAm_2_2),
                 names(DNAm_4_1),names(DNAm_4_2),names(DNAm_4_3),names(DNAm_4_4),names(DNAm_4_5),
                 names(DNAm_5_1),names(DNAm_5_2),names(DNAm_6),names(DNAm_7),
                 names(DNAm_8_1),names(DNAm_8_2),names(DNAm_8_3), names(DNAm_8_4),
                 names(DNAm_9_1),names(DNAm_9_2),names(DNAm_9_3),
                 names(DNAm_10_1),names(DNAm_10_2),names(DNAm_10_3),
                 names(DNAm_par_8_1), names(DNAm_par_8_2), names(DNAm_par_8_3), names(DNAm_par_8_4), names(DNAm_par_8_5), names(DNAm_par_8_6), names(DNAm_par_8_7),
                 names(DNAm_par_4_1), names(DNAm_par_4_2), names(DNAm_par_4_3), names(DNAm_par_4_4), names(DNAm_par_4_5), 
                 names(DNAm_par_4_6), names(DNAm_par_4_7), names(DNAm_par_4_8), names(DNAm_par_4_9), names(DNAm_par_4_10),
                 names(DNAm_par_9))

save(Sample_Names, file = "L:/people/Leo/ADHD/data/Sample_Names_QCed.RData")
load("L:/people/Leo/ADHD/data/merged_NPR_DNAm.RData")

# merge QCed and DNAm data with pheno
merged_QCed <- merged_no_dups[merged_no_dups$Sample_Name %in% Sample_Names, ]


unique_id_child <- merged_QCed %>% filter(ROLE == "Child")
unique_id_mother <- merged_QCed %>% filter(ROLE == "Mother")
unique_id_father <- merged_QCed %>% filter(ROLE == "Father")

table(unique_id_child$received_dx_F42_npr_child) # yes 63, no 6973   
table(unique_id_father$received_dx_F42_npr_father) # yes 11, no 3254
table(unique_id_mother$received_dx_F42_npr_mother) # yes 12, 3762


# total = 86 cases, 13989 controls

table(unique_id_child$Sample_Group) # 450k = 1825, EPICv1 = 5215
table(unique_id_father$Sample_Group) # EPICv1 = 3267 
table(unique_id_mother$Sample_Group) # EPICv1 = 3775 

# total 450k = 1825, EPICv1 = 12257  
gc()
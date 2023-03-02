library(dplyr)
library(GSVA)
library(GSEABase)

# saveRDS(GSE49710_GSVA_Pheno_Add, 
#         file="../leave one out/data/GSE49710_MES_Low_High.rds")

GSE49710_MES_Low_High <- readRDS(file="../data/GSE49710_MES_Low_High.rds")

MES_DF_1 <- GSE49710_MES_Low_High %>%
group_by(SlidingColumn) %>%
summarise_at(vars(List_RNA_Normoxia_ADRN_vs_MES_Down.txt), 
             list(List_RNA_Normoxia_ADRN_vs_MES_Down= mean))

MES_DF_2 <- GSE49710_MES_Low_High %>% group_by(SlidingColumn) %>% 
  dplyr::summarize(
    n=n(), 
    median_List_RNA_Normoxia_ADRN_vs_MES_Down=median(List_RNA_Normoxia_ADRN_vs_MES_Down.txt),
    mean_List_RNA_Normoxia_ADRN_vs_MES_Down=mean(List_RNA_Normoxia_ADRN_vs_MES_Down.txt),
    sd_List_RNA_Normoxia_ADRN_vs_MES_Down=sd(List_RNA_Normoxia_ADRN_vs_MES_Down.txt)) 


# Select and rename columns with MES bins in HR cohort

GSE49710_MES_selected <- GSE49710_MES_Low_High %>% 
  dplyr::select(title, List_RNA_Normoxia_ADRN_vs_MES_Down.txt, SlidingColumn)

GSE49710_MES_renamed <- GSE49710_MES_selected %>% 
  dplyr::rename(MES_GSVA_values=List_RNA_Normoxia_ADRN_vs_MES_Down.txt)

GSE49710_MES_renamed <- GSE49710_MES_renamed %>% 
  dplyr::rename(MES_GSVA_groups=SlidingColumn)


library(data.table)
library(stringr)
library(tibble)
# Evaluate how long this takes
start_time <- Sys.time()

# Read gene expression df as in "Survival MES.Rmd"
GSE49710 <- readRDS("../data/GSE49710.rds")
GSE49710_GeneMatrix <- as.matrix(GSE49710)

# p_1234 <- readRDS("../data/p_1234.rds")

##############################
###### 1) Read lists from file
##############################

RNA_ADRN_Norm_vs_Hypo_Down_All_1234 <- readRDS("../data/RNA_ADRN_Norm_vs_Hypo_Down_All_1234.rds")
intersect_mes_GSE49711_rownames     <- readRDS("../data/intersect_mes_GSE49711_rownames.rds")
intersect_mes_GSE49711_Hypo_Up      <- readRDS("../data/intersect_mes_GSE49711_Hypo_Up.rds")
intersect_mes_hypo_up_94            <- readRDS("../data/intersect_mes_hypo_up_94.rds")

# strict_mes_1375            <- readRDS("../data/strict_mes_1375.rds")
# strict_hypo_up_460           <- readRDS("../data/strict_hypo_up_460.rds")
intersect_mes_hypo_up_94   <- readRDS("../data/intersect_mes_hypo_up_94.rds")

# Pick genes at random
set.seed(1)
gene_list_random <- sample(intersect_mes_GSE49711_rownames, 6)

################################################
###### 2) Construct vectors to store information
################################################

gene_vector <- c()
gene_number = 0

# The p-value vectors will only be empty after entering the gene list for loop
pval_bucket_smaller <- c()
pval_bucket_bigger <- c()

gene_vector_smaller <- c()
gene_vector_bigger <- c()

pval.vec = c()
rank.vec <- c()

p_val_bonferroni_vector <- c()
my_gene_vector <- c()

# Vectors New GSVA scores
scores_high_group_gene_name             <- c()
scores_high_group_old_score_high        <- c()
scores_high_group_old_score_low         <- c()
scores_high_group_new                   <- c()

scores_low_group_gene_name              <- c()
scores_low_group_old_score_high         <- c()
scores_low_group_old_score_low          <- c()
scores_low_group_new                    <- c()

scores_both_groups_gene_name            <- c()
scores_both_groups_old_score_high       <- c()
scores_both_groups_old_score_low        <- c()
scores_both_groups_new_high             <- c()
scores_both_groups_new_low              <- c()


###########################################
###### 3) Vectors for variables of interest
###########################################

# Vectors Gene Expression Death.from.disease
gene_expression_vector_deceased        <- c() 
my_gene_vector_deceased                <- c() 

# Vectors Gene Expression Progression
gene_expression_vector_progression     <- c() 
my_gene_vector_progression             <- c() 

# Vectors Gene Expression os_bin_ch1 
gene_expression_vector_os_bin_ch1      <- c()
my_gene_vector_os_bin_ch1              <- c()

# Vectors Gene Expression efs_bin_ch1 
gene_expression_vector_efs_bin_ch1     <- c() 
my_gene_vector_efs_bin_ch1             <- c()

# MES groups 1
gene_expression_vector_mes_gsva_groups <- c() 
my_gene_vector_mes_gsva_groups         <- c() 

# MES groups 2
gene_expression_vector_mes             <- c() 
my_gene_vector_mes                     <- c() 



################################################################
###### 4) Decide list and elements in the list to enter for loop
################################################################

for (my_gene in intersect_mes_hypo_up_94) {
  
  # In case the loop breaks for any gene,
  # tryCatch continues the loop iteration after it breaks
  tryCatch({
  
  gene_number = gene_number + 1
  print (paste0("Starting with gene_number as: ", gene_number))
  print (paste0("Gene is: ", my_gene))
  
  gene_vector <- c(gene_vector, my_gene)
  
  #######################################################
  ###### 5.1) Construct GSVA w/o current gene in gene set
  #######################################################
  
  # Remove current gene from list
  partial_list <- intersect_mes_hypo_up_94[!grepl(my_gene, unlist(intersect_mes_hypo_up_94))]
  
  # Name gene set with gene name that was removed
  # This gene set name goes as column name in DF
  partial_list_appended <- append(partial_list, paste0("partial_list_", my_gene), 0)
  
  # Write gene set to read, in the next step, for GSVA calculation
  # Interesting that transposing the vector forced printing horizontally instead of vertically.
  write.table(t(partial_list_appended), file = paste0("../data/", "leave_", my_gene, "_out.txt"),
              sep = "\t",
              row.names = F,
              col.names = F,
              quote = F)
  
  print(paste0("Read gene set for GSVA calculation, from the 'data' folder for gene ", my_gene))
  leave_gene_out_gmt <- getGmt(paste0("../data/", "leave_", my_gene, "_out.txt"))
  
  ##############################
  ###### 5.2) Call GSVA function
  ##############################
  
  print(paste0("Calculate GSVA for samples excluding gene ", my_gene, " from the current gene set"))
  GSE49710_GSVA <- gsva(GSE49710_GeneMatrix,
                          leave_gene_out_gmt,
                          min.sz=1, max.sz=Inf,
                          verbose=TRUE)
  
  print(paste0("Bind GSVA and the gene matrix for gene ", my_gene))
  GSVA_GeneExpression_Matrix <- rbind(GSE49710_GSVA, GSE49710_GeneMatrix)
  
  
  #####################################################
  ###### 6) Transpose and Name Column with Sample Names
  #####################################################

  GSVA_GSE49710_DF_transpose <- t(GSVA_GeneExpression_Matrix)
  GSVA_GSE49710_DF_transpose <- as.data.frame(GSVA_GSE49710_DF_transpose)
  GSVA_GSE49710_DF_transpose["Sample_ID"] <- rownames(GSVA_GSE49710_DF_transpose)
  
  #################################################
  ###### 7) Load Metadata and Phenotype Information
  #################################################

  phenoGSE49710 <-read.delim(
  file ="../data/GSE49710_3.txt")

  dt1 <- data.table(phenoGSE49710, key = "Group_Accession") 
  dt2 <- data.table(GSVA_GSE49710_DF_transpose, key = "Sample_ID")

  GSE49710_GSVA_Pheno <- dt1[dt2]
  
  ######################################################
  ###### 8) Fix Sample Name Information from Pheno File
  ######################################################

  # Split Each Patient ID
  df_Title <- str_split_fixed(GSE49710_GSVA_Pheno$Title, " ", 4)
  # Construct Vector Maintaining parts of interest
  title_vector <- paste0(df_Title[,1], "_", df_Title[,2], df_Title[,4])
  # Add vector to DF, after column of interest
  GSE49710_GSVA_Pheno <- add_column(GSE49710_GSVA_Pheno, Title_Compatible = title_vector, .after = "Title")
  
  ############################################
  ###### 9) Metadata and Phenotype Information
  ############################################

  ####################
  ###### 9.1) Metadata
  ####################

  sampleInfo_GSE62564 <-read.delim(
  file ="../data/sampleInfo.txt")

  dt1 <- data.table(sampleInfo_GSE62564, key = "title") 
  dt2 <- data.table(GSE49710_GSVA_Pheno, key = "Title_Compatible")

  GSE49710_GSVA_Metadata_Complete <- dt1[dt2]
  
  
  #################################
  ############ 10) Select HR Samples
  #################################

  GSE49710_GSVA_Pheno_Subset <- GSE49710_GSVA_Metadata_Complete[ which(
  # GSE49710_GSVA_Metadata_Complete$High.risk == '0' |
  GSE49710_GSVA_Metadata_Complete$High.risk == '1'  ), ]
  
  print("Select gene and geneset of interest")
  GSE49710_GSVA_Pheno_Subset <- GSE49710_GSVA_Pheno_Subset %>% dplyr::select(1:78, all_of(my_gene))
  
  
  ############################################
  ############ 11) Merge MES partitions with DF
  ############################################
  
  dt1 <- data.table(GSE49710_GSVA_Pheno_Subset, key = "title") 
  dt2 <- data.table(GSE49710_MES_renamed, key = "title")

  GSE49710_HR_MES_Groups <- dt1[dt2]
  
  #####################################################
  ############ 12) Current means of the groups of tumors
  #####################################################
  
  print(paste0("Current Mean GSVA score, including gene ", my_gene, 
               " in list for High MES is: ", MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[1]))
  
  print(paste0("Current Mean GSVA score, including gene ", my_gene, 
             " in list for Low MES is: ", MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[2]))
  
  ###################################################################
  ############ 13) Calculate New MES Mean for partial_list w/o my_gene
  ###################################################################
  
  print("Constructing DFs for partial list calculations")
  my_gene_mean <- GSE49710_HR_MES_Groups %>%
  group_by(MES_GSVA_groups) %>%
  summarise_at(vars(!!(paste0("partial_list_", my_gene))), list(partial_list_my_gene = mean))

  print(paste0("Mean GSVA score of ", my_gene, " for High MES is: ", my_gene_mean$partial_list_my_gene[1]))
  print(paste0("Mean GSVA score of ", my_gene, " for Low MES is: ", my_gene_mean$partial_list_my_gene[2]))
  
  
  ##########################################################################
  ############ 14 ) If statement create vectors considering the partial list
  ##########################################################################
  
  # # If gene decreases MES scores of High MES group
  # if (  my_gene_mean$partial_list_my_gene[1] < MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[1]  ) {
  # print(paste0("Gene ", my_gene, " decreases High Group score"))
  # new_score_vector_high <- c(new_score_vector_high, my_gene_mean$partial_list_my_gene[1])
  # old_score_vector_high <- c(old_score_vector_high, MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[1])
  # old_score_vector_low  <- c(old_score_vector_low, MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[2])
  # my_gene_vector_high   <- c(my_gene_vector_high, my_gene)
  # }
  # 
  # # Try and make statement based on High or Low MES group
  # # ifelse(my_gene_mean$MES_GSVA_groups=="High", 
  # #        print("It is true"), 
  # #        print("It is false"))
  # 
  # # If gene increases MES scores of Low MES group
  # if (  my_gene_mean$partial_list_my_gene[2] > MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[2]  ) {
  # print(paste0("Gene ",my_gene, " increases Low Group Score"))
  # new_score_vector_low  <- c(new_score_vector_low, my_gene_mean$partial_list_my_gene[2])
  # my_gene_vector_low    <- c(my_gene_vector_low, my_gene)
  # }
  # 
  # # If gene decreases MES scores of High MES group and increases MES scores of Low MES group
  # if (  (my_gene_mean$partial_list_my_gene[1] < MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[1])  &
  #       (  my_gene_mean$partial_list_my_gene[2] > MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[2] )) {
  # print(paste0("Gene ",my_gene, " decreases High Group score and increases Low Group score"))
  # gene_vector_decrease_high_increase_low <- c(vector_decrease_high_increase_low, my_gene)
  # }
  
  
  
  # 1) Prepare DF genes scores_high_group
  
  if (my_gene_mean$partial_list_my_gene[1] < MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[1]) {
    
    # A) Write gene name
    scores_high_group_gene_name    <- c(scores_high_group_gene_name, my_gene)
    # B) Write old score of high group (0.0029)
    scores_high_group_old_score_high <- c(scores_high_group_old_score_high, MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[1])
    # C) Write old score of low group (-0.47)
    scores_high_group_old_score_low <- c(scores_high_group_old_score_low, MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[2])
    # D) Write new lower score of high group
    scores_high_group_new <- c(scores_high_group_new, my_gene_mean$partial_list_my_gene[1])
    
  }
  
  
  
  # 2) Prepare DF genes scores_low_group
  
  if (my_gene_mean$partial_list_my_gene[1] > MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[2]) {
    
    # A) Write gene name
    scores_low_group_gene_name    <- c(scores_low_group_gene_name, my_gene)
    # B) Write old score of high group (0.0029)
    scores_low_group_old_score_high <- c(scores_low_group_old_score_high, MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[1])
    # C) Write old score of low group (-0.47)
    scores_low_group_old_score_low <- c(scores_low_group_old_score_low,   MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[2])
    # D) Write new higher score of low group
    scores_low_group_new <- c(scores_low_group_new, my_gene_mean$partial_list_my_gene[2])
    
  }
  
  
   
  # 3) Prepare DF genes scores_both_groups
  
  if ((my_gene_mean$partial_list_my_gene[1] < MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[1])&
      
      (my_gene_mean$partial_list_my_gene[1] > MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[2])) {
    
    # A) Write gene name
    scores_both_groups_gene_name                <- c(scores_both_groups_gene_name, my_gene)
    # B) Write old score of high group (0.0029)
    scores_both_groups_old_score_high           <- c(scores_both_groups_old_score_high, MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[1])
    # C) Write old score of low group (-0.47)
    scores_both_groups_old_score_low            <- c(scores_both_groups_old_score_low, MES_DF_1$List_RNA_Normoxia_ADRN_vs_MES_Down[2])
    # D) Write new high scores
    scores_both_groups_new_high                 <- c(scores_both_groups_new_high, my_gene_mean$partial_list_my_gene[1])
    # E) Write new low scores
    scores_both_groups_new_low                  <- c(scores_both_groups_new_low, my_gene_mean$partial_list_my_gene[2])

  }

  
  
  ###############################################################################
  ############ 15) If statement create vectors considering expression of the gene
  ###############################################################################
  
  #####################################
  ############ 15.1) Death.from.disease
  #####################################
  
  print("Construct DFs to evaluate gene in Death.from.disease calculations")
  death_from_disease_df <- GSE49710_HR_MES_Groups %>%
  group_by(Death.from.disease) %>%
    mutate(Count = n()) %>%
  summarise_at(vars(!!(my_gene)), list(my_gene = mean))
  
  death_from_disease_renamed_df <- death_from_disease_df
  
  print("Rename my_gene in DF, in the hope to improve clarity")
  names(death_from_disease_renamed_df)[names(death_from_disease_renamed_df) == 'my_gene'] <- my_gene

  print("Evaluate if statement for death_from_disease")
  if (  death_from_disease_renamed_df[1,my_gene] > death_from_disease_renamed_df[2,my_gene]  ) {
  print(paste0("Gene ", my_gene, " expression is higher in deceased"))
  gene_expression_vector_deceased <- c(gene_expression_vector_deceased, death_from_disease_renamed_df[1,my_gene])
  my_gene_vector_deceased <- c(my_gene_vector_deceased, my_gene)
  }
  
  ##############################
  ############ 15.2) Progression
  ##############################
  
  print("Construct DFs to evaluate gene in Progression calculations")
  progression_df <- GSE49710_HR_MES_Groups %>%
  group_by(Progression) %>%
    mutate(Count = n()) %>%
  summarise_at(vars(!!(my_gene)), list(my_gene = mean))
  
  progression_renamed_df <- progression_df
  
  print("Rename my_gene in DF, in the hope to improve clarity")
  names(progression_renamed_df)[names(progression_renamed_df) == 'my_gene'] <- my_gene

  print("Evaluate if statement for progression")
  if (  progression_renamed_df[1,my_gene] > progression_renamed_df[2,my_gene]  ) {
  print(paste0("Gene ", my_gene, " expression is higher in Progression"))
  gene_expression_vector_progression <- c(gene_expression_vector_progression, progression_renamed_df[1,my_gene])
  my_gene_vector_progression <- c(my_gene_vector_progression, my_gene)
  }
  
  #############################
  ############ 15.3) os_bin_ch1
  #############################
  
  print("Construct DFs to evaluate gene in OS calculations")
  os_bin_ch1_df <- GSE49710_HR_MES_Groups %>%
  group_by(os.bin.ch1) %>%
    mutate(Count = n()) %>%
  summarise_at(vars(!!(my_gene)), list(my_gene = mean))
  
  os_bin_ch1_renamed_df <- os_bin_ch1_df
  
  print("Rename my_gene in DF, in the hope to improve clarity")
  names(os_bin_ch1_renamed_df)[names(os_bin_ch1_renamed_df) == 'my_gene'] <- my_gene

  print("Evaluate if statement for os_bin_ch1 overall survival")
  if (  os_bin_ch1_renamed_df[1,my_gene] > os_bin_ch1_renamed_df[2,my_gene]  ) {
  print(paste0("Gene ", my_gene, " expression is higher in 0"))
  gene_expression_vector_os_bin_ch1 <- c(gene_expression_vector_os_bin_ch1, os_bin_ch1_renamed_df[1,my_gene])
  my_gene_vector_os_bin_ch1 <- c(my_gene_vector_os_bin_ch1, my_gene)
  }
  
  ##############################
  ############ 15.4) efs_bin_ch1
  ##############################
  
  print("Construct DFs to evaluate gene in OS calculations")
  efs_bin_ch1_df <- GSE49710_HR_MES_Groups %>%
  group_by(efs.bin.ch1) %>%
    mutate(Count = n()) %>%
  summarise_at(vars(!!(my_gene)), list(my_gene = mean))
  
  efs_bin_ch1_renamed_df <- efs_bin_ch1_df
  
  print("Rename my_gene in DF, in the hope to improve clarity")
  names(efs_bin_ch1_renamed_df)[names(efs_bin_ch1_renamed_df) == 'my_gene'] <- my_gene

  print("Evaluate if statement for efs_bin_ch1 event free survival")
  if (  efs_bin_ch1_renamed_df[2,my_gene] > efs_bin_ch1_renamed_df[1,my_gene]  ) {
  print(paste0("Gene ", my_gene, " expression is higher in 1"))
  gene_expression_vector_efs_bin_ch1 <- c(gene_expression_vector_efs_bin_ch1, efs_bin_ch1_renamed_df[2,my_gene])
  my_gene_vector_efs_bin_ch1 <- c(my_gene_vector_efs_bin_ch1, my_gene)
  }
  
  ##################################
  ############ 15.5) MES_GSVA_Groups
  ##################################
  
  print("Construct DFs to evaluate gene in OS calculations")
  MES_GSVA_Groups_df <- GSE49710_HR_MES_Groups %>%
  group_by(MES_GSVA_groups) %>%
    mutate(Count = n()) %>%
  summarise_at(vars(!!(my_gene)), list(my_gene = mean))
  
  MES_GSVA_Groups_renamed_df <- MES_GSVA_Groups_df
  
  print("Rename my_gene in DF, in the hope to improve clarity")
  names(MES_GSVA_Groups_renamed_df)[names(MES_GSVA_Groups_renamed_df) == 'my_gene'] <- my_gene

  print("Evaluate if statement for MES_GSVA_Groups")
  if (  MES_GSVA_Groups_renamed_df[2,my_gene] > MES_GSVA_Groups_renamed_df[1,my_gene]  ) {
  print(paste0("Gene ", my_gene, " expression is higher in 1"))
  gene_expression_vector_mes_gsva_groups <- c(gene_expression_vector_mes_gsva_groups, MES_GSVA_Groups_renamed_df[2, my_gene])
  my_gene_vector_mes_gsva_groups <- c(my_gene_vector_mes_gsva_groups, my_gene)
  }
  
  print("Evaluate if statement for MES Groups Low and High")
  # If gene expression in High is lower than in Low
  if (  MES_GSVA_Groups_renamed_df[rlang::sym(my_gene)][1,] < MES_GSVA_Groups_renamed_df[rlang::sym(my_gene)][2,]  ) {
  print(paste0("Gene ", my_gene, " expression in High is lower than Low"))
  gene_expression_vector_mes <- c(gene_expression_vector_mes, MES_GSVA_Groups_renamed_df[rlang::sym(my_gene)][1,])
  my_gene_vector_mes <- c(my_gene_vector_mes, my_gene)
  }
  
  print(paste0("Ending for loop for gene ", my_gene))
  
  # Close tryCatch:
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  # End of gene list for loop:
}

# Evaluate how long this takes
end_time <- Sys.time()
end_time - start_time

################################################
############ 16) Construct vectors post for loop
################################################

# scores_high_group
tmp_matrix_1                   <- cbind(scores_high_group_gene_name, scores_high_group_old_score_high)
tmp_matrix_2                   <- cbind(scores_high_group_old_score_low, scores_high_group_new)
scores_high_group_matrix       <- cbind(tmp_matrix_1, tmp_matrix_2)
write.table(scores_high_group_matrix, 
            file = paste0("../results/", "scores_high_group_matrix.txt"),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F)

# scores_low_group
tmp_matrix_1                   <- cbind(scores_low_group_gene_name, scores_low_group_old_score_high)
tmp_matrix_2                   <- cbind(scores_low_group_old_score_low, scores_low_group_new)
scores_low_group_matrix        <- cbind(tmp_matrix_1, tmp_matrix_2)
write.table(scores_low_group_matrix, 
            file = paste0("../results/", "scores_low_group_matrix.txt"),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F)

# scores_both_groups
tmp_matrix_1                   <- cbind(scores_both_groups_gene_name, scores_both_groups_old_score_high)
tmp_matrix_2                   <- cbind(scores_both_groups_new_high, scores_both_groups_old_score_low )
tmp_matrix_3                   <- cbind(tmp_matrix_2, scores_both_groups_new_low)

scores_both_groups_matrix      <- cbind(tmp_matrix_1, tmp_matrix_3)
write.table(scores_both_groups_matrix, 
            file = paste0("../results/", "scores_both_groups_matrix.txt"),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F)

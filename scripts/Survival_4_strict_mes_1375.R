library("data.table")
library("stringr")
library("tibble")
library("survival")
library("survminer")
library("dplyr")

start_time <- Sys.time()

GSE49710_GSVA_Metadata_Complete <- readRDS("../data/GSE49710_GSVA_Metadata_Complete.rds")

##############################
############ 1) Select HR Samples
##############################

GSE49710_GSVA_Pheno_Subset <- GSE49710_GSVA_Metadata_Complete[ which(
  # GSE49710_GSVA_Metadata_Complete$High.risk == '0' |
  GSE49710_GSVA_Metadata_Complete$High.risk == '1'  ), ]


############################################################################
############ y is needed for the ranking of scores
############ i will be the dataframe index: position of the column in the DF
############################################################################

gene_number = 0

gene_expression_vector_surv_prob_smaller  <- c()
my_gene_vector_surv_prob_smaller          <- c() 
gene_expression_vector_surv_prob_bigger   <- c()
my_gene_vector_surv_prob_bigger           <- c()

partition_vector_low                      <- c()
partition_vector_high                     <- c()
gene_vector                               <- c()

smallest_pval_vector                      <- c()
smallest_pval_vector_assign               <- c()

surv_prob_1_500_vector                    <- c()
surv_prob_2_500_vector                    <- c()
surv_prob_1_1000_vector                   <- c() 
surv_prob_2_1000_vector                   <- c()
surv_prob_1_2000_vector                   <- c()
surv_prob_2_2000_vector                   <- c()

surv_prob_1_500_vector_high               <- c()
surv_prob_1_500_vector_low                <- c()

high_surv_prob_vec_500                    <- c()
low_surv_prob_vec_500                     <- c()

adjusted_pval_vector                      <- c()

surv_prob_string_vec                      <- c()

strict_mes_1375            <- readRDS("../data/strict_mes_1375.rds")
# strict_hypo_up_460         <- readRDS("../data/strict_hypo_up_460.rds")
# intersect_mes_hypo_up_94   <- readRDS("../data/intersect_mes_hypo_up_94.rds")

for (my_gene in strict_mes_1375) {
  
  tryCatch({
    
    gene_number = gene_number + 1
    
    # GSE49710_GSVA_Pheno_Subset <- GSE49710_GSVA_Pheno_Subset %>% dplyr::select(1:77, my_gene)
    
    # View(GSE49710_GSVA_Pheno_Subset)
    
    print(paste0("Starting first for loop for gene list and gene is: ", my_gene, " . Gene number in line is ", gene_number))
    
    y = 0
    pval.vec = c()
    rank.vec <- c()
    
    for (i in 1:175) {
      
      print(paste0("Starting second for loop for number of samples in lowest expression for gene: ", my_gene, "; y is ", y, " and i is ", i))
      
      y = y + 1
      
      
      #################################################
      ############ 1) Order Column / Add Rank column to DF
      #################################################
      
      GSE49710_GSVA_Pheno_Add = GSE49710_GSVA_Pheno_Subset[order(GSE49710_GSVA_Pheno_Subset[[my_gene]]),]
      GSE49710_GSVA_Pheno_Add$Rank2 <-  rank(GSE49710_GSVA_Pheno_Add[[my_gene]])
      
      
      ###########################################################################################
      ############ Using Rank column, create sliding window: Split of MES High and MES Low Tumors
      ###########################################################################################
      
      SlidingVector <- as.factor(
        ifelse(GSE49710_GSVA_Pheno_Add$Rank2 <= y, SlidingVector <- "Low", SlidingVector <- "High"))
      
      
      ##################################################################
      ############ 2) Add partition of High and Low MES scores, after Rank2
      ##################################################################
      
      GSE49710_GSVA_Pheno_Add <- add_column(GSE49710_GSVA_Pheno_Add, 
                                            SlidingColumn = SlidingVector, .after = "Rank2")
      
      
      #########################################################################
      ############ 3) Fit Survival, based on the Split of Low and High MES Samples
      #########################################################################
      
      fit <- survfit(Surv(efs.day.ch1, efs.bin.ch1) ~ GSE49710_GSVA_Pheno_Add$SlidingColumn,
                     data = GSE49710_GSVA_Pheno_Add)
      
      
      #####################################
      ############ 4) p-value and rank vectors
      #####################################
      
      pval.value <- surv_pvalue(fit)$pval
      pval.vec <- c(pval.vec, pval.value)
      
      rank.value <- print(GSE49710_GSVA_Pheno_Add$Rank2[y])
      rank.vec <- c(rank.vec, rank.value)
      
      print(paste0("Ending HR for loop. Number of elements in Low was: ", i, ". Gene is: ", my_gene, " and its number in line is ", gene_number))
      
      # End of HR for loop:
    }
    
    
    #######################################
    ############ 5) P-values and Ranks into DF
    #######################################
    
    p.valRank <- cbind(rank.vec, pval.vec)
    p.valRank <- as.data.frame(p.valRank)
    
    p.valRank$Bonferroni = 
      p.adjust(p.valRank$pval.vec, 
               method = "bonferroni")
    
    p.valRank$BH = 
      p.adjust(p.valRank$pval.vec, 
               method = "BH")
    
    p.valRank$hochberg = 
      p.adjust(p.valRank$pval.vec, 
               method = "hochberg")
    
    p.valRank$holm = 
      p.adjust(p.valRank$pval.vec, 
               method = "holm")
    
    p.valRank$hommel = 
      p.adjust(p.valRank$pval.vec, 
               method = "hommel")
    
    p.valRank$BY = 
      p.adjust(p.valRank$pval.vec, 
               method = "BY")
    
    
    p.valRank_ordered = p.valRank[order(p.valRank["BH"]),]
    p.valRank_ordered_my_gene <- assign(paste0("p.valRank_ordered_",  my_gene), p.valRank[order(p.valRank["BH"]),])
    
    
    #####################################################################
    ############ 6) Number of samples in lowest group based on smallest padj
    #####################################################################
    
    SlidingVectorOutside <- as.factor(
      ifelse(GSE49710_GSVA_Pheno_Add$Rank2 <= p.valRank_ordered[1, "rank.vec"], SlidingVectorOutside <- "Low",
             SlidingVectorOutside <- "High"))
    
    adjusted_pval_vector      <- c(adjusted_pval_vector, p.valRank_ordered[1, "BH"])
    adjusted_pval_vector_gene <- assign(paste0("adjusted_pval_vector_", my_gene), p.valRank_ordered[1, "BH"])
    
    
    ##################################################################
    #### 7) Add partition of High and Low MES scores, after channel_count
    ##################################################################
    
    GSE49710_GSVA_Pheno_Add <- add_column(GSE49710_GSVA_Pheno_Add, 
                                          SlidingColumnOutside = SlidingVectorOutside, 
                                          .after = "channel_count")
    
    GSE49710_GSVA_Pheno_Add$SlidingColumnOutside <- as.factor(GSE49710_GSVA_Pheno_Add$SlidingColumnOutside)
    
    
    
    #########################################################################
    ############ 8) Fit Survival, based on the Split of Low and High MES Samples
    #########################################################################
    
    print(paste0("fit KM for gene ", my_gene))
    fit_my_gene_smallest_pval      <- survfit(Surv(efs.day.ch1, efs.bin.ch1) ~ SlidingVectorOutside, data = GSE49710_GSVA_Pheno_Add)
    fit_smallest_pval_gene         <- assign(paste0("fit_smallest_pval_", my_gene),
                                             survfit(Surv(efs.day.ch1, efs.bin.ch1) ~ SlidingVectorOutside, data = GSE49710_GSVA_Pheno_Add))
    
    survfit_results_df <- broom::tidy(fit_my_gene_smallest_pval)
    
    survfit_table <- summary(fit_my_gene_smallest_pval)$table
    
    print("Assign to object")
    fit_gene_pval <- assign(paste0("fit_", my_gene,"_smallest_pval"), fit_my_gene_smallest_pval)
    
    print(paste0("Assign to entire function. Gene is ", my_gene, " Number ", gene_number, " in line."))
    print("format to call p-val for gene is fit_CYP4X1_smallest_pval")
    fit_gene_pval <- assign(paste0("fit_", my_gene,"_smallest_pval"), 
                            survfit(Surv(efs.day.ch1, efs.bin.ch1) ~ SlidingVectorOutside, data = GSE49710_GSVA_Pheno_Add))
    
    smallest_pval_value <- surv_pvalue(fit_my_gene_smallest_pval)$pval
    smallest_pval_vector <- c(smallest_pval_vector, smallest_pval_value)
    print("format to call smallest pvalue is smallest_pval_my_gene")
    smallest_pval_value <- assign(paste0("smallest_pval_", my_gene), surv_pvalue(fit_my_gene_smallest_pval)$pval)
    smallest_pval_vector_assign <- c(smallest_pval_vector, 
                                     assign(paste0("smallest_pval_", my_gene), surv_pvalue(fit_my_gene_smallest_pval)$pval))
    
    print("Assign partitions of gene expression as High and Low")
    partitions_gene <- assign(paste0("partitions_", my_gene), 
                              GSE49710_GSVA_Pheno_Add %>% group_by(!!SlidingVectorOutside) %>% count())
    # partitions_GBP4[2,2] # Low
    # partitions_GBP4[1,2] # High
    partition_vector_low  <- c(partition_vector_low, as.numeric(partitions_gene[2,2]))
    partition_vector_high <- c(partition_vector_high, as.numeric(partitions_gene[1,2]))
    gene_vector           <- c(gene_vector, my_gene)
    
    
    # Get survival probability of groups at a given time point
    print(paste0("summary.survfit_object_500 for gene ", my_gene))
    summary.survfit_object_500      <- summary(fit_my_gene_smallest_pval, times = 500)
    
    mes_high_surv_prob_500       <- summary_survfit_object_500$surv[1]
    mes_low_surv_prob_500        <- summary_survfit_object_500$surv[2]
    
    mes_high_surv_prob     <-  summary.survfit_object_500$surv[1] # Survival probability of MES High at the specified time point
    mes_low_surv_prob      <-  summary.survfit_object_500$surv[2] # Survival probability of MES Low at the specified time point
    
    high_surv_prob_vec_500 <- c(high_surv_prob_vec_500, mes_high_surv_prob_500)
    low_surv_prob_vec_500  <- c(low_surv_prob_vec_500,  mes_low_surv_prob_500)
    
    assign(paste0("surv_prob_high_500_", my_gene), summary_survfit_object_500$surv[1])
    assign(paste0("surv_prob_low_500_",  my_gene), summary_survfit_object_500$surv[2])
    
    # Print gene and survival probabilities together
    surv_prob_string <- paste0(my_gene, 
                               " ", 
                               get(paste0("surv_prob_high_500_", my_gene)), # This allows to get the value
                               " " ,
                               get(paste0("surv_prob_low_500_", my_gene)))
    
    surv_prob_string_vec <- c(surv_prob_string, surv_prob_string_vec)
    
    # Close tryCatch:
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  print(paste0("Ending gene list for loop. Gene was: ", my_gene))
  
  # End of gene list for loop:
}

print(paste0("It took this long to run code for ", gene_number, " genes", ":"))

end_time <- Sys.time()
end_time - start_time

print(my_gene_vector_surv_prob_smaller)
print(my_gene_vector_surv_prob_bigger)

# rbind(str_split_fixed(surv_prob_string_vec[1], " ", 3),(str_split_fixed(surv_prob_string_vec[2], " ", 3)))

#################################################
############ 9) Plot Survival Curve Outside for loop
#################################################

# ggsurvplot(fit, 
#            data = GSE49710_GSVA_Pheno_Add, 
#            risk.table = TRUE, 
#            pval = TRUE, 
#            pval.coord = c(2000, 1),
#            ylab = c(paste0("Survival Probability ", my_gene)))


######################################
############ 10) qplots Ranks vs. p-values
######################################

qplot (rank.vec, BH, data = p.valRank, ylab = paste0("BH ", my_gene))

qplot (rank.vec, BH, data = p.valRank, ylab = paste0("BH ", my_gene))+
  scale_y_continuous(
    breaks = c(0.05, .1, 0.25, 0.5, 1),
    labels = c("0.05", ".1", "0.25", "0.5", "1")
  )

qplot_BH <- qplot (rank.vec, BH, data = p.valRank, ylab = paste0("BH ", my_gene))+
  scale_y_continuous(
    breaks = c(0.05, .1, 0.25, 0.5, 1),
    labels = c("0.05", ".1", "0.25", "0.5", "1")
  )

# par(mar = c(5, 4, 4, 4) + 0.3)                               # Additional space for second y-axis
# plot(rank.vec, p.valRank$pval.vec, pch = 16, col = 2)        # Create first plot
# par(new = TRUE)                                              # Add new plot
# plot(rank.vec, p.valRank$BH, pch = 17, col = 3,              # Create second plot without axes
#      axes = FALSE, xlab = "", ylab = "")
# axis(side = 4, at = pretty(range(p.valRank$BH)))             # Add second axis
# mtext("BH", side = 4, line = 3)          


qplot (rank.vec, Bonferroni, data = p.valRank, ylab = paste0("Bonferroni ", my_gene))
qplot (rank.vec, hochberg,   data = p.valRank, ylab = paste0("hochberg "  , my_gene))
qplot (rank.vec, holm,       data = p.valRank, ylab = paste0("holm "      , my_gene))
qplot (rank.vec, hommel,     data = p.valRank, ylab = paste0("hommel "    , my_gene))
qplot (rank.vec, BY,         data = p.valRank, ylab = paste0("BY "        , my_gene))


###########################################
############ 11) Gene, Partition and padj DF
###########################################

partition_matrix                    <- cbind(gene_vector, partition_vector_low)
# tmp_pvalue                        <- cbind(smallest_pval_vector, smallest_pval_vector_assign)
partition_matrix                    <- cbind(partition_matrix, partition_vector_high)
smallest_pval_vector_matrix         <- cbind(partition_matrix, smallest_pval_vector)
adjusted_pval_matrix                <- cbind(smallest_pval_vector_matrix, adjusted_pval_vector)

surv_prob_matrix_500                <- cbind(high_surv_prob_vec_500, low_surv_prob_vec_500)

final_surv_prob_matrix_500          <- cbind(adjusted_pval_matrix, surv_prob_matrix_500)

write.table(adjusted_pval_matrix, file = paste0("../results/", "adjusted_pval_matrix_strict_mes_1375.txt"),
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F)


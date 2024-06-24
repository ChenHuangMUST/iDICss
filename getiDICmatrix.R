# Constructing a gene mutation matrix
# Gene mutation matrix, 1 is mutated, 0 is wild
get_mut_matrix <- function (maffile,nonsynonymous = TRUE) {#,percent = 0.01
  maf <- maffile
  mafneed <- maf[, c("Hugo_Symbol", "Tumor_Sample_Barcode", 
                     "Variant_Classification")]
  if (nonsynonymous == TRUE) {
    mafmissense <- mafneed[which(mafneed$Variant_Classification == 
                                   "Missense_Mutation" | mafneed$Variant_Classification == 
                                   "Frame_Shift_Del" | mafneed$Variant_Classification == 
                                   "Frame_Shift_Ins" | mafneed$Variant_Classification == 
                                   "In_Frame_Del" | mafneed$Variant_Classification == 
                                   "Nonsense_Mutation" | mafneed$Variant_Classification == 
                                   "In_Frame_Ins" | mafneed$Variant_Classification == 
                                   "Splice_Site" | mafneed$Variant_Classification == 
                                   "Nonstop_Mutation" | mafneed$Variant_Classification == 
                                   "Translation_Start_Site"), ]
    mafselect <- mafmissense
  }
  else {
    mafselect <- mafneed
  }
  sample_id <- unique(mafselect$Tumor_Sample_Barcode)
  gene_id <- unique(mafselect$Hugo_Symbol)
  mafselectmatrix <- matrix(0, nrow = length(gene_id), ncol = length(sample_id))
  rownames(mafselectmatrix) <- gene_id
  colnames(mafselectmatrix) <- sample_id
  for (i in 1:nrow(mafselect)) {
    gene_one <- mafselect$Hugo_Symbol[i]
    sample_one <- mafselect$Tumor_Sample_Barcode[i]
    mafselectmatrix[gene_one, sample_one] <- 1
  }
  sample_id <- substr(sample_id, 1, 16)
  colnames(mafselectmatrix) <- sample_id
  mafmatrix <- mafselectmatrix
  return(mafmatrix)
}

# Calculating the iDC matrix
calculate_iDC_matrix <- function(S_matrix, W_matrix,mut_matrix, driver_gene, lamda = lamda) {
  S_matrix_filtered <- S_matrix[driver_gene, ]
  mut_matrix_filtered <- mut_matrix[driver_gene, ]
  # Calculating the iDC matrix
  iDC_matrix <- lamda * W_matrix * (t(abs(S_matrix_filtered)) %*% as.matrix(mut_matrix_filtered)) +
    (1 - lamda) * W_matrix
  return(iDC_matrix)
}
# S_matrix: the original S matrix.
# mut_matrix: the original mutation matrix.
# W_matrix: the weight matrix.
# driver_gene: index vector used to select a specific gene from the S matrix and mutation matrix.
# lamda: coefficient used in the calculation, default value is 0.7.
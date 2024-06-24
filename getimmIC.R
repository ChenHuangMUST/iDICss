# Selects genes that are beyond the given number of standard deviations from the mean
select_genes <- function(S_matrix, num_std_dev) {
  # Calculate mean and standard deviation
  mean_weight <- apply(S_matrix, 2, mean)
  sd_weight <- apply(S_matrix, 2, sd)
  # Calculate thresholds
  lower_threshold <- mean_weight - num_std_dev * sd_weight
  upper_threshold <- mean_weight + num_std_dev * sd_weight
  # Initialize a list to hold the selected genes for each independent-component
  IC_top_genes <- list()
  # Iterate over each independent-component
  ICs <- colnames(S_matrix)
  for (IC in ICs) {
    # Select genes where weight is beyond the thresholds for this independent-component
    selected_index <- which(S_matrix[, IC] > upper_threshold[IC] | S_matrix[, IC] < lower_threshold[IC])
    selected_genes <- rownames(S_matrix)[selected_index]
    IC_top_genes[[IC]] <- selected_genes
    # Optionally print the independent-component being processed
    print(IC)  
  }
  return(IC_top_genes)
}

# Calculate Jaccard Coefficient between immune-pathway and independent-component
calculate_jaccard_scores <- function(imm_path_list, IC_top_genes) {
  jaccard_coef <- function(data_A, data_B) {
    intersect_A_B <- intersect(data_A, data_B)
    union_A_B <- union(data_A, data_B)
    length(intersect_A_B) / length(union_A_B)
  }
  # Initialize matrix for results
  path_size <- length(imm_path_list)
  IC_size <- length(IC_top_genes)
  path_IC_gene <- matrix(0, nrow = path_size, ncol = IC_size, dimnames = list(names(imm_path_list), names(IC_top_genes)))
  jaccard <- matrix(0, nrow = path_size, ncol = IC_size, dimnames = list(names(imm_path_list), names(IC_top_genes)))
  # Calculate Jaccard scores and intersection genes for each pair
  for (i in 1:path_size) {
    for (j in 1:IC_size) {
      datai <- imm_path_list[[i]]
      dataj <- IC_top_genes[[j]]
      gene <- intersect(datai, dataj)
      con_gene <- paste(gene, collapse = ",")
      path_IC_gene[i, j] <- con_gene
      jaccard[i, j] <- jaccard_coef(datai, dataj)
    }
  }
  # Return list containing both Jaccard scores and intersection genes
  list(JaccardScore = jaccard, PathICGene = path_IC_gene)
}



# Calculate the centrality score to select the immune-related IC
library(igraph)
calculate_centrality_scores <- function(JaccardScore, num_iterations = num_iterations) {
  # Constructing IC-IC crosstalk networks
  edge <- as.matrix(JaccardScore)
  edget <- t(edge)
  IC_IC <- edget%*%edge
  adj_matrix <- as.matrix(IC_IC)
  # Creating graph objects
  graph <- graph.adjacency(adj_matrix, mode = "undirected", weighted = TRUE, add.rownames = TRUE)
  # Calculating PageRank
  temp <- page.rank(graph, vids = V(graph), directed = FALSE, damping = 0.90, weights = NULL)
  real_centra <- as.matrix(temp$vector)
  # Initialization matrix for storing true centrality and random centrality scores
  centrality_scores <- matrix(nrow = ncol(adj_matrix), ncol = num_iterations + 1)
  centrality_scores[,1] <- real_centra
  # Generating random networks and computing centrality scores
  for (i in 1:num_iterations) {
    perm_IC <- sample(colnames(adj_matrix), replace = FALSE)
    perm_adj_matrix <- adj_matrix[perm_IC, perm_IC]
    graph_perm <- graph.adjacency(perm_adj_matrix, mode = "undirected", weighted = TRUE, add.rownames = TRUE)
    perm_centra <- page.rank(graph_perm, vids = V(graph_perm), directed = FALSE, damping = 0.90, weights = NULL)$vector
    centrality_scores[,i + 1] <- perm_centra
  }
  # Calculating p-value
  pval <- apply(centrality_scores[, -1, drop = FALSE], 1, function(x) mean(x >= centrality_scores[, 1]))
  # Adjusting p-value
  p_adjust <- p.adjust(pval, method = "BH")
  # Integrating results
  results <- data.frame(
    IC = colnames(adj_matrix),
    P_value = pval,
    FDR = p_adjust)
  # Ranking and selecting significant IC
  significant_ICs <- results[results$P_value < 0.05, "IC"]
  return(list(centrality_scores = centrality_scores, results = results, significant_ICs = significant_ICs))
}

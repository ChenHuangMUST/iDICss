library(caret)
library(pROC)
library(Hmisc)
calculate_auc_metrics <- function(data, markers, dataset_name) {
  results <- data.frame(marker = character(), AUC = numeric(), AUC_CI = character(), 
                        F1 = numeric(), C_index = numeric(), Sensitivity = numeric(), 
                        Specificity = numeric(), Accuracy = numeric(), PPV = numeric(), 
                        NPV = numeric(), dataset = character())
  
  for (marker in markers) {
    # Calculate ROC and AUC
    roc_result <- roc(data$Response, data[[marker]], ci = TRUE, ci.alpha = 0.95)
    auc_value <- round(roc_result$auc, 3)
    auc_ci <- paste("(", round(roc_result$ci[1], 2), "-", round(roc_result$ci[2], 2), ")", sep = "")
    
    # Calculate confusion matrix and other metrics
    if(marker == "RS"){
      predictions <- ifelse(data[[marker]] < median(data[[marker]]), 1, 0) 
    }else{
      predictions <- ifelse(data[[marker]] > median(data[[marker]]), 1, 0) 
    }
    cm <- confusionMatrix(as.factor(predictions), as.factor(data$Response), positive = "1")
    
    # Extract metrics from confusion matrix
    sensitivity <- as.numeric(cm$byClass['Sensitivity'])
    specificity <- as.numeric(cm$byClass['Specificity'])
    ppv <- as.numeric(cm$byClass['Pos Pred Value'])
    npv <- as.numeric(cm$byClass['Neg Pred Value'])
    accuracy <- as.numeric(cm$overall['Accuracy'])
    
    # Compute F1 Score manually
    f1 <- 2 * (ppv * sensitivity) / (ppv + sensitivity)
    
    # C-index 
    c_index <- round(somers2(as.numeric(predictions), as.numeric(as.character(data$Response)))[1],2)
    # Append results to dataframe
    results <- rbind(results, data.frame(marker = marker, AUC = auc_value, AUC_CI = auc_ci, 
                                         F1 = f1, C_index = c_index, Sensitivity = sensitivity, 
                                         Specificity = specificity, Accuracy = accuracy, 
                                         PPV = ppv, NPV = npv, dataset = dataset_name))
  }
  return(results)
}

# Example markers and dataset name
markers <- c("RS")
roc_results <- lapply(risk_score1,function(x){
  y <- rbind(calculate_auc_metrics(x[["dataset1"]], markers,"dataset1"),
             calculate_auc_metrics(x[["dataset2"]], markers,"dataset2"),
             calculate_auc_metrics(x[["dataset3"]], markers,"dataset3"))
  return(y)
})
roc_results1 <- do.call(rbind,roc_results)
roc_results1$Model <- rownames(roc_results1)
roc_results1$Model <- gsub("\\.C\\d*$", "",roc_results1$Model)
library(reshape2)
# Melting the data for visualization
roc_results2 <- melt(roc_results1, id.vars = c("Model", "dataset"), 
                     measure.vars = c("AUC", "F1", "C_index", "Sensitivity", 
                                      "Specificity", "Accuracy", "PPV", "NPV"))
write.csv(roc_results2,"perf_MLcomb_ORR.csv",quote = F,row.names = F)
perf_iDICss_ORR <- read.csv("performance_iDICss_ORR.csv",
                            header = T,sep = ",")
colnames(perf_iDICss_ORR)[1] <- "Model"
perf_iDICss_ORR$Model <- ifelse(perf_iDICss_ORR$Model == "risk_score","iDICss",perf_iDICss_ORR$Model)
perf_iDICss_ORR <- perf_iDICss_ORR[which(perf_iDICss_ORR$Model=="iDICss"),]

perf_integrate <- rbind(perf_iDICss_ORR[which(perf_iDICss_ORR$Model=="iDICss"),],roc_results2)

dd2 <- perf_integrate[which(perf_integrate$variable=="C_index"),]
dd2 <- cbind.data.frame("ID" = dd2$dataset,
                        "Cindex" = dd2$value,
                        "Model" = dd2$Model)
###plot the C-index heatmap of all the combinations 
library(tidyr)
library(tibble)
dd2 <- pivot_wider(dd2,names_from = 'ID',values_from = 'Cindex')%>%as.data.frame()
dd2[,-1] <- apply(dd2[,-1],2,as.numeric)

dd2<- dd2 %>%
  column_to_rownames("Model")
########visualization########
library(ComplexHeatmap)
library(BART)
library(snowfall)
library(RColorBrewer)

Cindex_mat <- dd2
avg_Cindex <- apply(Cindex_mat, 1, mean)          
avg_Cindex <- sort(avg_Cindex, decreasing = T)    
Cindex_mat <- Cindex_mat[names(avg_Cindex), ]     
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) 
row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,stat = 'identity',
                                          gp = gpar(fill = "steelblue", col = NA),
                                          add_numbers = T, numbers_offset = unit(-10, "mm"),
                                          axis_param = list("labels_rot" = 0),
                                          numbers_gp = gpar(fontsize = 11, col = "white"),
                                          width = unit(3, "cm")),
                       show_annotation_name = F)

#CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") 
#names(CohortCol) <- colnames(Cindex_mat)
CohortCol <- c("dataset1" = "#1F78B4",
               "dataset2" = "#B2DF8A",
               "dataset3" = "#33A02C")

col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                          col = list("Cohort" = CohortCol),
                          show_annotation_name = F)
#"#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C"
cellwidth = 1
cellheight = 0.5
hm <- Heatmap(as.matrix(Cindex_mat), name = "C-index",
              right_annotation = row_ha, 
              top_annotation = col_ha,
              col = c("#4195C1", "#FFFFFF", "#CB5746"), 
              rect_gp = gpar(col = "black", lwd = 1), 
              cluster_columns = FALSE, cluster_rows = FALSE, 
              show_column_names = FALSE, 
              show_row_names = TRUE,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 14),
              
              width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
              height = unit(cellheight * nrow(Cindex_mat), "cm"),
              column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)), 
              column_title = NULL,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                          x, y, gp = gpar(fontsize = 11))
              }
)

pdf("machine_Cindex_ORR.pdf", width = cellwidth * ncol(Cindex_mat) + 7, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm)
invisible(dev.off())
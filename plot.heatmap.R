setwd("")
load("res.RData")
data_table2 <- res[["Cindex.res"]]

###plot the C-index heatmap of all the combinations 
library(tidyr)
dd2 <- pivot_wider(data_table2,names_from = 'ID',values_from = 'Cindex')%>%as.data.frame()
dd2[,-1] <- apply(dd2[,-1],2,as.numeric)

dd2<- dd2 %>%
  column_to_rownames("Model")
save(dd2 ,file = "dd2.RData")

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

CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") 
names(CohortCol) <- colnames(Cindex_mat)
col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                          col = list("Cohort" = CohortCol),
                          show_annotation_name = F)

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

pdf("machine_Cindex.pdf", width = cellwidth * ncol(Cindex_mat) + 7, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm)
invisible(dev.off())
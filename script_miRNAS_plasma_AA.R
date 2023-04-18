install.packages("reshape")
install.packages("ggplot2")
install.packages("plotly")
install.packages("heatmaply")

library(tidyr)
library(reshape)
library(heatmaply)


########  1. Construcción del heatmap ------------------

miRNAs_Reactome <- read_csv("~/Desktop/miRNAs_inflammasome_AA_severas_vs_leves/miRNAs_Reactome.csv")

df2.long <- pivot_longer(miRNAs_Reactome, cols = c(MiRNA_1, MiRNA_2, MiRNA_3, MiRNA_4, MiRNA_5, MiRNA_6, MiRNA_7, MiRNA_8, MiRNA_9, MiRNA_10, MiRNA_11, MiRNA_12, MiRNA_13, MiRNA_14, MiRNA_15, MiRNA_16, MiRNA_17, MiRNA_18), 
                         names_to = c("column"), 
                         values_to = "miRNA")

data <- df2.long[c(1,8,5)]
colnames(data)[3]<-"p_value"

###### 1.1 obtener la matriz para el heatmap -----------
library(reshape2)
acast(data, pathway~miRNA, value.var=p_value)
data_matrix<-xtabs(p_value~pathway+miRNA, data=data)

######## 1.2. plot the heatmap ---------------
library(circlize)
col_fun = colorRamp2(c(0, 1, 2, 3, 4), c("cornflowerblue", "yellow", "orange", "red", "purple"))
col_fun(seq(-3, 3))
data_heatmap                       <- ComplexHeatmap::Heatmap(data_matrix, 
                                            cluster_columns = TRUE, 
                                            cluster_rows = TRUE, 
                                            col=col_fun,
                                            clustering_distance_rows = "euclidean",
                                            # clustering_method_rows = "complete",
                                            show_row_names = TRUE, 
                                            show_column_names = TRUE,
                                            width = unit(20, "mm"),
                                            height = unit(180, "mm"),
                                            border_gp = gpar(col = "black", lty = 2),                                            km = 4,
                                            column_names_gp = grid::gpar(fontsize = 2),
                                            row_names_gp = grid::gpar(fontsize = 1),
                                            #row_order = 1:164,
                                            show_heatmap_legend = FALSE)

postscript("allHeatMaps.eps", 
           width = 200, 
           height = 1000)
draw(data_heatmap)
graphics.off()

###### miRNAS_to_genes heatmap
miRNAS_to_gene <- read_delim("miRNAS_to_gene.csv",  delim = ";", escape_double = FALSE, trim_ws = TRUE)
library(dplyr)
library(magrittr)
miRNAS_to_gene_table<-miRNAS_to_gene %>%
	group_by(miRNA, gene) %>%
	summarise(freq = n()) %>%
	arrange(desc(freq))
write.csv2(miRNAS_to_gene_table, file="miRNAS_to_gene_table.csv", sep=";")


library(curl)       # read file from google drive
library(gplots)     # heatmap.2
library(dendextend) # make and color dendrogram
library(colorspace) 
library(reshape2)
table<-table(miRNAS_to_gene_table$gene, miRNAS_to_gene_table$miRNA)
hclust_rows <- as.dendrogram(hclust(dist(table)))  # Calculate hclust dendrograms
hclust_cols <- as.dendrogram(hclust(dist(t(table))))
colnames (table) <- c("1","2","3","4", "5", "6")
heatmap_miRNAs_gene<-heatmap(table, scale="row",cexCol = 0.7, Rowv = hclust_rows, col= colorRampPalette(brewer.pal(8, "Blues"))(25))
cairo_ps("heatmap_miRNAs_gene.eps")
heatmap(table, scale="row",cexCol = 0.7, Rowv = hclust_rows, col= colorRampPalette(brewer.pal(8, "Blues"))(25))
dev.off()
postscript("heatmap_miRNAs_gene", 
           width = 200, 
           height = 1000)
heatmap(table, scale="row",cexCol = 0.7, Rowv = hclust_rows, col= colorRampPalette(brewer.pal(8, "Blues"))(25))
graphics.off()

tiff("heatmap_miRNAs_gene.tiff", width = 7, height = 7, units = 'in', res = 300, compression = 'lzw')
heatmap(table, cexCol = 0.7, cexRow = 0.05,Colv=hclust_cols,Rowv = hclust_rows, col= colorRampPalette(brewer.pal(8, "Blues"))(25))
dev.off()



# dend1 <- as.dendrogram(hclust(dist(table)))
# c_group <- 8
# dend1 <- color_branches(dend1, k = c_group, col = rainbow_hcl) # add color to the lines
# dend1 <- color_labels(dend1, k = c_group, col = rainbow_hcl)   # add color to the labels
# col_labels <- get_leaves_branches_col(dend1)
# col_labels <- col_labels[order(order.dendrogram(dend1))]
# dend1 <- set(dend1, "labels_cex", 0.8)
# par(mar = c(1,1,1,14))
# plot_horiz.dendrogram(dend1, side = F) # use side = T to horiz mirror if needed
# 
# my_palette <- colorRampPalette(c("lightgrey", "black", "forestgreen"))(n = 1000)
# par(cex.main=0.3)                   # adjust font size of titles
# heatmap.2(table, main = 'Fasting to Lose Weight',
#           # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
#                                      # order by branch mean so the deepest color is at the top
#           dendrogram = "row",        # no dendrogram for columns
#           Rowv = dend1,              # * use self-made dendrogram
#           Colv = "NA",               # make sure the columns follow data's order
#           col = my_palette,         # color pattern of the heatmap
#           
#           trace="none",              # hide trace
#           density.info="none",       # hide histogram
#           
#           margins = c(5,18),         # margin on top(bottom) and left(right) side.
#           cexRow=0.5, cexCol = 0.8,      # size of row / column labels
#           xlab = "Year",
#           srtRow=0, adjRow = c(0.1,1),
#           srtCol=90, adjCol = c(1,0.1),# adjust the direction of row label to be horizontal
#           # margin for the color key
#           # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
#           key.par=list(mar=c(5,1,3,1)),
#           RowSideColors = col_labels, # to add nice colored strips        
#           colRow = col_labels         # add color to label
#           )

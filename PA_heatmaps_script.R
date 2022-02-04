#getting the libraries that I need
#install.packages("ape")
#install.packages("RColorBrewer")
#install.packages("colorRamps")
#install.packages("textshape")
library(ape)
library(RColorBrewer)
library(colorRamps)
library(gplots)
library(heatmap3)
library(textshape)


#set working directory
setwd(getwd())

#args are  TreeFile, PA_File, ColorsFile, imageName
args = commandArgs(trailingOnly = TRUE)
if (length(args)!=4){
  stop("args are  TreeFile, PA_File, ColorsFile, imageName")
}

#getting the tree into R with ape
core_tree <- read.tree(args[1])
#core_tree <- read.tree('inputs/core_gene_alignment.aln.treefile')

core_tree_r <- root(core_tree, outgroup = "ehirae", resolve.root = T)
core_tree_r$edge.length[which(core_tree_r$edge.length == 0)] <- 0.00001

core_tree_um <- chronopl(core_tree_r,
                         lambda = 0.1,
                         tol = 0)


core_tree_d <- as.dendrogram(as.hclust.phylo(core_tree_um))
plot(core_tree_d)

#getting PA table and tree to have the correct labels
#getting label names sorted properly
combined_matrix  <- read.csv2(args[2], sep = ',', na.strings="" )

name_order <- order.dendrogram(core_tree_d)
name_name <- labels(core_tree_d)
name_position <- data.frame(name_name, name_order)
name_position <- name_position[order(name_position$name_order),]
new_order <- match(name_position$name_name, combined_matrix$Genome_ID)
combined_ordered_matrix <- combined_matrix[new_order,]
names_for_rows <- combined_ordered_matrix$Isolate_Pathname

#making the heatmap
color <- colorRampPalette(c('#35193e'))(1)

print(length(combined_ordered_matrix))
combined_order_matrix_numerical <- apply(as.matrix(combined_ordered_matrix[, 2:length(combined_ordered_matrix)]), 2, as.numeric)
combined_order_matrix_numerical[combined_order_matrix_numerical==""]<-NA
combined_order_matrix_numerical[is.na(combined_order_matrix_numerical)] <- 0

#row of colours for habitat
##write.csv(names_for_rows,"/home/alex/Documents/arete/niche_data/haley_trees/AMR_tree_heatmap/row_names_efaecium_tree_heatmap_alexrun.csv", row.names = FALSE)
row_col_table <- read.csv(args[3], as.is = T)
row_col <- row_col_table$Env_Color
row_colors <- cbind(Habitat=c(row_col_table$Env_Color),Geography=c(row_col_table$Geo_Color),Clade=c(row_col_table$clade_col))

#heirarchal clustering
coltrans = t(combined_order_matrix_numerical)
colDistance = dist(coltrans, method = "manhattan")
colDistance2 = t(colDistance)
colCluster = hclust(colDistance, method = "complete")
colDend = as.dendrogram(colCluster)
colDend = reorder(colDend, colMeans(coltrans))
combined_order_matrix_numerical[combined_order_matrix_numerical == 0] <-NA

#saving heatmap
pdf(NULL)
png(args[4], width = 1000, height = 1000)
#labCol=NA,
heatmap3(as.matrix(combined_order_matrix_numerical),Colv = colDend, Rowv = as.dendrogram(core_tree_d), col = color, margins = c(20, 20), RowSideColors = row_colors, cexCol=2, legendfun=function()showLegend(legend=c(""), col=c("#FFFFFF")), labRow=NA,  scale = 'none')
legend(x="right", legend=c("Agriculture", "Clinical", "Wastewater Mun.", "Wastewater Agr.", "Natural Water", "", "Alberta,Canada", "United Kingdom", "", "Clade A", "Clade B"),
       fill=c("#33a02c", "#e31a1c", "#cab2d6", "#1f78b4", "#a6cee3", "white", "#e36951", "#9bbdff", "white", "#841e5a", "#f6c19f"), border=FALSE)
dev.off()

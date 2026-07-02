## 1. Quality control

# LIBRARIES
library(Seurat)
library(R.utils)
library(dplyr)
library(Matrix)
library(ggplot2)
library(scDblFinder)
library(patchwork)

###UPLOAD OBJECTS

patients <- c("", "", "", "")
dir_zip <- "/zipped/"
dir_data <- "/RNA/"
dir_muts <- "/MUTS/"
dir_plots <- "/"
dir_objects <- "/"

# unzip objects
unzip_data <- function(patients, input_dir, output_dir) {
  files <- c("barcodes.tsv", "features.tsv", "matrix.mtx")
  for (i in patients) {
    dir.create(paste0(output_dir, i), recursive = TRUE, showWarnings = FALSE)
    for (j in files) { 
      file <- paste0(input_dir, i, "/", j, ".gz")
      dest <- paste0(output_dir, i, "/", j)
      R.utils::gunzip(file, dest, remove = FALSE)
    }
  }
  
}
unzip_data(patients, dir_zip, dir_data)

# upload into seurat objects
create_seurat <- function(patients, dir_data) {
  for (i in patients) {
    dir <- paste0(dir_data, i)
    mat <- Read10X(dir, strip.suffix = TRUE)
    name <- paste0("s_", sub("-", "_", i))
    obj <- CreateSeuratObject(mat)
    assign(name, obj, envir = .GlobalEnv)
  }
}

create_seurat(patients, dir_data)

objects <- c("", "", "", "")

## ADD MUTATIONS

# upload data
mut_data_1 <- read.delim(paste0(dir_muts, "1_filtered_annot_matrix.tsv"), header = TRUE, sep = "\t")
mut_data_2 <- read.delim(paste0(dir_muts, "2_filtered_annot_matrix.tsv"), header = TRUE, sep = "\t")
mut_data_3 <- read.delim(paste0(dir_muts, "3_filtered_annot_matrix.tsv"), header = TRUE, sep = "\t")
mut_data_3 <- read.delim(paste0(dir_muts, "4_filtered_annot_matrix.tsv"), header = TRUE, sep = "\t")

# put the cells barcode as rownames
rownames(mut_data_1) <- mut_data_1$CellBarcode
mut_data_1$CellBarcode <- NULL

rownames(mut_data_2) <- mut_data_2$CellBarcode
mut_data_2$CellBarcode <- NULL

rownames(mut_data_3) <- mut_data_3$CellBarcode
mut_data_3$CellBarcode <- NULL

rownames(mut_data_3) <- mut_data_3$CellBarcode
mut_data_3$CellBarcode <- NULL

# filter the mutations with at least 5% cells mutated
results_1 <- list()
for (i in 1:ncol(mut_data_1)){
  counts <- table(mut_data_1[[i]])
  perc_genotyped <- (sum(counts[c("Wt", "Het", "Hom")], na.rm = TRUE)/(sum(counts[c("Wt", "Het", "Hom", "NotAvailable")], na.rm = TRUE))) *100
  perc_mutated <- (sum(counts[c("Het", "Hom")],na.rm = TRUE)/(sum(counts[c("Wt", "Het", "Hom")],na.rm = TRUE))) *100
  results_1[[i]] <- data.frame(
    mutation = colnames(mut_data_1)[i],
    perc_genotyped = perc_genotyped,
    perc_mutated = perc_mutated)
}
perc_mut_1 <- do.call(rbind, results_1)

results_2 <- list()
for (i in 1:ncol(mut_data_2)){
  counts <- table(mut_data_2[[i]])
  perc_genotyped <- (sum(counts[c("Wt", "Het", "Hom")], na.rm = TRUE)/(sum(counts[c("Wt", "Het", "Hom", "NotAvailable")], na.rm = TRUE))) *100
  perc_mutated <- (sum(counts[c("Het", "Hom")],na.rm = TRUE)/(sum(counts[c("Wt", "Het", "Hom")],na.rm = TRUE))) *100
  results_2[[i]] <- data.frame(
    mutation = colnames(mut_data_2)[i],
    perc_genotyped = perc_genotyped,
    perc_mutated = perc_mutated)
}
perc_mut_2 <- do.call(rbind, results_2)

results_3 <- list()
for (i in 1:ncol(mut_data_3)){
  counts <- table(mut_data_3[[i]])
  perc_genotyped <- (sum(counts[c("Wt", "Het", "Hom")], na.rm = TRUE)/(sum(counts[c("Wt", "Het", "Hom", "NotAvailable")], na.rm = TRUE))) *100
  perc_mutated <- (sum(counts[c("Het", "Hom")],na.rm = TRUE)/(sum(counts[c("Wt", "Het", "Hom")],na.rm = TRUE))) *100
  results_3[[i]] <- data.frame(
    mutation = colnames(mut_data_3)[i],
    perc_genotyped = perc_genotyped,
    perc_mutated = perc_mutated)
}
perc_mut_3 <- do.call(rbind, results_3)


results_4 <- list()
for (i in 1:ncol(mut_data_3)){
  counts <- table(mut_data_3[[i]])
  perc_genotyped <- (sum(counts[c("Wt", "Het", "Hom")], na.rm = TRUE)/(sum(counts[c("Wt", "Het", "Hom", "NotAvailable")], na.rm = TRUE))) *100
  perc_mutated <- (sum(counts[c("Het", "Hom")],na.rm = TRUE)/(sum(counts[c("Wt", "Het", "Hom")],na.rm = TRUE))) *100
  results_4[[i]] <- data.frame(
    mutation = colnames(mut_data_3)[i],
    perc_genotyped = perc_genotyped,
    perc_mutated = perc_mutated)
}
perc_mut_4 <- do.call(rbind, results_4)

# mutations identified by NGS in each patient
mut_1 <- c("BTK_c.1442G.C_p.Cys481Ser", "PLCG2_c.2126A.T_p.Tyr709Phe", "PLCG2_c.2120C.T_p.Ser707Phe")
mut_2 <- c("BTK_c.1442G.C_p.Cys481Ser", "BTK_c.1441T.A_p.Cys481Ser")
mut_3 <- c("BTK_c.1441T.A_p.Cys481Ser")
mut_4 <- c("BTK_c.1442G.C_p.Cys481Ser")

# table with only identified mutations info
sig_mut_1 <- mut_data_1[,colnames(mut_data_1) %in% mut_1, drop = FALSE]
sig_mut_2 <- mut_data_2[,colnames(mut_data_2) %in% mut_2, drop = FALSE]
sig_mut_3 <- mut_data_3[,colnames(mut_data_3) %in% mut_3, drop = FALSE]
sig_mut_4 <- mut_data_3[,colnames(mut_data_3) %in% mut_4, drop = FALSE]

# check if all cells in the seurat object have info in the mutations table before merging
setdiff(colnames(s_1), rownames(sig_mut_1))
# GTGCAAGGT_AACGCTAGT_AACAAGTGG
sig_mut_1["GTGCAAGGT_AACGCTAGT_AACAAGTGG", ] <- "NotAvailable"

setdiff(colnames(s_2), rownames(sig_mut_2))
# ATGACAGCA_AACGCTAGT_AACAAGTGG
sig_mut_2["ATGACAGCA_AACGCTAGT_AACAAGTGG", ] <- "NotAvailable"

setdiff(colnames(s_3), rownames(sig_mut_3))
# AGGACTCAC_AACGCTAGT_AACAAGTGG
sig_mut_3["AGGACTCAC_AACGCTAGT_AACAAGTGG", ] <- "NotAvailable"

setdiff(colnames(s_4), rownames(sig_mut_4))
# GGATAGATG_AACGTCCAA_AACAAGTGG
sig_mut_4["GGATAGATG_AACGTCCAA_AACAAGTGG", ] <- "NotAvailable"


# match significant mutations with cell_ids in seurat object
sig_mut_1 <- sig_mut_1[colnames(s_1), , drop = FALSE]
sig_mut_2 <- sig_mut_2[colnames(s_2), , drop = FALSE]
sig_mut_3 <- sig_mut_3[colnames(s_3), , drop = FALSE]
sig_mut_4 <- sig_mut_4[colnames(s_4), , drop = FALSE]

s_1 <- AddMetaData(s_1, metadata = sig_mut_1)
s_2 <- AddMetaData(s_2, metadata = sig_mut_2)
s_3 <- AddMetaData(s_3, metadata = sig_mut_3)
s_4 <- AddMetaData(s_4, metadata = sig_mut_4)


# filter heterozygous cells for BTK mutation (male patients)
table(s_3$BTK_c.1441T.A_p.Cys481Ser)
s_3 <- subset(s_3, subset = BTK_c.1441T.A_p.Cys481Ser != "Het")

table(s_4$BTK_c.1442G.C_p.Cys481Ser)
s_4 <- subset(s_4, subset =  BTK_c.1442G.C_p.Cys481Ser != "Het")

table(s_2$BTK_c.1442G.C_p.Cys481Ser, s_2$BTK_c.1441T.A_p.Cys481Ser)
s_2 <- subset(s_2, subset =  BTK_c.1442G.C_p.Cys481Ser != "Het")
s_2 <- subset(s_2, subset =  BTK_c.1441T.A_p.Cys481Ser != "Het")


## change the active identity + add % of mt genes --> returns a list with all the objects
add_percent.mt <- function(patients) {
  obj <- list() 
  for (i in patients) {
    pat <- gsub("-", "_", i)
    name <- paste0("s_", pat)
    sample <- get(name)
    sample@meta.data$patient <- factor(rep(pat, nrow(sample@meta.data)))
    Idents(sample) <- sample@meta.data$patient
    sample <- PercentageFeatureSet(sample, pattern = "^MT-", col.name = "percent.mt")
    obj[[pat]] <- sample
  }
  assign("obj_list", obj, envir = .GlobalEnv)
}

add_percent.mt(patients)

plot_metrics <- function(obj_list) {
  plots <- list()
  
  for (name in names(obj_list)) {
    obj <- obj_list[[name]]
    
    par(mfrow = c(3, 1))
    
    plot(density(obj@meta.data$nCount_RNA),
         main = paste(name, "- nCount_RNA"))
    
    plot(density(obj@meta.data$nFeature_RNA),
         main = paste(name, "- nFeature_RNA"))
    
    plot(density(obj@meta.data$percent.mt),
         main = paste(name, "- percent.mt"))
    plots[[name]] <- recordPlot()
  }
  return(plots)
}

metrics_prefilt <- plot_metrics(obj_list)

save_replay_plots <- function(plot_list, dir_plots, filename) {
  if (!dir.exists(dir_plots)) {
    dir.create(dir_plots, recursive = TRUE)
  }
  pdf(file.path(dir_plots, filename))
  for (name in names(plot_list)) {
    replayPlot(plot_list[[name]])
  }
  dev.off()
}

save_replay_plots(metrics_prefilt, dir_plots, "metrics_prefilter.pdf")

vln_seurat <- function(obj_list, feature) {
  meta_obj <- lapply(names(obj_list), function(i) {
    df <- obj_list[[i]]@meta.data
    df$sample <- i
    return(df)
  })
  meta_obj <- bind_rows(meta_obj)
  
  ggplot(meta_obj, aes(x = sample, y = .data[[feature]], fill = sample)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(x = NULL, y = feature)
}

plot_counts_obj   <- vln_seurat(obj_list, "nCount_RNA")
plot_features_obj <- vln_seurat(obj_list, "nFeature_RNA")
plot_mt_obj       <- vln_seurat(obj_list, "percent.mt")

plot_counts_obj + plot_features_obj + plot_mt_obj
ggsave(paste0(dir_plots,"count_feat_mt_prefilter.png"), plot_counts_obj + plot_features_obj + plot_mt_obj, width = 12, height = 5)



### QUALITY CONTROL
## hard filtering
hard_filter <- list()
for (i in names(obj_list)) {
  sample <- obj_list[[i]]
  
  # hard filtering (1%)
  sample <- subset(sample, subset = 
                     nFeature_RNA > quantile(nFeature_RNA, 0.01) &
                     nFeature_RNA < quantile(nFeature_RNA, 0.99) &
                     nCount_RNA > quantile(nCount_RNA, 0.01) &
                     nCount_RNA < quantile(nCount_RNA, 0.99) &
                     percent.mt > quantile(percent.mt, 0.01) &
                     percent.mt < quantile(percent.mt, 0.99)
  )
  
  # remove extreme values
  sample <- subset(sample, subset = 
                     nFeature_RNA > 50 & 
                     nCount_RNA > 200 & 
                     percent.mt < 60
  )
  
  hard_filter[[i]] <- sample
}

metrics_hardfilter <- plot_metrics(hard_filter)
metrics_hardfilter
save_replay_plots(metrics_hardfilter, dir_plots, "metrics_hard_filter.pdf")


## find doublets --> scDblFinder
db_finder <- list()
for (i in names(hard_filter)) {
  sample <- hard_filter[[i]]
  set.seed(1991)
  doublets <- scDblFinder(as.SingleCellExperiment(sample), clusters=FALSE)
  table(rownames(colData(doublets)) == rownames(sample@meta.data))
  cols_to_add <- grep("scDblFinder", colnames(colData(doublets)), value = TRUE)
  sample@meta.data <- cbind(sample@meta.data, colData(doublets)[match(rownames(sample@meta.data), rownames(colData(doublets))), cols_to_add])
  # en el match x[match(y,x),]
  db_finder[[i]] <- sample
}

par(mfrow = c(2,1))
db_scores_plots <- list()
for (i in names(db_finder)) {
  obj <- db_finder[[i]]
  plot(density(obj@meta.data$scDblFinder.score), main = paste0(i, " - scDblFinder.score"))
  db_scores_plots[[i]] <- recordPlot()
}

save_replay_plots(db_scores_plots, dir_plots, "scDblFinder_scores.pdf")

# table doublets
table_db <- data.frame()
for (i in names(db_finder)) {
  sample <- db_finder[[i]]
  singlet = sum(sample@meta.data$scDblFinder.class == "singlet")
  doublet = sum(sample@meta.data$scDblFinder.class == "doublet")
  percentage = (doublet/(singlet+doublet))*100
  df_temp <- data.frame(
    sample = i,
    singlet = singlet,
    doublet = doublet,
    perc_doublets = percentage
  )
  
  table_db <- rbind(table_db, df_temp)
}

table_db


##FILTERING
# define filter thresholds 
## for percent.mt
filters_mt <- list()
for (i in names(obj_list)) {
  sample <- obj_list[[i]]
  max_mt <- median(log10(sample@meta.data$percent.mt)) +  2*mad(log10(sample@meta.data$percent.mt), constant = 1)
  min_mt <- median(log10(sample@meta.data$percent.mt)) -  2*mad(log10(sample@meta.data$percent.mt), constant = 1)
  filter_mt <- c(10^min_mt, 10^max_mt)
  filters_mt[[i]] <- filter_mt
}

## for counts --> 
filters_counts <- list()
medians_counts <- list()
for (i in names(obj_list)) {
  sample <- obj_list[[i]]
  median <- 10^(median(log10(sample@meta.data$nCount_RNA)))
  max_counts <- median(log10(sample@meta.data$nCount_RNA)) +  2*mad(log10(sample@meta.data$nCount_RNA), constant = 1)
  min_counts <- median(log10(sample@meta.data$nCount_RNA)) -  2*mad(log10(sample@meta.data$nCount_RNA), constant = 1)
  filter_counts <- c(10^min_counts, 10^max_counts)
  filters_counts[[i]] <- filter_counts
  medians_counts[[i]] <- median
}

## for doublets --> scDblFinder.class == "singlet"

## plot the counts thresholds 
for (i in names(db_finder)) {
  
  obj <- db_finder[[i]]
  d <- density(obj@meta.data$nCount_RNA)
  
  plot(d, main = paste0(i, " - nCount_RNA"))
  
  abline(v = filters_counts[[i]], col = "red", lwd = 2, lty = 2)
  abline(v = medians_counts[[i]], col = "blue", lwd = 2, lty = 2)
}
## plot the percent.mt thresholds 
for (i in names(db_finder)) {
  
  obj <- db_finder[[i]]
  d <- density(obj@meta.data$percent.mt)
  
  plot(d, main = paste0(i, " - percent.mt"))
  
  abline(v = filters_mt[[i]], col = "red", lwd = 2, lty = 2)
}


# Filtering
filtered <- list()
for (i in names(db_finder)) {
  sample <- db_finder[[i]]
  # filtering
  sample <- subset(sample, subset = 
                     percent.mt < filters_mt[[i]][[2]] &
                     scDblFinder.class == "singlet"&
                     nCount_RNA > filters_counts[[i]][[1]] & 
                     nCount_RNA < filters_counts[[i]][[2]])
  filtered[[i]] <- sample
}

metrics_filtered <- plot_metrics(filtered)

save_replay_plots(metrics_filtered, dir_plots, "metrics_post_filtering.pdf")

par(mfrow = c(2,1))
db_scores_plots2 <- list()
for (i in names(filtered)) {
  obj <- filtered[[i]]
  plot(density(obj@meta.data$scDblFinder.score), main = paste0(i, " - scDblFinder.score"))
  db_scores_plots2[[i]] <- recordPlot()
}
save_replay_plots(db_scores_plots2, dir_plots, "scDblFinder_scores_post_filter.pdf")

f_1 <- filtered[[1]]

f_1$mutated_BTK <- ifelse(f_1$BTK_c.1442G.C_p.Cys481Ser == "Hom","YES",
                                 ifelse(f_1$BTK_c.1442G.C_p.Cys481Ser == "NotAvailable","NA","NO"))

f_1$mutated_PLCG2 <- ifelse(f_1$PLCG2_c.2120C.T_p.Ser707Phe %in% c("Hom", "Het") 
                                   | f_1$PLCG2_c.2126A.T_p.Tyr709Phe %in% c("Hom", "Het"),"YES",
                                   ifelse(f_1$PLCG2_c.2120C.T_p.Ser707Phe == "NotAvailable"
                                          | f_1$PLCG2_c.2126A.T_p.Tyr709Phe == "NotAvailable","NA","NO"))

f_1$mutation <- ifelse(f_1$mutated_BTK == "NO" & f_1$mutated_PLCG2 != "YES", "WT",
                              ifelse(f_1$mutated_PLCG2 == "YES", "PLCG2",
                                     ifelse(f_1$mutated_BTK == "YES", "BTK", "NA")))

f_1$mutated <- ifelse(f_1$mutation %in% c("BTK", "PLCG2") , "YES",
                             ifelse(f_1$mutation == "WT", "NO", "NA"))
## heatmap co-ocurrence of mutations
mat <- table(
  f_1$mutated_BTK,
  f_1$mutated_PLCG2)
df <- as.data.frame(mat)
colnames(df) <- c("BTK", "PLCG2", "Count")
df <- df %>%
  group_by(BTK) %>%
  mutate(prop = Count / sum(Count))

ggplot(df, aes(x = PLCG2, y = BTK, fill = prop)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = Count), size = 5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#A6BDDB") +
  labs(
    x = "PLCG2 mutation status",
    y = "BTK mutation status",
    title = "BTK/PLCG2 mutation co-occurrence (Patient 1)",
    fill = "Proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


## 2. INTEGRATION 

# LIBRARIES
library(Seurat)
library(R.utils)
library(dplyr)
library(Matrix)
library(ggplot2)
library(scDblFinder)
library(patchwork)
library(clustree)

p_1 <- filtered[[1]]
p_2 <- filtered[[2]]
p_3 <- filtered[[3]]
p_4 <- filtered[[4]]

rm(filtered)

# data normalization
p_1 <- NormalizeData(p_1)
p_1 <- FindVariableFeatures(p_1)
p_1 <- ScaleData(p_1)
p_1 <- RunPCA(p_1)
ElbowPlot(p_1, ndims = 50)
p_1 <- FindNeighbors(p_1, dims = 1:30)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (r in resolutions) {
  p_1 <- FindClusters(p_1, resolution = r)
}
clustree_uninteg <- clustree(p_1@meta.data, prefix = "RNA_snn_res.")
clustree_uninteg
p_1 <- RunUMAP(p_1, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(p_1, reduction = "umap.unintegrated")

# data normalization
p_2 <- NormalizeData(p_2)
p_2 <- FindVariableFeatures(p_2)
p_2 <- ScaleData(p_2)
p_2 <- RunPCA(p_2)
ElbowPlot(p_2, ndims = 50)
p_2 <- FindNeighbors(p_2, dims = 1:30)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (r in resolutions) {
  p_2 <- FindClusters(p_2, resolution = r)
}
clustree_uninteg <- clustree(p_2@meta.data, prefix = "RNA_snn_res.")
clustree_uninteg
p_2 <- RunUMAP(p_2, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(p_2, reduction = "umap.unintegrated")


# data normalization
p_3 <- NormalizeData(p_3)
p_3 <- FindVariableFeatures(p_3)
p_3 <- ScaleData(p_3)
p_3 <- RunPCA(p_3)
ElbowPlot(p_3, ndims = 50)
p_3 <- FindNeighbors(p_3, dims = 1:30)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (r in resolutions) {
  p_3 <- FindClusters(p_3, resolution = r)
}
clustree_uninteg <- clustree(p_3@meta.data, prefix = "RNA_snn_res.")
clustree_uninteg
p_3 <- RunUMAP(p_3, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(p_3, reduction = "umap.unintegrated")

# data normalization
p_4 <- NormalizeData(p_4)
p_4 <- FindVariableFeatures(p_4)
p_4 <- ScaleData(p_4)
p_4 <- RunPCA(p_4)
ElbowPlot(p_4, ndims = 50)
p_4 <- FindNeighbors(p_4, dims = 1:30)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (r in resolutions) {
  p_4 <- FindClusters(p_4, resolution = r)
}
clustree_uninteg <- clustree(p_4@meta.data, prefix = "RNA_snn_res.")
clustree_uninteg
p_4 <- RunUMAP(p_4, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(p_4, reduction = "umap.unintegrated")


## CELLTYPIST
for (i in patients) {
  cells <- get(paste0("p_", gsub("-", "_", i)))
  
  # Define path to celltypist output
  pathToCelltypist <- "/"
  
  # Specify celltypist annotation DIR
  dir_celltypist <- "/"
  
  library(Matrix)
  library(data.table)
  
  rna_assay <- cells[["RNA"]]
  
  # Agafem totes les layers de counts
  count_layers <- grep("^counts", Layers(rna_assay), value = TRUE)
  
  # Extraiem cada layer i assegurem que totes tinguin els mateixos gens
  counts_list <- lapply(count_layers, function(layer) {
    GetAssayData(rna_assay, layer = layer)
  })
  
  # Obtenim tots els gens presents en qualsevol layer
  all_genes <- unique(unlist(lapply(counts_list, rownames)))
  
  # Reindexem cada matriu per tenir tots els gens (afegint zeros on no existeixen)
  counts_list_aligned <- lapply(counts_list, function(mat) {
    # Creem matriu buida amb tots els gens
    missing_genes <- setdiff(all_genes, rownames(mat))
    if(length(missing_genes) > 0) {
      # Matriu esparsa de zeros per als gens que falten
      mat_extra <- Matrix(0, nrow = length(missing_genes), ncol = ncol(mat),
                          dimnames = list(missing_genes, colnames(mat)))
      # Combinem amb la matriu original
      mat <- rbind(mat, mat_extra)
    }
    # Ordenem les files com all_genes
    mat[all_genes, , drop = FALSE]
  })
  
  # Ara podem combinar-les en una sola matriu
  all_counts <- do.call(cbind, counts_list_aligned) # 44915 (genes) 65169 (cells)
  all_counts_t <- t(all_counts)
  dim(all_counts_t)
  
  # Escriure la matriu esparsa en format Matrix Market
  writeMM(all_counts_t, paste0(pathToCelltypist, "/celltypist_input.mtx"))
  
  # Escriure els noms de gens i cèl·lules
  write.table(rownames(all_counts),
              paste0(pathToCelltypist, "/genes.txt"),
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
  
  write.table(colnames(all_counts),
              paste0(pathToCelltypist, "/cells.txt"),
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
  
  rm(rna_assay)
  rm(all_counts)
  rm(all_counts_t)
  
  # Run CellTypist
  system(command = paste0(
    "celltypist --indata ", pathToCelltypist, "/celltypist_input.mtx",
    " --gene-file ", pathToCelltypist, "/genes.txt",
    " --cell-file ", pathToCelltypist, "/cells.txt",
    " --model Immune_All_High.pkl --mode best_match",
    " --outdir ", pathToCelltypist, "/model_high"
  ),intern = TRUE)
  
  system(command = paste0(
    "celltypist --indata ", pathToCelltypist, "/celltypist_input.mtx",
    " --gene-file ", pathToCelltypist, "/genes.txt",
    " --cell-file ", pathToCelltypist, "/cells.txt",
    " --model Immune_All_Low.pkl --mode best_match",
    " --outdir ", pathToCelltypist, "/model_low"
  ),intern = TRUE)
  
  # FOR LOW MODEL - prob 90
  library(readxl)
  library(stringr)
  # Load celltypist dictionary
  wiki <- read_xlsx(file.path(dir_celltypist, "encyclopedia_table.xlsx"))
  
  # Define b-cell, t-cell lineages
  blin <- wiki$`Low-hierarchy cell types`[wiki$`High-hierarchy cell types` %in% c("B cells","B-cell lineage", "Plasma cells") | grepl("Cycling B cells", wiki$`Low-hierarchy cell types`)]
  tlin <- wiki$`Low-hierarchy cell types`[wiki$`High-hierarchy cell types` %in% c("T cells","Double-negative thymocytes", "Double-positive thymocytes", "ETP") | grepl("Cycling gamma-delta T cells|Cycling T cells", wiki$`Low-hierarchy cell types`)]
  
  # Get predicted max label
  celltypist_low <- fread(file.path(pathToCelltypist, "/model_low/predicted_labels.csv"), data.table = F, header = T)
  
  # Define annotations by custom probability (prob = 0.9 by default)
  prob <- 0.9
  
  df <- fread(file.path(pathToCelltypist, "/model_low/probability_matrix.csv"), data.table = F)
  df <- data.frame(
    ID = df$V1,
    celltypist_low_prob90 = apply(df[,-1],1,function(y){
      paste(colnames(df)[-1][y > prob],collapse="|")
    })
  )
  df$celltypist_low_prob90[df$celltypist_low_prob90 == ""] <- "Unassigned"
  
  # Split each celltype label
  splits <- str_split(df$celltypist_low_prob90, pattern = "\\|")
  names(splits) <- df$ID
  
  # Define mixed cells
  mixing <- sapply(splits, function(y){
    if(length(y) <= 1){
      "No_mixed"
    }else if(length(y) > 1){
      if(sum(y %in% blin) == length(y)){
        "Mixed_B"
      }else if(sum(y %in% tlin) == length(y)){
        "Mixed_T"
      }else{
        "Mixed_other"
      }
    }
  })
  
  # Define Bcell or Tcell lineage
  lineage <- case_when(
    df$celltypist_low_prob90 %in% blin | mixing == "Mixed_B" ~ "B_lineage",
    df$celltypist_low_prob90 %in% tlin | mixing == "Mixed_T" ~ "T_lineage",
    TRUE ~ "Other"
  )
  
  # Incorporate mixing and lineage
  df$celltypist_low_prob90_mixing <- mixing
  df$celltypist_low_prob90_lineage <- lineage
  
  # Specify unassigned cells
  df$celltypist_low_prob90_mixing[df$celltypist_low_prob90 == "Unassigned"] <- "Unassigned"
  df$celltypist_low_prob90_lineage[df$celltypist_low_prob90 == "Unassigned"] <- "Unassigned"
  
  # Merge all celltypist info (max, 0.9 prob, and all probabilities)
  colnames(celltypist_low)[2] <- "celltypist_max_low"
  probs <- fread(file.path(pathToCelltypist, "/model_low/probability_matrix.csv"), data.table = F)
  colnames(probs)[2:ncol(probs)] <- paste0("celltypist_", gsub(" ", ".", colnames(probs)[2:ncol(probs)]))
  celltypist_low$celltypist_max_prob_low <- apply(probs[,2:ncol(probs)], 1, max)
  celltypist_low <- merge.data.frame(celltypist_low, df, by.x = "V1", by.y = "ID", all = T)
  celltypist_low <- merge.data.frame(celltypist_low, probs, by = "V1", all = T)
  cells@meta.data <- cbind(cells@meta.data, celltypist_low[match(rownames(cells@meta.data), celltypist_low$V1), 2:6])
  
  
  # FOR HIGH MODEL - prob 90
  # Get predicted max label
  celltypist_high <- fread(file.path(pathToCelltypist, "/model_high/predicted_labels.csv"), data.table = F, header = T)
  
  # Define annotations by custom probability (prob = 0.9 by default)
  prob <- 0.9
  
  df <- fread(file.path(pathToCelltypist, "/model_high/probability_matrix.csv"), data.table = F)
  df <- data.frame(
    ID = df$V1,
    celltypist_high_prob90 = apply(df[,-1],1,function(y){
      paste(colnames(df)[-1][y > prob],collapse="|")
    })
  )
  df$celltypist_high_prob90[df$celltypist_high_prob90 == ""] <- "Unassigned"
  
  # Split each celltype label
  splits <- str_split(df$celltypist_high_prob90, pattern = "\\|")
  names(splits) <- df$ID
  
  # Define mixed cells
  mixing <- sapply(splits, function(y){
    if(length(y) <= 1){
      "No_mixed"
    }else if(length(y) > 1){
      if(sum(y %in% blin) == length(y)){
        "Mixed_B"
      }else if(sum(y %in% tlin) == length(y)){
        "Mixed_T"
      }else{
        "Mixed_other"
      }
    }
  })
  
  # Define Bcell or Tcell lineage
  lineage <- case_when(
    df$celltypist_high_prob90 %in% blin | mixing == "Mixed_B" ~ "B_lineage",
    df$celltypist_high_prob90 %in% tlin | mixing == "Mixed_T" ~ "T_lineage",
    TRUE ~ "Other"
  )
  
  # Incorporate mixing and lineage
  df$celltypist_high_prob90_mixing <- mixing
  df$celltypist_high_prob90_lineage <- lineage
  
  # Specify unassigned cells
  df$celltypist_high_prob90_mixing[df$celltypist_high_prob90 == "Unassigned"] <- "Unassigned"
  df$celltypist_high_prob90_lineage[df$celltypist_high_prob90 == "Unassigned"] <- "Unassigned"
  
  # Merge all celltypist info (max, 0.9 prob, and all probabilities)
  colnames(celltypist_high)[2] <- "celltypist_max_high"
  probs <- fread(file.path(pathToCelltypist, "/model_high/probability_matrix.csv"), data.table = F)
  colnames(probs)[2:ncol(probs)] <- paste0("celltypist_", gsub(" ", ".", colnames(probs)[2:ncol(probs)]))
  celltypist_high$celltypist_max_prob_high <- apply(probs[,2:ncol(probs)], 1, max)
  celltypist_high <- merge.data.frame(celltypist_high, df, by.x = "V1", by.y = "ID", all = T)
  celltypist_high <- merge.data.frame(celltypist_high, probs, by = "V1", all = T)
  cells@meta.data <- cbind(cells@meta.data, celltypist_high[match(rownames(cells@meta.data), celltypist_high$V1), 2:6])
  
  assign(paste0("p_", gsub("-", "_", i)), cells, envir = .GlobalEnv)
}


# plots 
DimPlot(p_1, reduction = "umap.unintegrated", group.by = "celltypist_max_high")
DimPlot(p_2, reduction = "umap.unintegrated", group.by = "celltypist_max_high")
DimPlot(p_3, reduction = "umap.unintegrated", group.by = "celltypist_max_high")
DimPlot(p_4, reduction = "umap.unintegrated", group.by = "celltypist_max_high")


## FILTER B-CELLS
DimPlot(p_1, group.by = "RNA_snn_res.0.4")
bcell_ids_1 <- WhichCells(p_1, expression = celltypist_max_high== "B cells" & RNA_snn_res.0.4 %in% c(0,2))
b_1 <- subset(p_1,  cells = bcell_ids_1)
DimPlot(b_1, group.by = "RNA_snn_res.0.4")

umap_coords_1 <- Embeddings(b_1, "umap.unintegrated")
b_1@meta.data$UMAP_1 <- umap_coords_1[,1]
b_1@meta.data$UMAP_2 <- umap_coords_1[,2]

b_1 <- subset(b_1,  subset = UMAP_1 < 5 & UMAP_2 < 5)
DimPlot(b_1, group.by = "RNA_snn_res.0.4")


DimPlot(p_2, group.by = "RNA_snn_res.0.2")
bcell_ids_2 <- WhichCells(p_2, expression = celltypist_max_high== "B cells" & RNA_snn_res.0.2 %in% c(0,2,4))
b_2 <- subset(p_2, cells = bcell_ids_2)
DimPlot(b_2, group.by = "RNA_snn_res.0.2")

umap_coords_2 <- Embeddings(b_2, "umap.unintegrated")
b_2@meta.data$UMAP_1 <- umap_coords_2[,1]
b_2@meta.data$UMAP_2 <- umap_coords_2[,2]

b_2 <- subset(b_2,  subset = UMAP_1 < 5)
DimPlot(b_2, group.by = "RNA_snn_res.0.4")



DimPlot(p_3, group.by = "RNA_snn_res.0.2")
bcell_ids_3 <- WhichCells(p_3, expression = celltypist_max_high== "B cells" & RNA_snn_res.0.2 %in% c(0,2))
b_3 <- subset(p_3, cells = bcell_ids_3)
DimPlot(b_3, group.by = "RNA_snn_res.0.2")

umap_coords_2 <- Embeddings(b_3, "umap.unintegrated")
b_3@meta.data$UMAP_1 <- umap_coords_2[,1]
b_3@meta.data$UMAP_2 <- umap_coords_2[,2]

b_3 <- subset(b_3,  subset = UMAP_1 < 5)
DimPlot(b_3, group.by = "RNA_snn_res.0.4")



DimPlot(p_4, group.by = "RNA_snn_res.0.2")
bcell_ids_4 <- WhichCells(p_4, expression = celltypist_max_high== "B cells" & RNA_snn_res.0.2 %in% c(0,1,2))
b_4 <- subset(p_4, cells = bcell_ids_4)
DimPlot(b_4, group.by = "RNA_snn_res.0.2")

umap_coords_2 <- Embeddings(b_4, "umap.unintegrated")
b_4@meta.data$UMAP_1 <- umap_coords_2[,1]
b_4@meta.data$UMAP_2 <- umap_coords_2[,2]

b_4 <- subset(b_4,  subset = UMAP_1 < 4)
DimPlot(b_4, group.by = "RNA_snn_res.0.4")



DimPlot(b_1, reduction = "umap.unintegrated", group.by = "celltypist_max_high")
DimPlot(b_2, reduction = "umap.unintegrated", group.by = "celltypist_max_high")
DimPlot(b_3, reduction = "umap.unintegrated", group.by = "celltypist_max_high")
DimPlot(b_4, reduction = "umap.unintegrated", group.by = "celltypist_max_high")



# data normalization
b_1 <- NormalizeData(b_1)
b_1 <- FindVariableFeatures(b_1)
b_1 <- ScaleData(b_1)
b_1 <- RunPCA(b_1)
ElbowPlot(b_1, ndims = 50)
b_1 <- FindNeighbors(b_1, dims = 1:30)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (r in resolutions) {
  b_1 <- FindClusters(b_1, resolution = r)
}
clustree_uninteg <- clustree(b_1@meta.data, prefix = "RNA_snn_res.")
clustree_uninteg
b_1 <- RunUMAP(b_1, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(b_1, reduction = "umap.unintegrated")

# data normalization
b_2 <- NormalizeData(b_2)
b_2 <- FindVariableFeatures(b_2)
b_2 <- ScaleData(b_2)
b_2 <- RunPCA(b_2)
ElbowPlot(b_2, ndims = 50)
b_2 <- FindNeighbors(b_2, dims = 1:30)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (r in resolutions) {
  b_2 <- FindClusters(b_2, resolution = r)
}
clustree_uninteg <- clustree(b_2@meta.data, prefix = "RNA_snn_res.")
clustree_uninteg
b_2 <- RunUMAP(b_2, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(b_2, reduction = "umap.unintegrated")


# data normalization
b_3 <- NormalizeData(b_3)
b_3 <- FindVariableFeatures(b_3)
b_3 <- ScaleData(b_3)
b_3 <- RunPCA(b_3)
ElbowPlot(b_3, ndims = 50)
b_3 <- FindNeighbors(b_3, dims = 1:30)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (r in resolutions) {
  b_3 <- FindClusters(b_3, resolution = r)
}
clustree_uninteg <- clustree(b_3@meta.data, prefix = "RNA_snn_res.")
clustree_uninteg
b_3 <- RunUMAP(b_3, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(b_3, reduction = "umap.unintegrated")

# data normalization
b_4 <- NormalizeData(b_4)
b_4 <- FindVariableFeatures(b_4)
b_4 <- ScaleData(b_4)
b_4 <- RunPCA(b_4)
ElbowPlot(b_4, ndims = 50)
b_4 <- FindNeighbors(b_4, dims = 1:30)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (r in resolutions) {
  b_4 <- FindClusters(b_4, resolution = r)
}
clustree_uninteg <- clustree(b_4@meta.data, prefix = "RNA_snn_res.")
clustree_uninteg
b_4 <- RunUMAP(b_4, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(b_4, reduction = "umap.unintegrated")


## 3. ANALYSIS 


### Functions for pathway enrichment analysis
plot_gsea <- function(data,
                      p_cutoff = 0.05,
                      top_n = 20) {
  
  df <- data@result
  
  df <- df %>%
    dplyr::filter(!is.na(NES), !is.na(p.adjust)) %>%
    dplyr::filter(p.adjust < p_cutoff) %>%
    dplyr::arrange(p.adjust) %>%
    head(top_n)
  
  ggplot(df, aes(
    x = NES,
    y = reorder(Description, NES),
    color = p.adjust,
    size = setSize
  )) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = "GO pathway",
      title = deparse(substitute(data)),
      color = "Adjusted p-value",
      size = "Gene set size"
    )
}

plot_kegg <- function(data,
                      p_cutoff = 0.05,
                      top_n = 20) {
  
  df <- data@result
  
  df <- df %>%
    dplyr::filter(!is.na(NES), !is.na(p.adjust)) %>%
    dplyr::filter(p.adjust < p_cutoff) %>%
    dplyr::arrange(p.adjust) %>%
    head(top_n)
  
  ggplot(df, aes(
    x = NES,
    y = reorder(Description, NES),
    color = p.adjust,
    size = setSize
  )) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = "KEGG pathway",
      title = deparse(substitute(data)),
      color = "Adjusted p-value",
      size = "Gene set size"
    )
}


run_GSEA_pipeline <- function(markers_df,
                              orgdb = org.Hs.eg.db,
                              pct_cutoff = 0.1,
                              seed = 1991) {
  
  library(dplyr)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  
  set.seed(seed)
  
  filt_df <- markers_df %>%
    filter(pct.1 > pct_cutoff | pct.2 > pct_cutoff) %>%
    arrange(desc(avg_log2FC))
  
  if ("cluster" %in% colnames(filt_df)) {
    markers_list <- split(filt_df, filt_df$cluster)
    names(markers_list) <- paste0("cl_", names(markers_list))
  } else {
    markers_list <- list(all_genes = filt_df)
  }
  
  msig_h <- msigdbr(
    species = "Homo sapiens",
    collection = "H"
  ) %>%
    dplyr::select(gs_name, gene_symbol)
  
  run_go <- function(df) {
    df$gene <- rownames(df)
    gene_list <- df$avg_log2FC
    names(gene_list) <- df$gene
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    tryCatch(
      gseGO(
        geneList = gene_list,
        OrgDb = orgdb,
        ont = "BP",
        keyType = "SYMBOL",
        minGSSize = 10,
        maxGSSize = 250,
        pvalueCutoff = 0.05,
        eps = 0,
        nPermSimple = 10000,
        seed = TRUE
      ),
      error = function(e) NULL
    )
  }
  
  run_kegg <- function(df, orgdb) {
    df$gene <- rownames(df)
    df$gene <- ifelse(!"gene" %in% colnames(df),
                      rownames(df),
                      df$gene)
    
    gene_df <- tryCatch(
      bitr(
        df$gene,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = orgdb
      ),
      error = function(e) NULL
    )
    
    if (is.null(gene_df)) return(NULL)
    
    df <- inner_join(df, gene_df, by = c("gene" = "SYMBOL"))
    
    gene_list <- df$avg_log2FC
    names(gene_list) <- df$ENTREZID
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    tryCatch(
      gseKEGG(
        geneList = gene_list,
        organism = "hsa",
        minGSSize = 10,
        maxGSSize = 250,
        pvalueCutoff = 0.05,
        eps = 0,
        nPermSimple = 10000,
        seed = TRUE
      ),
      error = function(e) NULL
    )
  }
  
  run_hallmark <- function(df) {
    df$gene <- rownames(df)
    gene_list <- df$avg_log2FC
    names(gene_list) <- df$gene
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    tryCatch(
      GSEA(
        geneList = gene_list,
        TERM2GENE = msig_h,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        eps = 0,
        pAdjustMethod = "BH",
        seed = TRUE
      ),
      error = function(e) NULL
    )
  }
  
  
  GO_results <- lapply(markers_list, run_go)
  KEGG_results <- lapply(markers_list, function(x) run_kegg(x, orgdb))
  Hallmark_results <- lapply(markers_list, run_hallmark)
  
  return(list(
    filtered_markers = filt_df,
    markers_list = markers_list,
    GO = GO_results,
    KEGG = KEGG_results,
    Hallmark = Hallmark_results
  ))
}

plot_multi_patient_GSEA <- function(
    results_list,
    analysis = c("GO", "Hallmark", "KEGG", "Reactome"),
    padj_cutoff = 0.05,
    top_n = 20,
    split_direction = FALSE
) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(pheatmap)
  
  analysis <- match.arg(analysis)
  
  all_df <- lapply(names(results_list), function(patient){
    
    patient_res <- results_list[[patient]][[analysis]]
    
    if (is.null(patient_res))
      return(NULL)
    
    # per si hi ha clusters
    subdfs <- lapply(names(patient_res), function(cl){
      
      gsea <- patient_res[[cl]]
      
      if (is.null(gsea))
        return(NULL)
      
      df <- as.data.frame(gsea)
      
      if (nrow(df) == 0)
        return(NULL)
      
      df$patient <- patient
      df$cluster <- cl
      
      return(df)
    })
    
    bind_rows(subdfs)
  })
  
  df <- bind_rows(all_df)
  
  if (nrow(df) == 0) {
    message("No GSEA results found")
    return(NULL)
  }
  
  df <- df %>%
    filter(p.adjust < padj_cutoff)
  
  if (nrow(df) == 0) {
    message("No significant pathways")
    return(NULL)
  }
  
  df <- df %>%
    mutate(direction = ifelse(NES > 0, "UP", "DOWN"))
  
  if (split_direction) {
    
    top_up <- df %>%
      filter(NES > 0) %>%
      group_by(Description) %>%
      summarise(score = max(NES)) %>%
      arrange(desc(score)) %>%
      slice_head(n = top_n)
    
    top_down <- df %>%
      filter(NES < 0) %>%
      group_by(Description) %>%
      summarise(score = min(NES)) %>%
      arrange(score) %>%
      slice_head(n = top_n)
    
    top_terms <- c(top_up$Description,
                   top_down$Description)
    
  } else {
    
    top_terms <- df %>%
      group_by(Description) %>%
      summarise(score = max(abs(NES))) %>%
      arrange(desc(score)) %>%
      slice_head(n = top_n) %>%
      pull(Description)
  }
  
  df_top <- df %>%
    filter(Description %in% top_terms)
  
  p_dot <- ggplot(
    df_top,
    aes(
      x = patient,
      y = fct_reorder(Description, NES),
      color = NES,
      size = -log10(p.adjust)
    )
  ) +
    geom_point(alpha = 0.9) +
    facet_wrap(~cluster, scales = "free_y") +
    scale_color_gradient2(
      low = "blue",
      mid = "white",
      high = "#BE4158",
      midpoint = 0
    ) +
    theme_bw() +
    labs(
      title = paste0(analysis, " enrichment"),
      x = "",
      y = ""
    )
  
  heat_df <- df_top %>%
    dplyr::select(patient, Description, NES) %>%
    group_by(patient, Description) %>%
    summarise(NES = mean(NES), .groups = "drop") %>%
    pivot_wider(
      names_from = patient,
      values_from = NES
    )
  
  mat <- as.data.frame(heat_df)
  
  rownames(mat) <- mat$Description
  mat$Description <- NULL
  
  mat <- as.matrix(mat)
  
  mat[is.na(mat)] <- 0
  
  mat <- mat[
    apply(mat, 1, function(x) sd(x, na.rm = TRUE) > 0),
    ,
    drop = FALSE
  ]
  
  mat <- mat[
    ,
    apply(mat, 2, function(x) sd(x, na.rm = TRUE) > 0),
    drop = FALSE
  ]
  
  mat[is.infinite(mat)] <- 0
  
  pheat <- pheatmap(
    mat,
    clustering_method = "complete",
    main = paste0(analysis, " NES heatmap")
  )
  
  
  p_violin <- ggplot(
    df,
    aes(
      x = patient,
      y = NES,
      fill = patient
    )
  ) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    theme_bw() +
    labs(
      title = paste0(analysis, " NES distribution")
    )
  
  return(list(
    results = df_top,
    dotplot = p_dot,
    heatmap = pheat,
    violin = p_violin
  ))
}


# LIBRARIES
library(Seurat)
library(R.utils)
library(dplyr)
library(Matrix)
library(ggplot2)
library(scDblFinder)
library(patchwork)
library(clustree)
library(SeuratObject)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GeneNMF)
library(msigdbr)
library(pheatmap)

# Hallmark lists
msig_h <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_list <- msig_h %>%
  dplyr::select(gs_name, gene_symbol)

patients <- list(
  P1 = b_1,
  P2 = b_2,
  P3 = b_3,
  P4 = b_4
)

genes <- c("EBF1", "IGHM")

df_list <- lapply(names(patients), function(p) {
  obj <- patients[[p]]
  expr <- FetchData(obj, vars = genes)
  expr$patient <- p
  expr
})

df <- bind_rows(df_list)
dot_df <- df %>%
  group_by(patient) %>%
  summarise(
    EBF1_avg = mean(EBF1),
    IGHM_avg = mean(IGHM),
    EBF1_pct = mean(EBF1 > 0) * 100,
    IGHM_pct = mean(IGHM > 0) * 100
  ) %>%
  tidyr::pivot_longer(
    cols = -patient,
    names_to = c("gene", ".value"),
    names_pattern = "(EBF1|IGHM)_(avg|pct)"
  )

ggplot(dot_df, aes(x = patient, y = gene)) +
  geom_point(aes(size = pct, color = avg)) +
  scale_color_gradient(low = "lightgrey", high = "#BE4158") +
  scale_size(range = c(2, 10)) +
  labs(
    x = "Patient",
    y = "Gene",
    color = "Avg expression",
    size = "% cells expressing"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) + guides(
    color = guide_colorbar(title.position = "top", title.hjust = 0.5),
    size  = guide_legend(title.position = "top")
  ) +
  theme(
    legend.position = "right",
    legend.box = "horizontal"
  )

################################
# 20-7073A #
################################

table(b_1$BTK_c.1442G.C_p.Cys481Ser)
table(b_1$PLCG2_c.2120C.T_p.Ser707Phe, b_1$PLCG2_c.2126A.T_p.Tyr709Phe)
ncol(b_1) 

b_1$mutated_BTK <- ifelse(b_1$BTK_c.1442G.C_p.Cys481Ser == "Hom","YES",
                                 ifelse(b_1$BTK_c.1442G.C_p.Cys481Ser == "NotAvailable","NA","NO"))

b_1$mutated_PLCG2 <- ifelse(b_1$PLCG2_c.2120C.T_p.Ser707Phe %in% c("Hom", "Het") 
                                   | b_1$PLCG2_c.2126A.T_p.Tyr709Phe %in% c("Hom", "Het"),"YES",
                                   ifelse(b_1$PLCG2_c.2120C.T_p.Ser707Phe == "NotAvailable"
                                          | b_1$PLCG2_c.2126A.T_p.Tyr709Phe == "NotAvailable","NA","NO"))

b_1$mutation <- ifelse(b_1$mutated_BTK == "NO" & b_1$mutated_PLCG2 != "YES", "WT",
                              ifelse(b_1$mutated_PLCG2 == "YES", "PLCG2",
                                     ifelse(b_1$mutated_BTK == "YES", "BTK", "NA")))

b_1$mutated <- ifelse(b_1$mutation %in% c("BTK", "PLCG2") , "YES",
                             ifelse(b_1$mutation == "WT", "NO", "NA"))


DimPlot(b_1, group.by = "mutated_BTK", cols = c("NA" = "grey90", "NO" = "#A6BDDB", "YES" = "#BE4158"))
DimPlot(b_1, group.by = "mutated_PLCG2", cols = c("NA" = "grey90", "NO" = "#A6BDDB", "YES" = "#238B45"))

DimPlot(b_1, group.by = "mutation", cols = c("NA" = "grey90", "WT" = "#2C7BB6","PLCG2" = "#238B45", "BTK" = "#BE4158"))
DimPlot(b_1, group.by = "mutated", cols = c("NA" = "grey90", "NO" = "#2C7BB6", "YES" = "#BE4158"))

m_1 <- subset(b_1, subset = mutated != "NA")
DimPlot(m_1, group.by = "mutated", cols = c("NO" = "#2C7BB6", "YES" = "#BE4158"))

table(m_1$mutation)


## heatmap co-ocurrence of mutations
mat <- table(
  b_1$mutated_BTK,
  b_1$mutated_PLCG2)
df <- as.data.frame(mat)
colnames(df) <- c("BTK", "PLCG2", "Count")
df <- df %>%
  group_by(BTK) %>%
  mutate(prop = Count / sum(Count))

ggplot(df, aes(x = PLCG2, y = BTK, fill = prop)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = Count), size = 5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#A6BDDB") +
  labs(
    x = "PLCG2 mutation status",
    y = "BTK mutation status",
    title = "BTK/PLCG2 mutation co-occurrence (Patient 1)",
    fill = "Proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


##### MARKERS BY MUTATION with all cells

Idents(m_1) <- m_1@meta.data$mutated
mark.1 <- FindMarkers(m_1, ident.1 = "NO", ident.2 = "YES", min.pct = 0.1)

res.1 <- run_GSEA_pipeline(mark.1)


##### MARKERS BY CLUSTERS with all cells
FeaturePlot(b_1, features = c("CXCR4", "CD27", "MIR155HG"))
FeaturePlot(m_1, features = c("CXCR4", "CD27", "MIR155HG"))

DimPlot(b_1, group.by = "RNA_snn_res.0.6")
table(m_1$RNA_snn_res.0.6, m_1$mutation)
prop.table(table(m_1$RNA_snn_res.0.6, m_1$mutation), margin = 1)

Idents(b_1) <- b_1@meta.data$RNA_snn_res.0.6
mark.cl.1 <- FindMarkers(b_1, ident.1 = "2")
res.cl.1 <- run_GSEA_pipeline(mark.cl.1)


## plot markers
DimPlot(b_1, group.by = "mutated", cols = c("NA" = "grey90", "NO" = "#2C7BB6", "YES" = "#BE4158"))
DimPlot(b_1, group.by = "RNA_snn_res.0.6")
prop.table(table(m_1$mutated, m_1$RNA_snn_res.0.6), margin = 2)

p1 <- DimPlot(b_1, group.by = "RNA_snn_res.0.6", cols=c(
  "#8E2F3D",
  "#BE4158",
  "#E86D86",
  "#F6C1C8",
  "#6BAED6",
  "#2171B5"
)) +  labs(x = "UMAP 1", y = "UMAP 2",title = "Patient 1") +
  guides(color = guide_legend(title = "Cluster identity",override.aes = list(size = 4)))

p2 <- DimPlot(b_1, group.by = "mutation", order = c("PLCG2", "WT","BTK", "NA"),
              cols = c("NA" = "grey90", "WT" = "#A6BDDB", "BTK" = "#BE4158", "PLCG2" = "forestgreen")) +
  labs(x = "UMAP 1", y = "UMAP 2") +ggtitle(NULL) +
  guides(color = guide_legend(title = "Mutation",override.aes = list(size = 4)))


df <- as.data.frame(prop.table(
  table(m_1$RNA_snn_res.0.6, m_1$mutation),
  margin = 1
))

colnames(df) <- c("cluster", "mutation", "proportion")

# Barplot
p3 <- ggplot(df, aes(x = cluster, y = proportion, fill = mutation)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c(
    "NA" = "grey90",
    "WT" = "#A6BDDB",
    "BTK" = "#BE4158",
    "PLCG2" = "#1B9E77"
  )) + 
  labs(
    x = "Cluster",
    y = "Proportion of cells",
    fill = "Mutation"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )


wrap_plots(p2,p1, ncol=1)
p1 / (p2 | p3)  + plot_layout(heights = c(1.5, 1))


################################
# 2 #
################################

table(b_2$BTK_c.1441T.A_p.Cys481Ser, b_2$BTK_c.1442G.C_p.Cys481Ser)

b_2$mutated_BTK<- ifelse(b_2$BTK_c.1441T.A_p.Cys481Ser %in% c("Hom", "Het") 
                               | b_2$BTK_c.1442G.C_p.Cys481Ser %in% c("Hom", "Het"),"YES",
                               ifelse(b_2$BTK_c.1441T.A_p.Cys481Ser == "NotAvailable"
                                      | b_2$BTK_c.1442G.C_p.Cys481Ser == "NotAvailable","NA","NO"))

DimPlot(b_2, group.by = "mutated_BTK", cols = c("NA" = "grey90", "NO" = "#2C7BB6", "YES" = "#BE4158"))

m_2 <- subset(b_2, subset = mutated_BTK != "NA")
DimPlot(m_2, group.by = "mutated_BTK", cols = c("NO" = "#2C7BB6", "YES" = "#BE4158"))

## heatmap co-ocurrence of BTK mutations
mat <- table(
  b_2$BTK_c.1441T.A_p.Cys481Ser,
  b_2$BTK_c.1442G.C_p.Cys481Ser)
df <- as.data.frame(mat)
colnames(df) <- c("BTK_c.1441T.A.", "BTK_c.1442G.C", "Count")
df[] <- lapply(df, as.character)
df[df == "Wt"] <- "NO"
df[df == "Hom"] <- "YES"
df[df == "NotAvailable"] <- "NA"
df$Count <- as.numeric(df$Count)
df <- df %>%
  group_by(BTK_c.1441T.A.) %>%
  mutate(prop = Count / sum(Count))

ggplot(df, aes(x = BTK_c.1442G.C, y = BTK_c.1441T.A., fill = prop)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = Count), size = 5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#A6BDDB") +
  labs(
    x = "BTK_c.1442G.C mutation status",
    y = "BTK_c.1441T.A. mutation status",
    title = "BTK mutations co-occurrence (Patient 2)",
    fill = "Proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )



##### MARKERS BY MUTATION with all cells

Idents(m_2) <- m_2@meta.data$mutated_BTK
mark.2 <- FindMarkers(m_2, ident.1 = "NO", ident.2 = "YES",min.pct = 0.1)
res.2 <- run_GSEA_pipeline(mark.2)


##### MARKERS BY CLUSTERS with all cells
FeaturePlot(b_2, features = c("CXCR4", "CD27", "MIR155HG"))
FeaturePlot(m_2, features = c("CXCR4", "CD27", "MIR155HG"))

DimPlot(b_2, group.by = "RNA_snn_res.0.6")
table(m_2$RNA_snn_res.0.6, m_2$mutated_BTK)

table(m_2$RNA_snn_res.0.6, m_2$mutated_BTK)
prop.table(table(m_2$RNA_snn_res.0.6, m_2$mutated_BTK), margin = 1)

Idents(b_2) <- b_2@meta.data$RNA_snn_res.0.6
mark.cl.2 <- FindMarkers(b_2, ident.1 = "0")
res.cl.2 <- run_GSEA_pipeline(mark.cl.2)


################################
# 3 #
################################

table(b_3$BTK_c.1441T.A_p.Cys481Ser)

b_3$mutated_BTK<- ifelse(b_3$BTK_c.1441T.A_p.Cys481Ser %in% c("Hom", "Het") ,"YES",
                               ifelse(b_3$BTK_c.1441T.A_p.Cys481Ser == "Wt", "NO","NA"))

DimPlot(b_3, group.by = "mutated_BTK", cols = c("NA" = "grey90", "NO" = "#2C7BB6", "YES" = "#BE4158"))
m_3 <- subset(b_3, subset = mutated_BTK != "NA")
DimPlot(m_3, group.by = "mutated_BTK", cols = c("NO" = "#2C7BB6", "YES" = "#BE4158"))

mat <- table(b_3$BTK_c.1441T.A_p.Cys481Ser)
df <- as.data.frame(mat)
colnames(df) <- c("BTK_c.1441T.A.", "Count")

df$BTK_c.1441T.A. <- recode(df$BTK_c.1441T.A.,
                            "Wt" = "NO",
                            "Hom" = "YES",
                            "NotAvailable" = "NA")
df$BTK_c.1441T.A. <- factor(df$BTK_c.1441T.A.,
                            levels = c("NA", "NO", "YES"))
df$prop <- df$Count / sum(df$Count)
ggplot(df, aes(x = BTK_c.1441T.A.,
               y = 1,
               fill = prop)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = Count), size = 5,fontface = "bold") +
  scale_y_continuous(NULL, breaks = NULL) +
  scale_fill_gradient(low = "white", high = "#A6BDDB") +
  labs(
    x = "BTK_c.1442G.C mutation status",
    title = "BTK mutation (Patient 3)",
    fill = "Proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

##### MARKERS BY MUTATION with all cells

Idents(m_3) <- m_3@meta.data$mutated_BTK
mark.3 <- FindMarkers(m_3, ident.1 = "NO", ident.2 = "YES", min.pct = 0.1)

res.3 <- run_GSEA_pipeline(mark.3)


##### MARKERS BY CLUSTERS with all cells
FeaturePlot(b_3, features = c("CXCR4", "CD27", "MIR155HG"))
FeaturePlot(m_3, features = c("CXCR4", "CD27", "MIR155HG"))
p1 <- DimPlot(b_3, group.by = "mutated_BTK", 
              cols = c("NA" = "grey90", "NO" = "#2C7BB6", "YES" = "#BE4158")) +
  labs(x = "UMAP 1", y = "UMAP 2") + ggtitle("Patient 3")
p2 <- DimPlot(b_3, group.by = "RNA_snn_res.0.6", cols = c(
  "#C6DBEF",
  "#6BAED6",
  "#08306B",
  "#4292C6",
  "#2171B5"
)) + labs(x = "UMAP 1", y = "UMAP 2") + ggtitle(NULL)
wrap_plots(p1, p2, ncol=1)


DimPlot(b_3, group.by = "RNA_snn_res.0.6", label = TRUE)
table(m_3$RNA_snn_res.0.6, m_3$mutated_BTK)
prop.table(table(m_3$RNA_snn_res.0.6, m_3$mutated_BTK), margin=1)

Idents(b_3) <- b_3@meta.data$RNA_snn_res.0.6
mark.cl.3 <- FindMarkers(b_3, ident.1 = "2")
res.cl.3 <- run_GSEA_pipeline(mark.cl.3)

################################
# 4 #
################################

table(b_4$BTK_c.1442G.C_p.Cys481Ser)

b_4$mutated_BTK<- ifelse(b_4$BTK_c.1442G.C_p.Cys481Ser %in% c("Hom", "Het") ,"YES",
                               ifelse(b_4$BTK_c.1442G.C_p.Cys481Ser == "NotAvailable", "NA","NO"))

DimPlot(b_4, group.by = "mutated_BTK", cols = c("NA" = "grey90", "NO" = "#2C7BB6", "YES" = "#BE4158"))
m_4 <- subset(b_4, subset = mutated_BTK != "NA")
DimPlot(m_4, group.by = "mutated_BTK", cols = c("NA" = "grey90", "NO" = "#2C7BB6", "YES" = "#BE4158"))


mat <- table(b_4$BTK_c.1442G.C_p.Cys481Ser)
df <- as.data.frame(mat)
colnames(df) <- c("BTK_c.1442G.C.", "Count")

df$BTK_c.1442G.C. <- recode(df$BTK_c.1442G.C.,
                            "Wt" = "NO",
                            "Hom" = "YES",
                            "NotAvailable" = "NA")
df$BTK_c.1442G.C. <- factor(df$BTK_c.1442G.C.,
                            levels = c("NA", "NO", "YES"))
df$prop <- df$Count / sum(df$Count)
ggplot(df, aes(x = BTK_c.1442G.C.,
               y = 1,
               fill = prop)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = Count), size = 5,fontface = "bold") +
  scale_y_continuous(NULL, breaks = NULL) +
  scale_fill_gradient(low = "white", high = "#A6BDDB") +
  labs(
    x = "BTK_c.1442G.C mutation status",
    title = "BTK mutation (Patient 4)",
    fill = "Proportion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

##### MARKERS BY MUTATION with all cells

Idents(m_4) <- m_4@meta.data$mutated_BTK
mark.4 <- FindMarkers(m_4, ident.1 = "NO", ident.2 = "YES", min.pct = 0.1)
res.4 <- run_GSEA_pipeline(mark.4)


##### MARKERS BY CLUSTERS with all cells
FeaturePlot(b_4, features = c("CXCR4", "CD27", "MIR155HG"))
FeaturePlot(m_4, features = c("CXCR4", "CD27", "MIR155HG"))

DimPlot(b_4, group.by = "RNA_snn_res.0.6")
prop.table(table(m_4$RNA_snn_res.0.6, m_4$mutated_BTK), margin=1)

Idents(b_4) <- b_4@meta.data$RNA_snn_res.0.6
mark.cl.4 <- FindMarkers(b_4, ident.1 = "3")

res.cl.4 <- run_GSEA_pipeline(mark.cl.4)


### ANALYSIS OF COMMON GENES FOR MUTATED / UNMUTATED CELLS

markers_all <- bind_rows(
  mark.1 %>% mutate(patient = "p_20_7079A"),
  mark.2 %>% mutate(patient = "p_2"),
  mark.4 %>% mutate(patient = "p_4")
)
head(markers_all)


genes_WT <- markers_all %>%
  filter(cluster == "NO") %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>%
  distinct(patient, gene) %>%
  group_by(gene) %>%
  summarise(
    n_patients = n(),
    patients = paste(patient, collapse = ", ")
  ) %>%
  filter(n_patients > 1) %>%
  arrange(desc(n_patients))

genes_WT


## % genotyped cells per cluster
patients <- c("20-7073A", "23-0790", "23-6687", "1668-01")
pat <- gsub("-", "_", patients)

plots_genotyped_cells <- list()

patient_labels <- c(
  "1" = "Patient 1",
  "2"  = "Patient 2",
  "3"  = "Patient 3",
  "4"  = "Patient 4"
)

for (i in pat) {
  
  sample <- get(paste0("b_", i))
  
  if (i == "1"){
    genotyped_prop <- prop.table(table(sample$mutated, sample$RNA_snn_res.0.6), margin = 2)
  } else {
    genotyped_prop <- prop.table(table(sample$mutated_BTK, sample$RNA_snn_res.0.6), margin = 2)
  }
  
  df <- as.data.frame(genotyped_prop)
  colnames(df) <- c("mutation", "cluster", "proportion")
  df$cluster <- factor(as.numeric(as.character(df$cluster)),
                       levels = sort(unique(as.numeric(as.character(df$cluster)))))
  df$genotyped <- ifelse(df$mutation %in% c("YES", "NO"), "YES", "NO")
  
  plot <- ggplot(df, aes(cluster, proportion, fill = genotyped)) +
    geom_col() +
    labs(
      title = patient_labels[[i]],
      x = "Cluster",
      y = "Proportion"
    ) +
    scale_fill_manual(values = c(
      "NO" = "#A6BDDB",
      "YES" = "#BE4158"
    )) +
    theme_linedraw() +
    theme(legend.position = "right")
  
  plots_genotyped_cells[[i]] <- plot
}

p9  <- plots_genotyped_cells[["1"]]
p10 <- plots_genotyped_cells[["2"]]
p11 <- plots_genotyped_cells[["3"]]
p12 <- plots_genotyped_cells[["4"]]

library(patchwork)

(p9 + p10 + p11 + p12) +
  plot_layout(ncol = 4, guides = "collect") &
  theme(legend.position = "right")


## volcano plots
library(EnhancedVolcano)
EnhancedVolcano(mark.1, lab = rownames(mark.1), x = "avg_log2FC", y="p_val_adj",
                pCutoff = 0.05, FCcutoff = 0.5, pointSize = 3, labSize = 5,
                drawConnectors = TRUE, title = "Patient 1")

EnhancedVolcano(mark.2, lab = rownames(mark.2), x = "avg_log2FC", y="p_val_adj",
                pCutoff = 0.05, FCcutoff = 0.5, pointSize = 3, labSize = 5,
                drawConnectors = TRUE, title = "Patient 2")

EnhancedVolcano(mark.3, lab = rownames(mark.3), x = "avg_log2FC", y="p_val_adj",
                pCutoff = 0.05, FCcutoff = 0.5, pointSize = 3, labSize = 5,
                drawConnectors = TRUE, title = "3 - WT", max.overlaps = 15)

EnhancedVolcano(mark.4, lab = rownames(mark.4), x = "avg_log2FC", y="p_val_adj",
                pCutoff = 0.05, FCcutoff = 0.1, pointSize = 3, labSize = 5,
                drawConnectors = TRUE, title = "4 - WT")


# scatter plot WT - CXCR4
scatter_plot <- function(obj, patient_name, cluster_col = "RNA_snn_res.0.6", features) {
  
  # WT fraction per cluster
  df_wt <- obj@meta.data %>%
    group_by(.data[[cluster_col]]) %>%
    summarise(
      wt_fraction = mean(mutated_BTK == "NO"),
      .groups = "drop"
    ) %>%
    rename(cluster = all_of(cluster_col))
  
  feat_expr <- AverageExpression(
    obj,
    features = features,
    group.by = cluster_col,
    assays = "RNA"
  )$RNA
  
  df_expr <- data.frame(
    cluster = gsub("^g", "", colnames(feat_expr)),
    feat_expr = as.numeric(feat_expr[1,])
  )
  
  # merge
  df <- left_join(df_wt, df_expr, by = "cluster")
  
  # ensure numeric cluster ordering
  df$cluster <- as.numeric(as.character(df$cluster))
  
  # plot
  p <- ggplot(df, aes(x = wt_fraction, y = feat_expr)) +
    
    geom_point(size = 4, alpha = 0.85) +
    
    geom_smooth(method = "lm", se = FALSE,
                color = "black", linewidth = 0.8) +
    
    ggrepel::geom_text_repel(
      aes(label = cluster),
      size = 4,
      box.padding = 0.5,
      point.padding = 0.8,
      max.overlaps = Inf
    ) +
    labs(
      title = patient_name,
      x = "Fraction of wild-type cells per cluster",
      y = paste0(features, " expression (mean)")
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black")
    )
  
  return(p)
}

p1 <- scatter_plot(m_1, patient_name = "Patient 1", features = "CXCR4")
p2 <- scatter_plot(m_2, patient_name = "Patient 2", features = "CXCR4")
p3 <- scatter_plot(m_4, patient_name = "Patient 4", features = "CXCR4")

(p1 + p2 + p3) +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Association between WT cells proportion and CXCR4 expression",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold")
    ))



### test NF-KB and BCR signatures

hallmark_pathways <- msigdbr(
  species = "Homo sapiens",
  collection = "H")

tnfa_genes <- hallmark_pathways %>%
  dplyr::filter(gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB") %>%
  dplyr::pull(gene_symbol)

b_1 <- AddModuleScore(b_1, list(tnfa_genes), name = "NF-KB signature")
b_2 <- AddModuleScore(b_2, list(tnfa_genes), name = "NF-KB signature")
b_3 <- AddModuleScore(b_3, list(tnfa_genes), name = "NF-KB signature")
b_4 <- AddModuleScore(b_4, list(tnfa_genes), name = "NF-KB signature")


# ridgeplot --> mut vermell, WT blau (order de + a - mut)
order <- c("3", "2", "1", "4", "0", "5")
b_1$RNA_snn_res.0.6 <- factor(
  b_1$RNA_snn_res.0.6,
  levels = order)
Idents(b_1) <- b_1@meta.data$RNA_snn_res.0.6
plot6 <- RidgePlot(b_1, features = "NF-KB signature1", 
                   cols = c("#BE4158","#BE4158","#BE4158","#BE4158","#A6BDDB","#A6BDDB")) +
  ggtitle("Patient 1")

order <- c("3", "1", "2", "4","0")
b_2$RNA_snn_res.0.6 <- factor(
  b_2$RNA_snn_res.0.6,
  levels = order)
Idents(b_2) <- b_2@meta.data$RNA_snn_res.0.6
plot7 <- RidgePlot(b_2, features = "NF-KB signature1", 
                   cols = c("#BE4158","#BE4158","#BE4158","#BE4158","#A6BDDB")) +
  ggtitle("Patient 2")


order <- c("0", "3","2", "4", "1")
b_3$RNA_snn_res.0.6 <- factor(
  b_3$RNA_snn_res.0.6,
  levels = order)
Idents(b_3) <- b_3@meta.data$RNA_snn_res.0.6
plot8 <- RidgePlot(b_3, features = "NF-KB signature1", 
                   cols = c("#A6BDDB","#A6BDDB","#A6BDDB","#A6BDDB","#A6BDDB")) +
  ggtitle("Patient 3")


order <- c( "0", "2", "1", "4", "3")
b_4$RNA_snn_res.0.6 <- factor(
  b_4$RNA_snn_res.0.6,
  levels = order)
Idents(b_4) <- b_4@meta.data$RNA_snn_res.0.6
plot9 <- RidgePlot(b_4, features = "NF-KB signature1", 
                   cols = c("#BE4158","#BE4158","#BE4158","#BE4158","#A6BDDB")) +
  ggtitle("Patient 4")

plot6 <- plot6 + labs(x = NULL, y = NULL)
plot7 <- plot7 + labs(x = NULL, y = NULL)
plot8 <- plot8 + labs(x = NULL, y = NULL)
plot9 <- plot9 + labs(x = NULL, y = NULL)

library(cowplot)

p <- (plot6 + plot7  + plot9) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(
    axis.title = element_blank()
  )

ggdraw() +
  draw_plot(p, x = 0.1, y = 0.1, width = 0.88, height = 0.88) +

  draw_label(
    "Cluster color:",
    x = 0.78, y = 0.18,
    size = 12,
    fontface = "bold"
  ) +
  annotate("rect",
           xmin = 0.78, xmax = 0.80,
           ymin = 0.12, ymax = 0.14,
           fill = "#A6BDDB", color = "#A6BDDB") +
  draw_label(
    "WT-enriched",
    x = 0.81, y = 0.13,
    size = 11,
    hjust = 0
  ) +
  annotate("rect",
           xmin = 0.78, xmax = 0.80,
           ymin = 0.07, ymax = 0.09,
           fill = "#BE4158", color = "#BE4158") +
  draw_label(
    "Mutated-enriched",
    x = 0.81, y = 0.08,
    size = 11,
    hjust = 0
  ) +
  
  draw_label(
    "NF-κB signature score",
    x = 0.55, y = 0.1,
    size = 14
  ) +
  draw_label(
    "Cluster identity",
    x = 0.09, y = 0.5,
    angle = 90,
    size = 14
  )


## plot pathways per clusters
df_1 <- res.cl.1$Hallmark$all_genes@result %>%
  mutate(patient = "Patient 1")

df_2 <- res.cl.2$Hallmark$all_genes@result %>%
  mutate(patient = "Patient 2")

df_3 <- res.cl.4$Hallmark$all_genes@result %>%
  mutate(patient = "Patient 4")

df_all <- bind_rows(df_1, df_2, df_3)

df_all <- df_all %>%
  rename(
    pathway = ID, 
    nes = NES,
    padj = p.adjust,
    size = setSize
  )

df_all <- df_all %>% filter(padj < 0.05)

df_all$pathway <- factor(df_all$pathway,
                         levels = unique(df_all$pathway[order(df_all$patient)]))

ggplot(df_all, aes(x = nes,
                   y = reorder(pathway, nes),
                   color = padj,
                   size = size)) +
  geom_point(alpha = 0.9) +
  facet_grid(patient ~ ., scales = "free_y", space = "free_y") +
  scale_color_gradient(low = "#d73027", high = "#4575b4") +
  theme_bw() +
  labs(x = "Normalized Enrichment Score (NES)",
       y = "Hallmark pathways",
       color = "Adjusted p-value",
       size = "Gene set size") +
  theme(
    strip.text.y = element_text(angle = 0),
    panel.grid.major.y = element_blank()
  )



df_1 <- res.wt.1$Hallmark$all_genes@result %>%
  mutate(patient = "1")

df_2 <- res.wt.2$Hallmark$all_genes@result %>%
  mutate(patient = "2")

df_3 <- res.wt.4$Hallmark$all_genes@result %>%
  mutate(patient = "4")

df_all <- bind_rows(df_1, df_2, df_3)

df_all <- df_all %>%
  rename(
    pathway = ID, 
    nes = NES,
    padj = p.adjust,
    size = setSize
  )

df_all <- df_all %>% filter(padj < 0.05)

df_all$pathway <- factor(df_all$pathway,
                         levels = unique(df_all$pathway[order(df_all$patient)]))

ggplot(df_all, aes(x = nes,
                   y = reorder(pathway, nes),
                   color = padj,
                   size = size)) +
  geom_point(alpha = 0.9) +
  facet_grid(patient ~ ., scales = "free_y", space = "free_y") +
  scale_color_gradient(low = "#BE4158", high = "#F6C1C8") +
  theme_bw() +
  labs(x = "Normalized Enrichment Score (NES)",
       y = "Hallmark pathways",
       color = "Adjusted p-value",
       size = "Gene set size") +
  theme(
    strip.text.y = element_text(angle = 0),
    panel.grid.major.y = element_blank()
  )

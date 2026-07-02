#### QUALITY CONTROL #### 

# libraries
library(Seurat)
library(SeuratObject)
library(R.utils)
library(dplyr)
library(Matrix)
library(ggplot2)
library(scDblFinder)
library(patchwork)

# set working directory + sample path
setwd("")
sample_path <- ""

# name of the samples in files
samples_files <-  c("")
# name of the samples in R
samples <- c("")


### UPLOAD OBJECTS

## function to read the files 
create_seurat_objects <- function(samples_files) {
  seurat_list <- list()
  for (i in samples_files) {
    path <- paste0(
      sample_path,
      i,"/DGE_unfiltered")
    
    i_data <- ReadParseBio(path)
    
    # Check to see if empty gene names are present, add name if so.
    table(rownames(i_data) == "")
    rownames(i_data)[rownames(i_data) == ""] <- "unknown"
    
    # Read in cell meta data
    cell_meta <- read.csv(paste0(path, "/cell_metadata.csv"), row.names = 1)
    
    seurat_i <- CreateSeuratObject(
      counts = i_data,
      min.features = 1,
      min.cells = 1,
      meta.data = cell_meta
    )
    assign(paste0("obj_",gsub("-","_",i)), seurat_i,envir = .GlobalEnv)
  }
}

create_seurat_objects(samples_files)

## change the active identity + add % of mt genes 
obj <- list() 
for (i in samples) {
  name <- paste0("obj_", i)
  sample <- get(name)
  sample@meta.data$orig.ident <- factor(rep(i, nrow(sample@meta.data)))
  # change the cell identity slot to the sample (before: plate well numbers)
  Idents(sample) <- sample@meta.data$orig.ident
  
  # add percentage of mitochondrial genes expressed in each cell
  sample <- PercentageFeatureSet(sample, pattern = "^MT-", col.name = "percent.mt")
  obj[[i]] <- sample
}


# violin plots for counts, features and percent.mt after hard_filtering
meta_obj <- lapply(names(obj), function(i) {
  df <- obj[[i]]@meta.data
  df$sample <- i
  return(df)
})
meta_obj <- bind_rows(meta_obj)

vln_seurat_style <- function(df, feature) {
  ggplot(df, aes(x = sample, y = .data[[feature]], fill = sample)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(x = NULL, y = feature)
}
plot_counts_obj   <- vln_seurat_style(meta_obj, "nCount_RNA")
plot_features_obj <- vln_seurat_style(meta_obj, "nFeature_RNA")
plot_mt_obj       <- vln_seurat_style(meta_obj, "percent.mt")

plot_counts_obj + plot_features_obj + plot_mt_obj
ggsave("", plot_counts_obj + plot_features_obj + plot_mt_obj, width = 12, height = 5)


#######################
### QUALITY CONTROL ###
#######################

## filtering with the trailmaker cell size distribution filter (what are considered as cells)
# minimum number of transcripts per cell
filter_tm <- c(1375.385, 1111.67, 628.21, 730.115, 1255.08, 363.729, 820.542, 958.785, 998.52, 1191.805, 568.906)
cell_filter <- data.frame(samples, filter_tm) 
cells <- list()
for (i in names(obj)) {
  sample <- obj[[i]]
  limit <- cell_filter$filter_tm[cell_filter$samples == i]
  sample <- subset(sample, subset = nCount_RNA > limit)
  
  cells[[i]] <- sample
}

## table of number of cells before and after filtering
table_cells <- data.frame()
for (i in names(obj)) {
  sample_pre <- obj[[i]]
  sample_post <- cells[[i]]
  cells_pre = (ncol(sample_pre))
  cells_post = (ncol(sample_post))
  df_temp <- data.frame(
    sample = i,
    n_cells_pre = cells_pre,
    n_cells_post = cells_post)
  table_cells <- rbind(table_cells, df_temp)
}
table_cells

# violin plots for counts, features and percent.mt after hard_filtering
meta_cells <- lapply(names(cells), function(i) {
  df <- cells[[i]]@meta.data
  df$sample <- i
  return(df)
})
meta_cells <- bind_rows(meta_cells)

plot_counts_cells   <- vln_seurat_style(meta_cells, "nCount_RNA")
plot_features_cells <- vln_seurat_style(meta_cells, "nFeature_RNA")
plot_mt_cells       <- vln_seurat_style(meta_cells, "percent.mt")

plot_counts_cells + plot_features_cells + plot_mt_cells
ggsave("", plot_counts_cells + plot_features_cells + plot_mt_cells, width = 12, height = 5)



## hard filtering
hard_filter <- list()
for (i in names(cells)) {
  sample <- cells[[i]]
  
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


# violin plots for counts, features and percent.mt after hard_filtering
meta_hard_filter <- lapply(names(hard_filter), function(i) {
  df <- hard_filter[[i]]@meta.data
  df$sample <- i
  return(df)
})
meta_hard_filter <- bind_rows(meta_hard_filter)

plot_counts   <- vln_seurat_style(meta_hard_filter, "nCount_RNA")
plot_features <- vln_seurat_style(meta_hard_filter, "nFeature_RNA")
plot_mt       <- vln_seurat_style(meta_hard_filter, "percent.mt")

plot_counts + plot_features + plot_mt
ggsave("", plot_counts + plot_features + plot_mt, width = 12, height = 5)

## hard-filter max counts before finding doublets
filters_db <- list()
for (i in names(cells)) {
  sample <- cells[[i]] # filter defined from original cells
  max_counts <- median(log10(sample@meta.data$nCount_RNA)) +  3*mad(log10(sample@meta.data$nCount_RNA), constant = 1)
  filters_db[[i]] <- 10^max_counts
}
hard_filter2 <- list()
for (i in names(hard_filter)) {
  sample <- hard_filter[[i]]
  sample <- subset(sample, subset = nCount_RNA < filters_db[[i]])
  hard_filter2[[i]] <- sample
}

# violin plots for counts, features and percent.mt after hard_filtering max counts
meta_hard_filter2 <- lapply(names(hard_filter2), function(i) {
  df <- hard_filter2[[i]]@meta.data
  df$sample <- i
  return(df)
})
meta_hard_filter2 <- bind_rows(meta_hard_filter2)

plot_counts2   <- vln_seurat_style(meta_hard_filter2, "nCount_RNA")
plot_features2 <- vln_seurat_style(meta_hard_filter2, "nFeature_RNA")
plot_mt2      <- vln_seurat_style(meta_hard_filter2, "percent.mt")
plot_counts2 + plot_features2 + plot_mt2
ggsave("", plot_counts2 + plot_features2 + plot_mt2, width = 12, height = 5)


## find doublets --> scDblFinder
hard_db <- list()
for (i in names(hard_filter2)) {
  sample <- hard_filter2[[i]]
  set.seed(1991)
  doublets <- scDblFinder(as.SingleCellExperiment(sample), clusters=FALSE)
  table(rownames(colData(doublets)) == rownames(sample@meta.data))
  sample@meta.data <- cbind(sample@meta.data, colData(doublets)[match(rownames(sample@meta.data), rownames(colData(doublets))), 17:20])
  # en el match x[match(y,x),]
  hard_db[[i]] <- sample
}

# table doublets
table_db <- data.frame()
for (i in names(hard_db)) {
  sample <- hard_db[[i]]
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

meta_hard_db <- lapply(names(hard_db), function(i) {
  df <- hard_db[[i]]@meta.data
  df$sample <- i
  return(df)
})
meta_hard_db <- bind_rows(meta_hard_db)
db_score <- vln_seurat_style(meta_hard_db, "scDblFinder.score")
ggsave("", db_score)


## FILTERING 

# define filter thresholds 
### ATENCIÓ: els thresholds s'han de calcular amb els counts abans de fer cap filtratge

## for percent.mt
filters_mt <- list()
for (i in names(cells)) {
  sample <- cells[[i]]
  max_mt <- median(log10(sample@meta.data$percent.mt)) +  2*mad(log10(sample@meta.data$percent.mt), constant = 1)
  min_mt <- median(log10(sample@meta.data$percent.mt)) -  2*mad(log10(sample@meta.data$percent.mt), constant = 1)
  filter_mt <- c(10^min_mt, 10^max_mt)
  filters_mt[[i]] <- filter_mt
}

## for counts 
filters_counts <- list()
medians_counts <- list()
for (i in names(cells)) {
  sample <- cells[[i]]
  median <- 10^(median(log10(sample@meta.data$nCount_RNA)))
  max_counts <- median(log10(sample@meta.data$nCount_RNA)) +  2*mad(log10(sample@meta.data$nCount_RNA), constant = 1)
  min_counts <- median(log10(sample@meta.data$nCount_RNA)) -  2*mad(log10(sample@meta.data$nCount_RNA), constant = 1)
  filter_counts <- c(10^min_counts, 10^max_counts)
  filters_counts[[i]] <- filter_counts
  medians_counts[[i]] <- median
}

## for doublets --> scDblFinder.class == "singlet"


## plot the counts thresholds 
for (i in names(hard_db)) {
  
  obj <- hard_db[[i]]
  d <- density(obj@meta.data$nCount_RNA)
  
  plot(d, main = paste0(i, " - nCount_RNA"))
  
  abline(v = filters_counts[[i]], col = "red", lwd = 2, lty = 2)
  abline(v = medians_counts[[i]], col = "blue", lwd = 2, lty = 2)
}

# Filtering
filtered <- list()
for (i in names(hard_db)) {
  sample <- hard_db[[i]]
  # filtering
  sample <- subset(sample, subset = 
                     percent.mt < filters_mt[[i]][[2]] &
                     scDblFinder.class == "singlet"&
                     nCount_RNA > filters_counts[[i]][[1]] & 
                     nCount_RNA < filters_counts[[i]][[2]])
  filtered[[i]] <- sample
}

# violin plots for counts, features and percent.mt after all filtering
meta_filtered <- lapply(names(filtered), function(i) {
  df <- filtered[[i]]@meta.data
  df$sample <- i
  return(df)
})
meta_filtered <- bind_rows(meta_filtered)

plot_counts_filt   <- vln_seurat_style(meta_filtered, "nCount_RNA")
plot_features_filt <- vln_seurat_style(meta_filtered, "nFeature_RNA")
plot_mt_filt       <- vln_seurat_style(meta_filtered, "percent.mt")
plot_counts_filt + plot_features_filt + plot_mt_filt

ggsave("", plot_counts_filt + plot_features_filt + plot_mt_filt
       , width = 12, height = 5)

db_score2 <- vln_seurat_style(meta_filtered, "scDblFinder.score")
ggsave("", db_score2)

### 2. CONTAMINATION FILTERING for samples of scRNA-seq from PARSE ###

# libraries
library(Seurat)
library(R.utils)
library(dplyr)
library(Matrix)
library(ggplot2)
library(scDblFinder)
library(SeuratObject)
library(patchwork)
library(clustree)


# merge objects
merged_obj <- merge(
  filtered[[1]],
  y = filtered[-1]
)


# check that the idents are saved correctly and join layers
all(sapply(filtered, ncol) == as.vector(table(merged_obj$orig.ident)))
merged_obj <- JoinLayers(merged_obj)

# add patient column to metadata and separate counts layers by patient
merged_obj@meta.data$patient <- sub("_.*", "", merged_obj@meta.data$orig.ident)
merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], f = merged_obj$patient)

# data normalization
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)
ElbowPlot(merged_obj, ndims = 50)
merged_obj <- FindNeighbors(merged_obj, dims = 1:10)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (r in resolutions) {
  merged_obj <- FindClusters(merged_obj, resolution = r)
}
clustree(merged_obj@meta.data, prefix = "RNA_snn_res.")

merged_obj <- RunUMAP(merged_obj, dims = 1:10, reduction = "pca", reduction.name = "umap.unintegrated")

# pre-integration plots
p1 <- DimPlot(merged_obj, reduction = "umap.unintegrated", label = TRUE)
p2 <- DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = "patient", label = TRUE)
wrap_plots(p1,p2)


## REMOVE CONTAMINATION

# in this plot patients are defined by only 1 or 2 clusters
p3 <- DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = "RNA_snn_res.0.2", label = TRUE)
p4 <- DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = "patient", label = TRUE)
wrap_plots(p3,p4)
ggsave("", p3+p4, width = 12, height = 5)

# define the cluster for each sample
table(merged_obj$patient, merged_obj$RNA_snn_res.0.2)
patients <- c("", "", "", "", "")
patients_clusters <- list(2, 3, 0, 1, c(5,4))
clusters_list <- setNames(patients_clusters, patients)
clusters_list

# delete cells in clusters from other samples
keep <- mapply(function(pat, clu) {
  (clu %in% clusters_list[[pat]]) | (clu %in% c(6,8,9))
}, merged_obj$patient, merged_obj$RNA_snn_res.0.2)

filtered_obj <- subset(merged_obj, cells = colnames(merged_obj)[keep])


# post-filtering plots
p5 <- DimPlot(filtered_obj, reduction = "umap.unintegrated", group.by = "patient", label = TRUE)
p6 <- DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = "RNA_snn_res.0.2", label = TRUE)
wrap_plots(p5,p6)
ggsave("", p5+p6, width = 12, height = 5)

# table of number of cells from each patient in each cluster
table_before <- table(
  merged_obj$patient,
  merged_obj$RNA_snn_res.0.2
)
table_before

table_after <- table(
  filtered_obj$patient,
  filtered_obj$RNA_snn_res.0.2
)
table_after

# visualize the removed cells
table_removed <- table_before - table_after
df_removed <- as.data.frame(table_removed)
colnames(df_removed) <- c("patient", "cluster", "cells_removed")
df_removed <- df_removed[df_removed$cells_removed > 0, ]
cells_per_patient <- table(merged_obj$patient)
percent_removed <- sweep(table_removed,1, # operar per files (pacients)
                         cells_per_patient, "/") * 100
df_percent <- as.data.frame(percent_removed)
colnames(df_percent) <- c("patient", "cluster", "percent_removed")

df_removed
df_percent
heatmap(as.matrix(table_removed))

p1 <- DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = "label_timepoint", shuffle=TRUE,
              cols =  c("#8FD3E0","#2E8FA3","#D4A6E8","#8A3FB0",
                        "#A8E6B5","#2FAE5A",
                        "#F0C48A","#D08A2E", "#F3A7B6","#D14B6A")) + 
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) + NoLegend()
p2 <-DimPlot(filtered_obj, reduction = "umap.unintegrated", group.by = "label_timepoint", shuffle=TRUE,
             cols =  c("#8FD3E0","#2E8FA3","#D4A6E8","#8A3FB0",
                       "#A8E6B5","#2FAE5A",
                       "#F0C48A","#D08A2E", "#F3A7B6","#D14B6A")) + 
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL)
wrap_plots(p1,p2, ncol = 2)


### 3. INTEGRATION 
# libraries
library(Seurat)
library(R.utils)
library(dplyr)
library(Matrix)
library(ggplot2)
library(scDblFinder)
library(SeuratObject)
library(patchwork)
library(clustree)


path_plots <- ""
path_objects <- ""

# name of the samples for the analysis
samples <- c("")

# merge objects
merged_obj <- merge(
  filtered_obj[[1]],
  y = filtered_obj[-1]
)

# check that the idents are saved correctly and join layers
all(sapply(filtered, ncol) == as.vector(table(merged_obj$orig.ident)))
merged_obj <- JoinLayers(merged_obj)

# add patient column to metadata and separate counts layers by patient
merged_obj@meta.data$patient <- sub("_.*", "", merged_obj@meta.data$orig.ident)
merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], f = merged_obj$patient)
cells <- merged_obj

# delete ig_genes
#ig_genes <- grep("^IG[HKL][VJCMDAG][0-9-]*", rownames(cells), value = TRUE)
#cells <- cells[!(rownames(cells) %in% ig_genes), ]

# data normalization
cells <- NormalizeData(cells)
cells <- FindVariableFeatures(cells)
cells <- ScaleData(cells)
cells <- RunPCA(cells)
ElbowPlot(cells, ndims = 50)
cells <- FindNeighbors(cells, dims = 1:30)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (r in resolutions) {
  cells <- FindClusters(cells, resolution = r)
}
clustree_uninteg <- clustree(cells@meta.data, prefix = "RNA_snn_res.")

cells <- RunUMAP(cells, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# pre-integration plots
p1 <- DimPlot(cells, reduction = "umap.unintegrated", label = TRUE)
p2 <- DimPlot(cells, reduction = "umap.unintegrated", group.by = "orig.ident", label = TRUE)
wrap_plots(p1,p2)

#ggsave(paste0(path_plots, "/clustree_uninteg.png"), clustree_uninteg)
#ggsave(paste0(path_plots, "/preintegration.png"),p2)


## INTEGRATION
cells <- IntegrateLayers(object = cells, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harmony", verbose = FALSE)
cells[["RNA"]] <- JoinLayers(cells[["RNA"]])

cells <- FindNeighbors(cells, reduction = "integrated.harmony", dims = 1:30)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (r in resolutions) {
  cells <- FindClusters(cells, resolution = r, cluster.name = paste0("Harmony_snn_res.",r))
}
clustree_harmony <- clustree(cells@meta.data, prefix = "Harmony_snn_res.")

cells <- RunUMAP(cells, dims = 1:30, reduction = "integrated.harmony", reduction.name = "umap.harmony")

#ggsave(paste0(path_plots, "/clustree_harmony.png"), clustree_harmony)

## post-integration plots
DimPlot(cells, reduction = "umap.harmony", group.by = "patient", label = TRUE)
harmony_plot <- DimPlot(cells, reduction = "umap.harmony", group.by = "orig.ident", label = TRUE)

plots <- DimPlot(cells, reduction = "umap.harmony", split.by = "patient", ncol = 3, label = TRUE)
wrap_plots(plots)

#ggsave(paste0(path_plots, "/harmony_integration.png"),harmony_plot)


## CELLTYPIST
# Define path to celltypist output
pathToCelltypist <- ""

# Specify celltypist annotation DIR
dir_celltypist <- ""

library(Matrix)
library(data.table)

rna_assay <- cells[["RNA"]]

count_layers <- grep("^counts", Layers(rna_assay), value = TRUE)

counts_list <- lapply(count_layers, function(layer) {
  GetAssayData(rna_assay, layer = layer)
})

all_genes <- unique(unlist(lapply(counts_list, rownames)))

counts_list_aligned <- lapply(counts_list, function(mat) {
  missing_genes <- setdiff(all_genes, rownames(mat))
  if(length(missing_genes) > 0) {
    mat_extra <- Matrix(0, nrow = length(missing_genes), ncol = ncol(mat),
                        dimnames = list(missing_genes, colnames(mat)))
    mat <- rbind(mat, mat_extra)
  }
  mat[all_genes, , drop = FALSE]
})

all_counts <- do.call(cbind, counts_list_aligned) 
all_counts_t <- t(all_counts)
dim(all_counts_t)

writeMM(all_counts_t, paste0(pathToCelltypist, "/celltypist_input.mtx"))

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
colnames(celltypist_low)[2] <- "celltypist_max_low_90"
probs <- fread(file.path(pathToCelltypist, "/model_low/probability_matrix.csv"), data.table = F)
colnames(probs)[2:ncol(probs)] <- paste0("celltypist_", gsub(" ", ".", colnames(probs)[2:ncol(probs)]))
celltypist_low$celltypist_max_prob_low_90 <- apply(probs[,2:ncol(probs)], 1, max)
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
colnames(celltypist_high)[2] <- "celltypist_max_high_90"
probs <- fread(file.path(pathToCelltypist, "/model_high/probability_matrix.csv"), data.table = F)
colnames(probs)[2:ncol(probs)] <- paste0("celltypist_", gsub(" ", ".", colnames(probs)[2:ncol(probs)]))
celltypist_high$celltypist_max_prob_high_90 <- apply(probs[,2:ncol(probs)], 1, max)
celltypist_high <- merge.data.frame(celltypist_high, df, by.x = "V1", by.y = "ID", all = T)
celltypist_high <- merge.data.frame(celltypist_high, probs, by = "V1", all = T)
cells@meta.data <- cbind(cells@meta.data, celltypist_high[match(rownames(cells@meta.data), celltypist_high$V1), 2:6])

# plots 
p8 <- DimPlot(cells, reduction = "umap.harmony", group.by = "celltypist_max_high_90")
p9 <- DimPlot(cells, reduction = "umap.cca", group.by = "celltypist_max_high_90")
wrap_plots(p8,p9)
p10 <- DimPlot(cells, reduction = "umap.harmony", group.by = "celltypist_max_high_90", split.by = "patient", ncol=3)
p11 <- DimPlot(cells, reduction = "umap.cca", group.by = "celltypist_max_high_90", split.by = "patient", ncol=3)
wrap_plots(p10)
wrap_plots(p11)

#ggsave(paste0(path_plots, "/celltypist_max_high_harmony.png"), p8)
#ggsave(paste0(path_plots, "/celltypist_max_high_cca.png"), p9)
#ggsave(paste0(path_plots, "/celltypist_max_high_harmony_split.png"), p10)
#ggsave(paste0(path_plots, "/celltypist_max_high_cca_split.png"), p11)


## FILTER B-CELLS
DimPlot(cells, group.by = "RNA_snn_res.0.2", label = TRUE)
DimPlot(cells, group.by = "patient")

b_cells <- subset(cells, subset = !(RNA_snn_res.0.2 %in% c(8,9,10,11)))

DimPlot(b_cells, reduction = "umap.unintegrated", group.by = "patient")


#ggsave(paste0(path_plots, "/bcells_celltypist_max_high_harmony.png"), p12)
#ggsave(paste0(path_plots, "/bcells_celltypist_max_high_cca.png"), p13)


## INTEGRATION B-CELLS
DimPlot(b_cells, reduction = "umap.harmony", group.by = "patient")

# data normalization
b_cells <- NormalizeData(b_cells)
b_cells <- FindVariableFeatures(b_cells)
b_cells <- ScaleData(b_cells)
b_cells <- RunPCA(b_cells)
ElbowPlot(b_cells, ndims = 50)
b_cells <- FindNeighbors(b_cells, dims = 1:30)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (r in resolutions) {
  b_cells <- FindClusters(b_cells, resolution = r)
}
clustree_uninteg <- clustree(b_cells@meta.data, prefix = "RNA_snn_res.")
clustree_uninteg

b_cells <- RunUMAP(b_cells, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# pre-integration plots
p1 <- DimPlot(b_cells, reduction = "umap.unintegrated", label = TRUE)
p2 <- DimPlot(b_cells, reduction = "umap.unintegrated", group.by = "orig.ident", label = TRUE)
wrap_plots(p1,p2)

#ggsave(paste0(path_plots, "/clustree_uninteg.png"), clustree_uninteg)
#ggsave(paste0(path_plots, "/preintegration.png"),p2)

# integration
b_cells[["RNA"]] <- split(b_cells[["RNA"]], f = b_cells$patient)
b_cells <- IntegrateLayers(object = b_cells, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harmony", verbose = FALSE)
b_cells[["RNA"]] <- JoinLayers(b_cells[["RNA"]])

b_cells <- FindNeighbors(b_cells, reduction = "integrated.harmony", dims = 1:30)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
for (r in resolutions) {
  b_cells <- FindClusters(b_cells, resolution = r, cluster.name = paste0("Harmony_snn_res.",r))
}
clustree_harmony <- clustree(b_cells@meta.data, prefix = "Harmony_snn_res.")
clustree_harmony

b_cells <- RunUMAP(b_cells, dims = 1:30, reduction = "integrated.harmony", reduction.name = "umap.harmony")


## post-integration plots
DimPlot(b_cells, reduction = "umap.harmony", group.by = "patient", label = TRUE)
harmony_plot <- DimPlot(b_cells, reduction = "umap.harmony", group.by = "orig.ident", label = TRUE)

plots <- DimPlot(b_cells, reduction = "umap.harmony", split.by = "patient", ncol = 3, label = TRUE)
wrap_plots(plots)

#ggsave(paste0(path_plots, "/harmony_integration.png"),harmony_plot)


### Functions for pathway enrichment analysis
plot_gsea <- function(data,
                      p_cutoff = 0.05,
                      top_n = 20) {
  
  df <- data@result
  
  # 🔥 FILTRATGE NET
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

### 4. ANALYSIS

# libraries
library(Seurat)
library(R.utils)
library(dplyr)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(scDblFinder)
library(SeuratObject)
library(patchwork)
library(clustree)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GeneNMF)
library(msigdbr)
library(openxlsx)
library(ggrepel)


plot1 <- DimPlot(ig_cells, group.by = "label_timepoint", shuffle = TRUE,
                 cols = c("#8FD3E0","#2E8FA3","#D4A6E8","#8A3FB0",
                          "#A8E6B5","#2FAE5A",
                          "#F0C48A","#D08A2E", "#F3A7B6","#D14B6A"),pt.size = 0.2) + 
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL)

plot2 <- FeaturePlot(ig_cells,features = "EBF1") +  labs(x = "UMAP 1", y = "UMAP 2")
plot3 <- FeaturePlot(ig_cells,features = "IGHM") +  labs(x = "UMAP 1", y = "UMAP 2")
plot4 <- VlnPlot(cells, group.by = "label_timepoint",features = "BTK", pt.size= 0, 
                 cols = c("#8FD3E0","#2E8FA3","#D4A6E8","#8A3FB0",
                          "#A8E6B5","#2FAE5A",
                          "#F0C48A","#D08A2E", "#F3A7B6","#D14B6A"))+ NoLegend()

wrap_plots(plot1, plot2, plot3,plot4, ncol=2)


## MARKERS T1 vs T2
markers.1344 <- FindMarkers(cells, ident.1 = "1344_T2", ident.2 = "1344_T1")
res.1344 <- run_GSEA_pipeline(markers.1344)

res.1344$plots$GO
res.1344$plots$KEGG
res.1344$plots$Hallmark
res.1344$plots$Reactome


markers.2913 <- FindMarkers(cells, ident.1 = "2913_T2", ident.2 = "2913_T1")
res.2913 <- run_GSEA_pipeline(markers.2913)

res.2913$plots$GO
res.2913$plots$KEGG
res.2913$plots$Hallmark
res.2913$plots$Reactome


markers.V1 <- FindMarkers(cells, ident.1 = "V1_T2", ident.2 = "V1_T1")
res.V1 <- run_GSEA_pipeline(markers.V1)

res.V1$plots$GO
res.V1$plots$KEGG
res.V1$plots$Hallmark
res.V1$plots$Reactome


markers.V2 <- FindMarkers(cells, ident.1 = "V2_T2", ident.2 = "V2_T1")
res.V2 <- run_GSEA_pipeline(markers.V2)

res.V2$plots$GO
res.V2$plots$KEGG
res.V2$plots$Hallmark
res.V2$plots$Reactome

markers.M1 <- FindMarkers(cells, ident.1 = "M1_T2", ident.2 = "M1_T1")
res.M1 <- run_GSEA_pipeline(markers.M1)

res.M1$plots$GO
res.M1$plots$KEGG
res.M1$plots$Hallmark
res.M1$plots$Reactome


res_list <- list(
  P2 = res.2913,
  P3 = res.M1,
  P4 = res.V1,
  P5 = res.V2
)

hallmark <- plot_multi_patient_GSEA(
  res_list,
  analysis = "Hallmark",
  top_n = 5
)
hallmark$dotplot


df <- data.frame(hallmark$results)
df_path <- df %>% dplyr::filter(Description == "HALLMARK_TNFA_SIGNALING_VIA_NFKB")

gene_list <- lapply(split(df_path, df_path$patient), function(x) {
  unique(unlist(strsplit(x$core_enrichment, "/")))
})

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

patients <- names(gene_list)

mat <- outer(
  patients,
  patients,
  Vectorize(function(p1, p2) {
    jaccard(gene_list[[p1]], gene_list[[p2]])
  })
)

rownames(mat) <- patients
colnames(mat) <- patients
mat_clean <- mat
mat_clean[is.na(mat_clean)] <- 0
pheatmap(
  mat_clean,
  color = colorRampPalette(c("white", "#BE4158"))(100),
  main = ("HALLMARK_TNFA_SIGNALING_VIA_NFKB genes similarity"),
  clustering_distance_rows = as.dist(1 - mat_clean),
  clustering_distance_cols = as.dist(1 - mat_clean)
)


df_path <- df %>% dplyr::filter(Description == "HALLMARK_UV_RESPONSE_UP")

gene_list <- lapply(split(df_path, df_path$patient), function(x) {
  unique(unlist(strsplit(x$core_enrichment, "/")))
})

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

patients <- names(gene_list)

mat <- outer(
  patients,
  patients,
  Vectorize(function(p1, p2) {
    jaccard(gene_list[[p1]], gene_list[[p2]])
  })
)

rownames(mat) <- patients
colnames(mat) <- patients
mat_clean <- mat
mat_clean[is.na(mat_clean)] <- 0
pheatmap(
  mat_clean,
  color = colorRampPalette(c("white", "#BE4158"))(100),
  main = ("HALLMARK_UV_RESPONSE_UP genes similarity"),
  clustering_distance_rows = as.dist(1 - mat_clean),
  clustering_distance_cols = as.dist(1 - mat_clean)
)




## save markers list
filt.1344$patient <- "1344"
filt.2913$patient <- "2913"
filt.M1$patient <- "M1"
filt.V1$patient <- "V1"
filt.V2$patient <- "V2"

filt.1344$gene <- rownames(filt.1344)
filt.2913$gene <- rownames(filt.2913)
filt.M1$gene   <- rownames(filt.M1)
filt.V1$gene   <- rownames(filt.V1)
filt.V2$gene   <- rownames(filt.V2)

all_markers <- rbind(filt.1344, filt.2913, filt.M1, filt.V1, filt.V2)
rownames(all_markers) <- NULL

write.csv(all_markers, "/all_markers_by_patient.csv", row.names = TRUE)


## COMMON GENES
filter.markers <- function(df) {
  df[df$p_val_adj < 0.05 &
       abs(df$avg_log2FC) > 0.5 &
       (df$pct.1 > 0.1 | df$pct.2 > 0.1), ]
}

filt.2913 <- filter.markers(markers.2913)
filt.M1 <- filter.markers(markers.M1)
filt.V1 <- filter.markers(markers.V1)
filt.V2 <- filter.markers(markers.V2)

common_genes <- Reduce(intersect, list(rownames(filt.2913), 
                                       rownames(filt.M1), rownames(filt.V1), rownames(filt.V2)))

s
df2_c <- filt.2913[common_genes, ]
df3_c <- filt.M1[common_genes, ]
df4_c <- filt.V1[common_genes, ]
df5_c <- filt.V2[common_genes, ]

mat <- cbind(
  df2_c$avg_log2FC,
  df3_c$avg_log2FC,
  df4_c$avg_log2FC,
  df5_c$avg_log2FC
)

rownames(mat) <- common_genes
colnames(mat) <- c("P2", "P3", "P4", "P5")

same_direction <- apply(mat, 1, function(x) {
  all(x > 0) | all(x < 0)
})

mat_filtered <- mat[same_direction, ]

genes_same_dir <- rownames(mat_filtered)
df1_c <- filt.1344[genes_same_dir, ]

stress <- c("FOSB", "PPP1R15A", "DDIT3", "TSPYL2")
signaling_up <- c("DUSP1", "NFKBIA", "ARRDC3", "IRF1", "HSH2D", "SAT1")
translation <- c("RRP12", "DDX39A", "EIF4A2", "EIF5", "NCL", "CDK11B", "DDX5", "BRD2", "POLR2A")
signaling_down <- c("PLCL2", "RASA1", "MIR155HG")
others <- c("DCLK2", "MAN1A1", "SNX10", "RHEX")
plus <- c("CXCR4", "CD5", "CD27")

gene_groups <- c(
  setNames(rep("Stress", length(stress)), stress),
  setNames(rep("Signaling up", length(signaling_up)), signaling_up),
  setNames(rep("Translation", length(translation)), translation),
  setNames(rep("Signaling down", length(signaling_down)), signaling_down),
  setNames(rep("Others", length(others)), others)
)

# Assignar grup a cada gen del heatmap
row_group <- gene_groups[rownames(mat_filtered)]

# Si hi ha gens sense grup assignat:
row_group[is.na(row_group)] <- "Unclassified"

row_group <- factor(
  row_group,
  levels = c(
    "Stress",
    "Signaling up",
    "Translation",
    "Signaling down",
    "Others",
    "Unclassified"
  )
)

library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(
  c(min(mat_filtered, na.rm = TRUE), 0, max(mat_filtered, na.rm = TRUE)),
  c("#4C6A92", "grey90", "#BE4158")
)

Heatmap(
  mat_filtered,
  name = "avg_log2FC",
  show_row_names = TRUE,
  col = col_fun,
  
  # Separar els gens per grup
  row_split = row_group,
  
  # Opcional
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  
  row_title_rot = 0
)



### ANALYSIS OF THE PATWHAYS PER CLUSTER

run_seurat_pipeline <- function(cells) {
  cells[["RNA"]] <- split(cells[["RNA"]], f = cells$timepoint)
  cells <- NormalizeData(cells)
  cells <- FindVariableFeatures(cells)
  cells <- ScaleData(cells)
  cells <- RunPCA(cells)
  cells <- FindNeighbors(cells, dims = 1:30)
  
  resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
  
  for (r in resolutions) {
    cells <- FindClusters(cells, resolution = r)
  }
  
  clustree_uninteg <- clustree::clustree(
    cells@meta.data,
    prefix = "RNA_snn_res."
  )
  
  cells <- RunUMAP(
    cells,
    dims = 1:30,
    reduction = "pca",
    reduction.name = "umap.unintegrated"
  )
  
  cells[["RNA"]] <- JoinLayers(cells[["RNA"]])
  
  return(list(
    seurat_object = cells,
    clustree_unintegrated = clustree_uninteg
  ))
  
}

p_1344 <- subset(cells, subset = patient == "1344")
p_2913 <- subset(cells, subset = patient == "2913")
p_V1 <- subset(cells, subset = patient == "V1")
p_V2 <- subset(cells, subset = patient == "V2")
p_M1 <- subset(cells, subset = patient == "M1")

res.1 <- run_seurat_pipeline(p_1344)
res.2 <- run_seurat_pipeline(p_2913)
res.3 <- run_seurat_pipeline(p_V1)
res.4 <- run_seurat_pipeline(p_V2)
res.5 <- run_seurat_pipeline(p_M1)

p_1344 <- res.1$seurat_object
p_2913 <- res.2$seurat_object
p_V1 <- res.3$seurat_object
p_V2 <- res.4$seurat_object
p_M1 <- res.5$seurat_object


# -------

DimPlot(p_1344, reduction = "umap.unintegrated", group.by = "timepoint", label = TRUE)
DimPlot(p_1344, reduction = "umap.unintegrated", group.by = "RNA_snn_res.0.6", label = TRUE)

table(p_1344$timepoint, p_1344$RNA_snn_res.0.6)
barplot(prop.table(table(p_1344$timepoint, p_1344$RNA_snn_res.0.6), margin = 2))

Idents(p_1344) <- p_1344@meta.data$RNA_snn_res.0.6
mark.timepoint.1344 <- FindAllMarkers(p_1344)
res.timepoint.1344 <- run_GSEA_pipeline(mark.timepoint.1344)

res.timepoint.1344$plots$Hallmark$cl_0


# -------

DimPlot(p_2913, reduction = "umap.unintegrated", group.by = "timepoint", label = TRUE)
DimPlot(p_2913, reduction = "umap.unintegrated", group.by = "RNA_snn_res.0.6", label = TRUE)

table(p_2913$timepoint, p_2913$RNA_snn_res.0.6)
barplot(prop.table(table(p_2913$timepoint, p_2913$RNA_snn_res.0.6), margin = 2))
# el cluster 0,3 són T2, 5 mig mig

Idents(p_2913) <- p_2913@meta.data$RNA_snn_res.0.6
mark.timepoint.2913 <- FindAllMarkers(p_2913)
res.timepoint.2913 <- run_GSEA_pipeline(mark.timepoint.2913)

res.timepoint.2913$plots$GO$cl_0


# -------

DimPlot(p_V1, reduction = "umap.unintegrated", group.by = "timepoint", label = TRUE)
DimPlot(p_V1, reduction = "umap.unintegrated", group.by = "RNA_snn_res.0.8", label = TRUE)

table(p_V1$timepoint, p_V1$RNA_snn_res.0.8)
barplot(prop.table(table(p_V1$timepoint, p_V1$RNA_snn_res.0.8), margin = 2))
# el cluster 2,5,9 són T2

Idents(p_V1) <- p_V1@meta.data$RNA_snn_res.0.8
mark.timepoint.V1 <- FindAllMarkers(p_V1)
res.timepoint.V1 <- run_GSEA_pipeline(mark.timepoint.V1)

res.timepoint.V1$plots$Hallmark$cl_2


# -------

DimPlot(p_V2, reduction = "umap.unintegrated", group.by = "timepoint", label = TRUE)
DimPlot(p_V2, reduction = "umap.unintegrated", group.by = "RNA_snn_res.0.4", label = TRUE)

table(p_V2$timepoint, p_V2$RNA_snn_res.0.4)
barplot(prop.table(table(p_V2$timepoint, p_V2$RNA_snn_res.0.4), margin = 2))

Idents(p_V2) <- p_V2@meta.data$RNA_snn_res.0.4
mark.timepoint.V2 <- FindAllMarkers(p_V2)
res.timepoint.V2 <- run_GSEA_pipeline(mark.timepoint.V2)

res.timepoint.V2$plots$Hallmark$cl_1

# -------

DimPlot(p_M1, reduction = "umap.unintegrated", group.by = "timepoint", label = TRUE)
DimPlot(p_M1, reduction = "umap.unintegrated", group.by = "RNA_snn_res.0.6", label = TRUE)

barplot(prop.table(table(p_M1$timepoint, p_M1$RNA_snn_res.0.6), margin = 2))

Idents(p_M1) <- p_M1@meta.data$RNA_snn_res.0.6
mark.timepoint.M1 <- FindAllMarkers(p_M1)
res.timepoint.M1 <- run_GSEA_pipeline(mark.timepoint.M1)

res.timepoint.M1$plots$Hallmark$cl_1

# save results
patients <- c("1344", "2913", "V1", "V2", "M1")
for (i in patients) {
  input <- get(paste0("res.timepoint.", i))
  output <- paste0("table_T1T2_", i, ".xlsx")
  all_tables <- list()
  for (db in c("GO", "KEGG", "Hallmark")) {
    clusters <- names(input[[db]])
    for (cl in clusters) {
      obj <- input[[db]][[cl]]
      if (is.null(obj)) next
      df <- obj@result
      if (nrow(df) == 0) next
      df$cluster <- cl
      df$database <- db
      all_tables[[length(all_tables) + 1]] <- df
    }
  }
  table <- do.call(rbind, all_tables)
  table <- na.omit(table)
  table <- table[table$p.adjust <= 0.05, ]
  rownames(table) <- NULL
  write.xlsx(
    table,
    file = paste0(
      "/Users/paula/Desktop/tables_pathways_parse/",
      output))}



# EXTRACT SIGNATURE AND PLOT

hallmark_pathways <- msigdbr(
  species = "Homo sapiens",
  category = "H")

tnfa_genes <- hallmark_pathways %>%
  dplyr::filter(gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB") %>%
  dplyr::pull(gene_symbol)

p_1344 <- AddModuleScore(p_1344, features = list(tnfa_genes), name = "NF-KB signature")
p_2913 <- AddModuleScore(p_2913, features = list(tnfa_genes), name = "NF-KB signature")
p_V1 <- AddModuleScore(p_V1, features = list(tnfa_genes), name = "NF-KB signature")
p_V2 <- AddModuleScore(p_V2, features = list(tnfa_genes), name = "NF-KB signature")
p_M1 <- AddModuleScore(p_M1, features = list(tnfa_genes), name = "NF-KB signature")
p_M2 <- AddModuleScore(p_M2, features = list(tnfa_genes), name = "NF-KB signature")


## violin plots for the signature by clusters

plot_vln_by_timepoint_colored <- function(seurat_obj, label,
                                          feature,
                                          cluster_col = "RNA_snn_res.0.6",
                                          timepoint_col = "timepoint",
                                          tp_interest = "T1") {
  
  library(dplyr)
  library(ggplot2)
  
  patient <- deparse(substitute(seurat_obj))
  patient <- gsub("^p_", "", patient)
  
  meta <- seurat_obj@meta.data
  
  prop <- prop.table(
    table(meta[[timepoint_col]], meta[[cluster_col]]),
    margin = 2
  )
  
  prop_T1 <- as.numeric(prop[tp_interest, ])
  names(prop_T1) <- colnames(prop)
  
  df_clusters <- data.frame(
    cluster = names(prop_T1),
    prop_T1 = prop_T1
  )
  
  df_clusters$tp_group <- ifelse(df_clusters$prop_T1 > 0.6, "T1",
                                 ifelse(df_clusters$prop_T1 < 0.4, "T2", "T1/T2"))
  
  cluster_order <- df_clusters %>%
    arrange(desc(prop_T1)) %>%
    pull(cluster)
  
  cluster_colors <- c(
    "T1" = "#A6BDDB",
    "T1/T2" = "grey90",
    "T2" = "#BE4158")
  
  plot_df <- meta %>%
    dplyr::select(all_of(c(cluster_col, feature))) %>%
    dplyr::rename(
      cluster = all_of(cluster_col),
      value = all_of(feature)
    )
  
  
  plot_df$cluster <- as.character(plot_df$cluster)
  plot_df$cluster <- factor(plot_df$cluster, levels = cluster_order)
  
  df_clusters$cluster <- as.character(df_clusters$cluster)
  df_clusters$cluster <- factor(df_clusters$cluster, levels = cluster_order)
  
  plot_df <- plot_df %>%
    left_join(df_clusters[, c("cluster", "tp_group")], by = "cluster")
  
  median_t1 <- median(plot_df$value[plot_df$tp_group == "T1"],na.rm = TRUE)
  
  p <- ggplot(plot_df, aes(x = cluster, y = value, fill = tp_group)) +
    geom_violin(trim = FALSE, scale = "width") +
    scale_fill_manual(values = cluster_colors, name = "Timepoint group") +
    geom_hline(yintercept = median_t1, linetype = "dashed", color = "black") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    labs(
      x = "Cluster",
      y = "NF-KB signature score",
      title = paste0("Patient ", label)
    )
  
  return(list(
    plot = p,
    cluster_table = df_clusters,
    patient = patient
  ))
}



vln.1344 <- plot_vln_by_timepoint_colored(
  p_1344, label = "1",
  feature = "NF-KB signature1")
vln.2913 <- plot_vln_by_timepoint_colored(
  p_2913, label = "2",
  feature = "NF-KB signature1")
vln.V1 <- plot_vln_by_timepoint_colored(
  p_V1,  label = "4",
  cluster_col = "RNA_snn_res.0.8",
  feature = "NF-KB signature1")
vln.V2 <- plot_vln_by_timepoint_colored(
  p_V2,  label = "5",
  cluster_col = "RNA_snn_res.0.4",
  feature = "NF-KB signature1")
vln.M1 <- plot_vln_by_timepoint_colored(
  p_M1,  label = "3",
  feature = "NF-KB signature1")


plot1 <- vln.1344$plot + NoLegend()
plot2 <- vln.2913$plot
plot3 <- vln.M1$plot+ NoLegend()
plot4 <- vln.V1$plot+ NoLegend()
plot5 <- vln.V2$plot+ NoLegend()


combined_plot <- (plot1 | plot2 | plot3) /
  (plot4   | plot5)
combined_plot +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")



# ridgeplot 
order_1344 <- c("1", "2", "3", "0","5","4")
p_1344$RNA_snn_res.0.6 <- factor(
  p_1344$RNA_snn_res.0.6,
  levels = order_1344)
Idents(p_1344) <- p_1344@meta.data$RNA_snn_res.0.6
plot6 <- RidgePlot(p_1344, features = "NF-KB signature1", 
                   cols = c("grey90", "grey90", "#BE4158","#BE4158","#BE4158","#BE4158")) +
  ggtitle("NF-KB signature - 1344")


order_2913 <- c("4", "1", "2","5","0", "3")
p_2913$RNA_snn_res.0.6 <- factor(
  p_2913$RNA_snn_res.0.6,
  levels = order_2913)
Idents(p_2913) <- p_2913@meta.data$RNA_snn_res.0.6
plot7 <- RidgePlot(p_2913, features = "NF-KB signature1", 
                   cols = c("#4C6A92", "#4C6A92", "#4C6A92","grey90","#BE4158","#BE4158")) +
  ggtitle("NF-KB signature - 2913")


order_V1 <- c("6", "8", "3", "4","0", "7", "1", "9", "5", "2")
p_V1$RNA_snn_res.0.8 <- factor(
  p_V1$RNA_snn_res.0.8,
  levels = order_V1)
Idents(p_V1) <- p_V1@meta.data$RNA_snn_res.0.8
plot8 <- RidgePlot(p_V1, features = "NF-KB signature1", 
                   cols = c("#4C6A92", "#4C6A92", "#4C6A92","#4C6A92","#4C6A92", "#4C6A92", "#4C6A92","#BE4158","#BE4158","#BE4158")) +
  ggtitle("NF-KB signature - V1")


order_V2 <- c("2", "5", "6", "7", "9", "1", "3", "4", "0", "8", "10")
p_V2$RNA_snn_res.0.4 <- factor(
  p_V2$RNA_snn_res.0.4,
  levels = order_V2)
Idents(p_V2) <- p_V2@meta.data$RNA_snn_res.0.4
plot9 <- RidgePlot(p_V2, features = "NF-KB signature1",
                   cols = c("#4C6A92", "#4C6A92","#4C6A92", "#4C6A92", "#4C6A92","#4C6A92","#4C6A92","#BE4158","#BE4158","#BE4158","#BE4158")) +
  ggtitle("NF-KB signature - V2")


order_M1 <- c("0", "3", "2", "6", "1", "4", "5")
p_M1$RNA_snn_res.0.6 <- factor(
  p_M1$RNA_snn_res.0.6,
  levels = order_M1)
Idents(p_M1) <- p_M1@meta.data$RNA_snn_res.0.6
plot10 <- RidgePlot(p_M1, features = "NF-KB signature1", 
                    cols = c("#4C6A92", "#4C6A92", "#4C6A92","#4C6A92","#BE4158","#BE4158","#BE4158")) +
  ggtitle("NF-KB signature - M1")

combined_plot <- (plot6 | plot7 | plot8) /
  (plot9   | plot10)
combined_plot +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


### PSEUDOBULK
library(DESeq2)
library(ggrepel)

pseudo_cells <- AggregateExpression(cells, assays = "RNA", return.seurat = TRUE, group.by = "label_timepoint")
length(Cells(pseudo_cells))

counts_matrix <- pseudo_cells$RNA$counts
counts_matrix <- as.matrix(counts_matrix)

sample_info <- data.frame(
  sample_id = colnames(counts_matrix),
  Sample = c("P1_T1", "P1_T2", "P2_T1", "P2_T2", "P3_T1", "P3_T2","P4_T1", "P4_T2", "P5_T1", "P5_T2")
)
rownames(sample_info) <- sample_info$sample_id

dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = sample_info,
                              design = ~ Sample)

# Transformació VST (crucial per a PCAs nets)
vsd <- vst(dds, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = "Sample", returnData = TRUE)

ggplot(pca_data, aes(x = PC1, y = PC2, color = Sample, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(
    size = 4,                  # Mida de la lletra
    box.padding = 0.5,           # Espai al voltant de cada etiqueta
    point.padding = 0.3,         # Espai respecte al punt del PCA
    max.overlaps = Inf,          # Força a que s'enllacin TOTES les etiquetes, encara que n'hi hagin moltes
    segment.color = 'grey50',    # Color de la línia connector
    segment.size = 0.2           # Gruix de la línia connector
  ) +
  scale_color_manual(
    values = c(
      "#8FD3E0","#2E8FA3",
      "#D4A6E8","#8A3FB0",
      "#A8E6B5","#2FAE5A",
      "#F0C48A","#D08A2E",
      "#F3A7B6","#D14B6A"
    )) +
  theme_minimal() +
  labs(title = "PCA of pseudobulk profiles",
       x = paste0("PC1: ", round(attr(pca_data, "percentVar")[1] * 100), "%"),
       y = paste0("PC2: ", round(attr(pca_data, "percentVar")[2] * 100), "%"))



### correlation CXCR4/MIR155HG with NF-KB signature
patient <- p_2913

meta <- patient@meta.data %>%
  dplyr::filter(timepoint %in% c("T1", "T2"))

df_cluster <- meta %>%
  dplyr::group_by(RNA_snn_res.0.6) %>%
  dplyr::summarise(
    T2_fraction = mean(timepoint == "T2")
  ) %>%
  dplyr::rename(cluster = RNA_snn_res.0.6)

MIR155HG_expr <- AverageExpression(
  patient,
  features = "MIR155HG",
  group.by = "RNA_snn_res.0.6",
  assays = "RNA"
)$RNA

df_expr <- data.frame(
  cluster = gsub("^g", "", colnames(MIR155HG_expr)),
  MIR155HG_expression = MIR155HG_expr[1,]
)

nfkb_cluster <- patient@meta.data %>%
  dplyr::group_by(RNA_snn_res.0.6) %>%
  dplyr::summarise(
    NFkB_signature = mean(`NF-KB signature1`, na.rm = TRUE)
  ) %>%
  dplyr::rename(cluster = RNA_snn_res.0.6)


df <- nfkb_cluster %>%
  left_join(df_expr, by = "cluster") %>%
  left_join(df_cluster, by = "cluster")

ggplot(df, aes(x = NFkB_signature, y = MIR155HG_expression)) +
  geom_point(aes(color = T2_fraction), size = 4, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, color = "#BE4158") +
  ggrepel::geom_text_repel(aes(label = cluster), size = 4) +
  theme_classic(base_size = 14) +
  labs(
    x = "NF-κB signature (mean per cluster)",
    y = "MIR155HG expression (mean per cluster)",
    color = "T2_fraction",
    title = "Association between NF-κB activity and MIR155HG across clusters"
  )


### identify markers by clusters
filter.markers <- function(df) {
  df[df$p_val_adj < 0.05 &
       abs(df$avg_log2FC) > 0.5 &
       (df$pct.1 > 0.1 | df$pct.2 > 0.1), ]
}

filt.1344 <- filter.markers(mark.timepoint.1344)
filt.2913 <- filter.markers(mark.timepoint.2913)
filt.V1 <- filter.markers(mark.timepoint.V1)
filt.V2 <- filter.markers(mark.timepoint.V2)
filt.M1 <- filter.markers(mark.timepoint.M1)


### figures clusters
library(patchwork)

# ----- 2913
DimPlot(p_2913, group.by = "RNA_snn_res.0.6", cols = c(
  "#9ECAE1",
  "#6BAED6",
  "#3182BD",
  "#8C1D2C",
  "#D14B6A",
  "#F3A7B6" 
), label = TRUE) +  labs(x = "UMAP 1", y = "UMAP 2", title = "Patient 2") +
  guides(color = guide_legend(title = "Cluster identity",override.aes = list(size = 4)))



# ----- V1
DimPlot(p_V1, group.by = "RNA_snn_res.0.8", cols = c(
  "#E6F2FA",
  "#CFE7F5",
  "#9ECAE1",
  "#6BAED6",
  "#4A90C2",
  "#3182BD",
  "#1F5F9E",
  "#8C1D2C",
  "#D14B6A",
  "#F3A7B6" 
), label = TRUE) +  labs(x = "UMAP 1", y = "UMAP 2", title = "Patient 5")+
  guides(color = guide_legend(title = "Cluster identity",override.aes = list(size = 4)))


# ----- V2
DimPlot(p_V2, group.by = "RNA_snn_res.0.4", cols = c(
  "#E6F2FA",
  "#CFE7F5",
  "#9ECAE1",
  "#6BAED6",
  "#4A90C2",
  "#3182BD",
  "#1F5F9E",
  "#8C1D2C",
  "#D14B6A",
  "#E6788F",
  "#F3A7B6" 
), label=TRUE) +  labs(x = "UMAP 1", y = "UMAP 2", title = "Patient 6")+
  guides(color = guide_legend(title = "Cluster identity",override.aes = list(size = 4)))


# ----- M1
DimPlot(p_M1, group.by = "RNA_snn_res.0.6", cols = c(
  "#CFE7F5",
  "#9ECAE1",
  "#6BAED6",
  "#3182BD",
  "#8C1D2C",
  "#D14B6A",
  "#F3A7B6" 
), label = TRUE) +  labs(x = "UMAP 1", y = "UMAP 2", title = "Patient 3")+
  guides(color = guide_legend(title = "Cluster identity",override.aes = list(size = 4)))


# barplot clusters T1 T2
df <- as.data.frame(prop.table(
  table(p_V2$timepoint, p_V2$RNA_snn_res.0.4),
  margin = 2
))
colnames(df) <- c("timepoint", "cluster", "proportion")
ggplot(df, aes(x = cluster, y = proportion, fill = timepoint)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c(
    "T1" = "#A6BDDB",
    "T2" = "#BE4158"
  )) +theme_classic() +
  labs(x = "Cluster", y = "Cell proportion", fill = "Timepoint") +
  ggtitle("Patient 5") +
  theme(plot.title = element_text(
    size = 20,
    face = "bold",
    hjust = 0.5
  ))+
  theme(
    axis.title = element_text(size = 14),  # títols dels eixos
    axis.text = element_text(size = 14)                   # etiquetes dels eixos
  ) +theme(
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18)
  )



## plots AXIS des de taules de GSEA clusters 

col_fun <- scale_color_gradient2(
  low = "#4C6A92",
  high = "#BE4158",
  limits = range(-2,2),
  oob = scales::squish
)

size_fun <- scale_size(
  range = c(2, 10),
  limits = range(table %>%
                   group_by(sample, cluster, axis) %>%
                   summarise(n_pathways = n()) %>%
                   pull(n_pathways),
                 na.rm = TRUE)
)

dotplot_axis_clusters <- function(patient_id,
                                  df,
                                  clusters_T1 = NULL,
                                  clusters_T1T2 = NULL,
                                  clusters_T2 = NULL,
                                  split_groups = TRUE) {
  
  ordre_clusters <- c(clusters_T1, clusters_T1T2, clusters_T2)
  
  df_patient <- df %>%
    filter(label == patient_id) %>%
    group_by(cluster, axis) %>%
    summarise(
      n_pathways = n(),
      mean_NES = mean(NES, na.rm = TRUE),
      mean_padj = mean(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      cluster = factor(cluster, levels = ordre_clusters)
    )
  
  if (split_groups) {
    
    df_patient <- df_patient %>%
      mutate(
        group = case_when(
          as.character(cluster) %in% clusters_T1 ~ "T1",
          as.character(cluster) %in% clusters_T1T2 ~ "T1/T2",
          as.character(cluster) %in% clusters_T2 ~ "T2",
          TRUE ~ NA_character_
        ),
        group = factor(group, levels = c("T1", "T1/T2", "T2"))
      )
    
    p <- ggplot(df_patient, aes(x = cluster, y = axis)) +
      geom_point(aes(size = n_pathways, color = mean_NES)) +
      col_fun +
      size_fun +
      facet_grid(~group, scales = "free_x", space = "free_x")
    
  } else {
    
    p <- ggplot(df_patient, aes(x = cluster, y = axis)) +
      geom_point(aes(size = n_pathways, color = mean_NES)) +
      col_fun +
      size_fun
  }
  
  p +
    theme_bw() +
    labs(
      title = paste0("Patient ", gsub("^P", "", patient_id)),
      x = "",
      y = "Axis",
      size = "Nº pathways",
      color = "Mean NES"
    ) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(size = 14, face = "bold")
    )
}

library(readxl)
library(readxl)
library(dplyr)
library(ggplot2)

sheets <- c("1344", "2913", "V1", "V2", "M1")

all_tables <- list()
for (s in sheets) {
  df <- read_xlsx(
    "",
    sheet = s
  )
  df$sample <- s
  all_tables[[s]] <- df
}

axis_plot <- c("bcr", "stress", "immune", "cell cycle")

table <- bind_rows(all_tables)
table <- table[!is.na(table$axis), ]
table <- table %>% dplyr::filter(axis %in% axis_plot)
table <- table %>%
  mutate(label = recode(sample,
                        "1344" = "P1",
                        "2913" = "P2",
                        "M1" = "P3",
                        "V1" = "P4",
                        "V2" = "P5"))

table <- table %>%
  mutate(axis = recode(axis, "bcr" = "downstream BCR"))


plot_1344 <- dotplot_axis_clusters("P1", table, clusters_T1T2 = c("cl_1", "cl_2"), clusters_T2 = c( "cl_3", "cl_0","cl_5","cl_4"))
plot_2913 <- dotplot_axis_clusters("P2", table, clusters_T1 = c("cl_4","cl_1", "cl_2"), clusters_T1T2 = c("cl_5"), clusters_T2 = c("cl_0", "cl_3"))
plot_V1 <- dotplot_axis_clusters("P4", table, clusters_T1 = c("cl_6","cl_8", "cl_3", "cl_4","cl_0", "cl_7", "cl_1"), clusters_T2 = c("cl_9","cl_5","cl_2"))
plot_V2 <- dotplot_axis_clusters("P5", table, clusters_T1 = c("cl_2","cl_5", "cl_6", "cl_7", "cl_9", "cl_1", "cl_3"), clusters_T2 = c( "cl_4", "cl_0", "cl_8", "cl_10"))
plot_M1 <- dotplot_axis_clusters("P3", table, clusters_T1 = c("cl_0", "cl_3", "cl_2", "cl_6"), clusters_T2 = c("cl_1", "cl_4", "cl_5"))

combined_plot <- (plot_1344 | plot_2913 | plot_M1) /
  (plot_V1   | plot_V2)

ggsave("combined_plot.png",
       combined_plot + plot_layout(guides = "collect") &
         theme(legend.position = "right"),
       width = 16,
       height = 8,
       dpi = 300)
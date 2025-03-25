# WGCNA Tutorial


# pip line ----------------------------------------------------------------

# read_mt > normalization(FPKM recommended) > 

# Start -------------------------------------------------------------------
library(WGCNA)
library(dplyr)
library(textshape)
library(DESeq2)
library(magrittr)
library(ggplot2)
library(ComplexHeatmap)

#Setting string not as factor
options(stringsAsFactors = FALSE)
#Enable multithread
enableWGCNAThreads()


# Preparing the Data ------------------------------------------------------
setwd(dir = "C:/Users/SY/Desktop/Work_NOW/Project/AMLE")
#Reading the raw data (rows are the sample and columns the genes)
# expressiondata = read.csv("00 Raw/Readcount.csv", header = T)
# 
# expressiondata<- distinct(expressiondata,probe_ID, .keep_all = T) # 중복값제거
# 
# colnames(expressiondata) <- c("geneID", "DMSO1","DMSO2","DMSO3", "AMLE1","AMLE2","AMLE3")
# expressiondata<- column_to_rownames(expressiondata,loc = "geneID")
# 
# expressiondata[is.na(expressiondata)] <- 0
# 
# expressiondata <- as.matrix(expressiondata)
# 
# na_indices <- which(is.na(expressiondata), arr.ind = TRUE)
# 
# any_na <- any(is.na(expressiondata))
# if (any_na) {
#   print("Has NA in dataset.")
# } else {
#   print("Not NA in dataset.")
# }
# 
# head(expressiondata)
# write.csv(expressiondata,"./DMSO_APP_Readcount.csv")

# ------------------------------------------------------------------
expressiondata = read.csv("00 Raw/DMSO_APP_Readcount.csv", header = T, row.names = 1)

expressiondata <- round(expressiondata)
head(expressiondata)
metadata <-  data.frame(
  id = c(c("DMSO1","DMSO2","DMSO3","AMLE1","AMLE2","AMLE3")),
  treatment = c("control","control","control","AMLE","AMLE","AMLE")
  )
# filtered_counts <- expressiondata[rownames(expressiondata) %in% protein_coding_genes, ] #AMLE rnaseq fileter vactor
# expressiondata <- filtered_counts

all.equal(colnames(expressiondata), metadata$id) 


dds <- DESeqDataSetFromMatrix(
  countData = expressiondata,
  colData = metadata,
  design = ~1
)



# DESeq2 ------------------------------------------------------
# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
# function from the `DESEq2` R package
dds_norm <- vst(dds)

# Retrieve the normalized data from the `DESeqDataSet`
normalized_counts <- assay(dds_norm) %>%
  t() # Transpose this data


str(dds_norm)

#  WGCNA ------------------------------------------------------


sft <- pickSoftThreshold(normalized_counts,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed"
)

# This sft object has a lot of information, we will want to plot some of it to figure out 
# what our power soft-threshold should be. We have to first calculate a measure of the model fit, 
# the signed R2, and make that a new variable.
sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.90, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()


# otherplot ---------------------------------------------------------------

par(mfrow = c(1,2))
cex1 = 0.9

#Index the scale free topology adjust as a function of the power soft thresholding.
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],cex=cex1,col="red")

#This line corresponds to use a cut-off R² of h
abline(h=0.90,col="red")

#Connectivity mean as a function of soft power thresholding
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5],cex=cex1,col="red")

#This line corresponds to use a cut-off R² of h
abline(h=0.90,col="red")

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# 
# Run WGCNA ! # it well take the few time
# defult
# bwnet <- blockwiseModules(normalized_counts,
#                           maxBlockSize = 5000, # What size chunks (how many genes) the calculations should be run in
#                           TOMType = "signed", # topological overlap matrix
#                           power = 7, # soft threshold for network construction
#                           numericLabels = TRUE, # Let's use numbers instead of colors for module labels
#                           randomSeed = 1234, # there's some randomness associated with this calculation
#                           # so we should set a seed
)

# 
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------#2000개 이하로 잡으면 hub 찾는 STRING 분석가능
picked_power = 20
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(normalized_counts,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 2000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

cor <- temp_cor
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# save WGCNA results
# readr::write_rds(netwk,
#                  file = file.path("./01. Data/WGCNA with DEGseq2 normalization_USE_Protein_coding_gene/", "AMLE_wgcna_results_2_sft20.RDS")
# )
# -------------------------------------------------------------------------
# bwnet<- readRDS("./RDS/AMLE_wgcna_results.RDS")
# -------------------------------------------------------------------------
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

# netwk$colors[netwk$blockGenes[[1]]]
# table(netwk$colors)

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:6,]
# gene_id    colors
# 1 0610010K14Rik turquoise
# 2 0610030E20Rik      cyan
# 3 0610040J01Rik turquoise
# 4 1110004F10Rik turquoise
# 5 1110032F04Rik      blue
# 6 1110038F14Rik turquoise

# write_delim(module_df,
#             file = "gene_modules.txt",
#             delim = "\t")


# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(normalized_counts, mergedColors)$eigengenes

# MEs0 <-MEs0[,-14]

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

# Cont the gene numbers in each module
# summary(module_df)
# category_counts <- table(module_df[, 2])
# # Print the summary
# view(category_counts)


# -------------------------------------------------------------------------

###selacte ##  turquoise, blue, grey

# -------------------------------------------------------------------------

# pick out a few modules of interest here
modules_of_interest = c("blue", "turquoise")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
# expr_normalized<- t(normalized_counts)

expr_normalized[1:5,1:6]
#               DMSO1     DMSO2     DMSO3     AMLE1     AMLE2     AMLE3
# 0610010K14Rik  9.936115 10.011708  9.995117  9.827369  9.779074  9.653157
# 0610030E20Rik 11.131256 11.143314 11.011290 11.087152 11.134534 11.076914
# 0610040J01Rik  6.716757  6.795473  6.855045  6.306746  6.408461  6.595648
# 1110004F10Rik 11.244453 11.285536 11.283922 11.171566 11.172929 11.169228
# 1110032F04Rik  6.507804  6.589257  6.803953  7.469468  7.477599  7.391119
subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id))) +
  geom_line(aes(color = module),
            alpha = 0.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")



# -------------------------------------------------------------------------
# Creat the natwork # run-time 20min
# -------------------------------------------------------------------------

genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = expr_normalized[genes_of_interest$gene_id,]
expr_of_interest[1:5,1:5]

TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)

row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)


edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

head(edge_list)


write_delim(edge_list,
            file = "edgelist.tsv",
            delim = "\t")

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

#WGCNA 결과
module_eigengenes <- netwk$MEs

# Print out a preview
head(MEs0)
# equal cols?
all.equal(metadata$id, rownames(MEs0))

# Create the design matrix from the `treatment(group)` variable
des_mat <- model.matrix(~ metadata$treatment)

# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(MEs0), design = des_mat) 

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(MEs0)) %>%
  tibble::rownames_to_column("module")
#> Removing intercept from test coefficients


View(stats_df)


# -------------------------------------------------------------------------



# Let’s make plot of module ---------------------------------------------
module_df <- MEs0 %>%
  tibble::rownames_to_column("accession_code") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(metadata %>%
                      dplyr::select(id, treatment),
                    by = c("accession_code" = "id")
  ) %>%
  # treatment 순서를 지정
  mutate(treatment = factor(treatment, levels = c("control", "AMLE")))  # control 먼저, AMLE 나중

# > colnames(module_df)

ggplot(
  module_df,
  aes(
    x = treatment,
    y = MEmidnightblue,  
    color = treatment
  )
) +

  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()


# temp-------------------------------------------------------------------------
heatmap_data <- module_df %>%
  select(treatment, starts_with("ME")) %>%
  pivot_longer(
    cols = starts_with("ME"), 
    names_to = "module", 
    values_to = "expression"
  ) %>%
  mutate(treatment = factor(treatment, levels = c("control", "AMLE"))) %>%
  group_by(module) %>%
  mutate(avg_expression = mean(expression)) %>%
  ungroup() %>%
  arrange(desc(avg_expression)) %>%
  mutate(module = factor(module, levels = unique(module)))

stats_df <- na.omit(stats_df)

sorted_modules <- stats_df %>%
  arrange(desc(logFC)) %>%  
  pull(module)             

heatmap_data <- heatmap_data %>%
  mutate(module = factor(module, levels = sorted_modules))

heatmap_data<- na.omit(heatmap_data)



# 2. heatmap create
ggplot(
  heatmap_data,
  aes(
    x = treatment,  # col: treatment
    y = module,     # row: module
    fill = expression 
  )
) +
  geom_tile(color = "white") + 
  scale_fill_gradientn(
    colors = c("blue", "white", "red"),
    name = "Expression"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 10), 
    panel.grid = element_blank() 
  ) +
  labs(
    x = "Treatment",
    y = "Module",
    title = "Module Expression Heatmap"
  ) 

# zscoreing ----------------------------------------------
heatmap_data2 <- module_df %>%
  select(treatment, starts_with("ME")) %>%
  pivot_longer(
    cols = starts_with("ME"),
    names_to = "module",
    values_to = "expression"
  ) %>%
  mutate(treatment = factor(treatment, levels = c("control", "AMLE"))) %>%
  group_by(module) %>%
  mutate(
    module_mean = mean(expression),
    module_sd = sd(expression),    
    z_score = (expression - module_mean) / module_sd 
  ) %>%
  ungroup() %>%
  # Z-score 기준으로 정렬
  arrange(desc(z_score)) %>%
  mutate(module = factor(module, levels = unique(module)))

summary_df <- heatmap_data2 %>%
  group_by(treatment, module) %>%
  summarize(
    avg_expression = mean(expression, na.rm = TRUE), 
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = treatment,
    values_from = avg_expression 
  )




gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  dplyr::mutate(module = paste0("ME", module))

gene_module_key %>%
  dplyr::filter(module == "ME1") -> MODUL
# MODUL keys
# readr::write_tsv(gene_module_key,
#                  file = file.path("./01. Data/", "AMLE_wgcna_gene_to_module.tsv")
# )

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
### set heatmap funtion
make_module_heatmap <- function(module_name,
                                expression_mat = normalized_counts,
                                metadata_df = metadata,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes) {

  
  
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("id")
  
  
  col_annot_df <- metadata_df %>%
    dplyr::select(id, treatment) %>%
    dplyr::inner_join(module_eigengene, by = "id") %>%
    dplyr::arrange(treatment) %>%
    tibble::column_to_rownames("id")
  
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    treatment = col_annot_df$treatment,
    # Add annotation barplot
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    # Pick colors for each experimental group in treatment
    col = list(treatment = c("recovering" = "#f1a340", "acute illness" = "#998ec3"))
  )
  
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    dplyr::filter(rownames(.) %in% module_genes) %>%
    dplyr::select(rownames(col_annot_df)) %>%
    as.matrix()
  
  mod_mat <- mod_mat %>%
    t() %>%
    scale() %>%
    t()
  
  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = module_name,
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = FALSE,
                                     show_column_names = FALSE
  )
  return(heatmap)
}
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
mod_heatmap <- make_module_heatmap(module_name = "ME1")
mod_heatmap


cor_matrix <- module_df %>%
  select(starts_with("ME")) %>% 
  cor(method = "pearson") 

# Correlation Matrix
print(cor_matrix)


library(pheatmap)

# Correlation Matrix Heatmap
pheatmap(
  cor_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  main = "Module Correlation Heatmap"
)
#option1==================================================================================
cor_data <- melt(cor_matrix)



ggplot() +
  geom_tile(data = cor_data_lower, aes(x = Var1, y = Var2, fill = value), color = "black") +
  geom_tile(data = cor_data_upper, aes(x = Var1, y = Var2, fill = value), color = "black") + 
  geom_text(data = cor_data_upper, aes(x = Var1, y = Var2, label = round(value, 2)), size = 3) +
  scale_fill_gradientn(colors = c("blue", "white", "red"), limits = c(-1, 1), name = "Correlation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Module",
    y = "Module",
    title = "Module Correlation Heatmap with Numbers"
  )

# -------------------------------------------------------------------------


head(expressiondata)
library(edgeR)
dataset_norm <- cpm(expressiondata)
dataset_norm = log2(1+dataset_norm)
head(dataset_norm[,1:5])

plot(dataset_norm[,"DMSO1"],dataset_norm[,"AMLE1"])
#remove 0 counts
gsg <- goodSamplesGenes(datExpr = dataset_norm, verbose = 3)
print(paste("All ok?", gsg$allOK))

sprintf("Removing %d features", ncol(dataset_norm) - sum(gsg$goodGenes))
sprintf("Removing %d samples", nrow(dataset_norm) - sum(gsg$goodSamples))
dataset_norm = dataset_norm[gsg$goodSamples, gsg$goodGenes]

samplesTree <- hclust(d = dist(dataset_norm), method = "average")


plot(samplesTree, main = "Samples dendrogram", sub = "", xlab = "", cex.lab = 1, cex.axis = 1, cex.main = 1, cex = 0.5)
abline(h = 10, col = "red")


clust <- cutreeStatic(samplesTree, cutHeight = 15)
table(clust)



keep <- clust == 1
dataset_norm <- dataset_norm[keep,]
nFeatures <- ncol(dataset_norm)
nSamples <- nrow(dataset_norm)
sprintf("%d features, %d samples",nFeatures, nSamples)
# [1] "6 features, 24390 samples"

rownames(metadata) <- metadata$id
metadata$id <- NULL 

head(metadata)

# remove na
metadata$treatment[which(is.na(metadata$treatment))] = sample(na.omit(metadata$treatment),1)

#
traitColors <- labels2colors(metadata)
traitColors[,1] <- numbers2colors(metadata[,1]) # The first column is numerical
head(traitColors)

samplesTree <- hclust(d = dist(dataset_norm), method = "average")
plotDendroAndColors(dendro = samplesTree, colors = traitColors, groupLabels = names(metadata), 
                    main = "Patients dendogram and trait heatmap",cex.dendroLabels = 0.5, cex.colorLabels = 0.5) 


plot(dataset_norm[,1],dataset_norm[,4])
cor(dataset_norm[,1:6],dataset_norm[,1:6],method ="pearson")#spearman

cor_matrix <- cor(dataset_norm[, 1:6], dataset_norm[, 1:6], method = "pearson")


library(reshape2)
cor_data <- melt(cor_matrix)


library(ggplot2)
ggplot(cor_data, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") + 
  scale_fill_gradient(
    low = "white", high = "red", 
    limits = c(0.98, 1), 
    name = "Correlation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.text.y = element_text(size = 10) 
  ) +
  labs(
    x = "Samples",
    y = "Samples",
    title = "Heatmap of Pearson Correlation"
  )  +  geom_text(aes(Var1, Var2, label = round(value, 3)), color = "black", size = 4) 


enableWGCNAThreads()

powers <- c(c(1:20), seq(from = 10, to = 100, by = 4))
powers


sft <- pickSoftThreshold(dataset_norm, powerVector = powers, verbose = 26)

sprintf("Optimal soft-power = %d", sft$powerEstimate)


corRaw = abs(cor(dataset_norm))
corSoft = abs(corRaw**sft$powerEstimate)
par(mfrow=c(1,2))
hist(corRaw, main="Raw correlations")
hist(corSoft, main="Soft-thresholded correlations")

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])* sft$fitIndices[,2],
     xlab = "Soft Threshold power", ylab = "Scale Free Topology Model Fit, R²",
     type = "n", main = "Scale independence", cex.lab = 1.3)
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])* sft$fitIndices[,2], labels = powers, cex = 1, col = "black")
abline(h = 0.80, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold power", ylab = "Mean connectivity", 
     type = "n", main = "Mean connectivity", cex.lab = 1.3)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 1, col = "black")

softPower = sft$powerEstimate
similarityMat <- adjacency(dataset_norm, power = softPower)
head(similarityMat[, c(1:6)])

TOM <- TOMsimilarity(similarityMat)

colnames(TOM) <- colnames(similarityMat)
rownames(TOM) <- rownames(similarityMat)
head(TOM[, c(1:6)])


dissTOM <- 1 - TOM
head(dissTOM[, c(1:6)])


geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, 
                             distM = dissTOM, 
                             deepSplit = 2, 
                             pamRespectsDendro = FALSE, 
                             minClusterSize = minModuleSize)

table(dynamicMods)

dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

plotDendroAndColors(dendro = geneTree, colors = dynamicColors, groupLabels = "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

MEList <- moduleEigengenes(dataset_norm, colors = dynamicColors)
MEs <- MEList$eigengenes
head(MEs)

MEs <- orderMEs(MEs)
head(MEs)

MEDiss <- 1 - cor(MEs)
head(MEDiss)


METree <- hclust(as.dist(MEDiss))

MEDissCut <- 0.25
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h = MEDissCut, col = "red")

# # in case you want to merge some modules
# merge <- mergeCloseModules(dataset, dynamicColors, cutHeight = MEDissCut, verbose = 3)
# mergedColors <- merge$colors
# mergedMEs <- merge$newMEs
# plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
#                     c("Dynamic Tree Cut", "Merged dynamic"),
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# 
# moduleColors <- mergedColors
# colorOrder <- c("grey", standardColors(50))
# moduleLabels <- match(moduleColors, colorOrder) - 1 
# MEs <- mergedMEs

# save the files
# moduleColors = dynamicColors # Or mergedColors, if some modules have been merged
# names(moduleColors) = colnames(dataset)
# head(moduleColors)

# save(moduleColors, file = "results_TCGA/modules_mrna.rda")
# Store the unique bagger used in the module detection 
# save(MEs, file = "results_TCGA/MEs_mrna.rda")




# -------------------------------------------------------------------------

metadataNum <- data.frame("AgeDiagnosis"=scale(metadata$AgeDiagnosis), row.names = rownames(metadata))

AnatomicNeoplasm = binarizeCategoricalVariable(metadata$AnatomicNeoplasm,includePairwise = FALSE, includeLevelVsAll = TRUE)
rownames(AnatomicNeoplasm) = rownames(metadata)
metadataNum = cbind(metadataNum, AnatomicNeoplasm)

moduleTraitCor <- cor(MEs, metadataNum)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)



# module -trait relationship
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
colnames(textMatrix) <- colnames(moduleTraitCor)
rownames(textMatrix) <- rownames(moduleTraitCor)
for(r in row.names(moduleTraitCor)){
  for(c in colnames(moduleTraitCor)){
    if(moduleTraitPvalue[r, c] > 0.05){
      moduleTraitCor[r, c] <- 0
      textMatrix[r, c] <- ""
    }
  }
}
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(metadataNum),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors =  colorRampPalette(c("darkgreen", "white", "darkmagenta"))(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "Module-trait relationship")

# identify the gene correlation
moduleMembership <- as.data.frame(cor(dataset, MEs))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(moduleMembership), nSamples))

modNames <- substring(names(MEs), 3)
names(moduleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

head(moduleMembership[,1:5])
head(MMPvalue[,1:5])

#opt: Attribute of Interest Node Significance
geneTraitSignificance <- as.data.frame(cor(dataset, metadataNum, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS", names(metadataNum), sep = "")
names(GSPvalue) <- paste("p.GS", names(metadataNum), sep = "")


# scattor plot
module='turquoise'
moduleGenes = names(moduleColors)[which(moduleColors==module)]
geneTraitSignifColumn = "GSSubtype"
verboseScatterplot(abs(moduleMembership[moduleGenes, paste('MM', module, sep="")]),
                   abs(geneTraitSignificance[moduleGenes, geneTraitSignifColumn]),
                   xlab = paste("Module Membership (MM) in", module, "module"),
                   ylab = paste("Gene Significance (GS) for", geneTraitSignifColumn),
                   main = paste("Module Membership (MM) vs Gene Significance (GS)\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# export to downstrem analysis 
resultsModuleMembership = cbind(moduleMembership, MMPvalue)
mmNames = names(resultsModuleMembership)
mmNames = mmNames[order(gsub("p\\.", "", mmNames))]
resultsModuleMembership = resultsModuleMembership[,mmNames]

resultsSignificance = cbind(geneTraitSignificance, GSPvalue)
gsNames = names(resultsSignificance)
gsNames = gsNames[order(gsub("p\\.", "", gsNames))]
resultsSignificance = resultsSignificance[,gsNames]


results = data.frame("Gene"=names(moduleColors), "Module"=unlist(moduleColors), 
                     resultsModuleMembership, resultsSignificance)
results = results[order(results$Module),]
write.table(results, file="results_TCGA/TCGA_mRNA_results.tsv", sep="\t", row.names=F, col.names=T, quote=F)
head(results[,1:6])



load("results_TCGA/MEs_mirna.rda")
mirna = MEs
names(mirna) = paste("mirna", names(mirna), sep="_")
load("results_TCGA/MEs_proteins.rda")
protein = MEs
names(protein) = paste("protein", names(protein), sep="_")
load("results_TCGA/MEs_mrna.rda")
mrna = MEs
names(mrna) = paste("mrna", names(mrna), sep="_")

samples = Reduce("intersect", list(rownames(mirna), rownames(protein), rownames(mrna)))
MEs = do.call("cbind", list(mirna[samples,], mrna[samples,], protein[samples,]))


softPower = 1 # this time, we don't need a scale-free topology, so we just set softPower to 1
similarityMat <- adjacency(MEs, power = softPower)

annot = data.frame("Modality"=gsub("_.*", "", colnames(similarityMat)))
rownames(annot) = colnames(similarityMat)
head(annot)

heatmap <- pheatmap(similarityMat, show_rownames = TRUE, show_colnames = TRUE, annotation_row = annot, annotation_col = annot,
                    legend = T, cluster_cols = T, cutree_cols = 5, cluster_rows = T, cutree_rows  = 5)

clust <- cbind(similarityMat, cluster = cutree(heatmap$tree_col, k = 5))
clust = clust[,"cluster"]
clust[order(clust)]





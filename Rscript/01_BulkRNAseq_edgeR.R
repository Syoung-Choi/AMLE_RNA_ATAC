getwd()
setwd(dir = "C:/Users/SY/Desktop/Work_NOW/Project/AMLE/")



# Packages load -----------------------------------------------------------


# BiocManager::install("Mus.musculus")

library(dplyr)
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(Mus.musculus)
library(tidyverse)
#local_files===============================================================================

AMLE <- read.csv("./00 Raw/DMSO_APP_Readcount.csv")

AMLE<- distinct(AMLE,geneID, .keep_all = T) 
AMLE<- column_to_rownames(AMLE,var = "geneID")


AMLE[is.na(AMLE)] <- 0

AMLE <- as.matrix(AMLE)

na_indices <- which(is.na(AMLE), arr.ind = TRUE)

any_na <- any(is.na(AMLE))
if (any_na) {
  print("have NA.")
} else {
  print("not NA.")
}

# write.csv(AMLE,"./DMSO_APP_Readcount.csv")

# data load # ---------------------------------------------------------------
AMLE <- read.csv("./00 Raw/DMSO_APP_Readcount.csv",header = T, row.names = 1)

dim(AMLE)
names(AMLE)
str(AMLE)

class(AMLE)
# Convert to DEGlist, Gene infomation annotation #-----------------------------------

dgeObj <-DGEList(AMLE)

dim(dgeObj)
names(dgeObj)
str(dgeObj)

class(dgeObj)

dgeObj$samples

samplenames <- substring(colnames(dgeObj),1,nchar(colnames(dgeObj)))
samplenames

colnames(dgeObj) <- samplenames

group <- as.factor(c("DMSO","DMSO","DMSO","AMLE","AMLE","AMLE"))

dgeObj$samples$group <- group

dgeObj$samples

head(dgeObj$counts)
dim(dgeObj$counts)

keytypes(org.Mm.eg.db)



geneid <- rownames(dgeObj)

genes <- AnnotationDbi::select(org.Mm.eg.db,
                keys = geneid,
                columns = c("GENENAME","GENETYPE","ENTREZID","MGI"),
                keytype = "SYMBOL")



head(genes)

genes <- genes[!duplicated(genes$SYMBOL), ]

dgeObj$genes<- genes

dgeObj$genes

#########################FT Protine-coding use for WGCNA################################
# protein_coding_count <- sum(genes$GENETYPE == "ncRNA", na.rm = TRUE)
# cat("Number of protein coding genes:", protein_coding_count, "\n")
# 
# # -------------------------------------------------------------------------
# 
# 
# protein_coding_genes <- genes[genes$GENETYPE == "protein-coding", "SYMBOL"]
# 
# # filter protein coding genes
# filtered_counts <- dgeObj[rownames(dgeObj$counts) %in% protein_coding_genes, ]
# str(filtered_counts)
# #[1] 24419 -> 17934
#
# output_file <- "filtered_cpm_protein_coding.csv"
# write.csv(filtered_counts, file = output_file, row.names = TRUE)
#########################FT Protine-coding use for WGCNA################################

# Filter Genes ------------------------------------------------------------

CPM <- cpm(AMLE)
lcpm <- cpm(AMLE, log = TRUE)


plotDensities(lcpm, legend = F, main = "Before filterging")
 
abline(v = 0, lty = 3)

# save plot
# dir.create("../01.QC/plots")
# png(filename = "../01.QC/plots/distribution of gene expression plot(unnomrailzed).png", width =8, height = 5, units = "in", res = 300)
# plotDensities(lcpm, legend = F, main = "Before filterging")
# abline(v = 0, lty = 3)
# title("BMDC : MT treatmented MDS")
# dev.off()


# opt : only keep genes which have CPM greater then 1 in least 2 samples
keep.exprs <- rowSums(CPM > 1) >= 2
dgeObj <- dgeObj[keep.exprs, , keep.lib.sizes = FALSE]

dim(dgeObj) 

lcpm <- cpm(dgeObj, log=TRUE)

plotDensities(lcpm, legend = FALSE, main = "After filtering")
abline(v = 0, lty = 3)


# png(filename = "../01.QC/plots/distribution of gene expression plot(Nomrailzed).png", width =8, height = 5, units = "in", res = 300)
plotDensities(lcpm, legend = F, main = "Before filterging")
abline(v = 0, lty = 3)
# dev.off()


#  Normalization ----------------------------------------------------------

dgeObj$samples$norm.factors 

dgeObj <- calcNormFactors(dgeObj, method = "TMM")



# Data Exploration --------------------------------------------------------

levels(col.group) <-c("green","black")

col.group <- as.character(col.group)

plotMDS(lcpm, labels = group, col = col.group,
        main = "group")


# Construct the linear model ----------------------------------------------

# First, we will make a design matrix. This holds information about our samples

design <- model.matrix(~0 + group)
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(AMLEvsDMSO = AMLE - DMSO,
                              levels = colnames(design))

# contr.matrix <- makeContrasts(BasalvsLP = Basal - LP,
#                               BasalvsML = Basal - ML,
#                               LPvsML = LP - ML,
#                               levels = colnames(design))

contr.matrix


v <- voom(dgeObj, design, plot = TRUE)

v


vfit <- lmFit(v, design)

View(vfit$coefficients)
View(vfit$stdev.unscaled)

# Caculate the statistics for our specific contrasts of interest

vfit <- contrasts.fit(vfit, contrasts = contr.matrix)

# View(vfit$coefficients)
# View(vfit$stdev.unscaled)

# Calculate the t statistics for each gene in each comparison (hypothesis)
# 
efit <- eBayes(vfit)
# View(efit$t)

# Explore the relationship between residual variation versus the mean gene expression level. 

plotSA(efit, main = "Final model: Mean−variance trend")




# Tabulate the results

summary(decideTests(efit))

tfit <- treat(vfit, lfc = 1)
dt <- decideTests(tfit)
summary(dt)

head(dt)

de.common <- which(dt[, 1] != 0)
length(de.common)

head(tfit$Amean[de.common], n = 10)

Test <- as.data.frame(dt)
#write.fit(tfit, dt, file = "results.txt")

# Identify the top DE genes
AMLEvsDMSO <- topTreat(tfit, coef = 1, n = Inf)
head(AMLEvsDMSO)

#  MA plot
plotMD(tfit, column = 1, status = dt[, 1], main = colnames(tfit)[1],
       xlim = c(-8, 13))

glMDPlot(tfit, coef = 1, status = dt, main = colnames(tfit)[1],
         side.main = "SYMBOL", counts = AMLE$counts, groups = group,
         launch = TRUE)

# volcanoplot
ggplot(AMLEvsDMSO, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color =P.Value > 0.05 | abs(logFC)< 2)) +  # p-value < 0.05 and |logFC| > 2
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10 P-value") +
  scale_color_manual(values = c("red", "black"), name = "Significant")


# interactive volcano plot
# if (!requireNamespace("plotly", quietly = TRUE))
#   install.packages("plotly")
# if (!requireNamespace("ggplot2", quietly = TRUE))
#   install.packages("ggplot2")

library(ggplot2)
library(plotly)

ggplotly(
  ggplot(AMLEvsDMSO, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(color = P.Value > 0.05 | abs(logFC)< 2, text = paste("Gene:", rownames(AMLEvsDMSO)))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "blue") +
    theme_minimal() +
    labs(title = "Interactive Volcano Plot",
         x = "log2 Fold Change",
         y = "-log10 P-value") +
    scale_color_manual(values = c("red", "black"), name = "Significant"),
  tooltip = "text"
)



# Heatmap -------------------------------------------------------------

# View heatmap of top 100 DE genes between Basal and LP cells

library("gplots")

AMLEvsDMSO_toptages <- AMLEvsDMSO$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% AMLEvsDMSO_toptages)

mycol <- colorpanel(1000, "blue", "white", "red")
heatmap.2(v$E[i, ], scale = "row",
          labRow = v$genes$SYMBOL[i], labCol = group,
          col = mycol, trace = "none", density.info = "none",)

# -------------------------------------------------------------------------

#                          GSEA analysis

# -------------------------------------------------------------------------

# BiocManager::install("goseq")


####################################
##ClusterProfiler
library(clusterProfiler)
library(enrichplot)
library(GOSemSim)
library(DOSE)
#  App Broc_PPS Broc_PS Carr ChaM CowP DMSO DMSO_1 DMSO004 DMSO004_1 Fuco GavL GreT KRGE1000 LPS Lute NonSaponin Oste RGC_L RGC_P Saponin Spir VD3
DGE.List <- estimateDisp(dgeObj,design)
fit <- glmFit(DGE.List)
lrt <-glmLRT(fit, coef = 2,contrast=c(1,-1))

DGE_list <- topTags(lrt, n = Inf, sort.by = "PValue", adjust.method = "BH")


DGE_list %>% head

row.names(DGE_list) -> DGE_list_CP$table$SYMBOL 

result_deseq <- DGE_list$table

# logFC 
result_deseq %>% 
  mutate(
    rank = rank(logFC, ties.method = "random")
  ) %>% 
  arrange(desc(rank)) %>% data.frame() -> gsea_data


gsea_data$logFC -> gsea_vector

gsea_data$SYMBOL -> names(gsea_vector)


####GO: BP,MF,CC
gseGO(gsea_vector,               
      OrgDb = org.Mm.eg.db,      
      keyType = "SYMBOL",    
      ont = "ALL",     
      pvalueCutoff = 0.05,       
      pAdjustMethod = "fdr",   
      minGSSize = 10,   
      maxGSSize = 500,
      verbose = T,
      nPermSimple = 100000,
      eps = 0) -> IgAN_gsea_GO


# Group similar GO terms
simplify(IgAN_gsea_GO, cutoff=0.7, by="p.adjust", select_fun=min) -> IgAN_gsea_GO_simp

#run this to draw a tree plot.
pairwise_termsim(IgAN_gsea_GO )  -> IgAN_gsea_GO

# Gene annotation type switch
# bitr(gsea_data$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db") -> ensembl_entrezid
# left_join(gsea_data, ensembl_entrezid, by = "SYMBOL") -> gsea_data_entrezid
# gsea_data_entrezid -> gsea_data


  gsea_data %>%
  filter(!is.na(SYMBOL)) -> gsea_data

gsea_data$logFC -> gsea_vector_data
gsea_data$ENTREZID -> names(gsea_vector_data)

# KEGG
gseKEGG(gsea_vector_data, 
        organism = "mmu",   
        keyType = "kegg",    
        pvalueCutoff = 0.05,
        pAdjustMethod = "fdr", 
        verbose = T,
        maxGSSize = 500,
        nPermSimple = 10000,
        eps = 0) -> IgAN_gsea_KEGG

IgAN_gsea_KEGG@result$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", IgAN_gsea_KEGG@result$Description)
head(IgAN_gsea_KEGG@result$Description)


# WikiPathways

gseWP(
  gsea_vector_data, 
  organism = "Mus musculus",
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr", 
  verbose = T,
  nPermSimple = 100000,
  eps = 0
) -> IgAN_gsea_WP

# #######################Visualization
# 
# select_top_GSEA <- function(GSEA_result, GO = "no") {
# 
#   if(GO %in%  c("MF", "BP", "CC") ) {      
# 
#     c((GSEA_result %>% data.frame() %>%
#          filter(ONTOLOGY == GO) %>%
#          filter(NES >= 0) %>%
#          dplyr::top_n(n = 10, wt = NES))$Description,    
#       (GSEA_result %>% data.frame() %>%
#          filter(ONTOLOGY == GO) %>%
#          filter(NES <= 0) %>%
#          dplyr::top_n(n = -10, wt = NES))$Description) -> top_pathway_GO  
S#     GSEA_result %>%
#       filter(Description %in% top_pathway_GO) %>%
#       data.frame() %>%
#       arrange(NES) -> top_pathway_data

#   } else {   
# 
#     c((GSEA_result %>% data.frame() %>%
#          filter(NES >= 0) %>%
#          dplyr::top_n(n = 10, wt = NES))$Description,
#       (GSEA_result %>% data.frame() %>%
#          filter(NES <= 0) %>%
#          dplyr::top_n(n = -10, wt = NES))$Description) -> top_pathway
# 
#     GSEA_result %>%
#       filter(Description %in% top_pathway) %>%
#       data.frame() %>%
#       arrange(NES) -> top_pathway_data
# 
#   }
# 
#   factor(top_pathway_data$Description, levels = c(top_pathway_data$Description)) -> top_pathway_data$Description
# 
#   return(top_pathway_data)
# 
# }
# 
# split_pathway_name <- function(gsea_top_data, char_number = 40) {
#   
#   as.character(gsea_top_data[["Description"]]) -> label_gsea
#   
#   data.frame(Position = which(nchar(label_gsea) >= char_number),
#              Description = label_gsea[which(nchar(label_gsea) >= char_number)] ) -> df_label_gsea
#   
#   if(!is.na(df_label_gsea[1,1])) {
#     
#     for(i in 1:NROW(df_label_gsea) ) {
#       
#       which.min( abs(
#         (nchar(df_label_gsea$Description)[[i]] /2 ) - str_locate_all(df_label_gsea$Description, " ")[[i]][, 1] ) ) -> position_min
#       
#       str_locate_all(df_label_gsea$Description, " ")[[i]][, 1][position_min] -> position_min
#       
#       str_sub(label_gsea[ df_label_gsea$Position[[i]] ], position_min, position_min) <- "\n"
#       
#     }
#     
#     return(label_gsea)
#     
#   } else {
#     
#     return(label_gsea)
#     
#   }
#   
# }
# ##############################################################################
GOBP_top <- select_top_GSEA(GSEA_result = IgAN_gsea_GO_simp,GO = "BP" )
GOMF_top <- select_top_GSEA(GSEA_result = IgAN_gsea_GO_simp,GO = "MF" )
GOCC_top <- select_top_GSEA(GSEA_result = IgAN_gsea_GO_simp,GO = "CC" )
KEGG_top <- select_top_GSEA(GSEA_result = IgAN_gsea_KEGG )
Wiki_top <- select_top_GSEA(GSEA_result = IgAN_gsea_WP)

TopPathways<- list()
TopPathways[["GOBP"]] <- GOBP_top
TopPathways[["GOMF"]] <- GOMF_top
TopPathways[["GOCC"]] <- GOCC_top
TopPathways[["KEGG"]] <- KEGG_top
TopPathways[["Wiki"]] <- Wiki_top

path<- lapply(TopPathways, function(x) x$Description)
path

Input_Data <-  TopPathways$KEGG ## edit this line to change the pathway
titlename <- "Applemango GSEA: KEGG"

top_data <- split_pathway_name(gsea_top_data = Input_Data,char_number = 40)

# p1<-ggplot(Input_Data %>% filter(ONTOLOGY == "BP")) +
p1<-ggplot(Input_Data) +
  geom_bar(aes(x = Description, y = NES, fill = p.adjust),
           stat='identity') +
  scale_fill_continuous(low = "red", high = "blue")+
  coord_flip(ylim = c(-3.2, 3.2))+
  scale_x_discrete(labels = top_data )+
  scale_y_continuous(breaks = seq(-3,3,1))+
  theme_bw(14) +
  labs(title=titlename) +  
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank())
p1
# CP_lrt

ggsave(filename = "Nonsaponin vs Saponin Wiki.png",plot = p1,path = "../240829_Main_analysis/plots/",width = 20,height = 18,units = "cm", dpi = 600)





############################################################


# -------------------------------------------------------------------------

#                                sessionInfo

# -------------------------------------------------------------------------
# 
# sessionInfo()
# R version 4.4.1 (2024-06-14 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 10 x64 (build 19043)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=Korean_Korea.utf8  LC_CTYPE=Korean_Korea.utf8    LC_MONETARY=Korean_Korea.utf8 LC_NUMERIC=C                 
# [5] LC_TIME=Korean_Korea.utf8    
# 
# time zone: Asia/Seoul
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] DOSE_3.30.5                               GOSemSim_2.30.2                           enrichplot_1.24.4                        
# [4] clusterProfiler_4.12.6                    goseq_1.56.0                              geneLenDataBase_1.40.1                   
# [7] BiasedUrn_2.0.12                          plotly_4.10.4                             lubridate_1.9.3                          
# [10] forcats_1.0.0                             stringr_1.5.1                             purrr_1.0.2                              
# [13] readr_2.1.5                               tidyr_1.3.1                               tibble_3.2.1                             
# [16] ggplot2_3.5.1                             tidyverse_2.0.0                           Mus.musculus_1.3.1                       
# [19] TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0 org.Mm.eg.db_3.19.1                       GO.db_3.19.1                             
# [22] OrganismDbi_1.46.0                        GenomicFeatures_1.56.0                    GenomicRanges_1.56.1                     
# [25] GenomeInfoDb_1.40.1                       AnnotationDbi_1.66.0                      IRanges_2.38.1                           
# [28] S4Vectors_0.42.1                          Biobase_2.64.0                            BiocGenerics_0.50.0                      
# [31] RColorBrewer_1.1-3                        gplots_3.1.3.1                            Glimma_2.14.0                            
# [34] edgeR_4.2.1                               limma_3.60.4                              dplyr_1.1.4                              
# 
# loaded via a namespace (and not attached):
#   [1] splines_4.4.1               BiocIO_1.14.0               ggplotify_0.1.2             bitops_1.0-8               
# [5] filelock_1.0.3              R.oo_1.26.0                 polyclip_1.10-6             graph_1.82.0               
# [9] XML_3.99-0.17               lifecycle_1.0.4             httr2_1.0.5                 MASS_7.3-61                
# [13] lattice_0.22-6              crosstalk_1.2.1             magrittr_2.0.3              yaml_2.3.10                
# [17] cowplot_1.1.3               DBI_1.2.3                   abind_1.4-8                 zlibbioc_1.50.0            
# [21] R.utils_2.12.3              ggraph_2.2.1                RCurl_1.98-1.16             yulab.utils_0.1.7          
# [25] tweenr_2.0.3                rappdirs_0.3.3              sva_3.35.2                  GenomeInfoDbData_1.2.12    
# [29] ggrepel_0.9.5               tidytree_0.4.6              genefilter_1.84.0           annotate_1.82.0            
# [33] codetools_0.2-20            DelayedArray_0.30.1         ggforce_0.4.2               xml2_1.3.6                 
# [37] tidyselect_1.2.1            aplot_0.2.3                 UCSC.utils_1.0.0            farver_2.1.2               
# [41] viridis_0.6.5               matrixStats_1.3.0           BiocFileCache_2.12.0        GenomicAlignments_1.40.0   
# [45] jsonlite_1.8.8              tidygraph_1.3.1             survival_3.7-0              tools_4.4.1                
# [49] progress_1.2.3              treeio_1.28.0               snow_0.4-4                  Rcpp_1.0.13                
# [53] glue_1.7.0                  gridExtra_2.3               SparseArray_1.4.8           mgcv_1.9-1                 
# [57] DESeq2_1.44.0               qvalue_2.36.0               MatrixGenerics_1.16.0       withr_3.0.1                
# [61] BiocManager_1.30.25         fastmap_1.2.0               fansi_1.0.6                 caTools_1.18.2             
# [65] digest_0.6.36               gridGraphics_0.5-1          timechange_0.3.0            R6_2.5.1                   
# [69] colorspace_2.1-0            gtools_3.9.5                biomaRt_2.60.1              RSQLite_2.3.7              
# [73] R.methodsS3_1.8.2           utf8_1.2.4                  generics_0.1.3              data.table_1.15.4          
# [77] rtracklayer_1.64.0          graphlayouts_1.2.0          prettyunits_1.2.0           httr_1.4.7                 
# [81] htmlwidgets_1.6.4           S4Arrays_1.4.1              scatterpie_0.2.4            pkgconfig_2.0.3            
# [85] gtable_0.3.5                blob_1.2.4                  XVector_0.44.0              shadowtext_0.1.4           
# [89] htmltools_0.5.8.1           fgsea_1.30.0                RBGL_1.80.0                 scales_1.3.0               
# [93] png_0.1-8                   ggfun_0.1.6                 rstudioapi_0.17.0           tzdb_0.4.0                 
# [97] reshape2_1.4.4              rjson_0.2.21                nlme_3.1-165                curl_5.2.3                 
# [101] cachem_1.1.0                KernSmooth_2.23-24          parallel_4.4.1              restfulr_0.0.15            
# [105] pillar_1.9.0                grid_4.4.1                  vctrs_0.6.5                 dbplyr_2.5.0               
# [109] xtable_1.8-4                cli_3.6.3                   locfit_1.5-9.10             compiler_4.4.1             
# [113] Rsamtools_2.20.0            rlang_1.1.4                 crayon_1.5.3                labeling_0.4.3             
# [117] plyr_1.8.9                  fs_1.6.4                    stringi_1.8.4               viridisLite_0.4.2          
# [121] BiocParallel_1.38.0         txdbmaker_1.0.1             munsell_0.5.1               Biostrings_2.72.1          
# [125] lazyeval_0.2.2              Matrix_1.7-0                patchwork_1.3.0             hms_1.1.3                  
# [129] bit64_4.0.5                 KEGGREST_1.44.1             statmod_1.5.0               SummarizedExperiment_1.34.0
# [133] igraph_2.0.3                memoise_2.0.1               ggtree_3.12.0               fastmatch_1.1-4            
# [137] bit_4.0.5                   gson_0.1.0                  ape_5.8     
# 

# -------------------------------------------------------------------------

# Quick start -------------------------------------------------------------
x <- read.csv("./DMSO_APP_Readcount.csv",header = T, row.names = 1)

group <- factor(c(1,1,1,2,2,2)) # in orignal tutorial c(1,1,2,2) our sample has 3 replicate 
y <- DGEList(counts=x,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- normLibSizes(y)
design <- model.matrix(~group)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)

topTags(qlf)
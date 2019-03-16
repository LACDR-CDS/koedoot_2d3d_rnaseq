library("ggplot2")
library("reshape2")
library("tidyr")
library("pheatmap")
library("readr")
library("dplyr")

# To avoid data.frame vs tibble issues (eg with deprecated rownames)
read_tsv <- function(...) {
  as.data.frame(suppressMessages(readr::read_tsv(...)))
}

parent.dir <- "Add working directory" 
setwd(parent.dir)
shared.files.dir <- file.path("..", "FilesShared") # Set location of FilesShared folder here

GenesDF <- read.table("MultipleGeneTable.txt", sep = "\t", header = TRUE)
SetName <- "Name Gene Set"
GeneIDs <- GenesDF$EntrezID
EnsemblIDs <- as.character(GenesDF$Ensembl)

# Expression different cell lines 2D

dir.create(paste(parent.dir, "CellLine2D", sep="/")) # create a new folder to put the files in
outputdir <- paste(parent.dir, "CellLine2D/", sep="/")

CellLineExpression_2D <- read_tsv(file.path(shared.files.dir, "LogNormalizedCounts_2DCellLines.txt.xz")) # Load expression data
rownames(CellLineExpression_2D) <- CellLineExpression_2D$EnsemblID # Put the ensemblID in the rownames
CellLineExpression_2D$EnsemblID <- NULL # Remove the ensemblID column now they are in the rownames
CellLineExpression_2D <- as.data.frame(t(CellLineExpression_2D)) # transform the data frame so the cell lines are in the row names
CellLineExpression_2D$CellLine <- rownames(CellLineExpression_2D) # Create a cell line column that can be used to merge data frames

CellLine_Subtype <- read_tsv(file.path(shared.files.dir,"Subtype_CellLines.txt")) # Load the file containing cell line subtype information

Subtype_Expression_DF <- merge(CellLine_Subtype, CellLineExpression_2D)

GeneDF <- subset(Subtype_Expression_DF, select = c("CellLine", "Subtype", EnsemblIDs)) # select the genes of interest
AnnotationDF <- subset(GeneDF, select = c("Subtype")) # Create an annotation data frame for the heatmap
rownames(AnnotationDF) <- GeneDF$CellLine

rownames(GeneDF) <- GeneDF$CellLine # put the cell lines in the row names
GeneDF <- as.data.frame(t(subset(GeneDF, select = -c(CellLine, Subtype)))) # remove cell line and subtype columns from the data frame and transform
GeneDF$Ensembl <- rownames(GeneDF) # create a column with ensemblIDs, this will be used to convert ensembl IDs to gene symbols
GeneDF <- merge(GeneDF, GenesDF) # add ensemblID and gene symbol information

write.table(GeneDF, file = paste(outputdir, SetName, "_Log2Expression.txt", sep = ""), sep = "\t", row.names = FALSE)

rownames(GeneDF) <- GeneDF$GeneSymbol
GeneDF <- subset(GeneDF, select = -c(GeneSymbol, EntrezID, Ensembl))

AnnotationColors <- list(Subtype = c(Luminal = "darkorchid3", BasalA = "orange1", BasalB = "dodgerblue3"))

pdf(paste(outputdir, SetName, "_Log2Expression.pdf")  , width = 7, height = 8, onefile=FALSE)
pheatmap(as.matrix(GeneDF),
         color = colorRampPalette(c("green", "black", "red")) ( 100 ),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         annotation_col = AnnotationDF,
         annotation_colors = AnnotationColors,
         fontsize_col = 4,
         fontsize_row = 1,
         cellheight = 2,
         cellwidth = 5)
dev.off()

GenesToRemove <- read.table("GenesToRemove_2D.txt", sep = "\t", header = FALSE)[, c(1)]
GeneDF2 <- GeneDF[!(row.names(GeneDF) %in% GenesToRemove), ]

pdf(paste(outputdir, SetName, "_Log2Expression_AfterSelection.pdf")  , width = 7, height = 8, onefile=FALSE)
pheatmap(as.matrix(GeneDF2),
         color = colorRampPalette(c("green", "black", "red")) ( 100 ),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         annotation_col = AnnotationDF,
         annotation_colors = AnnotationColors,
         fontsize_col = 4,
         fontsize_row = 1,
         cellheight = 2,
         cellwidth = 5)
dev.off()

GeneDF2_MedianNorm <- GeneDF2 - apply(GeneDF2, 1, median)

bk = unique(c(seq(min(GeneDF2_MedianNorm),-0.001, length=50), 0, seq(0.0001, max(GeneDF2_MedianNorm), length=50)))
pdf(paste(outputdir, SetName, "_Log2Expression_MedianNormalizedPerGene.pdf")  , width = 7, height = 8, onefile=FALSE)
pheatmap(as.matrix(GeneDF2_MedianNorm),
         color = colorRampPalette(c("green", "black", "red")) ( 100 ),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         annotation_col = AnnotationDF,
         annotation_colors = AnnotationColors,
         breaks = bk, 
         fontsize_col = 4,
         fontsize_row = 1,
         cellheight = 2,
         cellwidth = 5)
dev.off()


# Expression 2D vs 3D

dir.create(paste(parent.dir, "CellLine2Dvs3D", sep="/")) # create a new folder to put the files in
outputdir <- paste(parent.dir, "CellLine2Dvs3D/", sep="/")

NormalizedCounts <- read_tsv(file.path(shared.files.dir, "NormalizedCounts_2Dvs3D.txt.xz"))
MetaDataDF <- read_tsv(file.path(shared.files.dir, "MetaData_2Dvs3D.txt"))

NormalizedCounts_Genes <- NormalizedCounts[NormalizedCounts$ensembl %in% EnsemblIDs, ] # select genes of interest
rownames(NormalizedCounts_Genes) <- NormalizedCounts_Genes$symbol
NormalizedCounts_Genes <- subset(NormalizedCounts_Genes, select = -c(symbol, entrez, ensembl))
colnames(NormalizedCounts_Genes) <- gsub("_A","",colnames(NormalizedCounts_Genes))

AnnotationDF <- unique(subset(MetaDataDF, select = c("Subtype", "CellLine"))) # Create an annotation data frame
AnnotationDF$CellLine <- gsub("-", "", AnnotationDF$CellLine)
rownames(AnnotationDF) <- AnnotationDF$CellLine
AnnotationDF$CellLine <- NULL
AnnotationColors <- list(Subtype = c(Luminal = "darkorchid3", BasalA = "orange1", BasalB = "dodgerblue3"))

pdf(paste(outputdir, SetName, "_Log2Expression_2Dvs3DSamples.pdf")  , width = 6, height = 8, onefile=FALSE)
pheatmap(as.matrix(NormalizedCounts_Genes),
         color = colorRampPalette(c("green", "black", "red")) ( 100 ),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         fontsize_col = 4,
         fontsize_row = 1,
         cellheight = 2,
         cellwidth = 5)
dev.off()

NormalizedCounts_Genes$Symbol <- rownames(NormalizedCounts_Genes)
write.table(NormalizedCounts_Genes, file = paste(outputdir, SetName, "_2Dvs3D_Log2NormalizedCounts.txt", sep = ""), sep = "\t", row.names = FALSE)
NormalizedCounts_Genes$Symbol <- NULL

GenesToRemove <- read.table("GenesToRemove_2Dvs3D.txt", sep = "\t", header = FALSE)[, c(1)]
NormalizedCounts_Genes2 <- NormalizedCounts_Genes[!(row.names(NormalizedCounts_Genes) %in% GenesToRemove), ]

pdf(paste(outputdir, SetName, "_Log2Expression_2Dvs3DSamples_AfterSelection.pdf")  , width = 6, height = 8, onefile=FALSE)
pheatmap(as.matrix(NormalizedCounts_Genes2),
         color = colorRampPalette(c("green", "black", "red")) ( 100 ),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         fontsize_col = 4,
         fontsize_row = 1,
         cellheight = 2,
         cellwidth = 5)
dev.off()

NormalizedCounts_Genes3 <- NormalizedCounts_Genes2 # Reformat the data frame and calculate log fold changes
NormalizedCounts_Genes3$Symbol <- rownames(NormalizedCounts_Genes3)
NormalizedCounts_Long <- melt(NormalizedCounts_Genes3, id.vars = c("Symbol"), 
                              variable.name = "Sample", value.name = "Expression") # create a long data frame
NormalizedCounts_Long2 <- separate(data = NormalizedCounts_Long, col = Sample, into = c("CellLine", "X2D.3D"), sep = "\\_")
NormalizedCounts2 <- dcast(NormalizedCounts_Long2, Symbol + CellLine ~ X2D.3D, value.var = "Expression")
colnames(NormalizedCounts2) [3:4] <- c("Expr_2D", "Expr_3D")
NormalizedCounts2$FoldChange2Dvs3D <- NormalizedCounts2$Expr_2D - NormalizedCounts2$Expr_3D
NormalizedCounts3 <- subset(NormalizedCounts2, select = -c(Expr_2D, Expr_3D))
NormalizedCounts4 <- dcast(NormalizedCounts3, Symbol ~ CellLine, value.var = "FoldChange2Dvs3D")
write.table(NormalizedCounts4, file = paste(outputdir, SetName, "_2Dvs3DFoldChange.txt", sep = ""), sep = "\t", row.names = FALSE)

rownames(NormalizedCounts4) <- NormalizedCounts4$Symbol
NormalizedCounts4$Symbol <- NULL

NormalizedCounts4 <- NormalizedCounts4[, c("BT474", "MDAMB175", "SKBR3", "T47D", "UACC893", "ZR751",
                                           "HCC1143", "HCC1806", "HCC1954", "MDAMB468",
                                           "BT549", "Hs578T", "MDAMB231", "SUM149PT")]
bk = unique(c(seq(min(NormalizedCounts4),-0.001, length=50), 0, seq(0.0001, max(NormalizedCounts4), length=50)))

pdf(paste(outputdir, SetName, "_LogFoldChange_2Dvs3D.pdf")  , width = 6, height = 8, onefile=FALSE)
pheatmap(as.matrix(NormalizedCounts4),
         color = colorRampPalette(c("green", "black", "red")) ( 100 ),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         breaks = bk, 
         annotation_col = AnnotationDF,
         annotation_colors = AnnotationColors,
         fontsize_col = 4,
         fontsize_row = 1,
         cellheight = 2,
         cellwidth = 5)
dev.off()

# Expression in different subtypes TCGA

dir.create(paste(parent.dir, "TCGAsubtypes", sep="/")) # create a new folder to put the files in
outputdir <- paste(parent.dir, "TCGAsubtypes/", sep="/")

TCGA_Subtype_DF <- read_tsv(file.path(shared.files.dir, "TCGA_PrimaryTumor_Subtype_Expression.txt.xz")) # Load TCGA data
TCGA_Subtype_DF[1:5, 1:5]
GeneDF <- select(TCGA_Subtype_DF, one_of(c("TumorType", paste("ID", GeneIDs, sep = "_"))))
GeneDF$PatientID <- seq(from = 1, to = nrow(GeneDF), by = 1)

AnnotationDF <- subset(GeneDF, select = c("TumorType"))
rownames(AnnotationDF) <- GeneDF$PatientID
AnnotationColors <- list(TumorType = c(TNBC = "darkorchid3", ERpos = "orange1"))
AnnotationDF2 <- AnnotationDF
AnnotationDF2$PatientID <- rownames(AnnotationDF)
write.table(AnnotationDF2, file = paste(outputdir, SetName, "_PatientSubtypeAnnotation.txt", sep = ""), sep = "\t", row.names = FALSE)

rownames(GeneDF) <- GeneDF$PatientID
GeneDF2 <- subset(GeneDF, select = -c(TumorType, PatientID))
GeneDF3 <- as.data.frame(t(GeneDF2))
GeneDF3$EntrezID <- rownames(GeneDF3)
GeneDF3$EntrezID <- gsub("ID_", "", GeneDF3$EntrezID)
GeneDF4 <- merge(GeneDF3, GenesDF)
rownames(GeneDF4) <- GeneDF4$GeneSymbol

write.table(GeneDF4, file = paste(outputdir, SetName, "_TNBCvsERpos_TCGApatients.txt", sep = ""), sep = "\t", row.names = FALSE)

GeneDF5 <- subset(GeneDF4, select = -c(Ensembl, GeneSymbol, EntrezID))

pdf(paste(outputdir, SetName, "_Log2Expression_TCGApatients_Subtypes.pdf"), width = 12, height = 8, onefile = FALSE)
pheatmap(as.matrix(GeneDF5),
         color = colorRampPalette(c("green", "black", "red")) ( 100 ),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         annotation_col = AnnotationDF,
         annotation_colors = AnnotationColors,
         fontsize_col = 4,
         fontsize_row = 1,
         cellheight = 2,
         cellwidth = .5)
dev.off()

GenesToRemove <- read.table("GenesToRemove_TCGA_Subtype.txt", sep = "\t", header = FALSE)[, c(1)]
GeneDF6 <- GeneDF5[!(row.names(GeneDF5) %in% GenesToRemove), ]

pdf(paste(outputdir, SetName, "_Log2Expression_TCGApatients_Subtypes_AfterGeneRemoval.pdf"), width = 12, height = 8, onefile = FALSE)
pheatmap(as.matrix(GeneDF6),
         color = colorRampPalette(c("green", "black", "red")) ( 100 ),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         annotation_col = AnnotationDF,
         annotation_colors = AnnotationColors,
         fontsize_col = 4,
         fontsize_row = 1,
         cellheight = 2,
         cellwidth = .5)
dev.off()

GeneDF7 <- GeneDF6 - apply(GeneDF6, 1, median)

pdf(paste(outputdir, SetName, "_Log2Expression_TCGApatients_Subtypes_MedianNorm.pdf"), width = 12, height = 8, onefile = FALSE)
pheatmap(as.matrix(GeneDF7),
         color = colorRampPalette(c("green", "black", "red")) ( 100 ),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         annotation_col = AnnotationDF,
         annotation_colors = AnnotationColors,
         fontsize_col = 4,
         fontsize_row = 1,
         cellheight = 2,
         cellwidth = .5)
dev.off()

# Expression tumor vs normal TCGA

dir.create(paste(parent.dir, "TCGA_TumorVsNormal", sep="/")) # create a new folder to put the files in
outputdir <- paste(parent.dir, "TCGA_TumorVsNormal/", sep="/")

TumorData <- read_tsv(file.path(shared.files.dir, "BRCA__illuminahiseq_rnaseqv2__GeneExp_PrimarySolidTumor.txt.xz"))
colnames(TumorData) <- c("GeneSymbol", "EntrezID", gsub(".{16}$", "", colnames(TumorData)[-c(1,2)])) # change the column names to patient IDs by removing the last 16 characters
TumorData_Genes <- TumorData[TumorData$EntrezID %in% GeneIDs, ]  # Select genes of interest
TumorData_Genes2 <- TumorData_Genes # rename the data frame
colnames(TumorData_Genes2) [-c(1,2)] <- paste("T", colnames(TumorData_Genes2) [-c(1,2)], sep = "_") # create patient IDs refering to tissue type

NormalTissueData <- read_tsv(file.path(shared.files.dir, "BRCA__illuminahiseq_rnaseqv2__GeneExp_SolidTissueNormal.txt.xz"))
colnames(NormalTissueData) <- c("GeneSymbol", "EntrezID", gsub(".{16}$", "", colnames(NormalTissueData)[-c(1,2)])) # change the column names to patient IDs by removing the last 16 characters
NormalTissueData_Gene <- NormalTissueData[NormalTissueData$EntrezID %in% GeneIDs, ] # Select gene of interest
NormalTissueData_Gene2 <- NormalTissueData_Gene # rename the data frame
colnames(NormalTissueData_Gene2) [-c(1,2)] <- paste("N", colnames(NormalTissueData_Gene2) [-c(1,2)], sep = "_") # create patient IDs refering to tissue type

NormalTumor_Comb <- merge(TumorData_Genes2, NormalTissueData_Gene2, by.all = TRUE) # merge data frames
NormalTumor_Comb[, -c(1:2)] <- log2(NormalTumor_Comb[, -c(1:2)] + 1) # perform log normalization

write.table(NormalTumor_Comb, file = paste(outputdir, SetName, "_NormalVsTumor_TCGApatients_Log2Expression.txt", sep = ""), sep = "\t", row.names = FALSE)

NormalTumor_Comb2 <- NormalTumor_Comb # rename the data frame
rownames(NormalTumor_Comb2) <- NormalTumor_Comb2$GeneSymbol
NormalTumor_Comb2 <- subset(NormalTumor_Comb2, select = -c(GeneSymbol, EntrezID))

AnnotationDF_N <- data.frame(colnames(NormalTissueData_Gene2) [-c(1,2)]) # create an annotation data frame
AnnotationDF_N$TissueType <- "Normal"
colnames(AnnotationDF_N)[1] <- "PatientID"
AnnotationDF_T <- data.frame(colnames(TumorData_Genes2)[-c(1,2)])
AnnotationDF_T$TissueType <- "Tumor"
colnames(AnnotationDF_T)[1] <- "PatientID"
AnnotationDF <- rbind(AnnotationDF_N, AnnotationDF_T)
AnnotationDF2 <- subset(AnnotationDF, select = c("TissueType"))
rownames(AnnotationDF2) <- AnnotationDF$PatientID
AnnotationColors <- list(TissueType = c(Tumor = "darkorchid3", Normal = "orange1"))

pdf(paste(outputdir, SetName, "_Log2Expression_TCGApatients_TumorVsNormal.pdf"), width = 12, height = 8, onefile = FALSE)
pheatmap(as.matrix(NormalTumor_Comb2),
         color = colorRampPalette(c("green", "black", "red")) ( 100 ),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         annotation_col = AnnotationDF2,
         annotation_colors = AnnotationColors,
         fontsize_col = 4,
         fontsize_row = 1,
         cellheight = 2,
         cellwidth = .5)
dev.off()

GenesToRemove <- read.table("GenesToRemove_TCGA_TumorVsNormal.txt", sep = "\t", header = FALSE)[, c(1)]
NormalTumor_Comb3 <- NormalTumor_Comb2[!(row.names(NormalTumor_Comb2) %in% GenesToRemove), ]

pdf(paste(outputdir, SetName, "_Log2Expression_TCGApatients_TumorVsNormal_AfterGeneRemoval.pdf"), width = 12, height = 8, onefile = FALSE)
pheatmap(as.matrix(NormalTumor_Comb3),
         color = colorRampPalette(c("green", "black", "red")) ( 100 ),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         annotation_col = AnnotationDF2,
         annotation_colors = AnnotationColors,
         fontsize_col = 4,
         fontsize_row = 1,
         cellheight = 2,
         cellwidth = .5)
dev.off()


NormalTumor_Comb4 <- NormalTumor_Comb3 - apply(NormalTumor_Comb3, 1, median)

pdf(paste(outputdir, SetName, "_Log2Expression_TCGApatients_TumorVsNormal_MedianNorm.pdf"), width = 12, height = 8, onefile = FALSE)
pheatmap(as.matrix(NormalTumor_Comb4),
         color = colorRampPalette(c("green", "black", "red")) ( 100 ),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         annotation_col = AnnotationDF2,
         annotation_colors = AnnotationColors,
         fontsize_col = 4,
         fontsize_row = 1,
         cellheight = 2,
         cellwidth = .5)
dev.off()

NormalTumor_Comb5 <- as.data.frame(t(NormalTumor_Comb3)) # match patients and calculate log fold change
NormalTumor_Comb5$Patient_Tissue_ID <- rownames(NormalTumor_Comb5)
NormalTumor_Comb6 <- separate(data = NormalTumor_Comb5, col = Patient_Tissue_ID, into = c("Tissue", "PatientID"), sep = "\\_")
NormalTumor_Comb7 <- melt(NormalTumor_Comb6, id.vars = c("PatientID", "Tissue"), variable.name = "GeneSymbol", value.name = "Expression")
NormalTumor_Comb8 <- dcast(NormalTumor_Comb7, PatientID + GeneSymbol ~ Tissue, value.var = "Expression")
NormalTumor_Comb9 <- na.omit(NormalTumor_Comb8)
NormalTumor_Comb9$TumorVsNormal <- NormalTumor_Comb9$T - NormalTumor_Comb9$N
NormalTumor_Comb10 <- subset(NormalTumor_Comb9, select = -c(N, T))
NormalTumor_Comb11 <- dcast(NormalTumor_Comb10, GeneSymbol ~ PatientID, value.var = "TumorVsNormal")
rownames(NormalTumor_Comb11) <- NormalTumor_Comb11$GeneSymbol
NormalTumor_Comb12 <- subset(NormalTumor_Comb11, select = -c(GeneSymbol))

bk = unique(c(seq(min(NormalTumor_Comb12),-0.001, length=50), 0, seq(0.0001, max(NormalTumor_Comb12), length=50)))
pdf(paste(outputdir, SetName, "_Log2Expression_TCGApatients_TumorVsNormal_LogFoldChangeMatchedPatients.pdf"), width = 6, height = 8, onefile = FALSE)
pheatmap(as.matrix(NormalTumor_Comb12),
         color = colorRampPalette(c("green", "black", "red")) ( 100 ),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         breaks = bk, 
         fontsize_col = 4,
         fontsize_row = 1,
         cellheight = 2,
         cellwidth = 2)
dev.off()

# Correlation cell line models with cell line data

## All tumors

dir.create(paste(parent.dir, "TCGApatientVs2DVs3D", sep="/")) # create a new folder to put the files in
outputdir <- paste(parent.dir, "TCGApatientVs2DVs3D/", sep="/")

RNAseq_2Dvs3D <- read_tsv(file.path(shared.files.dir, "NormalizedCounts_2Dvs3D.txt.xz"))
TCGA_ReceptorStatus <- read.table("TCGApatient_ReceptorStatus.txt", sep = "\t", header = TRUE)
RNAseq_TCGA <- read_tsv("TCGARNAseq_LogNormalized.txt.xz")

RNAseq_TCGA2 <- RNAseq_TCGA # Adjust TCGA data frame
RNAseq_TCGA2$PatientID <- substr(RNAseq_TCGA2$PatientID, 1, 12)
RNAseq_TCGA2$PatientID <- gsub("\\.", "-", RNAseq_TCGA2$PatientID)
rownames(RNAseq_TCGA2) <- RNAseq_TCGA2$PatientID
RNAseq_TCGA2$PatientID <- NULL
RNAseq_TCGA3 <- as.data.frame(t(RNAseq_TCGA2))
RNAseq_TCGA3$entrez <- gsub("X", "", rownames(RNAseq_TCGA3))

GenesToRemove_2Dvs3D <- read.table("GenesToRemove_2Dvs3D.txt", sep = "\t", header= FALSE) [,1] # Select genes from the data frames and merge them
GenesToRemove_TCGA <- read.table("GenesToRemove_TCGA_Subtype.txt", sep = "\t", header = FALSE) [,1]
GenesToRemove <- intersect(GenesToRemove_2Dvs3D, GenesToRemove_TCGA)
GenesIDs_Sel <- GenesDF[!GenesDF$GeneSymbol %in% GenesToRemove, ]$EntrezID

RNAseq_TCGA_gene <- RNAseq_TCGA3[RNAseq_TCGA3$entrez %in% GenesIDs_Sel, ]
TCGA_CellLine_Comb <- merge(RNAseq_TCGA_gene, RNAseq_2Dvs3D)

CorDF <- cor(TCGA_CellLine_Comb[, 2:1098], TCGA_CellLine_Comb[, 1099:1126], method = "spearman") # Calculate spearman correlations and create a clustering
colnames(CorDF) <- gsub("_A", "", colnames(CorDF))
CorDF <- CorDF[, c("BT474_2D", "BT474_3D", "MDAMB175_2D", "MDAMB175_3D", "SKBR3_2D", "SKBR3_3D", 
                       "T47D_2D", "T47D_3D", "UACC893_2D", "UACC893_3D", "ZR751_2D", "ZR751_3D", 
                       "HCC1143_2D", "HCC1143_3D", "HCC1806_2D", "HCC1806_3D", "HCC1954_2D", "HCC1954_3D",
                       "MDAMB468_2D", "MDAMB468_3D",
                       "BT549_2D", "BT549_3D", "Hs578T_2D", "Hs578T_3D",
                       "MDAMB231_2D", "MDAMB231_3D", "SUM149PT_2D", "SUM149PT_3D")]

MetaDataDF <- read_tsv(file.path(shared.files.dir, "MetaData_2Dvs3D.txt"))
MetaDataDF <- unique(subset(MetaDataDF, select = c("CellLine", "X2D.3D", "Subtype")))
MetaDataDF$CellLine <- gsub("-", "", MetaDataDF$CellLine)
MetaDataDF$SampleID <- paste(MetaDataDF$CellLine, MetaDataDF$X2D.3D, sep = "_")
AnnotationDF_col <- subset(MetaDataDF, select = c("Subtype", "X2D.3D"))
rownames(AnnotationDF_col) <- MetaDataDF$SampleID 
AnnotationDF_col$X2D.3D <- paste("X", AnnotationDF_col$X2D.3D, sep = "")

AnnotationDF_row <- subset(TCGA_ReceptorStatus, select = -c(PatientID))
rownames(AnnotationDF_row) <- TCGA_ReceptorStatus$PatientID

AnnotationDF_row_write <- AnnotationDF_row
AnnotationDF_row_write$PatientID <- rownames(AnnotationDF_row_write)
write.table(AnnotationDF_row_write, paste(outputdir, "TCGApatient_ReceptorStatus.txt", sep = ""), sep = "\t", row.names = FALSE)

AnnotationColors <- list(ERstatus = c(Positive = "black", Negative="grey", Unknown = "white"),
                         PRstatus = c(Positive = "black", Negative="grey", Unknown = "white"),
                         HER2status = c(Positive = "black", Negative="grey", Unknown = "white"),
                         Subtype = c(Luminal = "darkorchid3", BasalA = "orange1", BasalB = "dodgerblue3"),
                         X2D.3D = c(X2D = "palegreen4", X3D = "palegreen"))

pdf(paste(outputdir, SetName, "_TCGA_CellLine_CorrelationHeatmap_AllTumors.pdf"), width = 8, height = 10, onefile = FALSE)
pheatmap(as.matrix(CorDF),
         color = colorRampPalette(c("green", "black", "red")) ( 100 ),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = TRUE,
         show_rownames = FALSE,
         annotation_row = AnnotationDF_row, 
         annotation_col = AnnotationDF_col, 
         annotation_colors = AnnotationColors,
         fontsize_col = 4,
         fontsize_row = 1,
         cellheight = 0.5,
         cellwidth = 10)
dev.off()

CorDF_write <- as.data.frame(CorDF)
CorDF_write$PatientID <- rownames(CorDF_write)
write.table(CorDF_write, file = paste(outputdir, "TCGA_CellLine_Correlation_AllTumors.txt", sep = ""), sep = "\t", row.names = FALSE)

CorDF2 <- as.data.frame(CorDF) # Create boxplots 
CorDF2$PatientID <- rownames(CorDF2)
CorDF3 <- melt(CorDF2, id.vars = c("PatientID"), variable.name = "SampleID", value.name = "Correlation")
CorDF4 <- separate(CorDF3, col = SampleID, into = c("CellLine", "X2D.3D"))
CorDF5 <- merge(CorDF4, MetaDataDF, by.all = TRUE)

CellLines <- c("BT474", "MDAMB175", "SKBR3", "UACC893", "T47D", "ZR751",
               "HCC1143", "HCC1806", "HCC1954", "MDAMB468",
               "BT549", "Hs578T", "MDAMB231", "SUM149PT")

pdf(paste(outputdir, SetName, "_TCGA_CellLine_Correlation_Boxplot_AllTumors.pdf"), width = 5, height = 4, onefile = FALSE)
ggplot(data = CorDF5, aes(x = CellLine, y = Correlation)) +
  geom_boxplot(aes(fill = X2D.3D, color = Subtype)) +
  scale_fill_manual(values = c("black", "grey")) +
  scale_color_manual(values = c("red3", "dodgerblue3", "palegreen4")) +
  scale_x_discrete(limits = CellLines) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Cell line") +
  theme_bw() +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()

CorDF6 <- CorDF5
CorDF7 <- dcast(CorDF6, CellLine + PatientID + Subtype  ~ X2D.3D, value.var = "Correlation")
CorDF7$Cor_3Dmin2D <- CorDF7$`3D` - CorDF7$`2D` 

pdf(paste(outputdir, SetName, "_TCGA_CellLine_Correlation_Boxplot_Delta2Dvs3D_AllTumors.pdf"), width = 5, height = 4, onefile = FALSE)
ggplot(data = CorDF7, aes(x = CellLine, y = Cor_3Dmin2D)) +
  geom_boxplot(aes(fill = Subtype)) +
  scale_fill_manual(values = c("red3", "dodgerblue3", "palegreen4")) +
  scale_x_discrete(limits = CellLines) +
  scale_y_continuous(expand = c(0,0)) +
  geom_abline(slope = 0) +
  xlab("Cell line") +
  ylab("Delta Correlation (3D vs 2D)") +
  theme_bw() +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

## Only for TNBC

TNBC_patients <- rownames(subset(AnnotationDF_row, ERstatus == "Negative" & PRstatus == "Negative" & HER2status == "Negative"))
CorDF <- as.data.frame(CorDF)
CorDF_TNBC <- CorDF[c(TNBC_patients), ]

pdf(paste(outputdir, SetName, "_TCGA_CellLine_CorrelationHeatmap_TNBC.pdf"), width = 8, height = 10, onefile = FALSE)
pheatmap(as.matrix(CorDF_TNBC),
         color = colorRampPalette(c("green", "black", "red")) ( 100 ),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = TRUE,
         show_rownames = FALSE,
         annotation_row = AnnotationDF_row, 
         annotation_col = AnnotationDF_col, 
         annotation_colors = AnnotationColors,
         fontsize_col = 4,
         fontsize_row = 1,
         cellheight = 2,
         cellwidth = 10)
dev.off()

CorDF5_TNBC <- CorDF5[CorDF5$PatientID %in% TNBC_patients, ]

pdf(paste(outputdir, SetName, "_TCGA_CellLine_Correlation_Boxplot_TNBC.pdf"), width = 5, height = 4, onefile = FALSE)
ggplot(data = CorDF5_TNBC, aes(x = CellLine, y = Correlation)) +
  geom_boxplot(aes(fill = X2D.3D, color = Subtype)) +
  scale_fill_manual(values = c("black", "grey")) +
  scale_color_manual(values = c("red3", "dodgerblue3", "palegreen4")) +
  scale_x_discrete(limits = CellLines) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Cell line") +
  theme_bw() +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

CorDF7_TNBC <- CorDF7[CorDF7$PatientID %in% TNBC_patients, ]

pdf(paste(outputdir, SetName, "_TCGA_CellLine_Correlation_Boxplot_Delta2Dvs3D_TNBC.pdf"), width = 5, height = 4, onefile = FALSE)
ggplot(data = CorDF7_TNBC, aes(x = CellLine, y = Cor_3Dmin2D)) +
  geom_boxplot(aes(fill = Subtype)) +
  scale_fill_manual(values = c("red3", "dodgerblue3", "palegreen4")) +
  scale_x_discrete(limits = CellLines) +
  scale_y_continuous(expand = c(0,0)) +
  geom_abline(slope = 0) +
  xlab("Cell line") +
  ylab("Delta Correlation (3D vs 2D)") +
  theme_bw() +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

## ER positive tumors

ERpos_patients <- rownames(subset(AnnotationDF_row, ERstatus == "Positive"))
CorDF <- as.data.frame(CorDF)
CorDF_ERpos <- CorDF[c(ERpos_patients), ]
CorDF_ERpos <- na.omit(CorDF_ERpos)

pdf(paste(outputdir, SetName, "_TCGA_CellLine_CorrelationHeatmap_ERpos.pdf"), width = 8, height = 10, onefile = FALSE)
pheatmap(as.matrix(CorDF_ERpos),
         color = colorRampPalette(c("green", "black", "red")) ( 100 ),
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = TRUE,
         show_rownames = FALSE,
         annotation_row = AnnotationDF_row, 
         annotation_col = AnnotationDF_col, 
         annotation_colors = AnnotationColors,
         fontsize_col = 4,
         fontsize_row = 1,
         cellheight = 0.5,
         cellwidth = 10)
dev.off()

CorDF5_ERpos <- CorDF5[CorDF5$PatientID %in% ERpos_patients, ]

pdf(paste(outputdir, SetName, "_TCGA_CellLine_Correlation_Boxplot_ERpos.pdf"), width = 5, height = 4, onefile = FALSE)
ggplot(data = CorDF5_ERpos, aes(x = CellLine, y = Correlation)) +
  geom_boxplot(aes(fill = X2D.3D, color = Subtype)) +
  scale_fill_manual(values = c("black", "grey")) +
  scale_color_manual(values = c("red3", "dodgerblue3", "palegreen4")) +
  scale_x_discrete(limits = CellLines) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Cell line") +
  theme_bw() +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

CorDF7_ERpos <- CorDF7[CorDF7$PatientID %in% ERpos_patients, ]

pdf(paste(outputdir, SetName, "_TCGA_CellLine_Correlation_Boxplot_Delta2Dvs3D_ERpos.pdf"), width = 5, height = 4, onefile = FALSE)
ggplot(data = CorDF7_ERpos, aes(x = CellLine, y = Cor_3Dmin2D)) +
  geom_boxplot(aes(fill = Subtype)) +
  scale_fill_manual(values = c("red3", "dodgerblue3", "palegreen4")) +
  scale_x_discrete(limits = CellLines) +
  scale_y_continuous(expand = c(0,0)) +
  geom_abline(slope = 0) +
  xlab("Cell line") +
  ylab("Delta Correlation (3D vs 2D)") +
  theme_bw() +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()



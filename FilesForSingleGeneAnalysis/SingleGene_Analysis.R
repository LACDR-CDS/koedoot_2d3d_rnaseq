library("ggplot2")
library("reshape2")
library("readr")

# To avoid data.frame vs tibble issues (eg with deprecated rownames)
read_tsv <- function(...) {
  as.data.frame(suppressMessages(readr::read_tsv(...)))
}

parent.dir <- "Add working directory" 
setwd(parent.dir)
shared.files.dir <- file.path("..", "FilesShared") # Set location of FilesShared folder here

SingleGeneDF <- read.table("SingleGene.txt", sep = "\t", header = TRUE)
GeneSymbol <- SingleGeneDF$GeneSymbol
GeneID <- SingleGeneDF$EntrezID
EnsemblID <- as.character(SingleGeneDF$Ensembl)

# Expression in different cell lines (subtypes)

dir.create(paste(parent.dir, "CellLine2D", sep="/")) # create a new folder to put the files in
outputdir <- paste(parent.dir, "CellLine2D/", sep="/")

CellLineExpression_2D <- read_tsv(file.path(shared.files.dir, "LogNormalizedCounts_2DCellLines.txt.xz")) # Load expression data
rownames(CellLineExpression_2D) <- CellLineExpression_2D$EnsemblID # Put the ensemblID in the rownames
CellLineExpression_2D$EnsemblID <- NULL # Remove the ensemblID column now they are in the rownames
CellLineExpression_2D <- as.data.frame(t(CellLineExpression_2D)) # transform the data frame so the cell lines are in the row names
CellLineExpression_2D$CellLine <- rownames(CellLineExpression_2D) # Create a cell line column that can be used to merge data frames

CellLine_Subtype <- read_tsv(file.path(shared.files.dir,"Subtype_CellLines.txt")) # Load the file containing cell line subtype information

Subtype_Expression_DF <- merge(CellLine_Subtype, CellLineExpression_2D)

GeneDF <- subset(Subtype_Expression_DF, select = c("CellLine", "Subtype", EnsemblID))

write.table(GeneDF, file = paste(outputdir, GeneSymbol, "_Log2Expression.txt", sep = ""), sep = "\t", row.names = FALSE)

pdf(paste(outputdir, GeneSymbol, "_CountsPerCellLine.pdf", sep = ""), width = 9, height = 6)
ggplot(data = Subtype_Expression_DF, aes(x = CellLine, y = GeneDF[,3])) +
  geom_bar(aes(fill = Subtype), stat = "identity") +
  ylab("Log2 Expression Level") +
  xlab("Cell line") +
  ggtitle(GeneSymbol) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("red3", "dodgerblue3", "green4")) +
  theme_bw() +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
        plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste(outputdir, GeneSymbol, "_CountsPerCellLine_Sorted.pdf", sep = ""), width = 9, height = 6)
ggplot(data = Subtype_Expression_DF, aes(x = reorder(CellLine, -GeneDF[,3]), y = GeneDF[,3])) +
  geom_bar(aes(fill = Subtype), stat = "identity") +
  ylab("Log2 Expression Level") +
  xlab("Cell line") +
  ggtitle(GeneSymbol) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("red3", "dodgerblue3", "green4")) +
  theme_bw() +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
        plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste(outputdir, GeneSymbol, "_CountsPerSubtype.pdf", sep = ""), width = 5, height = 5)
ggplot(data = Subtype_Expression_DF, aes(x = Subtype, y = GeneDF[,3])) +
  geom_boxplot(aes(fill = Subtype), alpha = 0.7) +
  geom_jitter(shape = 21, aes(fill = Subtype)) +
  ylab("Log2 Expression Level") +
  xlab("Cell line") +
  ggtitle(GeneSymbol) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("red3", "dodgerblue3", "green4")) +
  theme_bw() +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5))
dev.off()


# Expression 2D vs 3D

dir.create(paste(parent.dir, "CellLine2Dvs3D", sep="/")) # create a new folder to put the files in
outputdir <- paste(parent.dir, "CellLine2Dvs3D/", sep="/")

NormalizedCounts <- read_tsv(file.path(shared.files.dir, "NormalizedCounts_2Dvs3D.txt.xz"))
MetaDataDF <- read_tsv(file.path(shared.files.dir, "MetaData_2Dvs3D.txt"))

NormalizedCounts_Long <- melt(NormalizedCounts, id.vars = c("symbol", "entrez", "ensembl"), 
                                       variable.name = "Sample", value.name = "Expression") # create a long data frame
NormalizedCounts_Long_Gene <- subset(NormalizedCounts_Long, entrez == GeneID) # select gene of interest
NormalizedCounts_Long_Gene2 <- merge(NormalizedCounts_Long_Gene, MetaDataDF, by.x = "Sample", by.y = "YourSampleID2") # merge with metadata


CellLines <- c("BT474", "MDA-MB-175", "SK-BR-3", "UACC893", "T47D", "ZR75-1",
               "HCC1143", "HCC1806", "HCC1954", "MDA-MB-468",
               "BT549", "Hs578T", "MDA-MB-231", "SUM149PT")

pdf(file = paste(outputdir, GeneSymbol, "_2D3D_Log2Counts.pdf", sep = ""), width = 5, height = 4)
ggplot(data = NormalizedCounts_Long_Gene2, aes(x = CellLine, y = Expression)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = X2D.3D, color = Subtype), size = 0.5) +
  scale_fill_manual(values = c("black", "grey")) +
  scale_color_manual(values = c("red2", "dodgerblue3", "springgreen4")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(limits = CellLines) +
  xlab("Cell line") +
  ylab("Log2 Normalized Read Counts") +
  theme_bw() +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()


NormalizedCounts_Gene <- dcast(NormalizedCounts_Long_Gene2, symbol + entrez + ensembl + CellLine + Subtype ~ X2D.3D, value.var = "Expression") # select the columns of interest and create a 2D/3D column
colnames(NormalizedCounts_Gene)[6:7] <- c("X2D", "X3D") # change the column names so that you don't start with a number
NormalizedCounts_Gene$FoldChange2Dvs3D <- NormalizedCounts_Gene$X2D - NormalizedCounts_Gene$X3D # calculate the fold change between 2D and 3D
NormalizedCounts_Gene$MeanExpr <- (NormalizedCounts_Gene$X2D + NormalizedCounts_Gene$X3D)/2 # calculate the mean expression of a gene over the 2D and 3D cell lines

pdf(file = paste(outputdir, GeneSymbol, "_FoldChange2Dvs3D.pdf", sep = ""), width = 4, height = 4)
ggplot(data = NormalizedCounts_Gene, aes(x = reorder(CellLine, -FoldChange2Dvs3D), y = FoldChange2Dvs3D)) +
  geom_bar(stat = "identity", aes(fill = Subtype)) +
  scale_fill_manual(values = c("red2", "dodgerblue3", "springgreen4")) +
  xlab("Cell line") +
  coord_cartesian(ylim = c(-max(abs(NormalizedCounts_Gene$FoldChange2Dvs3D)),
                           max(abs(NormalizedCounts_Gene$FoldChange2Dvs3D)))) +
  geom_abline(slope = 0) +
  theme_bw() +
  theme(axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()

write.table(NormalizedCounts_Gene, file = paste(outputdir, GeneSymbol, "_FoldChange2Dvs3D.txt", sep = ""), sep = "\t", row.names = FALSE)

# Expression in different subtypes TCGA

dir.create(paste(parent.dir, "TCGAsubtypes", sep="/")) # create a new folder to put the files in
outputdir <- paste(parent.dir, "TCGAsubtypes/", sep="/")

TCGA_Subtype_DF <- read_tsv(file.path(shared.files.dir, "TCGA_PrimaryTumor_Subtype_Expression.txt.xz")) # Load TCGA data

GeneDF <- subset(TCGA_Subtype_DF, select = c("TumorType", paste("ID", GeneID, sep = "_")))
colnames(GeneDF) <- c("TumorType", "ExpressionLevel")

pdf(file = paste(outputdir, GeneSymbol, "_TNBCvsERpos_TCGAprimarytumors.pdf", sep = ""), width = 3, height = 4)
ggplot(data = GeneDF, aes(x = TumorType, y = ExpressionLevel)) +
  geom_boxplot(aes(fill = TumorType)) +
  ylab("Log2 Expression level") +
  xlab("Tumor type") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("firebrick2", "dodgerblue3")) +
  ggtitle(GeneSymbol) +
  guides(fill = FALSE) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(hjust = 0.5))
dev.off()

write.table(GeneDF, file = paste(outputdir, GeneSymbol, "_TNBCvsERpos_TCGAprimarytumors.txt", sep = ""), row.names = FALSE, sep = "\t")

# Expression tumor vs normal TCGA

dir.create(paste(parent.dir, "TCGA_TumorVsNormal", sep="/")) # create a new folder to put the files in
outputdir <- paste(parent.dir, "TCGA_TumorVsNormal/", sep="/")

TumorData <- read_tsv(file.path(shared.files.dir, "BRCA__illuminahiseq_rnaseqv2__GeneExp_PrimarySolidTumor.txt.xz"))
colnames(TumorData) <- c("GeneSymbol", "EntrezID", gsub(".{16}$", "", colnames(TumorData)[-c(1,2)])) # change the column names to patient IDs by removing the last 16 characters
TumorPatients <- colnames(TumorData)[-c(1,2)] # vector of patients in tumor data
TumorData_Gene <- subset(TumorData, EntrezID == GeneID) # Select gene of interest

NormalTissueData <- read_tsv(file.path(shared.files.dir, "BRCA__illuminahiseq_rnaseqv2__GeneExp_SolidTissueNormal.txt.xz"))
colnames(NormalTissueData) <- c("GeneSymbol", "EntrezID", gsub(".{16}$", "", colnames(NormalTissueData)[-c(1,2)])) # change the column names to patient IDs by removing the last 16 characters
NormalTissuePatients <- colnames(NormalTissueData)[-c(1,2)] # vector of patients in normal data
NormalTissueData_Gene <- subset(NormalTissueData, EntrezID == GeneID) # Select gene of interest

TumorData_Long <- melt(TumorData_Gene, id.vars = c("GeneSymbol", "EntrezID"), variable.name = "PatientID", value.name = "Expression") # Create a long data frame
TumorData_Long$Expression <- log2(TumorData_Long$Expression+1)

NormalTissueData_Long <- melt(NormalTissueData_Gene, id.vars = c("GeneSymbol", "EntrezID"), variable.name = "PatientID", value.name = "Expression")
NormalTissueData_Long$Expression <- log2(NormalTissueData_Long$Expression + 1)

## Matched

TumorData_Long2 <- TumorData_Long # Rename the data frame
colnames(TumorData_Long2) [4] <- "Tumor"

NormalTissueData_Long2 <- NormalTissueData_Long
colnames(NormalTissueData_Long2) [4] <- "Normal"

TumorVsNormal_DF <- merge(TumorData_Long2, NormalTissueData_Long2, by.all = TRUE) # Combine the 2 data frames
TumorVsNormal_DF$TumorMinNormal <- TumorVsNormal_DF$Tumor - TumorVsNormal_DF$Normal

pdf(file = paste(outputdir, GeneSymbol, "_PrimaryTumorvsNormalTissue_Matched_DensityPlot.pdf", sep = ""), width = 4, height = 4)
ggplot(data = TumorVsNormal_DF, aes (x = TumorMinNormal)) +
  geom_density(fill = "grey80") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  ylab("Density") +
  xlab("Log2 Primary Tumor vs Normal Tissue") +
  ggtitle (GeneSymbol) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=10, color = "black"),
        axis.title=element_text(size=11,face="bold"),
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(hjust = 0.5))
dev.off()

write.table(TumorVsNormal_DF, file = paste(outputdir, GeneSymbol, "_PrimaryTumorvsNormalTissue_Matched.txt", sep = ""), sep = "\t", row.names = FALSE)

TumorVsNormal_DF2 <- subset(TumorVsNormal_DF, select = -c(TumorMinNormal)) # remove the delta column
TumorVsNormal_DF3 <- melt(TumorVsNormal_DF2, id.vars = c("GeneSymbol", "EntrezID", "PatientID"), variable.name = "Tumor_Normal", value.name = "Expression") # Create a long data frame

tumor_normal <- c("Normal", "Tumor")

pdf(file = paste(outputdir, GeneSymbol, "_PrimaryTumorvsNormalTissue_Matched_Boxplot.pdf", sep = ""), width = 3, height = 4)
ggplot(data = TumorVsNormal_DF3, aes(x = Tumor_Normal, y = Expression)) +
  geom_boxplot(aes(fill = Tumor_Normal)) +
  ylab("Log2 Expression level") +
  xlab("Tissue type") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("dodgerblue3", "firebrick2")) +
  scale_x_discrete(limits = tumor_normal) +
  ggtitle(GeneSymbol) +
  guides(fill = FALSE) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=10, color = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(hjust = 0.5))
dev.off()

## Non-matched

TumorData_Long3 <- TumorData_Long # Rename the tumor data frame
TumorData_Long3$Type <- "Tumor"

NormalTissueData_Long3 <- NormalTissueData_Long
NormalTissueData_Long3$Type <- "Normal"

TumorNormalComb <- rbind(NormalTissueData_Long3, TumorData_Long3)

write.table(TumorNormalComb, file = paste(outputdir, GeneSymbol, "_PrimaryTumorvsNormalTissue_NonMatched.txt", sep = ""), sep ="\t", row.names = FALSE)

pdf(file = paste(outputdir, GeneSymbol, "_PrimaryTumorvsNormalTissue_NonMatched_Boxplot.pdf", sep = ""), width = 3, height = 4)
ggplot(data = TumorNormalComb, aes(x = Type, y = Expression)) +
  geom_boxplot(aes(fill = Type)) +
  ylab("Log2 Expression level") +
  xlab("Tissue type") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("firebrick2", "dodgerblue3")) +
  ggtitle(GeneSymbol) +
  guides(fill = FALSE) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=10, color = "black"),
        axis.title=element_text(size=12,face="bold"),
        axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
        plot.title = element_text(hjust = 0.5))
dev.off()






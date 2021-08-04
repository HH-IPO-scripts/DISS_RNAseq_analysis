#Skeleton RNA-seq pipeline provided by Simon Tomlinson and modified to fit the needs of my diss project

#loading the required libraries l
library(Rsubread)
library(limma)
library(edgeR)
library(gplots)
library(DESeq2)
library(affy)
library(QuasR)
library(annotate)

#generating input files
#for files actually ending in ".bam"
alignment_files <- list.files(path = './sorted_bam', pattern = '*.bam',full.names = T)

#for renamed files not ending in ".bam"
#alignment_files <- list.files(path = './sorted_bam_2', full.names = T)

#use feature counts to load the BAM files 
# MQS = 30
my_table_features_gene <- featureCounts(files = alignment_files, annot.ext="./annotation/genome.gtf",isGTFAnnotationFile=TRUE, isPairedEnd=TRUE, GTF.featureType="gene", countMultiMappingReads = FALSE, minMQS = 30, nthreads = 30, minOverlap = 10)
# MQS = 15
# my_table_features_gene <- featureCounts(files = alignment_files, annot.ext="./annotation/genome.gtf",isGTFAnnotationFile=TRUE, GTF.featureType="gene", countMultiMappingReads = FALSE, minMQS = 15, nthreads = 10, minOverlap = 10)
# MQS = 0 (default)
#my_table_features_gene <- featureCounts(files = alignment_files, annot.ext="./annotation/genome.gtf",isGTFAnnotationFile=TRUE, GTF.featureType="gene", countMultiMappingReads = FALSE, minMQS = 0, nthreads = 10, minOverlap = 10)


#Some further processing to access count data and column isnames.ul 
table_rnaseq_gene <- (my_table_features_gene)$counts
head(table_rnaseq_gene)
colnames_rnaseq_gene <- colnames(table_rnaseq_gene)

#making annotated dataframe
#Simons way that didn't work when wanting colours, but works for DESEQ
adf <-read.table("samples.txt",sep='\t', row.names=1,fill=T,header=T)
rownames(adf) <- gsub("_", ".", rownames(adf))
rownames(adf) <- gsub("-", ".", rownames(adf))

#order the elements in the table
idx<-match(colnames_rnaseq_gene,rownames(adf))
colnames(table_rnaseq_gene) <-adf$ShortName[idx]
adf <-adf[idx,]
rownames(adf) <-colnames(table_rnaseq_gene)
adf

#create object for DESeq
#NOTE the design for differential analysis is fed into this object- the adf annotates each sample
#need to change adf to allow for DESEQ to work
dds_rnaseq <- DESeqDataSetFromMatrix(countData = table_rnaseq_gene, colData = adf, design = ~ Type)

#check the object dimensions
dim(dds_rnaseq)
head(rownames(dds_rnaseq))
colnames(dds_rnaseq)

#check the object type
class(dds_rnaseq)
typeof(dds_rnaseq)

#check what this object contains
slotNames(dds_rnaseq)
colnames(colData(dds_rnaseq))
colnames(assay(dds_rnaseq))

#apply filters that remove genes that do not pass quality filters - this is optional. 
#first filter based upon counts
keep <- rowSums(counts(dds_rnaseq)) > 1
dds_rnaseq2 <- dds_rnaseq[keep,]
nrow(dds_rnaseq2)

#second filter-require >10 in >= 3 samples
keep <- rowSums(counts(dds_rnaseq2) >= 10) >= 3
dds_rnaseq3 <- dds_rnaseq2[keep,]
nrow(dds_rnaseq3)

#size factor estimation for normalization- these are used to scale by library size when data is exported.
dds_rnaseq_gene <- estimateSizeFactors(dds_rnaseq3)
sizeFactors(dds_rnaseq_gene)

#DEseq2 also provides methods that outputs data with various normalisations

#normalised log counts
counts_rnaseq <-log2(counts(dds_rnaseq_gene, normalized=TRUE))
#fpm is similar to spm. here the data is log transformed
fpm_rnaseq <-log(fpm(dds_rnaseq_gene)) #this is also normalised by the sizeFactors
vsd_rnaseq <- vst(dds_rnaseq_gene, blind = T)
head(vsd_rnaseq, 3)
rld_rnaseq <- rlog(dds_rnaseq_gene, blind = T)

#alternatively, if the dataset is larger
#rld_rnaseq <- vst(dds_rnaseq, blind = T)
head(rld_rnaseq, 3)

#making new adf that works to add colours to plots
adf1 <- read.AnnotatedDataFrame('samples.txt', header=TRUE, row.names = 1, as.is=TRUE)
rownames(adf1) <- gsub("_", ".", rownames(adf1))
rownames(adf1) <- gsub("-", ".", rownames(adf1))

idx<-match(colnames_rnaseq_gene,rownames(adf1))
colnames(table_rnaseq_gene) <-adf1$ShortName[idx]
adf1 <-adf1[idx,]
rownames(adf1) <-colnames(table_rnaseq_gene)
adf1

#rlog will be used from this point onwards, notes that differential expression uses counts -rlog is just for visualisation
# time to make some plots with the normalised data!
png('counts_rnaseq_boxplot_gene.png')
par(mar=c(7,5,3,3))
boxplot(counts_rnaseq, col = adf1$Colour, las=2, names=adf1$ShortName, horizontal=TRUE, main = "GSE125285 Gene-level RNAseq counts boxplot")
dev.off()

png('fpm_rnaseq_boxplot_gene.png')
par(mar=c(7,5,3,3))
boxplot(fpm_rnaseq, col=adf1$Colour, las=2, names=adf1$ShortName, horizontal=TRUE, main = "GSE125285 Gene-level fpm_rnaseq boxplot")
dev.off()

png('assay_vsd_rnaseq_boxplot_gene.png')
par(mar=c(7,5,3,3))
boxplot(assay(vsd_rnaseq),col=adf1$Colour, las=2, names=adf1$ShortName, horizontal=TRUE, main = "GSE125285 Gene-level vsd_rnaseq boxplot")
dev.off()

png('assay_rld_rnaseq_boxplot_gene.png')
par(mar=c(7,5,3,3))
boxplot(assay(rld_rnaseq),col=adf1$Colour, las=2, names=adf1$ShortName, horizontal=TRUE, main = "GSE125285 Gene-level rld_rnaseq boxplot")
dev.off()

#but suppose we look at non-normalised data
na.rm=T
counts_un <-counts(dds_rnaseq_gene, normalized=FALSE)
png('unnormalised_counts_boxplot_gene.png')
par(mar=c(7,5,3,3))
boxplot(log2(counts_un),col = adf1$Colour, las=2, names=adf1$ShortName, horizontal = TRUE, main = "GSE125285 Gene-level unnormalised counts boxplot")
dev.off()

counts_un_rlog<-rlog(counts_un)
#OR
#counts_un_rlog<-vst(counts_un)
png(filename="mva.pairs.png",width=2000, height=2000)
mva.pairs(counts_un_rlog, main = "GSE125285 Gene-level mva_pairs plots") # (figure margins are too large)
dev.off()

#making PCA plots -  loading required libraries
library(scatterplot3d)
library(ggplot2)

#perform PCA
pca <- prcomp(t(na.omit(assay(rld_rnaseq))), scale=T)

#Plot the PCA results 
png('3D_PCA_plot_gene.png')
s3d<-scatterplot3d(pca$x[,1:3], pch=19, color=colData(dds_rnaseq_gene)$Colour, main = "GSE125285 Gene-level 3D PCA plot")
s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
text(s3d.coords$x, s3d.coords$y, labels = colnames(rld_rnaseq),pos = 3,offset = 0.5,cex=0.5)
dev.off()

#we can also output a distance matrix
sampleDists <- dist(t(assay(rld_rnaseq)))
sampleDists

#Generating further heatmaps with the data
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
png('sample_distance_heatmap_gene.png')
pheatmap(sampleDistMatrix, main = "GSE125285 Gene-level sample similarity heatmap")
dev.off()

#easier PCA plot
png('Easier_PCA_plot_gene.png')
plotPCA(rld_rnaseq, intgroup = c("Type"))
dev.off()

#Differential Gene Expression
dds_rnaseq_gene <- DESeq(dds_rnaseq3)
head(dds_rnaseq_gene)

# Gene IDs for genes of interest to this study
genes <- c("ENSG00000067225", "ENSG00000143627", "ENSG00000067057", "ENSG00000134333", "ENSG00000160211", "ENSG00000130313", "ENSG00000142657", "ENSG00000153574", "ENSG00000197713", "ENSG00000163931", "ENSG00000177156", "ENSG00000138413", "ENSG00000105953", "ENSG00000137764", "ENSG00000110092")

#This is how you specify the contrast
#all results
results_genes <- results(dds_rnaseq_gene, contrast=c("Type","Cancer","Normal"), parallel = TRUE)

#Storing the results of the genes of interest in a single variable
myresults_3 <- results_genes[c("ENSG00000067225", "ENSG00000143627", "ENSG00000067057", "ENSG00000134333", "ENSG00000160211", "ENSG00000130313", "ENSG00000142657", "ENSG00000153574", "ENSG00000197713", "ENSG00000163931", "ENSG00000177156", "ENSG00000138413", "ENSG00000105953", "ENSG00000137764", "ENSG00000110092", "ENSG00000143627"),]

#Download Ensembl annotation using BiomaRt and rename the samples
library(biomaRt)
ensembl_host <-"uswest.ensembl.org"
head(biomaRt::listMarts(host = ensembl_host), 15)
head(biomaRt::listAttributes(biomaRt::useDataset(dataset = "hsapiens_gene_ensembl",mart = useMart("ENSEMBL_MART_ENSEMBL",host = ensembl_host))), 40)
mart <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl",mart = useMart("ENSEMBL_MART_ENSEMBL", host=ensembl_host))
ResultAnnot <- biomaRt::getBM(values=rownames(dds_rnaseq), attributes=c('ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','description','strand'), filters='ensembl_gene_id', mart=mart, useCache=F)

#merge with input data
names <- ResultAnnot[,1]
ResultAnnot <- as.data.frame(ResultAnnot)
rownames(ResultAnnot) = names
idx<-match(rownames(dds_rnaseq_gene),rownames(ResultAnnot))

#make sure annotation is in the same order
all(rownames(dds_rnaseq_gene) == rownames(ResultAnnot))
grr<-ResultAnnot[match(rownames(dds_rnaseq_gene), ResultAnnot$ensembl_gene_id),]
all(rownames(dds_rnaseq_gene) == rownames(grr))
ResultAnnot <-grr
all(rownames(dds_rnaseq_gene) == rownames(ResultAnnot))

#make the names nice
nice_names<- paste(ResultAnnot$ensembl_gene_id,ResultAnnot$external_gene_name, sep = '_')
ResultAnnot$nice_names <-nice_names
head(ResultAnnot)
all(rownames(dds_rnaseq_gene) == rownames(ResultAnnot))

#adding annotation to genes
idx <- match(genes, rownames(dds_rnaseq_gene))
genes <- ResultAnnot$nice_names[idx]
genes

#adding annotation to myresults (key genes)
idx2 <- match(rownames(myresults_3), rownames(dds_rnaseq_gene))
rownames(myresults_3) <- ResultAnnot$nice_names[idx2]
myresults_3

#adding annotation to all results
idx3 <- match(rownames(results_genes), rownames(dds_rnaseq_gene))
rownames(results_genes) <- ResultAnnot$nice_names[idx3]
head(results_genes)

#adding annotation to dds_rnaseq
rownames(dds_rnaseq3) <- ResultAnnot$nice_names

# simple practice plots for users viewing. 
#Plotting the PKM counts in cancer vs normal
plotCounts(dds_rnaseq, gene = "ENSG00000067225_PKM", intgroup="Type", col = adf1$Colour, main = "GSE87096 PKM cancer vs normal gene counts") 
#Plotting the PFKP counts in cancer vs normal
plotCounts(dds_rnaseq, gene = "ENSG00000067057_PFKP", intgroup="Type", col = adf1$Colour, main = "GSE87096 PFKP cancer vs normal gene counts") 
#Plotting LDHA counts in cancer vs normal
plotCounts(dds_rnaseq, gene = "ENSG00000134333_LDHA", intgroup="Type", col = adf1$Colour, main = "GSE87096 LDHA cancer vs normal gene counts") 

#plotting of read counts for all genes of interest in different sample types. 
# loop for plotting gene counts 
dev.new()
par(mfrow = c(3,5))
for(g in 1:length(genes))
  {
  plotCounts(dds_rnaseq, gene = genes[g], intgroup = "Type", normalized=FALSE, transform=FALSE, las=2, col = adf1$Colour, pch=19)
   }

# PLOT DA VOLCANO
library(EnhancedVolcano)
#Potting all results 
EnhancedVolcano(results_genes, lab = rownames(results_genes), x = "log2FoldChange", y = 'pvalue', title='GSE125285 Cancer vs Normal DESEQ results', labSize=2, pointSize=2.0)

#making a MA plot for differential gene expression
png('Main_DESeq2_plot_gene.png')
plotMA(results_genes, main="GSE125285 DESeq2 plot of unfiltered results", ylim=c(-2,2))
dev.off()

#Selection of results with adjusted p value < 0.05 and sorting in by Log2foldchange
table(results_genes$padj < 0.01)
table(results_genes$padj < 0.05)
results_selected <- subset(results_genes, padj < 0.05)
results_selected <- results_selected[order(abs(results_selected$log2FoldChange), decreasing=TRUE), ]
head(results_selected)
top50 <-rownames(results_selected)[1:50]

#making a MA plot for filtered differential gene expression
png('Main_DESeq2_plot_filtered_exon.png')
plotMA(results_selected, main="GSE125285 DESeq2 plot of filtered results", ylim=c(-2,2))
dev.off()

#Searching filtered featurecount results for genes of interest (gene level)
#filtered results
r2 <- data.frame(results_selected)
# Pyruvate Kinase M1/M2
r2[which(rownames(r2) == "ENSG00000067225_PKM"),]
# Phosphofructokinase (platelet)
r2[which(rownames(r2) == "ENSG00000067057_PFKP"),]
# LDHA
r2[which(rownames(r2) == "ENSG00000134333_LDHA"),]

#Download Ensembl annotation using BiomaRt and rename the samples
library(biomaRt)
ensembl_host <-"uswest.ensembl.org"
head(biomaRt::listMarts(host = ensembl_host), 15)
head(biomaRt::listAttributes(biomaRt::useDataset(dataset = "hsapiens_gene_ensembl",mart = useMart("ENSEMBL_MART_ENSEMBL",host = ensembl_host))), 40)
mart <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl",mart = useMart("ENSEMBL_MART_ENSEMBL", host=ensembl_host))
ResultAnnot <- biomaRt::getBM(values=rownames(dds_rnaseq), attributes=c('ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','description','strand'), filters='ensembl_gene_id', mart=mart, useCache=F)

#merge with input data
names <- ResultAnnot[,1]
ResultAnnot <- as.data.frame(ResultAnnot)
rownames(ResultAnnot) = names
idx<-match(rownames(dds_rnaseq),rownames(ResultAnnot))

#make sure annotation is in the same order
all(rownames(dds_rnaseq) == rownames(ResultAnnot))
grr<-ResultAnnot[match(rownames(dds_rnaseq), ResultAnnot$ensembl_gene_id),]
all(rownames(dds_rnaseq) == rownames(grr))
ResultAnnot <-grr
all(rownames(dds_rnaseq) == rownames(ResultAnnot))

#make the names nice
nice_names<- paste(ResultAnnot$ensembl_gene_id,ResultAnnot$external_gene_name, sep = '_')
ResultAnnot$nice_names <-nice_names
head(ResultAnnot)
all(rownames(dds_rnaseq) == rownames(ResultAnnot))

#check the names
rld_rnaseq <- vst(dds_rnaseq, blind = TRUE)
idxx <- match(rownames(results_selected)[1:50],rownames(dds_rnaseq)) #for top 50 differentially expressed genes
plotme <-(rld_rnaseq)[rownames(results_selected)[1:50],]
rownames(plotme)<-ResultAnnot$nice_names[idxx]

#make a heatmap with the candidate genes
#plotting the top 50 differentiated genes between the cancer Vs normal samples
png(filename="heatmap_top50_cancervsnormal_exon.png",width=800, height=1000)
pheatmap(assay(plotme),scale="row",fontsize_row = 10,cellheight =12, cellwidth =12,treeheight_row = 40, treeheight_col = 40, main = "GSE125285 top 50 differentiated genes between Cancer and Normal sample")
dev.off()


#Exon level featurecounts section!!
#PKM1/M2 alternative splicing 
# below are the partiular exon IDs that are to be used to distinguish PKM1 and PKM2 levels
# NCBI(PKM) transcipt was ENST00000335181.10
#PKM Exon 9 ID : ENSE00003618390
#PKM Exon 10 ID : ENSE00003632285

adf <-read.table("samples.txt",sep='\t', row.names=1,fill=T,header=T)

#exon level feauturecounts
mytable_features_exon <- featureCounts(files = alignment_files, annot.ext="./annotation/genome.gtf",isGTFAnnotationFile=TRUE, isPairedEnd=TRUE, GTF.featureType="exon", GTF.attrType="exon_id", useMetaFeatures = FALSE, countMultiMappingReads = TRUE, allowMultiOverlap = TRUE, minMQS = 0, nthreads = 20, minOverlap = 10)

#check featurecount stats
mytable_features_exon$stat

#check dimenstions of featurecounts
dim(mytable_features_exon$counts)

#filtering 
keep <- rowSums(mytable_features_exon$counts >= 10) >= 3
length(keep)
keep <- which(rowSums(mytable_features_exon$counts)>0)
length(keep)

#Some further processing to access count data and column isnames.ul
table_rnaseq_exon <- (mytable_features_exon)$counts
head(table_rnaseq_exon)
colnames_rnaseq_exon <- colnames(table_rnaseq_exon)

#Simons way that didn't work when wanting colours, but works for DESEQ
adf <-read.table("samples.txt",sep='\t', row.names=1,fill=T,header=T)
rownames(adf) <- gsub("_", ".", rownames(adf))
rownames(adf) <- gsub("-", ".", rownames(adf))
adf

#order the elements in the table
idx<-match(colnames_rnaseq_exon,rownames(adf))
colnames(table_rnaseq_exon) <-adf$ShortName[idx]
adf <-adf[idx,]
rownames(adf) <-colnames(table_rnaseq_exon)
adf

#create object for DESeq
#NOTE the design for differential analysis is fed into this object- the adf annotates each sample
#need to change adf to allow for DESEQ to work
dds_rnaseq <- DESeqDataSetFromMatrix(countData = table_rnaseq_exon, colData = adf, design = ~ Type)

#check the object dimensions
dim(dds_rnaseq)
head(rownames(dds_rnaseq))
colnames(dds_rnaseq)

#check the object type
class(dds_rnaseq)
typeof(dds_rnaseq)

#check what this object contains
slotNames(dds_rnaseq)
colnames(colData(dds_rnaseq))
colnames(assay(dds_rnaseq))

#check number of dds_rnaseq rows before filtering
nrow(dds_rnaseq)

#apply filters that remove genes that do not pass quality filters
#first filter based upon counts
keep <- rowSums(counts(dds_rnaseq)) > 1
dds_rnaseq1 <- dds_rnaseq[keep,]
nrow(dds_rnaseq1)

#second filter-require >10 in >= 3 samples
keep <- rowSums(counts(dds_rnaseq1) >= 10) >= 3
dds_rnaseq2 <- dds_rnaseq1[keep,]
nrow(dds_rnaseq2)

#size factor estimation for normalisation- these are used to scale by library size when data is exported.
dds_rnaseq3 <- estimateSizeFactors(dds_rnaseq1)
sizeFactors(dds_rnaseq3)

#DEseq2 also provides methods that outputs data with various normalisations
#normalised log counts
counts_rnaseq <-log2(counts(dds_rnaseq3, normalized=TRUE))
#fpm is similar to spm. here the data is log transformed
fpm_rnaseq <-log(fpm(dds_rnaseq3)) #this is also normalised by the sizeFactors
vsd_rnaseq <- vst(dds_rnaseq3, blind = T)
head(vsd_rnaseq, 3)
rld_rnaseq <- vst(dds_rnaseq3, blind = T)
head(rld_rnaseq, 3)

#making new adf that works to add colours to plots
adf1 <- read.AnnotatedDataFrame('samples.txt', header=TRUE, row.names = 1, as.is=TRUE)
rownames(adf1) <- gsub("_", ".", rownames(adf1))
rownames(adf1) <- gsub("-", ".", rownames(adf1))
adf1

idx<-match(colnames_rnaseq_exon,rownames(adf1))
colnames(table_rnaseq_exon) <-adf1$ShortName[idx]
adf1 <-adf1[idx,]
rownames(adf1) <-colnames(table_rnaseq_exon)
adf1

#Perform DESEQ, creating dds_rnaseq2
library("BiocParallel")
register(MulticoreParam(60))
dds_rnaseq3 <- DESeq(dds_rnaseq3, parallel=TRUE, minReplicatesForReplace=12)
dds_rnaseq3 #up to here

# Obtaining contrast results
# cancer relative to normal
# exon_results <- results(dds_rnaseq2, contrast=c("Type","Cancer","Normal"),alpha=alpha, parallel=TRUE)

# normal relative to cancer
# exon_results <- results(dds_rnaseq2, contrast=c("Type","Normal","Cancer"),alpha=alpha, parallel=TRUE)

#exons of interest 
these_exons <- c("ENSE00003618390","ENSE00003632285") # PKM1/M2

#get results of these particular exons
exon_results[match(these_exons, mytable_features_exon$annotation$GeneID),]
exon_results2 <- results(dds_rnaseq3, contrast=c("Type","Cancer","Normal"), parallel=TRUE)
exon_GB_normal_results <- results(dds_rnaseq3, contrast=c("Type","GB","Normal"), parallel=TRUE)
exon_GBnormal <- as.data.frame(exon_GB_normal_results)
exon_GBnormal[these_exons,]

#plot counts of particular exons
dev.new()
par(mfrow = c(1,2))
for(g in 1:length(these_exons))
{
  plotCounts(dds_rnaseq3, gene = these_exons[g], intgroup="Type", normalized=FALSE, transform=FALSE, ylim=c(0,2500), las=2, col = adf1$Colour ,pch=19)
}

#getting the number of read counts for each exon for comparison and ratio determination
exon_featurecounts <- as.data.frame(table_rnaseq_exon)
exon_featurecounts["ENSE00003618390",] #exon 9 counts
exon_featurecounts["ENSE00003632285",]# exon 10 counts

#Exons of interest 
exonsIwant <- c("ENSE00003618390","ENSE00003632285")

#getting index for these exons
myexonIndex <- match(exonsIwant, rownames(exon_featurecounts))
exon_data <- exon_featurecounts[myexonIndex,]

#make sure its the right data
head(exon_data)

#checking - these should be equal
length(exonsIwant)
dim(exon_data)[1]

#exon ratios
exonratios <- exon_data[1,] / exon_data[2,] #ratio of exon 9 over exon 10
# exonratios2 <- exon_data[2,] / exon_data[1,] - ratio of exon 10 over exon 9 if preferred 

#exon proportions
exonprops <- exon_data[1,]/(exon_data[1,] + exon_data[2,])

#getting ratios for the normal samples
normals <- which(dds_rnaseq3$Type=="Normal")
normal_ratios <- exonratios[normals]

#getting ratios for cancer samples
cancers <- which(dds_rnaseq3$Type=="Cancer")
cancer_ratios <- exonratios[cancers]

#plotting exon 9 ratios in cancer and normal samples
boxplot(
  cbind(
    unlist(as.list(normal_ratios)),
    unlist(as.list(cancer_ratios))
  ),
  col=c("blue","red"),
  main="Ratios of reads aligning to exon 9 vs exon 10"
)

#plotting exon 9 ratios of normal, glioblastoma and lower grade glioblastoma samples. 
boxplot(
  cbind(
    unlist(as.list(normal_ratios)),
    unlist(as.list(GB_ratios)),
    unlist(as.list(LGGB_ratios))
  ),
  col=c("blue","red","yellow"),
  main="Ratios of reads aligning to exon 9 vs exon 10"
)


#DONE


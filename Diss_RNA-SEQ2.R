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
my_table_features_gene <- featureCounts(files = alignment_files, annot.ext="./annotation/genome.gtf",isGTFAnnotationFile=TRUE, isPairedEnd=TRUE, GTF.featureType="gene", countMultiMappingReads = FALSE, minMQS = 30, nthreads = 20, minOverlap = 10)
# MQS = 15
# my_table_features_gene <- featureCounts(files = alignment_files, annot.ext="./annotation/genome.gtf",isGTFAnnotationFile=TRUE, GTF.featureType="gene", countMultiMappingReads = FALSE, minMQS = 15, nthreads = 10, minOverlap = 10)
# MQS = 0 (default)
#my_table_features_gene2 <- featureCounts(files = alignment_files, annot.ext="./annotation/genome.gtf",isGTFAnnotationFile=TRUE, is.PairedEnd=TRUE, GTF.featureType="gene", countMultiMappingReads = FALSE, minMQS = 0, nthreads = 10, minOverlap = 10)


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
#rld_rnaseq <- vst(dds_rnaseq_gene, blind = T)
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
boxplot(counts_rnaseq, col = adf1$Colour, las=2, names=adf1$ShortName, horizontal=TRUE, main = "GSE142083 Gene-level RNAseq counts boxplot")
dev.off()

png('fpm_rnaseq_boxplot_gene.png')
par(mar=c(7,5,3,3))
boxplot(fpm_rnaseq, col=adf1$Colour, las=2, names=adf1$ShortName, horizontal=TRUE, main = "GSE142083 Gene-level fpm_rnaseq boxplot")
dev.off()

png('assay_vsd_rnaseq_boxplot_gene.png')
par(mar=c(7,5,3,3))
boxplot(assay(vsd_rnaseq),col=adf1$Colour, las=2, names=adf1$ShortName, horizontal=TRUE, main = "GSE142083 Gene-level vsd_rnaseq boxplot")
dev.off()

png('assay_rld_rnaseq_boxplot_gene.png')
par(mar=c(7,5,3,3))
boxplot(assay(rld_rnaseq),col=adf1$Colour, las=2, names=adf1$ShortName, horizontal=TRUE, main = "GSE142083 Gene-level rld_rnaseq boxplot")
dev.off()

#but suppose we look at non-normalised data
na.rm=T
counts_un <-counts(dds_rnaseq_gene, normalized=FALSE)
png('unnormalised_counts_boxplot_gene.png')
par(mar=c(7,5,3,3))
boxplot(log2(counts_un),col = adf1$Colour, las=2, names=adf1$ShortName, horizontal = TRUE, main = "GSE142083 Gene-level unnormalised counts boxplot")
dev.off()

counts_un_rlog<-rlog(counts_un)
#OR
#counts_un_rlog<-vst(counts_un)
png(filename="mva.pairs.png",width=2000, height=2000)
mva.pairs(counts_un_rlog, main = "GSE142083 Gene-level mva_pairs plots") # (figure margins are too large)
dev.off()

#making PCA plots -  loading required libraries
library(scatterplot3d)
library(ggplot2)

#perform PCA
pca <- prcomp(t(na.omit(assay(rld_rnaseq))), scale=T)

#Plot the PCA results 
png('3D_PCA_plot_gene.png')
s3d<-scatterplot3d(pca$x[,1:3], pch=19, color=colData(dds_rnaseq_gene)$Colour, main = "GSE142083 Gene-level 3D PCA plot")
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
pheatmap(sampleDistMatrix, main = "GSE142083 Gene-level sample similarity heatmap")
dev.off()

#easier PCA plot
png('Easier_PCA_plot_gene.png')
plotPCA(rld_rnaseq, intgroup = c("Type"))
dev.off()

#Differential Gene Expression
dds_rnaseq_gene <- DESeq(dds_rnaseq)
head(dds_rnaseq_gene)

# Gene IDs for genes of interest to this study
genes <- c("ENSG00000067225", "ENSG00000143627", "ENSG00000067057", "ENSG00000134333", "ENSG00000160211", "ENSG00000130313", "ENSG00000142657", "ENSG00000153574", "ENSG00000197713", "ENSG00000163931", "ENSG00000177156", "ENSG00000138413", "ENSG00000105953", "ENSG00000137764", "ENSG00000110092")

#This is how you specify the contrast
#all results
results_genes <- results(dds_rnaseq_gene, contrast=c("Type","Cancer","Normal"), parallel = TRUE)
results_genes2 <- results(dds_rnaseq_gene, contrast=c("Type","Cancer","Normal"), parallel = TRUE)
Results <- as.data.frame(results_genes)

#Storing the results of the genes of interest in a single variable
myresults <- Results[c("ENSG00000067225", "ENSG00000143627", "ENSG00000067057", "ENSG00000134333", "ENSG00000160211", "ENSG00000130313", "ENSG00000142657", "ENSG00000153574", "ENSG00000197713", "ENSG00000163931", "ENSG00000177156", "ENSG00000138413", "ENSG00000105953", "ENSG00000137764", "ENSG00000110092"),]
myresults_ <- Results[c("ENSG00000067225", "ENSG00000067057", "ENSG00000134333", "ENSG00000160211", "ENSG00000130313", "ENSG00000142657", "ENSG00000153574", "ENSG00000197713", "ENSG00000163931", "ENSG00000177156", "ENSG00000138413", "ENSG00000105953", "ENSG00000137764", "ENSG00000110092"),]

# myresults <- results_genes[c("ENSG00000067225", "ENSG00000143627", "ENSG00000067057", "ENSG00000134333", "ENSG00000160211", "ENSG00000130313", "ENSG00000142657", "ENSG00000153574", "ENSG00000197713", "ENSG00000163931", "ENSG00000177156", "ENSG00000138413", "ENSG00000105953", "ENSG00000137764", "ENSG00000110092"),]

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
idx2 <- match(rownames(myresults), rownames(dds_rnaseq_gene))
rownames(myresults) <- ResultAnnot$nice_names[idx2]
myresults

#adding annotation to all results
idx3 <- match(rownames(results_genes), rownames(dds_rnaseq_gene))
rownames(results_genes) <- ResultAnnot$nice_names[idx3]
head(results_genes)

#adding annotation to dds_rnaseq
rownames(dds_rnaseq3) <- ResultAnnot$nice_names
rownames(dds_rnaseq) <- ResultAnnot$nice_names
# simple practice plots for users viewing. 
#Plotting the PKM counts in cancer vs normal
plotCounts(dds_rnaseq, gene = "ENSG00000067225", intgroup="Type", col = adf1$Colour, main = "GSE140343 PKM cancer vs normal gene counts") 
plotCounts(dds_rnaseq, gene = "ENSG00000143627", intgroup="Type", col = adf1$Colour, main = "GSE140343 PKLR cancer vs normal gene counts") 

#Plotting the PFKP counts in cancer vs normal
plotCounts(dds_rnaseq3, gene = "ENSG00000067057_PFKP", intgroup="Type", col = adf1$Colour, main = "GSE140343 PFKP cancer vs normal gene counts") 
#Plotting LDHA counts in cancer vs normal
plotCounts(dds_rnaseq3, gene = "ENSG00000134333_LDHA", intgroup="Type", col = adf1$Colour, main = "GSE140343 LDHA cancer vs normal gene counts") 


plotCounts(dds_rnaseq, gene = "ENSG00000143627", intgroup="Type", col = adf1$Colour, main = "GSE140343 PK cancer vs normal gene counts") 

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
dev.new()
EnhancedVolcano(Results, lab=rownames(results_genes), x = "log2FoldChange", y = 'pvalue', selectLab=c("ENSG00000067057_PFKP", "ENSG00000067225_PKM", "ENSG00000134333_LDHA"),  title='GSE140343 cancer vs normal gene-level  DESEQ results', labSize=5, pointSize=2.0, legendPosition = 'bottom', legendLabSize=15, drawConnectors=TRUE) + coord_flip()
EnhancedVolcano(Results, lab=rownames(results_genes), x = "log2FoldChange", y = 'pvalue', title='GSE142083 cancer vs normal gene-level  DESEQ results', labSize=5, pointSize=2.0, legendPosition = 'bottom', legendLabSize=15, drawConnectors=TRUE) + coord_flip()
EnhancedVolcano(Results, lab=rownames(results_genes), x = "log2FoldChange", y = 'pvalue', title='GSE142083 cancer vs normal gene-level  DESEQ results', labSize=5, pointSize=2.0, legendPosition = 'bottom', legendLabSize=15) + coord_flip()


#making a MA plot for differential gene expression
png('Main_DESeq2_plot_gene.png')
plotMA(results_genes, main="GSE142083 DESeq2 plot of unfiltered results", ylim=c(-2,2))
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
#here

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
dds_rnaseq_exon <- dds_rnaseq[keep,]
nrow(dds_rnaseq_exon)

#second filter-require >10 in >= 3 samples
keep <- rowSums(counts(dds_rnaseq_exon) >= 10) >= 3
dds_rnaseq_exon <- dds_rnaseq_exon[keep,]
nrow(dds_rnaseq_exon)

#size factor estimation for normalisation- these are used to scale by library size when data is exported.
dds_rnaseq_exon <- estimateSizeFactors(dds_rnaseq_exon)
sizeFactors(dds_rnaseq_exon)

#DEseq2 also provides methods that outputs data with various normalisations
#normalised log counts
counts_rnaseq <-log2(counts(dds_rnaseq_exon, normalized=TRUE))
#fpm is similar to spm. here the data is log transformed
fpm_rnaseq <-log(fpm(dds_rnaseq_exon)) #this is also normalised by the sizeFactors
vsd_rnaseq <- vst(dds_rnaseq_exon, blind = T)
head(vsd_rnaseq, 3)
rld_rnaseq <- vst(dds_rnaseq_exon, blind = T)
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
dds_rnaseq_exon_results <- DESeq(dds_rnaseq_exon, parallel=TRUE, minReplicatesForReplace=12)
dds_rnaseq_exon_results

# Obtaining contrast results
# cancer relative to normal
exon_results <- results(dds_rnaseq_exon_results, contrast=c("Type","Cancer","Normal"), alpha=alpha, parallel=TRUE)
exon_results <- results(dds_rnaseq_exon_results, contrast=c("Type","Cancer","Normal"), parallel=TRUE)

# normal relative to cancer
# exon_results <- results(dds_rnaseq2, contrast=c("Type","Normal","Cancer"),alpha=alpha, parallel=TRUE)

#exons of interest 
these_exons <- c("ENSE00001339202", "ENSE00003618390") # PKM1/M2

#get results of these particular exons
exon_results[match(these_exons, mytable_features_exon$annotation$GeneID),]

exon_results2 <- results(dds_rnaseq3, contrast=c("Type","Cancer","Normal"), parallel=TRUE)

exon_GB_normal_results <- results(dds_rnaseq3, contrast=c("Type","GB","Normal"), parallel=TRUE)

exon_GBnormal <- as.data.frame(exon_GB_normal_results)
exon_GBnormal[these_exons,]

#key_exon_results <- exon_results[match(these_exons, mytable_features_exon$annotation$GeneID),]

Results_exons <- as.data.frame(exon_results)
Results_exons[these_exons, ]
#plot counts of particular exons
dev.new()
par(mfrow = c(1,2))
for(g in 1:length(these_exons))
{
  plotCounts(dds_rnaseq_exon, gene = these_exons[g], intgroup="Type", normalized=FALSE, transform=FALSE, ylim=c(0,11000), las=2, col = adf1$Colour ,pch=19)
}

dev.new()
par(mfrow = c(1,2))
for(g in 1:length(these_exons))
{
  plotCounts(dds_rnaseq_exon, gene = these_exons[g], intgroup="Type", normalized=FALSE, transform=FALSE, las=2, col = adf1$Colour ,pch=19)
}

dev.new()
par(mfrow = c(1,2))
for(g in 1:length(these_exons))
{
  plotCounts(dds_rnaseq3, gene = these_exons[g], intgroup="Type", normalized=FALSE, transform=FALSE, las=2, col = adf1$Colour ,pch=19)
}

#getting the number of read counts for each exon for comparison and ratio determination
exon_featurecounts <- as.data.frame(table_rnaseq_exon)
exon_featurecounts["ENSE00003618390",] #exon 10 counts
exon_featurecounts["ENSE00001339202",]# exon 9 counts

#Exons of interest 
exonsIwant <- c("ENSE00001339202", "ENSE00003618390")

#getting index for these exons
myexonIndex <- match(exonsIwant, rownames(exon_featurecounts))
exon_data <- exon_featurecounts[myexonIndex,]

#turn the exon data into a table for supevisors
df <- as.data.frame(t(as.data.frame(exon_data)))
write.csv(df, file="PKM1_M2_counts.csv", sep='\t')

#adding columns for sample type
df_new <- cbind(dds_rnaseq2$Type, df)
names(df_new)
names(df_new)[1] <- "SampleType"

#sorting these new dfs by sample type
df_new <- df_new[order(df_new$SampleType),]
#saving new dfs has a file
write.csv(df_new, file="GSE147352_new_PKM1_M2.csv", sep = '\t')

#adding exon length information to exon_counts dft
exon_info <- read.table("Homo_sapiens.104.exon.excel.xls", header=FALSE, sep='\t')
exon_info <- as.data.frame(exon_info)
exon_lengths <- subset(exon_info, select = c(10,4))
names(exon_lengths)[2] <- "exon_length"
#remove column name for the exon Ids and remove row numbers
names(exon_lengths)[1] <- "ExonID"
typeof(exon_lengths)
head(exon_lengths)
#exon_lengths_df <- as.data.frame(t(as.data.frame(exon_lengths)))


# which ones do we want: this should give a vector of positions in the exon_lengths
exonswanted <- match(names(df_new),exon_lengths$ExonID)
#now just get the values into its own df
these_exon_lengths <- exon_lengths[exonswanted,]
# now r bind the size columns, eg put the sizes in the first row of current datafile
df_new_size_row <- rbind(t(these_exon_lengths)[2,], df_new)

#now divide each of the rows in the table by the first row
df_normalised1 <- as.matrix(df_new_size_row)
x1 <- df_normalised1[,-1]
for(i in 1:ncol(x1)){
  x1[, i] <- (t(t(as.numeric(x1[,i]))/as.numeric(x1[1,i])))
}
normalised_df_1 <- as.data.frame(unlist(x1))
# re add the sample type column
normalised_df_1$SampleType <- df_new_size_row$SampleType
normalised_df1 <- normalised_df_1[c(3,1:2)]


#do the same with all PKLR exons...
exonsIwant2 <- c("ENSE00001624144", "ENSE00001791302", "ENSE00003512685", 
                 "ENSE00003658804", "ENSE00000961246", "ENSE00000961247", 
                 "ENSE00000961248", "ENSE00000961249", "ENSE00000961250", 
                 "ENSE00000961251", "ENSE00001373175", "ENSE00002299405", 
                 "ENSE00001835019", "ENSE00002163140", "ENSE00002149851", 
                 "ENSE00002200075")

#PKLR_exons <- PKLR_exons <- read.table("PKLR_exons", header=FALSE)
#cat(sprintf('"%s",', PKLR_exons$V1))
#exonsIwant2 <- c("ENSE00002299405", "ENSE00001791302", "ENSE00003512685", 
                 "ENSE00003658804", "ENSE00000961246", "ENSE00000961247", 
                 "ENSE00000961248", "ENSE00000961249", "ENSE00000961250", 
                 "ENSE00000961251", "ENSE00001835019", "ENSE00001624144", 
                 "ENSE00001791302", "ENSE00003512685", "ENSE00003658804",
                 "ENSE00000961246", "ENSE00000961247", "ENSE00000961248", 
                 "ENSE00000961249", "ENSE00000961250", "ENSE00000961251", 
                 "ENSE00001373175", "ENSE00002163140", "ENSE00002149851", 
                 "ENSE00003512685", "ENSE00003658804", "ENSE00002200075")

myexonIndex2 <- match(exonsIwant2, rownames(exon_featurecounts))
exon_data2 <- exon_featurecounts[myexonIndex2,]
df2 <- as.data.frame(t(as.data.frame(exon_data2)))
write.csv(df2, file="GSE147352_PKLR_all_exon_counts.csv", sep='\t')

#adding columns for sample type
df_new2 <- cbind(dds_rnaseq2$Type, df2)
names(df_new2)
names(df_new2)[1] <- "SampleType"

#sorting these new dfs by sample type
df_new2 <- df_new2[order(df_new2$SampleType),]
#saving new dfs has a file
write.csv(df_new2, file="GSE147352_new_PKLR_exons.csv", sep = '\t')

#adding exon lengths to the dfs
exonswanted2 <- match(names(df_new2),exon_lengths$ExonID)

these_exon_lengths2 <- exon_lengths[exonswanted2,]

# now r bind the size columns, eg put the sizes in the first row of current datafile
df_new_size_row2 <- rbind(t(these_exon_lengths2)[2,], df_new2)

#getting the normalised results
df_normalised2 <- as.matrix(df_new_size_row2)
x2 <- df_normalised2[,-1]
for(i in 1:ncol(x2)){
  x2[, i] <- (t(t(as.numeric(x2[,i]))/as.numeric(x2[1,i])))
}
normalised_df_2 <- as.data.frame(unlist(x2))
# re add the sample type column
normalised_df_2$SampleType <- df_new_size_row2$SampleType
dim(normalised_df_2)
normalised_df2 <- normalised_df_2[c(17,1:16)]

#and then PKM from my particular transcript
exonsIwant3 <- c("ENSE00001546931", "ENSE00003522795", "ENSE00003790251", 
                 "ENSE00003566465", "ENSE00003535056", "ENSE00003590685", 
                 "ENSE00003475296", "ENSE00003647308", "ENSE00003618390", 
                 "ENSE00003632285", "ENSE00001560298")


myexonIndex3 <- match(exonsIwant3, rownames(exon_featurecounts))
exon_data3 <- exon_featurecounts[myexonIndex3,]
df3 <- as.data.frame(t(as.data.frame(exon_data3)))
write.csv(df3, file="GSE147352_PKM_all_transcript_exon_counts.csv", sep='\t')

#adding columns for sample type
df_new3 <- cbind(dds_rnaseq2$Type, df3)
names(df_new3)
names(df_new3)[1] <- "SampleType"

#sorting these new dfs by sample type
df_new3 <- df_new3[order(df_new3$SampleType),]
#saving new dfs has a file
write.csv(df_new3, file="GSE147352_new_PKM_transcript.csv", sep = '\t')

# adding the exon lengths to counts df
exonswanted3 <- match(names(df_new3),exon_lengths$ExonID)
these_exon_lengths3 <- exon_lengths[exonswanted3,]

# now r bind the size columns, eg put the sizes in the first row of current datafile
df_new_size_row3 <- rbind(t(these_exon_lengths3)[2,], df_new3)

#Normalise the data!
df_normalised3 <- as.matrix(df_new_size_row3)
x3 <- df_normalised3[,-1]
for(i in 1:ncol(x3)){
  x3[, i] <- (t(t(as.numeric(x3[,i]))/as.numeric(x3[1,i])))
}
normalised_df_3 <- as.data.frame(unlist(x3))
# re add the sample type column
normalised_df_3$SampleType <- df_new_size_row3$SampleType
dim(normalised_df_3)
normalised_df3 <- normalised_df_3[c(12,1:11)]

# and then all PKM exons from all PKM transcripts, 
exonsIwant4 <- c("ENSE00002615002", "ENSE00003535770", "ENSE00003548841",
                 "ENSE00003715228", "ENSE00003659292", "ENSE00003559627", 
                 "ENSE00003566465", "ENSE00003535056", "ENSE00003590685", 
                 "ENSE00003475296", "ENSE00003647308", "ENSE00003618390", 
                 "ENSE00003632285", "ENSE00003563330", "ENSE00001339295", 
                 "ENSE00003522795", "ENSE00003790251", "ENSE00003566465", 
                 "ENSE00003535056", "ENSE00003590685", "ENSE00003475296", 
                 "ENSE00003647308", "ENSE00001339202", "ENSE00003632285", 
                 "ENSE00003563330", "ENSE00003715228", "ENSE00003522795", 
                 "ENSE00003747482", "ENSE00003752520", "ENSE00003590685", 
                 "ENSE00003475296", "ENSE00003647308", "ENSE00003618390", 
                 "ENSE00003632285", "ENSE00003563330", "ENSE00002618127", 
                 "ENSE00003522795", "ENSE00003790251", "ENSE00003566465", 
                 "ENSE00003535056", "ENSE00003590685", "ENSE00003475296", 
                 "ENSE00003647308", "ENSE00001339202", "ENSE00003632285", 
                 "ENSE00001560298", "ENSE00003522795", 
                 "ENSE00003790251", "ENSE00003566465", "ENSE00003535056", 
                 "ENSE00003590685", "ENSE00003475296", "ENSE00003647308", 
                 "ENSE00003618390", "ENSE00003632285", "ENSE00001560298", 
                 "ENSE00002586217", "ENSE00002582218", "ENSE00002580368", 
                 "ENSE00003475296", "ENSE00003647308", "ENSE00001339202", 
                 "ENSE00003632285", "ENSE00001560298", "ENSE00002601471", 
                 "ENSE00002590935", "ENSE00001546931", "ENSE00003467662",
                 "ENSE00003522795", "ENSE00003534760", "ENSE00003497716", 
                 "ENSE00003618174", "ENSE00003509895", "ENSE00003646231", 
                 "ENSE00003647862", "ENSE00003535770", "ENSE00003651040", 
                 "ENSE00002600860", "ENSE00003467662", "ENSE00003522795", 
                 "ENSE00003790251", "ENSE00003566465", "ENSE00003535056", 
                 "ENSE00003590685", "ENSE00003475296", "ENSE00003647308", 
                 "ENSE00001339202", "ENSE00003632285", "ENSE00002601115", 
                 "ENSE00002591894", "ENSE00003522795", "ENSE00003790251", 
                 "ENSE00003566465", "ENSE00003535056", "ENSE00003590685", 
                 "ENSE00003475296", "ENSE00003647308", "ENSE00001339202", 
                 "ENSE00003632285", "ENSE00002614041", "ENSE00002592781", 
                 "ENSE00003522795", "ENSE00003790251", "ENSE00003566465", 
                 "ENSE00003535056", "ENSE00003590685", "ENSE00003475296", 
                 "ENSE00003647308", "ENSE00001339202", "ENSE00002599842", 
                 "ENSE00002616730", "ENSE00002604788", "ENSE00001360078", 
                 "ENSE00003522795", "ENSE00003790251", "ENSE00003566465", 
                 "ENSE00003535056", "ENSE00002595435", "ENSE00003558214", 
                 "ENSE00003646231", "ENSE00002606051", "ENSE00002628200", 
                 "ENSE00002577098", "ENSE00002598286", "ENSE00002619480", 
                 "ENSE00002596208", "ENSE00003509895", "ENSE00002607472", 
                 "ENSE00002625069", "ENSE00002623632", "ENSE00002589475", 
                 "ENSE00003567910", "ENSE00002587987", "ENSE00002629928", 
                 "ENSE00003522795", "ENSE00003790251", "ENSE00003566465", 
                 "ENSE00002617601", "ENSE00002603073", "ENSE00003522795", 
                 "ENSE00003790251", "ENSE00003566465", "ENSE00002625427", 
                 "ENSE00002576229", "ENSE00003522795", "ENSE00003790251", 
                 "ENSE00003566465", "ENSE00002595493", "ENSE00001546931", 
                 "ENSE00003489416", "ENSE00003597522", "ENSE00002618575", 
                 "ENSE00002601215", "ENSE00003522795", "ENSE00003790251", 
                 "ENSE00002586023", "ENSE00002614691", "ENSE00002576672")


exonsIwant4_unique <- unique(exonsIwant4)
myexonIndex4 <- match(exonsIwant4_unique, rownames(exon_featurecounts))
exon_data4 <- exon_featurecounts[myexonIndex4,]
df4 <- as.data.frame(t(as.data.frame(exon_data4)))
write.csv(df4, file="GSE147352_PKM_all_exon_counts.csv", sep='\t')

#adding columns for sample type
df_new4 <- cbind(dds_rnaseq2$Type, df4)
names(df_new4)
names(df_new4)[1] <- "SampleType"

#sorting these new dfs by sample type
df_new4 <- df_new4[order(df_new4$SampleType),]
#saving new dfs has a file
write.csv(df_new4, file="GSE147352_new_PKM_all.csv", sep = '\t')

#adding exon lengths row to the exon_read counts dfs

exonswanted4 <- match(names(df_new4),exon_lengths$ExonID)
these_exon_lengths4 <- exon_lengths[exonswanted4,]

# now r bind the size columns, eg put the sizes in the first row of current datafile
df_new_size_row4 <- rbind(t(these_exon_lengths4)[2,], df_new4)

#Normalise the data!
df_normalised4 <- as.matrix(df_new_size_row4)
x4 <- df_normalised4[,-1]
for(i in 1:ncol(x4)){
  x4[, i] <- (t(t(as.numeric(x4[,i]))/as.numeric(x4[1,i])))
}
normalised_df_4 <- as.data.frame(unlist(x4))
# re add the sample type column
normalised_df_4$SampleType <- df_new_size_row4$SampleType
dim(normalised_df_4)
normalised_df4 <- normalised_df_4[c(73,1:72)]

#make sure its the right datae
head(exon_data)

#checking - these should be equal
length(exonsIwant)
dim(exon_data)[1]

#exon ratios
exonratios <- exon_data[1,] / exon_data[2,] #ratio of exon 9 over exon 10
#exonratios2 <- exon_data[2,] / exon_data[1,] - ratio of exon 10 over exon 9 if preferred 

#exon proportions
exonprops <- exon_data[1,]/(exon_data[1,] + exon_data[2,])

#getting ratios for the normal samples
normals <- which(dds_rnaseq_exon_results$Type=="Normal")
normal_ratios <- exonratios[normals]

#getting ratios for cancer samples
cancers <- which(dds_rnaseq_exon_results$Type=="Cancer")
cancer_ratios <- exonratios[cancers]

#plotting exon 9 ratios in cancer and normal samples
boxplot(
  cbind(
    unlist(as.list(normal_ratios)),
    unlist(as.list(cancer_ratios))
  ),
  col=c("blue","red"),
  main="GSE142083 ENSE00001339202 (exon 9)  to ENSE00003618390 (exon 10) ratio",
  ylab="Ratio of PKM1 (exon 9) to PKM2 (exon 10)",
  cex.axis=1.25,
  cex.labels=1.5
)

#plotting exon 9 ratios of normal, glioblastoma and lower grade glioblastoma samples. 
boxplot(
  cbind(
    unlist(as.list(normal_ratios)),
    unlist(as.list(GB_ratios)),
    unlist(as.list(LGGB_ratios))
  ),
  col=c("blue","red","yellow"),
  main="GSE147352 ENSE00001339202 (exon 9) to ENSE00003618390 (exon 10) ratio",
  ylab="Ratio of PKM1 (exon 9) to PKM2 (exon 10)",
  cex.axis=1.25,
  cex.labels=1.5
)


#DONE


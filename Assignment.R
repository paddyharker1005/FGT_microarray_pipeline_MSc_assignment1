##the code in this script is adapted from the BasicAnalysisScript provided in class

setwd("~/FGT/Assignment")
#load the required libraries for the workflow
library(affy)
library(limma)
library(annotate)
library(scatterplot3d)
library(mouse4302.db)# load chip-specific annotation
library(EnhancedVolcano)

#read in a targets file
targets <- read.table("targets.txt",header=T,as.is=T)
#load all CEL files in the R working directory
mydata <- ReadAffy(filenames = targets$Filename)
#assign the sample names as the mousenames from the targets file - easier to understand
sampleNames(mydata) <-targets$MouseName

#Quality control plot - save to file
png("plots/QC_histogram.png")
hist(mydata, main = "Diagnostic Histogram of Data set", cex.main = 2, cex.lab=1.3)
dev.off()

# Boxplot with different colour per sample group
png("plots/QC_boxplot.png")
boxplot(mydata, col=targets$Colour, las=1, cex.axis=0.8, cex.lab=1.3, cex.main = 2, xlab = "Sample", ylab = "Log Intensity", main = "Diagnostic Boxplot of Data set")
dev.off()

# Normalise data using RMA
eset <- rma(mydata)

# Obtain a matrix of the expression values
values <- exprs(eset)

# Boxplot with normalized data
png("plots/QC_boxplot_normalised.png")
boxplot(values, col=targets$Colour, las=1, cex.axis=0.8, cex.lab=1.3, cex.main = 2, xlab = "Sample", ylab = "Log Intensity", main = "Boxplot of RMA normalised Data")
dev.off()

# MA plot of the samples 1 and 4
png("plots/MA_plot_normalised.png", width = 800, height = 600)
mva.pairs(values[, c(1,2,3,4,5,6)], cex=1.6, main = "MvA (MvsA) Diagnostic Plot Normalised Data")
dev.off()

# The same plot for the non-normalised raw data
png("plots/MA_plot_nonnormalised.png", width = 800, height = 600)
mva.pairs(pm(mydata)[, c(1,2,3,4,5,6)],cex=1.6, main = "MvA (MvsA) Diagnostic Plot Non-normalised Data", cex.main=2)
dev.off()



# Performs hierarchical clustering with average linkage based on
# Pearson’s Correlation Coefficient and save to file

hc<-hclust(as.dist(1-cor(values, method="pearson")), method="average")
png("plots/HC_dendogram.png", width = 800, height = 400)
plot(hc, main = "Hierarchical clustering of normalised data based on \n Pearson's Correlation Coefficient", xlab = NULL, cex.main =2, cex.lab = 1.5, font =2)
dev.off()

## Perform PCA

pca <- prcomp(t(values), scale=T)

# Plot the PCA results and save to file
png("plots/PCA_plot.png")
s3d<-scatterplot3d(pca$x[,1:3], pch=19, color=c("yellow", "yellow 2", "yellow 3", "red", "red 3", "red 4"), main = "Principal Components Analysis (PCA) \n of normalised data", cex.lab=1.5, cex.main =1.5)
s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
text(s3d.coords$x, s3d.coords$y, labels = colnames(values),pos = 3,offset = 0.5, cex=0.8, font=2)
dev.off()

## Perform fold filtering
#obtaining a matrix of expression values
exprsvals <- exprs(eset)
#RMA outputs log2 data while MAS5 outputs linear data
#To convert from log…
exprsvals10 <-2^exprsvals
#check conversion
exprsvals[1:10,]
#converted
exprsvals10[1:10,]

#More fold filtering
#check order of sample names and probe sets
mysamples <- sampleNames(eset)
mysamples
probesets <- probeNames(mydata)
probesets[1:50]

#Build final fold table
#Calculate the means - removed DNMTs-Tg_plus.1, but at reduced power
#Note mean of the log is not the same as the log of the mean!!
WT.mean <- apply(exprsvals10[,c("WT.1", "WT.2", "WT.3")],1,mean)
DNMTs_TG_min.mean <- apply(exprsvals10[,c("DNMTs-Tg_min.1", "DNMTs-Tg_min.2")],1,mean)

#calculate fold changes

WT_DNMT_Tg_min <-  DNMTs_TG_min.mean/WT.mean

#build a summary table to hold all the data

all.data= cbind(WT.mean, DNMTs_TG_min.mean, WT_DNMT_Tg_min)
colnames(all.data)

write.table(all.data,file="group_means.txt", quote=F, sep="\t",col.names=NA)


## Beginning statistical analysis


#Check original sample order
sampleNames(eset)
##Building annotation for differential gene identification

eset@annotation

#packages in the annotation package
ls("package:mouse4302.db")

#build an annotation table
ID <- featureNames(eset)
Symbol <- getSYMBOL(ID, "mouse4302.db")
Name <- as.character(lookUp(ID, "mouse4302.db", "GENENAME"))
tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name, stringsAsFactors=F)
tmp[tmp=="NA"] <- NA #fix padding with NA characters 
#assign as feature data of the current Eset
fData(eset) <- tmp


##Shiny expression data - extract 50 most variable probe sets between samples
variances <- apply(exprsvals10, 1, var)

expression <- exprsvals10[order(variances, decreasing=TRUE)[1:50],]

gene_names <- tmp$Symbol[match(rownames(expression), tmp$ID)]

rownames(expression)[! is.na(gene_names)] <- paste(gene_names[! is.na(gene_names)], rownames(expression)[! is.na(gene_names)], sep = ' / ')

save(expression,file="expression.Rdata")

pheatmap(expression, fontsize_row = 4, fontsize_col=10, show_rownames=TRUE, scale="row")

## Statistical analysis using Limma

design <- model.matrix(~-1+factor(c(1,1,1,2,2,2)))
colnames(design) <- c("DNMT_Tg", "WT")
#output design matrix
design

#This instructs Limma which comparisons to make
contrastmatrix <- makeContrasts(DNMT_Tg-WT, levels=design)
contrastmatrix

#issue these commands to fit the model
#and make the contrasts
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrastmatrix)

#this last part essentially moderates the t-statistic using 
#the borrowed variance approach described in class
fit2 <- eBayes(fit2)
#get the results
topTable(fit2,coef=1,adjust="fdr")
myresults <-topTable(fit2,coef=1, adjust="fdr", number=nrow(eset))
write.table(myresults,"myresults.txt", sep="\t")

#table for top 10 based on Pval
myresultsPval10 <- head(myresults[order(abs(myresults$adj.P.Val), decreasing = FALSE),], 10)
write.table(myresultsPval10,"myresultsPval10.txt", sep="\t")

#table for top 10 based on FC
myresultsFC10 <- head(myresults[order(abs(myresults$logFC), decreasing = TRUE),], 10)
write.table(myresultsFC10,"myresultsFC10.txt", sep="\t")

##volcano plot - adapted code from https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
#assign colours to upregulated and downregulated genes with custom colour scheme using key-value pairs
keyvals <- ifelse(
  (myresults$logFC < 0 & myresults$adj.P.Val < 0.05), 'blue',
  ifelse((myresults$logFC > 0 & myresults$adj.P.Val < 0.05), 'red', 'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red'] <- 'Up'
names(keyvals)[keyvals == 'grey'] <- 'No change'
names(keyvals)[keyvals == 'blue'] <- 'Down'

pdf("plots/Volcano.pdf")
EnhancedVolcano(myresults, lab=myresults$Symbol, colCustom = keyvals, x = 'logFC', y = 'adj.P.Val' , pCutoff=0.05, FCcutoff=0, pointSize = 1, xlim = c(-9, 9), ylim = c(0, 8),  pLabellingCutoff = 1e-03, labSize = 2)
dev.off()

##shiny differential expression data
save(myresults,file="differential_expression.Rdata")


## Carry out Functional Enrichment analysis

Mm.H <- readRDS("/shared_files/MSigDB/Mm.h.all.v7.1.entrez.rds") 

#Select from the annotation a number of keys with the primary key being PROBEID
res <- select(mouse4302.db, keys = rownames(eset), columns = c("ENTREZID", "ENSEMBL","SYMBOL"), keytype="PROBEID")
#View the top of the table
head(res)
#find the index of each row of the expression set in the #annotation object res
idx <- match(rownames(eset), res$PROBEID)
#Use the index to set the phenotypic data in the ExpressionSet
fData(eset) <- res[idx, ]
head(fData(eset), 10)
#Find all rows that don’t have an EntrezID and remove then
eset_t<-eset[is.na(fData(eset)$ENTREZID)==0,]


## Functional Enrichment Analysis

#convert to indexes
H.indices <- ids2indices(Mm.H,fData(eset_t)$ENTREZID)

#run mroast
results <-mroast(eset_t,index=H.indices,design=design,contrast=contrastmatrix[,1],adjust.method = "BH")
#run camera
#this is the tool used in the analysis
results2 <-camera(eset_t,index=H.indices,design=design,contrast=contrastmatrix[,1])
#run romer
results3 <-romer(eset_t,index=H.indices,design=design,contrast=contrastmatrix[,1],adjust.method = "BH")


#we can also use exactly the same model as our differential gene analysis for
#the enrichment analysis- in this case we can extract it from
#the fit
sv <- squeezeVar(fit$sigma^2,df=fit$df.residual)
results_prior <-mroast(eset_t,index=H.indices,design=design,contrast=contrastmatrix[,1], var.prior = sv$var.prior, df.prior = sv$df.prior, adjust.method = "BH")

#write out the results of the functional enrichment analysis to individual files
write.table(results,"enrichment_mroast.txt",sep="\t")
write.table(results2,"enrichment_camera.txt",sep="\t")
write.table(results3,"enrichment_romer.txt",sep="\t")
write.table(results_prior,"enrichment_mroast_prior.txt",sep="\t")


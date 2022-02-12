############# RNA-seq analysis #############
### VERY important: https://hbctraining.github.io/main/
  
############# In Terminal #############

### Unzip raw data
#   $ gunzip *.gz


### Trim adapter or low quality sequences in local laptop using Terminal followed by Quality check by Fastqc
# https://github.com/FelixKrueger/TrimGalore
#   $ trim_galore -o ./1_trimming *.fastq --paired --fastqc --illumina


### Make an index genome by bowtie2
# https://github.com/BenLangmead/bowtie2
#   $ cp GCF_000007805.1_ASM780v1_genomic.fna GCF_000007805.1_ASM780v1_genomic.fa
#   $ bowtie2-build GCF_000007805.1_ASM780v1_genomic.fa ./Index/PsyDC3000


### Alignment of the sequence reads to the index genome by bowtie2
#   $ for f1 in *R1_val_1.fq ; do for f2 in ${f1%%_R1_val_1.fq}"_R2_val_2.fq" ; do bowtie2 -x ./Index/PsyDC3000 -1 $f1 -2 $f2 -S ./2_sam/$f1.sam ; done; done


### File format change by samtools
# https://github.com/samtools/samtools
# 1)  $ for file in ./*.sam; do samtools view -bS $file -o ./3_bam/$file.bam ; done
# 2)  $ for file in ./*.bam; do samtools sort $file -o ./4_bam_sort/$file.sort.bam ; done
# 3)  $ for file in ./*.bam; do samtools index $file ; done


### Read counts by htseq
# https://htseq.readthedocs.io/en/release_0.11.1/index.html
#   $ for file in ./*.bam ; do htseq-count -t gene -i Name --stranded=no -f bam -r pos ./$file ./GCF_000007805.1_ASM780v1_genomic.gff > ./5_count/$file.count ; done





############# In R or Rstudio #############

### Identification of differently expressed genes (DEG) by DESeq2
# http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
rm(list = ls()) 
 
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(dendextend)


### Create a DESeqDataSet object
# Read the count file containing all bio reps
mycount <- as.matrix(read.csv("Data_2_2021.csv", row.names="gene_id"))
dim(mycount)
class(mycount)
head(mycount, 2)

# Read the meta file 
coldata <- read.csv("Meta_2_2021.csv", row.names = 1)
dim(coldata)
class(coldata)
coldata <- coldata[,c("Strain", "Chemical")]

# Are they identical as they follow different orders?
all(rownames(coldata) == colnames(mycount))
# If needed, reorder the columns of the count matrix according to the order of Label in columnData
#mycount <- mycount[,rownames(coldata)]
#all(rownames(coldata) == colnames(mycount))


### Make DESeq object on the basis of the counts
dds.mydata <- DESeqDataSetFromMatrix(
  mycount,                      # Count matrix
  coldata,                      # metadata
  ~ Chemical)  # design formula

dds.mydata$Chemical <- relevel(dds.mydata$Chemical, ref = "none")
dds.mydata
slotNames(dds.mydata)
dds.mydata@design # Check your design for the comparison.
dds.mydata@colData

### The steps for a basic analysis are: 1) estimate size factors, 2) estimate dispersion parameters, and carry out 3) DE analysis.
### 1) Estimate Size Factors to inspect the raw count data which is generated from various lanes and bio reps. 
dds.mydata <- estimateSizeFactors(dds.mydata)
dds.mydata
sizeFactors(dds.mydata)
normalized_counts <- counts(dds.mydata, normalized=TRUE)
write.csv(normalized_counts, file="mydata1_norm.csv")

sizeFactors(dds.mydata) %>%
  as.data.frame %>%
  rownames_to_column -> mydf

colnames(mydf)[2] <- "sizefac"
mydf
ggplot(mydf, aes(rowname, sizefac)) + geom_point() + theme(axis.text.x = element_text(face="bold", color="blue", angle=45))

### 2) Dispersion Parameters
#Dispersion is a measure of spread or variability in the data.
#DESeq2 uses a specific measure of dispersion (α) related to the mean (μ) and variance of the data: Var = μ + α*μ^2. 
#DESeq2 dispersion estimates are inversely related to the mean and directly related to variance. Based on this relationship, the dispersion is higher for small mean counts and lower for large mean counts. The dispersion estimates for genes with the same mean will differ only based on their variance. Therefore, the dispersion estimates reflect the variance in gene expression for a given mean value.

dds.mydata <- estimateDispersions(dds.mydata)
dds.mydata
dds.mydata@dispersionFunction
alphas <- dispersions(dds.mydata) # Verify that the number of dispersion factors equals the number of genes
length(alphas) # The number of disperion factors
mcols(dds.mydata)
#baseMean:	  mean of normalized counts for all samples
#baseVar: 	  variance of normalized counts for all samples
#allZero: 	  all counts for a gene are zero
#dispGeneEst:	gene-wise estimates of dispersion
#dispFit: 	  fitted values of dispersion
#dispersion:	final estimate of dispersion
#dispIter:  	number of iterations
#dispOut:   	dispersion flagged as outlier
#dispMAP:   	maximum a posteriori estimate
boxplot(log(dispersions(dds.mydata))) # Summarize the dispersion factors using a box plot
plotDispEsts(dds.mydata) # Plot the dispersion estimates. The shrinkage method is applied to reduce false positives in the differential expression analysis.

### 3-1) Differential Expression Analysis
# We can now conduct a differential expression analysis using the DESeq() function. 
# Keep in mind that to get to this step, we first estimated the size factors and then the dispersion parameters.
ddsDE <- DESeq(dds.mydata)
ddsDE # Look at object
colSums(counts(ddsDE)) # Total number of raw counts per sample
raw.count <- colSums(counts(ddsDE)) %>%
  data.frame() %>%
  rownames_to_column(var="sample") %>%
  as_tibble()
ggplot(raw.count, aes(x=sample, y=.)) + geom_point() + theme(axis.text.x = element_text(face="bold", color="blue", angle=45))

colSums(counts(ddsDE, normalized = T)) # Total number of normalized counts per sample
norm.count <- colSums(counts(ddsDE, normalized = T)) %>%
  data.frame() %>%
  rownames_to_column(var="sample") %>%
  as_tibble()
ggplot(norm.count, aes(x=sample, y=.)) + geom_point() + theme(axis.text.x = element_text(face="bold", color="blue", angle=45))

results(ddsDE) # We can get the results for the differential expression analysis.
results(ddsDE, tidy = TRUE)
summary(results(ddsDE))

# Clustering: Regularized log transformation. The regularized log transform can be obtained using the rlog() function. Note that an important argument for this function is blind (TRUE by default). The default "blinds" the normalization to the design. This is very important so as to not bias the analyses (e.g. class discovery)
rld <- rlog(ddsDE, blind = TRUE)
rld.table <- assay(rld)
write.csv(rld.table, "mydata1_norm_log2.csv") #edit the first column name as "Locus" in log2_norm_count.csv from your local folder
sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(coldata$rownames)
#rownames(sampleDistMatrix) <- paste(rld$Strain, rld$Chemical, sep="-")
#colnames(sampleDistMatrix) <- paste(rld$Strain, rld$Chemical, sep="-")
#colnames(sampleDistMatrix) <- NULL # Use this, if you want to remove the column name.
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.print(pdf, 'Fig1_heatmap.pdf', width = 10, height =12)

# Dendrogram of samples: showing strain & media of each sample. Hierarchical clustering using rlog transformation.
options(repr.plot.width = 9, repr.plot.height = 5)
plot(hclust(sampleDists, method = "complete"))
dev.print(pdf, 'Fig2_dendrogram.pdf', width = 10, height =12)

# Principal Components Analysis (i.e., Color by condition): used to emphasize variation and bring out strong patterns in a dataset (dimensionality reduction).
# If two samples have similar levels of expression for the genes that contribute significantly to the variation represented by PC1, they will be plotted close together on the PC1 axis. 
plotPCA(rld, intgroup = c("Chemical"))
dev.print(pdf, 'Fig3_PCA.pdf', width = 10, height =12)

# Here, you can compare two group of samples specified by the contrast. 
# If not, the default contrast would be the last term in your additive model design(dds).
coldata # none, Hopea
myres_EryVSnone <- results(ddsDE, contrast = c("Chemical", "Hopea", "none"))
myres_EryVSnone
myres_table_myres_EryVSnone <- myres_EryVSnone %>% # Make a result table
  data.frame() %>%
  rownames_to_column(var="Locus") %>%
  as_tibble()
myres_table_myres_EryVSnone

# Gen Key
anno.data <- read.csv("gene_key.csv")
myres_table_myres_EryVSnone$Old_locus <- anno.data$old_locus_tag[match(myres_table_myres_EryVSnone$Locus, anno.data$locus_tag)]
myres_table_myres_EryVSnone$Annotation <- anno.data$annotation[match(myres_table_myres_EryVSnone$Locus, anno.data$locus_tag)]
# log2FC > 0 suggests that "Hopea" is associated with increased expression. 
# log2FC < 0 suggests that "Hopea" is associated with lower expression. 
# Therefore, "none" is the reference condition. You can change this contrast by flipping like this: c("Chemical", "none", "Hopea")

myres_table_myres_EryVSnone <- myres_table_myres_EryVSnone[order(myres_table_myres_EryVSnone$padj),]   # Sorted in ascending order by adjusted p-value
dim(myres_table_myres_EryVSnone)
head(myres_table_myres_EryVSnone, 5)
write.csv(myres_table_myres_EryVSnone, file="Res_0_all_total.csv")


myres_table_myres_EryVSnone_sig <- myres_table_myres_EryVSnone %>%    # Threshold applied (significant genes)
  filter(abs(log2FoldChange) > 1 , padj < 0.05)
dim(myres_table_myres_EryVSnone_sig)
write.csv(myres_table_myres_EryVSnone_sig, file="Res_1_Hopea_vs_none_log2FC1_padj05_all.csv")

# up-regulated genes
myres_table_myres_EryVSnone_sig_up <- myres_table_myres_EryVSnone_sig %>% filter(log2FoldChange > 0)
dim(myres_table_myres_EryVSnone_sig_up)
write.csv(myres_table_myres_EryVSnone_sig_up, file="Res_1_Hopea_vs_none_log2FC1_padj05_up.csv")

# down-regulated genes
myres_table_myres_EryVSnone_sig_down <- myres_table_myres_EryVSnone_sig %>% filter(log2FoldChange < 0)
dim(myres_table_myres_EryVSnone_sig_down)
write.csv(myres_table_myres_EryVSnone_sig_down, file="Res_1_Hopea_vs_none_log2FC1_padj05_down.csv")



# Volcano  plot (red if both padj < 0.05 and log2FC > 1)
#https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

library(ggplot2)
x <- myres_table_myres_EryVSnone
x$expression = ifelse(x$padj < 0.05 & abs(x$log2FoldChange) > 1, 
                      ifelse(x$log2FoldChange > 1 ,'Up','Down'),'Stable')

p <- ggplot(data = x, 
            aes(x = log2FoldChange, 
                y = -log10(x$padj), 
                colour=expression,
                label = x$Locus)) +
  geom_point(alpha=0.4, size=2.5) +
  scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-10, 10)) +
  geom_vline(xintercept=c(-1,1),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = 1.301,lty=2,col="black",lwd=0.5) +
  labs(x="log2(fold change)",
       y="-log10 (adjusted P)",
       title="Differential expression")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        text = element_text(size = 20))
p
dev.print(pdf, 'Res_2.pdf', width = 7, height =9)

##------------------------
# Functional enrichment
library(tidyverse)
library(reutils)
library(plyr)
library(ggplot2)

rm(list = ls()) 
`%ni%` = Negate(`%in%`)
#this function converts genomes annotation files (.gff) from NCBI to a readable table
gff.parse <- function(x) {
  ids <- list()
  x %>% filter(type == "gene" | type == "CDS" | type == "tRNA" | type == "rRNA" ) -> x
  
  #if ("old_locus_tag" %in% colnames(x)) {
  if ("locus_tag" %in% colnames(x)) {
    for (i in 1:(nrow(x)-1)) {
      if (all(x[i,2:4] == x[i+1,2:4])) {
        vec <- data.frame("chr"=x[["seqnames"]][i],
                          "acc"= x[["protein_id"]][(i+1)], 
                          "locus_tag" = x[["locus_tag"]][i], 
                          "old_locus_tag" = x[["old_locus_tag"]][i], 
                          "length_nt"= x[["width"]][i], 
                          "start"= x[["start"]][i], 
                          "end"= x[["end"]][i],
                          "strand" = x[["strand"]][i],
                          "type"= x[["gene_biotype"]][i],
                          "annotation"= x[["product"]][(i+1)],
                          stringsAsFactors = F)
        ids[[i]] <- vec
      }}} 
  else if ("gene" %in% colnames(x)) {
    for (i in 1:(nrow(x)-1)) {
      if (all(x[i,2:4] == x[i+1,2:4])) {
        vec <- data.frame("chr"=x[["seqnames"]][i],
                          "acc"= x[["protein_id"]][(i+1)], 
                          "locus_tag" = x[["locus_tag"]][i], 
                          "gene_name" = x[["gene"]][i], 
                          "length_nt"= x[["width"]][i], 
                          "start"= x[["start"]][i], 
                          "end"= x[["end"]][i],
                          "strand" = x[["strand"]][i],
                          "type"= x[["gene_biotype"]][i],
                          "annotation"= x[["product"]][(i+1)],
                          stringsAsFactors = F)
        ids[[i]] <- vec
      }}}
  
  df <- bind_rows(ids)
  if (length(na.omit(unique(x[["protein_id"]]))) == length(na.omit(unique(df[["acc"]])))) {
    print("all protein ids accounted for")
  }
  return(df)
}
getNOG <- function(x, colname = "eggNOGs", sep = " ", factors = FALSE) {
  lis <- list()
  for (i in 1:length(x[[colname]])) {
    #vec <- c(capture.output(cat(unlist(str_extract_all(x[[colname]][i], "[a-zA-Z0-9]*@NOG")), sep = sep)))
    vec <- c(capture.output(cat(unlist(str_extract_all(x[[colname]][i], "[a-zA-Z0-9]*@")), sep = sep)))
    lis[[i]] <- vec
    column.names <- "NOGs"
  } 
  df <- as.data.frame(do.call(rbind, lis), stringsAsFactors = FALSE)
  colnames(df) <- column.names
  df <- cbind("query" = x[["query"]], df, "COG_category" = x[["COG Cat."]], "HMM_description" = x[["eggNOG HMM Desc."]])
  return(df)
}

#This function gets sequences based on your acc list
getseqs <- function(x, filename) {
  #x must be a vector of acc or gi numbers.
  require(reutils)
  path <- file.path(getwd(), filename)
  uid1 <- esearch(x, "protein", rettype = "uilist", usehistory = F)
  efetch(uid1, "protein", "fasta", "text", outfile = filename)
  print(paste("sequences can be found at ", path))
}

#TESTS FOR FUNCTIONAL ENRICHMENT
#performs hypergeometric test on provided subset of genes/proteins relative to the genome
#code adapted from Keely Dulmage (https://journals.asm.org/doi/10.1128/mBio.00649-15)
nogtest <- function(namelist,nogfile,pvalue, cutoff = 5) {
  #namelist is a vector of protein on gene names you wnat to test for enrichment
  #nogfile is the genome-wide GETNOG output
  #p-value is significance threshold desired
  #cutoff preset prevents functional categories with less than the designated number of genes/proteins being displayed 
  
  nogs <- nogfile[nogfile[["query"]] %in% namelist,]
  clust <-  table(nogs[["COG_category"]])
  resm <- matrix(0,length(clust),3) #create 0 matrix
  res <- data.frame(resm)  #make 0 df
  rownames(res) <- names(clust)
  colnames(res) <- c("probability", "expected","count")
  all <- table(nogfile[["COG_category"]][nogfile[["COG_category"]] %in% nogs[["COG_category"]]])
  tot <- sum(table(nogfile[["COG_category"]]))
  #print(tot); print(all); print(clust)
  for (i in 1:length(clust)){   #calc expected frequencies and pval by hypergeo and append to DF
    
    res[i,1] <- signif(phyper(clust[i], all[i], tot-all[i], nrow(nogs),lower.tail=F), digits = 4)
    res[i,2] <- signif(all[i]*(nrow(nogs)/tot), digits = 4)
    res[i,3] <- clust[i]
  }
  fin <- subset(res, probability <= pvalue & count >= cutoff)
  fin$COG_category <- rownames(fin)
  fin <- fin[, c("COG_category", "probability", "expected", "count")]
  return(fin)
}
#Use the following function to look at the genes in your cluster associated with a particular COG.
nogset= function(namelist,nogfile, cog.category) {
  subset(nogfile, is.element(nogfile$query, namelist) & is.element(nogfile$COG_category, cog.category)==TRUE)
}



### ===================== Input your DEG data ===================== ###
#read in your .gff file and list of significant genes here
my.genes <- read.csv("4_Functional_enrichment/Res2_1_pex10DGA1_vs_po1fEmpty_None_sig_up.csv")
my.gff <- as.data.frame(rtracklayer::import.gff("GCF_000002525.2_ASM252v1_genomic.gff"))
#execute function 
gff.table <- gff.parse(my.gff)
#convert from locus tags to GI accession numbers
my.accs <- gff.table[gff.table$locus_tag %in% my.genes$gene,]$acc
#write.csv(my.accs, "my.accs.csv")

#check length
length(my.accs) == length(my.genes$gene)
my.genes[my.genes$gene %ni% gff.table$locus_tag,]

#These differentially expressed genes do not exist in the gff table
my.genes$gene[my.genes$gene %ni% gff.table$locus_tag]

#get sequences : you can download at here too, https://www.ncbi.nlm.nih.gov/sites/batchentrez
#my.accs <- read.csv("my.accs_1.csv", header = F)

#getseqs(my.accs, paste(Sys.Date(), "_seqs.for.eggnog.faa", sep = ""))



### ===================== Parsing the eggNOG output file ===================== ###
my.eggNOG <- read_csv("GCF_000002525.2_ASM252v1_emapper.csv")
# How to make "emapper.csv": Input the ".faa file" from NCBI (Download sequences in FASTA format for protein) into EggNOG (http://eggnogdb.embl.de/#/app/emapper)
my.genome.NOGs <- getNOG(my.eggNOG)

#these genes were identified as differentially expressed, but are not present in the eggNOG database
gff.table[gff.table$acc %in% my.accs[my.accs %ni% my.genome.NOGs$query],]


#append unknown function annotation to proteins with missing COG classifications
my.genome.NOGs$COG_category[is.na(my.genome.NOGs$COG_category)] <- "S"
my.genome.NOGs$HMM_description[is.na(my.genome.NOGs$HMM_description)] <- "function unknown, manually assigned"
#calculate enrichment of all functional categories in the DEGs
my.genes.hypg <-  nogtest(my.accs, my.genome.NOGs, 1, cutoff = 1)
#correct for multiple testing, using FDR. 
my.genes.hypg$p.adj <- p.adjust(my.genes.hypg$probability, method = "fdr")
my.genes.hypg


#get proteins with specific function, where "X" is the single letter functional code
ls <- list()
counter <- 1
for (i in my.genes.hypg$COG_category){
  ls[[counter]] <- nogset(my.accs, my.genome.NOGs, i)
  counter <- counter +1
}
my.sig.proteins <- bind_rows(ls)
p.vals <- vector()
for (i in 1:length(my.sig.proteins$COG_category)) {
  cog <- my.sig.proteins$COG_category[i]
  p.vals[i] <- filter(my.genes.hypg, my.genes.hypg$COG_category == cog)$p.adj
}

#convert those proteins from acc back to locus_tags
my.sig.proteins$locus_tag <- gff.table$locus_tag[match(my.sig.proteins$query, gff.table$acc)]
my.sig.proteins$old_locus_tag <- gff.table$old_locus_tag[match(my.sig.proteins$locus_tag, gff.table$locus_tag)]
my.sig.proteins$log2FoldChange <- my.genes$log2FoldChange[match(my.sig.proteins$locus_tag, my.genes$gene)]
my.sig.proteins$COG.padj <- p.vals
(my.sig.proteins<- my.sig.proteins[c(5,6,1,4,2,3,7,8)])

#write out files
write_csv(my.sig.proteins, "4_Functional_enrichment/output/_func.csv")

my.sig.proteins.filter <- my.sig.proteins %>%     # Sorted in ascending order by adjusted p-value
  filter(COG.padj < 0.05)
#write out files
write_csv(my.sig.proteins.filter, "4_Functional_enrichment/output/_func_sig.csv")


#Bar graph
func <- data.frame(my.sig.proteins.filter)
query <- count(func, "COG_category")
query.data <- data.frame(query)
query.data$COG_category <- factor(query.data$COG_category, 
                                  level = query.data$COG_category[order(query.data$freq)])
myplot <- ggplot(query.data, aes(x=COG_category, y=freq)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label=freq), vjust=0)
myplot

ggsave("4_Functional_enrichment/output/_func_sig.png")

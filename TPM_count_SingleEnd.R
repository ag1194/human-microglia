rm(list = ls())

library("GenomicAlignments")
library("GenomicRanges")
library("rtracklayer")
library("Rsamtools")
library("GenomicFeatures")
library("org.Mm.eg.db")
library("org.Hs.eg.db")
library("BiocParallel")
library("matrixStats")
library("magrittr")
library("readr")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
register(MulticoreParam(workers = 12))
args = commandArgs(trailingOnly=TRUE)

### genome
geno =args[2] 
setwd(args[1])
bamDir =(paste0(getwd(),"/rawsam/"))


###read strand###
mystring <- read_file(paste0(getwd(),"/strand.txt"))  #create this file forward or reverse
mystring=gsub("\n","",mystring)
print(mystring)

### bam files from STAR results
files <- BamFileList(dir(bamDir, ".bam$", full=TRUE)) ##add 'asMates' for pair-end
names(files) <- basename(names(files))
names(files) =gsub(".Aligned.sortedByCoord.DeDup.bam","",names(files))
samplename=names(files)

#####################################calculate CDS TPM######################################
if(geno=="hg19")
{
  ####build genome annotation###
 	gtffile <- file.path("Homo_sapiens.GRCh38.97.gtf")
	txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character())
	gtf <- rtracklayer::import('Homo_sapiens.GRCh38.97.gtf')
	gtf_df=as.data.frame(gtf)
    IDDB <- org.Hs.eg.db

}

if(geno=="mm9")
{
  ####build genome annotation###
  gtffile <- file.path("gencode.vM20.primary_assembly.annotation.gtf")
	txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character())
	gtf <- rtracklayer::import('encode.vM20.primary_assembly.annotation.gtf')
	gtf_df=as.data.frame(gtf)
	IDDB <- org.Mm.eg.db

}

CDSbygene <- exonsBy(txdb, by="gene")

x = unlist(CDSbygene)
x1=as.data.frame(names(x))
colnames(x1)[1]="gene_id"
x1$gene_name= gtf_df$gene_name[match(x1$gene_id, gtf_df$gene_id)]
#x1[[1]] = gsub("\\..*","",x1[[1]])
names(x)=x1$gene_name


CDSbygene = split(x, names(x))
CDSbygene = reduce(CDSbygene)
cdslength=as.data.frame(sum(width(CDSbygene)))
names(cdslength)[1]="length"
cdslength$gene_symbol=rownames(cdslength)

if(mystring=="forward")
{
  CDStbl<- summarizeOverlaps(CDSbygene, files, mode="IntersectionNotEmpty", singleEnd=T, ignore.strand=F) ##added 'singleEnd=F' for pair-end;'ignore.strand=F' for strand_specific
} else if(mystring=="reverse"){
  CDStbl<- summarizeOverlaps(CDSbygene, files, mode="IntersectionNotEmpty", singleEnd=T, ignore.strand=F, preprocess.reads=invertStrand)
} else {
  CDStbl<- summarizeOverlaps(CDSbygene, files, mode="IntersectionNotEmpty", singleEnd=T, ignore.strand=T)
}


finalDF=as.data.frame(assays(CDStbl)$counts)

names(finalDF)=paste0(samplename,"_CDSreads")
readscols=names(finalDF)
finalDF=merge(finalDF, cdslength, by="row.names")
RPKcols=paste0(samplename,"_RPK")
finalDF[,RPKcols]=finalDF[,readscols]/(finalDF[,"length"]/1000)
TPMcols=paste0(samplename,"_TPM")
finalDF[TPMcols]=finalDF[RPKcols]/colSums(finalDF[RPKcols])*1000000
for (i in seq(1,length(TPMcols),1)){
finalDF[TPMcols][i]= finalDF[RPKcols][i]/sum(finalDF[RPKcols][i])*1000000
}
finalDF=finalDF[,c("gene_symbol", readscols, TPMcols,"length")]
write.table(finalDF, paste0(getwd(),'/rawout/',"all.Gene_reads_TPM.txt"), sep="\t",row.names=FALSE,quote=FALSE)




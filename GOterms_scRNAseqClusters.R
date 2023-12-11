
############################################### GOTERMS ################################################

table1<- read.csv("MSADPrkMsNew_allClusters.markers.csv")

library(org.Hs.eg.db)
hs <- org.Hs.eg.db

my.symbols <- as.character(table1 $gene)
myEntrez<- select(hs, 
       keys = my.symbols,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")
       
             

table1 $ ENTREZID<- myEntrez$ENTREZID[match(table1 $gene, myEntrez$SYMBOL)]


c0 = subset(table1, cluster==0)
c1 = subset(table1, cluster==1)
c2 = subset(table1, cluster==2)
c3 = subset(table1, cluster==3)
c4 = subset(table1, cluster==4)
c5 = subset(table1, cluster==5)
c6 = subset(table1, cluster==6)
c7 = subset(table1, cluster==7)
c8 = subset(table1, cluster==8)
c9 = subset(table1, cluster==9)
c10 = subset(table1, cluster==10)
c11 = subset(table1, cluster==11)
c12 = subset(table1, cluster==12)
c13 = subset(table1, cluster==13)

 x=merge(stabLess,stabLess1, by="gene_symbol")
 y=merge(stabMore, stabMore1, by="gene_symbol")
 q=merge(stabLess, stabMore1, by="gene_symbol")
 z=merge(stabMore, stabLess1, by="gene_symbol")
### 241 common for Less 270 common for More

library(GOstats)
#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db")
require(org.Hs.eg.db)

datatables=c("c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13")

#selectedEntrezIds = unique(stabMore1 $gene_id)   ###foreground
#entrezUniverse = unique(hesc2 $gene_id)  ###background


for (k in 1:length(datatables)) {
	data<-get(datatables[k])
	#dataBG<-get(datatablesBG[k])
	
selectedEntrezIds = unique(data $ENTREZID)   ###foreground
entrezUniverse = unique(table1 $ENTREZID)  ###background

print(length(selectedEntrezIds))
print(length(entrezUniverse))

###set parameters
annotation ="org.Hs.eg.db"
hgCutoff = 0.01
conditional = FALSE
categorySize=5
params.bp <- new("GOHyperGParams",
                   geneIds = as.character(selectedEntrezIds),
                   universeGeneIds = as.character(entrezUniverse),
                   annotation = annotation, ######################## was "hgu95av2.db" in the demo
                   ontology="BP",
                   pvalueCutoff=hgCutoff,
                   conditional=conditional, ######################## read the manual
                   testDirection="over")
#params.cc <- params.mf <- params.bp
#ontology(params.cc) = "CC"
#ontology(params.mf) = "MF"
params.cc <- new("GOHyperGParams",
                   geneIds = as.character(selectedEntrezIds),
                   universeGeneIds = as.character(entrezUniverse),
                   annotation = annotation, ######################## was "hgu95av2.db" in the demo
                   ontology="CC",
                   pvalueCutoff=hgCutoff,
                   conditional=conditional, ######################## read the manual
                   testDirection="over")
params.mf <- new("GOHyperGParams",
                   geneIds = as.character(selectedEntrezIds),
                   universeGeneIds = as.character(entrezUniverse),
                   annotation = annotation, ######################## was "hgu95av2.db" in the demo
                   ontology="MF",
                   pvalueCutoff=hgCutoff,
                   conditional=conditional, ######################## read the manual
                   testDirection="over")


############# BP
hgOver.bp <- hyperGTest(params.bp)
df4 = summary(hgOver.bp, categorySize=categorySize) # , categorySize is the minimum category size
if(nrow(df4) == 0){
	df4 = summary(hgOver.bp, categorySize=2)
if(nrow(df4) == 0){
	pvalueCutoff(params.bp) = 1
    hgOver.bp <- hyperGTest(params.bp)
    df4 = summary(hgOver.bp, categorySize=5)
  }
}

# filter out generic terms
df4 <- df4[df4$Size < 1000,]

  
names(df4)[1] = "GOID"
df4$Ontology = "BP"
# add gene symbols
df4$GenesFound = ""
geneIdList = geneIdsByCategory(hgOver.bp)[sigCategories(hgOver.bp)]

if(annotation == "org.Hs.eg.db"){
  library(org.Hs.eg.db)
  for(i in 1:nrow(df4)){
    df4$GenesFound[i] = paste(select(org.Hs.eg.db, keys = as.character(geneIdList[[df4$GOID[i]]]), columns = "SYMBOL")$SYMBOL, collapse = ", ")
    }
  }
  
if(annotation == "org.Mm.eg.db"){
  library(org.Mm.eg.db)
  for(i in 1:nrow(df4)){
    df4$GenesFound[i] = paste(select(org.Mm.eg.db, keys = as.character(geneIdList[[df4$GOID[i]]]), columns = "SYMBOL")$SYMBOL, collapse = ", ")
  }
}

# merge GO terms with significant overlaps, Wencheng's method
if(nrow(df4) > 1){
  df4$filter = 0
  for(i in 1:(nrow(df4)-1)){
    if(df4$filter[i]){next}
    for(j in (i+1):nrow(df4)){
      if(df4$filter[j]){next}
      # calculate fraction of overlapped genes between GOs
      overlap_i = length(intersect(geneIdList[[df4$GOID[i]]], geneIdList[[df4$GOID[j]]]))/length(geneIdList[[df4$GOID[i]]])
      overlap_j = length(intersect(geneIdList[[df4$GOID[i]]], geneIdList[[df4$GOID[j]]]))/length(geneIdList[[df4$GOID[j]]])
      # if the overlap is significant: 
      if(max(overlap_i, overlap_j) > 0.75 & min(overlap_i, overlap_j) > 0.25){
        # if p values are similar, keep the GO with smaller universe gene list
        if(abs(log2(df4$Pvalue[i]/df4$Pvalue[j])) < 1){
          index = ifelse(length(geneIdUniverse(hgOver.bp)[[df4$GOID[i]]]) < length(geneIdUniverse(hgOver.bp)[[df4$GOID[j]]]), j, i)
          # if GO universe gene list sizes are the same, keep the more significant one.
          if(length(geneIdUniverse(hgOver.bp)[[df4$GOID[i]]]) == length(geneIdUniverse(hgOver.bp)[[df4$GOID[j]]])){
            index = ifelse(df4$Pvalue[i] < df4$Pvalue[j], j, i)
          }
        }else{
          index = ifelse(df4$Pvalue[i] < df4$Pvalue[j], j, i) 
        }
        df4$filter[index] = 1
      }  
    }
  }
  df4 = df4[df4$filter == 0, ]
  df4$filter = NULL
}



###CC

hgOver.cc <- hyperGTest(params.cc)
df5 = summary(hgOver.cc, categorySize=categorySize) # , categorySize is the minimum category size
if(nrow(df5) == 0){
	df5 = summary(hgOver.cc, categorySize=2)
if(nrow(df5) == 0){
	pvalueCutoff(params.cc) = 1
    hgOver.cc <- hyperGTest(params.cc)
    df5 = summary(hgOver.cc, categorySize=5)
  }
}

# filter out generic terms
df5 <- df5[df5 $Size < 1000,]

  
names(df5)[1] = "GOID"
df5$Ontology = "CC"
# add gene symbols
df5$GenesFound = ""
geneIdList = geneIdsByCategory(hgOver.cc)[sigCategories(hgOver.cc)]

if(annotation == "org.Hs.eg.db"){
  library(org.Hs.eg.db)
  for(i in 1:nrow(df5)){
    df5$GenesFound[i] = paste(select(org.Hs.eg.db, keys = as.character(geneIdList[[df5$GOID[i]]]), columns = "SYMBOL")$SYMBOL, collapse = ", ")
    }
  }
  
if(annotation == "org.Mm.eg.db"){
  library(org.Mm.eg.db)
  for(i in 1:nrow(df5)){
    df5$GenesFound[i] = paste(select(org.Mm.eg.db, keys = as.character(geneIdList[[df5$GOID[i]]]), columns = "SYMBOL")$SYMBOL, collapse = ", ")
  }
}

# merge GO terms with significant overlaps, Wencheng's method
if(nrow(df5) > 1){
  df5$filter = 0
  for(i in 1:(nrow(df5)-1)){
    if(df5$filter[i]){next}
    for(j in (i+1):nrow(df5)){
      if(df5$filter[j]){next}
      # calculate fraction of overlapped genes between GOs
      overlap_i = length(intersect(geneIdList[[df5$GOID[i]]], geneIdList[[df5$GOID[j]]]))/length(geneIdList[[df5$GOID[i]]])
      overlap_j = length(intersect(geneIdList[[df5$GOID[i]]], geneIdList[[df5$GOID[j]]]))/length(geneIdList[[df5$GOID[j]]])
      # if the overlap is significant: 
      if(max(overlap_i, overlap_j) > 0.75 & min(overlap_i, overlap_j) > 0.25){
        # if p values are similar, keep the GO with smaller universe gene list
        if(abs(log2(df5$Pvalue[i]/df5$Pvalue[j])) < 1){
          index = ifelse(length(geneIdUniverse(hgOver.cc)[[df5$GOID[i]]]) < length(geneIdUniverse(hgOver.cc)[[df5$GOID[j]]]), j, i)
          # if GO universe gene list sizes are the same, keep the more significant one.
          if(length(geneIdUniverse(hgOver.cc)[[df5$GOID[i]]]) == length(geneIdUniverse(hgOver.cc)[[df5$GOID[j]]])){
            index = ifelse(df5$Pvalue[i] < df5$Pvalue[j], j, i)
          }
        }else{
          index = ifelse(df5$Pvalue[i] < df5$Pvalue[j], j, i) 
        }
        df5$filter[index] = 1
      }  
    }
  }
  df5 = df5[df5$filter == 0, ]
  df5$filter = NULL
}



###MF
hgOver.mf <- hyperGTest(params.mf)
df6 = summary(hgOver.mf, categorySize=categorySize) # , categorySize is the minimum category size
if(nrow(df6) == 0){
	df6 = summary(hgOver.mf, categorySize=2)
if(nrow(df6) == 0){
	pvalueCutoff(params.bp) = 1
    hgOver.mf <- hyperGTest(params.bp)
    df6 = summary(hgOver.mf, categorySize=5)
  }
}
 
# filter out generic terms
df6 <- df6[df6$Size < 1000,]
 
 
names(df6)[1] = "GOID"
df6$Ontology = "MF"
# add gene symbols
df6$GenesFound = ""
geneIdList = geneIdsByCategory(hgOver.mf)[sigCategories(hgOver.mf)]
if(annotation == "org.Hs.eg.db"){
  library(org.Hs.eg.db)
  for(i in 1:nrow(df6)){
    df6$GenesFound[i] = paste(select(org.Hs.eg.db, keys = as.character(geneIdList[[df6$GOID[i]]]), columns = "SYMBOL")$SYMBOL, collapse = ", ")
    }
  }
  
if(annotation == "org.Mm.eg.db"){
  library(org.Mm.eg.db)
  for(i in 1:nrow(df6)){
    df6$GenesFound[i] = paste(select(org.Mm.eg.db, keys = as.character(geneIdList[[df6$GOID[i]]]), columns = "SYMBOL")$SYMBOL, collapse = ", ")
  }
}

# merge GO terms with significant overlaps, Wencheng's method
if(nrow(df6) > 1){
  df6$filter = 0
  for(i in 1:(nrow(df6)-1)){
    if(df6$filter[i]){next}
    for(j in (i+1):nrow(df6)){
      if(df6$filter[j]){next}
      # calculate fraction of overlapped genes between GOs
      overlap_i = length(intersect(geneIdList[[df6$GOID[i]]], geneIdList[[df6$GOID[j]]]))/length(geneIdList[[df6$GOID[i]]])
      overlap_j = length(intersect(geneIdList[[df6$GOID[i]]], geneIdList[[df6$GOID[j]]]))/length(geneIdList[[df6$GOID[j]]])
      # if the overlap is significant: 
      if(max(overlap_i, overlap_j) > 0.75 & min(overlap_i, overlap_j) > 0.25){
        # if p values are similar, keep the GO with smaller universe gene list
        if(abs(log2(df6$Pvalue[i]/df6$Pvalue[j])) < 1){
          index = ifelse(length(geneIdUniverse(hgOver.mf)[[df6$GOID[i]]]) < length(geneIdUniverse(hgOver.mf)[[df6$GOID[j]]]), j, i)
          # if GO universe gene list sizes are the same, keep the more significant one.
          if(length(geneIdUniverse(hgOver.mf)[[df6$GOID[i]]]) == length(geneIdUniverse(hgOver.mf)[[df6$GOID[j]]])){
            index = ifelse(df6$Pvalue[i] < df6$Pvalue[j], j, i)
          }
        }else{
          index = ifelse(df6$Pvalue[i] < df6$Pvalue[j], j, i) 
        }
        df6$filter[index] = 1
      }  
    }
  }
  df6 = df6[df6$filter == 0, ]
  df6$filter = NULL
}


                               
#Now Merge all three GO

df_f = rbind(df4, df5, df6)
df_f$BgCount = df_f$Size - df_f $Count
df_2 = df_f[, c("Ontology", "GOID", "Pvalue", "OddsRatio", "Count", "BgCount", "GenesFound", "Term")]
colnames(df_2)[5] = "FgCount"



#write.csv(df_2, paste("GOTerms_",datatables[i],"APAstabMarkers_D>P.csv", sep=""), row.names = FALSE)
#write.csv(df_2, paste("GOTerms_",datatables[k],"APAstabMarkers_P>D.csv", sep=""), row.names = FALSE)

#write.csv(df_2, paste("GOTerms_",datatables[k],"3UTR_APA.csv", sep=""), row.names = FALSE)
write.csv(df_2, paste("./GOterms/GOTerms_",datatables[k],"cluster.csv", sep=""), row.names = FALSE)

}

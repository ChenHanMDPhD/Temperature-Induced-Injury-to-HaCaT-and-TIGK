#### Packages ####
#### Packages for DESeq2, DEA, and Enrichment 
install.packages("ggVennDiagram")
install.packages("gplots")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ReactomePA")
BiocManager::install("limma")
BiocManager::install("GEOquery")
BiocManager::install("umap")
BiocManager::install("illuminaHumanv4.db")
BiocManager::install("biomaRt")
BiocManager::install("MethReg")
devtools::install_github("stephenturner/annotables")
install_github("wjawaid/enrichR", force=TRUE)
install.packages("cowplot")
install.packages("scales")
install.packages("corrr")
install.packages("ggcorrplot")
install.packages("FactoMineR")

library(scales)
library(cowplot)
library(devtools)
library(GEOquery)
library(limma)
library(umap)
library(DESeq2)
library("ggVennDiagram")
library(gplots)
library(BiocManager)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library("ReactomePA")
library(AnnotationDbi)
library(VennDiagram)
library(ggplot2)
library(openxlsx)
library(readxl)
library(stringr)
library(tidyverse)
library("illuminaHumanv4.db") #Get this library if you don't have - converts illumina probes into gene
library(biomaRt)
library("MethReg")
library(annotables)
library(enrichR)
library(ComplexHeatmap)
library("org.Mm.eg.db")
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
dbs <- c("GO_Biological_Process_2023","GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "Reactome_2022")
library(stringr)
library(dplyr)
##### correlation matrix packages
library("FactoMineR")
library(ggcorrplot)
library('corrr')
library("reshape2")
##### ggplot packages for PCA plots
# install.packages("factoextra")
library(ggplot2)
library(ggplotify)
library(ggrepel)
library(factoextra)
library(ggforce)

#### Reading in the data and filtering only "protein_coding" ####
Thermal.counts<-read.csv("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/gene_counts_protein coding.csv")
#Only reading mRNA/protein coding
Thermal.counts<-Thermal.counts%>%filter(gene_biotype == "protein_coding")
rownames(Thermal.counts)<-Thermal.counts[,1]
colnames(Thermal.counts)

#####PCA plot with ggplot ####
coldata.Thermal.PCA<-matrix(c("H60C_1","H60C_2","H60C_3","T60C_1","T60C_2","T60C_3",
                              "HCtrl_1","HCtrl_3","HCtrl_3","TCtrl_1","TCtrl_2","TCtrl_3",
                              "Hn25C_1","Hn25C_2","Hn25C_3","Tn25C_1","Tn25C_2","Tn25C_3",
                              "H_Heat","H_Heat","H_Heat","T_Heat","T_Heat","T_Heat",
                              "H_ctrl","H_ctrl","H_ctrl","T_ctrl","T_ctrl","T_ctrl",
                              "H_cold","H_cold","H_cold","T_cold","T_cold","T_cold",
                              "Heat","Heat","Heat","Heat","Heat","Heat",
                              "Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl", 
                              "Cold","Cold","Cold","Cold","Cold","Cold",
                              "HaCaT","HaCaT","HaCaT","TIGK","TIGK","TIGK",
                              "HaCaT","HaCaT","HaCaT","TIGK","TIGK","TIGK",
                              "HaCaT","HaCaT","HaCaT","TIGK","TIGK","TIGK"),
                            nrow = 18,
                            ncol = 4,
                            byrow = FALSE)
colnames(coldata.Thermal.PCA)<-c("sample","SampleTemp","Temp","Group")

Thermal.counts2<-Thermal.counts[, c(1,2,3,4,5,6,7,14,15,16,17,18,19,8,9,10,11,12,13,20,21,22,23,24,25,26,27,28)] #reordering so Ctrl is in center
colnames(Thermal.counts2)
dds.Thermal.PCA2<-DESeqDataSetFromMatrix(countData = Thermal.counts2[,2:19],
                                        colData = coldata.Thermal.PCA,
                                        design = ~ SampleTemp)

#filtering all datasets: keep only avg counts >=10
keep<-rowMeans(counts(dds.Thermal.PCA2))>=10 #could do rowsums too
dds.Thermal.PCA2<-dds.Thermal.PCA2[keep,] 

#transformation of data
rld.Thermal.PCA2<-vst(dds.Thermal.PCA2, blind=FALSE)

pcaData2<- plotPCA(rld.Thermal.PCA2, intgroup=c("sample","SampleTemp","Temp","Group"), returnData=TRUE)

#Making the plots
pcaData2.reoder <-  factor(pcaData2$SampleTemp, levels = c("H_Heat", "H_ctrl", "H_cold","T_Heat", "T_ctrl", "T_cold"))
scale_map <- data.frame(name = c("HaCaT_Heat", "HaCaT_Control","HaCaT_Cold", "TIGK_Heat","TIGK_Control", "TIGK_Cold"),
                        color = c("red", "chartreuse4", "blue", "tomato1","chartreuse","cyan"),
                        shape = c(15, 16, 17,15, 16, 17))
ggplot(pcaData2, aes(x = PC1, y = PC2,shape = pcaData2.reoder, colour = pcaData2.reoder)) + 
  geom_point(size=4) +
  scale_color_manual(values = scale_map$color, labels = scale_map$name,
                   name = "Samples") + 
  scale_shape_manual(values = scale_map$shape, labels = scale_map$name,
                     name = "Samples")

#####Generating averaged correlation matrix ####
DF<-Thermal.counts[,2:19]
keep<-rowMeans((DF))>=10 #could do rowsums too
DF<-DF[keep,] 
nrow(DF)
colnames(DF)
prefix <- unique(unlist(strsplit(names(DF), "//_[0-9]")))
DF<-sapply(prefix, function(i) rowMeans(DF[, grepl(i, names(DF))]))

corr.avg <- round(cor(DF[,1:6]), 2) #now we have a correlation matrix where we have only groups and indvidual samples
upper_tri.avg<-get_upper_tri(corr.avg)
upper_tri.avg
melted_corr.avg<-melt(upper_tri.avg)

ggplot(data = melted_corr.avg, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+geom_text(aes(label=value), size = 8)+
  scale_fill_gradient(
    name = "Cor", # changes legend title
    low = "blue",
    high = "red",
    limit = c(min(melted_corr.avg$value), 1),
    space = "Lab",
    guide = "colourbar",
    na.value="white"
  )+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

#####Comparing HaCaT and TIGK to primary cell lines (NHEK and NHOK) ####
#Counts for HaCaT and TIGK
Thermal.counts<-Thermal.counts%>%filter(gene_biotype == "protein_coding")
rownames(Thermal.counts)<-Thermal.counts[,1]
colnames(Thermal.counts)
Thermal.counts.HaCaT.Ctrl<-Thermal.counts[,c(14:16,20)]
Thermal.counts.TIGK.Ctrl<-Thermal.counts[,17:20]
Thermal.counts.HaCaT.Ctrl$Average<-rowMeans(Thermal.counts.HaCaT.Ctrl[,1:3])
Thermal.counts.TIGK.Ctrl$Average<-rowMeans(Thermal.counts.TIGK.Ctrl[,1:3])

Thermal.counts.HaCaT.Ctrl<-Thermal.counts.HaCaT.Ctrl%>%filter(Average >= 10)
Thermal.counts.HaCaT.Ctrl<-Thermal.counts.HaCaT.Ctrl[order(-Thermal.counts.HaCaT.Ctrl$Average),]
head(Thermal.counts.HaCaT.Ctrl)

Thermal.counts.TIGK.Ctrl<-Thermal.counts.TIGK.Ctrl%>%filter(Average >= 10)
Thermal.counts.TIGK.Ctrl<-Thermal.counts.TIGK.Ctrl[order(-Thermal.counts.TIGK.Ctrl$Average),]
head(Thermal.counts.TIGK.Ctrl)

Thermal.counts.HaCaT.Ctrl<-Thermal.counts.HaCaT.Ctrl[1:5000,]
Thermal.counts.TIGK.Ctrl<-Thermal.counts.TIGK.Ctrl[1:5000,]

nrow(Thermal.counts.HaCaT.Ctrl)
nrow(Thermal.counts.TIGK.Ctrl)

#Counts for NHEK (GSE184119)
NHEK.GSE184119<-read.csv("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/HaCaT vs TIGK manuscript/Wietecha Related Work/GSE184119-raw counts_NHEK.csv")
NHEK.GSE184119<-NHEK.GSE184119%>%filter(Average >= 10)
NHEK.GSE184119<- bitr(NHEK.GSE184119$GeneID, fromType = "ENTREZID",
            toType = c("SYMBOL"),
            OrgDb = org.Hs.eg.db) #not every gene mapped but that's ok
NHEK.GSE184119<-na.omit(NHEK.GSE184119)
NHEK.GSE184119<-NHEK.GSE184119[1:5000,] #top5000

#Counts for NHEK (GSE185309)
NHEK.GSE185309<-read.csv("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/GSE185309_NHEK_IH6000/GSE185309.csv")
NHEK.GSE185309<-NHEK.GSE185309%>%filter(Average >= 10)
NHEK.GSE185309<- bitr(NHEK.GSE185309$GeneID, fromType = "ENTREZID",
                      toType = c("SYMBOL"),
                      OrgDb = org.Hs.eg.db) #not every gene mapped but that's ok
NHEK.GSE185309<-na.omit(NHEK.GSE185309)
NHEK.GSE185309<-NHEK.GSE185309[1:5000,] #top5000

#Counts for NHOK-1 (GSE262505)
NHOK.GSE262505<-read.csv("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/HaCaT vs TIGK manuscript/Wietecha Related Work/GSE262505_raw counts_NHOK.csv")
NHOK.GSE262505<-NHOK.GSE262505%>%filter(Average >= 10)
NHOK.GSE262505<- bitr(NHOK.GSE262505$GeneID, fromType = "ENTREZID",
            toType = c("SYMBOL"),
            OrgDb = org.Hs.eg.db) #not every gene mapped but that's ok
NHOK.GSE262505<-na.omit(NHOK.GSE262505)
NHOK.GSE262505<-NHOK.GSE262505[1:5000,] #top5000

#Counts for NHOK (GSE121627)
NHOK.GSE121627<-read.csv("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/GSE121627_NHOK_IH2000/GSE121627.csv")
NHOK.GSE121627<-NHOK.GSE121627%>%filter(Average >= 10)
NHOK.GSE121627<- bitr(NHOK.GSE121627$GeneID, fromType = "ENTREZID",
                      toType = c("SYMBOL"),
                      OrgDb = org.Hs.eg.db) #not every gene mapped but that's ok
NHOK.GSE121627<-na.omit(NHOK.GSE121627)
NHOK.GSE121627<-NHOK.GSE121627[1:5000,] #top5000

#saving the sheets
write.csv(NHEK.GSE184119, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/NHEK_NHOK_HaCaT_TIGK_Raw Count/NHEK.GSE184119.T5000.csv" )
write.csv(NHOK.GSE262505, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/NHEK_NHOK_HaCaT_TIGK_Raw Count/NHOK.GSE262505.T5000.csv" )
write.csv(Thermal.counts.HaCaT.Ctrl, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/NHEK_NHOK_HaCaT_TIGK_Raw Count/HaCaT.T5000.csv" )
write.csv(Thermal.counts.TIGK.Ctrl, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/NHEK_NHOK_HaCaT_TIGK_Raw Count/TIGK.T5000.csv" )
write.csv(NHEK.GSE185309, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/NHEK_NHOK_HaCaT_TIGK_Raw Count/NHEK.GSE185309.T5000.csv")
write.csv(NHOK.GSE121627, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/NHEK_NHOK_HaCaT_TIGK_Raw Count/NHOK.GSE121627.T5000.csv")


#Enrichment of top 1000 genes
enriched.TIGK.Top1000 <- enrichr(Thermal.counts.TIGK.Ctrl[1:1000,]$gene_name, dbs)
enriched.HaCaT.Top1000 <- enrichr(Thermal.counts.HaCaT.Ctrl[1:1000,]$gene_name, dbs)
enriched.NHEK.GSE184119.Top1000  <- enrichr(NHEK.GSE184119[1:1000,]$SYMBOL, dbs)
enriched.NHOK.GSE262505.Top1000  <- enrichr(NHOK.GSE262505[1:1000,]$SYMBOL, dbs)
enriched.NHOK.GSE121627.Top1000  <- enrichr(NHOK.GSE121627[1:1000,]$SYMBOL, dbs)
enriched.NHEK.GSE185309.Top1000  <- enrichr(NHEK.GSE185309[1:1000,]$SYMBOL, dbs)


Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP TIGK")
addWorksheet(Enrichments, "Reactome TIGK")
addWorksheet(Enrichments, "BP NHOKGSE262505")
addWorksheet(Enrichments, "Reactome NHOKGSE262505")
addWorksheet(Enrichments, "BP NHOK.GSE121627")
addWorksheet(Enrichments, "Reactome NHOK.GSE121627")

addWorksheet(Enrichments, "BP HaCaT")
addWorksheet(Enrichments, "Reactome HaCaT")
addWorksheet(Enrichments, "BP NHEKGSE184119")
addWorksheet(Enrichments, "Reactome NHEKGSE184119")
addWorksheet(Enrichments, "BP NHEK.GSE185309")
addWorksheet(Enrichments, "Reactome NHEK.GSE185309")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP TIGK", x =(as.data.frame(enriched.TIGK.Top1000$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome TIGK", x =(as.data.frame(enriched.TIGK.Top1000$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "BP NHOKGSE262505", x =(as.data.frame(enriched.NHOK.GSE262505.Top1000$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome NHOKGSE262505", x =(as.data.frame(enriched.NHOK.GSE262505.Top1000$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "BP NHOK.GSE121627",  x =(as.data.frame(enriched.NHOK.GSE121627.Top1000$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome NHOK.GSE121627", x =(as.data.frame(enriched.NHOK.GSE121627.Top1000$Reactome))%>%filter(Adjusted.P.value < 0.05))

writeData(Enrichments, sheet = "BP HaCaT", x =(as.data.frame(enriched.HaCaT.Top1000$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome HaCaT", x =(as.data.frame(enriched.HaCaT.Top1000$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "BP NHEKGSE184119", x =(as.data.frame(enriched.NHEK.GSE184119.Top1000$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome NHEKGSE184119", x =(as.data.frame(enriched.NHEK.GSE184119.Top1000$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "BP NHEK.GSE185309", x =(as.data.frame(enriched.NHEK.GSE185309.Top1000$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome NHEK.GSE185309", x =(as.data.frame(enriched.NHEK.GSE185309.Top1000$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
# Export the file
saveWorkbook(Enrichments,  file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/EnrichR Terms_NHEK_NHOK_HaCaT_TIGK_T1000.xlsx")

#Enrichment of the different genes between TIGK and NHOK and HaCaT and NHEK
#TIGK vs NHOK GSE262505
TvsN262505<-read_excel("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/NHEK_NHOK_HaCaT_TIGK_Raw Count/Gene Comparisons/Gene Comparisons for manuscript.xlsx",
           sheet = "TIGK vs NHOK (GSE262505)")
TvsN262505<-na.omit(TvsN262505)
EnrichR.TvsN262505.TIGKOnly<-enrichr(TvsN262505$A, dbs)

#TIGK vs NHOK GSE121627
TvsN121627<-read_excel("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/NHEK_NHOK_HaCaT_TIGK_Raw Count/Gene Comparisons/Gene Comparisons for manuscript.xlsx",
                       sheet = "TIGK vs NHOK (GSE121627)")
TvsN121627<-na.omit(TvsN121627)
EnrichR.TvsN121627.TIGKOnly<-enrichr(TvsN121627$A, dbs)

#HaCaT vs NHEK GSE184119
HvsN184119<-read_excel("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/NHEK_NHOK_HaCaT_TIGK_Raw Count/Gene Comparisons/Gene Comparisons for manuscript.xlsx",
                       sheet = "HaCaT vs NHEK (GSE184119)")
HvsN184119<-na.omit(HvsN184119)
EnrichR.HvsN184119.HaCaTOnly<-enrichr(HvsN184119$A, dbs)

#HaCaT vs NHEK GSE185309
HvsN185309<-read_excel("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/NHEK_NHOK_HaCaT_TIGK_Raw Count/Gene Comparisons/Gene Comparisons for manuscript.xlsx",
                       sheet = "HaCaT vs NHEK (GSE185309)")
HvsN185309<-na.omit(HvsN185309)
EnrichR.HvsN185309.HaCaTOnly<-enrichr(HvsN185309$A, dbs)

Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP TvsN262505")
addWorksheet(Enrichments, "Reactome TvsN262505")
addWorksheet(Enrichments, "BP TvsN121627")
addWorksheet(Enrichments, "Reactome TvsN121627")

addWorksheet(Enrichments, "BP TvsN184119")
addWorksheet(Enrichments, "Reactome TvsN184119")
addWorksheet(Enrichments, "BP Tvs185309")
addWorksheet(Enrichments, "Reactome Tvs185309")

# Write the data to the sheets
writeData(Enrichments, sheet = "BP TvsN262505", x =(as.data.frame(EnrichR.TvsN262505.TIGKOnly$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome TvsN262505", x =(as.data.frame(EnrichR.TvsN262505.TIGKOnly$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "BP TvsN121627", x =(as.data.frame(EnrichR.TvsN121627.TIGKOnly$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome TvsN121627", x =(as.data.frame(EnrichR.TvsN121627.TIGKOnly$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))

writeData(Enrichments, sheet = "BP HvsN184119", x =(as.data.frame(EnrichR.HvsN184119.HaCaTOnly$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome HvsN184119", x =(as.data.frame(EnrichR.HvsN184119.HaCaTOnly$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "BP Hvs185309", x =(as.data.frame(EnrichR.HvsN185309.HaCaTOnly$GO_Biological_Process_2023))%>%filter(Adjusted.P.value < 0.05))
writeData(Enrichments, sheet = "Reactome Hvs185309", x =(as.data.frame(EnrichR.HvsN185309.HaCaTOnly$Reactome_2022))%>%filter(Adjusted.P.value < 0.05))

# Export the file
saveWorkbook(Enrichments,  file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/NHEK_NHOK_HaCaT_TIGK_Raw Count/Enrichment for TIGK or HaCaT only.xlsx")

#Visualization for the enrichment of the different genes between TIGK and NHOK and HaCaT and NHEK
#TIGK vs GSE262505
EnrichR.TvsN262505.TIGKOnly.Reactome<-as.data.frame(EnrichR.TvsN262505.TIGKOnly$Reactome_2022)[1:10,]
frac <- EnrichR.TvsN262505.TIGKOnly.Reactome$Overlap
EnrichR.TvsN262505.TIGKOnly.Reactome$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))
p1<-ggplot(EnrichR.TvsN262505.TIGKOnly.Reactome, aes(x=EnrichR.TvsN262505.TIGKOnly.Reactome$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.TvsN262505.TIGKOnly.Reactome$Adjusted.P.value, size=(EnrichR.TvsN262505.TIGKOnly.Reactome$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png( "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/TIGK vs GSE262505_TIGK Only Reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

#TIGK vs GSE121627
EnrichR.TvsN121627.TIGKOnly.reactome<-as.data.frame(EnrichR.TvsN121627.TIGKOnly$Reactome_2022)[1:10,]
frac <- EnrichR.TvsN121627.TIGKOnly.reactome$Overlap
EnrichR.TvsN121627.TIGKOnly.reactome$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))
p1<-ggplot(EnrichR.TvsN121627.TIGKOnly.reactome, aes(x=EnrichR.TvsN121627.TIGKOnly.reactome$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.TvsN121627.TIGKOnly.reactome$Adjusted.P.value, size=(EnrichR.TvsN121627.TIGKOnly.reactome$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png( "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/TIGK vs GSE121627_TIGK Only Reactome.png",
     width = 1600,
     height = 2000,
     units = "px",
     res = 250)
p1
dev.off()


#TIGK vs GSE184119
EnrichR.HvsN184119.HaCaTOnly.reactome<-as.data.frame(EnrichR.HvsN184119.HaCaTOnly$Reactome_2022)[1:10,]
frac <- EnrichR.HvsN184119.HaCaTOnly.reactome$Overlap
EnrichR.HvsN184119.HaCaTOnly.reactome$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))
p1<-ggplot(EnrichR.HvsN184119.HaCaTOnly.reactome, aes(x=EnrichR.HvsN184119.HaCaTOnly.reactome$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.HvsN184119.HaCaTOnly.reactome$Adjusted.P.value, size=(EnrichR.HvsN184119.HaCaTOnly.reactome$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png( "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/HaCaT vs GSE184119_HaCaT Only Reactome.png",
     width = 1600,
     height = 2000,
     units = "px",
     res = 250)
p1
dev.off()

#HaCaT vs GSE185309
EnrichR.HvsN185309.HaCaTOnly.reactome<-as.data.frame(EnrichR.TvsN121627.TIGKOnly$Reactome_2022)[1:10,]
frac <- EnrichR.HvsN185309.HaCaTOnly.reactome$Overlap
EnrichR.HvsN185309.HaCaTOnly.reactome$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))
p1<-ggplot(EnrichR.HvsN185309.HaCaTOnly.reactome, aes(x=EnrichR.HvsN185309.HaCaTOnly.reactome$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=EnrichR.HvsN185309.HaCaTOnly.reactome$Adjusted.P.value, size=(EnrichR.HvsN185309.HaCaTOnly.reactome$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Overlap', y=NULL,
       color='Adjusted.P.value',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png( "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/HaCaT vs GSE185309_HaCaT Only Reactome.png",
     width = 1600,
     height = 2000,
     units = "px",
     res = 250)
p1
dev.off()




###### Differential Expression Analysis using DESeq2 and count matrix - Same time point comparisons ####
#setting the columns to measure for DESeq2
colnames(Thermal.counts)
comp.Heat <- c(2:7)
comp.Cold <- c(8:13)
comp.Ctrl <- c(14:19)
Comp.HaCaT.HeatvsCtrl <- c(2:4, 14:16)
Comp.TIGK.HeatvsCtrl <- c(5:7, 17:19)
Comp.HaCaT.ColdvsCtrl <- c(8:10, 14:16)
Comp.TIGK.ColdvsCtrl <- c(11:13, 17:19)

#column data already set from above - will repost it here
coldata.Thermal.PCA<-matrix(c("H60C_1","H60C_2","H60C_3","T60C_1","T60C_2","T60C_3",
                              "Hn25C_1","Hn25C_2","Hn25C_3","Tn25C_1","Tn25C_2","Tn25C_3",
                              "HCtrl_1","HCtrl_3","HCtrl_3","TCtrl_1","TCtrl_2","TCtrl_3",
                              "H_Heat","H_Heat","H_Heat","T_Heat","T_Heat","T_Heat",
                              "H_cold","H_cold","H_cold","T_cold","T_cold","T_cold",
                              "H_ctrl","H_ctrl","H_ctrl","T_ctrl","T_ctrl","T_ctrl",
                              "Heat","Heat","Heat","Heat","Heat","Heat",
                              "Cold","Cold","Cold","Cold","Cold","Cold",
                              "Ctrl","Ctrl","Ctrl","Ctrl","Ctrl","Ctrl"),
                            nrow = 18,
                            ncol = 3,
                            byrow = FALSE)
colnames(coldata.Thermal.PCA)<-c("sample","SampleTemp","Temp")

#DESeq2 analysis
dds.Heat.Comparison<-DESeqDataSetFromMatrix(countData = Thermal.counts[,comp.Heat],
                                    colData = coldata.Thermal.PCA[1:6,],
                                    design = ~ SampleTemp)

dds.Cold.Comparison<-DESeqDataSetFromMatrix(countData = Thermal.counts[,comp.Cold],
                                            colData = coldata.Thermal.PCA[7:12,],
                                            design = ~ SampleTemp)

dds.Ctrl.Comparison<-DESeqDataSetFromMatrix(countData = Thermal.counts[,comp.Ctrl],
                                            colData = coldata.Thermal.PCA[13:18,],
                                            design = ~ SampleTemp)

HaCaT.HeatvsCtrl<-c(1:3, 13:15)
dds.HaCaT.HeatvsCtrl.Comparison<-DESeqDataSetFromMatrix(countData = Thermal.counts[,Comp.HaCaT.HeatvsCtrl],
                                            colData = coldata.Thermal.PCA[HaCaT.HeatvsCtrl,],
                                            design = ~ SampleTemp)

HaCaT.ColdvsCtrl<-c(7:9, 13:15)
dds.HaCaT.ColdvsCtrl.Comparison<-DESeqDataSetFromMatrix(countData = Thermal.counts[,Comp.HaCaT.ColdvsCtrl],
                                                        colData = coldata.Thermal.PCA[HaCaT.ColdvsCtrl,],
                                                        design = ~ SampleTemp)

TIGK.HeatvsCtrl<-c(4:6, 16:18)
dds.TIGK.HeatvsCtrl.Comparison<-DESeqDataSetFromMatrix(countData = Thermal.counts[,Comp.TIGK.HeatvsCtrl],
                                                        colData = coldata.Thermal.PCA[TIGK.HeatvsCtrl,],
                                                        design = ~ SampleTemp)

TIGK.ColdvsCtrl<-c(10:12, 16:18)
dds.TIGK.ColdvsCtrl.Comparison<-DESeqDataSetFromMatrix(countData = Thermal.counts[,Comp.TIGK.ColdvsCtrl],
                                                        colData = coldata.Thermal.PCA[TIGK.ColdvsCtrl,],
                                                        design = ~ SampleTemp)
#filtering all datasets: keep only avg counts >=10
keep<-rowMeans(counts(dds.Heat.Comparison))>=10 
dds.Heat.Comparison<-dds.Heat.Comparison[keep,] 

keep<-rowMeans(counts(dds.Cold.Comparison))>=10
dds.Cold.Comparison<-dds.Cold.Comparison[keep,]

keep<-rowMeans(counts(dds.Ctrl.Comparison))>=10 
dds.Ctrl.Comparison<-dds.Ctrl.Comparison[keep,] 

keep<-rowMeans(counts(dds.HaCaT.HeatvsCtrl.Comparison))>=10 
dds.HaCaT.HeatvsCtrl.Comparison<-dds.HaCaT.HeatvsCtrl.Comparison[keep,] 

keep<-rowMeans(counts(dds.HaCaT.ColdvsCtrl.Comparison))>=10
dds.HaCaT.ColdvsCtrl.Comparison<-dds.HaCaT.ColdvsCtrl.Comparison[keep,]

keep<-rowMeans(counts(dds.TIGK.HeatvsCtrl.Comparison))>=10 
dds.TIGK.HeatvsCtrl.Comparison<-dds.TIGK.HeatvsCtrl.Comparison[keep,] 

keep<-rowMeans(counts(dds.TIGK.ColdvsCtrl.Comparison))>=10
dds.TIGK.ColdvsCtrl.Comparison<-dds.TIGK.ColdvsCtrl.Comparison[keep,]

#set the factor level (reference)
dds.Heat.Comparison$SampleTemp<-relevel(dds.Heat.Comparison$SampleTemp, ref = "H_Heat")
dds.Cold.Comparison$SampleTemp<-relevel(dds.Cold.Comparison$SampleTemp, ref = "H_cold")
dds.Ctrl.Comparison$SampleTemp<-relevel(dds.Ctrl.Comparison$SampleTemp, ref = "H_ctrl")

dds.HaCaT.HeatvsCtrl.Comparison$SampleTemp<-relevel(dds.HaCaT.HeatvsCtrl.Comparison$SampleTemp, ref = "H_ctrl")
dds.HaCaT.ColdvsCtrl.Comparison$SampleTemp<-relevel(dds.HaCaT.ColdvsCtrl.Comparison$SampleTemp, ref = "H_ctrl")
dds.TIGK.HeatvsCtrl.Comparison$SampleTemp<-relevel(dds.TIGK.HeatvsCtrl.Comparison$SampleTemp, ref = "T_ctrl")
dds.TIGK.ColdvsCtrl.Comparison$SampleTemp<-relevel(dds.TIGK.ColdvsCtrl.Comparison$SampleTemp, ref = "T_ctrl")

#####Running DESeq and annotating gene names
#running DESeq Heat
dds.Heat.Comparison<- DESeq(dds.Heat.Comparison)
res.Heat.comparison<-results(dds.Heat.Comparison)

write.csv(res.Heat.comparison, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/HaCaTvsTIGK.60C.DEGs.csv")
Heat.comparison<-read.csv("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/HaCaTvsTIGK.60C.DEGs.csv")

Heat.comparison<-Heat.comparison %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(grch38, by = c("X" = "ensgene"))
Heat.comparison <- Heat.comparison[!duplicated(Heat.comparison$symbol),] 
na.omit(Heat.comparison)
head(Heat.comparison)

#running DESeq Cold
dds.Cold.Comparison<- DESeq(dds.Cold.Comparison)
res.Cold.comparison<-results(dds.Cold.Comparison)

write.csv(res.Cold.comparison, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/HaCaTvsTIGK.-25C.DEGs.csv")
Cold.comparison<-read.csv("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/HaCaTvsTIGK.-25C.DEGs.csv")

Cold.comparison<-Cold.comparison %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(grch38, by = c("X" = "ensgene"))
Cold.comparison <- Cold.comparison[!duplicated(Cold.comparison$symbol),] 
na.omit(Cold.comparison)
head(Cold.comparison)

#running DESeq Ctrl
dds.Ctrl.Comparison<- DESeq(dds.Ctrl.Comparison)
res.Ctrl.Comparison<-results(dds.Ctrl.Comparison)

write.csv(res.Ctrl.Comparison, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/HaCaTvsTIGK.37C.DEGs.csv")
Ctrl.comparison<-read.csv("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/HaCaTvsTIGK.37C.DEGs.csv")

Ctrl.comparison<-Ctrl.comparison %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(grch38, by = c("X" = "ensgene"))
Ctrl.comparison <- Ctrl.comparison[!duplicated(Ctrl.comparison$symbol),] 
na.omit(Ctrl.comparison)
head(Ctrl.comparison)

#running DESeq HaCaT: Heat vs Ctrl
dds.HaCaT.HeatvsCtrl.Comparison<- DESeq(dds.HaCaT.HeatvsCtrl.Comparison)
res.HaCaT.HeatvsCtrl.Comparison<-results(dds.HaCaT.HeatvsCtrl.Comparison)

write.csv(res.HaCaT.HeatvsCtrl.Comparison, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/HaCaT.HeatvsCtrl.DEGs.csv")
HaCaT.HeatvsCtrl<-read.csv("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/HaCaT.HeatvsCtrl.DEGs.csv")

HaCaT.HeatvsCtrl<-HaCaT.HeatvsCtrl %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(grch38, by = c("X" = "ensgene"))
HaCaT.HeatvsCtrl <- HaCaT.HeatvsCtrl[!duplicated(HaCaT.HeatvsCtrl$symbol),] 
na.omit(HaCaT.HeatvsCtrl)
head(HaCaT.HeatvsCtrl)

#running DESeq HaCaT: Cold vs Ctrl
dds.HaCaT.ColdvsCtrl.Comparison<- DESeq(dds.HaCaT.ColdvsCtrl.Comparison)
res.HaCaT.ColdvsCtrl.Comparison<-results(dds.HaCaT.ColdvsCtrl.Comparison)

write.csv(res.HaCaT.ColdvsCtrl.Comparison, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/HaCaT.ColdvsCtrl.DEGs.csv")
HaCaT.ColdvsCtrl<-read.csv("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/HaCaT.ColdvsCtrl.DEGs.csv")

HaCaT.ColdvsCtrl<-HaCaT.ColdvsCtrl %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(grch38, by = c("X" = "ensgene"))
HaCaT.ColdvsCtrl <- HaCaT.ColdvsCtrl[!duplicated(HaCaT.ColdvsCtrl$symbol),] 
na.omit(HaCaT.ColdvsCtrl)
head(HaCaT.ColdvsCtrl)

#running DESeq TIGK: Heat vs Ctrl
dds.TIGK.HeatvsCtrl.Comparison<- DESeq(dds.TIGK.HeatvsCtrl.Comparison)
res.TIGK.HeatvsCtrl.Comparison<-results(dds.TIGK.HeatvsCtrl.Comparison)

write.csv(res.TIGK.HeatvsCtrl.Comparison, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/TIGK.HeatvsCtrl.DEGs.csv")
TIGK.HeatvsCtrl<-read.csv("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/TIGK.HeatvsCtrl.DEGs.csv")

TIGK.HeatvsCtrl<-TIGK.HeatvsCtrl %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(grch38, by = c("X" = "ensgene"))
TIGK.HeatvsCtrl <- TIGK.HeatvsCtrl[!duplicated(TIGK.HeatvsCtrl$symbol),] 
na.omit(TIGK.HeatvsCtrl)
head(TIGK.HeatvsCtrl)

#running DESeq TIGK: Cold vs Ctrl
dds.TIGK.ColdvsCtrl.Comparison<- DESeq(dds.TIGK.ColdvsCtrl.Comparison)
res.TIGK.ColdvsCtrl.Comparison<-results(dds.TIGK.ColdvsCtrl.Comparison)

write.csv(res.TIGK.ColdvsCtrl.Comparison, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/TIGK.ColdvsCtrl.DEGs.csv")
TIGK.ColdvsCtrl<-read.csv("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/TIGK.ColdvsCtrl.DEGs.csv")

TIGK.ColdvsCtrl<-TIGK.ColdvsCtrl %>% 
  dplyr::arrange(padj) %>% 
  dplyr::inner_join(grch38, by = c("X" = "ensgene"))
TIGK.ColdvsCtrl <- TIGK.ColdvsCtrl[!duplicated(TIGK.ColdvsCtrl$symbol),] 
na.omit(TIGK.ColdvsCtrl)
head(TIGK.ColdvsCtrl)

#### Acquiring list of DEGs using p.adj < 0.01 as cutoff
Heat.comparison.C1<-Heat.comparison%>%filter(padj <  0.01)
Cold.comparison.C1<-Cold.comparison%>%filter(padj <  0.01)
Ctrl.comparison.C1<-Ctrl.comparison%>%filter(padj <  0.01)
HaCaT.HeatvsCtrl.C1<-HaCaT.HeatvsCtrl%>%filter(padj <  0.01)
HaCaT.ColdvsCtrl.C1<-HaCaT.ColdvsCtrl%>%filter(padj <  0.01)
TIGK.HeatvsCtrl.C1<-TIGK.HeatvsCtrl%>%filter(padj <  0.01)
TIGK.ColdvsCtrl.C1<-TIGK.ColdvsCtrl%>%filter(padj <  0.01)

nrow(TIGK.ColdvsCtrl)

write.csv(Heat.comparison.C1, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/HaCaTvsTIGK.60C.DEGs.csv")
write.csv(Cold.comparison.C1, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/HaCaTvsTIGK.-25C.DEGs.csv")
write.csv(Ctrl.comparison.C1, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/HaCaTvsTIGK.37C.DEGs.csv")
write.csv(HaCaT.HeatvsCtrl.C1, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/HaCaT.HeatvsCtrl.DEGs.csv")
write.csv(HaCaT.ColdvsCtrl.C1, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/HaCaT.ColdvsCtrl.DEGs.csv")
write.csv(TIGK.HeatvsCtrl.C1, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/TIGK.HeatvsCtrl.DEGs.csv")
write.csv(TIGK.ColdvsCtrl.C1, file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Differential Gene Expression Analyses/TIGK.ColdvsCtrl.DEGs.csv")
#### Volcanoe Plots for all DEGs ####
#Heat Comparisons
Heat.comparison$diffexpressed <- "No differential expression"
Heat.comparison$diffexpressed[Heat.comparison$log2FoldChange > 0 & Heat.comparison$padj < 0.01] <- "Upregulated in TIGK"
Heat.comparison$diffexpressed[Heat.comparison$log2FoldChange < 0 & Heat.comparison$padj < 0.01] <- "Downregulated in HaCaT"

ggplot(data=Heat.comparison, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) + 
  theme(axis.title.x=element_text(size = 30), axis.title.y=element_text(size = 30), axis.text.x=element_text(size = 30), axis.text.y=element_text(size = 30), 
        legend.text = element_text(size=25), legend.title = element_blank(), panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = "grey", linewidth = 0.25),panel.grid.minor = element_line(color = "grey", linewidth = 0.25) )+
  geom_point(aes(colour=diffexpressed, fill=diffexpressed), size = 7, alpha = 1, shape = 21,colour = "black") + 
  geom_vline(xintercept=c(0), col="red", size = 2) +   geom_hline(yintercept=-log10(0.01), col="red", size = 2) 

#Cold Comparisons
Cold.comparison$diffexpressed <- "No differential expression"
Cold.comparison$diffexpressed[Cold.comparison$log2FoldChange > 0 & Cold.comparison$padj < 0.01] <- "Upregulated in TIGK"
Cold.comparison$diffexpressed[Cold.comparison$log2FoldChange < 0 & Cold.comparison$padj < 0.01] <- "Downregulated in HaCaT"

ggplot(data=Cold.comparison, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) + 
  theme(axis.title.x=element_text(size = 30), axis.title.y=element_text(size = 30), axis.text.x=element_text(size = 30), axis.text.y=element_text(size = 30), 
        legend.text = element_text(size=25), legend.title = element_blank(), panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = "grey", linewidth = 0.25),panel.grid.minor = element_line(color = "grey", linewidth = 0.25) )+
  geom_point(aes(colour=diffexpressed, fill=diffexpressed), size = 7, alpha = 1, shape = 21,colour = "black") + 
  geom_vline(xintercept=c(0), col="red", size = 2) +   geom_hline(yintercept=-log10(0.01), col="red", size = 2) 

#Control Comparisons
Ctrl.comparison$diffexpressed <- "No differential expression"
Ctrl.comparison$diffexpressed[Ctrl.comparison$log2FoldChange > 0 & Ctrl.comparison$padj < 0.01] <- "Upregulated in TIGK"
Ctrl.comparison$diffexpressed[Ctrl.comparison$log2FoldChange < 0 & Ctrl.comparison$padj < 0.01] <- "Downregulated in HaCaT"

ggplot(data=Ctrl.comparison, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) + 
  theme(axis.title.x=element_text(size = 30), axis.title.y=element_text(size = 30), axis.text.x=element_text(size = 30), axis.text.y=element_text(size = 30), 
        legend.text = element_text(size=25), legend.title = element_blank(), panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = "grey", linewidth = 0.25),panel.grid.minor = element_line(color = "grey", linewidth = 0.25) )+
  geom_point(aes(colour=diffexpressed, fill=diffexpressed), size = 7, alpha = 1, shape = 21,colour = "black") + 
  geom_vline(xintercept=c(0), col="red", size = 2) +   geom_hline(yintercept=-log10(0.01), col="red", size = 2) 


#HaCaT: 60C vs 37C Comparisons
HaCaT.HeatvsCtrl$diffexpressed <- "No differential expression"
HaCaT.HeatvsCtrl$diffexpressed[HaCaT.HeatvsCtrl$log2FoldChange > 0 & HaCaT.HeatvsCtrl$padj < 0.01] <- "Upregulated in 60C"
HaCaT.HeatvsCtrl$diffexpressed[HaCaT.HeatvsCtrl$log2FoldChange < 0 & HaCaT.HeatvsCtrl$padj < 0.01] <- "Downregulated in 60C"

ggplot(data=HaCaT.HeatvsCtrl, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) + 
  theme(axis.title.x=element_text(size = 30), axis.title.y=element_text(size = 30), axis.text.x=element_text(size = 30), axis.text.y=element_text(size = 30), 
        legend.text = element_text(size=25), legend.title = element_blank(), panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = "grey", linewidth = 0.25),panel.grid.minor = element_line(color = "grey", linewidth = 0.25) )+
  geom_point(aes(colour=diffexpressed, fill=diffexpressed), size = 7, alpha = 1, shape = 21,colour = "black") + 
  geom_vline(xintercept=c(0), col="red", size = 2) +   geom_hline(yintercept=-log10(0.01), col="red", size = 2) 



#Cold Comparisons
HaCaT.ColdvsCtrl$diffexpressed <- "No differential expression"
HaCaT.ColdvsCtrl$diffexpressed[HaCaT.ColdvsCtrl$log2FoldChange > 0 & HaCaT.ColdvsCtrl$padj < 0.01] <- "Upregulated in -25C"
HaCaT.ColdvsCtrl$diffexpressed[HaCaT.ColdvsCtrl$log2FoldChange < 0 & HaCaT.ColdvsCtrl$padj < 0.01] <- "Downregulated in -25C"

ggplot(data=HaCaT.ColdvsCtrl, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) + 
  theme(axis.title.x=element_text(size = 30), axis.title.y=element_text(size = 30), axis.text.x=element_text(size = 30), axis.text.y=element_text(size = 30), 
        legend.text = element_text(size=25), legend.title = element_blank(), panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = "grey", linewidth = 0.25),panel.grid.minor = element_line(color = "grey", linewidth = 0.25) )+
  geom_point(aes(colour=diffexpressed, fill=diffexpressed), size = 7, alpha = 1, shape = 21,colour = "black") + 
  geom_vline(xintercept=c(0), col="red", size = 2) +   geom_hline(yintercept=-log10(0.01), col="red", size = 2) 

#TIGK: 60C vs 37C Comparisons
TIGK.HeatvsCtrl$diffexpressed <- "No differential expression"
TIGK.HeatvsCtrl$diffexpressed[TIGK.HeatvsCtrl$log2FoldChange > 0 & TIGK.HeatvsCtrl$padj < 0.01] <- "Upregulated in 60C"
TIGK.HeatvsCtrl$diffexpressed[TIGK.HeatvsCtrl$log2FoldChange < 0 & TIGK.HeatvsCtrl$padj < 0.01] <- "Downregulated in 60C"

ggplot(data=TIGK.HeatvsCtrl, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) + 
  theme(axis.title.x=element_text(size = 30), axis.title.y=element_text(size = 30), axis.text.x=element_text(size = 30), axis.text.y=element_text(size = 30), 
        legend.text = element_text(size=25), legend.title = element_blank(), panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = "grey", linewidth = 0.25),panel.grid.minor = element_line(color = "grey", linewidth = 0.25) )+
  geom_point(aes(colour=diffexpressed, fill=diffexpressed), size = 7, alpha = 1, shape = 21,colour = "black") + 
  geom_vline(xintercept=c(0), col="red", size = 2) +   geom_hline(yintercept=-log10(0.01), col="red", size = 2) 


#Cold Comparisons
TIGK.ColdvsCtrl$diffexpressed <- "No differential expression"
TIGK.ColdvsCtrl$diffexpressed[TIGK.ColdvsCtrl$log2FoldChange > 0 & TIGK.ColdvsCtrl$padj < 0.01] <- "Upregulated in -25C"
TIGK.ColdvsCtrl$diffexpressed[TIGK.ColdvsCtrl$log2FoldChange < 0 & TIGK.ColdvsCtrl$padj < 0.01] <- "Downregulated in -25C"

ggplot(data=TIGK.ColdvsCtrl, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) + 
  theme(axis.title.x=element_text(size = 30), axis.title.y=element_text(size = 30), axis.text.x=element_text(size = 30), axis.text.y=element_text(size = 30), 
        legend.text = element_text(size=25), legend.title = element_blank(), panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(color = "grey", linewidth = 0.25),panel.grid.minor = element_line(color = "grey", linewidth = 0.25) )+
  geom_point(aes(colour=diffexpressed, fill=diffexpressed), size = 7, alpha = 1, shape = 21,colour = "black") + 
  geom_vline(xintercept=c(0), col="red", size = 2) +   geom_hline(yintercept=-log10(0.01), col="red", size = 2) 


##### EnrichR Gene Enrichment Analysis using DEGs from comparisons####
Heat.comparison.C1.upTIGK<-Heat.comparison.C1%>%filter(padj <  0.01 & log2FoldChange > 0)
Heat.comparison.C1.upHaCaT<-Heat.comparison.C1%>%filter(padj <  0.01 & log2FoldChange < 0)

Cold.comparison.C1.upTIGK<-Cold.comparison.C1%>%filter(padj <  0.01 & log2FoldChange > 0)
Cold.comparison.C1.upHaCaT<-Cold.comparison.C1%>%filter(padj <  0.01 & log2FoldChange < 0)

Ctrl.comparison.C1.upTIGK<-Ctrl.comparison.C1%>%filter(padj <  0.01 & log2FoldChange > 0)
Ctrl.comparison.C1.upHaCaT<-Ctrl.comparison.C1%>%filter(padj <  0.01 & log2FoldChange < 0)

HaCaT.HeatvsCtrl.C1.upregHeat<-HaCaT.HeatvsCtrl%>%filter(padj <  0.01 & log2FoldChange > 0)
HaCaT.HeatvsCtrl.C1.downregHeat<-HaCaT.HeatvsCtrl%>%filter(padj <  0.01 & log2FoldChange < 0)

HaCaT.ColdvsCtrl.C1.upregCold<-HaCaT.ColdvsCtrl%>%filter(padj <  0.01 & log2FoldChange > 0)
HaCaT.ColdvsCtrl.C1.downregCold<-HaCaT.ColdvsCtrl%>%filter(padj <  0.01 & log2FoldChange < 0)

TIGK.HeatvsCtrl.C1.upregHeat<-TIGK.HeatvsCtrl%>%filter(padj <  0.01 & log2FoldChange > 0)
TIGK.HeatvsCtrl.C1.downregHeat<-TIGK.HeatvsCtrl%>%filter(padj <  0.01 & log2FoldChange < 0)

TIGK.ColdvsCtrl.C1.upregCold<-TIGK.ColdvsCtrl%>%filter(padj <  0.01 & log2FoldChange > 0)
TIGK.ColdvsCtrl.C1.downregCold<-TIGK.ColdvsCtrl%>%filter(padj <  0.01 & log2FoldChange < 0)

#Heat comparisons
enriched.HvsT.Heat.HaCaT <- enrichr(Heat.comparison.C1.upHaCaT$symbol, dbs)
enriched.HvsT.Heat.TIGK <- enrichr(Heat.comparison.C1.upTIGK$symbol, dbs)
#Cold comparisons
enriched.HvsT.Cold.HaCaT <- enrichr(Cold.comparison.C1.upHaCaT$symbol, dbs)
enriched.HvsT.Cold.TIGK <- enrichr(Cold.comparison.C1.upTIGK$symbol, dbs)
#Ctrl comparisons
enriched.HvsT.Ctrl.HaCaT <- enrichr(Ctrl.comparison.C1.upHaCaT$symbol, dbs)
enriched.HvsT.Ctrl.TIGK <- enrichr(Ctrl.comparison.C1.upTIGK$symbol, dbs)

#HaCaT:60C vs 37C comparisons
enriched.H60vsH37.upH60 <- enrichr(HaCaT.HeatvsCtrl.C1.upregHeat$symbol, dbs)
enriched.H60vsH37.downH60 <- enrichr(HaCaT.HeatvsCtrl.C1.downregHeat$symbol, dbs)

#HaCaT:-25C vs 37C comparisons
enriched.Hneg25vsH37.upHneg25 <- enrichr(HaCaT.ColdvsCtrl.C1.upregCold$symbol, dbs)
enriched.Hneg25vsH37.downHneg25 <- enrichr(HaCaT.ColdvsCtrl.C1.downregCold$symbol, dbs)

#TIGK:60C vs 37C comparisons
enriched.T60vsT37.upT60 <- enrichr(TIGK.HeatvsCtrl.C1.upregHeat$symbol, dbs)
enriched.T60vsT37.downT60 <- enrichr(TIGK.HeatvsCtrl.C1.downregHeat$symbol, dbs)

#TIGK:-25C vs 37C comparisons
enriched.Tneg25vsT37.upTneg25 <- enrichr(TIGK.ColdvsCtrl.C1.upregCold$symbol, dbs)
enriched.Tneg25vsT37.downTneg25 <- enrichr(TIGK.ColdvsCtrl.C1.downregCold$symbol, dbs)

# Heat comparisons
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up HaCaT 60C")
addWorksheet(Enrichments, "MF Up HaCaT 60C")
addWorksheet(Enrichments, "Reactome Up HaCaT 60C")
addWorksheet(Enrichments, "BP Up TIGK 60C")
addWorksheet(Enrichments, "MF Up TIGK 60C")
addWorksheet(Enrichments, "Reactome Up TIGK 60C")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up HaCaT 60C", x = as.data.frame(enriched.HvsT.Heat.HaCaT$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up HaCaT 60C", x = as.data.frame(enriched.HvsT.Heat.HaCaT$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up HaCaT 60C", x = as.data.frame(enriched.HvsT.Heat.HaCaT$Reactome_2022))
writeData(Enrichments, sheet = "BP Up TIGK 60C",  x = as.data.frame(enriched.HvsT.Heat.TIGK$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up TIGK 60C",  x = as.data.frame(enriched.HvsT.Heat.TIGK$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up TIGK 60C", x = as.data.frame(enriched.HvsT.Heat.TIGK$Reactome_2022))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of all DEGs/H60vsT60.xlsx")

# Cold comparisons
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up HaCaT -25C")
addWorksheet(Enrichments, "MF Up HaCaT -25C")
addWorksheet(Enrichments, "Reactome Up HaCaT -25C")
addWorksheet(Enrichments, "BP Up TIGK -25C")
addWorksheet(Enrichments, "MF Up TIGK -25C")
addWorksheet(Enrichments, "Reactome Up TIGK -25C")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up HaCaT -25C", x = as.data.frame(enriched.HvsT.Cold.HaCaT$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up HaCaT -25C", x = as.data.frame(enriched.HvsT.Cold.HaCaT$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up HaCaT -25C", x = as.data.frame(enriched.HvsT.Cold.HaCaT$Reactome_2022))
writeData(Enrichments, sheet = "BP Up TIGK -25C",  x = as.data.frame(enriched.HvsT.Cold.TIGK$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up TIGK -25C",  x = as.data.frame(enriched.HvsT.Cold.TIGK$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up TIGK -25C", x = as.data.frame(enriched.HvsT.Cold.TIGK$Reactome_2022))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of all DEGs/H-25vsT-25.xlsx")

# Control comparisons
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up HaCaT 37C")
addWorksheet(Enrichments, "MF Up HaCaT 37C")
addWorksheet(Enrichments, "Reactome Up HaCaT 37C")
addWorksheet(Enrichments, "BP Up TIGK 37C")
addWorksheet(Enrichments, "MF Up TIGK 37C")
addWorksheet(Enrichments, "Reactome Up TIGK 37C")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up HaCaT 37C", x = as.data.frame(enriched.HvsT.Ctrl.HaCaT$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up HaCaT 37C", x = as.data.frame(enriched.HvsT.Ctrl.HaCaT$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up HaCaT 37C", x = as.data.frame(enriched.HvsT.Ctrl.HaCaT$Reactome_2022))
writeData(Enrichments, sheet = "BP Up TIGK 37C",  x = as.data.frame(enriched.HvsT.Ctrl.TIGK$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up TIGK 37C",  x = as.data.frame(enriched.HvsT.Ctrl.TIGK$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up TIGK 37C", x = as.data.frame(enriched.HvsT.Ctrl.TIGK$Reactome_2022))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of all DEGs/H37vsT37.xlsx")

#HaCaT:60C vs 37C comparisons
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up HaCaT 60C")
addWorksheet(Enrichments, "MF Up HaCaT 60C")
addWorksheet(Enrichments, "Reactome Up HaCaT 60C")
addWorksheet(Enrichments, "BP Down HaCaT 60C")
addWorksheet(Enrichments, "MF Down HaCaT 60C")
addWorksheet(Enrichments, "Reactome Down HaCaT 60C")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up HaCaT 60C", x = as.data.frame(enriched.H60vsH37.upH60$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up HaCaT 60C", x = as.data.frame(enriched.H60vsH37.upH60$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up HaCaT 60C", x = as.data.frame(enriched.H60vsH37.upH60$Reactome_2022))
writeData(Enrichments, sheet = "BP Down HaCaT 60C",  x = as.data.frame(enriched.H60vsH37.downH60$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Down HaCaT 60C",  x = as.data.frame(enriched.H60vsH37.downH60$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Down HaCaT 60C", x = as.data.frame(enriched.H60vsH37.downH60$Reactome_2022))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of all DEGs/H60vsH37.xlsx")

#HaCaT:-25C vs 37C comparisons
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up HaCaT -25C")
addWorksheet(Enrichments, "MF Up HaCaT -25C")
addWorksheet(Enrichments, "Reactome Up HaCaT -25C")
addWorksheet(Enrichments, "BP Down HaCaT -25C")
addWorksheet(Enrichments, "MF Down HaCaT -25C")
addWorksheet(Enrichments, "Reactome Down HaCaT -25C")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up HaCaT -25C", x = as.data.frame(enriched.Hneg25vsH37.upHneg25$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up HaCaT -25C", x = as.data.frame(enriched.Hneg25vsH37.upHneg25$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up HaCaT -25C", x = as.data.frame(enriched.Hneg25vsH37.upHneg25$Reactome_2022))
writeData(Enrichments, sheet = "BP Down HaCaT -25C",  x = as.data.frame(enriched.Hneg25vsH37.downHneg25$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Down HaCaT -25C",  x = as.data.frame(enriched.Hneg25vsH37.downHneg25$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Down HaCaT -25C", x = as.data.frame(enriched.Hneg25vsH37.downHneg25$Reactome_2022))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of all DEGs/H-25vsH37.xlsx")

#TIGK:60C vs 37C comparisons
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up TIGK 60C")
addWorksheet(Enrichments, "MF Up TIGK 60C")
addWorksheet(Enrichments, "Reactome Up TIGK 60C")
addWorksheet(Enrichments, "BP Down TIGK 60C")
addWorksheet(Enrichments, "MF Down TIGK 60C")
addWorksheet(Enrichments, "Reactome Down TIGK 60C")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up TIGK 60C", x = as.data.frame(enriched.H60vsH37.upH60$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up TIGK 60C", x = as.data.frame(enriched.H60vsH37.upH60$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up TIGK 60C", x = as.data.frame(enriched.H60vsH37.upH60$Reactome_2022))
writeData(Enrichments, sheet = "BP Down TIGK 60C",  x = as.data.frame(enriched.H60vsH37.downH60$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Down TIGK 60C",  x = as.data.frame(enriched.H60vsH37.downH60$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Down TIGK 60C", x = as.data.frame(enriched.H60vsH37.downH60$Reactome_2022))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of all DEGs/T60vsT37.xlsx")

#TIGK:-25C vs 37C comparisons
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up TIGK -25C")
addWorksheet(Enrichments, "MF Up TIGK -25C")
addWorksheet(Enrichments, "Reactome Up TIGK -25C")
addWorksheet(Enrichments, "BP Down TIGK -25C")
addWorksheet(Enrichments, "MF Down TIGK -25C")
addWorksheet(Enrichments, "Reactome Down TIGK -25C")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up TIGK -25C", x = as.data.frame(enriched.Hneg25vsH37.upHneg25$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up TIGK -25C", x = as.data.frame(enriched.Hneg25vsH37.upHneg25$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up TIGK -25C", x = as.data.frame(enriched.Hneg25vsH37.upHneg25$Reactome_2022))
writeData(Enrichments, sheet = "BP Down TIGK -25C",  x = as.data.frame(enriched.Hneg25vsH37.downHneg25$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Down TIGK -25C",  x = as.data.frame(enriched.Hneg25vsH37.downHneg25$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Down TIGK -25C", x = as.data.frame(enriched.Hneg25vsH37.downHneg25$Reactome_2022))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of all DEGs/T-25vsT37.xlsx")


#### Extracting out the HSPs and making Heatmaps HSPs ####
#loading in dataframe for HGNC list of approved gene names for heat shock proteins (HSPs)
Chaperonins<-read.csv(file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Heat Shock Gene Name_HGNC/Chaperonins.csv")
DNAJ<-read.csv(file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Heat Shock Gene Name_HGNC/DNAJ.csv")
HSP70kda<-read.csv(file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Heat Shock Gene Name_HGNC/Heatshock70kDA.csv")
HSP90kda<-read.csv(file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Heat Shock Gene Name_HGNC/Heatshock90kDA.csv")
SmallHSP<-read.csv(file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Heat Shock Gene Name_HGNC/Small heat shock proteins.csv")


#Pulling out all HGNC list of approved HSP gene names - heat injury response HaCaT
HaCaT.HeatvsCtrl.C1.chaperonins<-  subset(HaCaT.HeatvsCtrl.C1, symbol %in% Chaperonins$Approved.symbol)
HaCaT.HeatvsCtrl.C1.DNAJ<-  subset(HaCaT.HeatvsCtrl.C1, symbol %in% DNAJ$Approved.symbol)
HaCaT.HeatvsCtrl.C1.HSP70kDA<-  subset(HaCaT.HeatvsCtrl.C1, symbol %in% HSP70kda$Approved.symbol)
HaCaT.HeatvsCtrl.C1.HSP90kDA<-  subset(HaCaT.HeatvsCtrl.C1, symbol %in% HSP90kda$Approved.symbol)
HaCaT.HeatvsCtrl.C1.SmallHSP<-  subset(HaCaT.HeatvsCtrl.C1, symbol %in% SmallHSP$Approved.symbol)

HaCaT.HeatResponse.HSPs<-rbind(HaCaT.HeatvsCtrl.C1.chaperonins,HaCaT.HeatvsCtrl.C1.DNAJ, HaCaT.HeatvsCtrl.C1.HSP70kDA,
                               HaCaT.HeatvsCtrl.C1.HSP90kDA, HaCaT.HeatvsCtrl.C1.SmallHSP)

mat.HaCaT.HeatResponse.HSPs<-counts(dds.HaCaT.HeatvsCtrl.Comparison, normalized = T)[HaCaT.HeatResponse.HSPs$X,] #this allows you to get counts data for only the Ensembl IDs correlating to Diff. Expressed Collagens
mat.HaCaT.HeatResponse.HSPs.z<-t(apply(mat.HaCaT.HeatResponse.HSPs, 1, scale))
colnames(mat.HaCaT.HeatResponse.HSPs.z)<-colnames(mat.HaCaT.HeatResponse.HSPs)
pheatmap(mat.HaCaT.HeatResponse.HSPs.z, labels_row = HaCaT.HeatResponse.HSPs$symbol,
         fontsize_row = 10,
         fontsize_col = 10)

#Pulling out all HGNC list of approved HSP gene names - heat injury response TIGK
TIGK.HeatvsCtrl.C1.chaperonins<-  subset(TIGK.HeatvsCtrl.C1, symbol %in% Chaperonins$Approved.symbol)
TIGK.HeatvsCtrl.C1.DNAJ<-  subset(TIGK.HeatvsCtrl.C1, symbol %in% DNAJ$Approved.symbol)
TIGK.HeatvsCtrl.C1.HSP70kDA<-  subset(TIGK.HeatvsCtrl.C1, symbol %in% HSP70kda$Approved.symbol)
TIGK.HeatvsCtrl.C1.HSP90kDA<-  subset(TIGK.HeatvsCtrl.C1, symbol %in% HSP90kda$Approved.symbol)
TIGK.HeatvsCtrl.C1.SmallHSP<-  subset(TIGK.HeatvsCtrl.C1, symbol %in% SmallHSP$Approved.symbol)

TIGK.HeatResponse.HSPs<-rbind(TIGK.HeatvsCtrl.C1.chaperonins,TIGK.HeatvsCtrl.C1.DNAJ, TIGK.HeatvsCtrl.C1.HSP70kDA,
                              TIGK.HeatvsCtrl.C1.HSP90kDA, TIGK.HeatvsCtrl.C1.SmallHSP)

mat.TIGK.HeatResponse.HSPs<-counts(dds.TIGK.HeatvsCtrl.Comparison, normalized = T)[TIGK.HeatResponse.HSPs$X,] #this allows you to get counts data for only the Ensembl IDs correlating to Diff. Expressed Collagens
mat.TIGK.HeatResponse.HSPs.z<-t(apply(mat.TIGK.HeatResponse.HSPs, 1, scale))
colnames(mat.TIGK.HeatResponse.HSPs.z)<-colnames(mat.TIGK.HeatResponse.HSPs)
pheatmap(mat.TIGK.HeatResponse.HSPs.z, labels_row = TIGK.HeatResponse.HSPs$symbol,
         fontsize_row = 10,
         fontsize_col = 10)

#Pulling out all HGNC list of approved HSP gene names - Cold injury response HaCaT
HaCaT.ColdvsCtrl.C1.chaperonins<-  subset(HaCaT.ColdvsCtrl.C1, symbol %in% Chaperonins$Approved.symbol)
HaCaT.ColdvsCtrl.C1.DNAJ<-  subset(HaCaT.ColdvsCtrl.C1, symbol %in% DNAJ$Approved.symbol)
HaCaT.ColdvsCtrl.C1.HSP70kDA<-  subset(HaCaT.ColdvsCtrl.C1, symbol %in% HSP70kda$Approved.symbol)
HaCaT.ColdvsCtrl.C1.HSP90kDA<-  subset(HaCaT.ColdvsCtrl.C1, symbol %in% HSP90kda$Approved.symbol)
HaCaT.ColdvsCtrl.C1.SmallHSP<-  subset(HaCaT.ColdvsCtrl.C1, symbol %in% SmallHSP$Approved.symbol)

HaCaT.ColdResponse.HSPs<-rbind(HaCaT.ColdvsCtrl.C1.chaperonins,HaCaT.ColdvsCtrl.C1.DNAJ, HaCaT.ColdvsCtrl.C1.HSP70kDA,
                               HaCaT.ColdvsCtrl.C1.HSP90kDA, HaCaT.ColdvsCtrl.C1.SmallHSP)

mat.HaCaT.ColdResponse.HSPs<-counts(dds.HaCaT.ColdvsCtrl.Comparison, normalized = T)[HaCaT.ColdResponse.HSPs$X,] #this allows you to get counts data for only the Ensembl IDs correlating to Diff. Expressed Collagens
mat.HaCaT.ColdResponse.HSPs.z<-t(apply(mat.HaCaT.ColdResponse.HSPs, 1, scale))
colnames(mat.HaCaT.ColdResponse.HSPs.z)<-colnames(mat.HaCaT.ColdResponse.HSPs)
pheatmap(mat.HaCaT.ColdResponse.HSPs.z, labels_row = HaCaT.ColdResponse.HSPs$symbol,
         fontsize_row = 10,
         fontsize_col = 10)

#Pulling out all HGNC list of approved HSP gene names - Cold injury response TIGK
TIGK.ColdvsCtrl.C1.chaperonins<-  subset(TIGK.ColdvsCtrl.C1, symbol %in% Chaperonins$Approved.symbol)
TIGK.ColdvsCtrl.C1.DNAJ<-  subset(TIGK.ColdvsCtrl.C1, symbol %in% DNAJ$Approved.symbol)
TIGK.ColdvsCtrl.C1.HSP70kDA<-  subset(TIGK.ColdvsCtrl.C1, symbol %in% HSP70kda$Approved.symbol)
TIGK.ColdvsCtrl.C1.HSP90kDA<-  subset(TIGK.ColdvsCtrl.C1, symbol %in% HSP90kda$Approved.symbol)
TIGK.ColdvsCtrl.C1.SmallHSP<-  subset(TIGK.ColdvsCtrl.C1, symbol %in% SmallHSP$Approved.symbol)

TIGK.ColdResponse.HSPs<-rbind(TIGK.ColdvsCtrl.C1.chaperonins,TIGK.ColdvsCtrl.C1.DNAJ, TIGK.ColdvsCtrl.C1.HSP70kDA,
                              TIGK.ColdvsCtrl.C1.HSP90kDA, TIGK.ColdvsCtrl.C1.SmallHSP)

mat.TIGK.ColdResponse.HSPs<-counts(dds.TIGK.ColdvsCtrl.Comparison, normalized = T)[TIGK.ColdResponse.HSPs$X,] #this allows you to get counts data for only the Ensembl IDs correlating to Diff. Expressed Collagens
mat.TIGK.ColdResponse.HSPs.z<-t(apply(mat.TIGK.ColdResponse.HSPs, 1, scale))
colnames(mat.TIGK.ColdResponse.HSPs.z)<-colnames(mat.TIGK.ColdResponse.HSPs)
pheatmap(mat.TIGK.ColdResponse.HSPs.z, labels_row = TIGK.ColdResponse.HSPs$symbol,
         fontsize_row = 10,
         fontsize_col = 10)


#Pulling out all HGNC list of approved HSP gene names - 60C (heat) comparison
HaCaTvsTIGK.Heat.chaperonins<-  subset(Heat.comparison.C1, symbol %in% Chaperonins$Approved.symbol)
HaCaTvsTIGK.Heat.DNAJ<-  subset(Heat.comparison.C1, symbol %in% DNAJ$Approved.symbol)
HaCaTvsTIGK.Heat.HSP70kDA<-  subset(Heat.comparison.C1, symbol %in% HSP70kda$Approved.symbol)
HaCaTvsTIGK.Heat.HSP90kDA<-  subset(Heat.comparison.C1, symbol %in% HSP90kda$Approved.symbol)
HaCaTvsTIGK.Heat.SmallHSP<-  subset(Heat.comparison.C1, symbol %in% SmallHSP$Approved.symbol)

Heat.comparison.HSPs<-rbind(HaCaTvsTIGK.Heat.chaperonins,HaCaTvsTIGK.Heat.DNAJ, HaCaTvsTIGK.Heat.HSP70kDA,
                            HaCaTvsTIGK.Heat.HSP90kDA, HaCaTvsTIGK.Heat.SmallHSP)

mat.Heat.comparison.HSPs<-counts(dds.Heat.Comparison, normalized = T)[Heat.comparison.HSPs$X,] #this allows you to get counts data for only the Ensembl IDs correlating to Diff. Expressed Collagens
mat.Heat.comparison.HSPs.z<-t(apply(mat.Heat.comparison.HSPs, 1, scale))
colnames(mat.Heat.comparison.HSPs.z)<-colnames(mat.Heat.comparison.HSPs)
pheatmap(mat.Heat.comparison.HSPs.z, labels_row = Heat.comparison.HSPs$symbol, cluster_cols = T, 
         fontsize_row = 10,
         fontsize_col = 10)

#Pulling out all HGNC list of approved HSP gene names - -25C (cold) comparison
HaCaTvsTIGK.Cold.chaperonins<-  subset(Cold.comparison.C1, symbol %in% Chaperonins$Approved.symbol)
HaCaTvsTIGK.Cold.DNAJ<-  subset(Cold.comparison.C1, symbol %in% DNAJ$Approved.symbol)
HaCaTvsTIGK.Cold.HSP70kDA<-  subset(Cold.comparison.C1, symbol %in% HSP70kda$Approved.symbol)
HaCaTvsTIGK.Cold.HSP90kDA<-  subset(Cold.comparison.C1, symbol %in% HSP90kda$Approved.symbol)
HaCaTvsTIGK.Cold.SmallHSP<-  subset(Cold.comparison.C1, symbol %in% SmallHSP$Approved.symbol)

Cold.comparison.HSPs<-rbind(HaCaTvsTIGK.Cold.chaperonins,HaCaTvsTIGK.Cold.DNAJ, HaCaTvsTIGK.Cold.HSP70kDA,
                            HaCaTvsTIGK.Cold.HSP90kDA, HaCaTvsTIGK.Cold.SmallHSP)

mat.Cold.comparison.HSPs<-counts(dds.Cold.Comparison, normalized = T)[Cold.comparison.HSPs$X,] #this allows you to get counts data for only the Ensembl IDs correlating to Diff. Expressed Collagens
mat.Cold.comparison.HSPs.z<-t(apply(mat.Cold.comparison.HSPs, 1, scale))
colnames(mat.Cold.comparison.HSPs.z)<-colnames(mat.Cold.comparison.HSPs)
pheatmap(mat.Cold.comparison.HSPs.z, labels_row = Cold.comparison.HSPs$symbol, cluster_cols = T, 
         fontsize_row = 10,
         fontsize_col = 10)

#Pulling out all HGNC list of approved HSP gene names - baseline (37C) comparison
HaCaTvsTIGK.Ctrl.chaperonins<-  subset(Ctrl.comparison.C1, symbol %in% Chaperonins$Approved.symbol)
HaCaTvsTIGK.Ctrl.DNAJ<-  subset(Ctrl.comparison.C1, symbol %in% DNAJ$Approved.symbol)
HaCaTvsTIGK.Ctrl.HSP70kDA<-  subset(Ctrl.comparison.C1, symbol %in% HSP70kda$Approved.symbol)
HaCaTvsTIGK.Ctrl.HSP90kDA<-  subset(Ctrl.comparison.C1, symbol %in% HSP90kda$Approved.symbol)
HaCaTvsTIGK.Ctrl.SmallHSP<-  subset(Ctrl.comparison.C1, symbol %in% SmallHSP$Approved.symbol)

Ctrl.comparison.HSPs<-rbind(HaCaTvsTIGK.Ctrl.chaperonins,HaCaTvsTIGK.Ctrl.DNAJ, HaCaTvsTIGK.Ctrl.HSP70kDA,
                            HaCaTvsTIGK.Ctrl.HSP90kDA, HaCaTvsTIGK.Ctrl.SmallHSP)

mat.Ctrl.comparison.HSPs<-counts(dds.Ctrl.Comparison, normalized = T)[Ctrl.comparison.HSPs$X,] #this allows you to get counts data for only the Ensembl IDs correlating to Diff. Expressed Collagens
mat.Ctrl.comparison.HSPs.z<-t(apply(mat.Ctrl.comparison.HSPs, 1, scale))
colnames(mat.Ctrl.comparison.HSPs.z)<-colnames(mat.Ctrl.comparison.HSPs)
pheatmap(mat.Ctrl.comparison.HSPs.z, labels_row = Ctrl.comparison.HSPs$symbol, cluster_cols = T, 
         fontsize_row = 10,
         fontsize_col = 10)

write.csv(Heat.comparison.HSPs,  
          file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/DE of HSPs/HvsT.60C.HSPs.csv")
write.csv(Cold.comparison.HSPs,  
          file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/DE of HSPs/HvsT.-25C.HSPs.csv")
write.csv(Ctrl.comparison.HSPs,  
          file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/DE of HSPs/HvsT.37C.HSPs.csv")
write.csv(HaCaT.HeatResponse.HSPs,  
          file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/DE of HSPs/H60vsH37.HSPs.csv")
write.csv(HaCaT.ColdResponse.HSPs,  
          file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/DE of HSPs/H-25vsH37.HSPs.csv")
write.csv(TIGK.HeatResponse.HSPs,  
          file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/DE of HSPs/T60vsT37.HSPs.csv")
write.csv(TIGK.ColdResponse.HSPs,  
          file = "C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/DE of HSPs/T-25vsT37.HSPs.csv")

#### Enrichment of upregulated HSPs ####
Heat.comparison.HSPs.HaCaT<-Heat.comparison.HSPs%>%filter(log2FoldChange       <  0)
Heat.comparison.HSPs.TIGK<-Heat.comparison.HSPs%>%filter(log2FoldChange       >  0)
Cold.comparison.HSPs.HaCaT<-Cold.comparison.HSPs%>%filter(log2FoldChange       <  0)
Cold.comparison.HSPs.TIGK<-Cold.comparison.HSPs%>%filter(log2FoldChange       >  0)
Ctrl.comparison.HSPs.TIGK<-Ctrl.comparison.HSPs%>%filter(log2FoldChange       >  0)
Ctrl.comparison.HSPs.HaCaT<-Ctrl.comparison.HSPs%>%filter(log2FoldChange       <  0)

HaCaT.HeatResponse.HSPs.upregulated <- HaCaT.HeatResponse.HSPs%>%filter(log2FoldChange       >  0)
HaCaT.ColdResponse.HSPs.upregulated <- HaCaT.ColdResponse.HSPs%>%filter(log2FoldChange       >  0)
TIGK.HeatResponse.HSPs.upregulated  <- TIGK.HeatResponse.HSPs%>%filter(log2FoldChange       >  0)
TIGK.ColdResponse.HSPs.upregulated  <- TIGK.ColdResponse.HSPs%>%filter(log2FoldChange       >  0)

#Heat comparisons
enriched.Heat.comparison.HSPs.HaCaT <- enrichr(Heat.comparison.HSPs.HaCaT$symbol, dbs)
enriched.Heat.comparison.HSPs.TIGK <- enrichr(Heat.comparison.HSPs.TIGK$symbol, dbs)
#Cold comparisons
enriched.Cold.comparison.HSPs.HaCaT <- enrichr(Cold.comparison.HSPs.HaCaT$symbol, dbs)
enriched.Cold.comparison.HSPs.TIGK <- enrichr(Cold.comparison.HSPs.TIGK$symbol, dbs)
#Ctrl comparisons
enriched.Ctrl.comparison.HSPs.HaCaT <- enrichr(Ctrl.comparison.HSPs.HaCaT$symbol, dbs)
enriched.Ctrl.comparison.HSPs.TIGK <- enrichr(Ctrl.comparison.HSPs.TIGK$symbol, dbs)

#HaCaT:60C vs 37C comparisons
enriched.HaCaT.HeatResponse.HSPs.upregulated <- enrichr(HaCaT.HeatResponse.HSPs.upregulated$symbol, dbs)

#HaCaT:-25C vs 37C comparisons
enriched.HaCaT.ColdResponse.HSPs.upregulated <- enrichr(HaCaT.ColdResponse.HSPs.upregulated$symbol, dbs)

#TIGK:60C vs 37C comparisons
enriched.TIGK.HeatResponse.HSPs.upregulated <- enrichr(TIGK.HeatResponse.HSPs.upregulated$symbol, dbs)

#TIGK:-25C vs 37C comparisons
enriched.TIGK.ColdResponse.HSPs.upregulated <- enrichr(TIGK.ColdResponse.HSPs.upregulated$symbol, dbs)

# Heat comparisons
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up HaCaT 60C")
addWorksheet(Enrichments, "MF Up HaCaT 60C")
addWorksheet(Enrichments, "Reactome Up HaCaT 60C")
addWorksheet(Enrichments, "BP Up TIGK 60C")
addWorksheet(Enrichments, "MF Up TIGK 60C")
addWorksheet(Enrichments, "Reactome Up TIGK 60C")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up HaCaT 60C", x = as.data.frame(enriched.Heat.comparison.HSPs.HaCaT$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up HaCaT 60C", x = as.data.frame(enriched.Heat.comparison.HSPs.HaCaT$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up HaCaT 60C", x = as.data.frame(enriched.Heat.comparison.HSPs.HaCaT$Reactome_2022))
writeData(Enrichments, sheet = "BP Up TIGK 60C",  x = as.data.frame(enriched.Heat.comparison.HSPs.TIGK$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up TIGK 60C",  x = as.data.frame(enriched.Heat.comparison.HSPs.TIGK$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up TIGK 60C", x = as.data.frame(enriched.Heat.comparison.HSPs.TIGK$Reactome_2022))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs/H60vsT60.xlsx")

# Cold comparisons
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up HaCaT -25")
addWorksheet(Enrichments, "MF Up HaCaT -25")
addWorksheet(Enrichments, "Reactome Up HaCaT -25")
addWorksheet(Enrichments, "BP Up TIGK -25")
addWorksheet(Enrichments, "MF Up TIGK -25")
addWorksheet(Enrichments, "Reactome Up TIGK -25")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up HaCaT -25", x = as.data.frame(enriched.Cold.comparison.HSPs.HaCaT$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up HaCaT -25", x = as.data.frame(enriched.Cold.comparison.HSPs.HaCaT$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up HaCaT -25", x = as.data.frame(enriched.Cold.comparison.HSPs.HaCaT$Reactome_2022))
writeData(Enrichments, sheet = "BP Up TIGK -25",  x = as.data.frame(enriched.Cold.comparison.HSPs.TIGK$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up TIGK -25",  x = as.data.frame(enriched.Cold.comparison.HSPs.TIGK$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up TIGK -25", x = as.data.frame(enriched.Cold.comparison.HSPs.TIGK$Reactome_2022))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs/H-25vsT-25.xlsx")

# Ctrl comparisons
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up HaCaT 37C")
addWorksheet(Enrichments, "MF Up HaCaT 37C")
addWorksheet(Enrichments, "Reactome Up HaCaT 37C")
addWorksheet(Enrichments, "BP Up TIGK 37C")
addWorksheet(Enrichments, "MF Up TIGK 37C")
addWorksheet(Enrichments, "Reactome Up TIGK 37C")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up HaCaT 37C", x = as.data.frame(enriched.Ctrl.comparison.HSPs.HaCaT$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up HaCaT 37C", x = as.data.frame(enriched.Ctrl.comparison.HSPs.HaCaT$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up HaCaT 37C", x = as.data.frame(enriched.Ctrl.comparison.HSPs.HaCaT$Reactome_2022))
writeData(Enrichments, sheet = "BP Up TIGK 37C",  x = as.data.frame(enriched.Ctrl.comparison.HSPs.TIGK$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up TIGK 37C",  x = as.data.frame(enriched.Ctrl.comparison.HSPs.TIGK$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up TIGK 37C", x = as.data.frame(enriched.Ctrl.comparison.HSPs.TIGK$Reactome_2022))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs/H37CvsT37C.xlsx")

# Heat Injury comparisons
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up HaCaT 60C")
addWorksheet(Enrichments, "MF Up HaCaT 60C")
addWorksheet(Enrichments, "Reactome Up HaCaT 60C")
addWorksheet(Enrichments, "BP Up TIGK 60C")
addWorksheet(Enrichments, "MF Up TIGK 60C")
addWorksheet(Enrichments, "Reactome Up TIGK 60C")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up HaCaT 60C", x = as.data.frame(enriched.HaCaT.HeatResponse.HSPs.upregulated$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up HaCaT 60C", x = as.data.frame(enriched.HaCaT.HeatResponse.HSPs.upregulated$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up HaCaT 60C", x = as.data.frame(enriched.HaCaT.HeatResponse.HSPs.upregulated$Reactome_2022))
writeData(Enrichments, sheet = "BP Up TIGK 60C",  x = as.data.frame(enriched.TIGK.HeatResponse.HSPs.upregulated$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up TIGK 60C",  x = as.data.frame(enriched.TIGK.HeatResponse.HSPs.upregulated$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up TIGK 60C", x = as.data.frame(enriched.TIGK.HeatResponse.HSPs.upregulated$Reactome_2022))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs/Thermal Injury Response.HSP.xlsx")

# Cold Injury comparisons
Enrichments <- createWorkbook()
# Add some sheets to the workbook
addWorksheet(Enrichments, "BP Up HaCaT -25C")
addWorksheet(Enrichments, "MF Up HaCaT -25C")
addWorksheet(Enrichments, "Reactome Up HaCaT -25C")
addWorksheet(Enrichments, "BP Up TIGK -25C")
addWorksheet(Enrichments, "MF Up TIGK -25C")
addWorksheet(Enrichments, "Reactome Up TIGK -25C")
# Write the data to the sheets
writeData(Enrichments, sheet = "BP Up HaCaT -25C", x = as.data.frame(enriched.HaCaT.ColdResponse.HSPs.upregulated$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up HaCaT -25C", x = as.data.frame(enriched.HaCaT.ColdResponse.HSPs.upregulated$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up HaCaT -25C", x = as.data.frame(enriched.HaCaT.ColdResponse.HSPs.upregulated$Reactome_2022))
writeData(Enrichments, sheet = "BP Up TIGK -25C",  x = as.data.frame(enriched.TIGK.ColdResponse.HSPs.upregulated$GO_Biological_Process_2023))
writeData(Enrichments, sheet = "MF Up TIGK -25C",  x = as.data.frame(enriched.TIGK.ColdResponse.HSPs.upregulated$GO_Molecular_Function_2023))
writeData(Enrichments, sheet = "Reactome Up TIGK -25C", x = as.data.frame(enriched.TIGK.ColdResponse.HSPs.upregulated$Reactome_2022))
# Export the file
saveWorkbook(Enrichments, file="C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs/Cold Injury Response.HSP.xlsx")

########################Data visualization for Enrichment for DEGs: Bubble plots ####
#HvsT_37C Control Comparisons ####
enriched.HvsT.Ctrl.HaCaT <- enrichr(Ctrl.comparison.C1.upHaCaT$symbol, dbs)

Ctrl.HaCaT.DEGs.BP<-as.data.frame(enriched.HvsT.Ctrl.HaCaT$GO_Biological_Process_2023)
Ctrl.HaCaT.DEGs.MF<-as.data.frame(enriched.HvsT.Ctrl.HaCaT$GO_Molecular_Function_2023)
Ctrl.HaCaT.DEGs.Reactome<-as.data.frame(enriched.HvsT.Ctrl.HaCaT$Reactome_2022)


Ctrl.HaCaT.DEGs.BP.Top10<-Ctrl.HaCaT.DEGs.BP[1:10,]
Ctrl.HaCaT.DEGs.BP.Top10
frac <- Ctrl.HaCaT.DEGs.BP.Top10$Overlap
Ctrl.HaCaT.DEGs.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(Ctrl.HaCaT.DEGs.BP.Top10, aes(x=Ctrl.HaCaT.DEGs.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Ctrl.HaCaT.DEGs.BP.Top10$Adjusted.P.value, size=(Ctrl.HaCaT.DEGs.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/37C_HaCaT_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

Ctrl.HaCaT.DEGs.MF.Top10<-Ctrl.HaCaT.DEGs.MF[1:10,]
Ctrl.HaCaT.DEGs.MF.Top10
frac <- Ctrl.HaCaT.DEGs.MF.Top10$Overlap
Ctrl.HaCaT.DEGs.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(Ctrl.HaCaT.DEGs.MF.Top10, aes(x=Ctrl.HaCaT.DEGs.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=Ctrl.HaCaT.DEGs.MF.Top10$`Gene Ratio`, colour=Ctrl.HaCaT.DEGs.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/37C_HaCaT_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

Ctrl.HaCaT.DEGs.Reactome.Top10<-Ctrl.HaCaT.DEGs.Reactome[1:10,]
Ctrl.HaCaT.DEGs.Reactome.Top10
frac <- Ctrl.HaCaT.DEGs.Reactome.Top10$Overlap
Ctrl.HaCaT.DEGs.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(Ctrl.HaCaT.DEGs.Reactome.Top10, aes(x=Ctrl.HaCaT.DEGs.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Ctrl.HaCaT.DEGs.Reactome.Top10$Adjusted.P.value, size=(Ctrl.HaCaT.DEGs.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/37C_HaCaT_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/37C_HaCaT_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, labels = "AUTO", nrow = 1)
dev.off()



enriched.HvsT.Ctrl.TIGK <- enrichr(Ctrl.comparison.C1.upTIGK$symbol, dbs)
Ctrl.TIGK.DEGs.BP<-as.data.frame(enriched.HvsT.Ctrl.TIGK$GO_Biological_Process_2023)
Ctrl.TIGK.DEGs.MF<-as.data.frame(enriched.HvsT.Ctrl.TIGK$GO_Molecular_Function_2023)
Ctrl.TIGK.DEGs.Reactome<-as.data.frame(enriched.HvsT.Ctrl.TIGK$Reactome_2022)

Ctrl.TIGK.DEGs.BP.Top10<-Ctrl.TIGK.DEGs.BP[1:10,]
Ctrl.TIGK.DEGs.BP.Top10
frac <- Ctrl.TIGK.DEGs.BP.Top10$Overlap
Ctrl.TIGK.DEGs.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(Ctrl.TIGK.DEGs.BP.Top10, aes(x=Ctrl.TIGK.DEGs.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Ctrl.TIGK.DEGs.BP.Top10$Adjusted.P.value, size=(Ctrl.TIGK.DEGs.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/37C_TIGK_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

Ctrl.TIGK.DEGs.MF.Top10<-Ctrl.TIGK.DEGs.MF[1:10,]
Ctrl.TIGK.DEGs.MF.Top10
frac <- Ctrl.TIGK.DEGs.MF.Top10$Overlap
Ctrl.TIGK.DEGs.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p5<-ggplot(Ctrl.TIGK.DEGs.MF.Top10, aes(x=Ctrl.TIGK.DEGs.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Ctrl.TIGK.DEGs.MF.Top10$Adjusted.P.value, size=(Ctrl.TIGK.DEGs.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/37C_TIGK_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p5
dev.off()

Ctrl.TIGK.DEGs.Reactome.Top10<-Ctrl.TIGK.DEGs.Reactome[1:10,]
Ctrl.TIGK.DEGs.Reactome.Top10
frac <- Ctrl.TIGK.DEGs.Reactome.Top10$Overlap
Ctrl.TIGK.DEGs.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p6<-ggplot(Ctrl.TIGK.DEGs.Reactome.Top10, aes(x=Ctrl.TIGK.DEGs.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Ctrl.TIGK.DEGs.Reactome.Top10$Adjusted.P.value, size=(Ctrl.TIGK.DEGs.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/37C_TIGK_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p6
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/37C_TIGK_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p4,p5,p6, labels = "AUTO", nrow = 1)
dev.off()

#HvsT_60C comparisons ####
enriched.HvsT.Heat.HaCaT <- enrichr(Heat.comparison.C1.upHaCaT$symbol, dbs)

Heat.HaCaT.DEGs.BP<-as.data.frame(enriched.HvsT.Heat.HaCaT$GO_Biological_Process_2023)
Heat.HaCaT.DEGs.MF<-as.data.frame(enriched.HvsT.Heat.HaCaT$GO_Molecular_Function_2023)
Heat.HaCaT.DEGs.Reactome<-as.data.frame(enriched.HvsT.Heat.HaCaT$Reactome_2022)


Heat.HaCaT.DEGs.BP.Top10<-Heat.HaCaT.DEGs.BP[1:10,]
Heat.HaCaT.DEGs.BP.Top10
frac <- Heat.HaCaT.DEGs.BP.Top10$Overlap
Heat.HaCaT.DEGs.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(Heat.HaCaT.DEGs.BP.Top10, aes(x=Heat.HaCaT.DEGs.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Heat.HaCaT.DEGs.BP.Top10$Adjusted.P.value, size=(Heat.HaCaT.DEGs.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60CHaCaT_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

Heat.HaCaT.DEGs.MF.Top10<-Heat.HaCaT.DEGs.MF[1:10,]
Heat.HaCaT.DEGs.MF.Top10
frac <- Heat.HaCaT.DEGs.MF.Top10$Overlap
Heat.HaCaT.DEGs.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(Heat.HaCaT.DEGs.MF.Top10, aes(x=Heat.HaCaT.DEGs.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=Heat.HaCaT.DEGs.MF.Top10$`Gene Ratio`, colour=Heat.HaCaT.DEGs.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60CHaCaT_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

Heat.HaCaT.DEGs.Reactome.Top10<-Heat.HaCaT.DEGs.Reactome[1:10,]
Heat.HaCaT.DEGs.Reactome.Top10
frac <- Heat.HaCaT.DEGs.Reactome.Top10$Overlap
Heat.HaCaT.DEGs.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(Heat.HaCaT.DEGs.Reactome.Top10, aes(x=Heat.HaCaT.DEGs.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Heat.HaCaT.DEGs.Reactome.Top10$Adjusted.P.value, size=(Heat.HaCaT.DEGs.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60CHaCaT_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60CHaCaT_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, labels = "AUTO", nrow = 1)
dev.off()



enriched.HvsT.Heat.TIGK <- enrichr(Heat.comparison.C1.upTIGK$symbol, dbs)
Heat.TIGK.DEGs.BP<-as.data.frame(enriched.HvsT.Heat.TIGK$GO_Biological_Process_2023)
Heat.TIGK.DEGs.MF<-as.data.frame(enriched.HvsT.Heat.TIGK$GO_Molecular_Function_2023)
Heat.TIGK.DEGs.Reactome<-as.data.frame(enriched.HvsT.Heat.TIGK$Reactome_2022)

Heat.TIGK.DEGs.BP.Top10<-Heat.TIGK.DEGs.BP[1:10,]
Heat.TIGK.DEGs.BP.Top10
frac <- Heat.TIGK.DEGs.BP.Top10$Overlap
Heat.TIGK.DEGs.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(Heat.TIGK.DEGs.BP.Top10, aes(x=Heat.TIGK.DEGs.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Heat.TIGK.DEGs.BP.Top10$Adjusted.P.value, size=(Heat.TIGK.DEGs.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60CTIGK_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

Heat.TIGK.DEGs.Reactome.Top10<-Heat.TIGK.DEGs.Reactome[1:10,]
Heat.TIGK.DEGs.Reactome.Top10
frac <- Heat.TIGK.DEGs.Reactome.Top10$Overlap
Heat.TIGK.DEGs.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))


Heat.TIGK.DEGs.MF.Top10<-Heat.TIGK.DEGs.MF[1:10,]
Heat.TIGK.DEGs.MF.Top10
frac <- Heat.TIGK.DEGs.MF.Top10$Overlap
Heat.TIGK.DEGs.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p5<-ggplot(Heat.TIGK.DEGs.MF.Top10, aes(x=Heat.TIGK.DEGs.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Heat.TIGK.DEGs.MF.Top10$Adjusted.P.value, size=(Heat.TIGK.DEGs.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60CTIGK_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p5
dev.off()

p6<-ggplot(Heat.TIGK.DEGs.Reactome.Top10, aes(x=Heat.TIGK.DEGs.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Heat.TIGK.DEGs.Reactome.Top10$Adjusted.P.value, size=(Heat.TIGK.DEGs.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60CTIGK_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p6
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60CTIGK_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p4,p5,p6, labels = "AUTO", nrow = 1)
dev.off()


#HvsT_-25C comparisons ####
enriched.HvsT.Cold.HaCaT <- enrichr(Cold.comparison.C1.upHaCaT$symbol, dbs)

Cold.HaCaT.DEGs.BP<-as.data.frame(enriched.HvsT.Cold.HaCaT$GO_Biological_Process_2023)
Cold.HaCaT.DEGs.MF<-as.data.frame(enriched.HvsT.Cold.HaCaT$GO_Molecular_Function_2023)
Cold.HaCaT.DEGs.Reactome<-as.data.frame(enriched.HvsT.Cold.HaCaT$Reactome_2022)


Cold.HaCaT.DEGs.BP.Top10<-Cold.HaCaT.DEGs.BP[1:10,]
Cold.HaCaT.DEGs.BP.Top10
frac <- Cold.HaCaT.DEGs.BP.Top10$Overlap
Cold.HaCaT.DEGs.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(Cold.HaCaT.DEGs.BP.Top10, aes(x=Cold.HaCaT.DEGs.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Cold.HaCaT.DEGs.BP.Top10$Adjusted.P.value, size=(Cold.HaCaT.DEGs.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25CHaCaT_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

Cold.HaCaT.DEGs.MF.Top10<-Cold.HaCaT.DEGs.MF[1:10,]
Cold.HaCaT.DEGs.MF.Top10
frac <- Cold.HaCaT.DEGs.MF.Top10$Overlap
Cold.HaCaT.DEGs.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(Cold.HaCaT.DEGs.MF.Top10, aes(x=Cold.HaCaT.DEGs.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=Cold.HaCaT.DEGs.MF.Top10$`Gene Ratio`, colour=Cold.HaCaT.DEGs.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25CHaCaT_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

Cold.HaCaT.DEGs.Reactome.Top10<-Cold.HaCaT.DEGs.Reactome[1:10,]
Cold.HaCaT.DEGs.Reactome.Top10
frac <- Cold.HaCaT.DEGs.Reactome.Top10$Overlap
Cold.HaCaT.DEGs.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(Cold.HaCaT.DEGs.Reactome.Top10, aes(x=Cold.HaCaT.DEGs.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Cold.HaCaT.DEGs.Reactome.Top10$Adjusted.P.value, size=(Cold.HaCaT.DEGs.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25CHaCaT_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25CHaCaT_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, labels = "AUTO", nrow = 1)
dev.off()



enriched.HvsT.Cold.TIGK <- enrichr(Cold.comparison.C1.upTIGK$symbol, dbs)
Cold.TIGK.DEGs.BP<-as.data.frame(enriched.HvsT.Cold.TIGK$GO_Biological_Process_2023)
Cold.TIGK.DEGs.MF<-as.data.frame(enriched.HvsT.Cold.TIGK$GO_Molecular_Function_2023)
Cold.TIGK.DEGs.Reactome<-as.data.frame(enriched.HvsT.Cold.TIGK$Reactome_2022)

Cold.TIGK.DEGs.BP.Top10<-Cold.TIGK.DEGs.BP[1:10,]
Cold.TIGK.DEGs.BP.Top10
frac <- Cold.TIGK.DEGs.BP.Top10$Overlap
Cold.TIGK.DEGs.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(Cold.TIGK.DEGs.BP.Top10, aes(x=Cold.TIGK.DEGs.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Cold.TIGK.DEGs.BP.Top10$Adjusted.P.value, size=(Cold.TIGK.DEGs.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25CTIGK_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

Cold.TIGK.DEGs.MF.Top10<-Cold.TIGK.DEGs.MF[1:10,]
Cold.TIGK.DEGs.MF.Top10
frac <- Cold.TIGK.DEGs.MF.Top10$Overlap
Cold.TIGK.DEGs.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p5<-ggplot(Cold.TIGK.DEGs.MF.Top10, aes(x=Cold.TIGK.DEGs.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Cold.TIGK.DEGs.MF.Top10$Adjusted.P.value, size=(Cold.TIGK.DEGs.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25CTIGK_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p5
dev.off()

Cold.TIGK.DEGs.Reactome.Top10<-Cold.TIGK.DEGs.Reactome[1:10,]
Cold.TIGK.DEGs.Reactome.Top10
frac <- Cold.TIGK.DEGs.Reactome.Top10$Overlap
Cold.TIGK.DEGs.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p6<-ggplot(Cold.TIGK.DEGs.Reactome.Top10, aes(x=Cold.TIGK.DEGs.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Cold.TIGK.DEGs.Reactome.Top10$Adjusted.P.value, size=(Cold.TIGK.DEGs.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25CTIGK_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p6
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25CTIGK_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p4,p5,p6, labels = "AUTO", nrow = 1)
dev.off()




#HaCaT:60C vs 37C comparisons ####
enriched.H60vsH37.upH60 <- enrichr(HaCaT.HeatvsCtrl.C1.upregHeat$symbol, dbs)
HeatResponse.UpDEGs.HaCaT.BP<-as.data.frame(enriched.H60vsH37.upH60$GO_Biological_Process_2023)
HeatResponse.UpDEGs.HaCaT.MF<-as.data.frame(enriched.H60vsH37.upH60$GO_Molecular_Function_2023)
HeatResponse.UpDEGs.HaCaT.Reactome<-as.data.frame(enriched.H60vsH37.upH60$Reactome_2022)


HeatResponse.UpDEGs.HaCaT.BP.Top10<-HeatResponse.UpDEGs.HaCaT.BP[1:10,]
HeatResponse.UpDEGs.HaCaT.BP.Top10
frac <- HeatResponse.UpDEGs.HaCaT.BP.Top10$Overlap
HeatResponse.UpDEGs.HaCaT.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(HeatResponse.UpDEGs.HaCaT.BP.Top10, aes(x=HeatResponse.UpDEGs.HaCaT.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=HeatResponse.UpDEGs.HaCaT.BP.Top10$Adjusted.P.value, size=(HeatResponse.UpDEGs.HaCaT.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60Cvs37C_UpDEGs_HaCaT_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

HeatResponse.UpDEGs.HaCaT.MF.Top10<-HeatResponse.UpDEGs.HaCaT.MF[1:10,]
HeatResponse.UpDEGs.HaCaT.MF.Top10
frac <- HeatResponse.UpDEGs.HaCaT.MF.Top10$Overlap
HeatResponse.UpDEGs.HaCaT.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(HeatResponse.UpDEGs.HaCaT.MF.Top10, aes(x=HeatResponse.UpDEGs.HaCaT.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=HeatResponse.UpDEGs.HaCaT.MF.Top10$`Gene Ratio`, colour=HeatResponse.UpDEGs.HaCaT.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60Cvs37C_UpDEGs_HaCaT_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

HeatResponse.UpDEGs.HaCaT.Reactome.Top10<-HeatResponse.UpDEGs.HaCaT.Reactome[1:10,]
HeatResponse.UpDEGs.HaCaT.Reactome.Top10
frac <- HeatResponse.UpDEGs.HaCaT.Reactome.Top10$Overlap
HeatResponse.UpDEGs.HaCaT.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(HeatResponse.UpDEGs.HaCaT.Reactome.Top10, aes(x=HeatResponse.UpDEGs.HaCaT.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=HeatResponse.UpDEGs.HaCaT.Reactome.Top10$Adjusted.P.value, size=(HeatResponse.UpDEGs.HaCaT.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60Cvs37C_UpDEGs_HaCaT_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60Cvs37C_UpDEGs_HaCaT_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, labels = "AUTO", nrow = 1)
dev.off()


enriched.H60vsH37.downH60 <- enrichr(HaCaT.HeatvsCtrl.C1.downregHeat$symbol, dbs)
HeatResponse.DownDEGs.HaCaT.BP<-as.data.frame(enriched.H60vsH37.downH60$GO_Biological_Process_2023)
HeatResponse.DownDEGs.HaCaT.MF<-as.data.frame(enriched.H60vsH37.downH60$GO_Molecular_Function_2023)
HeatResponse.DownDEGs.HaCaT.Reactome<-as.data.frame(enriched.H60vsH37.downH60$Reactome_2022)

HeatResponse.DownDEGs.HaCaT.BP.Top10<-HeatResponse.DownDEGs.HaCaT.BP[1:10,]
HeatResponse.DownDEGs.HaCaT.BP.Top10
frac <- HeatResponse.DownDEGs.HaCaT.BP.Top10$Overlap
HeatResponse.DownDEGs.HaCaT.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(HeatResponse.DownDEGs.HaCaT.BP.Top10, aes(x=HeatResponse.DownDEGs.HaCaT.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=HeatResponse.DownDEGs.HaCaT.BP.Top10$Adjusted.P.value, size=(HeatResponse.DownDEGs.HaCaT.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60Cvs37C_DownDEGs_HaCaT_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

HeatResponse.DownDEGs.HaCaT.MF.Top10<-HeatResponse.DownDEGs.HaCaT.MF[1:10,]
HeatResponse.DownDEGs.HaCaT.MF.Top10
frac <- HeatResponse.DownDEGs.HaCaT.MF.Top10$Overlap
HeatResponse.DownDEGs.HaCaT.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p5<-ggplot(HeatResponse.DownDEGs.HaCaT.MF.Top10, aes(x=HeatResponse.DownDEGs.HaCaT.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=HeatResponse.DownDEGs.HaCaT.MF.Top10$`Gene Ratio`, colour=HeatResponse.DownDEGs.HaCaT.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60Cvs37C_DownDEGs_HaCaT_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p5
dev.off()

HeatResponse.DownDEGs.HaCaT.Reactome.Top10<-HeatResponse.DownDEGs.HaCaT.Reactome[1:10,]
HeatResponse.DownDEGs.HaCaT.Reactome.Top10
frac <- HeatResponse.DownDEGs.HaCaT.Reactome.Top10$Overlap
HeatResponse.DownDEGs.HaCaT.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p6<-ggplot(HeatResponse.DownDEGs.HaCaT.Reactome.Top10, aes(x=HeatResponse.DownDEGs.HaCaT.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=HeatResponse.DownDEGs.HaCaT.Reactome.Top10$Adjusted.P.value, size=(HeatResponse.DownDEGs.HaCaT.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60Cvs37C_DownDEGs_HaCaT_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p6
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60Cvs37C_DownDEGs_HaCaT_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p4,p5,p6, labels = "AUTO", nrow = 1)
dev.off()

#HaCaT:-25C vs 37C comparisons ####
enriched.Hneg25vsH37.upHneg25 <- enrichr(HaCaT.ColdvsCtrl.C1.upregCold$symbol, dbs)
ColdResponse.UpDEGs.HaCaT.BP<-as.data.frame(enriched.Hneg25vsH37.upHneg25$GO_Biological_Process_2023)
ColdResponse.UpDEGs.HaCaT.MF<-as.data.frame(enriched.Hneg25vsH37.upHneg25$GO_Molecular_Function_2023)
ColdResponse.UpDEGs.HaCaT.Reactome<-as.data.frame(enriched.Hneg25vsH37.upHneg25$Reactome_2022)


ColdResponse.UpDEGs.HaCaT.BP.Top10<-ColdResponse.UpDEGs.HaCaT.BP[1:10,]
ColdResponse.UpDEGs.HaCaT.BP.Top10
frac <- ColdResponse.UpDEGs.HaCaT.BP.Top10$Overlap
ColdResponse.UpDEGs.HaCaT.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(ColdResponse.UpDEGs.HaCaT.BP.Top10, aes(x=ColdResponse.UpDEGs.HaCaT.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=ColdResponse.UpDEGs.HaCaT.BP.Top10$Adjusted.P.value, size=(ColdResponse.UpDEGs.HaCaT.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25Cvs37C_UpDEGs_HaCaT_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

ColdResponse.UpDEGs.HaCaT.MF.Top10<-ColdResponse.UpDEGs.HaCaT.MF[1:10,]
ColdResponse.UpDEGs.HaCaT.MF.Top10
frac <- ColdResponse.UpDEGs.HaCaT.MF.Top10$Overlap
ColdResponse.UpDEGs.HaCaT.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(ColdResponse.UpDEGs.HaCaT.MF.Top10, aes(x=ColdResponse.UpDEGs.HaCaT.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=ColdResponse.UpDEGs.HaCaT.MF.Top10$`Gene Ratio`, colour=ColdResponse.UpDEGs.HaCaT.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25Cvs37C_UpDEGs_HaCaT_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

ColdResponse.UpDEGs.HaCaT.Reactome.Top10<-ColdResponse.UpDEGs.HaCaT.Reactome[1:10,]
ColdResponse.UpDEGs.HaCaT.Reactome.Top10
frac <- ColdResponse.UpDEGs.HaCaT.Reactome.Top10$Overlap
ColdResponse.UpDEGs.HaCaT.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(ColdResponse.UpDEGs.HaCaT.Reactome.Top10, aes(x=ColdResponse.UpDEGs.HaCaT.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=ColdResponse.UpDEGs.HaCaT.Reactome.Top10$Adjusted.P.value, size=(ColdResponse.UpDEGs.HaCaT.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25Cvs37C_UpDEGs_HaCaT_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25Cvs37C_UpDEGs_HaCaT_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, labels = "AUTO", nrow = 1)
dev.off()


enriched.Hneg25vsH37.downHneg25 <- enrichr(HaCaT.ColdvsCtrl.C1.downregCold$symbol, dbs)
ColdResponse.DownDEGs.HaCaT.BP<-as.data.frame(enriched.Hneg25vsH37.downHneg25$GO_Biological_Process_2023)
ColdResponse.DownDEGs.HaCaT.MF<-as.data.frame(enriched.Hneg25vsH37.downHneg25$GO_Molecular_Function_2023)
ColdResponse.DownDEGs.HaCaT.Reactome<-as.data.frame(enriched.Hneg25vsH37.downHneg25$Reactome_2022)

ColdResponse.DownDEGs.HaCaT.BP.Top10<-ColdResponse.DownDEGs.HaCaT.BP[1:10,]
ColdResponse.DownDEGs.HaCaT.BP.Top10
frac <- ColdResponse.DownDEGs.HaCaT.BP.Top10$Overlap
ColdResponse.DownDEGs.HaCaT.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(ColdResponse.DownDEGs.HaCaT.BP.Top10, aes(x=ColdResponse.DownDEGs.HaCaT.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=ColdResponse.DownDEGs.HaCaT.BP.Top10$Adjusted.P.value, size=(ColdResponse.DownDEGs.HaCaT.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25Cvs37C_DownDEGs_HaCaT_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

ColdResponse.DownDEGs.HaCaT.MF.Top10<-ColdResponse.DownDEGs.HaCaT.MF[1:10,]
ColdResponse.DownDEGs.HaCaT.MF.Top10
frac <- ColdResponse.DownDEGs.HaCaT.MF.Top10$Overlap
ColdResponse.DownDEGs.HaCaT.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p5<-ggplot(ColdResponse.DownDEGs.HaCaT.MF.Top10, aes(x=ColdResponse.DownDEGs.HaCaT.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=ColdResponse.DownDEGs.HaCaT.MF.Top10$`Gene Ratio`, colour=ColdResponse.DownDEGs.HaCaT.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25Cvs37C_DownDEGs_HaCaT_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p5
dev.off()

ColdResponse.DownDEGs.HaCaT.Reactome.Top10<-ColdResponse.DownDEGs.HaCaT.Reactome[1:10,]
ColdResponse.DownDEGs.HaCaT.Reactome.Top10
frac <- ColdResponse.DownDEGs.HaCaT.Reactome.Top10$Overlap
ColdResponse.DownDEGs.HaCaT.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p6<-ggplot(ColdResponse.DownDEGs.HaCaT.Reactome.Top10, aes(x=ColdResponse.DownDEGs.HaCaT.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=ColdResponse.DownDEGs.HaCaT.Reactome.Top10$Adjusted.P.value, size=(ColdResponse.DownDEGs.HaCaT.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25Cvs37C_DownDEGs_HaCaT_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p6
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25Cvs37C_DownDEGs_HaCaT_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p4,p5,p6, labels = "AUTO", nrow = 1)
dev.off()

#TIGK:60C vs 37C comparisons ####
enriched.T60vsT37.upT60 <- enrichr(TIGK.HeatvsCtrl.C1.upregHeat$symbol, dbs)
HeatResponse.UpDEGs.TIGK.BP<-as.data.frame(enriched.T60vsT37.upT60$GO_Biological_Process_2023)
HeatResponse.UpDEGs.TIGK.MF<-as.data.frame(enriched.T60vsT37.upT60$GO_Molecular_Function_2023)
HeatResponse.UpDEGs.TIGK.Reactome<-as.data.frame(enriched.T60vsT37.upT60$Reactome_2022)


HeatResponse.UpDEGs.TIGK.BP.Top10<-HeatResponse.UpDEGs.TIGK.BP[1:10,]
HeatResponse.UpDEGs.TIGK.BP.Top10
frac <- HeatResponse.UpDEGs.TIGK.BP.Top10$Overlap
HeatResponse.UpDEGs.TIGK.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(HeatResponse.UpDEGs.TIGK.BP.Top10, aes(x=HeatResponse.UpDEGs.TIGK.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=HeatResponse.UpDEGs.TIGK.BP.Top10$Adjusted.P.value, size=(HeatResponse.UpDEGs.TIGK.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60Cvs37C_UpDEGs_TIGK_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

HeatResponse.UpDEGs.TIGK.MF.Top10<-HeatResponse.UpDEGs.TIGK.MF[1:10,]
HeatResponse.UpDEGs.TIGK.MF.Top10
frac <- HeatResponse.UpDEGs.TIGK.MF.Top10$Overlap
HeatResponse.UpDEGs.TIGK.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(HeatResponse.UpDEGs.TIGK.MF.Top10, aes(x=HeatResponse.UpDEGs.TIGK.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=HeatResponse.UpDEGs.TIGK.MF.Top10$`Gene Ratio`, colour=HeatResponse.UpDEGs.TIGK.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60Cvs37C_UpDEGs_TIGK_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

HeatResponse.UpDEGs.TIGK.Reactome.Top10<-HeatResponse.UpDEGs.TIGK.Reactome[1:10,]
HeatResponse.UpDEGs.TIGK.Reactome.Top10
frac <- HeatResponse.UpDEGs.TIGK.Reactome.Top10$Overlap
HeatResponse.UpDEGs.TIGK.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(HeatResponse.UpDEGs.TIGK.Reactome.Top10, aes(x=HeatResponse.UpDEGs.TIGK.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=HeatResponse.UpDEGs.TIGK.Reactome.Top10$Adjusted.P.value, size=(HeatResponse.UpDEGs.TIGK.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60Cvs37C_UpDEGs_TIGK_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60Cvs37C_UpDEGs_TIGK_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, labels = "AUTO", nrow = 1)
dev.off()


enriched.T60vsT37.downT60 <- enrichr(TIGK.HeatvsCtrl.C1.downregHeat$symbol, dbs)
HeatResponse.DownDEGs.TIGK.BP<-as.data.frame(enriched.T60vsT37.downT60$GO_Biological_Process_2023)
HeatResponse.DownDEGs.TIGK.MF<-as.data.frame(enriched.T60vsT37.downT60$GO_Molecular_Function_2023)
HeatResponse.DownDEGs.TIGK.Reactome<-as.data.frame(enriched.T60vsT37.downT60$Reactome_2022)

HeatResponse.DownDEGs.TIGK.BP.Top10<-HeatResponse.DownDEGs.TIGK.BP[1:10,]
HeatResponse.DownDEGs.TIGK.BP.Top10
frac <- HeatResponse.DownDEGs.TIGK.BP.Top10$Overlap
HeatResponse.DownDEGs.TIGK.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(HeatResponse.DownDEGs.TIGK.BP.Top10, aes(x=HeatResponse.DownDEGs.TIGK.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=HeatResponse.DownDEGs.TIGK.BP.Top10$Adjusted.P.value, size=(HeatResponse.DownDEGs.TIGK.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60Cvs37C_DownDEGs_TIGK_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

HeatResponse.DownDEGs.TIGK.MF.Top10<-HeatResponse.DownDEGs.TIGK.MF[1:10,]
HeatResponse.DownDEGs.TIGK.MF.Top10
frac <- HeatResponse.DownDEGs.TIGK.MF.Top10$Overlap
HeatResponse.DownDEGs.TIGK.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p5<-ggplot(HeatResponse.DownDEGs.TIGK.MF.Top10, aes(x=HeatResponse.DownDEGs.TIGK.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=HeatResponse.DownDEGs.TIGK.MF.Top10$`Gene Ratio`, colour=HeatResponse.DownDEGs.TIGK.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60Cvs37C_DownDEGs_TIGK_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p5
dev.off()

HeatResponse.DownDEGs.TIGK.Reactome.Top10<-HeatResponse.DownDEGs.TIGK.Reactome[1:10,]
HeatResponse.DownDEGs.TIGK.Reactome.Top10
frac <- HeatResponse.DownDEGs.TIGK.Reactome.Top10$Overlap
HeatResponse.DownDEGs.TIGK.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p6<-ggplot(HeatResponse.DownDEGs.TIGK.Reactome.Top10, aes(x=HeatResponse.DownDEGs.TIGK.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=HeatResponse.DownDEGs.TIGK.Reactome.Top10$Adjusted.P.value, size=(HeatResponse.DownDEGs.TIGK.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60Cvs37C_DownDEGs_TIGK_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p6
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/60Cvs37C_DownDEGs_TIGK_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p4,p5,p6, labels = "AUTO", nrow = 1)
dev.off()

#TIGK:-25C vs 37C comparisons ####
enriched.Tneg25vsT37.upTneg25 <- enrichr(TIGK.ColdvsCtrl.C1.upregCold$symbol, dbs)
ColdResponse.UpDEGs.TIGK.BP<-as.data.frame(enriched.Tneg25vsT37.upTneg25$GO_Biological_Process_2023)
ColdResponse.UpDEGs.TIGK.MF<-as.data.frame(enriched.Tneg25vsT37.upTneg25$GO_Molecular_Function_2023)
ColdResponse.UpDEGs.TIGK.Reactome<-as.data.frame(enriched.Tneg25vsT37.upTneg25$Reactome_2022)


ColdResponse.UpDEGs.TIGK.BP.Top10<-ColdResponse.UpDEGs.TIGK.BP[1:10,]
ColdResponse.UpDEGs.TIGK.BP.Top10
frac <- ColdResponse.UpDEGs.TIGK.BP.Top10$Overlap
ColdResponse.UpDEGs.TIGK.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(ColdResponse.UpDEGs.TIGK.BP.Top10, aes(x=ColdResponse.UpDEGs.TIGK.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=ColdResponse.UpDEGs.TIGK.BP.Top10$Adjusted.P.value, size=(ColdResponse.UpDEGs.TIGK.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25Cvs37C_UpDEGs_TIGK_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

ColdResponse.UpDEGs.TIGK.MF.Top10<-ColdResponse.UpDEGs.TIGK.MF[1:10,]
ColdResponse.UpDEGs.TIGK.MF.Top10
frac <- ColdResponse.UpDEGs.TIGK.MF.Top10$Overlap
ColdResponse.UpDEGs.TIGK.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(ColdResponse.UpDEGs.TIGK.MF.Top10, aes(x=ColdResponse.UpDEGs.TIGK.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=ColdResponse.UpDEGs.TIGK.MF.Top10$`Gene Ratio`, colour=ColdResponse.UpDEGs.TIGK.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25Cvs37C_UpDEGs_TIGK_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

ColdResponse.UpDEGs.TIGK.Reactome.Top10<-ColdResponse.UpDEGs.TIGK.Reactome[1:10,]
ColdResponse.UpDEGs.TIGK.Reactome.Top10
frac <- ColdResponse.UpDEGs.TIGK.Reactome.Top10$Overlap
ColdResponse.UpDEGs.TIGK.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(ColdResponse.UpDEGs.TIGK.Reactome.Top10, aes(x=ColdResponse.UpDEGs.TIGK.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=ColdResponse.UpDEGs.TIGK.Reactome.Top10$Adjusted.P.value, size=(ColdResponse.UpDEGs.TIGK.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25Cvs37C_UpDEGs_TIGK_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25Cvs37C_UpDEGs_TIGK_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, labels = "AUTO", nrow = 1)
dev.off()


enriched.Tneg25vsT37.downTneg25 <- enrichr(TIGK.ColdvsCtrl.C1.downregCold$symbol, dbs)
ColdResponse.DownDEGs.TIGK.BP<-as.data.frame(enriched.Tneg25vsT37.downTneg25$GO_Biological_Process_2023)
ColdResponse.DownDEGs.TIGK.MF<-as.data.frame(enriched.Tneg25vsT37.downTneg25$GO_Molecular_Function_2023)
ColdResponse.DownDEGs.TIGK.Reactome<-as.data.frame(enriched.Tneg25vsT37.downTneg25$Reactome_2022)

ColdResponse.DownDEGs.TIGK.BP.Top10<-ColdResponse.DownDEGs.TIGK.BP[1:10,]
ColdResponse.DownDEGs.TIGK.BP.Top10
frac <- ColdResponse.DownDEGs.TIGK.BP.Top10$Overlap
ColdResponse.DownDEGs.TIGK.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(ColdResponse.DownDEGs.TIGK.BP.Top10, aes(x=ColdResponse.DownDEGs.TIGK.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=ColdResponse.DownDEGs.TIGK.BP.Top10$Adjusted.P.value, size=(ColdResponse.DownDEGs.TIGK.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25Cvs37C_DownDEGs_TIGK_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

ColdResponse.DownDEGs.TIGK.MF.Top10<-ColdResponse.DownDEGs.TIGK.MF[1:10,]
ColdResponse.DownDEGs.TIGK.MF.Top10
frac <- ColdResponse.DownDEGs.TIGK.MF.Top10$Overlap
ColdResponse.DownDEGs.TIGK.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p5<-ggplot(ColdResponse.DownDEGs.TIGK.MF.Top10, aes(x=ColdResponse.DownDEGs.TIGK.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=ColdResponse.DownDEGs.TIGK.MF.Top10$`Gene Ratio`, colour=ColdResponse.DownDEGs.TIGK.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25Cvs37C_DownDEGs_TIGK_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p5
dev.off()

ColdResponse.DownDEGs.TIGK.Reactome.Top10<-ColdResponse.DownDEGs.TIGK.Reactome[1:10,]
ColdResponse.DownDEGs.TIGK.Reactome.Top10
frac <- ColdResponse.DownDEGs.TIGK.Reactome.Top10$Overlap
ColdResponse.DownDEGs.TIGK.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p6<-ggplot(ColdResponse.DownDEGs.TIGK.Reactome.Top10, aes(x=ColdResponse.DownDEGs.TIGK.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=ColdResponse.DownDEGs.TIGK.Reactome.Top10$Adjusted.P.value, size=(ColdResponse.DownDEGs.TIGK.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25Cvs37C_DownDEGs_TIGK_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p6
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of DEGs_Bubble Plots/-25Cvs37C_DownDEGs_TIGK_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p4,p5,p6, labels = "AUTO", nrow = 1)
dev.off()

###############################Data visualization for Enrichment for HSPs: Bubble plots ####
#HvsT_37C Control Comparisons ####
enriched.Ctrl.comparison.HSPs.HaCaT <- enrichr(Ctrl.comparison.HSPs.HaCaT$symbol, dbs)
enriched.Ctrl.comparison.HSPs.TIGK <- enrichr(Ctrl.comparison.HSPs.TIGK$symbol, dbs)

Ctrl.HaCaT.BP<-as.data.frame(enriched.Ctrl.comparison.HSPs.HaCaT$GO_Biological_Process_2023)
Ctrl.HaCaT.MF<-as.data.frame(enriched.Ctrl.comparison.HSPs.HaCaT$GO_Molecular_Function_2023)
Ctrl.HaCaT.Reactome<-as.data.frame(enriched.Ctrl.comparison.HSPs.HaCaT$Reactome_2022)


Ctrl.HaCaT.BP.Top10<-Ctrl.HaCaT.BP[1:10,]
Ctrl.HaCaT.BP.Top10
frac <- Ctrl.HaCaT.BP.Top10$Overlap
Ctrl.HaCaT.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(Ctrl.HaCaT.BP.Top10, aes(x=Ctrl.HaCaT.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Ctrl.HaCaT.BP.Top10$Adjusted.P.value, size=(Ctrl.HaCaT.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/37C_HaCaT_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

Ctrl.HaCaT.MF.Top10<-Ctrl.HaCaT.MF[1:10,]
Ctrl.HaCaT.MF.Top10
frac <- Ctrl.HaCaT.MF.Top10$Overlap
Ctrl.HaCaT.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(Ctrl.HaCaT.MF.Top10, aes(x=Ctrl.HaCaT.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=Ctrl.HaCaT.MF.Top10$`Gene Ratio`, colour=Ctrl.HaCaT.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/37C_HaCaT_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

Ctrl.HaCaT.Reactome.Top10<-Ctrl.HaCaT.Reactome[1:10,]
Ctrl.HaCaT.Reactome.Top10
frac <- Ctrl.HaCaT.Reactome.Top10$Overlap
Ctrl.HaCaT.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(Ctrl.HaCaT.Reactome.Top10, aes(x=Ctrl.HaCaT.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Ctrl.HaCaT.Reactome.Top10$Adjusted.P.value, size=(Ctrl.HaCaT.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/37C_HaCaT_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/37C_HaCaT_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, labels = "AUTO", nrow = 1)
dev.off()



enriched.Ctrl.comparison.HSPs.TIGK <- enrichr(Ctrl.comparison.HSPs.TIGK$symbol, dbs)
Ctrl.TIGK.BP<-as.data.frame(enriched.Ctrl.comparison.HSPs.TIGK$GO_Biological_Process_2023)
Ctrl.TIGK.MF<-as.data.frame(enriched.Ctrl.comparison.HSPs.TIGK$GO_Molecular_Function_2023)
Ctrl.TIGK.Reactome<-as.data.frame(enriched.Ctrl.comparison.HSPs.TIGK$Reactome_2022)


Ctrl.TIGK.BP.Top10<-Ctrl.TIGK.BP[1:10,]
Ctrl.TIGK.BP.Top10
frac <- Ctrl.TIGK.BP.Top10$Overlap
Ctrl.TIGK.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(Ctrl.TIGK.BP.Top10, aes(x=Ctrl.TIGK.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Ctrl.TIGK.BP.Top10$Adjusted.P.value, size=(Ctrl.TIGK.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/37C_TIGK_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

Ctrl.TIGK.MF.Top10<-Ctrl.TIGK.MF[1:10,]
Ctrl.TIGK.MF.Top10
frac <- Ctrl.TIGK.MF.Top10$Overlap
Ctrl.TIGK.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p5<-ggplot(Ctrl.TIGK.MF.Top10, aes(x=Ctrl.TIGK.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Ctrl.TIGK.MF.Top10$Adjusted.P.value, size=(Ctrl.TIGK.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/37C_TIGK_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p5
dev.off()

Ctrl.TIGK.MF.Significant<-Ctrl.TIGK.MF[1:2,]
Ctrl.TIGK.MF.Significant
frac <- Ctrl.TIGK.MF.Significant$Overlap
Ctrl.TIGK.MF.Significant$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))
p5<-ggplot(Ctrl.TIGK.MF.Significant, aes(x=Ctrl.TIGK.MF.Significant$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Ctrl.TIGK.MF.Significant$Adjusted.P.value, size=(Ctrl.TIGK.MF.Significant$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/37C_TIGK_MF_significant.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p5
dev.off()

Ctrl.TIGK.Reactome.Top10<-Ctrl.TIGK.Reactome[1:10,]
Ctrl.TIGK.Reactome.Top10
frac <- Ctrl.TIGK.Reactome.Top10$Overlap
Ctrl.TIGK.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p6<-ggplot(Ctrl.TIGK.Reactome.Top10, aes(x=Ctrl.TIGK.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Ctrl.TIGK.Reactome.Top10$Adjusted.P.value, size=(Ctrl.TIGK.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/37C_TIGK_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p6
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/37C_TIGK_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p4,p5,p6, labels = "AUTO", nrow = 1)
dev.off()

#HvsT_60C comparisons ####
enriched.Heat.comparison.HSPs.HaCaT <- enrichr(Heat.comparison.HSPs.HaCaT$symbol, dbs)

Heat.HaCaT.BP<-as.data.frame(enriched.Heat.comparison.HSPs.HaCaT$GO_Biological_Process_2023)
Heat.HaCaT.MF<-as.data.frame(enriched.Heat.comparison.HSPs.HaCaT$GO_Molecular_Function_2023)
Heat.HaCaT.Reactome<-as.data.frame(enriched.Heat.comparison.HSPs.HaCaT$Reactome_2022)


Heat.HaCaT.BP.Top10<-Heat.HaCaT.BP[1:10,]
Heat.HaCaT.BP.Top10
frac <- Heat.HaCaT.BP.Top10$Overlap
Heat.HaCaT.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(Heat.HaCaT.BP.Top10, aes(x=Heat.HaCaT.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Heat.HaCaT.BP.Top10$Adjusted.P.value, size=(Heat.HaCaT.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60C_HaCaT_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

Heat.HaCaT.MF.Top10<-Heat.HaCaT.MF[1:10,]
Heat.HaCaT.MF.Top10
frac <- Heat.HaCaT.MF.Top10$Overlap
Heat.HaCaT.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(Heat.HaCaT.MF.Top10, aes(x=Heat.HaCaT.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=Heat.HaCaT.MF.Top10$`Gene Ratio`, colour=Heat.HaCaT.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60C_HaCaT_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

Heat.HaCaT.Reactome.Top10<-Heat.HaCaT.Reactome[1:10,]
Heat.HaCaT.Reactome.Top10
frac <- Heat.HaCaT.Reactome.Top10$Overlap
Heat.HaCaT.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(Heat.HaCaT.Reactome.Top10, aes(x=Heat.HaCaT.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Heat.HaCaT.Reactome.Top10$Adjusted.P.value, size=(Heat.HaCaT.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60C_HaCaT_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60C_HaCaT_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, labels = "AUTO", nrow = 1)
dev.off()



enriched.Heat.comparison.HSPs.TIGK <- enrichr(Heat.comparison.HSPs.TIGK$symbol, dbs)

Heat.TIGK.BP<-as.data.frame(enriched.Heat.comparison.HSPs.TIGK$GO_Biological_Process_2023)
Heat.TIGK.MF<-as.data.frame(enriched.Heat.comparison.HSPs.TIGK$GO_Molecular_Function_2023)
Heat.TIGK.Reactome<-as.data.frame(enriched.Heat.comparison.HSPs.TIGK$Reactome_2022)


Heat.TIGK.BP.Top10<-Heat.TIGK.BP[1:10,]
Heat.TIGK.BP.Top10
frac <- Heat.TIGK.BP.Top10$Overlap
Heat.TIGK.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(Heat.TIGK.BP.Top10, aes(x=Heat.TIGK.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Heat.TIGK.BP.Top10$Adjusted.P.value, size=(Heat.TIGK.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60C_TIGK_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

Heat.TIGK.MF.Top10<-Heat.TIGK.MF[1:10,]
Heat.TIGK.MF.Top10
frac <- Heat.TIGK.MF.Top10$Overlap
Heat.TIGK.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p5<-ggplot(Heat.TIGK.MF.Top10, aes(x=Heat.TIGK.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Heat.TIGK.MF.Top10$Adjusted.P.value, size=(Heat.TIGK.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60C_TIGK_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p5
dev.off()

Heat.TIGK.MF.Significant<-Heat.TIGK.MF[1:4,]
Heat.TIGK.MF.Significant
frac <- Heat.TIGK.MF.Significant$Overlap
Heat.TIGK.MF.Significant$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p5<-ggplot(Heat.TIGK.MF.Significant, aes(x=Heat.TIGK.MF.Significant$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Heat.TIGK.MF.Significant$Adjusted.P.value, size=(Heat.TIGK.MF.Significant$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60C_TIGK_MF_significant.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p5
dev.off()


Heat.TIGK.Reactome.Top10<-Heat.TIGK.Reactome[1:10,]
Heat.TIGK.Reactome.Top10
frac <- Heat.TIGK.Reactome.Top10$Overlap
Heat.TIGK.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p6<-ggplot(Heat.TIGK.Reactome.Top10, aes(x=Heat.TIGK.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Heat.TIGK.Reactome.Top10$Adjusted.P.value, size=(Heat.TIGK.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60C_TIGK_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p6
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60C_TIGK_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p4,p5,p6, labels = "AUTO", nrow = 1)
dev.off()



#HvsT_-25C comparisons ####
enriched.Cold.comparison.HSPs.HaCaT <- enrichr(Cold.comparison.HSPs.HaCaT$symbol, dbs)

Cold.HaCaT.BP<-as.data.frame(enriched.Cold.comparison.HSPs.HaCaT$GO_Biological_Process_2023)
Cold.HaCaT.MF<-as.data.frame(enriched.Cold.comparison.HSPs.HaCaT$GO_Molecular_Function_2023)
Cold.HaCaT.Reactome<-as.data.frame(enriched.Cold.comparison.HSPs.HaCaT$Reactome_2022)


Cold.HaCaT.BP.Top10<-Cold.HaCaT.BP[1:10,]
Cold.HaCaT.BP.Top10
frac <- Cold.HaCaT.BP.Top10$Overlap
Cold.HaCaT.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(Cold.HaCaT.BP.Top10, aes(x=Cold.HaCaT.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Cold.HaCaT.BP.Top10$Adjusted.P.value, size=(Cold.HaCaT.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25C_HaCaT_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

Cold.HaCaT.MF.Top10<-Cold.HaCaT.MF[1:10,]
Cold.HaCaT.MF.Top10
frac <- Cold.HaCaT.MF.Top10$Overlap
Cold.HaCaT.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(Cold.HaCaT.MF.Top10, aes(x=Cold.HaCaT.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=Cold.HaCaT.MF.Top10$`Gene Ratio`, colour=Cold.HaCaT.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25C_HaCaT_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

Cold.HaCaT.Reactome.Top10<-Cold.HaCaT.Reactome[1:10,]
Cold.HaCaT.Reactome.Top10
frac <- Cold.HaCaT.Reactome.Top10$Overlap
Cold.HaCaT.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(Cold.HaCaT.Reactome.Top10, aes(x=Cold.HaCaT.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Cold.HaCaT.Reactome.Top10$Adjusted.P.value, size=(Cold.HaCaT.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25C_HaCaT_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25C_HaCaT_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, labels = "AUTO", nrow = 1)
dev.off()



enriched.Cold.comparison.HSPs.TIGK <- enrichr(Cold.comparison.HSPs.TIGK$symbol, dbs)

Cold.TIGK.BP<-as.data.frame(enriched.Cold.comparison.HSPs.TIGK$GO_Biological_Process_2023)
Cold.TIGK.MF<-as.data.frame(enriched.Cold.comparison.HSPs.TIGK$GO_Molecular_Function_2023)
Cold.TIGK.Reactome<-as.data.frame(enriched.Cold.comparison.HSPs.TIGK$Reactome_2022)


Cold.TIGK.BP.Top10<-Cold.TIGK.BP[1:10,]
Cold.TIGK.BP.Top10
frac <- Cold.TIGK.BP.Top10$Overlap
Cold.TIGK.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p4<-ggplot(Cold.TIGK.BP.Top10, aes(x=Cold.TIGK.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Cold.TIGK.BP.Top10$Adjusted.P.value, size=(Cold.TIGK.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25C_TIGK_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p4
dev.off()

Cold.TIGK.MF.Top10<-Cold.TIGK.MF[1:10,]
Cold.TIGK.MF.Top10
frac <- Cold.TIGK.MF.Top10$Overlap
Cold.TIGK.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p5<-ggplot(Cold.TIGK.MF.Top10, aes(x=Cold.TIGK.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Cold.TIGK.MF.Top10$Adjusted.P.value, size=(Cold.TIGK.MF.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25C_TIGK_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p5
dev.off()

Cold.TIGK.MF.Significant<-Cold.TIGK.MF[1:3,]
Cold.TIGK.MF.Significant
frac <- Cold.TIGK.MF.Significant$Overlap
Cold.TIGK.MF.Significant$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p5<-ggplot(Cold.TIGK.MF.Significant, aes(x=Cold.TIGK.MF.Significant$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Cold.TIGK.MF.Significant$Adjusted.P.value, size=(Cold.TIGK.MF.Significant$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25C_TIGK_MF_signficant.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p5
dev.off()

Cold.TIGK.Reactome.Top10<-Cold.TIGK.Reactome[1:10,]
Cold.TIGK.Reactome.Top10
frac <- Cold.TIGK.Reactome.Top10$Overlap
Cold.TIGK.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p6<-ggplot(Cold.TIGK.Reactome.Top10, aes(x=Cold.TIGK.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=Cold.TIGK.Reactome.Top10$Adjusted.P.value, size=(Cold.TIGK.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25C_TIGK_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p6
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25C_TIGK_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p4,p5,p6, labels = "AUTO", nrow = 1)
dev.off()


#HaCaT:60C vs 37C comparisons ####
enriched.HaCaT.HeatResponse.HSPs.upregulated <- enrichr(HaCaT.HeatResponse.HSPs.upregulated$symbol, dbs)

HeatResponse.HaCaT.BP<-as.data.frame(enriched.HaCaT.HeatResponse.HSPs.upregulated$GO_Biological_Process_2023)
HeatResponse.HaCaT.MF<-as.data.frame(enriched.HaCaT.HeatResponse.HSPs.upregulated$GO_Molecular_Function_2023)
HeatResponse.HaCaT.Reactome<-as.data.frame(enriched.HaCaT.HeatResponse.HSPs.upregulated$Reactome_2022)


HeatResponse.HaCaT.BP.Top10<-HeatResponse.HaCaT.BP[1:10,]
HeatResponse.HaCaT.BP.Top10
frac <- HeatResponse.HaCaT.BP.Top10$Overlap
HeatResponse.HaCaT.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(HeatResponse.HaCaT.BP.Top10, aes(x=HeatResponse.HaCaT.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=HeatResponse.HaCaT.BP.Top10$Adjusted.P.value, size=(HeatResponse.HaCaT.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60Cvs37C_HaCaT_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

HeatResponse.HaCaT.MF.Top10<-HeatResponse.HaCaT.MF[1:10,]
HeatResponse.HaCaT.MF.Top10
frac <- HeatResponse.HaCaT.MF.Top10$Overlap
HeatResponse.HaCaT.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(HeatResponse.HaCaT.MF.Top10, aes(x=HeatResponse.HaCaT.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=HeatResponse.HaCaT.MF.Top10$`Gene Ratio`, colour=HeatResponse.HaCaT.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60Cvs37C_HaCaT_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

HeatResponse.HaCaT.Reactome.Top10<-HeatResponse.HaCaT.Reactome[1:10,]
HeatResponse.HaCaT.Reactome.Top10
frac <- HeatResponse.HaCaT.Reactome.Top10$Overlap
HeatResponse.HaCaT.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(HeatResponse.HaCaT.Reactome.Top10, aes(x=HeatResponse.HaCaT.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=HeatResponse.HaCaT.Reactome.Top10$Adjusted.P.value, size=(HeatResponse.HaCaT.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60Cvs37C_HaCaT_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60Cvs37C_HaCaT_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, labels = "AUTO", nrow = 1)
dev.off()


#HaCaT:-25C vs 37C comparisons ####
enriched.HaCaT.ColdResponse.HSPs.upregulated <- enrichr(HaCaT.ColdResponse.HSPs.upregulated$symbol, dbs)

ColdResponse.HaCaT.BP<-as.data.frame(enriched.HaCaT.ColdResponse.HSPs.upregulated$GO_Biological_Process_2023)
ColdResponse.HaCaT.MF<-as.data.frame(enriched.HaCaT.ColdResponse.HSPs.upregulated$GO_Molecular_Function_2023)
ColdResponse.HaCaT.Reactome<-as.data.frame(enriched.HaCaT.ColdResponse.HSPs.upregulated$Reactome_2022)


ColdResponse.HaCaT.BP.Top10<-ColdResponse.HaCaT.BP[1:10,]
ColdResponse.HaCaT.BP.Top10
frac <- ColdResponse.HaCaT.BP.Top10$Overlap
ColdResponse.HaCaT.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(ColdResponse.HaCaT.BP.Top10, aes(x=ColdResponse.HaCaT.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=ColdResponse.HaCaT.BP.Top10$Adjusted.P.value, size=(ColdResponse.HaCaT.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25Cvs37C_HaCaT_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

ColdResponse.HaCaT.MF.Top10<-ColdResponse.HaCaT.MF[1:10,]
ColdResponse.HaCaT.MF.Top10
frac <- ColdResponse.HaCaT.MF.Top10$Overlap
ColdResponse.HaCaT.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(ColdResponse.HaCaT.MF.Top10, aes(x=ColdResponse.HaCaT.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=ColdResponse.HaCaT.MF.Top10$`Gene Ratio`, colour=ColdResponse.HaCaT.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25Cvs37C_HaCaT_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

ColdResponse.HaCaT.Reactome.Top10<-ColdResponse.HaCaT.Reactome[1:10,]
ColdResponse.HaCaT.Reactome.Top10
frac <- ColdResponse.HaCaT.Reactome.Top10$Overlap
ColdResponse.HaCaT.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(ColdResponse.HaCaT.Reactome.Top10, aes(x=ColdResponse.HaCaT.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=ColdResponse.HaCaT.Reactome.Top10$Adjusted.P.value, size=(ColdResponse.HaCaT.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25Cvs37C_HaCaT_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25Cvs37C_HaCaT_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, labels = "AUTO", nrow = 1)
dev.off()



#TIGK:60C vs 37C comparisons ####
enriched.TIGK.HeatResponse.HSPs.upregulated <- enrichr(TIGK.HeatResponse.HSPs.upregulated$symbol, dbs)

HeatResponse.TIGK.BP<-as.data.frame(enriched.TIGK.HeatResponse.HSPs.upregulated$GO_Biological_Process_2023)
HeatResponse.TIGK.MF<-as.data.frame(enriched.TIGK.HeatResponse.HSPs.upregulated$GO_Molecular_Function_2023)
HeatResponse.TIGK.Reactome<-as.data.frame(enriched.TIGK.HeatResponse.HSPs.upregulated$Reactome_2022)


HeatResponse.TIGK.BP.Top10<-HeatResponse.TIGK.BP[1:10,]
HeatResponse.TIGK.BP.Top10
frac <- HeatResponse.TIGK.BP.Top10$Overlap
HeatResponse.TIGK.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(HeatResponse.TIGK.BP.Top10, aes(x=HeatResponse.TIGK.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=HeatResponse.TIGK.BP.Top10$Adjusted.P.value, size=(HeatResponse.TIGK.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60Cvs37C_TIGK_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

HeatResponse.TIGK.MF.Top10<-HeatResponse.TIGK.MF[1:10,]
HeatResponse.TIGK.MF.Top10
frac <- HeatResponse.TIGK.MF.Top10$Overlap
HeatResponse.TIGK.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(HeatResponse.TIGK.MF.Top10, aes(x=HeatResponse.TIGK.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=HeatResponse.TIGK.MF.Top10$`Gene Ratio`, colour=HeatResponse.TIGK.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60Cvs37C_TIGK_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

HeatResponse.TIGK.MF.Signficant<-HeatResponse.TIGK.MF[1:5,]
HeatResponse.TIGK.MF.Signficant
frac <- HeatResponse.TIGK.MF.Signficant$Overlap
HeatResponse.TIGK.MF.Signficant$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(HeatResponse.TIGK.MF.Signficant, aes(x=HeatResponse.TIGK.MF.Signficant$`Gene Ratio`,y=Term))+
  geom_point(aes(size=HeatResponse.TIGK.MF.Signficant$`Gene Ratio`, colour=HeatResponse.TIGK.MF.Signficant$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60Cvs37C_TIGK_MF_Signficant.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

HeatResponse.TIGK.Reactome.Top10<-HeatResponse.TIGK.Reactome[1:10,]
HeatResponse.TIGK.Reactome.Top10
frac <- HeatResponse.TIGK.Reactome.Top10$Overlap
HeatResponse.TIGK.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(HeatResponse.TIGK.Reactome.Top10, aes(x=HeatResponse.TIGK.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=HeatResponse.TIGK.Reactome.Top10$Adjusted.P.value, size=(HeatResponse.TIGK.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60Cvs37C_TIGK_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

HeatResponse.TIGK.Reactome.Signficant<-HeatResponse.TIGK.Reactome[1:7,]
HeatResponse.TIGK.Reactome.Signficant
frac <- HeatResponse.TIGK.Reactome.Signficant$Overlap
HeatResponse.TIGK.Reactome.Signficant$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(HeatResponse.TIGK.Reactome.Signficant, aes(x=HeatResponse.TIGK.Reactome.Signficant$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=HeatResponse.TIGK.Reactome.Signficant$Adjusted.P.value, size=(HeatResponse.TIGK.Reactome.Signficant$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60Cvs37C_TIGK_reactome_Significant.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()


png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/60Cvs37C_TIGK_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, labels = "AUTO", nrow = 1)
dev.off()


#TIGK:-25C vs 37C comparisons ####
enriched.TIGK.ColdResponse.HSPs.upregulated <- enrichr(TIGK.ColdResponse.HSPs.upregulated$symbol, dbs)

ColdResponse.TIGK.BP<-as.data.frame(enriched.TIGK.ColdResponse.HSPs.upregulated$GO_Biological_Process_2023)
ColdResponse.TIGK.MF<-as.data.frame(enriched.TIGK.ColdResponse.HSPs.upregulated$GO_Molecular_Function_2023)
ColdResponse.TIGK.Reactome<-as.data.frame(enriched.TIGK.ColdResponse.HSPs.upregulated$Reactome_2022)


ColdResponse.TIGK.BP.Top10<-ColdResponse.TIGK.BP[1:10,]
ColdResponse.TIGK.BP.Top10
frac <- ColdResponse.TIGK.BP.Top10$Overlap
ColdResponse.TIGK.BP.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p1<-ggplot(ColdResponse.TIGK.BP.Top10, aes(x=ColdResponse.TIGK.BP.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=ColdResponse.TIGK.BP.Top10$Adjusted.P.value, size=(ColdResponse.TIGK.BP.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25Cvs37C_TIGK_BP.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p1
dev.off()

ColdResponse.TIGK.MF.Top10<-ColdResponse.TIGK.MF[1:10,]
ColdResponse.TIGK.MF.Top10
frac <- ColdResponse.TIGK.MF.Top10$Overlap
ColdResponse.TIGK.MF.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(ColdResponse.TIGK.MF.Top10, aes(x=ColdResponse.TIGK.MF.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(size=ColdResponse.TIGK.MF.Top10$`Gene Ratio`, colour=ColdResponse.TIGK.MF.Top10$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25Cvs37C_TIGK_MF.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()

ColdResponse.TIGK.MF.Significant<-ColdResponse.TIGK.MF[1:7,]
ColdResponse.TIGK.MF.Significant
frac <- ColdResponse.TIGK.MF.Significant$Overlap
ColdResponse.TIGK.MF.Significant$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p2<-ggplot(ColdResponse.TIGK.MF.Significant, aes(x=ColdResponse.TIGK.MF.Significant$`Gene Ratio`,y=Term))+
  geom_point(aes(size=ColdResponse.TIGK.MF.Significant$`Gene Ratio`, colour=ColdResponse.TIGK.MF.Significant$Adjusted.P.value))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))

png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25Cvs37C_TIGK_MF_significant.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p2
dev.off()


ColdResponse.TIGK.Reactome.Top10<-ColdResponse.TIGK.Reactome[1:10,]
ColdResponse.TIGK.Reactome.Top10
frac <- ColdResponse.TIGK.Reactome.Top10$Overlap
ColdResponse.TIGK.Reactome.Top10$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(ColdResponse.TIGK.Reactome.Top10, aes(x=ColdResponse.TIGK.Reactome.Top10$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=ColdResponse.TIGK.Reactome.Top10$Adjusted.P.value, size=(ColdResponse.TIGK.Reactome.Top10$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25Cvs37C_TIGK_reactome.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()

ColdResponse.TIGK.Reactome.Signficant<-ColdResponse.TIGK.Reactome[1:7,]
ColdResponse.TIGK.Reactome.Signficant
frac <- ColdResponse.TIGK.Reactome.Signficant$Overlap
ColdResponse.TIGK.Reactome.Signficant$"Gene Ratio"<-as.numeric(lapply(sapply(frac, function(x) eval(parse(text=x))), round, 3))

p3<-ggplot(ColdResponse.TIGK.Reactome.Signficant, aes(x=ColdResponse.TIGK.Reactome.Signficant$`Gene Ratio`,y=Term))+
  geom_point(aes(colour=ColdResponse.TIGK.Reactome.Signficant$Adjusted.P.value, size=(ColdResponse.TIGK.Reactome.Signficant$`Gene Ratio`)))+
  scale_color_gradientn(colours = rainbow(5))+
  labs(x='Gene Ratio', y=NULL,
       color='P.adjusted',size='Gene Ratio')+
  theme(axis.title = element_text(face='bold'),
        axis.text = element_text(face='bold', size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border =element_rect(color = "black", linewidth = 0.2, fill="NA"))+
  scale_y_discrete(labels = label_wrap(30))
png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25Cvs37C_TIGK_reactome_signficant.png",
    width = 1600,
    height = 2000,
    units = "px",
    res = 250)
p3
dev.off()


png("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/Enrichment of HSPs_Bubble Plots/-25Cvs37C_TIGK_BP+MF+Reactome.png",
    width = 4800,
    height = 2000,
    units = "px",
    res = 250)
plot_grid(p1,p2,p3, labels = "AUTO", nrow = 1)
dev.off()




#### Making heatmaps for the gene set 1-6 (genes upregulated in TIGK at baseline, compared to upregulated genes in HaCaT) ####
####Reordering the DESeq2 matrix 
#Reading in the data
Thermal.counts.reordered<-read.csv("C:/Users/chenh/Desktop/DiPietro Lab_Manuscripts/Thermal Injury to HaCaT and TIGK/gene_counts_protein coding_reordered.csv")

#Only reading mRNA/protein coding
Thermal.counts.reordered<-Thermal.counts.reordered%>%filter(gene_biotype == "protein_coding")
rownames(Thermal.counts.reordered)<-Thermal.counts.reordered[,1]
colnames(Thermal.counts.reordered)

#getting DESeq2 normalized matrix
coldata.Thermal.reordered<-matrix(c("TCtrl_1","TCtrl_2", "TCtrl_3","T60C_1","T60C_2","T60C_3","Tn25C_1","Tn25C_2","Tn25C_3",
                                    "HCtrl_1","HCtrl_2","HCtrl_3","H60C_1","H60C_2","H60C_3","Hn25C_1","Hn25C_2","Hn25C_3",
                                    "T_ctrl","T_ctrl","T_ctrl","T_Heat","T_Heat","T_Heat","T_cold","T_cold","T_cold",
                                    "H_ctrl","H_ctrl","H_ctrl","H_Heat","H_Heat","H_Heat","H_cold","H_cold","H_cold",
                                    "Ctrl","Ctrl","Ctrl", "Heat","Heat","Heat", "Cold","Cold","Cold",
                                    "Ctrl","Ctrl","Ctrl", "Heat","Heat","Heat", "Cold","Cold","Cold"),
                                  nrow = 18,
                                  ncol = 3,
                                  byrow = FALSE)
colnames(coldata.Thermal.reordered)<-c("sample","SampleTemp","Temp")

dds.Thermal.reordered.reordered<-DESeqDataSetFromMatrix(countData = Thermal.counts.reordered[,2:19],
                                                        colData = coldata.Thermal.reordered,
                                                        design = ~ SampleTemp)

#filtering all datasets: keep only avg counts >=10
keep<-rowMeans(counts(dds.Thermal.reordered))>=10 #could do rowsums too
dds.Thermal.reordered<-dds.Thermal.reordered[keep,] 

dds.Thermal.reordered<- DESeq(dds.Thermal.reordered)

#####P1 - Pulling out the HSPs upregulated in TIGK at baseline compared to genes upregulated in HaCaT after Heat Injury
TIGK.P1<-data.frame(Symbols = c("DNAJC12", "DNAJC22", "HSPA4L", "CRYAB", "BBS12", "DNAJA1", "DNAJB5", "CCT5", "DNAJB4", "DNAJB1", "HSPA13", "DNAJC18", "DNAJB2", "CCT8", "HSPA9", "DNAJC21", "DNAJA2",
                                "HSPB8", "DNAJC4", "HSPA2", "HSPB1", 
                                "HSPB3", "HSPA5", "HSP90B1", "HYOU1", "DNAJB11"))

TIGK.P1<-TIGK.P1 %>% 
  dplyr::inner_join(grch38, by = c("Symbols" = "symbol"))
TIGK.P1 <- TIGK.P1[!duplicated(TIGK.P1$ensgene),] 
na.omit(TIGK.P1)
TIGK.P1<-TIGK.P1[-13,] #remove duplicate
TIGK.P1<-TIGK.P1[-26,] #remove duplicate

dds.Thermal.reordered.count<-data.frame(counts(dds.Thermal.reordered, normalized = T))
dds.Thermal.reordered.count$ensembl<-row.names(dds.Thermal.reordered.count)
mat.TIGK.P1.HSPs<-subset(dds.Thermal.reordered.count, ensembl %in% TIGK.P1$ensgene)

#reordering to fit the right order for genes
reorder_idx <- match(TIGK.P1$ensgene,mat.TIGK.P1.HSPs$ensembl)
mat.TIGK.P1.HSPs<-mat.TIGK.P1.HSPs[reorder_idx,]

#deleting the ensembl column
colnames(mat.TIGK.P1.HSPs)
mat.TIGK.P1.HSPs<-mat.TIGK.P1.HSPs[,-19] #remove ensembl column
colnames(mat.TIGK.P1.HSPs)

#Making the heatmap for the gene sets 1-3 (hot)
mat.TIGK.P1.HSPs<-data.matrix(mat.TIGK.P1.HSPs)
mat.TIGK.P1.HSPs.z<-t(apply(mat.TIGK.P1.HSPs, 1, scale))

colnames(mat.TIGK.P1.HSPs.z)<-colnames(mat.TIGK.P1.HSPs)
Gene.symbols<-rownames(mat.TIGK.P1.HSPs.z)
Gene.symbols<-as.data.frame(Gene.symbols)
Gene.symbols<-Gene.symbols %>% 
  dplyr::inner_join(grch38, by = c("Gene.symbols" = "ensgene"))
na.omit(Gene.symbols)

pheatmap(mat.TIGK.P1.HSPs.z, color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         legend_breaks = -4:4, legend_labels = c("-4","-3", "-2", "-1", "0", "1","2","3","4"),
         border_color = "black",
         labels_row = Gene.symbols$symbol,
         cluster_cols = FALSE, cluster_row = FALSE,
         gaps_col = c(3, 6,9,12,15,18), gaps_row = c(17,21),
         fontsize_row = 20,
         fontsize_col = 10)

pheatmap(mat.TIGK.P1.HSPs.z,
         legend_breaks = -4:4, legend_labels = c("-4","-3", "-2", "-1", "0", "1","2","3","4"),
         border_color = "black",
         labels_row = Gene.symbols$symbol,
         cluster_cols = FALSE, cluster_row = FALSE,
         gaps_col = c(3, 6,9,12,15,18), gaps_row = c(17,21),
         fontsize_row = 10,
         fontsize_col = 10)


######P2 - Pulling out the HSPs upregulated in TIGK at baseline compared to genes upregulated in HaCaT after Cold Injury
TIGK.P2<-data.frame(Symbols = c("DNAJC12", "DNAJC22", "HSPA4L", "CRYAB", "BBS12", "DNAJA1", "DNAJB5", "CCT5", "DNAJB4", "DNAJB1", "HSPA13", "DNAJC18", "DNAJB2", "CCT8", "HSPA9", "DNAJC21", "DNAJA2",
                                "HSPB8", "DNAJC4", "HSPA2", "HSPB1", 
                                "DNAJC3", "HSPA5", "HSP90B1", "HYOU1"))

TIGK.P2<-TIGK.P2 %>% 
  dplyr::inner_join(grch38, by = c("Symbols" = "symbol"))
TIGK.P2 <- TIGK.P2[!duplicated(TIGK.P2$ensgene),] 
na.omit(TIGK.P2)
TIGK.P2<-TIGK.P2[-13,] #remove duplicate
TIGK.P2<-TIGK.P2[-26,] #remove duplicate

dds.Thermal.reordered.count<-data.frame(counts(dds.Thermal.reordered, normalized = T))
dds.Thermal.reordered.count$ensembl<-row.names(dds.Thermal.reordered.count)
mat.TIGK.P2.HSPs<-subset(dds.Thermal.reordered.count, ensembl %in% TIGK.P2$ensgene)

#reordering to fit the right order for genes
reorder_idx <- match(TIGK.P2$ensgene,mat.TIGK.P2.HSPs$ensembl)
mat.TIGK.P2.HSPs<-mat.TIGK.P2.HSPs[reorder_idx,]

#deleting the ensembl column
colnames(mat.TIGK.P2.HSPs)
mat.TIGK.P2.HSPs<-mat.TIGK.P2.HSPs[,-19] #remove ensembl column
colnames(mat.TIGK.P2.HSPs)

#Making the heatmap for the gene sets 4-6 (cold)
mat.TIGK.P2.HSPs<-data.matrix(mat.TIGK.P2.HSPs)
mat.TIGK.P2.HSPs.z<-t(apply(mat.TIGK.P2.HSPs, 1, scale))

colnames(mat.TIGK.P2.HSPs.z)<-colnames(mat.TIGK.P2.HSPs)
pheatmap(mat.TIGK.P2.HSPs.z, color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         legend_breaks = -4:4, legend_labels = c("-4","-3", "-2", "-1", "0", "1","2","3","4"),
         border_color = "black",
         labels_row = TIGK.P2$Symbols,
         cluster_cols = FALSE, cluster_row = FALSE,
         gaps_col = c(3, 6,9,12,15,18), gaps_row = c(17,21),
         fontsize_row = 20,
         fontsize_col = 10)

pheatmap(mat.TIGK.P2.HSPs.z, 
         legend_breaks = -4:4, legend_labels = c("-4","-3", "-2", "-1", "0", "1","2","3","4"),
         border_color = "black",
         labels_row = TIGK.P2$Symbols,
         cluster_cols = FALSE, cluster_row = FALSE,
         gaps_col = c(3,6,9,12,15,18), gaps_row = c(17,21),
         fontsize_row = 10,
         fontsize_col = 10)


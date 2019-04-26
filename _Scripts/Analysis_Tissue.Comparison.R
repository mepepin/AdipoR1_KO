# Developed 10/20/17
# Developer: Mark E. Pepin
Proj_Name="Koentges_Proteomics"
##Set working directory (relative path to saved .r file location)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#Load the libraries needed for this particular analysis
library(dplyr)
library(readxl)
library(pheatmap)
#Import the Tissue Index File
Mt.FXN_Liver<-read_xlsx("IPA_Analysis/Canonical Pathway Comparison/Mitochondrial.Dysfunction.xlsx", sheet="Liver")
Mt.FXN_Heart<-read_xlsx("IPA_Analysis/Canonical Pathway Comparison/Mitochondrial.Dysfunction.xlsx", sheet="Heart")
Mt.FXN_Kidney<-read_xlsx("IPA_Analysis/Canonical Pathway Comparison/Mitochondrial.Dysfunction.xlsx", sheet="Kidney")
#join the individual tissue lists to create an index table
Index.Row<-full_join(Mt.FXN_Liver, Mt.FXN_Heart, by="GeneSymbol")
Index.Row<-full_join(Index.Row, Mt.FXN_Kidney, by="GeneSymbol")
Index.Row[is.na(Index.Row)]<-""
#Convert the human gene name from IPA into the mouse gene (as described in the proteomics dataset)
Orthologous.Genes<-read.csv("ensembl_Human.Mouse.Orthologues.txt")
Index.Row_mouse<-merge(Index.Row, Orthologous.Genes, by.x="GeneSymbol", by.y="Gene.name")
Index.Row_mouse<-dplyr::select(Index.Row_mouse, GeneSymbol=Mouse.gene.name, Liver, Heart, Kidney)
Index.Row_mouse<-distinct(Index.Row_mouse)
#import expression dataset to determine directionality (heatmap-style)
  #Liver
Data.raw_Liver<-read_xlsx("../_Data/Liver - PN-0064-1 Liver - mep.xlsx")
Liver_Expression<-Data.raw_Liver %>% dplyr::select(GeneSymbol=`Gene names`, log2FC_liver=`log2 Ratio`, log.pval_liver=`-log P-Value`)
Liver_Expression<-dplyr::select(Liver_Expression, GeneSymbol, log2FC_liver, log.pval_liver)
Liver_Expression<-dplyr::filter(Liver_Expression, GeneSymbol!=is.na(GeneSymbol))
  #heart
Data.raw_Heart<-read_xlsx("../_Data/Heart - PN-0064 Heart - mep.xlsx")
Heart_Expression<-dplyr::select(Data.raw_Heart, GeneSymbol=`Gene names`, log2FC_heart=`log2 Ratio`, log.pval_heart=`-log P-Value`)
Heart_Expression<-dplyr::select(Heart_Expression, GeneSymbol, log2FC_heart, log.pval_heart)
Heart_Expression<-dplyr::filter(Heart_Expression, GeneSymbol!=is.na(GeneSymbol))
  #Kidney
Data.raw_Kidney<-read_xlsx("../_Data/Kidney - PN-0064-3 Ergebnisse - mep.xlsx")
Kidney_Expression<-dplyr::select(Data.raw_Kidney, GeneSymbol=`Gene names`, log2FC_kidney=`log2 Ratio`, log.pval_kidney=`-log P-Value`)
Kidney_Expression<-dplyr::select(Kidney_Expression, GeneSymbol, log2FC_kidney, log.pval_kidney)
Kidney_Expression<-dplyr::filter(Kidney_Expression, GeneSymbol!=is.na(GeneSymbol))
#Merge the tissues
Tissue.Expression<-full_join(Liver_Expression, Heart_Expression, by="GeneSymbol")
Tissue.Expression<-full_join(Tissue.Expression, Kidney_Expression, by="GeneSymbol")

#Create a heatmap of the changes
#filter the Tissue.Expression dataset by the Annotating Genes to pathway
HM.Data<-merge(Tissue.Expression, Index.Row_mouse, by="GeneSymbol")
HM.Data_significant<-HM.Data %>% filter(log.pval_heart>3 | log.pval_kidney>3 | log.pval_liver>3)
HM.Data_significant$GeneSymbol<-sub('[.]', '_', make.names(HM.Data_significant$GeneSymbol, unique=TRUE))
rownames(HM.Data_significant)<-HM.Data_significant$GeneSymbol
HM.input<-select(HM.Data_significant, log2FC_liver, log2FC_heart, log2FC_kidney)
HM.input[is.na(HM.input)]<-"0"
HM.input<-data.matrix(HM.input)
##Annotation above the Heat Map:
library(RColorBrewer)
color.hm<-rev(brewer.pal(11, "RdBu"))
pdf(file = paste0("Proteomics_Tissue.Heatmap_P<0..05.pdf"), width=10, height=7.5, onefile = FALSE)
pheatmap(HM.input, cluster_cols = FALSE, cluster_rows=FALSE, scale='row', color = color.hm, show_colnames=FALSE, show_rownames = TRUE, border_color = NA)
dev.off()


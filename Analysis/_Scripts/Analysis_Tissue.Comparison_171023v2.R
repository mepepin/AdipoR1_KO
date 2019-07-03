# Developed 10/20/17
# Developer: Mark E. Pepin
Proj_Name="Koentges_Proteomics"
##Set working directory (relative path to saved .r file location)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#Load the libraries needed for this particular analysis
library(dplyr)
library(readxl)
library(pheatmap)
library(RColorBrewer)
#Import the Tissue Index File
Mt.FXN_Liver<-read_xlsx("Input/Pathways/Mitochondrial_Dysfunction.xlsx", sheet="Liver")
Mt.FXN_Heart<-read_xlsx("Input/Pathways/Mitochondrial_Dysfunction.xlsx", sheet="Heart")
Mt.FXN_Kidney<-read_xlsx("Input/Pathways/Mitochondrial_Dysfunction.xlsx", sheet="Kidney")
#join the individual tissue lists to create an index table
Mt.FXN_Liver$Liver<-"Liver"
Mt.FXN_Heart$Heart<-"Heart"
Mt.FXN_Kidney$Kidney<-"Kidney"
Mt.FXN_Liver.Index<-select(Mt.FXN_Liver, Symbol, Liver)
Mt.FXN_Heart.Index<-select(Mt.FXN_Heart, Symbol, Heart)
Mt.FXN_Kidney.Index<-select(Mt.FXN_Kidney, Symbol, Kidney)
Index.Row<-full_join(Mt.FXN_Liver.Index, Mt.FXN_Heart.Index, by="Symbol")
Index.Row<-full_join(Index.Row, Mt.FXN_Kidney.Index, by="Symbol")
Index.Row[is.na(Index.Row)]<-0

## Venn Diagram

#Create a Heatmap of data
Mt.FXN_Heart.FC<-select(Mt.FXN_Heart, Symbol, Heart_FC=`Expr Log Ratio`)
Mt.FXN_Liver.FC<-select(Mt.FXN_Liver, Symbol, Liver_FC=`Expr Log Ratio`)
Mt.FXN_Kidney.FC<-select(Mt.FXN_Kidney, Symbol, Kidney_FC=`Expr Log Ratio`)
HM.Data_MtFXN<-full_join(Mt.FXN_Heart.FC, Mt.FXN_Liver.FC, by="Symbol")
HM.Data_MtFXN<-full_join(HM.Data_MtFXN, Mt.FXN_Kidney.FC, by="Symbol")  
rownames(HM.Data_MtFXN)<-HM.Data_MtFXN$Symbol
HM.Data_MtFXN$Heart_FC<-as.numeric(HM.Data_MtFXN$Heart_FC)
HM.Data_MtFXN$Liver_FC<-as.numeric(HM.Data_MtFXN$Liver_FC)
HM.Data_MtFXN$Kidney_FC<-as.numeric(HM.Data_MtFXN$Kidney_FC)
HM.Data_MtFXN[HM.Data_MtFXN==100]<-5
HM.Data_MtFXN[HM.Data_MtFXN==-100]<--5
HM.Input<-data.matrix(select(HM.Data_MtFXN, -Symbol))
HM.Input[is.na(HM.Input)]<-0
paletteLength <- 100
myColor <- colorRampPalette(c("purple", "dodgerblue4", "white", "brown4", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(min(HM.Input), seq(-4, 0, length.out=ceiling(paletteLength/2)), 
              seq(2/paletteLength, 4, length.out=floor(paletteLength/2)-1), max(HM.Input))
pdf(file = paste0("Proteomics_Tissue.Heatmap_P<0..05.pdf"), width=5, height=15, onefile = FALSE)
pheatmap(HM.Input, cluster_cols = TRUE, cluster_rows = TRUE, breaks = myBreaks, color = myColor, show_rownames = TRUE, border_color = NA)
dev.off()

#set the working directory
setwd("~/Box/Work/4_PhD/_Papers/To Do/Goetnes Manuscript/Pepin_Analysis")
#load necessary packages
library(dplyr)    #useful for cleaning/"tidying" data
library(pheatmap)    #best heatmap package I can find
library(readxl)   #import data from xlsx format
#Import Datasets (must preserve the relative file structure of the data/analysis folder for this to work)
Heart_raw <- read_excel("../_Data/Heart - PN-0064 Heart - mep.xlsx", sheet = "IPA_Import")
Kidney_raw <- read_excel("../_Data/Kidney - PN-0064-3 Ergebnisse - mep.xlsx", sheet = "IPA_Import")
Liver_raw <- read_excel("../_Data/Liver - PN-0064-1 Liver - mep.xlsx", sheet = "IPA_Import")
#full-join these 3 datasets (missing values are given an empty cell)
IPA_raw<-full_join(Heart_raw, Kidney_raw, by="Uniprot_ID")
IPA_raw<-select(IPA_raw, Uniprot_ID, Gene=Gene.x, qvalue_heart=qvalue.x, logpvalue_heart=minuslogpvalue.x, log2FC_heart=log2FC.x, qvalue_kidney=qvalue.y, logpvalue_kidney=minuslogpvalue.y, log2FC_kidney=log2FC.y)
IPA_raw<-full_join(IPA_raw, Liver_raw, by="Uniprot_ID")
IPA_raw<-select(IPA_raw, Uniprot_ID, Gene=Gene.x, -Gene.y, qvalue_heart, logpvalue_heart, log2FC_heart, qvalue_kidney, logpvalue_kidney, log2FC_kidney, qvalue_liver=qvalue, logpvalue_liver=minuslogpvalue, log2FC_liver=log2FC)
IPA_raw[is.na(IPA_raw)]<-""
##condense the dataset into only necessary columns
write.csv(IPA_raw, "IPA.upload_combined.proteomics_170922.csv")


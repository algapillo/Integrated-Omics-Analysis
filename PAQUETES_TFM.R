

## Descargas de paquetes ##
if(!require("BIocmanager", quietly =TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
BiocManager::install("iClusterPlus")
BiocManager::install("SummarizedExperiment")
BiocManager::install("sesame")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("EDASeq")
BiocManager::install("ggrpel")
BiocManager::install("Enmix")

#Paquete MOVIC y maxConway
devtools::install_github("xlucpu/MOVICS")
devtools::install_github("maxconway/SNFtool")

#

#Paquetes revision sistematica
install.packages("rvest")


###Paquetes analisis de Supervivencia
install.packages("survminer")
install.packages("survival")
install.packages("SummarirezExperiment")
install.packages("tidyverse")

#Paquete NOrmalization


#Paquetes Expresion diferencial
install.packages("limma")
BiocManager::install("edgeR")

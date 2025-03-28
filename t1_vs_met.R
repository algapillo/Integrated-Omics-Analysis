
HOMOGEINIZACION

trans= t(datos_trans_filtrados)





###########################################################
#Obtención de datos

##Visualizamos los datos metilómicos 
datos_met <- GDCquery(
  project = "TCGA-SKCM",
  data.category = "DNA Methylation",
  data.type="Methylation Beta Value", 
  platform = "Illumina Human Methylation 450",
  sample.type = c("Primary Tumor", "Metastatic")
)


##Descargamos los datos

GDCdownload(query= datos_met)

###Transformamos los datos en Summarized Experiment

library("sesame")
datos_met_se <- GDCprepare(
  query = datos_met)


##Visualizamos el numero de genes y muestras en la matriz de expresion
library(SummarizedExperiment)
mat_exp_datos_met_se = assay(datos_met_se)


#NUmero genes y muestras
dim(mat_exp_datos_met_se)


#Descartamos los valores nulos

datos_met_se_clean <- datos_met_se[complete.cases(assay(datos_met_se)), ]

matriz_datos_met_se_clean <- assay(datos_met_se_clean)




## Transponemos las matriz

met <-t(matriz_datos_met_se_clean)

dim(met)
#######################################################

################### PROTEINA ######################


##Visualizamos los datos proteómicos 
datos_prot <- GDCquery(
  project = "TCGA-SKCM",
  data.category = "Proteome Profiling"
)

GDCdownload(datos_prot)

###  Preprocesamiento datos prot ##
datos_prot=na.omit(datos_prot)

# Creamos Summarize Experiment de datos proteomicos
datos_prot_se <- GDCprepare(query = datos_prot
                          ,save =TRUE,
                          summarizedExperiment= TRUE) 

#Eliminamos Missings values
datos_prot_se<- datos_prot_se[complete.cases(datos_prot_se), ]


#Matriz
mat_exp_datos_prot_se  <- subset(datos_prot_se_clean, select = -c(1:5))

dim(mat_exp_datos_prot_se)





##################################

#####################################################

###################################################


trans= t(datos_trans_filtrados)

met= t(matriz_datos_met_se_clean)

prot <-t(mat_exp_datos_prot_se)
  
# Transformación datos trans
rownames(trans) <- substr(rownames(trans),1,16)
rownames(trans)[1]

#Transformación datos met
rownames(met) <- substr(rownames(met),1,16)
rownames(met)[1]

#Transformación datos prot
rownames(prot) <- substr(rownames(prot),1,16)
rownames(prot)[1]

#Muestras comunes prot-trans

prot_comun <- prot[rownames(prot) %in% rownames(trans),]
prot_comun <- prot_comun[order(rownames(prot_comun)),]
dim(prot_comun)

#Muestras comunes trans-prot

trans_comun <- trans[rownames(trans)  %in% rownames(prot_comun),]
trans_comun <- trans_comun[order(rownames(trans_comun)),]
dim(trans_comun)



#Muestras comunes met-prot

met_comun <- met[rownames(met)  %in% rownames(trans_comun),]
met_comun <- met_comun[order(rownames(met_comun)),]
dim(met_comun)

# Identificar los barcodes duplicados
duplicated_barcodes <- rownames(prot_comun)[duplicated(rownames(prot_comun))]


# Eliminar las filas con barcodes duplicados
prot_comun <- prot_comun[!rownames(prot_comun) %in% duplicated_barcodes, ]
trans_comun <- trans_comun[!rownames(trans_comun) %in% duplicated_barcodes, ]
met_comun <- met_comun[!rownames(met_comun) %in% duplicated_barcodes, ]

dim(trans_comun)
dim(met_comun)
dim(prot_comun)


#######################################################
#####################################################
######################################################
####################################################
################################################
###############################################
###############################################

#########   Integracion de datos    ############

#######################################################
#####################################################
######################################################
####################################################
################################################
###############################################
###############################################

####   Reduccion dimensional


## TRANS
tmp.a <-t(trans_comun)
elite.trans_comun.tmp <- getElites(dat = tmp.a,
                                      method    = "sd",
                                      elite.pct = 0.1,
                                      elite.num= 150
)




# MET 
tmp.a <-t(met_comun)
elite.met_comun.tmp <- getElites(dat       = tmp.a,
                                    method    = "pca",
                                    pca.ratio = 0.95,
                                    elite.num= 100
)




#PROT
tmp.a <-t(prot_comun)
elite.prot_comun.tmp <- getElites(dat       = tmp.a,
                                     method    = "pca",
                                     pca.ratio = 0.95,
                                     elite.num= 100
)


#TAMAÑO DE LAS MUESTRAS

dim(elite.trans_comun.tmp$elite.dat)
dim(elite.met_comun.tmp$elite.dat)
dim(elite.prot_comun.tmp$elite.dat)



mo.a <- list(omics1 = elite.trans_comun.tmp$elite.dat,
                omics2 = elite.met_comun.tmp$elite.dat,
                omics3=elite.prot_comun.tmp$elite.dat)

omics1 = elite.trans$elite.dat
omics2 = elite.met$elite.dat
omics3= elite.prot$elite.dat

###############################
###############################

##### Numero optimo de clustering

###########################
#########################

optk.brca <- getClustNum(data        = mo.a,
                         is.binary   = c(F,F,F),
                         try.N.clust = 2:8, # try cluster number from 2 to 8
                         fig.name    = "CLUSTER NUMBER OF TCGA-SKCM")



moic.res.list.k2.a <- getMOIC(data        = mo.a,
                                 methodslist = list("SNF","NEMO", "MoCluster",
                                                    "LRAcluster", "IntNMF", "iClusterBayes"),
                                 N.clust     = 2,
                                 type        = c("gaussian", "gaussian", "gaussian"))

cmoic.skcm.k2.a <- getConsensusMOIC(moic.res.list = moic.res.list.k2.a,
                                  fig.name      = "CONSENSUS HEATMAP k2",
                                  distance      = "euclidean",
                                  linkage       = "average",
)


#SILHOUETTE
sil_k2.a= getSilhouette(sil      = cmoic.skcm.k2.a$sil, # a sil object returned by getConsensusMOIC()
                      fig.path = getwd(),
                      fig.name = "SILHOUETTE_k2",
                      height   = 5.5,
                      width    = 5)

##############################
#################################


# Assuming datos_clinicos is your data frame


rownames(trans_comun) <- substr(rownames(trans_comun), 1,12)

rownames(trans_comun)[1]

#Extraemos tumores primarios de datos clinicos

rownames(datos_clinicos)
datos_clinicos_subid = datos_clinicos$submitter_id


tumor= datos_clinicos[datos_clinicos_subid %in% rownames(trans_comun),]
dim(tumor)

tumor <- tumor[order(rownames(tumor)),]


all(rownames(tumor) == rownames(trans_comun))

table(tumor$ajcc_pathologic_stage)


#Añadimos las columnas de los clusters

tumor$cluster = as.character(cmoic.skcm.k2.a$clust.res[,"clust"])

head(tumor$cluster)

tumor <- subset(tumor, select =c("ajcc_pathologic_stage","cluster"))

#Excluimos los valores nulos

tumor= na.omit(tumor)
dim(tumor)

library(caret)
library(stringr)

#Hacemos las agrupaciones en los 4 estadios

i <- 1
while (i <= dim(tumor)[1]) {
  if (str_detect(tumor$ajcc_pathologic_stage[i], "Stage IV")) {
    tumor$ajcc_pathologic_stage[i] <- "Stage II"
  } else {
    if (str_detect(tumor$ajcc_pathologic_stage[i], "Stage III")) {
      tumor$ajcc_pathologic_stage[i] <- "Stage II"
    } else {
      if (str_detect(tumor$ajcc_pathologic_stage[i], "Stage II")) {
        tumor$ajcc_pathologic_stage[i] <- "Stage I"
      } else {
        if (str_detect(tumor$ajcc_pathologic_stage[i], "Stage I")) {
          tumor$ajcc_pathologic_stage[i] <- "Stage I"
        }
      }
    }
  }
  i <- i + 1
}



i <- 1
while (i <= dim(tumor)[1]) {
  if (str_detect(tumor$cluster[i], "1")) {
    tumor$cluster[i] <- "Stage II"
  }
  else
    if (str_detect(tumor$cluster[i], "2")) {
      tumor$cluster[i] <- "Stage I"
    }
  i <- i + 1
}




tumor_filtered <- tumor[tumor$ajcc_pathologic_stage != "Not Reported", ]

# Recomputa la matriz de confusión
confusionMatrix(factor(tumor_filtered$ajcc_pathologic_stage), factor(tumor_filtered$cluster))


# Recomputa la matriz de confusión
confusionMatrix(factor(tumor$ajcc_pathologic_stage), factor(tumor$cluster))

mo.data_tp <- list(omics1 = elite.trans_t1_comun.tmp$elite.dat,
                 omics3=elite.prot_t1_comun.tmp$elite.dat)

moic.res.list.k2_tp <- getMOIC(data        = mo.data_tp,
                            methodslist = list("SNF", "CIMLR", "PINSPlus", "NEMO", "COCA", "MoCluster",
                                               "LRAcluster", "ConsensusClustering", "IntNMF", "iClusterBayes"),
                            N.clust     = 2,
                            type        = c("gaussian", "gaussian"))


cmoic.skcm.k2_tp <- getConsensusMOIC(moic.res.list = moic.res.list.k2_tp,
                                  fig.name      = "CONSENSUS HEATMAP k2",
                                  distance      = "euclidean",
                                  linkage       = "average",
)

rownames(prot_t1_comun) <- substr(rownames(prot_t1_comun), 1, 12)

rownames(prot_t1_comun)[1]

#Extraemos tumores primarios de datos clinicos

rownames(datos_clinicos)
datos_clinicos_subid = datos_clinicos$submitter_id


tumor_prim_tp= datos_clinicos[datos_clinicos_subid %in% rownames(prot_t1_comun),]
dim(tumor_prim_tp)


tumor_prim_tp <- tumor_prim_tp[order(rownames(tumor_prim_tp)),]

all(rownames(tumor_prim_tp) == rownames(prot_t1_comun))

table(tumor_prim_tp$ajcc_pathologic_stage)


#Añadimos las columnas de los clusters

tumor_prim_tp$cluster = as.character(cmoic.skcm.k2_tp$clust.res[,"clust"])

head(tumor_prim_tp$cluster)

tumor_prim_tp <- subset(tumor_prim_tp, select =c("ajcc_pathologic_stage","cluster"))

#Excluimos los valores nulos

tumor_prim_tp= na.omit(tumor_prim_tp)
dim(tumor_prim_tp)



#Hacemos las agrupaciones en los 4 estadios

i <- 1
while (i <= dim(tumor_prim_tp)[1]) {
  if (str_detect(tumor_prim_tp$ajcc_pathologic_stage[i], "Stage IV")) {
    tumor_prim_tp$ajcc_pathologic_stage[i] <- "Stage II"
  } else {
    if (str_detect(tumor_prim_tp$ajcc_pathologic_stage[i], "Stage III")) {
      tumor_prim_tp$ajcc_pathologic_stage[i] <- "Stage II"
    } else {
      if (str_detect(tumor_prim_tp$ajcc_pathologic_stage[i], "Stage II")) {
        tumor_prim_tp$ajcc_pathologic_stage[i] <- "Stage I"
      } else {
        if (str_detect(tumor_prim_tp$ajcc_pathologic_stage[i], "Stage I")) {
          tumor_prim_tp$ajcc_pathologic_stage[i] <- "Stage I"
        }
      }
    }
  }
  i <- i + 1
}



i <- 1
while (i <= dim(tumor_prim_tp)[1]) {
  if (str_detect(tumor_prim_tp$cluster[i], "1")) {
    tumor_prim_tp$cluster[i] <- "Stage II"
  }
  else
    if (str_detect(tumor_prim_tp$cluster[i], "2")) {
      tumor_prim_tp$cluster[i] <- "Stage I"
    }
  i <- i + 1
}




tumor_prim_tp_filtered <- tumor_prim_tp[tumor_prim_tp$ajcc_pathologic_stage != "Not Reported", ]

# Recomputa la matriz de confusión
confusionMatrix(factor(tumor_prim_tp_filtered$ajcc_pathologic_stage), factor(tumor_prim_tp_filtered$cluster))
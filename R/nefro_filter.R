nefro_filter <- function(df){
  #Comprobamos que esten todos los campos correspondientes a loci de receptores y donantes, aquellos que no esten se añaden como NA
  loci_HLA <- c("PacHLA1", "PacHLA2", "PacHLB1", "PacHLB2", "PacHLC1", "PacHLC2", "PacHLDr1", "PacHLDr2", "PacHLDQ1", "PacHLDQ2", "PacHLDQA1",
                "PacHLDQA2", "Don HLA A1", "Don HLA A2", "Don HLA B1", "Don HLA B2", "Don HLA C1", "Don HLA C2", "Don HLA Dr1", "Don HLA Dr2",
                "Don HLA DQ1", "Don HLA DQ2", "Don HLA DQA1", "Don HLA DQA2")
  for (locus in loci_HLA){
    if (!(locus %in% names(df))){
      df[[locus]] <- NA
    }
  }
  #Eliminamos las columnas correspondientes a DP, el algoritmo de imputación no las acepta como input ni los devuelve como output
  df <- df[, !grepl("DP", names(df))]

  #Dividimos en dos dataframes distintos los datos de donantes y receptores
  receptors <- df[, which(names(df) == "NumTx" | names(df) == "NumHistoria" | grepl("^Pac", names(df)))]
  receptors$Sample <- "Receptor"
  receptors %<>% select(NumTx, NumHistoria, Sample, everything())


  donors <- df[, which(names(df) == "NumTx" | names(df) == "NumHistoria" | grepl("^Don", names(df)))]
  donors$Sample <- "Donor"
  donors <- donors %>%
    mutate(NumHistoria = row_number())
  donors %<>% select(NumTx, NumHistoria, Sample, everything())

  #Comprobación del número de columnas de ambos dataframes
  df_names <-  c("ID_Tx", "Número", "Sample", "A1", "A2", "B1", "B2", "C1", "C2", "DRB11", "DRB12"
                 , "DQB11", "DQB12", "DQA1", "DQA2")
  if (length(df_names) != ncol(donors) | length(df_names) != ncol(receptors)){
    warning("Número de columnas de df_donors o df_receptors erroneo, comprobar las columnas de la tabla original.")
  }
  #Se renombran las columnas de ambos y se fusionan
  names(receptors) <- df_names
  names(donors) <- df_names
  df_filtered <- arrange(rbind(receptors, donors), ID_Tx)

  prefixes <- c("A*", "A*", "B*", "B*", "C*", "C*", "DRB1*", "DRB1*", "DQB1*", "DQB1*", "DQA1*", "DQA1*")
  for (i in 1:length(prefixes)){
    df_filtered[[i+3]][!is.na(df_filtered[[i+3]])] <- paste0(prefixes[i], df_filtered[[i+3]][!is.na(df_filtered[[i+3]])])
  }

  #Se añaden columnas DRB345 y se reordena el dataframe final
  df_filtered$DRB3451 <- NA
  df_filtered$DRB3452 <- NA
  df_filtered %<>% select(Número, ID_Tx, Sample, A1, A2, B1, B2, C1, C2, DRB11, DRB12, DRB3451, DRB3452, everything())
}



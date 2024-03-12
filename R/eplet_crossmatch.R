eplet_crossmatch <- function(eplet_df){
  regex <- c("^A[0-9]_eplet$", "^B[0-9]_eplet$", "^C[0-9]_eplet$", "^DRB1[0-9]_eplet$","^DRB345[0-9]_eplet$", "^DQB1[0-9]_eplet$", "^DQA[0-9]_eplet$")
  new_cols <- c("A_eplets", "B_eplets", "C_eplets", "DRB1_eplets", "DRB345_eplets", "DQB1_eplets", "DQA1_eplets", "HLA_I_eplets", "HLA_II_eplets", "HLA_all_eplets")

  #Funcion que obtiene los valores de varias celdas, los pega, separa todos los elementos que esten delimitados por un espacio en una lista,
  #obtiene los elementos únicos, los ordena y los vuelve a concatenar en un string
  split_unique_collapse <- function(cell) {
    all_values <- paste(cell, collapse = " ")
    unique_vals <- unique(strsplit(all_values, " ")[[1]])
    non_na_unique_vals <- sort(unique_vals[unique_vals != "NA"])
    collapsed_str <- paste(non_na_unique_vals, collapse = " ")
    return(collapsed_str)
  }

  #Loop que itera por todas las filas, evalua todas las regex extrayendo las columnas que las cumplan, almacena los indices, obtiene los eplets unicos
  #y los almacena en una nueva columna (p.e. obtiene eplets de A1_eplet y A2_eplet y almacena los eplets unicos en una nueva columna A_eplets)
  for (i in 1:nrow(eplet_df)){
    for (j in 1:length(regex)){
      indexes <- which(grepl(regex[j], names(eplet_df)))
      eplet_df[[new_cols[j]]] <- apply(eplet_df[indexes], 1, split_unique_collapse)

    }
  }

  #Codigo para la comparación de eplets únicos entre receptor y donante para cada loci: obtendremos solo los eplets que esten en el donante pero no en el receptor
  #Iteramos las nuevas columnas que hemos creado
  for (col in new_cols) {
    #Creamos nuevas columnas con el sufijo "_MM" que almacenaran los missmatches para cada loci
    MM_col <- gsub("_eplets", "_MM", col)

    #Itera por los valores úniocos de ID_Tx
    for (id_tx in unique(eplet_df$ID_Tx)) {
      #Para donantes y receptores obtiene los eplets de la columna correspondiente y los almacena en una nueva variable
      receptor_vals <- unlist(strsplit(as.character(eplet_df[eplet_df$ID_Tx == id_tx & eplet_df$Sample == "Receptor", col]), " "))
      donante_vals <- unlist(strsplit(as.character(eplet_df[eplet_df$ID_Tx == id_tx & eplet_df$Sample == "Donor", col]), " "))

      #Obtenemos los eplets que contiene el donante y que no estan en el receptor, si es 0 se asigna NA y si no se concatena el resultado en un string
      #diferencia_vals <- setdiff(donante_vals, receptor_vals)
      MM_val <- ifelse(length(setdiff(donante_vals, receptor_vals)) > 0, paste(setdiff(donante_vals, receptor_vals), collapse = " "), NA)

      #Asignamos el valor obtenido a la columna correspondiente del receptor
      eplet_df[eplet_df$ID_Tx == id_tx & eplet_df$Sample == "Receptor", MM_col] <- MM_val
    }
  }

  #Creamos un nuevo dataframe solo con los datos de interés
  eplet_MM <- eplet_df[eplet_df$Sample == "Receptor", c("Número", "ID_Tx", "Cw_ignored", "min_loci_allowed", grep("_MM$", colnames(eplet_df), value = TRUE))]

  MM_cols <- c("A_MM", "B_MM", "C_MM", "DRB1_MM", "DRB345_MM", "DQB1_MM", "DQA1_MM", "HLA_I_MM", "HLA_II_MM", "HLA_all_MM")

  #Bucle para crear las columnas que contabilizan el numero de MM por locus, por loci de clase I, de clase II y en total
  for (col in MM_cols) {
    count_col <- paste0(col, " (count)")
    eplet_MM[[count_col]] <- sapply(eplet_MM[[col]], function(x) ifelse(is.na(x), 0, length(strsplit(as.character(x), " ")[[1]])))
  }

  #Creación de la columna Risk_Score y calculo del riesgo en base al MM de eplets en DRB1, DRB345, DQB1, DQA1 segñún el algoritmo de la librería hlaR
  eplet_MM$Risk_Score <- NA
  eplet_MM %<>% select(Número, ID_Tx, Cw_ignored, min_loci_allowed, Risk_Score, grep("count", names(.), value = TRUE), everything())

  for (i in 1:nrow(eplet_MM)){
    DQ <- sum(eplet_MM[i, "DQB1_MM (count)"], eplet_MM[i, "DQA1_MM (count)"])
    DR <- sum(eplet_MM[i, "DRB1_MM (count)"], eplet_MM[i, "DRB345_MM (count)"])
    eplet_MM[i, "Risk_Score"] <- ifelse(between(DQ, 15, 31), "High",
                                        ifelse((DR >= 7 & DQ <= 14) | (DR < 7 & between(DQ, 9, 15)), "Moderate",
                                               ifelse( DR < 7 & DQ < 9, "Low", "Out of bond")))

  }

  return(eplet_MM)
}

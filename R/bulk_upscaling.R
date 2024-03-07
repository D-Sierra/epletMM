bulk_upscaling <- function(low_list = "low_list", df = "subdf1", haplotypes_path = getwd()){
  high_list <- list()
  path = paste0(haplotypes_path,"/A~C~B~DRB3-4-5~DRB1~DQB1.xlsx")
  #Loop para la imputacion usando hlapro::upscale_typings()
  for (i in seq_along(low_list)){
    highHLA <- hlapro::upscale_typings(filepath = path,
                               typing = low_list[[i]],
                               loci_input = c("A", "B", "C", "DRB1", "DRB.", "DQB1"))
    #Para algunos tipados la función hlapro::upscale_typings() no devuelve resultado si usa un DQ que ha sido previamente
    #pasado a resolucion serologica usando hlapro::get_serology() en estos casos se vuelve a realizar la imputación eliminando el locus DQ

    #Si tras la primera imputacion no se obtienen resultados se elimina el locus DQ y realiza una nueva imputacion
    if (is.na(highHLA[1,"unphased_geno"])){
      highHLA <- hlapro::upscale_typings(filepath = path,
                                 typing = gsub("DQ[0-9]{1,2}", "", low_list[[i]]),
                                 loci_input = c("A", "B", "C", "DRB1", "DRB.", "DQB1"))
      #Si tras la segunda imputacion no tenemos resultados se añade el resultado como NA y se imprime un mensaje de advertencia
      if (is.na(highHLA[1,"unphased_geno"])){
        high_list <- append(high_list, list(highHLA[1, "unphased_geno"]))
        names(high_list[[i]]) <- names(low_list[i])
        warning(paste(names(low_list[i]), "no pudo ser imputado."))
        #Si la segunda imputacion obtiene un resultado entonces se añade a la lista y se imprime un mensaje indicando que se ha realizado sin usar DQ
      } else {
        high_list <- append(high_list, list(highHLA[1, "unphased_geno"]))
        names(high_list[[i]]) <- names(low_list[i])
        warning(paste(names(low_list[i]), "fue imputado ignorando el locus DQ."))
      }
      #Si la primera imputacion es correcta se añade directamente a la posicion correspondiente de high_list
    } else {
      #Asignación del primer genotipo a la lista high_list
      high_list <- append(high_list, list(highHLA[1, "unphased_geno"]))
      names(high_list[[i]]) <- names(low_list[i])
    }
  }
  high_list <- high_list[!sapply(high_list, is.na)]

  #Creación de un nuevo dataframe para poblar con los resultados de la lista high_list
  imputed_df <- as.data.frame(t(as.data.frame(high_list, row.names = "high_typing"))) %>%
    mutate(Número = sub("^X", "", rownames(.))) %>% #Pasar la lista a dataframe añade X al principio del nombre de cada elemento, se borra y se asigna a una columna
    relocate(high_typing, high_typing, .after = Número) %>%
    separate(high_typing, into = c("A1" , "A2", "B1", "B2", "C1", "C2", "DQB11", "DQB12", "DRB11", "DRB12", "DRB3451", "DRB3452"), sep = " ") %>%
    mutate(across(everything(), ~sub("g", "", .))) #Se eliminan las "g" añadidas al final de algunos alelos por upscale_typings()
  row.names(imputed_df) <- NULL #Se eliminan los rownames

  for (i in seq_along(imputed_df$Número)) {
    #Para cada fila de imputed_df encuentra la fila en df con el mismo valor para "Número" y obtiene el índice
    index_df <- which(df$Número == imputed_df$Número[i])

    #Si el índice existe (>0) itera todas las columnas de imputed_df
    if (length(index_df) > 0) {
      for (col_name in names(imputed_df)[-1]) {
        #Para no sobreescribir tipajes ya en int/alta si pone la condición de que solo se actualizan las celdas que no contengan ":" en df
        if (!any(grepl(":", df[index_df, col_name]))) {
          df[index_df, col_name] <- imputed_df[i, col_name]
        }
      }
    }
  }

  #Tabla de equivalencias para alelos DQB1, DRB1, DQA1
  equivalence_table <- data.frame(
    DQB1 = c(".*\\*02:01.*", ".*\\*02:02.*", ".*\\*02:02.*", ".*\\*03:01.*", ".*\\*03:01.*", ".*\\*03:02.*", ".*\\*03:03.*", ".*\\*03:03.*", ".*\\*03:19.*", ".*\\*04.*", ".*\\*04.*",
             ".*\\*05:01.*", ".*\\*05:02.*", ".*\\*05:03.*", ".*\\*06:01.*", ".*\\*06:02.*", ".*\\*06:03.*", ".*\\*06:04.*", ".*\\*06:08.*", ".*\\*06:09.*"),
    DRB1 = c(NA, ".*\\*04:05.*", "^(?!.*\\*04:05).*$", ".*\\*04:01.*", "^(?!.*\\*04:01).*$", NA, ".*\\*07:01.*", "^(?!.*\\*07:01).*$", NA,  ".*\\*08.*",
             "^(?!.*\\*08).*$", NA, NA,NA,NA,NA,NA,NA,NA,NA),
    DQA1 = c("DQA1*05:01", "DQA1*03:01", "DQA1*02:01", "DQA1*03:01", "DQA1*05:01", "DQA1*03:01", "DQA1*02:01","DQA1*03:01", "DQA1*05:05", "DQA1*04:01", "DQA1*03:01", "DQA1*01:01", "DQA1*01:02",
             "DQA1*01:01", "DQA1*01:01", "DQA1*01:02", "DQA1*01:03", "DQA1*01:02", "DQA1*01:01", "DQA1*01:01")
  )

  #Asignación de DQA1 en base al tipaje de DQB1 y DRB1
  for (i in 1:nrow(df)) {
    for (j in 1:2) {
      dqb_value <- ifelse(j == 1, df[i, "DQB11"], df[i, "DQB12"])
      drb_value <- paste(df[i, "DRB11"], ", ", df[i, "DRB12"])

      index <- which(sapply(equivalence_table[, "DQB1"], function(x) any(grepl(x, dqb_value, perl = TRUE))))

      if (length(index) == 1) {
        df[i, ifelse(j == 1, "DQA1", "DQA2")] <- equivalence_table[index, "DQA1"]
      } else if (length(index) > 1) {
        index2 <- sapply(equivalence_table[index, "DRB1"], function(regex) any(grepl(regex, drb_value, perl = TRUE)))
        if (any(index2)) {
          df[i, ifelse(j == 1, "DQA1", "DQA2")] <- equivalence_table[index[index2], "DQA1"]
        } else {
          next
        }
      } else {
        next
      }
    }
  }
  return(df)
}

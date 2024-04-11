preimputation_HLA <- function(df = "mixed_reso"){

  #Loop para rebajar la resolución de alelos de intermedia a baja
  df_serology <- df
  for (i in 1:nrow(df_serology)) {
    for (j in 4:ncol(df_serology)){
      if (grepl(":", df_serology[i, j])){
        df_serology[i,j] <- hlapro::get_broad(df[i,j])
      }
    }
  }

  #Evaluamos un criterio de minimo de loci tipados para poder hacer la imputacion, que seran A, B y DR
  df$min_loci_allowed <- TRUE
  for (i in 1:nrow(df_serology)){
    if(any(is.na(df_serology[i,c("A1", "A2", "B1", "B2", "DRB11", "DRB12")]))){
      df[i, "min_loci_allowed"] <- FALSE
    }
  }
  #Filtramos la matriz serológica solo cuando cumpla el criterio de minimos loci tipados
  df_serology <- subset(df_serology, df$min_loci_allowed == TRUE)

  #Añadimos a cada celda el prefijo de locus correspondiente (pe. A, B, C...)
  # Iterar sobre las filas y columnas del dataframe
  for (i in 1:nrow(df_serology)) {
    for (j in 4:ncol(df_serology)) {

      # Obtener el valor de la celda actual
      valor_celda <- df_serology[i, j]

      # Verificar si el valor no es nulo
      if (!is.na(valor_celda)) {
        # Aplicar las sustituciones usando gsub
        df_serology[i, j] <- gsub("^C\\*|^Cw", "Cw*", valor_celda)
        df_serology[i, j] <- gsub("^DRB1\\*|^DR", "DR*", df_serology[i, j])
        df_serology[i, j] <- gsub("^DQB1\\*|^DQ", "DQ*", df_serology[i, j])
        df_serology[i, j] <- gsub("^B(\\d+)", "B*\\1", df_serology[i, j])
        df_serology[i, j] <- gsub("^A(\\d+)", "A*\\1", df_serology[i, j])
      }
    }
  }

  #Algunos alelos serólogicos del locus C no son reconocidos correctamente por hlapro::upscale_typings() y devuelve un resultado nulo
  #Por esta razón es mejor eliminar estos alelos (Cw8, Cw11, Cw13, Cw16 y Cw17) antes de la imputación
  #Creamos un nueva columna en df para guardar los tipajes que se han obtenido sin tener en cuenta el locus C
  df$Cw_ignored <- FALSE
  #Se evalua si alguno de estos alelos se encuentra en las columnas C1 o C2, y en caso afirmativo se asigna NA en df_serology y TRUE al registro en df
  for (i in 1:nrow(df_serology)){
    if (any(c("Cw*8", "Cw*11", "Cw*13", "Cw*16", "Cw*17") %in% df_serology[i, "C1"] |
            c("Cw*8", "Cw*11", "Cw*13", "Cw*16", "Cw*17") %in% df_serology[i, "C2"])){
      df_serology[i, "C1"] <- NA
      df_serology[i, "C2"] <- NA
      df[i,"Cw_ignored"] <- TRUE
    } else {
      next
    }
  }

  loci_test <- data.frame(matrix(rep(NA,14), nrow = 1))
  names(loci_test) <- names(df_serology[,-c(1, 2, 3)])
  checks <- c("^A\\*\\d{1,2}$", "^A\\*\\d{1,2}$",
              "^B\\*\\d{1,2}$", "^B\\*\\d{1,2}$",
              "^Cw\\*\\d{1,2}$", "^Cw\\*\\d{1,2}$",
              "^DR\\*\\d{1,2}$", "^DR\\*\\d{1,2}$",
              "^DRY\\*\\d{1,2}$", "^DRY\\*\\d{1,2}$",
              "^DQ\\*\\d{1,2}$", "^DQ\\*\\d{1,2}$",
              "^DQA\\*\\d{1,2}$", "^DQA\\*\\d{1,2}$")

  failed_elements_info <- character(0)

  #Comprueba para cada locus si todos los alelos de df_serology tienen una estrutura compatible con un alelo correcto
  for (i in seq_along(loci_test)) {
    column_name <- names(loci_test)[i]
    pattern <- checks[i]

    condition_result <- grepl(pattern, df_serology[[column_name]]) | is.na(df_serology[[column_name]])

    if (all(condition_result)) {
      loci_test[[column_name]] <- "PASS"
    } else {
      loci_test[[column_name]] <- "FAIL"
      warning(paste("Loci test failed for column ", column_name))

      #Registra los elementos que no han pasado el control de calidad en el vector failed_elements
      failed_elements <- df_serology[!condition_result, c("Número", column_name)]
      for (j in seq_len(nrow(failed_elements))) {
        failed_elements_info <- c(
          failed_elements_info,
          paste0(failed_elements[j, "Número"]," (",failed_elements[j, column_name], ")")
        )
      }
    }
  }

  #Imprimir la información de los elementos que eno han pasado el control de calidad
  if (length(failed_elements_info) > 0) {
    failed_elements_info <- paste(failed_elements_info, collapse = ", ")
    warning(paste("Elements that failed quality test:", (failed_elements_info)))
  } else {
    message("Quality test passed for all loci.")
  }

  #Reordenamos las columnas de df
  df %<>% select(Número, ID_Tx, Sample, Cw_ignored, min_loci_allowed, everything())

  return(curated = list(df, df_serology, loci_test))
}

modulab_HLA_cleaner <- function(x = "HLA_baja"){
  #Se evalua las columnas para todos los loci (pe.LOCUS A) existe en la tabla, de no ser así se crea
  loci_HLA <- c("LOCUS A", "LOCUS B", "LOCUS C", "LOCUS DRB1", "LOCUS DRB345","LOCUS DQB1", "LOCUS DQA1")

  for (locus in loci_HLA){
    if (!(locus %in% names(HLA_baja))){
      HLA_baja[[locus]] <- NA
    }
  }

  #Creamos una tabla intermedia sobre la que trabajaremos y que solo contiene el numero de muestra y los valores de todos los LOCUS
  df <- subset(HLA_baja, select = c("Número", "LOCUS A", "LOCUS B", "LOCUS C", "LOCUS DRB1", "LOCUS DRB345","LOCUS DQB1", "LOCUS DQA1"))
  #Eliminamos todas las muestras (filas) que no contengan ningun valor de tipaje al hacer la exportacion de datos
  df <- subset(df, !apply(df[, -1], 1, function(row) all(is.na(row))))
  warning(paste(nrow(HLA_baja %>% anti_join(df, by = "Número") %>% select("Número")),"muestras no tienen datos de tipaje y fueron descartadas."))  #Aviso con el número de muestras que se han ignorado por no tener tipado ningun loci

  #La información importada desde Modulab puede contener anotaciones manuales que no son admitidas por algunos de los pasos de preprocesamiento
  #Los siguientes filtros eliminan todas las anotaciones manuales o erratas encontradas hasta el momento

  #Eliminamos todos los valores contenidos entre paréntesis en todas las columnas menos la primera (numero de muestra)
  #Asi eliminamos las anotaciones sobre tipaje serologico de algunas celdas
  df[, -1] <- lapply(df[, -1], function(x) gsub("\\([^\\)]+\\)", "", x))
  #Sustituye "*-" encontrado en algún caso para indicar homocigoto por "-"
  df[, -1] <- lapply(df[, -1], function(x) gsub("\\*-", "-", x))
  #Sustituye ":-" encontrado por ", -"
  df[, -1] <- lapply(df[, -1], function(x) gsub(":\\s*-", ", -", x))
  #Sustituye "Homocigoto" por ", -"
  df[, -1] <- lapply(df[, -1], function(x) gsub(" Homocigoto", ", -", x))
  #Sustituye todos los "." por "," para separar alelos
  df[, -1] <- lapply(df[, -1], function(x) gsub("\\.", ",", x))
  #Sustituye ", ," por ", -" ya que en algún caso el homocigoto se estaba señalando originalmente por "." en vez de "-" y este fue cambiado a "-" en la línea anterior
  df[, -1] <- lapply(df[, -1], function(x) gsub(",\\s*\\,", ", -", x))
  #Encuentra patrones correspondientes a alelos que esten separados por un espacio y los separa por ", "
  df[, -1] <- lapply(df[, -1], function(x) gsub("([0-9]{1,2})\\s+\\*", "\\1, \\*", x))
  #Encuentra patrones correspondientes a tipajes serologicos del tipo "*XX, *YY" y elimina todo lo que haya a continuación
  df[, -1] <- lapply(df[, -1], function(x) gsub("(\\*\\d{1,2}, \\*\\d{1,2}).*", "\\1", x))
  #Encuentra patrones XXY donde XX son 1 o 2 cifras e Y es una letra cualquier y elimina la letras, se eliminan así anotaciones de expresion alternativa de alelos (p.e. A*03:35N)
  df[, -1] <- lapply(df[, -1], function(x) gsub("([0-9]{1,2})[A-Za-z]", "\\1", x))
  #Renombramos las columnas para eliminar el termino "LOCUS " y usar la nomenclatura que usa hlapro::upscale_typings()
  names(df) <- c("Número", "A", "B", "C", "DRB1", "DRB345","DQB1", "DQA")

  #Modificamos el dataframe para obtener el formato final
  for (locus in colnames(df[, -1])) {
    locus1 <- paste0(locus, "1")
    locus2 <- paste0(locus, "2")
    new_locus <- c(locus1, locus2)
    df[[locus]] <- gsub("\\s", "", df[[locus]]) #Elimina los espacios en blanco en cada celda para evitar que el segundo alelo tenga un espacio (pe A* 01)
    df <- separate(df, col = locus, into = c(locus1, locus2), sep = "\\,")
    #Conversión de "-" que indican homocigotos
    #Si el formato de las celdas originales es correcto (pe. "*xx, *xx" o "*xx, -") no debería de ocurrir que una columna resultante sea NA y otra no lo sea
    #Por eso solo especificamos estas condiciones
    for (i in 1:nrow(df)) {
      if (is.na(df[i, locus1]) && is.na(df[i, locus2])){ #No hacer nada si ambas columnas son NA, columna original NA (locus no tipado)
      } else if (df[i, locus1] == "-") {
        df[i, locus1] <- df[i, locus2]
      } else if (df[i, locus2] == "-") {
        df[i, locus2] <- df[i, locus1]
      }
    }

    for(col in new_locus){
      #Sustituye "*0" por "*" solo si no contiene ":", es decir, si ya esta en alta resolución
      contains_colon <- grepl(":", df[[col]])
      df[[col]] <- ifelse(contains_colon, df[[col]], gsub("\\*0", "*", df[[col]]))

      df[[col]][!is.na(df[[col]])] <- paste0(locus, df[[col]][!is.na(df[[col]])])
    }
  }

  #Algunos loci pueden tener un alelo tipado en resolución intermedia y otros en resolución baja, esto imposibilita su uso para la imputación
  #Loop para rebajar la resolución de alelos de intermedia a baja
  df_serology <- df
  for (i in 1:nrow(df_serology)) {
    for (j in 2:ncol(df_serology)){
      if (grepl(":", df_serology[i, j])){
        df_serology[i,j] <- hlapro::get_serology(df[i,j])
      }
    }
  }

  #Algunos alelos serólogicos del locus C no son reconocidos correctamente por hlapro::upscale_typings() y devuelve un resultado nulo
  #Por esta razón es mejor eliminar estos alelos (Cw8, Cw11, Cw13, Cw16 y Cw17) antes de la imputación
  #Creamos un nueva columna en df para guardar los tipajes que se han obtenido sin tener en cuenta el locus C
  df$Cw_ignored <- FALSE
  #Se evalua si alguno de estos alelos se encuentra en las columnas C1 o C2, y en caso afirmativo se asigna NA en df_serology y TRUE al registro en df
  for (i in 1:nrow(df_serology)){
    if (any(c("C*8", "C*11", "C*13", "C*16", "C*17") %in% df_serology[i, "C1"] |
        c("C*8", "C*11", "C*13", "C*16", "C*17") %in% df_serology[i, "C2"])){
      df_serology[i, "C1"] <- NA
      df_serology[i, "C2"] <- NA
      df[i,"Cw_ignored"] <- TRUE
    } else {
      next
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

  # Iterar sobre las filas y columnas del dataframe
  for (i in 1:nrow(df_serology)) {
    for (j in 2:ncol(df_serology)) {

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

  #Añadimos a cada celda el prefijo de locus correspondiente (pe. A, B, C...)

  return(curated = list(df, df_serology))
}

vectorization_HLA <- function(df){
  low_list <- list()
  #Loop que itera las filas del dataframe df
  for (i in 1:nrow(df)){
    HLAlist <- list()

    #Loop anidado que itera por las celdas de la fila i, comprueba si cada celda contiene el simbolo : y si NO lo tiene lo añade a HLAlist
    #Ignoramos de esta manera los alelos que ya estan en resolución intermedia
    for(j in 2:ncol(df)){
      if (!grepl(":", df[i,j], fixed = TRUE)){
        cell <- df[i,j]
        #Todas las celdas (j) que cumplan la condición para la fila i se añaden a una lista temporal denominada HLAlist
        HLAlist <- c(HLAlist, cell)
      }
    }
    #La lista temporal se delista, se ignoran NAs y se concatena generando un string que se añade como nuevo elemento a la lista low_list
    x <- HLAlist %>%  unlist() %>%  na.omit() %>%  paste0(collapse = " ") %>% gsub("\\*", "", .)
    names(x) <- df[i,1]

    low_list <- c(low_list, x)
  }
  low_list <- low_list[nzchar(low_list)]

  # #Crea un nuevo df con los indices de las filas que no han sido añadidas a lowlist, aquellas que tienen todos
  # df2 <- df[which(!(df$Número %in% names(low_list))), ]

  return(low_list)
}

eplet_assignment <- function(hres_df, verified = "Yes"){
  #Evalua que el argumento verief tome los valores "Yes" solo para eplets verificados por anticuerpo o "All" que incluye a todos
  allowed_values <- c("Yes", "All")
  verify <- match.arg(verified, allowed_values)

  #Inicialización del data frame
  eplet_assignment <- data.frame(Número = hres_df$Número,
                                 setNames(data.frame(matrix(NA, nrow = nrow(hres_df), ncol = ncol(hres_df) - 1)), paste0(names(hres_df)[-1], "_eplets")))

  df_eplets <- hlapro::load_eplet_registry() # Cargar la base de datos de eplets de HLA eplet registry
  if (verified == "Yes"){
    df_eplets <- subset(df_eplets, df_eplets$confirmation == "Yes")
  } else {
    next
  }
  # Loop para asignar eplets a cada celda de eplet_assignment
  for (i in 1:nrow(hres_df)){
    for (j in 2:ncol(hres_df)){
      eplet_assignment[i, j] <- paste(hlapro::lookup_eplets(df_eplets, hres_df[i, j])[[1]], collapse = " ")
    }
  }
  #Convertimos todas las celdas vacias a NA
  eplet_assignment[eplet_assignment == ""] <- NA

  #Obtentener los eplets únicos entre clase I o clase II o todos
  #Obtener los índices de las columnas que corresponden a clase I o clase II
  HLA_I_indices <- which(grepl("^[ABC]", names(eplet_assignment)))
  HLA_II_indices <- which(grepl("^D", names(eplet_assignment)))
  HLA_all_indices <- which(grepl("^[ABCD]", names(eplet_assignment)))

  #Función para dividir, obtener valores únicos y colapsar
  split_unique_collapse <- function(cell) {
    all_values <- paste(cell, collapse = " ")
    unique_vals <- unique(strsplit(all_values, " ")[[1]])
    non_na_unique_vals <- unique_vals[unique_vals != "NA"]
    collapsed_str <- paste(non_na_unique_vals, collapse = " ")
    return(collapsed_str)
  }

  #Aplicar la función split_unique_collapse a las columnas HLA_I_eplets y HLA_II_eplets
  eplet_assignment$HLA_I_eplets <- apply(eplet_assignment[HLA_I_indices], 1, split_unique_collapse)
  eplet_assignment$HLA_II_eplets <- apply(eplet_assignment[HLA_II_indices], 1, split_unique_collapse)
  eplet_assignment$HLA_all_eplets <- apply(eplet_assignment[HLA_all_indices], 1, split_unique_collapse)

  #Fusionamos las tablas hres_df y eplet_assignment por la columna Número
  HLA_final <- merge(hres_df, eplet_assignment, by = "Número", all = TRUE)
  return(HLA_final)
}

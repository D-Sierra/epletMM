
# epletMM <a href="https://D-sierra.github.io/epletMM/"><img src="man/figures/hex.png" align="right" height="250" alt="hlapro website" /></a>

<!-- badges: start -->
  [![R-CMD-check](https://github.com/D-Sierra/epletMM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/D-Sierra/epletMM/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

Esta librería proporciona una serie de herramientas para trabajar con datos de tipajes HLA.

Las funciones contenidas en este paquete están pensadas para realizar las siguientes operaciones:

-   Preprocesar un archivo excel que contenga información de tipajes HLA
-   Automatizar la imputación a alta resolución de los loci A, B, C, DRB1, DR345, DQB1 y DQA1.
-   Automatizar la asignación de eplets de tipajes HLA de alta resolución
-   Realizar el crossmatch virtual entre parejas donante/receptor y calcular el missmatch (MM) a nivel de eplets

## Información importante

1.  Actualmente este paquete depende directamente de las funciones `upscale_typings()`, `get_serology()`, `load_eplet_registry()` y `lookup_eplets()` de la librería [hlapro](https://github.com/lcreteig/hlapro/tree/main).

    -   `upscale_typings()` requiere descargar previamente la información sobre la frecuencia de haplotipos proporcionada por NMDP de la siguiente [web](https://frequency.nmdp.org/). Una vez registrado veremos un listado de posibles archivos para descargar, es necesario seleccionar el formato .xlsx que contiene la información para los loci A, B, C, DRB1, DR345, DQB1 (HLA\-A\~C\~B\~DRB3\-4-5\~DRB1\~DQB1).
    -   `load_eplet_registry()` realiza un scrapping de la base de datos HLA Eplet Registry y descarga una tabla con el listado de eplets descritos hasta el momento si esta no existe ya de forma local. Cualquier actualización en la web HLA Eplet Registry que altere su estructura puede imposibilitar la descargar de la información de eplets necesaria en la primera ejecución de la función en un nuevo equipo. Para estos casos este paquete contiene la última versión actualizada del registro de eplets a fecha 03/03/2024, este archivo (eplets.rds) debe copiarse a la carpeta cache de la librería hlapro, normalmente localizada el directorio \~/AppData/Local/R/cache/R/hlapro/ dentro del usuario correspondiente.

2.  En su versión actual (0.1.3) las funciones disponibles para el preprocesamiento de datos con tipajes HLA son `modulab_filter()` y `preimputation()`. Estas funciones trabajan en serie, la primera toma como argumento un dataframe obtenido a partir de una exportación de datos del LIS Modulab. La tabla cargada en este dataframe puede contener cualquier número de columnas, pero debe de cumplir ciertos requisitos:

    -   Debe de existir una columna con un identificador único para cada individuo y su nombre debe de ser "Número", por defecto una exportación de Modulab debería de proporcionar el número de historía clínica en una columna usando este nombre.
    -   No es necesario que los tipajes contengan información sobre todos los loci, ni que la tabla contenga las columnas para todos ellos. Sin embargo, los nombres de las columnas de los loci que deseemos incluir en el análisis deben de ser los siguiente: LOCUS A, LOCUS B, LOCUS C, LOCUS DPB1, LOCUS DQA1, LOCUS DQB1, LOCUS DRB1.

3.  La función `preimputation()`establece ciertas condiciones para realizar la imputación:

    -  Hay ciertos alelos serológicos del locus C que no son admitidos por `hlapro::upscale_typings()` (Cw8, Cw11, Cw13, Cw16 y Cw17) si un genotipo contiene 1 o 2 de estos alelos la información del locus C es eliminada para evitar que falle la imputación de alta resolución. Algunos de estos alelos serológicos tampoco son admitidos por Haplostats.
    -  La falta de información suficiente de tipaje hace que la imputación no sea precisa y sus resultados no sean fiables, por ello se establece que las muestras tengan tipado por lo menos los loci A, B y DRB1 para realizar la imputación. Las muestras que no cumplan este requisito serán eliminadas del análisis.
    -  La función añade columnas adicionales al dataframe procesado indicando en que muestras la imputación no se ha realizado por no alcanzar el mínimo de información disponible y cuales se han hecho ignorando el tipaje del locus C.
  
## Esquema de flujo de trabajo

<a href="https://D-sierra.github.io/epletMM/"><img src="man/figures/Pipeline de análisis.png" height="450" alt="pipeline" /></a>
## Instalación

```{r}
# install.packages("devtools")
devtools::install_github("D-Sierra/epletMM")
```

## Uso

```{r}
library(epletMM)
```

### Preprocesamiento de los datos

El preprocesamiento dependerá del archivo de partida que vayamos a utilizar. En su versión actual (0.1.3) esta librería contiene funciones para el preprocesamiento de archivo procedentes de dos fuentes distintas:

1.  `modulab_filter()`: esta función admite archivos de tipajes exportados desde Modulab. Estas tablas deben de tener siempre una tabla con IDs únicos por individuo llamada "Número", además los nombres de las columnas con los tipajes deben de ser los siguientes: "LOCUS A", "LOCUS B", "LOCUS C", "LOCUS DRB1", "LOCUS DRB345","LOCUS DQB1", "LOCUS DQA1". No es necesario que todas las columnas correspondientes a los loci estén presentes en los datos originales, aquellas que no lo esten serán creadas, asignadas el valor **NA** y más adelante serán imputadas.

2.  `nefro_filter()`:esta función admite archivos de tipajes exportados desde la base de datos de nefrología. Los datos importados de esta manera estan pensados para realizar el crossmatch entre pacientes y receptores, por lo que deberán tener el tipaje de ambos. La función requiere que la tabla original contenga las columnas "NumTx" con identificador del trasplante, y "NumHistoria" con el ID único de cada receptor. Se considera que en los datos originales existe una única fila por cada par donante/receptor con las siguientes columnas de tipaje: PacHLA1", "PacHLA2", "PacHLB1", "PacHLB2", "PacHLC1", "PacHLC2", "PacHLDr1", "PacHLDr2", "PacHLDQ1", "PacHLDQ2", "PacHLDQA1", "PacHLDQA2", "Don HLA A1", "Don HLA A2", "Don HLA B1", "Don HLA B2", "Don HLA C1", "Don HLA C2", "Don HLA Dr1", "Don HLA Dr2", "Don HLA DQ1", "Don HLA DQ2", "Don HLA DQA1", "Don HLA DQA2". Como en el caso anterior, cualquier columna de tipaje que no exista será creada, poblada con el valor **NA** e imputada más adelante.

```{r}
#Procedimiento para datos obtenidos de exportación de modulab
HLA_baja_modulab <- read_xls("HLA_baja_modulab.xls")
filtered_df_modulab <- modulab_filter(HLA_baja_modulab)
processed_data <- preimputation(filtered_df_modulab)

#Procedimiento para datos obtenidos de exportación de la base de datos de nefro
HLA_baja_nefro <- read_xls("HLA_baja_nefro.xls")
filtered_df_nefro <- nefro_filter(HLA_baja_nefro)
processed_data <- preimputation(filtered_df_nefro)

mixed_reso <- processed_data[[1]]
head(mixed_reso)
# # A tibble: 6 × 15
#   Número   A1      A2      B1      B2      C1      C2      DRB11      DRB12      DRB3451 DRB3452 DQB11      DQB12      DQA1  DQA2 
#   <chr>    <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>      <chr>      <chr>   <chr>   <chr>      <chr>      <chr> <chr>
# 1 40145308 A*26    A*31    B*44    B*57    C*5     C*7     DRB1*4     DRB1*7     NA      NA      DQB1*03:01 DQB1*03:02 NA    NA   
# 2 40145307 A*1     A*32    B*40    B*57    C*2     C*6     DRB1*7     DRB1*11    NA      NA      DQB1*03:01 DQB1*03:03 NA    NA   
# 3 40145305 A*1     A*2     B*40    B*41    C*7     C*15    DRB1*13    DRB1*15    NA      NA      DQB1*03:01 DQB1*6     NA    NA   
# 4 10205567 A*01:01 A*02:01 B*08:01 B*50:01 C*06:02 C*07:01 DRB1*03:01 DRB1*07:01 NA      NA      DQB1*02:01 DQB1*02:01 NA    NA   

low_reso <- processed_data[[2]]
head(low_reso)
# # A tibble: 6 × 15
#   Número   A1    A2    B1    B2    C1    C2    DRB11 DRB12 DRB3451 DRB3452 DQB11 DQB12 DQA1  DQA2 
#   <chr>    <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr>   <chr>   <chr> <chr> <chr> <chr>
# 1 40145308 A*26  A*31  B*44  B*57  Cw*5  Cw*7  DR*4  DR*7  NA      NA      DQ*7  DQ*8  NA    NA   
# 2 40145307 A*1   A*32  B*40  B*57  Cw*2  Cw*6  DR*7  DR*11 NA      NA      DQ*7  DQ*9  NA    NA   
# 3 40145305 A*1   A*2   B*40  B*41  Cw*7  Cw*15 DR*13 DR*15 NA      NA      DQ*7  DQ*6  NA    NA   
# 4 10205567 A*1   A*2   B*8   B*50  Cw*6  Cw*7  DR*17 DR*7  NA      NA      DQ*2  DQ*2  NA    NA 
```

La función `preimputation()` devuelve una lista que contiene tres dataframes distintos, para los siguientes pasos es recomendable asignar cada uno de los dos primeros a una nueva variable. En este caso:

-   processed_data[[1]]: es un dataframe que contiene la información del dataframe de partida manteniendo la resolución original para los tipajes HLA. Este dataframe será el actualizado con los nuevos tipajes en alta y la información de eplets.
-   processed_data[[2]]: tiene la misma estructura, pero todos los tipajes en resolución intermedia o alta han sido transformados en los equivalentes serologicos usando `hlapro:get_serology()`. Este dataframe será utilizado para realizar la imputación de alta resolución, ya que la función `hlapro::upscale_typings()`, al contrario que [HaploStats](https://www.haplostats.org/) no tiene en cuenta aquellos alelos proporcionados ya en alta resolución para realizar la imputación.
-   processed_data[[3]]: almacena la información del control de calidad.
  
**CONTROL DE CALIDAD**

`preimputation()` ejecuta un control de calidad de los tipajes de baja resolución contenidos en processed_data[[2]], esto incluye aquellos filtrados desde el df original y los que han sidos imputados mediante `hlapro::getserology()` a partir de datos de resolución alta/media. La función revisa para cada columna correspondiente a un alelo que todas las filas de processed_data[[2]] contengan datos con estructura correcta (p.e. alelos HLA-A deben de tener estructura "^A\\*\\d{1,2}$", que indica que deben de seguir un patrón exacto *A\*XX* donde XX son uno o dos caracteres numericos.
Si todas las columnas pasan correctamente el test se mostrará un mensaje indicándolo, en caso contrario se imprimirán mensajes de error que informaran que columnas no han pasado el test en su conjunto, y además se mostrará el valor o valores que no cumplen las condiciones de estructura y el número de ID ("Número") asociado para que pueda ser revisado el preprocesamiento.

```{r}
typings_vector <- vectorization_HLA(low_reso)

vec[1:2]
# $`40145308`
# [1] "A26 A31 B44 B57 Cw5 Cw7 DR4 DR7 DQ7 DQ8"

# $`40145307`
# [1] "A1 A32 B40 B57 Cw2 Cw6 DR7 DR11 DQ7 DQ9"
```

HLA_vec convierte el dataframe de tipajes en una lista de vectores formado por un elemento por cada individuo, estos elementos contienen el tipaje completo concatenado en una única cadena de texto y están nombrados según el ID contenido en la columna "Número" del dataframe original. Este paso es necesario ya que este será el formato de tipajes que reconocerá la función `hlapro::upscale_typings()` utilizada por `bulk_upscaling()`.

### Imputación HLA alta resolución

En primer lugar, debemos de conocer la localización del archivo que contiene la información de haplotipos y frecuencias descargado de NMDP, y asegurarnos de que su nombre sea "A\~C\~B\~DRB3-4-5\~DRB1\~DQB1.xlsx". Por defecto la función `bulk_upscaling()` marcará como ruta del archivo la correspondiente a la carpeta sobre la que estemos trabajando `getwd()`. Si no tuviesemos una copia de dicho archivo en esa dirección, deberemos de especificarlo mediante el argumento **haplotypes_path**.

Esta función requiere otros dos argumentos:

-   low_list: lista de vectores de baja resolución para realizar la imputación.
-   df: dataframe con información de tipajes obtenida durante el preprocesamiento. Por defecto la función bulk_upscaling sobreescribe los tipajes de baja resolución con la imputación de alta resolución, pero conserva todos los tipajes que ya consten en el dataframe proporcionado con alta resolución. Por esta razón, si queremos conservar la información en alta resolución que ya tuviesemos debemos de seleccionar el df processed_data[[1]] obtenido incialmente. Por el contrario, si queremos sobreescribir toda la información, debemos de utilizar el df processed_data[[2]], *hay que tener en cuenta que esta última opción puede producir el cambio de alelos ya tipados en alta resolución por otros distintos obtenidos mediante la imputación*.

```{r}
hres_df <- bulk_upscaling(low_list = "low_list", df = "subdf1", haplotypes_path = getwd())
```

`bulk_upscaling()` realiza la imputación en dos pasos, en primer lugar se sirve de la función `hlapro::upscale_typings()` para realizar la imputación en alta resolución de los loci A, B, C, DRB1, DRB345 y DQB1. En segundo lugar crea una tabla de equivalencias y asociación de haplotipos para los loci DRB1, DQB1 y DQA1 y realiza la impitación de alta resolución del locus DQA1.

La función devuelve dos tipos de mensaje de advertencia que no impiden que el código termine su ejecución:

-   Al aplicar la función `hlapro::upscale_typings` puede no devolver haplotipos compatibles con los datos proporcionados. Esto es más probable con haplotipos raros o si algún alelo, especialmente DQ, ha sido previamente pasado a resolución serológica. En estos casos el programa elimina el tipaje correspondiente a DQB1 y repite la imputación. Si la imputación obtiene algún resultado este se guarda y se devuelve un mensaje de advertencia que indica el ID del individuo cuyo tipaje ha sido imputado ignorando el locus DQB1.
-   En el caso de que al eliminar el locus DQB1 se siga sin obtener ningún resultado tras la imputación se devuelve un mensaje que indica el ID del individuo cuyos haplotipos no hayan podido ser imputados.

### Asignación de eplets

```{r}
hres_eplets <- eplet_assigment(hres_df, verified = "Yes")
```

El dataframe proporcionado debe de contener información HLA de resolución intermedia/alta (p.e. A\*01:01), cualquier alelo en resolución baja no se tendrá en cuenta para el cálculo de eplets. La función devolverá un nuevo df que contendrá la información original más una nueva columna por cada alelo en el dataframe proporcionado, y que contendrá la colección de eplets asociados al alelo correspondiente. Además contendrá tres nuevas columnas que mostrarán todos los eplets únicos para los loci HLA de clase I y de clase II, así como todos los eplets únicos totales para cada uno de los individuos.

El argumento verified indica si deseamos tener en cuenta solo aquellos eplets que HLA Eplet Registry asegure que han sido confirmados por anticuerpos (verified = "Yes") o todos los eplets del registro incluidos los teóricos (verified = "All")

### Virtual eplet crossmatch

```{r}
eplets_MM <- eplet_crossmatch(hres_eplets)
```

El dataframe de entrada debe de ser obtenido mediante la función `eplet_assigment()` o tener la misma estructura de datos, con una columna por cada alelo A, B, C, DRB1, DRB345, DQB1 y DQA1, con los eplets únicos correspondientes a cada uno de ellos. La función devuelve por cada loci (A, B, C... ) dos columnas, la primera contiene el número de eplets presentes en el donante que no se encuentran en el receptor, y la segunda el listado de esos eplets. Esto mismo se repite para el total de eplets en los loci HLA de clase I, de clase II y para la totalidad del tipaje. El dataframe de salida contiene además una puntuación de riesgo del MM calculada con la fórmula de la librería `hlaR` y que utiliza el número de MM de eplets en los alelos DRB1, DRB345, DQB1 y DQA1.

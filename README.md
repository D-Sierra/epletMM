
# epletMM <a href="https://D-sierra.github.io/epletMM/"><img src="man/figures/hex.png" align="right" height="250" alt="hlapro website" /></a>

Esta librería proporciona una serie de herramientas para trabajar con datos de tipajes HLA.

Las funciones contenidas en este paquete están pensadas para realizar las siguientes operaciones:

-   Preprocesar un archivo excel que contenga información de tipajes HLA
-   Automatizar la imputación a alta resolución de los loci A, B, C, DRB1, DR345 y DQB
-   Automatizar la asignación de eplets de tipajes HLA de alta resolución
-   Realizar el crossmatch virtual entre parejas donante/receptor y calcular el missmatch (MM) a nivel de eplets

## Información importante

1.  Actualmente este paquete depende directamente de las funciones `upscale_typings()`, `get_serology()`, `load_eplet_registry()` y `lookup_eplets()` de la librería [hlapro](https://github.com/lcreteig/hlapro/tree/main).

    -   `upscale_typings()` requiere descargar previamente la información sobre la frecuencia de haplotipos proporcionada por NMDP de la siguiente [web](https://frequency.nmdp.org/). Una vez registrado veremos un listado de posibles archivos para descargar, es necesario seleccionar el formato .xlsx que contiene la información para los loci A, B, C, DRB1, DR345, DQB1 (HLA\-A\~C\~B\~DRB3\-4-5\~DRB1\~DQB1).
    -   `load_eplet_registry()` realiza un scrapping de la base de datos HLA Eplet Registry y descarga una tabla con el listado de eplets descritos hasta el momento si esta no existe ya de forma local. Cualquier actualización en la web HLA Eplet Registry que altere su estructura puede imposibilitar la descargar de la información de eplets necesaria en la primera ejecución de la función en un nuevo equipo. Para estos casos este paquete contiene la última versión actualizada del registro de eplets a fecha 03/03/2024, este archivo (eplets.rds) debe copiarse a la carpeta cache de la librería hlapro, normalmente localizada el directorio \~/AppData/Local/R/cache/R/hlapro/ dentro del usuario correspondiente.

2.  En su versión actual (0.1.0) la única función disponible para el preprocesamiento de datos con tipajes HLA es `modulab_HLA_cleaner()`, esta función toma como argumento un dataframe obtenido a partir de una exportación de datos del LIS Modulab. La tabla cargada en este dataframe puede contener cualquier número de columnas, pero debe de cumplir ciertos requisitos:

    -   Debe de existir una columna con un identificador único para cada individuo y su nombre debe de ser "Número", por defecto una exportación de Modulab debería de proporcionar el número de historía clínica en una columna usando este nombre.
    -   No es necesario que los tipajes contengan información sobre todos los loci, ni que la tabla contenga las columnas para todos ellos. Sin embargo, los nombres de las columnas de los loci que deseemos incluir en el análisis deben de ser los siguiente: LOCUS A, LOCUS B, LOCUS C, LOCUS DPB1, LOCUS DQA1, LOCUS DQB1, LOCUS DRB1.

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

En primer lugar nos aseguramos que nuestro archivo (.xlsx o .csv) cumple los requisitos necesarios. Cargamos el archivo lo modificamos para prepararlo para la imputación.

```{r}
HLA_baja <- read_xls("HLA_baja.xls")

processed_data <- modulab_HLA_cleaner(HLA_baja)

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

La función `modulab_HLA_cleaner()` devuelve una lista que contiene dos dataframes distintas, para los siguientes pasos es recomendable asignar cada una de ellas a una nueva variable. En este caso:

-   processed_data[[1]]: es un dataframe que contiene la información del dataframe de partida manteniendo la resolución original para los tipajes HLA. Este dataframe será el actualizado con los nuevos tipajes en alta y la información de eplets.
-   processed_data[[2]]: tiene la misma estructura, pero todos los tipajes en resolución intermedia o alta han sido transformados en los equivalentes serologicos usando `hlapro:get_serology()`. Este dataframe será utilizado para realizar la imputación de alta resolución, ya que la función `hlapro::upscale_typings()`, al contrario que [HaploStats](https://www.haplostats.org/) no tiene en cuenta aquellos alelos proporcionados ya en alta resolución para realizar la imputación.

```{r}
typings_vector <- HLA_vec(low_reso)

vec[1:2]
# $`40145308`
# [1] "A26 A31 B44 B57 Cw5 Cw7 DR4 DR7 DQ7 DQ8"

# $`40145307`
# [1] "A1 A32 B40 B57 Cw2 Cw6 DR7 DR11 DQ7 DQ9"
```

HLA_vec convierte el dataframe de tipajes en una colección una lista de vectores formado por un elemento por cada individuo, estos elementos contienen el tipaje completo concatenado en una única cadena de texto y están nombrados según el ID contenido en la columna "Número" del dataframe original. Este paso es necesario ya que este será el formato de tipajes que reconocerá la función `hlapro::upscale_typings()` utilizada por `bulk_upscaling()`.

### Imputación HLA alta resolución

En primer lugar, debemos de conocer la localización del archivo que contiene la información de haplotipos y frecuencias descargado de NMDP, y asegurarnos de que su nombre sea "A~C~B~DRB3-4-5~DRB1\~DQB1.xlsx". Por defecto la función `bulk_upscaling()` marcará como ruta del archivo la correspondiente a la carpeta sobre la que estemos trabajando `getwd()`. Si no tuviesemos una copia de dicho archivo en esa dirección, deberemos de especificarlo meidante el argumento **haplotypes_path**.

Esta función requiere otros dos argumentos:

-   low_list: lista de vectores de baja resolución para realizar la imputación.
-   df: dataframe con información de tipajes obtenida durante el preprocesamiento. Por defecto la función bulk_upscaling sobreescribe los tipajes de caja resolución con la imputación de alta resolución, pero conserva todos los tipajes que ya consten en el dataframe proporcionado con alta resolución. Por esta razón, si queremos conservar la información en alta resolución que ya tuviesemos debemos de seleccionar el df processed_data[[1]] obtenido incialmente. Por el contrario, si queremos sobreescribir toda la información, debemos de utilizar el df processed_data[[2]], *hay que tener en cuenta que esta última opción puede producir el cambio de alelos ya tipados en alta resolución por otros distintos obtenidos mediante la imputación*.

```{r}
hres_df <- bulk_upscaling(low_list = "low_list", df = "subdf1", haplotypes_path = getwd())
```

La función `bulk_upscaling()` devuelve dos tipos de mensaje de advertencia que no impiden que el código se ejecute por completo:

-   Al aplicar la función `hlapro::upscale_typings` el tipaje de algunos individuos puede no devolver resultados si contiene información de muchos loci o si alguno de ellos, especialmente DQ, ha sido previamente pasado a resolución serológica. En estos casos el programa elimina el tipaje correspondiente a DQB1 y repite la imputación. Si la imputación obtiene algún resultado este se guarda y se devuelve un mensaje de advertencia que indica el ID del individuo cuyo tipaje ha sido imputado ignorando el locus DQB1.
-   En el caso de que al eliminar el locus DQB1 se siga sin obtener ningún resultado tras la imputación se devuelve un mensaje que indica el ID del individuo cuyos haplotipos no hayan podido ser imputados.

### Asignación de eplets

```{r}
hres_eplets <- assign_eplets(hres_df, verified = "Yes")
```

El dataframe proporcionado debe de contener información HLA de resolución intermedia/alta (p.e. A\*01:01), cualquier alelo en resolución baja no se tendrá en cuenta para el cálculo de eplets. La función devolverá un nuevo df que contendrá la información original más una nueva cada columna por cada alelo en el dataframe proporcionado, y que contendrá la colección de eplets asociados al alelo correspondiente. Además contendrá tres nuevas columnas que mostrarán todos los eplets únicos para los loci HLA de clase I y de clase II, así como todos los eplets únicos totales para cada uno de los individuos.

El argumento verified indica si deseamos tener en cuenta solo aquellos eplets que HLA Eplet Registry asegure que han sido confirmados por anticuerpos (verified = "Yes") o todos los eplets del registro incluidos los teóricos (verified = "All")

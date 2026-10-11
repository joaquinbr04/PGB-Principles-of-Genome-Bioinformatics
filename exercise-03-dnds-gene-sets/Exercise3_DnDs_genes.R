# Exercise 3: comparison of dN/dS ratios among gene sets
# Course: Principles of Genome Bioinformatics


# =============================================================================
# 1. PAQUETES NECESARIOS 
# =============================================================================

# install.packages(c("dplyr", "ggplot2", "viridis"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install(c("biomaRt", "GO.db", "AnnotationDbi"))

# Cargamos los paquetes necesarios para el análisis.

library(biomaRt)
library(GO.db)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(viridis)


# =============================================================================
# 2. LEER EL ARCHIVO CON LOS VALORES DE dN, dS Y dN/dS
# =============================================================================

# HumanDnDsW.txt debe estar en la carpeta principal del proyecto.
# Usamos una ruta relativa para que el script funcione también cuando
# otra persona descargue el repositorio en un ordenador diferente.

input_file <- "HumanDnDsW.txt"

if (!file.exists(input_file)) {
  stop(
    "No se encuentra HumanDnDsW.txt. ",
    "Comprueba que el archivo esté en la carpeta del proyecto."
  )
}

dnds_data <- read.delim(
  file = input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Columnas que esperamos encontrar en el archivo.
# Conservamos Homo_sapiens.dw porque forma parte del formato original,
# aunque no se utilice directamente en los pasos posteriores.

required_columns <- c(
  "GeneID",
  "Homo_sapiens.dN",
  "Homo_sapiens.dS",
  "Homo_sapiens.dw",
  "Homo_sapiens.w"
)

missing_columns <- setdiff(
  required_columns,
  colnames(dnds_data)
)

if (length(missing_columns) > 0) {
  stop(
    "Faltan las siguientes columnas en HumanDnDsW.txt: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Eliminamos posibles versiones añadidas a los identificadores Ensembl.
# Por ejemplo, ENSG00000123456.12 pasa a ENSG00000123456.

dnds_data$GeneID <- sub(
  pattern = "\\.[0-9]+$",
  replacement = "",
  x = dnds_data$GeneID
)


# =============================================================================
# 3. CONECTARSE A ENSEMBL BIOMART
# =============================================================================

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  GRCh = 37
)


# =============================================================================
# 4. OBTENER INFORMACION DE LOS GENES
# =============================================================================

gene_info <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol",
    "chromosome_name",
    "start_position",
    "end_position"
  ),
  filters = "ensembl_gene_id",
  values = unique(dnds_data$GeneID),
  mart = ensembl
)

# Queremos una sola fila por identificador Ensembl antes de hacer el join.

gene_info <- gene_info %>%
  distinct(
    ensembl_gene_id,
    .keep_all = TRUE
  )

# GeneID en dnds_data corresponde a ensembl_gene_id en gene_info.

dnds_data <- dnds_data %>%
  left_join(
    gene_info,
    by = c("GeneID" = "ensembl_gene_id")
  )

# Contamos los genes que no han recibido un simbolo HGNC.

genes_without_hgnc <- sum(
  is.na(dnds_data$hgnc_symbol) |
    dnds_data$hgnc_symbol == ""
)

message("Genes without an HGNC symbol: ", genes_without_hgnc)


# =============================================================================
# 5. OBTENER LOS TERMINOS DE GENE ONTOLOGY
# =============================================================================

# GOTERM contiene los terminos de Gene Ontology. Lo convertimos en una lista
# y extraemos el nombre de cada termino, conservando los identificadores GO
# como nombres del vector.

go_objects <- as.list(GOTERM)

go_terms <- vapply(
  go_objects,
  function(x) slot(x, "Term"),
  FUN.VALUE = character(1)
)


# =============================================================================
# 6. BUSCAR TERMINOS RELACIONADOS CON COMPLEMENT
# =============================================================================

complement_terms <- grep(
  pattern = "complement",
  x = go_terms,
  ignore.case = TRUE,
  value = TRUE
)

complement_go_ids <- names(complement_terms)

complement_terms_table <- data.frame(
  GO_ID = complement_go_ids,
  GO_term = unname(complement_terms)
)

# print(complement_terms_table)


# =============================================================================
# 7. BUSCAR TERMINOS RELACIONADOS CON SPLICING
# =============================================================================

splicing_terms <- grep(
  pattern = "splicing",
  x = go_terms,
  ignore.case = TRUE,
  value = TRUE
)

splicing_go_ids <- names(splicing_terms)

splicing_terms_table <- data.frame(
  GO_ID = splicing_go_ids,
  GO_term = unname(splicing_terms)
)

# print(splicing_terms_table)


# =============================================================================
# 8. BUSCAR TERMINOS RELACIONADOS CON NEUROGENESIS
# =============================================================================

neurogenesis_terms <- grep(
  pattern = "neurogenesis",
  x = go_terms,
  ignore.case = TRUE,
  value = TRUE
)

neurogenesis_go_ids <- names(neurogenesis_terms)

neurogenesis_terms_table <- data.frame(
  GO_ID = neurogenesis_go_ids,
  GO_term = unname(neurogenesis_terms)
)

# print(neurogenesis_terms_table)


# =============================================================================
# 9. CREAR UNA FUNCION PARA OBTENER GENES A PARTIR DE TERMINOS GO
# =============================================================================

get_genes_from_go <- function(go_ids, mart) {

  # La funcion consulta en Ensembl los genes humanos asociados a los
  # identificadores GO, elimina identificadores vacios y duplicados,
  # y devuelve una tabla.

  if (length(go_ids) == 0) {
    stop("No se han encontrado identificadores GO.")
  }

  gene_table <- getBM(
    attributes = c(
      "ensembl_gene_id",
      "hgnc_symbol",
      "go_id"
    ),
    filters = "go",
    values = go_ids,
    mart = mart
  )

  gene_table <- gene_table %>%
    filter(ensembl_gene_id != "") %>%
    distinct(
      ensembl_gene_id,
      .keep_all = TRUE
    )

  return(gene_table)
}


# =============================================================================
# 10. OBTENER LOS GENES DE CADA CATEGORIA
# =============================================================================

complement_genes <- get_genes_from_go(
  go_ids = complement_go_ids,
  mart = ensembl
)

splicing_genes <- get_genes_from_go(
  go_ids = splicing_go_ids,
  mart = ensembl
)

neurogenesis_genes <- get_genes_from_go(
  go_ids = neurogenesis_go_ids,
  mart = ensembl
)

# Estas comprobaciones son opcionales:

# length(unique(complement_genes$ensembl_gene_id))
# length(unique(splicing_genes$ensembl_gene_id))
# length(unique(neurogenesis_genes$ensembl_gene_id))


# =============================================================================
# 11. CREAR UNA FUNCION PARA CONSTRUIR CADA GENE SET
# =============================================================================

make_gene_set <- function(data, gene_ids, label) {

  # Selecciona de la tabla principal los genes incluidos en gene_ids,
  # conserva sus valores evolutivos y añade la etiqueta del grupo.

  data %>%
    filter(GeneID %in% unique(gene_ids)) %>%
    select(
      GeneID,
      Homo_sapiens.dN,
      Homo_sapiens.dS,
      Homo_sapiens.w,
      hgnc_symbol
    ) %>%
    mutate(Geneset = label)
}


# =============================================================================
# 12. CREAR LOS TRES GENE SETS
# =============================================================================

df_complement <- make_gene_set(
  data = dnds_data,
  gene_ids = complement_genes$ensembl_gene_id,
  label = "Complement"
)

df_splicing <- make_gene_set(
  data = dnds_data,
  gene_ids = splicing_genes$ensembl_gene_id,
  label = "Splicing"
)

df_neurogenesis <- make_gene_set(
  data = dnds_data,
  gene_ids = neurogenesis_genes$ensembl_gene_id,
  label = "Neurogenesis"
)


# =============================================================================
# 13. UNIR LOS TRES GENE SETS
# =============================================================================

gene_sets <- bind_rows(
  df_complement,
  df_splicing,
  df_neurogenesis
)

# Convertimos las columnas evolutivas en variables numericas.

gene_sets <- gene_sets %>%
  mutate(
    across(
      c(
        Homo_sapiens.dN,
        Homo_sapiens.dS,
        Homo_sapiens.w
      ),
      as.numeric
    )
  ) %>%
  mutate(
    Geneset = factor(
      Geneset,
      levels = c(
        "Complement",
        "Splicing",
        "Neurogenesis"
      )
    )
  )


# =============================================================================
# 14. ELIMINAR VALORES DE OMEGA AUSENTES O INFINITOS
# =============================================================================

# Eliminamos valores no finitos de omega (NA, NaN, Inf o -Inf) porque no
# pueden representarse correctamente. Pueden surgir por datos ausentes,
# errores de conversion o por calcular omega = dN/dS cuando dS = 0.

plot_data <- gene_sets %>%
  filter(is.finite(Homo_sapiens.w))

# Esta comprobación es opcional:

# genes_per_set <- plot_data %>%
#   count(
#     Geneset,
#     name = "number_of_genes"
#   )
#
# print(genes_per_set)


# =============================================================================
# 15. REPRESENTAR LOS RESULTADOS
# =============================================================================

dnds_plot <- ggplot(
  data = plot_data,
  mapping = aes(
    x = Geneset,
    y = Homo_sapiens.w,
    fill = Geneset
  )
) +
  geom_violin(
    width = 0.9,
    trim = FALSE,
    alpha = 0.85,
    na.rm = TRUE
  ) +
  geom_boxplot(
    width = 0.12,
    color = "grey20",
    fill = "white",
    alpha = 0.70,
    outlier.shape = NA,
    na.rm = TRUE
  ) +
  scale_fill_viridis_d(
    option = "D",
    end = 0.85
  ) +
  # Este limite amplia visualmente el intervalo de omega entre 0 y 1.
  # Los valores superiores a 1 permanecen en plot_data y en los resumenes,
  # pero no aparecen en esta vista de la grafica.
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  labs(
    title = "dN/dS ratios among human gene sets",
    x = "Gene set",
    y = expression(omega == d[N] / d[S])
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      size = 12,
      face = "bold"
    ),
    panel.grid.major.x = element_blank()
  )

print(dnds_plot)


# =============================================================================
# 16. GUARDAR LA GRAFICA
# =============================================================================

# Descomenta estas lineas si quieres guardar la figura.

# ggsave(
#   filename = "DnDs_gene_sets.png",
#   plot = dnds_plot,
#   width = 7,
#   height = 5,
#   dpi = 300
# )


# =============================================================================
# 17. CALCULAR UN RESUMEN NUMERICO
# =============================================================================

summary_table <- plot_data %>%
  group_by(Geneset) %>%
  summarise(
    n = n_distinct(GeneID),
    median_omega = median(
      Homo_sapiens.w,
      na.rm = TRUE
    ),
    mean_omega = mean(
      Homo_sapiens.w,
      na.rm = TRUE
    ),
    IQR_omega = IQR(
      Homo_sapiens.w,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(summary_table)

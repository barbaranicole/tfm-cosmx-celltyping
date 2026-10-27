# Identificación de tipos celulares en transcriptómica espacial: comparativa de métodos de machine learning

Código del Trabajo de Fin de Máster **"Estudio comparativo de métodos de machine learning para la identificación del tipo celular mediante single-cell RNAseq con resolución espacial"**.

- **Autora:** Bárbara Bernal Gómez
- **Máster:** Máster en Bioestadística, Universidad Complutense de Madrid
- **Tutores:** Adrián Portillo, Melina Peressini, Silvia Pineda Sanjuan
- **Memoria:** [Repositorio Docta UCM](https://docta.ucm.es/entities/publication/0bf5fe30-981f-46d3-9db2-77e20d72b17b)

## Descripción

Se comparan dos métodos de identificación de tipos celulares sobre datos de transcriptómica espacial a resolución de célula única (CosMx SMI, cáncer de pulmón no microcítico). Se evalúa la robustez de cada método frente a la degradación de los datos, comparando los resultados sobre los datos originales y sobre una versión degradada (`thin_75`).

| Método | Tipo | Paquete |
|---|---|---|
| Leiden | No supervisado | Seurat |
| InSituType | Semisupervisado | InSituType |


## Datos

Los datos son públicos y deben descargarse desde Bruker Spatial Biology (antes NanoString):

**[CosMx SMI Human NSCLC FFPE Dataset](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/)**

> **Nota:** la versión degradada de los datos (`thin_75`) se generó mediante un procedimiento interno del laboratorio que no forma parte de este repositorio. El pipeline parte de los objetos ya degradados.

## Preprocesamiento (fuera del pipeline)

Por requerimientos de memoria, los siguientes pasos se ejecutaron en un dispositivo
del laboratorio y no forman parte de este repositorio:

- Unión de las muestras de los 5 pacientes (Lung5_Rep3, Lung6, Lung9_Rep1,
  Lung12, Lung13) en un único objeto Seurat: 960 genes, 416.087 células
  en el escenario original y 398.813 en el degradado (`thin_75`).
- Control de calidad, normalización, PCA (50 componentes) y UMAP, realizados con
  el código del laboratorio.
- Ejecución de Leiden e InSituType mediante los scripts de `codigo/Metodos/`,
  preparados para recibir parámetros por línea de comandos.

Leiden se ejecutó directamente sobre el objeto Seurat integrado. Para InSituType
se usó el mismo objeto con las capas de conteos unidas (`JoinLayers`), formato que
requiere el método. **Ambos métodos se aplicaron sobre las mismas células y genes,
y a partir de sus resultados todo el análisis (anotación, marcadores, comparación
y figuras) es idéntico para los dos.**

Los objetos resultantes (`*_joinmerged.rds`) y los datos previos al control de
calidad (`datos_histograma.rds`, usados solo para las figuras de QC) se recibieron
ya procesados y son el punto de partida del pipeline local.

### Perfil de referencia

InSituType utiliza como referencia el perfil **IO_1k** (panel de 960 genes), disponible en el repositorio [Nanostring-Biostats/CosMx-Cell-Profiles](https://github.com/Nanostring-Biostats/CosMx-Cell-Profiles/tree/main/Human/IO_1k). El archivo `IO_1k.profiles.csv` debe descargarse y ubicarse en `codigo/referencia_anotacion/`.

## Estructura del repositorio

```
codigo/
├── Metodos/                     # Pipeline: un script por paso (se ejecutan con Rscript)
│   ├── 01_QC.R                  # Control de calidad, sondas negativas y normalización
│   ├── 02_Seurat.R              # Clustering Leiden
│   └── 03_IST.R                 # InSituType
├── 01_Resultados_clustering.R   # Une los resultados de ambos métodos en un único archivo
├── 02_Clustering_joinmerged.R   # Añade los resultados de clustering a los objetos Seurat
├── 03_Marcadores.R              # Genes marcadores por cluster (FindAllMarkers)
├── 04_Resultados_IST.qmd        # Reporte de resultados: InSituType
├── 04_Resultados_Leiden.qmd     # Reporte de resultados: Leiden
├── 05_Labels.R                  # Unificación de etiquetas de tipo celular
├── 06_Confusion.R               # Matrices de confusión y ARI entre métodos y escenarios
├── funciones.R                  # Funciones auxiliares
├── qc_plots.R                   # Figuras de control de calidad
└── resultados_pacientes.R       # Resultados por paciente
```

## Uso

Los scripts de `Metodos/` se ejecutan desde la terminal:

```bash
Rscript codigo/Metodos/01_QC.R <input_raw.rds> <out_subset.rds> <out_neg_vector.rds>
Rscript codigo/Metodos/02_Seurat.R <input_preprocesado.rds> <dims> <resolution> <directorio_salida>
Rscript codigo/Metodos/03_IST.R <input.rds> <neg_vector.rds> <reference.csv> <n_clusts> <directorio_salida>
```

Después se ejecutan en orden los scripts `01_` a `06_` de `codigo/` y se renderizan los reportes `.qmd`.

## Resultados

Los reportes renderizados con los resultados están disponibles en: *(pendiente — GitHub Pages)*

## Requisitos

- R 4.4.2
- Quarto 1.8.27
- Paquetes de R: `Seurat`, `InSituType`, `clustree`, `mclust`, `tidyverse`, `patchwork`, `pheatmap`, `DT`, `kableExtra`, `knitr`, `here`

## Licencia

Código bajo licencia [MIT](LICENSE).
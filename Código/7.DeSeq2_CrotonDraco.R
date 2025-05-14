#############################################################################
# Differential gene expression analysis 
# Workflow: [Analyzing RNA-seq data with DESeq2](https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts) |Mirna Vazquez-Rosas-Landa|Anahí Canedo-Téxon
############################################################################
# INSTALACION DE PAQUETERIAS

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

# EN TERMINAL(desde BASH):
# sudo apt-get installl ibxml2-dev
# sudo apt-get installl libcurl4-openssl-dev
# 3. sudo apt-get -y install libcurl4-gnutls-dev libxml2-dev libssl-dev 

# Desde R:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")

##############################################################
install.packages("pheatmap")
install.packages("ggplot2")
install.packages(" RColorBrewer")
install.packages(" tximport")
##############################################################
# CARGAR PAQUETERIAS
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("tximport")
search()
#############################################################################
getwd()
setwd("../DifExpGen/")

dir <- getwd()

#Cargar y explorar uno de los archivos de genes.results.

hojas <- read.csv("Hoja.genes.results", sep = "\t")
head(hojas)
dim(hojas)
##############################################################################
# 1. Cargar todos los archivos a partir de un directorio determinado.
# Cargar el archivo de settings, donde se establecen los nombres de las condiciones del experimento, y se les asigna un nombre
samples <- read.csv("2.SamplesDescribed.csv", header = TRUE, sep = "\t") 

# Leer el archivo de settings y cargar de manera automatizada todos los archivos al mismo tiempo
files <- file.path(dir, samples$Name)
names(files)
# Preguntamos si se cargaron todos los achivos establecidos en el settings
all(file.exists(files)) # ¿Existen los 6 archivos del objeto samples?
# Asignar nombres a los archivos 
print(files) # Hasta ahora el nombre de los archivos, es la ruta completa
names(files) <-gsub("sub.","", samples$Name, perl = T) # Vamos a modificar los nombres, utilizando la función gsub ()
print(names(files)) # ¿Qué paso?

##############################################################################
# 2. Convertir los reusltados de cada biblioteca, a una matriz
?tximport()
txi <- tximport(files, type = "rsem")

# Explora el contenido del objeto
names(txi)
head(txi$counts) #rawcounts
head(txi$abundance) #TPM
head(txi$length) #effective_length del gen 

##############################################################################
# Convertir la matriz de conteos a un objeto tipo DESeqDataSet para que se pueda reconocer por el siguiente programa/funcion
?DESeqDataSetFromTximport()
dds <- DESeqDataSetFromTximport(txi, samples, ~1) 
class(dds)
str(dds)
dds
colData(dds)

# En caso de tener datos faltantes convertir 0 a 1
# which(txi$length == 0)
# zero <- which(txi$length == 0)
# txi$length[zero] <- 1
# which(txi$length == 0)
#dds <- DESeqDataSetFromTximport(txi, samples, ~1)


##############################################################################
# 3.Core analysis para seleccionar genes expresados diferencialmente
?DESeq()
dds <- DESeq(dds) #Esta funcion tiene 3 funciones integradas, estimationSizeFactors, estimateDispersions, nbinom....
res <- results(dds)
##############################################################################
?vst()
vsd <- vst(dds, blind=FALSE) #Extraer información del objeto dds, para observar la tendecia de dispersión 

##############################################################################
# 4.PCA
plotPCA(vsd, intgroup=c("Tratamiento"))

##############################################################################
# 5. HEATMAP

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Tratamiento, vsd$Dia, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)





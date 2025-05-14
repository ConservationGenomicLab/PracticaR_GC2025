###############################################################################
# CARGAR BASES DE DATOS
###############################################################################

# Conoce tu directorio actual de trabajo, copia la estructura de la dirección y modificala de acuerdo al sitio dónde se encuentren tus archivos.
getwd() 
#Establece tu directorio de trabajo:
setwd("/home/anahi/Escritorio/Primeros-Pasos-en-R/DataBases/") #establecer la ruta al directorio en que se encuentran los archivos que deseas leer

###############################################################################  
# Archivos .txt
# Los archivos con extensión .txt pueden ser leídos empleando la función read.table().  
# Aunque no son tan comunes para leer en R.
###############################################################################  

arboles <- read.table("Arboles.txt", sep = " ", header = T)
dim(arboles)
head(arboles)             

############################################################################### 
# Archivos .csv 
# Los archivos con extensión .csv pueden ser leídos empleando la función read.csv().  
###############################################################################

?read.csv()

pob <- read.csv("Poblaciones.csv")
head(pob)          
pob <- read.csv("Poblaciones.csv", sep = "\t", header = T)
head(pob)     
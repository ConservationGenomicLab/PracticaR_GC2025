###############################################################################
## FFUNCIONES PARA EXPLORAR DATAFRAME
###############################################################################
# Vamos a usar como ejemplo una base de datos pre-cargada en R base.


data() #Explorar las opciones de las bases de datos precargadas en R base
data(PlantGrowth) #Seleccionar la base de datos "PlantGrowth": Results from an Experiment on Plant Growth
---

str(PlantGrowth) #obtener el tipo de objeto/vector de PlantGrowth
dim(PlantGrowth) #obtener las dimensiones del objetos PlantGrowth
names(PlantGrowth) #debido a qué es un data frame (compuesto de columnas y renglones), queremos conocer el nombre de las etiquetas de cada columna
summary(PlantGrowth) #nos muestra datos estadísticos de cada uno de los vectores en el data
frame

---
# Acceder a los elementos, o columnas del data frame

PlantGrowth$weight
PlantGrowth$group
?attach() # Permite acceder a los objetos en la base de datos simplemente dando sus nombres.
attach(PlantGrowth) 


weight
group
---
# Seleccionar elementos de un data frame
# nombre[renglones, columnas]  
  PlantGrowth[1,] #Seleccionar el renglón 1
  PlantGrowth[,1] #Seleccionar la columna 1
  PlantGrowth[1:5,] #Seleccionar un rango de renglones
  PlantGrowth[-1,] #Excluir el renglón 1
  PlantGrowth[-(1:10),] #Excluir los renglones del 1 al 10
---

ncol(PlantGrowth) #obtener el numero de columnas del data frame
nrow(PlantGrowth) #obtener el numero de renglones del data frame


#explorar las categorias o los diferentes niveles de los grupos de tratamientos en el experimento
levels(group)

#Contar el número de categorías que resultan de la función anterior:
#Funciones aninadas
length(levels(group)) #¿cuántos tratamientos se realizaron en el experimento?

#############################################################################
# CREANDO UN SUBSET
###############################################################################
  # ¿Cuántas plantas se utilizaron por cada tratamieto? 
  # Hacer un subset de cada uno de los grupos  de tratamientos en el experimento
?subset()
ctrl <- subset(PlantGrowth, group =="ctrl")
c1 <- nrow(ctrl)
trt1  <- subset(PlantGrowth, group =="trt1")
c2 <- nrow(trt1)
trt2  <- subset(PlantGrowth, group =="trt2")
c3 <-nrow(trt2)

---
#Guardar el conteo de plantas por cada tratamiento en una tabla:
results <- as.table(cbind(c1,c2,c3))

#Renombrarlas columnas y renglones:
colnames(results) <- c("control", "nitrogeno", "fosforo")
names(results)
row.names(results) <- "No.Plantas"
print(results)


---
#############################################################################
# GENERACIÓN DE GRÁFICOS
###############################################################################

barplot(PlantGrowth$weight) # Gráfico de barras, del peso de las plantas
hist(PlantGrowth$weight) #Histogramas
boxplot(PlantGrowth$weight~PlantGrowth$group) #Boxplot por tratamiento
---
#Gráfico de barras, del peso de las plantas por cada tratamiento
  
  
peso.mean <- as.table(cbind(mean(ctrl$weight),mean(trt1$weight),mean(trt2$weight)))

colnames(peso.mean) <- c("control", "nitrogeno", "fosforo")
peso.mean

# Enchular el gŕafico 
barplot(peso.mean, beside= T, main = "Promedio de Crecimiento", ylab = "Peso(gr.)", xlab = "Tratamientos", col = c(4,26,39))

# Exportar figuras
png("PlantGrowth.png")
barplot(peso.mean, beside= T, main = "Promedio de Crecimiento", ylab = "Peso(gr.)", xlab = "Tratamientos", col = c(4,26,39))
dev.off()


              
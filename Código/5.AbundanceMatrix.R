#######################################################
#
# Explorando matriz de abundancias generada por RSEM 
# para realizar preparar el análisis de expresión génica con
# datos provenientes de RNA-Seq
# 
#######################################################
library(dplyr)
#######################################################
data <-  read.table("../DifExpGen/AbundanceMatrix.isoform.counts.matrix" , sep = "\t", header = T)
names(data)[1] <-"GeneID"
names(data)
class(data)
head(data)
ncol(data)
dim(data)

#######################################################
#  Remover  unigenes con 0 cuentas en todas las columnas (condiciones)
#######################################################
data.clean <-  data %>%
  filter(!if_all(-GeneID, ~ . == 0))   

#######################################################
# 2.1 Obtener promedio por réplicas para Flores
#######################################################
averageTx1 <-  data.clean %>%
  head() %>%
  rowwise() %>% # Trabajar fila por fila
  mutate(averageTx1 = mean(c_across(2:5))) 
dim(averageTx1)
#######################################################
# 2.2 Obtener promedio por réplicas de cada tejido
#######################################################
data.clean <-   data.clean %>%
  rowwise() %>% # Trabajar fila por fila
  mutate(averageTx1 = mean(c_across(2:4))) %>%
  mutate(averageTx2 = mean(c_across(5:7))) %>%
  mutate(averageTx3 = mean(c_across(8:10))) %>%
  mutate(averageControl = mean(c_across(11:13))) %>%
  ungroup() # Liberar el modo rowwise
names(data.clean) 

dim(data.clean)
data.clean$
#######################################################
# Crear un dataframe unicamente con las columnas de promedios entre réplicas
#######################################################
data.mean <- data.clean %>% 
select(GeneID, averageTx1, averageTx2, averageTx3, averageControl) 
names(data.mean)
dim(data.mean)
summary(data.mean)

#######################################################
# Filtrar unigenes con al menos 10 cuentas mapeadas en cada experimento
#######################################################
c10.Tx1 <-  data.mean %>% filter(averageTx1 > 10) 
c10.Tx2 <-  data.mean %>% filter(averageTx2 > 10)    
c10.Tx3 <-  data.mean %>% filter(averageTx3 > 10) 
c10.C <-  data.mean %>% filter(averageControl > 10)

c10.number <- data.frame(
  "Tx1" = nrow(c10.Tx1), 
  "Tx2" = nrow(c10.Tx2),
  "Tx3" = nrow(c10.Tx3),
  "Control" = nrow(c10.C))


#######################################################
# Obtener genes con conteos mayoritarios para cada condición
#######################################################
cm.Tx1 <-  data.mean %>%
  rowwise() %>%
  filter(all(averageTx1 > c(averageTx2 & averageTx3 & averageControl ))) %>% #verifica si averageTx1 es mayor que todas las otras columnas en la fila
  ungroup() 

cm.Tx2 <-  data.mean %>%
  rowwise() %>%
  filter(all(averageTx2 > c(averageTx1 & averageTx3& averageControl ))) %>%
  ungroup()

cm.Tx3 <-  data.mean %>%
  rowwise() %>%
  filter(all(averageTx3 > c(averageTx2 & averageTx1& averageControl ))) %>%
  ungroup()
 

cm.C <-  data.mean %>%
  rowwise() %>%  # Evalúa fila por fila
  filter( all(averageControl > c(averageTx2, averageTx3, averageTx1)))%>%  
  ungroup()


# Crear un vector para almacenar los resultados
cm.Genes <- data.frame(
  Control = nrow(cm.C),
  Tx1 = nrow(cm.Tx1),
  Tx2 = nrow(cm.Tx2),
  Tx3 = nrow(cm.Tx3) )

###############################
# rename()
##############################
Renombrar columna en el dataframe
data.mean <- rename(data.mean, pob1 = averageTx1 ) # renombrer una columna
names(data.mean)

data.mean <- data.mean %>% # renombrer varias columnas
  rename("pob1"= "averageTx1", "pob2" = "averageTx2", "pob3" = "averageTx3")

############################################################
# Obtener foldchange Control vs. Tratamiento1
############################################################
# FC = cuentas gen1 (Pob1) / cuentas gen1 (Control)
############################################################
data.mean %>% 
mutate(averageControl = ifelse(is.infinite(averageControl), min(data.mean$averageControl), averageControl))%>%
  head()

data.mean <- data.mean %>% 
  rowwise() %>% # Trabajar fila por fila
  mutate(FC = pob1 /averageControl) %>% 
  select(GeneID,pob1, averageControl, FC) 
dim(data.mean)
boxplot(log10(data.mean$averageControl), log10(data.mean$pob1))

FCmax <-  data.mean %>%
   filter(FC > 2) %>%
  mutate(FC = ifelse(is.infinite(FC), NA, FC)) %>%
  filter(!is.na(FC))

FCmin <-  data.mean %>%
  filter(FC < 2)  %>%
  mutate(FC = ifelse(is.infinite(FC), NA, FC)) %>%
  filter(!is.na(FC))
  

################################################################################
# Matriz de Abundancias Calculadas con datos de RNA-Seq del transcriptoma de una especie arbórea, para hacer un heatmap y observar el perfil génico de cada condición 
################################################################################
setwd("DataBases/")
# Cargar y explorar los datos
data2 <- read.csv("../DifExpGen/AbundanceMatrix.TMM.fpkm.matrix", sep = "\t", header = T)
colnames(data2) <- c("Unigenes" , "Hoja", "Fruto" ,"Tallo" ,"Raíz"  ,"Inf_1" ,"Inf_2")
rownames(data2) <- data$Unigenes
dim(data2)
################################################################################
SUB <- data2[1:171474,]
SUB2 <- SUB[1:100000,2:7]

SUB3 <- scale(data2[1:10000 , 2:7])

library(RColorBrewer)
col <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
heatmap(SUB3, main = "DEG en Transcriptoma C.D.", col = col,  xlab= "Organos", ylab= "Unigenes")

### heatmap(as.matrix(data[1:10000,2:7]))  ### NO CORRER EN PC PERSONAL
#------------
# Guardar
png("Heatmap_all.png")
heatmap(SUB3, main = "DEG en Transcriptoma C.D.", col = col, xlab= "Organos", ylab= "Unigenes")
dev.off()
#------------


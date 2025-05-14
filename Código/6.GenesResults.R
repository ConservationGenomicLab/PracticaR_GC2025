#######################################################
#
# Explorando dataframe para valores de expresión génica provenientes de RNA-Seq
# 
#######################################################
library(dplyr)
#######################################################
datatx1 <-  read.csv("../DataBases/Tx1.rep1.genes.results.csv" , sep = ",", header = T, )
names(datatx1)

class(datatx1) # clase de objeto
dim(datatx1) #obtener las dimensiones del objeto 
ncol(datatx1) # numero de columnas
nrow(datatx1) # numero de genes en el datatx1 frame

---
  # Conocer cual es el valor mínimo y máximo de una variable:
  
min(datatx1$expected_count)
max(datatx1$expected_count)
---
  # Identificar el gen de mayor longitud:
  
max.Length <-row.names(datatx1[which.max(datatx1$length), ] )
x.Length <-row.names(datatx1)[datatx1$length=="615"]
y.Length <- (datatx1)[datatx1$length > 10000, ] 
---
  # Función subsets:
  
  y.Length <- (datatx1)[datatx1$length > 10000, ]
z.Length <- subset(datatx1, datatx1$length > 10000 ) 
max.Counts <- subset(datatx1, datatx1$expected_count== max(datatx1$expected_count))

---
# Obtener un subset de genes con un número de conteos determinado:
  
threshold.1 <- subset(datatx1, datatx1$expected_count > 1000)
length(threshold.1$transcript_id.s.)

---
  # Es tu turno:
  # Obten genes que con una longitud mayor a 1000 pb y que tegan más de 50 conteos asignados por RSEM.
  
  ---
############################################
#  dplyr::filter():
#############################################
library(dplyr)

datatx1 %>% 
  filter(expected_count > 1000) %>%
  count()


---
############################################
#  dplyr::summarize()
#############################################
#  ¿Cómo se distribuyen los valores de TPM?
datatx1 %>%
  summarise(TOTAL = n() , MEAN= mean(TPM) , SD= sd(TPM), VAR= var(TPM)) %>% 
  datatx1wizard::datatx1_rotate(colnames = T, rownames = "TOTAL")


############################################
#  dplyr::mutate()
#############################################
# Transformación logaritmica de los conteos observados para cada gen:
datatx1 <- datatx1 %>% 
  mutate(LOG = log10(expected_count))
datatx1$LOG

# Calcular manualmente el valor de FPKM
head(datatx1)

counts.Sample <- sum(datatx1$expected_count)
FPKM <-  datatx1 %>% 
  mutate(FPKM = expected_count*10^9 / (length*counts.Sample)) 

names(datatx1)

---
  #  Tu turno:
  # Crea un nuevo objeto donde puedas obtener genes con valores de FPKM > 15 considerados como "sobrerepresetnados"; y valores de FPKM < 15 considerados como "reprimidos".
  
  UpExpression <-   
  
  DownExpression <-
  ---

############################################
#  dplyr::select()
#############################################

datatx1 %>% 
  select(transcript_id.s.,length, expected_count)%>% 
  head()
---
  # ¿Podemos combinar funciones dplyr?
  #  filter() %>% select() %>% mutate() 
  
  ---
  # Obtener unigenes con valores de FPKM mayores a 5, y con una longitud mayor a 1000 pb 
  # Abundancia Relativa= Valor en cada celda/ Suma total de la columna
  
  AbRel <-  OverExpression %>% 
  select(FPKM, length) %>% 
  filter(FPKM > 5 & length >1000) %>% 
  mutate(AbunRel= FPKM/sum(datatx1$FPKM))

##############################################
# Vamos a comparar control vs tratamiento:
##############################################
library(tidyverse)

tx1 <-  read.table("../DataBases/Tx1.rep1.genes.results" , sep = "\t", header = T)
control <- read.table("../DataBases/Control.rep1.genes.results.csv" , sep = ",", header = T) 

# Calcula los valores de FPKM para el dataframe del control.




# Unir data frames, usando el nombre de los genes
merged <- inner_join(control, tx1, by = "transcript_id.s.", suffix = c(".ctrl", ".tx1"))

# Calcular el logaritmo de la tasa de cambio tx1/control: 
merged <- merged %>%
  mutate(log2FC_FPKM = log2(FPKM.tx1/FPKM.ctrl)) 

head(merged)

# Definir punto de corte para up/down regulated genes: 
upGenes <- merged %>%
  filter(is.finite(log2FC_FPKM),abs(log2FC_FPKM) > 5)                 # elimina Inf y Na  / # fold change log2 mayor a 1 o menor a -1

downGenes <- merged %>%
  filter(is.finite(log2FC_FPKM),abs(log2FC_FPKM) < 5) 

dim(upGenes)
dim(downGenes)
 

library(ggplot2)

# Crear columna de estado de expresión
diffexp <- merged %>%
  mutate(expression = case_when(
    log2FC_FPKM > 10  ~ "UpRegulated",
    log2FC_FPKM < -10 ~ "DownRegulated"
  ))

# Graficar
ggplot(diffexp, aes(x = log2FC_FPKM, y = 0, color = expression)) +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("UpRegulated" = "red", "DownRegulated" = "blue")) +
  labs(
    title = "Diferential Expression Genes",
    x = "log2 Fold Change (FPKM)",
    y = "",
    color = "DEG"
  ) +
  theme_minimal()

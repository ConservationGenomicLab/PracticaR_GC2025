###############################################################################
# ESTRUCTURAS DE DATOS 
###############################################################################

# Tipos de objetos en R 

###############################################################################
# Objetos simples
###############################################################################
  # - numeric: números reales o decimales
  # - integer: números enteros
  # - complex: números complejos
  # - character: carácteres
  # - factor: registrar datos cualitativos.
  # - logical: resultados lógicos (TRUE o FALSE)

---
# Para conocer la clase de cada objeto podemos emplear la función class() , por ejemplo
# class(objetos)

---

# NUMÉRICOS
###############################################################################
precios <- c(20,50,100)
class(precios)

---
# CARACTERES
###############################################################################
  
nombres <- c("Jose", "Pepe", "Mariana", "Daniel", "Oscar", "Paola")
class(nombres)

--- 
# FACTORES
###############################################################################
sexo <- factor(c("Femenino", "Masculino"))
class(sexo)
str(sexo)
levels(sexo)
---
# LÓGICOS
###############################################################################
  
x <- 1
y <- 2   
respuesta <- (x > y)  # ¿x es mayor que y?; guardar el resultado en el objeto r
respuesta             # imprimir el resultado

class(respuesta)       # clase del objeto respuesta
---
  
###############################################################################  
## Objetos compuestos
###############################################################################
  # - Listas
  # - Data frames
  # - Matrices
  # - Arreglos
  
---
# LISTAS
# Función list()
###############################################################################
  

Colecta <- list(especie = "Piper nigrum", muestras= 200,
                  sitio = c("Bosque Mesófilo", "Bosque Tropical","Selva") )
Colecta
class(Colecta)

---
# MATRICES
# Función matrix()
###############################################################################

matrix(1:15)
matrix(1:15, nrow = 3, byrow = TRUE)
matrix(1:15, nrow = 3, ncol = 3)
mi.MATRIX <- matrix(1:15, nrow = 4, ncol = 4)

dim(matrix(1:15, nrow = 4, ncol = 4))
################################################################################
# TO/DO: 
# Crea una matriz con  50 números.
# Identifica las diferencias obtenidas en cada uno de los siguientes casos:
# a) Con 5 columnas y 10 renglones.
# b) Con 10 columnas y 10 renglones
# c) Con sólo 2 columnas 


################################################################################

---
# ARREGLOS
# SINTAXIS: array(datos, dimensiones)
###############################################################################
miarray <- array(data=LETTERS[1:24], dim=c(3, 4, 2))
#Si queremos extraer la letra almacenada en la fila 1 y columna 3 de la segunda capa de miarray usamos el siguiente código:
miarray[1, 3, 2]
#Si queremos extraer la segunda capa completa usamos el siguiente código.
miarray[,, 2]

################################################################################
# TO/DO: Extrae la letra "k"
################################################################################


---
# DATA FRAME
# Función data.frame()
###############################################################################
df <- data.frame(  "id" = 1:4, 
  "profesion" = c("musico", "pintor", "enfermedo", "actor"), 
  "edad" = c(20, 30, 18, 32),
  "nombre" = as.character(c("Ana", "Pepe", "Maria", "Arturo")) )

df
class(df)


  
  
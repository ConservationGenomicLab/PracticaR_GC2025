###############################################################################
# CREACIÓN Y MANIPULACIÓN DE OBJETOS
###############################################################################

# Asignar objetos: 
a<-20
20 -> a
a=20

# Almacenar variables en un objeto:
  
f1 <- "manzana"
f2 <- "pera"

---
## R cómo calculadora 

8+9 #sumas
1000-84 #restas
200/9 #división
5*10 #multiplicación
10^2 #exponenciales

## Operaciones con objetos

a=10
b=50
c=2

(a+b)^c

---

---
# Enlistar los objetos en el entorno
  
ls()
objects()

# Borrar objetos
rm(a)
ls()

---
  
# Consultar el manual de ayuda en R de las funciones
  
?ls()
?sum()
?sqrt()
?mean()
?median()


---
################################################################################
# TO/DO: 
# 1.Crea 5 objetos numéricos (enteros) de 3 dígitos
# 2.Obten el promedio, meadina y la cuadrática de los objetos que creaste
# 3.Guarda el resultado de cada operación en un objeto

  
################################################################################
# CONCATENAR ELEMENTOS

#Crear un objeto que almacene el nombre de los estudiantes de un grupo

estudiantes <- c("Jose", "Pepe", "Mariana", "Daniel", "Oscar", "Paola")
print(estudiantes)

#Crear un objeto que almacene la edad de cada estudiante
edades <- c(10, 10, 9, 8, 9, 11)
print(edades)

#Crear un objeto con la estatura de los estudiantes
alturas <- c(119, 118, 110, 116, 111, 120)
print(alturas)
alturas/2
---
# Crear una tabla con la información
  
  
cbind(estudiantes,edades, alturas) #concatenamos los objetos por columnas
grupoA <- cbind(estudiantes,edades, alturas) # asignamos la salida de la función cbind, a un nuevo objeto
print(grupoA)
dim(grupoA) #dim() nos permite conocer las dimensiones de un objeto

---
# Extraer información de una tabla

print(grupoA)
grupoA[6,] #obtener la fila de la posición 6, correspondiente a la información de Paola
grupoA[,3] #obtener la columna de alturas de de todos los estudiantes
grupoA[6] #obtener el nombre del estudiante Paola

grupoA[6,3] #obtener la estatura de Paola 

################################################################################
# TO/DO: Obten la edad de Mariana
 
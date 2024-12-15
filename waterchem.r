rm(list = ls())

library(cglasso)
library(corrplot)

#cargamos la base de datos y la especificación de los nombres de las variables
waterchem <- read.csv(gzfile("nla22_waterchem_wide.csv.gz"))
waterchem_specs <-read.delim("nla22_waterchem_wide.txt")

#El objetivo de esta sección de código es dejar los nombres de las variables en forma prolija
#y generar un archivo: variables_description.csv que contiene el nombre de cada variable numérica
#y su correspondiente descripción

#nos quedamos con las columnas que tienen datos numéricos
results <- waterchem[, grep("RESULT", names(waterchem))]
data <- results[, -grep("UNITS", names(results))]
#le sacamos la palabra RESULT a esas columnas
var_names_orig <- colnames(data)
variables_names <- data.frame(PARAMETER = var_names_orig)
variables_description <- merge(variables_names,waterchem_specs)
variables_description <- variables_description[, colnames(variables_description) %in% c("PARAMETER", "LABEL")]
variables_description$PARAMETER <- gsub("_RESULT", "", variables_description$PARAMETER)
write.csv(variables_description, "variables_description.csv", row.names = FALSE)

colnames(data) <- gsub("_RESULT", "", colnames(data))

#sacamos una columna con solo seis valores no nulos: la de nitrato-nitrito
sum(!is.na(data$NITRATE_NITRITE_N))
data_complete <- data[, !colnames(data) %in% c("NITRATE_NITRITE_N")]

#realizamos un heatmap de la matriz de covarianzas sin estandarizar los datos
cov_matrix <- cov(data_complete, use='pairwise')
heatmap(cov_matrix)


#invertimos la matriz de covarianzas para obtener la matriz de precisión y la graficamos
library(matlib)

precision_matrix <-inv(cov_matrix)
rownames(precision_matrix) <- rownames(cov_matrix)
colnames(precision_matrix) <- colnames(cov_matrix)
heatmap(precision_matrix)

#estandarizamos los datos y lo gruardamos en el dataframe data_complete_std
#y realizamos gráficos de la matriz de covarianzas que es una matriz de correlación
#por estar los datos estandarizados
cor_matrix <- cor(data_complete, use = "pairwise.complete.obs")
data_complete_std <- scale(data_complete)
cov_matrix_std <- cov(data_complete_std, use='pairwise')
heatmap(cov_matrix_std)
corrplot(cov_matrix_std, method = "circle", tl.cex = 0.8, tl.col = "black")

#calculamos la matriz de precisión a partir de invertir la matriz de covarianzas y la graficamos
precision_matrix_std <-inv(cov_matrix_std)
rownames(precision_matrix_std) <- rownames(cov_matrix_std)
colnames(precision_matrix_std) <- colnames(cov_matrix_std)
heatmap(precision_matrix_std)
corrplot(precision_matrix_std, method = "circle", tl.cex = 0.8, tl.col = "black", is.corr=FALSE, 
         col.lim=c(min(precision_matrix_std),max(precision_matrix_std)))


#calculamos la matriz de correlaciones parciales a partir de la matriz de covarianzas empleando
#la función cov2pcor del paquete gRbase, pero la matriz resultante tiene la mayoría de sus coeficinentes nulos
library(gRbase)

pcor_matrix <- cov2pcor(cov_matrix_std)
corrplot(pcor_matrix, method = "circle", type = "upper", tl.cex = 0.8, tl.col = "black")

corrplot(cov_matrix_std, method = "circle", type = "upper", tl.cex = 0.8, tl.col = "black")
corrplot(cor_matrix, method = "circle", type = "upper", tl.cex = 0.8, tl.col = "black")

#este es un fragmento de código que intenta calcular las correlaciones parciales a través de un método
#de regresión, pero que no resulta computacionalmente adecuado

#library(generalCorr)
#intentamos emplear esta librería pero no terminaba de realizar los cómputos en tiempos razonables
#parcor_mtx <- parcorMtx(data_complete_std, verbo  = TRUE)
#parcor_mtx2 <- parcorBMany(data_complete_std, idep = 1, blksiz = 10, verbo = TRUE)
# <- parcorSilent(cov_matrix_std, idep = 1, verbo = TRUE)


#Esta es la parte central del código donde se genera la estructura de datos (cggm) que requiere
#el paquete cglasso y su función principal homónima para funcionar,
#y donde asimismo se realiza la aplicación del algoritmo propiamente dicho
#En los argumentos de la llamada a la función cglasso
#en primer lugar se especifica la fórmula que, en este caso, indica que el modelo a ser ajustado
#es de cada variable contra todas las demás
#maxit.bcd indica el máximo de iteriaciones a realizar, que nosotros seteamos a un valor alto por las dudas
#para que la cantidad de iteraciones no sea un impedimento para la convergencia del algoritmo
#rho es una secuencia descendente que indica los distintos pesos que se le da a la regularización L1 o Lasso
#finalmente, convertimos el resultado a un grafo
#como no especificamos en este llamado a to_graph, el parámetro rho.id, tomará el grafo correspondiente
#al valor de rho más bajo de nuestra secuencia

Z <- datacggm(data_complete_std)
system.time(out <- cglasso(. ~ ., data = Z, maxit.bcd = 2e9, rho = c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1)))
out.graph <- to_graph(out) #probar, por ejemplo rho.id = 5

#ploteamos el grafo
plot(out.graph, type = "Gyy")

#introducimos una curvatura en la graficación de los enlaces para que sean más
#fáciles de leer
out.igraph <- getGraph(out.graph, type = c("Gyy"))
plot(out.igraph, vertex.label.family="Times", edge.curved=0.1)

#mostramos el grafo pero con pesos en los enlaces, que afecten el grosor del enlace
out.graph_weight <- to_graph(out, weighted = TRUE, rho.id = 5)
E(out.graph_weight$Gyy)$weight = -E(out.graph_weight$Gyy)$weight
out.igraph_weight <- getGraph(out.graph_weight, type = c("Gyy"))
plot(out.igraph_weight, vertex.label.family="Times", edge.curved=0.1, edge.lty=c("solid"), 
     edge.width=(E(out.igraph_weight)$weight/min(E(out.igraph_weight)$weight)))


#en esta sección de código calculamos las correlaciones parciales empleando el
#paquete psych que permite trabajar con los datos originales, en lugar de
#realizar un cálculo a partir de la matriz de covarianzas
#Empleamos este paquete porque permite trabajar con valores nulos
#la inclusión de la variable nitrato-nitrito rompía esta parte del código
#probablemente por el exceso de valores nulos, así que, en la versión actual,
#esta columna se quita al inicio del código, como ya aclaramos.
library(psych)

#los coeficientes de correlación parcial exceden el rango [-1,1] por la presencia de datos faltantes
partial_cor_matrix <- partial.r(data=data_complete, method = "pearson")
corrplot(partial_cor_matrix, is.corr=FALSE, col.lim=c(min(partial_cor_matrix),max(partial_cor_matrix)),
         method = "circle", type = "upper", tl.cex = 0.8, tl.col = "black")



#En esta sección graficamos las correlaciones parciales a partir de las
#matrices de covarianzas empírica y estimadas
#que nos entrega la función cglasso: hay una matriz estimada para cada valor de 
#rho que hayamos especificado en la llamado original a la función cglasso
#aquí como se grafica una detrás de otra queda el foco puesto en la última,
#es decir para el valor de rho más chico
library(corpcor)

#working empirical covariance matrices
for (i in 1: dim(out$S[,,,])[3]){
pcor_matrix <- cor2pcor(cov2cor(out$S[,,,i]))
variable_names <- colnames(data_complete)
colnames(pcor_matrix) <- variable_names
rownames(pcor_matrix) <- variable_names
corrplot(pcor_matrix, method = "circle", type = "upper", tl.cex = 0.8, tl.col = "black")
}

#estimated covariance matrices
for (i in 1: dim(out$S[,,,])[3]){
  pcor_matrix <- cor2pcor(cov2cor(out$Sgm[,,,i]))
  variable_names <- colnames(data_complete)
  colnames(pcor_matrix) <- variable_names
  rownames(pcor_matrix) <- variable_names
  corrplot(pcor_matrix, method = "circle", type = "upper", tl.cex = 0.8, tl.col = "black")
}



#esta sección comentada de código intentaba emplear el paquete huge para realizar los análisis
#library(huge)
#out1 <- huge(cor_matrix, method = "glasso", lambda = seq(log(0.1), log(0.001), length = 30) )
#plot(out1)

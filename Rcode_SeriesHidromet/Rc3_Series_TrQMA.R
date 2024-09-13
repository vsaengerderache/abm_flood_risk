# ..................................................................................................
# Proyecto de sociohidrologia de inundaciones
# Analisis de Caudales Máximos Anuales (QMA)
#
# Autor: Jorge Hurtado
# Fecha: 22 de julio de 2024
# 
# Objetivos: 1) Generar una serie de QMA y TR 
#
# ..................................................................................................

# Preparing the environment

cat("\014")       #clean console   
rm(list = ls())   #clean environment
graphics.off()    #clean plots

# Libraries

{
  #library(data.table)
  library(hydroTSM)
  #library(dplyr)
  #library(ggplot2)
}

# ..................................................................................................

# 1. Datos CR2-CAMELS (descargados) ----

# Load data

setwd("C:/JorgeR/SocioHidrologiaUdeC/RespositorioInformacion/2_CR2_Camels/camels_cl_8530001")
flow <- read.csv("q_m3s_day.csv", header = TRUE)

# Preparing data

names(flow) <- c("date","flow")
flow$date2 <- as.Date(flow$date, format="%Y-%m-%d"); head(flow)
flow <- na.omit(flow); head(flow)

# Converting to Series 

flow.ts = zoo (flow$flow, order.by=flow$date2) # con la funcion zoo se convierte la serie de tiempo en objeto de serie de tiempo
head(flow.ts)
plot(flow.ts, xlab="", ylab="Daily Flow, m3/s")

flow.anual <- daily2annual(flow.ts, FUN=max, na.rm=TRUE) # caudales maximos anuales
plot (flow.anual, xlab="years",ylab="m3/s", main="CR2 Annual Maximum Flow")

# Serie QMA y barplot

n = 30 # numero de años de la serie

barplot(flow.anual,ylim=c(0,1000))

t <- n/length(flow.anual);t # numero de veces para repetir la serie en 50 años

serie.qma <- rep(flow.anual,t,length.out=n)
serie.qma <- data.frame(serie.qma)
years <- 1:n
#row.names(serie.qma) <- years

barplot(serie.qma$serie.qma, names.arg = years, ylim=c(0,1000), 
        xlab="Time Step (year)", ylab="QMA (m3/s)",main="Caudales Máximos Anuales")
abline(h=690, col="red") # tr10 caudal de diseño para los escenarios con proyecto

# Serie Tr y barplot

# periodos de retorno de los QMA de carampangue

qt1 <- 280; qt2 <- 487; qt5 <- 610; qt10 <- 690; qt25 <- 794; qt50 <- 869; qt100 <- 947; qt200 <- 1026

# Crear una función para asignar la etiqueta basada en el valor más cercano

closest_label <- function(x) {
  if (abs(x - qt1) < abs(x - qt2)) {
    return(1)
  } else if (abs(x - qt2) < abs(x - qt5)) {
    return(2)
  } else if (abs(x - qt5) < abs(x - qt10)) {
    return(5)
  } else if (abs(x - qt10) < abs(x - qt25)) {
    return(10)
  } else if (abs(x - qt25) < abs(x - qt50)) {
    return(25)
  } else if (abs(x - qt50) < abs(x - qt100)) {
    return(50)
  } else if (abs(x - qt100) < abs(x - qt200)) {
    return(100)
  } else {
    return(200)
  }
}

# Aplicar la función a la columna qq

serie.qma <- c(serie.qma$serie.qma)
tr <- sapply(serie.qma, closest_label)

# Crear el data frame con la columna qq y la columna de etiquetas basadas en la columna qq

df <- data.frame(
  year = 1:n,
  qma = serie.qma,
  qtr = tr
)

# Verificar el data frame
print(df)

barplot(tr, names.arg = years, ylim=c(0,30), 
        xlab="Time Step (year)", ylab="Return Period (year)",main="Caudales Máximos Anuales")
abline(h=10, col="red") # tr10 caudal de diseño para los escenarios con proyecto

# exportando df
nombre_archivo <- paste("serie", n, "years.csv", sep = "")
setwd("C:/JorgeR/SocioHidrologiaUdeC/AgentBasedModelling/Rcode_SeriesHidromet")
write.csv(df, nombre_archivo, row.names = FALSE)




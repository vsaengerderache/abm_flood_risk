# ..................................................................................................
# Proyecto de sociohidrologia de inundaciones
# Autor: Jorge Hurtado
# Fecha: 20 de junio de 2024
# Objetivos: Reformatear de ascii a tif los raster de nivel de agua obtenidos con IBER (escenarios con proyecto v3)
#    
# ..................................................................................................

# Preparing the environment

cat("\014")       #clean console   
rm(list = ls())   #clean environment
graphics.off()    #clean plots

# Libraries

library(sf)
library(dplyr)
library(terra)

# ..................................................................................................

# load raster

ruta_entrada <- "C:/JorgeR/SocioHidrologiaUdeC/IberHydraulicModelling/RasterDepth_ConProyecto_v3"
(archivos_raster <- list.files(path = ruta_entrada, pattern = "\\.asc$", full.names = TRUE))
(nraster <- length(archivos_raster))

rasterlist <- list()
rastername <- c()

for (i in 1:nraster){
  rasterlist[[i]] <- rast(archivos_raster[i])
  rastername[i] <- c(names(rasterlist[[i]]))
}

names(rasterlist) <- rastername

# names(rasterlist) # asigno nombres a los raster de la lista
# plot(rasterlist[[1]])

# resampling

# raster de referencia
ruta_referencia <- "C:/JorgeR/SocioHidrologiaUdeC/AgentBasedModelling/ReformatedRasterDepth/RasterDepth_SinProy10m/FloodQT1.tif"
r0 <- rast(ruta_referencia); r0

r1 <- rast(ext=r0, res=10, crs= crs(r0)); r1 # nuevo raster con valores vacios

ruta_salida <- "C:/JorgeR/SocioHidrologiaUdeC/AgentBasedModelling/ReformatedRasterDepth/RasterDepth_ConProy10m_v3"
resamplelist <- list()

for (i in 1:nraster){
  resamplelist[[i]] <- resample(rasterlist[[i]], r1)               # resampling
  nombre_archivo <- paste0(names(rasterlist[[i]]), ".tif")         # nombre de archivo con extension .tif
  ruta_archivo <- file.path(ruta_salida, nombre_archivo)           # ruta de archivo
  writeRaster(resamplelist[[i]], ruta_archivo, overwrite=TRUE)     # escritura de archivos raster resampleados
}

# nota: ahora los raster de "resamplelist", ubicados en "ruta_salida" tienen las mismas caracteristicas
# (dimensiones, resolucion, proyeccion, formato) que los raster usados de referencia (r0) 

# fin del script...

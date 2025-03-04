# model development (urban and development areas for ABM)

# Preparing the environment

cat("\014")       #clean console   
rm(list = ls())   #clean environment
graphics.off()    #clean plots

# libraries -----

#library(sf)
library(tidyverse)
library(terra)


# ..................................................................................................

# 1. Definir las coordenadas del área de estudio ----

study_area <- c(644000, 663000, 5865000, 5885000)

# Crear una función para generar y plotear el cuadrado

  area <- function(coords) {
  # Extraer las coordenadas
  x_min <- coords[1]
  x_max <- coords[2]
  y_min <- coords[3]
  y_max <- coords[4]
  
  # Crear las coordenadas del cuadrado
  square_coords <- matrix(c(x_min, y_min,
                            x_max, y_min,
                            x_max, y_max,
                            x_min, y_max,
                            x_min, y_min), 
                          ncol = 2, byrow = TRUE)
# Convertir a SpatVector
square <- vect(square_coords, type = "polygons")
return(square)

}

area <- area(study_area)
  
# Plotear el cuadrado
plot(area, main = "Cuadrado del Área de Estudio", col = "lightblue", border = "blue")


# 2. area urbana -----

setwd("C:/JorgeR/SocioHidrologiaUdeC/RespositorioInformacion/3_SH_Geoinfo/uso_de_suelo_2/")
cob <- vect("08_regi_n_de__biobio_actualizaci_n_2015.shp")
cob <- crop(cob, area) # cut by extend area 
plot(cob)

urban  <- subset(cob, cob$ID_USO == "01") # filter by urban use

ramad <- subset(urban,urban$ID==99652|urban$ID==105428)
plot(ramad)

caram <- subset(urban,urban$ID == 177372 | urban$ID == 1196420) 
plot(caram)

arau <- subset(urban, urban$ID == 176043 | urban$ID == 176810 | urban$ID == 176039 | urban$ID == 176040| urban$ID == 176033|urban$ID == 99653)
plot(arau)
 
urban3 <- union(ramad, caram) %>% union(.,arau)
plot(urban3)

# 3. area crecimiento ----

development <- cob %>% terra::subset(cob$ID_USO %in% c("01", "02", "03", "05")) #1 (areas urbanas e industriales); 2 (Terrenos agricolas); 3(Praderas y matorrales); 5(humedales) 
plot(development, col = 'lightgray')
plot(urban1,col="blue",add=T)

# guardando shapefiles
output_dir <- "C:/JorgeR/SocioHidrologiaUdeC/RespositorioInformacion/1_GeoinfoGenEdit/model_development3"

# Guardar 'area' como shapefile
writeVector(area, file.path(output_dir, "area3.shp"), overwrite = TRUE)

# Guardar 'urban' como shapefile
writeVector(urban3, file.path(output_dir, "urban3.shp"), overwrite = TRUE)

# Guardar 'development_area' como shapefile
writeVector(development, file.path(output_dir, "development3.shp"), overwrite = TRUE)

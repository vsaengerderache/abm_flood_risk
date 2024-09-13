# ..................................................................................................
# Proyecto de sociohidrologia de inundaciones
# Autor: Jorge Hurtado & Vicente Saenger
# Fecha: 12 de septiembre de 2024
#
# Objetivos: 1) Modelo de desarrollo de viviendas (se incorpora una flood_zone para seleccionar hogares en la zona de riesgo).
#            2) Modelo de reduccion de vulnerabilidad.
#            3) Modelo de reduccion de amenaza (se desactiva).
#            4) Calculo del riesgo (con 4 modelos de aprendizaje). 
#            5) Modelo Integrado.
#            6) Simulaciones Individuales.
#            7) Analisis de sensibilidad.
#            8) Implementacion de paralelizacion.
#    
# Notas: 1) Modelo con dos agentes: hogares (reduccion de vulnerabilidad) e instituciones (reduccion de amenaza).
#        2) Espacialmente, en cada periodo solo se agregan los agentes nuevos, los anteriores permanecen igual.
#        3) Para el mismo hogar se mantenie vc durante los periodos y se actualiza segun el modelo de aprendizaje.
#        4) Los escenarios de reduccion de amenaza son los raster CON proyecto (v3:modelacion en iber con defensas fluviales para Tr10). 
#        5) Zonas urbanas (urb) y de desarrollo (development area) obtenidos del uso de suelo 2015. Flood_zone es el flood map tr100.
#        6) Serie de QMA (Tr) obtenidos de caudales diarios (1970-1983) de DGA (cr2-camels). La serie que usa es temporal hasta definir la final.
#        7) Cambios respecto a la version anterior: 
#           I) se ajusta el area de los hogares con flood_zone; 
#           II) no se usa el escenario con proyecto (reduccion de amenaza); 
#           III) se usan 4 modelos de aprendizaje (MA);
#           IV) en una misma simulacion hay diferentes MA, sin embargo, el agente no cambia su MA.
#
# ..................................................................................................

# Preparing the environment

cat("\014")       #clean console   
rm(list = ls())   #clean environment
graphics.off()    #clean plots

# libraries -----

{
library(sf)
library(tidyverse)
library(terra)
library(multisensi)
library(furrr)
}

# ..................................................................................................

# series de caudal maximo anual ----

serie <- read.csv("Rcode_SeriesHidromet/Rc3_outputs/serie30years.csv", header = TRUE)
serie <- c(serie$qtr)

barplot(serie, main="Retorno QMA")
n <- length(serie)

# raster series
load_flood_raster_list <- function(serie, ruta_carpeta) {

  raster_depth_files <- list.files(path = ruta_carpeta, pattern = "\\.tif$", full.names = TRUE) # raster files
  
  floodlist <- raster_depth_files %>% 
    map(rast)
  floodname <- floodlist %>% 
    map_chr(names)
  names(floodlist) <- floodname
  
  # Definir los nombres de las capas de inundación
  flood_layers <- c("1" = "FloodQT1", "2" = "FloodQT2", "5" = "FloodQT5", 
                    "10" = "FloodQT10", "25" = "FloodQT25", "35" = "FloodQT35", 
                    "50" = "FloodQT50", "100" = "FloodQT100", "200" = "FloodQT200")
  
  # Asociar cada elemento de la serie con el raster correspondiente
floodserie <- map(serie, function(x) {
    layer_name <- flood_layers[as.character(x)]
    if (!is.null(layer_name) && layer_name %in% names(floodlist)) {
      return(floodlist[[layer_name]])
    } else {
      warning(paste("no raster for year:", x))
      return(NULL)
    }
  })
  return(floodserie)
}

floodserie_sp <- load_flood_raster_list(serie, "Rcode_Resampling/outputs_resampling/RasterDepth_SinProy10m_v1") # Load raster sin proyecto
floodserie_cp <- load_flood_raster_list(serie, "Rcode_Resampling/outputs_resampling/RasterDepth_ConProy10m_v3") # Load raster con proyecto

# study area

# poligonos de uso de suelo
area_urban <- terra::vect("Rcode_AreasDesarrollo/Rc3_outputs/urban3.shp")
area_development <-  terra::vect("Rcode_AreasDesarrollo/Rc3_outputs/development3.shp")

development_initial <- terra::aggregate(buffer(area_urban, 0)) ; #plot(development_initial)

# flood zone

{
flood_zone_rst <- load_flood_raster_list(100, "Rcode_Resampling/outputs_resampling/RasterDepth_SinProy10m_v1")
flood_zone_rst <- flood_zone_rst[[1]]
flood_zone_rst[] <- ifelse(flood_zone_rst[] > 0, 1, NA)
flood_zone_vct <- as.polygons(flood_zone_rst, values=FALSE, na.rm=TRUE)
flood_zone_vct <- simplifyGeom(flood_zone_vct, tolerance = 10)
}

# model: development ----

model_development <- function(development_initial, n, h, m, b) { 
  
  set.seed(1234)
  
  #n <- 30; h <- 1000; m <- 200; b <- 10

  list_hog_initial <- spatSample(development_initial, h, method = "random") ; #plot(list_hog_initial, pch=1,add= T)
  
  buf <- seq(b, (b*(n-1)), b); buf # crece b m por período
  list_buffer <- map(buf, ~ terra::aggregate(buffer(development_initial, .x)))

  # house generation per period within each polygon + buffer
  
  number_hog_new <- rep(m, n - 1)
  list_buffer <- map(list_buffer, ~ crop(.x, area_development))
  list_hog_new <- map2(list_buffer, number_hog_new, ~ spatSample(.x, .y, method = "random"))
  list_hog <- c(list(list_hog_initial), list_hog_new)  # join the initial list with the new lists.
  
  list_hog <- map(list_hog, ~ crop(.x, flood_zone_vct)) # uso de flood_zone para usar hogares con amenaza

  return(list_hog)
}

# model: reduction of vulnerability ----

model_reduction_vulnerability <- function(list_hog_model_development, pma, ua_1) { 
  
  set.seed(1234)
  
  #list_hog_model_development <- list_hog; ua_1 <- 0.5
  #pma <- c(0.97, 0.01, 0.01, 0.01)
  
  # add the will of accomplishment (vc) using map
  list_hog_model_development <- map(list_hog_model_development, ~ {
  #sf_hog <- st_as_sf(list_hog[[2]])
  sf_hog <- st_as_sf(.x)  # Convertir a sf
  len <- length(sf_hog[[1]])  # Número de hogares en cada período
  vc <- runif(len)  # Voluntad de cumplimiento (aleatorio entre 0 y 1)
  sf_hog <- data.frame(sf_hog)  # Convertir a data.frame
  sf_hog$vc <- vc  # Agregar el campo de voluntad de cumplimiento
  #sf_hog  # Devolver el dataframe actualizado
  
  # asignacion del modelo de aprendizaje: 1, 2, 3, 4
  # Verificar que las proporciones sumen 1
  if (sum(pma) != 1) {
    stop("Las proporciones en 'pma' deben sumar 1")
  }
  
  n_reps <- pma * len                 
  
  # Redondear hacia abajo y calcular la diferencia
  n_reps_floor <- floor(n_reps)
  sobra <- len - sum(n_reps_floor)
  
  # Distribuir la diferencia aleatoriamente entre los valores
  indices_no_cero <- which(pma != 0)
  ajuste <- sample(indices_no_cero, size = sobra) #solo usar donde pma no es 0
  
  # Aumentar en uno los valores seleccionados para ajustar el total
  n_reps_final <- n_reps_floor
  n_reps_final[ajuste] <- n_reps_final[ajuste] + 1
  non_zero_reps <- which(n_reps_final != 0)
  
  n_reps_final <- n_reps_final[non_zero_reps]
  
  # Crear el vector con valores repetidos
  #values <- rep(indices_no_cero, times = n_reps_final)  # vector con valores repetidos
  values <- rep(non_zero_reps, times = n_reps_final)  # vector con valores repetidos
  values <- sample(values)
  
  sf_hog$ma <- values 

  sf_hog  # Devolver el dataframe actualizado
  
  })
  
  # evaluating whether agents adopt measures
  list_hog_model_vulnerability <- map(list_hog_model_development, ~ mutate(.x, lev = ifelse(vc > ua_1, "si", "no")))  
  
  list_hog_model_vulnerability <- list_hog_model_vulnerability %>%
    accumulate(., rbind) # bind hog
  
  return(list_hog_model_vulnerability)
}

# model: reduction of hazard ----

model_reduction_hazard <- function(list_hog_model_development, ua_2) { 
  
  #list_hog_model_development <- list_hog_model_vulnerability; ua_2 <- 100
  
  alarm <- map_int(serie, ~ ifelse(.x >= ua_2, 1, 0)) # build alarm serie
  alarm_position <- which(alarm == 1)[1] # find alarm position
  
  # build serie with 1 from alarm position and 0 befor
  serie_infrastructure <- case_when(
    is.na(alarm_position) ~ rep(0, length(alarm)),
    TRUE ~ ifelse(seq_along(alarm) >= (alarm_position + 1), 1, 0)
  )
  
  # extract water level for each house
  list_hog_model_hazard <- map2(list_hog_model_development, seq_along(serie_infrastructure), function(hog, i) {
    depth <- if (serie_infrastructure[i] == 0) {
      terra::extract(floodserie_sp[[i]], vect(hog$geometry))
    } else {
      terra::extract(floodserie_cp[[i]], vect(hog$geometry))
    }
    depth[is.na(depth)] <- 0
    hog$depth <- depth
    return(hog)
  })
  
  return(list_hog_model_hazard)
}

# model: risk 1 status quo ----

model_risk_1 <- function(list_hog_model_vulnerability_hazard, ua_1, f) { 
  
  #list_hog_model_vulnerability_hazard <- list_hog_model_hazard; ua_1=0.5; f <- 0.5
  
  list_hog_risk1 <- list_hog_model_vulnerability_hazard
  
  #aqui debo filtrar por ma
  
  # Filtrar la lista y luego mutar para agregar la columna 'inun'
  list_hog_risk1 <- map(list_hog_risk1, function(df) {
    df %>%
      filter(ma == 1) %>%  # Filtrar donde 'ma' es igual a 1
      mutate(inun = case_when(
        depth$Flood > f ~ "si",                # Se inunda en todos los casos
        depth$Flood == 0 ~ "no",               # No se inunda en ningún caso
        depth$Flood < f & lev == "no" ~ "si",  # Si la inundación es menos de f y no se levanta, se inunda
        depth$Flood < f & lev == "si" ~ "no",  # Si la inundación es menos de f y si se levanta, no se inunda
        TRUE ~ "NAN"                           # En cualquier otro caso
      ))
  })
  
  return(list_hog_risk1)
  
  }

# model: risk 2 reactivo ----

model_risk_2 <- function(list_hog_model_vulnerability_hazard, ua_1, f, ap) { 
  
  #list_hog_model_vulnerability_hazard <- list_hog_model_hazard; ua_1=0.7; f <- 0.5; ap <- 0.2
  
  list_hog_risk2 <- list_hog_model_vulnerability_hazard
  
  # Filtrar la lista y luego mutar para agregar la columna 'inun'
  list_hog_risk2 <- map(list_hog_risk2, function(df) {
    df %>%
      filter(ma == 2) %>%  # Filtrar donde 'ma' es igual a 1
      mutate(inun = case_when(
        depth$Flood > f ~ "si",                # Se inunda en todos los casos
        depth$Flood == 0 ~ "no",               # No se inunda en ningún caso
        depth$Flood < f & lev == "no" ~ "si",  # Si la inundación es menos de f y no se levanta, se inunda
        depth$Flood < f & lev == "si" ~ "no",  # Si la inundación es menos de f y si se levanta, no se inunda
        TRUE ~ "NAN"                           # En cualquier otro caso
      ))
  })
  
  # incorporates learning
  
  # agregar aqui: si la lista no esta vacia ejecutar learning caso devolver lista vacia
  # nota esta parte es necesaria porque si pma tiene una fraccion pequeña y h y m son pequeños pueden devolver listas vacias para algun ma, esta lista vacia no ejecuta los modelos de riesgo de esta parte
  
  list_hog_risk2 <- list_hog_risk2 #??
  for (i in 2:n) {
    list_hog_risk2[[i]]$hog_exist <- rep(FALSE, nrow(list_hog_risk2[[i]])) 
    list_hog_risk2[[i]]$hog_exist[1:nrow(list_hog_risk2[[i-1]])] <- TRUE # Marcar como TRUE los hogares existentes
  }
  
  for (i in 2:n) {
    
    list_hog_risk2[[i]]$vc   <- ifelse(list_hog_risk2[[i]]$hog_exist, list_hog_risk2[[i-1]]$vc, list_hog_risk2[[i]]$vc)
    list_hog_risk2[[i]]$inun_prev <- ifelse(list_hog_risk2[[i]]$hog_exist, list_hog_risk2[[i-1]]$inun, "no") # un nuevo campo que tiene la condicion de inundacion del periodo anterior(i-1) para el mismo punto
    list_hog_risk2[[i]]$vc   <- ifelse(list_hog_risk2[[i]]$inun_prev =="si", list_hog_risk2[[i]]$vc + ap, 
                                      ifelse((list_hog_risk2[[i]]$inun_prev =="no"), list_hog_risk2[[i]]$vc, NaN))                              # si este nuevo campo dice que si hubo inundacion en el periodo anterior le suma ap a vc, en un nuevo campo
    
    list_hog_risk2[[i]]$vc <- pmax(0, pmin(1, list_hog_risk2[[i]]$vc)) # restringe vc dentro del rango 0 y 1
    
    list_hog_risk2[[i]]$lev  <- ifelse(list_hog_risk2[[i]]$vc > ua_1, "si", "no")
  }
  
  list_hog_risk2 <- map(list_hog_risk2, function(df) {
    df %>%
      mutate(inun = case_when(
        depth$Flood > f ~ "si",               # se inunda en todos los casos
        depth$Flood == 0 ~ "no",              # no se inunda en ningun caso
        depth$Flood < f & lev == "no" ~ "si", # si la inundacion es menos de f y no se levanta se inunda
        depth$Flood < f & lev == "si" ~ "no", # si la inundacion es menos de f y si se levanta no se inunda
        TRUE ~ "NAN"
      ))
  })
  
  return(list_hog_risk2)
}

# model: risk 3 olvido ----

model_risk_3 <- function(list_hog_model_vulnerability_hazard, ua_1, f, ap, r) { 
  
  #list_hog_model_vulnerability_hazard <- list_hog_model_hazard; ua_1=0.7; f <- 0.5; ap <- 0.2; r <- 10
  
  list_hog_risk3 <- list_hog_model_vulnerability_hazard
  
  # Filtrar la lista y luego mutar para agregar la columna 'inun'
  list_hog_risk3 <- map(list_hog_risk3, function(df) {
    df %>%
      filter(ma == 3) %>%  # Filtrar donde 'ma' es igual a 1
      mutate(inun = case_when(
        depth$Flood > f ~ "si",                # Se inunda en todos los casos
        depth$Flood == 0 ~ "no",               # No se inunda en ningún caso
        depth$Flood < f & lev == "no" ~ "si",  # Si la inundación es menos de f y no se levanta, se inunda
        depth$Flood < f & lev == "si" ~ "no",  # Si la inundación es menos de f y si se levanta, no se inunda
        TRUE ~ "NAN"                           # En cualquier otro caso
      ))
  })
  
  # incorporates learning
  #list_hog_risk <- list_hog_risk #??
  for (i in 2:n) {
    list_hog_risk3[[i]]$hog_exist <- rep(FALSE, nrow(list_hog_risk3[[i]])) 
    list_hog_risk3[[i]]$hog_exist[1:nrow(list_hog_risk3[[i-1]])] <- TRUE # Marcar como TRUE los hogares existentes
  }
  
  for (i in 2:n) {
    
    list_hog_risk3[[i]]$vc   <- ifelse(list_hog_risk3[[i]]$hog_exist, list_hog_risk3[[i-1]]$vc, list_hog_risk3[[i]]$vc)
    list_hog_risk3[[i]]$inun_prev <- ifelse(list_hog_risk3[[i]]$hog_exist, list_hog_risk3[[i-1]]$inun, "no") # un nuevo campo que tiene la condicion de inundacion del periodo anterior(i-1) para el mismo punto
    list_hog_risk3[[i]]$vc   <- ifelse(list_hog_risk3[[i]]$inun_prev =="si", list_hog_risk3[[i]]$vc + ap, 
                                      ifelse((list_hog_risk3[[i]]$hog_exist & list_hog_risk3[[i]]$inun_prev =="no"), 
                                             list_hog_risk3[[i]]$vc - (ap/r), list_hog_risk3[[i]]$vc))                              # si este nuevo campo dice que si hubo inundacion en el periodo anterior le suma ap a vc, en un nuevo campo
    
    list_hog_risk3[[i]]$vc <- pmax(0, pmin(1, list_hog_risk3[[i]]$vc)) # restringe vc dentro del rango 0 y 1
    
    list_hog_risk3[[i]]$lev  <- ifelse(list_hog_risk3[[i]]$vc > ua_1, "si", "no")
  }
  
  list_hog_risk3 <- map(list_hog_risk3, function(df) {
    df %>%
      mutate(inun = case_when(
        depth$Flood > f ~ "si",               # se inunda en todos los casos
        depth$Flood == 0 ~ "no",              # no se inunda en ningun caso
        depth$Flood < f & lev == "no" ~ "si", # si la inundacion es menos de f y no se levanta se inunda
        depth$Flood < f & lev == "si" ~ "no", # si la inundacion es menos de f y si se levanta no se inunda
        TRUE ~ "NAN"
      ))
  })
  
  return(list_hog_risk3)
}

# model: risk 4 proactivo ----

model_risk_4 <- function(list_hog_model_vulnerability_hazard, ua_1, f, ap, r) { 
  
  #list_hog_model_vulnerability_hazard <- list_hog_model_hazard; ua_1=0.7; f <- 0.5; ap <- 0.2; r <- 10
  
  list_hog_risk4 <- list_hog_model_vulnerability_hazard
  
  # Filtrar la lista y luego mutar para agregar la columna 'inun'
  list_hog_risk4 <- map(list_hog_risk4, function(df) {
    df %>%
      filter(ma == 4) %>%  # Filtrar donde 'ma' es igual a 1
      mutate(inun = case_when(
        depth$Flood > f ~ "si",                # Se inunda en todos los casos
        depth$Flood == 0 ~ "no",               # No se inunda en ningún caso
        depth$Flood < f & lev == "no" ~ "si",  # Si la inundación es menos de f y no se levanta, se inunda
        depth$Flood < f & lev == "si" ~ "no",  # Si la inundación es menos de f y si se levanta, no se inunda
        TRUE ~ "NAN"                           # En cualquier otro caso
      ))
  })
  
  # incorporates learning
  #list_hog_risk <- list_hog_risk #??
  for (i in 2:n) {
    list_hog_risk4[[i]]$hog_exist <- rep(FALSE, nrow(list_hog_risk4[[i]])) 
    list_hog_risk4[[i]]$hog_exist[1:nrow(list_hog_risk4[[i-1]])] <- TRUE # Marcar como TRUE los hogares existentes
  }
  
  for (i in 2:n) {
    
    list_hog_risk4[[i]]$vc   <- ifelse(list_hog_risk4[[i]]$hog_exist, list_hog_risk4[[i-1]]$vc, list_hog_risk4[[i]]$vc)
    list_hog_risk4[[i]]$inun_prev <- ifelse(list_hog_risk4[[i]]$hog_exist, list_hog_risk4[[i-1]]$inun, "no") # un nuevo campo que tiene la condicion de inundacion del periodo anterior(i-1) para el mismo punto
    list_hog_risk4[[i]]$vc   <- ifelse(list_hog_risk4[[i]]$inun_prev =="si", list_hog_risk4[[i]]$vc + ap, 
                                      ifelse((list_hog_risk4[[i]]$hog_exist & list_hog_risk4[[i]]$inun_prev =="no"), 
                                             list_hog_risk4[[i]]$vc + (ap/r), list_hog_risk4[[i]]$vc))  #                            # si este nuevo campo dice que si hubo inundacion en el periodo anterior le suma ap a vc, en un nuevo campo
    
    list_hog_risk4[[i]]$vc <- pmax(0, pmin(1, list_hog_risk4[[i]]$vc)) # restringe vc dentro del rango 0 y 1
    
    list_hog_risk4[[i]]$lev  <- ifelse(list_hog_risk4[[i]]$vc > ua_1, "si", "no")
  }
  
  list_hog_risk4 <- map(list_hog_risk4, function(df) {
    df %>%
      mutate(inun = case_when(
        depth$Flood > f ~ "si",               # se inunda en todos los casos
        depth$Flood == 0 ~ "no",              # no se inunda en ningun caso
        depth$Flood < f & lev == "no" ~ "si", # si la inundacion es menos de f y no se levanta se inunda
        depth$Flood < f & lev == "si" ~ "no", # si la inundacion es menos de f y si se levanta no se inunda
        TRUE ~ "NAN"
      ))
  })
  
  return(list_hog_risk4)
}

# model: integrated 1 ----

model_integrated <- function(development_initial, n, h, m, b, pma, ua_1, ua_2, f, ap, r) { 
  
  #n <- 30; h <- 500; m <- 100; b <- 10 ; ua_1 <- 0.7; ua_2 <- 100 ; f <- 0.5 ; ap <- 0.3 ; r <- 7 
  #pma <- c(0.97, 0.01, 0.01, 0.01)
  
  # Generar la lista base hasta el paso model_reduction_hazard
  list_base <- development_initial %>%
    model_development(., n, h, m, b) %>%
    model_reduction_vulnerability(., pma, ua_1) %>% 
    model_reduction_hazard(., ua_2)
  
  # Aplicar cada model_risk sobre la lista base
  # list_risk_1 <- model_risk_1(list_base, ua_1, f)
  # list_risk_2 <- model_risk_2(list_base, ua_1, f, ap)
  # list_risk_3 <- model_risk_3(list_base, ua_1, f, ap, r)
  # list_risk_4 <- model_risk_4(list_base, ua_1, f, ap, r)
  
  if (pma[1] > 0) {
    # Si es mayor que 0, crea list_risk_1 utilizando model_risk_1
    list_risk_1 <- model_risk_1(list_base, ua_1, f)
  } else {
    # Si no es mayor que 0, crea una lista vacía
    list_risk_1 <- vector("list", n)
  } 
  
  if (pma[2] > 0) {
    # Si es mayor que 0, crea list_risk_1 utilizando model_risk_1
    list_risk_2 <- model_risk_2(list_base, ua_1, f,ap)
  } else {
    # Si no es mayor que 0, crea una lista vacía
    list_risk_2 <- vector("list", n)
  }
  
  if (pma[3] > 0) {
    # Si es mayor que 0, crea list_risk_1 utilizando model_risk_1
    list_risk_3 <- model_risk_3(list_base, ua_1, f,ap,r)
  } else {
    # Si no es mayor que 0, crea una lista vacía
    list_risk_3 <- vector("list", n)
  }
  
  if (pma[4] > 0) {
    # Si es mayor que 0, crea list_risk_1 utilizando model_risk_1
    list_risk_4 <- model_risk_4(list_base, ua_1, f,ap,r)
  } else {
    # Si no es mayor que 0, crea una lista vacía
    list_risk_4 <- vector("list", n)
  }
  
  # Unificar las listas en una sola, intercalando los elementos
  #list_hog_inter <- pmap(list(list_risk_1, list_risk_2, list_risk_3, list_risk_4), list)
  
  list_hog_inter <- map(seq_len(n), function(i) {
    bind_rows(
      list_risk_1[[i]],
      list_risk_2[[i]],
      list_risk_3[[i]],
      list_risk_4[[i]]
    )
  })
  
  # Aplicar la función de riesgo basada en el valor de 'ma'
  
  # flooded houses per year
  hoginun <- map_dbl(list_hog_inter, ~ sum(.x$inun == "si"))
  lenlist <- map_dbl(list_hog_inter, ~ nrow(.x))
  hog_flooded <- t(data.frame((hoginun/lenlist)*100))
  
  return(hog_flooded)
}

# test 1 ----
# run one model integrated

n <- 30; h <- 1603; m <- 220; b <- 10 ; ua_1 <- 0.5; ua_2 <- 100 ; f <- 0.5 ; ap <- 0.3 ; r <- 7 
pma <- c(0.25, 0.25, 0.25, 0.25)
#pma <- c(0.97, 0.01, 0.01, 0.01)

#t_beg <- Sys.time()
test1 <- model_integrated(development_initial, n, h, m, b, pma, ua_1, ua_2, f, ap, r)
#t_end <- Sys.time()
#test_time <- t_end - t_beg; test_time
plot(seq(1:n), test1, type = "l", xlab = "Años", ylab = "% hogares afectados", ylim = c(0, 100), lty=2, lwd=2.5, col = "red")

# test 2

#n <- 30; h <- 1000; m <- 200; b <- 10 ; ua_1 <- 0.5; ua_2 <- 100 ; f <- 0.5 ; ap <- 0.3 ; r <- 7 
pma <- c(1, 0, 0, 0)
test2 <- model_integrated(development_initial, n, h, m, b, pma, ua_1, ua_2, f, ap, r)
lines(seq(1:n), test2, type = "l", xlab = "Años", ylab = "% hogares afectados", ylim = c(0, 100), lwd=2.5,col = "blue")

# test 3

#n <- 30; h <- 500; m <- 100; b <- 10 ; ua_1 <- 0.7; ua_2 <- 100 ; f <- 0.5 ; ap <- 0.3 ; r <- 7 
pma <- c(0, 1, 0, 0)
test3 <- model_integrated(development_initial, n, h, m, b, pma, ua_1, ua_2, f, ap, r)
lines(seq(1:n), test3, type = "l", xlab = "Años", ylab = "% hogares afectados", ylim = c(0, 100), lwd=2.5,col = "gray")

# test 4

#n <- 30; h <- 500; m <- 100; b <- 10 ; ua_1 <- 0.7; ua_2 <- 100 ; f <- 0.5 ; ap <- 0.3 ; r <- 7 
pma <- c(0, 0, 1, 0)
test4 <- model_integrated(development_initial, n, h, m, b, pma, ua_1, ua_2, f, ap, r)
lines(seq(1:n), test4, type = "l", xlab = "Años", ylab = "% hogares afectados", ylim = c(0, 100), lwd=2.5,col = "green")

# test 5

#n <- 30; h <- 500; m <- 100; b <- 10 ; ua_1 <- 0.7; ua_2 <- 100 ; f <- 0.5 ; ap <- 0.3 ; r <- 7 
pma <- c(0, 0, 0, 1)
test5 <- model_integrated(development_initial, n, h, m, b, pma, ua_1, ua_2, f, ap, r)
lines(seq(1:n), test5, type = "l", xlab = "Años", ylab = "% hogares afectados", ylim = c(0, 100), lwd=2.5,col = "black")

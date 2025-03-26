# ..................................................................................................
# Proyecto de sociohidrologia de inundaciones
# Autor: Jorge Hurtado & Vicente Saenger
# Fecha: 24 de febrero de 2025
#
# Objetivos: 1) Modelo de desarrollo de viviendas (se corta con flood_zone TR100 para seleccionar hogares en la zona de riesgo).
#            2) Modelo de reduccion de vulnerabilidad.
#            3) Modelo de reduccion de amenaza (se desactiva).
#            4) Calculo del riesgo (con 4 modelos de aprendizaje). 
#            5) Modelo Integrado.
#            6) Simulaciones Individuales.
#            7) Plots con estadisticos globales (preparados/afectados) de calibracion (2017) y validacion (2024)
#            8) Metricas de calibracion/validacion espacial.
#            
# Notas:  1) Espacialmente, en cada periodo solo se agregan los agentes nuevos, los anteriores permanecen igual.
#         2) Para el mismo hogar se mantenie vc durante los periodos y se actualiza segun el modelo de aprendizaje.
#         3) Flood_zone es el flood map tr100.
#         4) Serie de QMA rellenada con SWAT. El Tr se aproxima al inmediato superior.
#         5) Aunque esta implementado, NO se usa el escenario con proyecto (ua_2: umbral para reduccion de amenaza). 
#         6) Se implementaron 4 modelos de aprendizaje (MA). Sin embargo para calibracion/validacion solo se usa Status quo y Olvido
#         7) En una misma simulacion hay diferentes MA, sin embargo, el agente no cambia su MA durante la simulación.
#
# ..................................................................................................

# Preparing the environmentt

cat("\014")       #clean console   
rm(list = ls())   #clean environment
graphics.off()    #clean plots

# glossary ----

# n: number of years
# n_switch: year of switch area of development 
# h: number of houses initial 
# m: number of new houses per year
# pma: proportion of behaviors at household level
# ua_1: threshold at household level at which flood adaptive measures are taken 
# ua_2: threshold at state level in function of flood return period 
# f: threshold at household level that damage of flood depth can be reduced by adaptive measures
# ap: learning rate 
# r: years in which what has been learned is forgotten  

# libraries -----

{
library(tidyverse)
library(readxl)
library(sf)
library(terra)
library(multisensi)
library(furrr)
library(gganimate)
library(gifski)
library(tidyterra)
library(RColorBrewer)
library(ggrepel)
library(gridExtra) 
library(mapview)
library(tmap)
library(spatstat.geom) # as.ppp in NNI
library(dbscan) #clustering for NNI
library(spatstat.explore) # Hopkins-Skellam test
library(nngeo) # k-Nearest Neighbor Join for Spatial Data
library(caret) # confusion matrix
}

# ..................................................................................................
# path working directory -----
setwd("C:/JorgeR/SocioHidrologiaUdeC/AgentBasedModelling/data")

# set up ----
year_initial <- 1985
year_end <- 2024
# n <- year_end - year_initial + 1

# discharge series
data <- read.csv("./hydroseries/serie40years.csv", header = TRUE) #%>%
  # dplyr::filter(year %in% seq(year_initial, year_end)) %>%
  # dplyr::filter(serie == "serie_4")
serie <- data$qtr
years <- data$year
barplot(serie, main = "Retorno QMA con serie SWAT", names.arg = years, las = 2, cex.names = 0.8)

n <- length(serie)

# households growth
hog_initial <- read.csv("./hog_growth.csv", header = TRUE) %>%
  dplyr::filter(year == year_initial, localidad == "total") %>%
  dplyr::select(hogares) %>%
  pull() %>% 
  first() %>%
  as.numeric()  
hog_growth <- read.csv("./hog_growth.csv", header = TRUE) %>%
  dplyr::filter(year == year_initial, localidad == "total") %>%
  dplyr::select(growth_rate) %>%
  pull() %>% 
  first() %>%
  as.numeric()  

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


floodserie_sp <- load_flood_raster_list(serie, "RasterDepth_SinProy10m_v1") # Load raster sin proyecto
floodserie_cp <- load_flood_raster_list(serie, "RasterDepth_ConProy10m_v3") # Load raster con proyecto

# development area
development_area_list <- list(
  terra::vect("model_development4/1997_urban.shp") %>% terra::aggregate(),
  terra::vect("model_development4/2008_urban.shp") %>% terra::aggregate(),
  terra::vect("model_development4/2015_urban.shp") %>% terra::aggregate()
)

# flood zone
{
flood_zone_rst <- load_flood_raster_list(100, "RasterDepth_SinProy10m_v1")
flood_zone_rst <- flood_zone_rst[[1]]
flood_zone_rst[] <- base::ifelse(flood_zone_rst[] > 0, 1, NA)
flood_zone_vct <- terra::as.polygons(flood_zone_rst, values=FALSE, na.rm=TRUE)
flood_zone_vct <- terra::simplifyGeom(flood_zone_vct, tolerance = 10)
#flood_zone_vct <- buffer(flood_zone_vct, width = 50)
}

#setwd("C:\\010_r\\project_sociohydro_abm_flood_risk_r_model_integrated")
#flood_zone_vct <- terra::vect("survey_2017\\survey_polygon.shp") %>%
#  sf::st_as_sf()
flood_zone_vct <- sf::st_as_sf(flood_zone_vct)

# survey 2017 households
crs = "EPSG:32718"

survey_2017 <- read_excel("./surveys/survey2017_edit.xlsx", sheet = 3, col_names = T) %>%
  na.omit(.) %>%
  rename_with(tolower) %>%
  sf::st_as_sf(., coords = c("coord_x", "coord_y"), crs = crs) %>% 
  sf::st_intersection(flood_zone_vct)

sum_2017 <- survey_2017 %>%
  #dplyr::filter(afectados == "si") %>%
  dplyr::summarise(
    total = n(), # total no inun
    inun = sum(str_detect(bienes, "vivienda") & danos == "moderados" | danos == "graves"),
    
    porc_inun = (inun / total) * 100,
    med = sum(medidas_proteccion == "si"),
    porc_med = (med / total) * 100,
  )

sum_2017$porc_inun
sum_2017$porc_med

# survey 2024 households
survey_2024 <- read_excel("./surveys/survey2024_edit.xlsx", col_names = T) %>%
  rename_with(tolower) %>%
  filter(comuna == "Arauco") %>%
  sf::st_as_sf(., coords = c("lon_encuesta", "lat_encuesta"), crs = 4326) %>%
  sf::st_transform(crs) %>%
  sf::st_intersection(flood_zone_vct)

#mapview::mapView(survey_2024)

sum_2024 <- survey_2024 %>% 
  #dplyr::filter(c10=="Propia") %>% # can be actived
  dplyr::select(d19, d20,contains("g36_1"), contains("g36_2")) %>%       #D.19. Con relación al último evento de inundación ¿Entró agua a su hogar? . #D.20. Y respecto a ese último evento de inundación ¿Sufrió daño material en su vivienda? 
                     
  #mutate(measure = coalesce(g36_1_1, g36_1_2, g36_1_3, g36_1_6, g36_1_7,g36_1_11, g36_1_12)) %>%  #medidas en relacion a la ultima inundacion (2024), 
  mutate(measure = coalesce(g36_1_1, g36_1_2, g36_1_3, g36_1_6, g36_1_7,g36_1_11, g36_1_12,
                            g36_2_1, g36_2_2, g36_2_3, g36_2_6, g36_2_7,g36_2_11, g36_2_12)) %>%
  
  dplyr::summarise(
    total = n(), # total no inun
    inun = sum(d19 == "Sí" & d20 == "Sí"), # afected
    #inun = sum(d19 == "Sí"), # afected
    porc_inun = (inun / total) * 100, # percent
    med = sum(!is.na(measure)), # measure
    porc_med = (med / total) * 100 # percent
  ) #%>% .$porcentaje_inun

# model: development ----
model_development <- function(development_area_list, n_switch, n, h, m, pma) {
  #n <- n; n_switch <- c(10,25); h <- hog_initial; m <- hog_growth; b <- 10 ; ua_1 <- 0.7; ua_2 <- 100 ; f <- 0.5 ; ap <- 0.3 ; r <- 7 
  #pma <- c(1, 0, 0, 0)
  
  set.seed(1234)
  
  # intervals for development areas
  area_intervals <- c(n_switch-1, n) # Append n as the final limit
  land_use_layers <- data.frame(
    area = seq_along(development_area_list),
    start = c(1, area_intervals[-length(area_intervals)] + 1),
    end = area_intervals
  )
  
  # create random point in development area list
  list_hog <- purrr::map2(seq_len(n), c(h, rep(m, n - 1)), 
                          ~ {
                            current_layer <- if (.x == 1) {
                              1  
                            } else {
                              land_use_layers$area[.x >= land_use_layers$start & .x <= land_use_layers$end]
                            }
                            terra::spatSample(development_area_list[[current_layer]], .y, method = "random")
                          }
  )
  
  # filter by flood zone
  list_hog <- purrr::map(list_hog, ~ terra::crop(.x, vect(flood_zone_vct)))
  
  # create a data frame for houses
  list_hog <- purrr::imap(list_hog, ~ {
    sf_hog <- sf::st_as_sf(.x) %>%
      data.frame() %>%
      dplyr::mutate(id = seq_len(n()) + ifelse(.y == 1, 0, sum(sapply(list_hog[1:(.y-1)], nrow)))) %>%
      dplyr::mutate(year_entry = .y) %>%
      dplyr::mutate(ma = sample(c(1, 2, 3, 4), n(), replace = TRUE, prob = pma)) %>%
    return(sf_hog)
  })
  
  # accumulate houses
  list_hog <- list_hog %>%
    purrr::accumulate(., base::rbind)
  
  # add year and existing houses column
  list_hog <- purrr::imap(list_hog, ~ {
    .x %>%
      dplyr::mutate(year = .y) %>%
      dplyr::mutate(hog_exist = .y == 1 | dplyr::row_number() <= nrow(list_hog[[max(1, .y - 1)]])) %>%
      dplyr::select(id, year, year_entry, geometry, ma, hog_exist)
  })
  
  return(list_hog)
}

# model: reduction of vulnerability ----
model_reduction_vulnerability <- function(list_hog_model_development, ua_1) { 
  #list_hog_model_development <- list_hog_inter; ua_1 <- 0.5
  
  # add vc
  list_hog_model_development <- list_hog_model_development %>%
    map(as_tibble) %>%
    reduce(bind_rows) %>%
    group_by(id) %>%
    dplyr::mutate(vc = unique(runif(1, 0, 1))) %>%
    dplyr::select(id, year, year_entry, geometry, ma, vc, hog_exist) %>%
    base::split(.$year)
    
  # evaluating whether agents adopt measures
  list_hog_model_vulnerability <- purrr::map(list_hog_model_development, ~ dplyr::mutate(.x, lev = base::ifelse(vc > ua_1, TRUE, FALSE)))  
  
  return(list_hog_model_vulnerability)
}
model_reduction_vulnerability <- function(list_hog_model_development, ua_1) { 
  #list_hog_model_development <- list_hog_inter; ua_1 <- 0.5
  
  # add vc
  vc_values <- list_hog_model_development %>%
    map(as_tibble) %>%
    reduce(bind_rows) %>%
    distinct(id) %>%
    mutate(vc = runif(n(), 0, 1))  # Generar valores aleatorios únicos por id
  
  # Unir valores precomputados y procesar el resto
  list_hog_model_development <- list_hog_model_development %>%
    map(as_tibble) %>%
    reduce(bind_rows) %>%
    left_join(vc_values, by = "id") %>%  # Añadir valores aleatorios únicos
    dplyr::select(id, year, year_entry, geometry, ma, vc, hog_exist) %>%
    base::split(.$year)
  
  return(list_hog_model_development)
}

# model: reduction of hazard ----
model_reduction_hazard <- function(list_hog_model_development, ua_2) { 
  #list_hog_model_development <- list_hog_inter; ua_2 <- 100
  
  alarm <- map_int(serie, ~ ifelse(.x >= ua_2, 1, 0)) # build alarm serie
  alarm_position <- base::which(alarm == 1)[1] # find alarm position
  
  # build serie with 1 from alarm position and 0 befor
  serie_infrastructure <- dplyr::case_when(
    base::is.na(alarm_position) ~ base::rep(0, length(alarm)),
    TRUE ~ base::ifelse(seq_along(alarm) >= (alarm_position + 1), 1, 0))
  
  # extract water level for each house
  list_hog_model_hazard <- map2(list_hog_model_development, seq_along(serie_infrastructure), function(hog, i) {
    depth <- if (serie_infrastructure[i] == 0) {
      terra::extract(floodserie_sp[[i]], vect(hog$geometry), ID = FALSE)
    } else {
      terra::extract(floodserie_cp[[i]], vect(hog$geometry), ID = FALSE)
    }
    depth[is.na(depth)] <- 0
    hog$return_period <- serie[[i]]
    hog$depth <- depth[[1]]
    return(hog)
  })
  
  return(list_hog_model_hazard)
}

# model: risk ----
model_risk <- function(list_hog_model_vulnerability_hazard, ua_1, f, ap, r) { 
  #list_hog_model_vulnerability_hazard <- list_hog_inter

  df_hog_risk <- list_hog_model_vulnerability_hazard %>%
    bind_rows() 
  list_hog_risk <- list()
  for (id_unique in unique(df_hog_risk$id)) {
    
    df <- df_hog_risk %>%
      filter(id == id_unique) %>%
      arrange(year)
    
    for (i in seq_len(nrow(df))) {
      
      df$lev[i] = df$vc[i] > ua_1    # Define si el hogar se "levanta" o no
      df$inun[i] = case_when(
        df$depth[i]  > f ~ TRUE,
        df$depth[i]  == 0 ~ FALSE,
        df$depth[i]  < f & !df$lev[i]  ~ TRUE,
        df$depth[i]  < f & df$lev[i]  ~ FALSE
      )
      if (i <= nrow(df) - 1) { 
      df$vc[i + 1] = case_when(
        # status quo
        df$ma[i] == 1 ~ df$vc[i],                        
        # reactivo
        df$ma[i] == 2 &  df$inun[i] ~ df$vc[i] + ap, 
        df$ma[i] == 2 & !df$inun[i] ~ df$vc[i],                        
        # olvido
        df$ma[i] == 3 &  df$inun[i] ~ df$vc[i] + ap,                        
        df$ma[i] == 3 & !df$inun[i] ~ df$vc[i] - ap/r,  
        # apredizaje
        df$ma[i] == 4 &  df$inun[i] ~ df$vc[i] + ap,  
        df$ma[i] == 4 & !df$inun[i] ~ df$vc[i] + ap/r )
      df$vc[i + 1] = pmax(0, pmin(1, df$vc[i + 1]))
      }
    }
    
    df$inun_once <- cumany(df$inun)
    list_hog_risk <- c(list_hog_risk, list(df))
  }
  
  list_hog_risk <- list_hog_risk %>%
    bind_rows() %>%
    split(.$year)
  return(list_hog_risk)
}

# model: integrated ----
model_integrated <- function(development_area_list, n_switch, n, h, m, pma, ua_1, ua_2, f, ap, r) { 
  #n <- n; n_switch <- c(10,25); h <- hog_initial; m <- hog_growth; b <- 10 ; ua_1 <- 0.7; ua_2 <- 100 ; f <- 0.5 ; ap <- 0.3 ; r <- 7 
  #pma <- c(0.25, 0.25, 0.25, 0.25)
  
  # verification
  if (base::sum(pma) != 1) {stop("The proportions in 'pma' must sum 1")}
  if (ua_1 >= 1) {stop("Threshold at household level must be less than 1")}
  if (ua_1 <  0) {stop("Threshold at household level must be at less 0")}
  if (m >= 0.1 * h) {stop("Very high growth")}
  
  # run model_integrated
  list_hog_inter <- development_area_list %>%
    model_development(., n_switch, n, h, m, pma) %>%
    model_reduction_vulnerability(., ua_1) %>% 
    model_reduction_hazard(., ua_2) %>%
    model_risk(., ua_1, f, ap, r)
  
  flooded_houses <- unlist(map_dbl(list_hog_inter, ~ sum(.x$inun == TRUE)))
  hog_flooded <- (flooded_houses / unlist(map_dbl(list_hog_inter, ~ nrow(.x)))) * 100
  
  flooded_once <- unlist(map_dbl(list_hog_inter, ~ sum(.x$inun_once == TRUE)))
  hog_flooded_once <- (flooded_once / unlist(map_dbl(list_hog_inter, ~ nrow(.x)))) * 100
  
  measures <- unlist(map_dbl(list_hog_inter, ~ sum(.x$lev == TRUE)))
  hog_measures <- (measures / unlist(map_dbl(list_hog_inter, ~ nrow(.x)))) * 100
  
  return_list <- list(list_hog_inter, hog_flooded, hog_flooded_once, hog_measures)
  names(return_list) <- c('list_hog_inter', 'hog_flooded', 'hog_flooded_once', 'hog_measures')
  return(return_list)
}

# test individuales

# test1 for calibration and validation-----

# run one model integrated
n <- n; n_switch <- c(15,25); h <- hog_initial; m <- hog_growth; ua_1 <- 0.7; ua_2 <- 100 ; f <- 0.5 ; ap <- 0.1; r <- 7 
pma <- c(0.24, 0.76, 0, 0)
t_beg <- Sys.time()
test1 <- model_integrated(development_area_list, n_switch, n, h, m, pma, ua_1, ua_2, f, ap, r)
t_end <- Sys.time()
test_time <- t_end - t_beg; test_time

#survey 2017
obs.flood2017 <- sum_2017$porc_inun
obs.measu2017 <- sum_2017$porc_med

#survey 2024
obs.flood2024 <- sum_2024$porc_inun
obs.measu2024 <- sum_2024$porc_med

{
#plot
plot(seq(year_initial, year_end), test1$hog_flooded_once, type = "l", xlab = "Años", ylab = "% hogares afectados", ylim = c(0, 100), col = "red")
lines(seq(year_initial, year_end), test1$hog_measures, type = "l", xlab = "Años", ylab = "% hogares que toman medidas", ylim = c(0, 100), col = "black")
lines(seq(year_initial, year_end), test1$hog_flooded, type = "l", xlab = "Años", ylab = "% hogares que toman medidas", ylim = c(0, 100), col = "blue")

points(2017,obs.flood2017, pch=3, col="red")
points(2017,obs.measu2017)
points(2024,obs.flood2024, pch=3, col="blue")
points(2024,obs.measu2024)
}

# test 2. different flood behaviour models----

# parameters
# n <- n; n_switch <- c(15,25); h <- hog_initial; m <- hog_growth; ua_1 <- 0.5; ua_2 <- 100 ; f <- 0.5 ; ap <- 0.1 ; r <- 7 
# 
# status_quo <- model_integrated(development_area_list, n_switch, n, h, m, c(1,0,0,0), ua_1, ua_2, f, ap, r)
# reactivo   <- model_integrated(development_area_list, n_switch, n, h, m, c(0,1,0,0), ua_1, ua_2, f, ap, r)
# olvido     <- model_integrated(development_area_list, n_switch, n, h, m, c(0,0,1,0), ua_1, ua_2, f, ap, r)
# proactivo  <- model_integrated(development_area_list, n_switch, n, h, m, c(0,0,0,1), ua_1, ua_2, f, ap, r)
# 
# # plot of behaviour models
# 
# #flooded by year
# {
# plot(seq(year_initial, year_end), status_quo$hog_flooded, type = "l", 
#      xlab = "Años", ylab = "% hogares afectados", ylim = c(0, 100), col = "red")
# lines(seq(year_initial, year_end), reactivo$hog_flooded, type = "l", 
#       xlab = "Años", ylab = "", ylim = c(0, 100), col = "black")
# lines(seq(year_initial, year_end), olvido$hog_flooded, type = "l", 
#       xlab = "Años", ylab = "", ylim = c(0, 100), col = "green")
# lines(seq(year_initial, year_end), proactivo$hog_flooded, type = "l", 
#       xlab = "Años", ylab = "", ylim = c(0, 100), col = "blue")
# }
# 
# # flooded once
# {
# plot(seq(year_initial, year_end), status_quo$hog_flooded_once, type = "l", 
#      xlab = "Años", ylab = "% hogares afectados", ylim = c(0, 100), col = "red")
# lines(seq(year_initial, year_end), reactivo$hog_flooded_once, type = "l", 
#       xlab = "Años", ylab = "", ylim = c(0, 100), col = "black")
# lines(seq(year_initial, year_end), olvido$hog_flooded_once, type = "l", 
#       xlab = "Años", ylab = "", ylim = c(0, 100), col = "green")
# lines(seq(year_initial, year_end), proactivo$hog_flooded_once, type = "l", 
#       xlab = "Años", ylab = "", ylim = c(0, 100), col = "blue")
# }

# calibration and validation----

# a. global statistics----

# calibration
(obs.flood2017 <- sum_2017$porc_inun)
(obs.measu2017 <- sum_2017$porc_med)

(sim.flood2017 <- test1$hog_flooded_once[[33]]) 
(sim.measu2017 <- test1$hog_measures[[33]])

# validation
(obs.flood2024 <- sum_2024$porc_inun)
(obs.measu2024 <- sum_2024$porc_med)

(sim.flood2024 <- test1$hog_flooded[[40]]) 
(sim.measu2024 <- test1$hog_measures[[40]])

# Summary Table

summary_table1 <- matrix(
  c(obs.flood2017, sim.flood2017,
    obs.flood2024, sim.flood2024,
    obs.measu2017, sim.measu2017,
    obs.measu2024, sim.measu2024),
  nrow = 4, byrow = TRUE,
  dimnames = list(
    c("Calibration Flood 2017", "Validation Flood 2024", "Calibration Measures 2017", "Validation Measures 2024"),
    c("Observed", "Simulated")
  )
)

# Convert to data frame
(df_summary1 <- as.data.frame(summary_table1)%>% round(2))

# b. Nearest Neighbour Index (NNI)----

# b.1. NNI for calibration (2017)

# b.1.2. Envelope polygon

db=100 # buffer for envelope area (optional)

# observed calibration
obs_envelope2017 <- survey_2017 %>%
  group_by(nucleo) %>%
  summarise(geometry = st_union(geometry)) %>%  # combine points by group
  mutate(geometry = st_convex_hull(geometry))%>% # envelope polygon
  st_buffer(dist = db) %>%
  mutate(area = st_area(geometry)) # area

(obs.area2017 <- sum(obs_envelope2017$area)) #area [m2]

# simulated calibration

#clip simulated points with polygon of observed points
sim_point2017_sf <- st_as_sf(test1$list_hog_inter[[33]]) %>% #33 corresponde a 2017
sf::st_intersection(obs_envelope2017) 

# clustering with DBSCAN (equivalent to nucleo field in survey points)
coords2017 <- st_coordinates(sim_point2017_sf)
dbscan_result2017 <- dbscan(coords2017, eps = 1000, minPts = 5)
sim_point2017_sf$cluster <- dbscan_result2017$cluster

sim_envelope2017 <- sim_point2017_sf %>%
  group_by(cluster) %>%
  summarise(geometry = st_union(geometry)) %>%  # combine points by group
  mutate(geometry = st_convex_hull(geometry))%>% # envelope polygon
  st_buffer(dist = db) %>%
  mutate(area = st_area(geometry)) # area

(sim.area2017 <- sum(sim_envelope2017$area)) #area [m2]

# plots of calibration
#plot(st_geometry(sim_point2017_sf), col = sim_point2017_sf$cluster, pch = 16, main = "Agrupamiento de Puntos por Grupo")
plot(st_geometry(obs_envelope2017), main="calibracion")
points(survey_2017,cex = 0.5)
points(sim_point2017_sf, cex = 0.5,col="red")

# b.1.1. Mean distance to nearest neighbour

# observed calibration
obs_point2017 <- as.ppp(survey_2017)
obs_point2017$dist <- nndist(obs_point2017)
(obs_mean_dist2017 <- mean(obs_point2017$dist))

# simulated calibration
#sim_point2017_sf <- st_as_sf(test1$list_hog_inter[[33]]) #33 corresponde a 2017
sim_point2017 <- as.ppp(sim_point2017_sf)
sim_point2017$dist <- nndist(sim_point2017)
(sim_mean_dist2017 <- mean(sim_point2017$dist))

# b.1.3. Nearest Neighbour Index

(n.obs2017 <- length(survey_2017$geometry))
(n.sim2017 <- length(sim_point2017_sf$geometry))

(rsp.obs2017 <- 0.5*sqrt(obs.area2017/n.obs2017))
(rsp.sim2017 <- 0.5*sqrt(sim.area2017/n.sim2017))

(NNI.obs2017 <-  obs_mean_dist2017/rsp.obs2017)
(NNI.sim2017 <-  sim_mean_dist2017/rsp.sim2017)

# b.2. NNI for validation (2024)

# b.2.2. Envelope polygon

# observed validation
obs_envelope2024 <- survey_2024 %>%
  filter(folio != 2035) %>% # wrong location
  group_by(distrito) %>%
  summarise(geometry = st_union(geometry)) %>%  # combine points by group
  mutate(geometry = st_convex_hull(geometry))%>% # envelope polygon
  st_buffer(dist = db) %>%
  mutate(area = st_area(geometry)) # area

(obs.area2024 <- sum(obs_envelope2024$area)) #area [m2]

# simulated validation

#clip simulated points with polygon of observed points
sim_point2024_sf <- st_as_sf(test1$list_hog_inter[[40]]) %>% #2024
sf::st_intersection(obs_envelope2024) 

# clustering with DBSCAN (equivalent to nucleo/distrito field in survey points)

coords2024 <- st_coordinates(sim_point2024_sf)
dbscan_result2024 <- dbscan(coords2024, eps = 1000, minPts = 5)
sim_point2024_sf$cluster <- dbscan_result2024$cluster

sim_envelope2024 <- sim_point2024_sf %>%
  group_by(cluster) %>%
  summarise(geometry = st_union(geometry)) %>%  # combine points by group
  mutate(geometry = st_convex_hull(geometry))%>% # envelope polygon
  st_buffer(dist = db) %>%
  mutate(area = st_area(geometry)) # area

(sim.area2024 <- sum(sim_envelope2024$area)) #area [m2]

# plots of validation
#plot(st_geometry(sim_point2024_sf), col = sim_point2024_sf$cluster, pch = 16, main = "Agrupamiento de Puntos por Grupo")
plot(st_geometry(obs_envelope2024), main="validation")
points(survey_2024,cex = 0.5)
points(sim_point2024_sf, cex = 0.5,col="red")

# b.2.1. Mean distance to nearest neighbour

# observed validation
obs_point2024 <- as.ppp(survey_2024)
obs_point2024$dist <- nndist(obs_point2024)
(obs_mean_dist2024 <- mean(obs_point2024$dist))

# simulated validation
#sim_point2024_sf <- st_as_sf(test1$list_hog_inter[[40]]) #40 corresponde a 2024
sim_point2024 <- as.ppp(sim_point2024_sf)
sim_point2024$dist <- nndist(sim_point2024)
(sim_mean_dist2024 <- mean(sim_point2024$dist))

# b.2.3. Nearest Neighbour Index

(n.obs2024 <- length(survey_2024$geometry))
(n.sim2024 <- length(sim_point2024_sf$geometry))

(rsp.obs2024 <- 0.5*sqrt(obs.area2024/n.obs2024))
(rsp.sim2024 <- 0.5*sqrt(sim.area2024/n.sim2024))

(NNI.obs2024 <-  obs_mean_dist2024/rsp.obs2024)
(NNI.sim2024 <-  sim_mean_dist2024/rsp.sim2024)

# c. Hopkins-Skellam test statistic----

# calibration
#(hopskel.obs2017 <- hopskel.test(obs_point2017, alternative ="clustered"))
(hopskel.obs2017 <- hopskel(obs_point2017))
(hopskel.sim2017 <- hopskel(sim_point2017))

# validation
(hopskel.obs2024 <- hopskel(obs_point2024))
(hopskel.sim2024 <- hopskel(sim_point2024))

# Summary table of points spatial patterns indexes

summary_table2 <- matrix(
  c(as.numeric(NNI.obs2017), as.numeric(NNI.sim2017),
    as.numeric(NNI.obs2024), as.numeric(NNI.sim2024),
    hopskel.obs2017, hopskel.sim2017,
    hopskel.obs2024, hopskel.sim2024),
  nrow = 4, byrow = TRUE,
  dimnames = list(
    c("Calibration NNI 2017", "Validation NNI 2024", "Calibration hopskel 2017","Validation hopskel 2024"),
    c("Observed", "Simulated")
  )
)

# Convert to data frame
(df_summary2 <- as.data.frame(summary_table2) %>% round(7))

# Interpretation of NNI
# NNI ≈ 1: Random
# NNI < 1: Clustered
# NNI > 1: Uniform

# Interpretation of Hopskel test (H)
# H ≈ 0: Uniform
# H ≈ 0.5: Random
# H ≈ 1: Clustered

# d. kappa----
{
nb <- 3 # number of neasrest neighbor

# d.1. kappa for calibration (2017)

# nearest-neighbour
nearest.2017 <- st_nn(sim_point2017_sf, survey_2017, k = nb, returnDist = FALSE)

# list with attributes of NN for each point
attributes.2017 <- lapply(nearest.2017, function(neighbors) {
  survey_2017[neighbors, c("bienes", "danos", "medidas_proteccion")]
})

# most repeated value (mode) for each column in each element of the list
modes.list.2017 <- lapply(attributes.2017, function(subset_data) {
   col_modes <- sapply(subset_data, function(column) {
   most_frequent <- names(sort(table(column), decreasing = TRUE))[1]
   return(most_frequent)
  })
  return(col_modes)
})

# convert list to dataframe
modes.df.2017 <- do.call(rbind, lapply(modes.list.2017, function(x) {
  data.frame(bienes = x[1], danos = x[2], medidas_proteccion = x[3])
}))


# Table for metrics
sim.points.2017.mode <- cbind(sim_point2017_sf, modes.df.2017)

# edition to obtain metrics
sim.points.2017.mode <- sim.points.2017.mode %>%
  select(lev,inun_once,medidas_proteccion,bienes,danos) %>% #filter
  mutate(afectado = ifelse(bienes == "vivienda" & (danos == "moderados" | danos == "graves"), "si", "no")) %>%
  rename(sim_measu = lev,                        #rename
         sim_afect = inun_once,
         obs_measu = medidas_proteccion,
         obs_afect = afectado) %>%
  mutate(                                        #change values 
    sim_measu = if_else(sim_measu == "TRUE", 1, 0),
    sim_afect = if_else(sim_afect == "TRUE", 1, 0),
    obs_measu = if_else(obs_measu == "si", 1, 0),
    obs_afect = if_else(obs_afect == "si", 1, 0),
 )

# confusion matrix and kappa with caret

# measure
(conf.matrix.2017.measu <- confusionMatrix(factor(sim.points.2017.mode$sim_measu),
                                           factor(sim.points.2017.mode$obs_measu)))

kappa.measu.2017 <- conf.matrix.2017.measu$overall[["Kappa"]]

# afected
(conf.matrix.2017.afect <- confusionMatrix(factor(sim.points.2017.mode$sim_afect),
                               factor(sim.points.2017.mode$obs_afect)))

kappa.afect.2017 <- conf.matrix.2017.afect$overall[["Kappa"]]

# confusion matrix and kappa with functions

f.kappa <- function(sim,obs){
  # Calcular la matriz de confusión
  TP <- sum(sim == 1 & obs == 1)
  TN <- sum(sim == 0 & obs == 0)
  FP <- sum(sim == 1 & obs == 0)
  FN <- sum(sim == 0 & obs == 1)
  
  #Accuracy
  accuracy.o <- (TP + TN) / (TP + TN + FP + FN) # precision observada
  accuracy.e <- ((TP + FP) / (TP + TN + FP + FN)) * ((TP + FN) / (TP + TN + FP + FN)) +
    ((TN + FN) / (TP + TN + FP + FN)) * ((TN + FP) / (TP + TN + FP + FN)) # precision esperada
  
  kappa <- (accuracy.o - accuracy.e) / (1- accuracy.e) 
  
  sensitivity <- TP / (TP + FN) #sensividad
  specificity <- TN / (TN + FP) #especificidad
  
  ppv <- TP / (TP+FP) # positive predictive value
  npv <- TN / (TN+FN) # negative predictive value
  
  return(c(kappa, accuracy.o, sensitivity, specificity, ppv,npv))
  
}

(f.kappa.measu.2017 <- f.kappa(sim.points.2017.mode$sim_measu, sim.points.2017.mode$obs_measu))

(f.kappa.afect.2017 <- f.kappa(sim.points.2017.mode$sim_afect, sim.points.2017.mode$obs_afect))

# d.2. kappa for validation (2024)

# nearest-neighbour
nearest.2024 <- st_nn(sim_point2024_sf, survey_2024, k = nb, returnDist = FALSE) 

# list with attributes of NN for each point
attributes.2024 <- lapply(nearest.2024, function(neighbors) {
  survey_2024[neighbors, ] %>%
    dplyr::select(d19, d20, contains("g36_1"), contains("g36_2")) %>%
    #dplyr::select(d19, d20, contains("g36_1")) %>%
    mutate(measure = coalesce(g36_1_1, g36_1_2, g36_1_3, g36_1_6, g36_1_7, g36_1_11, g36_1_12, 
                              g36_2_1, g36_2_2, g36_2_3, g36_2_6, g36_2_7, g36_2_11, g36_2_12)) # NA in neighbors without measures
    #mutate(measure = coalesce(g36_1_1, g36_1_2, g36_1_3,g36_1_6, g36_1_7,g36_1_11, g36_1_12))
  })

# most repeated value (mode) for each column in each element of the list. 
modes.list.2024 <- lapply(attributes.2024, function(subset_data) {
  col_modes <- sapply(subset_data, function(column) {
    most_frequent <- names(sort(table(column), decreasing = TRUE))[1]
    return(most_frequent)
  })
  return(col_modes)
})

# convert list to dataframe
# modes.df.2024 <- do.call(rbind, lapply(modes.list.2024, function(x) {
#   data.frame(x[c("d19", "d20", "measure"), drop = FALSE]) 
# }))

modes.df.2024 <- do.call(rbind, lapply(modes.list.2024, function(x) {
  # Convert each element to a data frame and transpose it
  x <- as.data.frame(t(x))  
  
  # Define the required columns
  vars <- c("d19", "d20", "measure")
  
  # Check which columns are missing and add them with NA
  missing_vars <- setdiff(vars, names(x))  # Find missing columns
  x[missing_vars] <- NA  # Add missing columns with NA values
  
  # Select only the required columns
  x <- x[, vars, drop = FALSE]
  
  return(x)
}))

#any(is.na(modes.df.2024))

# Table for metrics
sim.points.2024.mode <- cbind(sim_point2024_sf, modes.df.2024)

# edition to obtain metrics
sim.points.2024.mode <- sim.points.2024.mode %>%
  select(lev, inun, measure, d19, d20) %>%  # 
  rename(
    sim_measu = lev,   # rename
    sim_afect = inun   # related to last flood
  ) %>%
  mutate(  
    obs_measu = if_else(measure=="NULL", 0, 1),
    #obs_measu = if_else(!is.na(measure), 1, 0),  
    #obs_afect = if_else(d19 == "Sí" & d20 == "Sí", 1, 0),  
    obs_afect = if_else(d19 == "Sí", 1, 0),  
    sim_measu = if_else(sim_measu == TRUE, 1, 0),  
    sim_afect = if_else(sim_afect == TRUE, 1, 0)
    #obs_measu = if_else(obs_measu == "si", 1, 0),  
    #obs_afect = if_else(obs_afect == "si", 1, 0)  
  )

# confusion matrix and kappa with caret

# measure
(conf.matrix.2024.measu <- confusionMatrix(factor(sim.points.2024.mode$obs_measu),
                                           factor(sim.points.2024.mode$sim_measu)))

kappa.measu.2024 <- conf.matrix.2024.measu$overall[["Kappa"]]

# afected
(conf.matrix.2024.afect <- confusionMatrix(factor(sim.points.2024.mode$obs_afect),
                                           factor(sim.points.2024.mode$sim_afect)))

kappa.afect.2024 <- conf.matrix.2024.afect$overall[["Kappa"]]
#conf.matrix.2017.afect$overall

(f.kappa.measu.2024 <- f.kappa(sim.points.2024.mode$sim_measu, sim.points.2024.mode$obs_measu))

(f.kappa.afect.2024 <- f.kappa(sim.points.2024.mode$sim_afect, sim.points.2024.mode$obs_afect))

# e. chi squared test----

# calibration

# measures
table.chisq.measu.2017 <- table(sim.points.2017.mode$obs_measu,sim.points.2017.mode$sim_measu)
(test.chisq.measu.2017 <- chisq.test(table.chisq.measu.2017))

# afected
table.chisq.afect.2017 <- table(sim.points.2017.mode$obs_afect,sim.points.2017.mode$sim_afect)
(test.chisq.afect.2017 <- chisq.test(table.chisq.afect.2017))

# validation

# measures
(table.chisq.measu.2024 <- table(sim.points.2024.mode$obs_measu,sim.points.2024.mode$sim_measu))
(test.chisq.measu.2024 <- chisq.test(table.chisq.measu.2024))


# afected
(table.chisq.afect.2024 <- table(sim.points.2024.mode$obs_afect,sim.points.2024.mode$sim_afect))
(test.chisq.afect.2024 <- chisq.test(table.chisq.afect.2024))

# plots bars

plot_chisq_bar <- function(chisq_table, main) {
  # Extract observed and simulated counts
  observed_0 <- sum(chisq_table[1, ])  # Sum of first row (0 observed)
  observed_1 <- sum(chisq_table[2, ])  # Sum of second row (1 observed)
  simulated_0 <- sum(chisq_table[, 1]) # Sum of first column (0 simulated)
  simulated_1 <- sum(chisq_table[, 2]) # Sum of second column (1 simulated)
  
  # Create a data frame for plotting
  bar_data <- data.frame(
    Category = c("Observed", "Observed", "Simulated", "Simulated"),
    Label = c("0", "1", "0", "1"),
    Count = c(observed_0, observed_1, simulated_0, simulated_1)
  )
  
  # Plot using ggplot2
  library(ggplot2)
  ggplot(bar_data, aes(x = Label, y = Count, fill = Category)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = main,
         x = "Category (0 = No, 1 = Yes)",
         y = "Count") +
    theme_minimal() +
    scale_fill_manual(values = c("Observed" = "blue", "Simulated" = "red"))
}

# Example usage
plot_chisq_bar(table.chisq.measu.2017, "cal-measure")
plot_chisq_bar(table.chisq.afect.2017, "cal-afected")
plot_chisq_bar(table.chisq.measu.2024, "val-measure")
plot_chisq_bar(table.chisq.afect.2024, "val-afected")

# summary table 

summary_table3 <- matrix(
  c(kappa.measu.2017, kappa.afect.2017,
    kappa.measu.2024, kappa.afect.2024,
    test.chisq.measu.2017$p.value,test.chisq.afect.2017$p.value,
    test.chisq.measu.2024$p.value,test.chisq.afect.2024$p.value),
  nrow = 4, byrow = TRUE,
  dimnames = list(
    c("Calibration.kappa.2017","Validation.kappa.2024","Calibration.Chi-Sq.2017", "Validation.Chi-Sq.2024"),
    c("Measures", "Afected")
  )
)

# Convert to data frame
(df_summary3 <- as.data.frame(summary_table3) %>% round(5))

table.chisq.measu.2024
table.chisq.afect.2024
}

# Interpretation of the Kappa value (k):
# k=1: Perfect agreement between predictions and observations.
# k=0: Agreement similar to that expected by chance.
# k<0: There is no agreement and the predictions are worse than chance.
# Values between 0 and 1 indicate the quality of the agreement, where values close to 1 indicate high precision and values close to 0 indicate low precision.

# Chi-squared
#H0= there are no differences between observed and simulated frequencies
#p-value = 0.12, H0 is accepted, there are no differences

# list of changes----

# L162: crop de los puntos de encuesta 2017 con el area de inundacion
# L189 se incorpora preparados futuros en la encuesta 2024 (g36_2)
# L358 aprendizaje cuando hay inundacion en el modelo olvido (+ ap)
# L414 cambio en los periodos para el uso del suelo urbano n_switch <- c(15,25)


# code trash-----

# plot of behaviour models

# data_combined <- bind_rows(
#   data.frame(year = seq(year_initial, year_end), hog_flooded = status_quo$hog_flooded, hog_flooded_once = status_quo$hog_flooded_once, comportamiento = "Status Quo"),
#   data.frame(year = seq(year_initial, year_end), hog_flooded = reactivo$hog_flooded, hog_flooded_once = reactivo$hog_flooded_once, comportamiento = "Reactivo"),
#   data.frame(year = seq(year_initial, year_end), hog_flooded = olvido$hog_flooded, hog_flooded_once = olvido$hog_flooded_once, comportamiento = "Olvido"),
#   data.frame(year = seq(year_initial, year_end), hog_flooded = proactivo$hog_flooded, hog_flooded_once = proactivo$hog_flooded_once, comportamiento = "Proactivo")
# )
# 
# data_labels <- data_combined %>%
#   group_by(comportamiento) %>%
#   filter(year == max(year))
# 
# p1 <- ggplot(data_combined, aes(x = year, y = hog_flooded, color = comportamiento)) +
#   geom_line() +
#   geom_text_repel(data = data_labels, aes(label = comportamiento), 
#                   nudge_x = 0.5, direction = "y", show.legend = FALSE) +
#   labs(
#     y = "% Hogares afectados por año"
#   ) +
#   scale_x_continuous(
#     breaks = seq(year_initial, year_end, by = 5),      
#     minor_breaks = seq(year_initial, year_end, by = 1) 
#   ) +
#   theme_classic() +
#   theme(
#     legend.position = "none",              
#     axis.ticks.length = unit(5, "pt"),      
#     axis.text.x = element_text(hjust = 0.5) 
#   )
# 
# p2 <- ggplot(data_combined, aes(x = year, y = hog_flooded_once, color = comportamiento)) +
#   geom_line() +
#   geom_text_repel(data = data_labels, aes(label = comportamiento), 
#                   nudge_x = 0.5, direction = "y", show.legend = FALSE) +
#   labs(
#     y = "% Hogares afectados acumulado"
#   ) +
#   scale_x_continuous(
#     breaks = seq(year_initial, year_end, by = 5),      
#     minor_breaks = seq(year_initial, year_end, by = 1) 
#   ) +
#   theme_classic() +
#   theme(
#     legend.position = "none",               
#     axis.ticks.length = unit(5, "pt"),      
#     axis.text.x = element_text(hjust = 0.5) 
#   )
# 
# grid.arrange(p1, p2, ncol = 1)
# ..................................................................................................
# Project: Climate Change and the sociohydrology of floods
# Authors: Jorge Hurtado & Vicente Saenger
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
#        6) Serie de QMA (Tr) obtenidos de caudales diarios (1979-1983) de DGA (cr2-camels). La serie_return_period que usa es temporal hasta definir la final.
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

# glossary ----

# n: number of years
# n_switch: year of switch area of development 
# h: number of initial households 
# m: number of new households per year
# pfb: proportion of flood behaviors at household level
# household_compliance_threshold: threshold at household level at which flood-proof-measures are taken 
# ua_2: threshold at state level in function of flood return period 
# f: threshold at household level that damage of flood depth can be reduced by adaptive measures
# learning: learning rate 
# r: years in which what has been learned is forgotten  

# assumptions ----


# libraries -----

{
  library(tidyverse)
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
  library(yardstick)
  library(leafpop)
  library(tools)
}
###
#setwd("G:\\My Drive\\010_code\\010_r\\project_sociohydro_abm_flood_risk_r_code\\abm_flood_risk")

# ..................................................................................................
# set up ----
year_initial <- 1985
year_end <- 2024
n <- year_end - year_initial + 1
crs = "EPSG:32718"

# domain hydraulic model 
domain_hydraulic_model <- st_read("domain_hydraulic_model/domain_hydraulic_model.shp", crs = crs)

# series discharge
df_discharge <- read.csv("series_discharge\\serie_4_1.csv", header = TRUE)
series_return_period <- df_discharge$qtr
year_simulation <- df_discharge$year
barplot(series_return_period, main = "Return period", names.arg = year_simulation, las = 2, cex.names = 0.8)

# series flood raster 
load_flood_vect_list <- function(series_return_period, folder_path) {
  
  #folder_path <- "series_flood_shp\\shp"
  vect_files <- list.files(path = folder_path, pattern = "\\.gpkg$", full.names = TRUE) # shp files
  
  floodlist <- vect_files %>% 
    map(st_read)
  floodname <- vect_files %>% 
    map_chr(~ file_path_sans_ext(basename(.x)))
  names(floodlist) <- floodname
  
  # define the names of the flood layers
  flood_layers <- c("1" = "FloodQT1", "2" = "FloodQT2", "5" = "FloodQT5", 
                    "10" = "FloodQT10", "25" = "FloodQT25", "35" = "FloodQT35", 
                    "50" = "FloodQT50", "100" = "FloodQT100", "200" = "FloodQT200")
  
  # associate each element of the series_return_period with the corresponding raster
  floodseries <- map(series_return_period, function(x) {
    layer_name <- flood_layers[as.character(x)]
    if (!is.null(layer_name) && layer_name %in% names(floodlist)) {
      return(floodlist[[layer_name]])
    } else {
      warning(paste("no raster for year:", x))
      return(NULL)
    }
  })
  return(floodseries)
}
list_floodseries <- load_flood_vect_list(series_return_period, "series_flood_vect\\vect") # load raster without project
#list_floodseries_wp <- load_flood_raster_list(series_return_period, "series_flood_raster\\raster_depth_10m_with_project") # load raster with
year_simulation <- seq(year_initial, year_end, by = 1)
names(list_floodseries) <- year_simulation
#names(list_floodseries_wp) <- year_simulation
list_floodseries <- list_floodseries %>%
  map(~ select(.x, "max_depth"))


# development area
list_development_area <- list.files("development_area", pattern = "_urban.shp$", full.names = TRUE) %>%
  map(., ~ st_read(.x) %>% st_union() %>% st_transform(crs))
list_development_area <- map(list_development_area, ~ st_intersection(.x, list_development_area[[3]]))
mapview(list_development_area[[1]])+mapview(list_development_area[[2]])+mapview(list_development_area[[3]])

# survey
df_villa <- st_sf(
  name_villa = c("villa esperanza", "villa molino el sol"),  
  geometry = st_sfc(
    st_polygon(list(matrix(c(
      655627.2, 655830.1, 655602.4, 655571.0, 655627.2,  
      5875653.8, 5875835.1, 5875899.1, 5875665.2, 5875653.8  
    ), ncol = 2, byrow = FALSE))),
    
    st_polygon(list(matrix(c(
      654685, 654767.4, 654645.8, 654613.6, 654685,  
      5869646.7, 5869884.1, 5869909.6, 5869721.6, 5869646.7  
    ), ncol = 2, byrow = FALSE)))
  ),
  crs = crs  # definir el sistema de referencia
)

# survey 
df_survey <- bind_rows(
  # survey 2017
  read.csv("survey/survey_2017.csv", skip = 1) %>%
    select(where(~ !all(is.na(.)))) %>%
    rename_with(tolower) %>%
    filter(!is.na(objectid)) %>%
    mutate_if(is.character, tolower) %>%
    sf::st_as_sf(., coords = c("coord_x", "coord_y"), crs = crs) %>%
    dplyr::mutate(household_id_survey = objectid,
                  household_address = dirección,
                  household_flood_once = afectados == "si",
                  household_measure = medidas_proteccion == "si",
                  household_flood_risk_perception = recode(nivel_riesgo,
                                                           "alto" = 4,
                                                           "medio" = 3,
                                                           "bajo" = 2,
                                                           "nulo" = 1)) %>%
    mutate(year_current = 2017),
  # survey 2024
  read.csv("survey/survey_2024.csv",  na.strings = c("", "-")) %>%
    rename_with(tolower) %>%
    filter(comuna == "Arauco") %>%
    mutate_if(is.character, tolower) %>%
    sf::st_as_sf(., coords = c("lon_encuesta", "lat_encuesta"), crs = 4326) %>%
    sf::st_transform(crs) %>%
    rowwise() %>%
    dplyr::mutate(household_id_survey = folio,
                  household_address = paste(direccion,numero,sep=" "),
                  #household_flood_year_current = ifelse(d19 == "sí" & d20 == "sí", TRUE, FALSE),
                  household_flood_year_current = case_when(
                    d19 == "sí" ~ TRUE,
                    d19 == "no" ~ FALSE,
                    d19 == "no aplica" ~ FALSE,
                    TRUE ~ NA),
                  household_measure = !is.na(coalesce(g36_1_1, g36_1_2, g36_1_3, g36_1_4, g36_1_5, g36_1_6, g36_1_7, 
                                                      g36_1_8, g36_1_9, g36_1_10, g36_1_11, g36_1_12, g36_1_13, g36_1_14, 
                                                      g36_1_15, g36_1_16, g36_1_17, g36_1_18, g36_1_19, g36_1_20, g36_1_21,
                                                      g36_1_22, g36_1_23, g36_1_24, g36_1_25 #,
                                                      #g36_2_1, g36_2_2, g36_2_3, g36_2_4, g36_2_5, g36_2_6, g36_2_7, 
                                                      #g36_2_8, g36_2_9, g36_2_10, g36_2_11, g36_2_12, g36_2_13, g36_2_14, 
                                                      #g36_2_15, g36_2_16, g36_2_17, g36_2_18, g36_2_19, g36_2_20, g36_2_21,
                                                      #g36_2_22, g36_2_23, g36_2_24, g36_2_25
                                                      )),
                  household_measure_resp = if_else(
                    household_measure,
                    str_trunc(paste(na.omit(c(
                      g36_1_1, g36_1_2, g36_1_3, g36_1_4, g36_1_5, g36_1_6, g36_1_7, g36_1_8, g36_1_9, 
                      g36_1_10, g36_1_11, g36_1_12, g36_1_13, g36_1_14, g36_1_15, g36_1_16, g36_1_17, 
                      g36_1_18, g36_1_19, g36_1_20, g36_1_21, g36_1_22, g36_1_23, g36_1_24, g36_1_25
                    )), collapse = ", "), 100),
                    NA_character_
                  ),
                  household_measure_flood_proof = !is.na(coalesce(g36_1_1, g36_1_2, g36_1_3, g36_1_7, g36_1_12 #, 
                                                                  #g36_2_1, g36_2_2, g36_2_3
                                                                  )),
                  household_measure_flood_proof_resp = ifelse( household_measure_flood_proof, 
                                                               str_trunc(paste(unique(na.omit(c(g36_1_1, g36_1_2, g36_1_3, g36_1_7, g36_1_12
                                                               ))), collapse = ", "), 100), NA),
                  household_flood_risk_perception = recode(e26,
                                                           "muy alto" = 5,
                                                           "alto" = 4,
                                                           "ni bajo ni alto" = 3,
                                                           "bajo" = 2,
                                                           "muy bajo" = 1),
                  household_year_existence = c9_1) %>%
    mutate(year_current = 2024,
           household_year_entry = year_current - household_year_existence)
) %>%
  mutate(
    household_address = str_replace_all(household_address, c(
      "," = "",
      "nº" = " ",
      "arturo perez canto|arturo perez" = "arturo perez canto",
      "avda" = "avenida",
      "carpites" = "carpinteros",
      "calle malloga" = "callejon malloga",
      "calle 1 df_villa esperanza" = "df_villa esperanza calle 1",
      "campolican|compolican" = "caupolican",
      "casa" = "",
      "freisa|freise|calle fresia" = "fresia",
      "fresia con o'higgins" = "fresia con o'higgins sn",
      "ignaio" = "ignacio",
      "hernan pelen pucheu|hernan pelen puchen|hernan pelen" = "hernan pelen pucheu",
      "millar|milan" = "millan",
      "montts" = "montt",
      "nolinos|molinos" = "molino",
      "villa molino del sol|villa el molino el sol|villa molino del sol|molino del sol" = "villa molino el sol", 
      "naitenez|naitez|los naitenes" = "maitenes",
      "pasaje 10 de julio|pasaje 10 julio|10 julio" = "10 de julio",
      "peaje|pasje" = "pasaje",
      "pesaje 61" = "pasaje 2 61",
      "punto" = "pinto",
      "(sin numero)|sin número|(sn)" = "sn",
      "extension manuel zarate|manuel zarrate|manuel zarate" = "samuel zarate",
      "ohiggins" = "o'higgins",
      "viecente" = "vicente",
      "villa esperanza pasaje 1" = "villa esperanza calle 1",
      "virginia fondo" = "virginia pardo",
      "tres" = "3",
      "  |   " = " "
    )), household_address = household_address) %>%
  select(household_id_survey, household_address, year_current, household_year_entry, household_year_existence, 
         household_measure, household_measure_resp, household_measure_flood_proof,household_measure_flood_proof_resp, household_flood_year_current, household_flood_once, 
         household_flood_risk_perception) %>%
  mutate(household_address = case_when(
    str_detect(household_address, "villa\\s+[\\w\\s]+") ~ 
      paste0(str_extract(household_address, "villa\\s+[\\w\\s]+"), " ", 
             str_remove(household_address, "villa\\s+[\\w\\s]+,?\\s*")),
    TRUE ~ household_address
  )) %>%
  mutate(household_address = case_when(
    # if it is a df_villa, leave the address as it is
    str_detect(household_address, paste(df_villa$name_villa, collapse = "|")) ~ household_address,
    # if the address is in a df_villa and dont have the df_villa name, add the df_villa name
    !is.na(sapply(st_within(geometry, df_villa), function(x) ifelse(length(x) > 0, x, NA))) ~ 
      paste0(df_villa$name_villa[sapply(st_within(geometry, df_villa), function(x) ifelse(length(x) > 0, x, NA))], 
             " ", household_address),
    # if the address is not in a df_villa do nothing
    TRUE ~ household_address
  )) %>%
  mutate(
    household_address = str_replace_all(household_address, c(
      "  " = " ", 
      "los boldos con samuel zarrate 16" = "samuel zarate 16",
      "o'higgins c fresia s n" = "fresia con o'higgins sn",
      "vicente millan pasaje los cuervos" = "pasaje los cuervos",
      "villa esperanza pasaje 1" = "villa esperanza calle 1"
    )), household_address = household_address) %>%
  mutate(household_address = str_trim(household_address))
mapview(df_survey)
# keep survey points that are within the domain of the hydraulic model
df_survey <- st_intersection(df_survey, domain_hydraulic_model)
  
# data for calibration: household interviewed in both surveys
df_household_calibration <- df_survey %>%
  group_by(household_address) %>%
  filter(n_distinct(year_current) > 1) %>%
  # for same address same geometry
  mutate(geometry = geometry[year_current == 2024][1]) %>% # 2017 or 2024 
  ungroup() #%>%
  #filter(grepl('villa', household_address))

# data for validation: household interviewed in only one survey
df_household_validation <- df_survey %>%
  group_by(household_address) %>%
  filter(n_distinct(year_current) == 1)

mapview(df_survey)


df_survey %>%
  mutate(decade = floor(household_year_entry / 10) * 10) %>%
  group_by(decade) %>%
  summarise(
    n = n(),
    household_with_flood_proof_measure = sum(household_measure_flood_proof, na.rm = TRUE),
    percentage_measures = (household_with_flood_proof_measure / n) * 100
  ) %>%
  arrange


df_survey %>%
  mutate(decade = floor(household_year_entry / 10) * 10) %>%
  filter(!is.na(decade)) %>%
  mutate(decade = factor(decade)) %>%
  mapview(zcol = "decade")

df_household <- df_household_calibration %>% 
  group_by(household_address) %>%
  filter(year_current == 2017 & household_measure == F | year_current == 2024 & household_measure == T)

# model: development ----
model_development <- function(list_development_area, n_switch, n, df_household, pfb) {
  #n <- n; n_switch <- c(10,25); h <- household_initial; m <- household_growth; b <- 10 ; household_compliance_threshold <- 0.7; ua_2 <- 100 ; f <- 0.5 ; learning <- 0.3 ; r <- 7 
  #pfb <- c(1, 0, 0, 0)
  #df_household <- df_household_validation
  
  set.seed(1234) 
  
  # range for development areas
  area_intervals <- c(n_switch - 1, n) 
  df_development_area_range <- data.frame(
    id_development_area = seq_along(list_development_area),
    start = c(1, area_intervals[-length(area_intervals)] + 1),
    end = area_intervals
  )
  
  # polygon of development areas
  polygon_development_area <- do.call(rbind, lapply(seq_along(list_development_area), function(i) {
    st_sf(id_development_area = i, geometry = st_geometry(list_development_area[[i]]))
  }))
  
  # assign polygons to households
  df_household <- st_join(df_household, polygon_development_area) %>% 
    group_by(geometry) %>% 
    slice_min(id_development_area, na_rm = F) %>%  
    ungroup() #%>% 
    #st_buffer(1)
  
  
  # assign time step  
  assign_time_step <- function(area_id) {
    area_info <- df_development_area_range %>% filter(id_development_area == area_id)
    
    if (nrow(area_info) == 0 || is.na(area_id)) {
      return(sample(seq(1, n), 1))  # random time step
    }
    
    sample(seq(area_info$start[1], area_info$end[1]), 1)   # random time step in the range of the first area of development in which the point is inside
  }
  
  # assign household year entry and create list of households
  list_household <- df_household %>%
    mutate(time_step = map_int(id_development_area, assign_time_step)) %>%
    mutate(
      # if household_year_entry = NA: assign random value 
      # else: preserve the household_year_entry value from the survey
      household_year_entry = if_else(is.na(household_year_entry), time_step + year_initial - 1, household_year_entry),
      # if household_year_entry < year_initial: assign year_initial
      household_year_entry = if_else(household_year_entry < year_initial, year_initial, household_year_entry)
    ) %>%
    complete(household_year_entry = 1985:2024) %>%
    split(.$household_year_entry)
  
  # set list as sf object
  list_household <- purrr::map(list_household, ~ sf::st_sf(.)) %>% 
    unname(.)
  
  # add flood behavior 
  list_household <- purrr::imap(list_household, ~ {
    .x %>%
      dplyr::mutate(
        household_id_model = seq_len(n()) + ifelse(.y == 1, 0, sum(sapply(list_household[1:(.y-1)], nrow))),
        household_fb = sample(1:4, n(), replace = TRUE, prob = pfb)
      )
  })
  
  # accumulate household
  list_household <- purrr::accumulate(list_household, base::rbind)
  
  # add time columns
  list_household <- purrr::imap(list_household, ~ {
    .x %>%
      dplyr::mutate(
        time_step = .y,
        year_current = year_initial + .y - 1,
        household_year_existence = year_current - household_year_entry
      ) %>%
      dplyr::select(household_id_survey, household_id_model, household_address, time_step, year_current, household_year_entry, household_year_existence, household_fb) %>%
      filter(!is.na(household_id_survey)) # delete artificial points in years with no new household
  })
  
  return(list_household)
}

# model: reduction of vulnerability ----
model_reduction_vulnerability <- function(list_household_model_development) { 
  #list_household_model_development <- list_household_inter;
  
  # add compliance rate
  household_compliance_rate <- list_household_model_development %>%
    map(as_tibble) %>%
    reduce(bind_rows) %>%
    distinct(household_id_model, household_fb) %>%
    mutate(household_compliance_rate = case_when(
      household_fb == 1 ~ 0,
      household_fb == 2 ~ 0,
      household_fb == 3 ~ 0,
      household_fb == 4 ~ 0,
      TRUE ~ 0 ))
  
  # bind
  list_household_model_development <- list_household_model_development %>%
    map(as_tibble) %>%
    reduce(bind_rows) %>%
    left_join(household_compliance_rate, by = c("household_id_model","household_fb")) %>%  
    base::split(.$year_current)
  
  return(list_household_model_development)
}

# model: reduction of hazard ----
model_reduction_hazard <- function(list_household_model_development, ua_2) { 
  #list_household_model_development <- list_household_inter; ua_2 <- 100
  
  alarm <- map_int(series_return_period, ~ ifelse(.x >= ua_2, 1, 0)) # build alarm series_return_period
  alarm_position <- base::which(alarm == 1)[1] # find alarm position
  
  # build series_return_period with 1 from alarm position and 0 befor
  series_infrastructure <- dplyr::case_when(
    base::is.na(alarm_position) ~ base::rep(0, length(alarm)),
    TRUE ~ base::ifelse(seq_along(alarm) >= (alarm_position + 1), 1, 0))
  
  # extract water level for each household
  list_household_model_hazard <- map2(list_household_model_development, seq_along(series_infrastructure), function(household, i) {
    depth <- if (series_infrastructure[i] == 0) {
      terra::extract(list_floodseries[[i]], vect(household$geometry), ID = FALSE)
    } else {
      terra::extract(list_floodseries_cp[[i]], vect(household$geometry), ID = FALSE)
    }
    depth[is.na(depth)] <- 0
    household$flood_return_period_year_current <- series_return_period[[i]]
    household$flood_depth_year_current <- depth[[1]]
    return(household)
  })
  
  return(list_household_model_hazard)
}
model_reduction_hazard <- function(list_household_model_development, ua_2) { 
  list_household_model_development <- list_household_inter; ua_2 <- 100
  
  alarm <- map_int(series_return_period, ~ ifelse(.x >= ua_2, 1, 0)) # build alarm series_return_period
  alarm_position <- base::which(alarm == 1)[1] # find alarm position
  
  # build series_return_period with 1 from alarm position and 0 befor
  series_infrastructure <- dplyr::case_when(
    base::is.na(alarm_position) ~ base::rep(0, length(alarm)),
    TRUE ~ base::ifelse(seq_along(alarm) >= (alarm_position + 1), 1, 0))
  
  
  df_household <- list_household_model_development[[1]] %>%
    st_as_sf() #%>%
    st_set_crs(32718)
  flood_layer <- list_floodseries[[1]] %>%
    st_set_crs(32718)
  depth_value <- st_join(df_household, flood_layer, left = TRUE) #%>%
   # pull(max_depth)
  
  # extract water level for each household
  list_household_model_hazard <- map2(list_household_model_development, seq_along(series_infrastructure), function(household, i) {
    depth <- if (series_infrastructure[i] == 0) {
      terra::extract(list_floodseries[[i]], vect(household$geometry), ID = FALSE)
    } else {
      terra::extract(list_floodseries_cp[[i]], vect(household$geometry), ID = FALSE)
    }
    depth[is.na(depth)] <- 0
    household$flood_return_period_year_current <- series_return_period[[i]]
    household$flood_depth_year_current <- depth[[1]]
    return(household)
  })
  
  return(list_household_model_hazard)
}
model_reduction_hazard <- function(list_household_model_development, ua_2 = 100) {
  #list_household_model_development <- list_household_inter; ua_2 <- 100
  
  # construir serie de alarma
  alarm <- map_int(series_return_period, ~ ifelse(.x >= ua_2, 1, 0))
  alarm_position <- base::which(alarm == 1)[1]
  
  # definir series_infrastructure con 1 desde alarm_position en adelante
  series_infrastructure <- dplyr::case_when(
    base::is.na(alarm_position) ~ base::rep(0, length(alarm)),
    TRUE ~ base::ifelse(seq_along(alarm) >= (alarm_position + 1), 1, 0)
  )
  
  # recorrer cada household y extraer max_depth
  list_household_model_hazard <- map2(list_household_model_development, seq_along(series_infrastructure), function(household, i) {
    
    # convertir household a sf si no lo es
    if (!inherits(household, "sf")) {
      household <- st_as_sf(household)
    }
    
    # seleccionar la capa de inundación
    flood_layer <- if (series_infrastructure[i] == 0) {
      list_floodseries[[i]]
    } else {
      list_floodseries_cp[[i]]
    }
    
    # asegurar que ambas capas tienen el mismo CRS
    household <- st_transform(household, st_crs(flood_layer))
    
    # unir household con flood_layer y extraer max_depth
    household <- st_join(household, flood_layer, left = TRUE)
    
    # reemplazar NA por 0 en max_depth
    household$flood_depth_year_current <- ifelse(is.na(household$max_depth), 0, household$max_depth)
    
    # agregar el periodo de retorno
    household$flood_return_period_year_current <- series_return_period[[i]]
    
    return(household)
  })
  
  return(list_household_model_hazard)
}

# model: risk ----
model_risk <- function(list_household_model_vulnerability_hazard, f, learning, r) { 
  #list_household_model_vulnerability_hazard <- list_household_inter
  
  df_household_risk <- list_household_model_vulnerability_hazard %>%
    bind_rows() 
  
  list_household_risk <- list()
  for (id_unique in unique(df_household_risk$household_id_model)) {
    
    df <- df_household_risk %>%
      filter(household_id_model == id_unique) %>%
      arrange(year_current) %>% 
      mutate(flood_return_period_max = flood_return_period_year_current,
             flood_depth_max = flood_depth_year_current)
    
    for (i in seq_len(nrow(df))) {
      
      df$household_measure_flood_proof[i] = df$household_compliance_rate[i] == 1   # Define si el hogar se "levanta" o no
      df$household_flood_year_current[i] = case_when(
        df$flood_depth_year_current[i]  > f ~ TRUE,
        df$flood_depth_year_current[i]  == 0 ~ FALSE,
        df$flood_depth_year_current[i]  < f & !df$household_measure_flood_proof[i] ~ TRUE,
        df$flood_depth_year_current[i]  < f & df$household_measure_flood_proof[i] ~ FALSE
      )
      
      if (i <= nrow(df) - 1) { 
        df$household_compliance_rate[i + 1] = case_when (
          # status quo
          df$household_fb[i] == 1 ~ df$household_compliance_rate[i],                        
          # learning effect
          df$household_fb[i] == 2 &  df$household_flood_year_current[i] ~ df$household_compliance_rate[i] + learning, 
          df$household_fb[i] == 2 & !df$household_flood_year_current[i] ~ df$household_compliance_rate[i],                        
          # levee effect
          df$household_fb[i] == 3 &  df$household_flood_year_current[i] ~ df$household_compliance_rate[i] + learning,                        
          df$household_fb[i] == 3 & !df$household_flood_year_current[i] ~ df$household_compliance_rate[i] - learning/r,  
          # good students
          df$household_fb[i] == 4 &  df$household_flood_year_current[i] ~ df$household_compliance_rate[i] + learning,  
          df$household_fb[i] == 4 & !df$household_flood_year_current[i] ~ df$household_compliance_rate[i] + learning/r )
        
        df$household_compliance_rate[i + 1] = pmax(0, pmin(1, df$household_compliance_rate[i + 1]))
        df$flood_depth_max <- cummax(df$flood_depth_year_current)
        df$flood_return_period_max <- cummax(df$flood_return_period_year_current)
        df$household_learning_year_current[i] <- df$household_compliance_rate[i+1] - df$household_compliance_rate[i]
        df$household_learning_accumulated[i] <- sum(df$household_learning_year_current[1:(i + 1)], na.rm = TRUE)
      }
    }
    
    df$household_flood_once <- cumany(df$household_flood_year_current)
    df$household_flood_counter <- cumsum(df$household_flood_year_current)
    list_household_risk <- c(list_household_risk, list(df))
  }
  
  list_household_risk <- list_household_risk %>%
    bind_rows() %>%
    split(.$year_current)
  return(list_household_risk)
}

# model: integrated ----
model_integrated <- function(list_development_area, n_switch, n, df_household, pfb, f, learning, r) {  
  #n <- n; n_switch <- c(10,25); #h <- household_initial; m <- household_growth; b <- 10 ;
  #household_compliance_threshold <- 0.7; ua_2 <- 100 ; f <- 0.5 ; learning <- 0.3 ; r <- 7 ; df_household <- df_household_calibration
  #n_switch <- c(10,25) ;   pfb <- c(0.25, 0.25, 0.25, 0.25) ;   df_household <- df_household_calibration
  #f <- 0.5 ; learning <- 0.3 ; r <- 7
  
  # verification
  if (base::sum(pfb) != 1) {stop("The proportions in 'pfb' must sum 1")}
  #if (household_compliance_threshold >  1) {stop("Threshold at household level must be less than 1")}
  #if (household_compliance_threshold <  0) {stop("Threshold at household level must be at less 0")}
  #if (m >= 0.1 * h) {stop("Very high growth")}
  
  # run model_integrated
  list_household_inter <- list_development_area %>%
    model_development(., n_switch, n, df_household, pfb) %>%
    model_reduction_vulnerability(.) %>% 
    model_reduction_hazard(.) %>%
    model_risk(., f, learning, r)
  
  household_flood_per_year <- (unlist(map_dbl(list_household_inter, ~ sum(.x$household_flood_year_current == TRUE))) / unlist(map_dbl(list_household_inter, ~ nrow(.x)))) * 100
  household_flood_accumulated <- (unlist(map_dbl(list_household_inter, ~ sum(.x$household_flood_once == TRUE))) / unlist(map_dbl(list_household_inter, ~ nrow(.x)))) * 100
  household_measure_flood_proof <- (unlist(map_dbl(list_household_inter, ~ sum(.x$household_measure_flood_proof == TRUE))) / unlist(map_dbl(list_household_inter, ~ nrow(.x)))) * 100
  
  return_list <- list(list_household_inter, household_flood_per_year, household_flood_accumulated, household_measure_flood_proof)
  names(return_list) <- c('list_household_inter', 'household_flood_per_year', 'household_flood_accumulated', 'household_measure_flood_proof')
  return(return_list)
}

# calibration ----
# parameters 
begin <- Sys.time()
n <- n; n_switch <- c(15,25) ; pfb=c(0,1,0,0) ; f <- 0.5 ; learning <- 0.2 ; r <- 7 
#household_compliance_threshold <- 1; ua_2 <- 100 ; 
list_model_calibration <- model_integrated(list_development_area, n_switch, n, df_household_calibration %>% filter(year_current==2024), 
                                           pfb, f, learning, r)
end <- Sys.time()
end-begin
n <- n; n_switch <- c(15,25) ; pfb=c(1,0,0,0) ; f <- 0.5 ; learning <- 0.2 ; r <- 7 
#household_compliance_threshold <- 1; ua_2 <- 100 ; 
list_model_calibration <- model_integrated(list_development_area, n_switch, n, df_household_calibration %>% filter(year_current==2024), 
                                           pfb, f, learning, r)

compare_observed_simulated <- function(df_observed, list_model, prefix) {
  #df_observed <- df_household_calibration; list_model <- list_model_calibration ;prefix <- "calibration"
  #df_observed <- df_household_validation; list_model <- list_model_validation ;prefix <- "validation"
  
  # observed simulated data frame -
  df_observed_simulated <- bind_rows(
    # observed
    df_observed %>% mutate(source = "observed"),
    # simulated
    map_dfr(c(2017, 2024), ~ list_model$list_household_inter[[as.character(.x)]] %>%
              st_as_sf() %>%
              mutate(source = "simulated"))
  ) %>%
    group_by(household_address, year_current) %>%
    filter(n() > 1) %>%
    ungroup()
  
  # observed simulated confusion matrix function
  # - generate confusion matrix plot
  # - generate confusion matrix map
  conf_matrix <- function(df_observed_simulated, year, variable, label) {
    #df_observed_simulated <- df_observed_simulated; year <- 2017; variable <- "household_flood_once";label <- "flood once"
    
    # filter data, adapt df to evaluate selected variable
    df <- df_observed_simulated %>%
      filter(year_current == year) %>% # filter by year
      distinct(geometry, source, .keep_all = TRUE) %>% # keep only unique geometry per source
      pivot_wider(
        id_cols = "geometry",
        names_from = "source", 
        # variable to see in popup in map view
        values_from = c(variable, "year_current", "household_year_entry", "household_compliance_rate",
                        "household_learning_accumulated", "household_flood_counter", "household_measure_flood_proof", 
                        "household_measure", "household_flood_year_current", "household_flood_once", "flood_depth_year_current",
                      "flood_depth_max", "flood_return_period_year_current", "flood_return_period_max")
      ) %>%
      mutate(!!sym(paste0(variable, "_observed")) := factor(!!sym(paste0(variable, "_observed")), levels = c(TRUE, FALSE))) %>%
      mutate(!!sym(paste0(variable, "_simulated")) := factor(!!sym(paste0(variable, "_simulated")), levels = c(TRUE, FALSE)))
    
    
    # confusion matrix
    cm <- conf_mat(df, truth = !!sym(paste0(variable, "_observed")), estimate = !!sym(paste0(variable, "_simulated")))
    
    # metrics to evaluate performance
    metrics <- cm %>%
      summary() %>%
      filter(.metric %in% c("accuracy", "f_meas","kap")) %>%
      select(.metric, .estimate) %>%
      pivot_wider(names_from = .metric, values_from = .estimate)
    
    # sub title plot
    title_text <- paste0(label, ": ", year, "\n",
                         "accuracy: ", round(metrics$accuracy, 2), " | ",
                         "f-1: ", round(metrics$f_meas, 2), " | ",
                         "kappa: ", round(metrics$kap, 2))
    
    # plot
    confusion_matrix <- 
      autoplot(cm, type = "heatmap") +
      scale_fill_gradient(low = "#D6EAF8", high = "dodgerblue4") +  
      labs(
        title = title_text,
        x = "Observed",
        y = "Simulated",
        fill = "Count"
      ) +
      theme_minimal() 
    
    
    # add category to df
    df <- df %>%
      mutate(
        conf_matrix_category = case_when(
          !!sym(paste0(variable, "_observed")) == TRUE  & !!sym(paste0(variable, "_simulated")) == TRUE  ~ "True Positive (TP)",
          !!sym(paste0(variable, "_observed")) == FALSE & !!sym(paste0(variable, "_simulated")) == FALSE ~ "True Negative (TN)",
          !!sym(paste0(variable, "_observed")) == FALSE & !!sym(paste0(variable, "_simulated")) == TRUE  ~ "False Positive (FP)",
          !!sym(paste0(variable, "_observed")) == TRUE  & !!sym(paste0(variable, "_simulated")) == FALSE ~ "False Negative (FN)",
          TRUE ~ NA
        )
      )
    
    # visualize map
    map <- mapview(df %>% 
                     filter(!is.na(conf_matrix_category)) %>%
                     mutate(across(matches("observed|simulated"), as.character)) %>%
                     mutate(conf_matrix_category = factor(conf_matrix_category, 
                                                          levels = c("True Positive (TP)", 
                                                                     "True Negative (TN)", 
                                                                     "False Positive (FP)", 
                                                                     "False Negative (FN)"))), 
                   zcol = "conf_matrix_category",
                   col.regions = c(
                     "True Positive (TP)" = "green4",
                     "True Negative (TN)" = "green1",
                     "False Positive (FP)" = "red4",
                     "False Negative (FN)" = "red1"),
                   na.color = "gray",
                   layer.name = variable,
                   popup = popupTable(df, 
                                      zcol = grep("observed|simulated", names(df), value = TRUE)))
    map
    return(list(plot = confusion_matrix, map = map))
  }
  
  # generate confusion matrix plots for different variables and years
  plot1 <- conf_matrix(df_observed_simulated, 2017, "household_flood_once", "flood once")$plot
  plot2 <- conf_matrix(df_observed_simulated, 2024, "household_flood_year_current", "flood current year")$plot
  plot3 <- conf_matrix(df_observed_simulated, 2024, "household_measure_flood_proof", "flood proof measures")$plot
  
  plot(plot1)
  plot(plot2)
  plot(plot3)
  
  # generate confusion matrix plots for different variables and years
  map1 <- conf_matrix(df_observed_simulated, 2017, "household_flood_once", "flood once")$map
  map2 <- conf_matrix(df_observed_simulated, 2024, "household_flood_year_current", "flood current year")$map
  map3 <- conf_matrix(df_observed_simulated, 2024, "household_measure_flood_proof", "flood proof measures")$map
  
  # save plot
  png(paste0("plot/", prefix, "_conf_matrix.png"), width = 2200, height = 3400, res = 300)
  grid.arrange(plot1, plot2, plot3, ncol = 1)
  dev.off()
  
  
  # observed simulated general answer -
  df_comparison <- df_observed_simulated %>%
    sf::st_drop_geometry() %>%  # delete geoemtry
    group_by(year_current, source) %>%  # group by year and source
    summarise(
      household_flood_once = mean(household_flood_once, na.rm = TRUE) * 100,
      household_flood_year_current = mean(household_flood_year_current, na.rm = TRUE) * 100,
      household_measure_flood_proof = mean(household_measure_flood_proof, na.rm = TRUE) * 100,
      household_measure = mean(household_measure, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(household_flood_once, household_flood_year_current, household_measure_flood_proof, household_measure),
      names_to = "variable",
      values_to = "percentage"
    ) %>%
    pivot_wider(names_from = source, values_from = percentage) %>%
    mutate(error = simulated - observed)
  
  # generate plot 
  data_combined <- bind_rows(
    data.frame(year_current = seq(year_initial, year_end), 
               household_flood_per_year = list_model_calibration$household_flood_per_year, 
               household_flood_accumulated = list_model_calibration$household_flood_accumulated, 
               household_measure_flood_proof = list_model_calibration$household_measure_flood_proof, 
               household_flood_behaviour = "true proportions"))
  
  data_labels <- data_combined %>%
    group_by(household_flood_behaviour) %>%
    filter(year_current == max(year_current))
  {
    x_axis <- scale_x_continuous(
      breaks = seq(year_initial, year_end, by = 5),
      minor_breaks = seq(year_initial, year_end, by = 1),
      expand = c(0.07,0) 
    )
    plot_flood_return_period <- 
      ggplot( data.frame(year_simulation = year_simulation, caudales = series_return_period), aes(x = year_simulation, y = caudales)) +
      geom_bar( stat = "identity") +
      labs( x = "Year", y = "Return period") +
      theme_classic() +
      theme(legend.position = "none", 
            axis.title.x = element_blank()) +
      x_axis + 
      scale_y_continuous(limits = c(0, 100), expand = c(0, 0))
    
    plot_household_flood_per_year <- 
      ggplot(data_combined, aes(x = year_current, y = household_flood_per_year, color = household_flood_behaviour)) +
      geom_line() +
      annotate("point", 
               x = 2024, 
               y = df_comparison %>% filter(year_current == 2024, variable == "household_flood_year_current") %>% 
                 pull(observed), 
               size = 2) +
      annotate("text", 
               x = 2024, 
               y = df_comparison %>% filter(year_current == 2024, variable == "household_flood_year_current") %>% 
                 pull(observed), 
               label = paste("survey 2024", "\nerror %:", round(df_comparison$error[df_comparison$year_current == 2024 & df_comparison$variable == "household_flood_year_current"], 0)),
               vjust = -0.5) +
      #geom_text_repel(data = data_labels,
      #aes(label = household_flood_behaviour), 
      #                nudge_x = 0.5, direction = "y", show.legend = FALSE) +
      labs(y = "% Flood per year") +
      x_axis +
      scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
      theme_classic() +
      theme(legend.position = "none", 
            axis.title.x = element_blank(), 
            #axis.text.x = element_blank(),       
            #axis.ticks.x = element_blank()
      )
    
    
    plot_household_flood_accumulated <- 
      ggplot(data_combined, aes(x = year_current, y = household_flood_accumulated, color = household_flood_behaviour)) +
      geom_line() +
      annotate("point", x = 2017, 
               y = df_comparison %>% filter(year_current == 2017, variable == "household_flood_once") %>% 
                 pull(observed)) +
      annotate("text", x = 2017, 
               y = df_comparison %>% filter(year_current == 2017, variable == "household_flood_once") %>% 
                 pull(observed), 
               label = paste("survey 2017", "\nerror %:", round(df_comparison$error[df_comparison$year_current == 2017 & df_comparison$variable == "household_flood_once"], 0)), 
               vjust = -0.5) +
      #geom_text_repel(data = data_labels, 
      #                #aes(label = household_flood_behaviour), 
      #                nudge_x = 0.5, direction = "y", show.legend = FALSE) +
      #geom_text_repel(data = data_labels, aes(label = household_flood_behaviour), 
      #                nudge_x = 0.5, direction = "y", show.legend = FALSE) +
      labs(y = "% Flood accumulated") +
      x_axis + 
      theme_classic() +
      scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
      theme(legend.position = "none",               
            axis.title.x = element_blank(),
            #axis.text.x = element_blank(),       
            #axis.ticks.x = element_blank()
      )
    
    plot_household_measures_flood_proof <- 
      ggplot(data_combined, aes(x = year_current, y = household_measure_flood_proof, color = household_flood_behaviour)) +
      geom_line() +
      annotate("point", 
               x = 2024, 
               y = df_comparison %>% filter(year_current == 2024, variable == "household_measure_flood_proof") %>% 
                 pull(observed)) +
      annotate("text", 
               x = 2024, 
               y = df_comparison %>% filter(year_current == 2024, variable == "household_measure_flood_proof") %>% 
                 pull(observed), 
               label = paste("survey 2024", "\nerror %:", round(df_comparison$error[df_comparison$year_current == 2024 & df_comparison$variable == "household_measure_flood_proof"], 0)), 
               vjust = -0.5) +
      #geom_text_repel(data = data_labels, 
      #                #aes(label = household_flood_behaviour), 
      #               nudge_x = 0.5, direction = "y", show.legend = FALSE) +
      labs(y = "% Flood-proof-measures") +
      x_axis +
      theme_classic() +
      scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
      theme(legend.position = "none",           
            axis.title.x = element_blank(),
            #axis.text.x = element_blank(),       
            #axis.ticks.x = element_blank()
      )
    
    plot_household_measures <- 
      ggplot(data_combined, aes(x = year_current, y = household_measure_flood_proof, color = household_flood_behaviour)) +
      geom_line() +
      annotate("point", x = 2017, 
               y = df_comparison %>% filter(year_current == 2017, variable == "household_measure") %>% 
                 pull(observed)) +
      annotate("text", x = 2017, 
               y = df_comparison %>% filter(year_current == 2017, variable == "household_measure") %>% 
                 pull(observed), 
               label = paste("survey 2017", "\nerror %:", round(df_comparison$error[df_comparison$year_current == 2017 & df_comparison$variable == "household_measure"], 0)), 
               vjust = -0.5) +
      annotate("point", 
               x = 2024, 
               y = df_comparison %>% filter(year_current == 2024, variable == "household_measure") %>% 
                 pull(observed)) +
      annotate("text", 
               x = 2024, 
               y = df_comparison %>% filter(year_current == 2024, variable == "household_measure") %>% 
                 pull(observed), 
               label = paste("survey 2024", "\nerror %:", round(df_comparison$error[df_comparison$year_current == 2024 & df_comparison$variable == "household_measure_flood_proof"], 0)), 
               vjust = -0.5) +
      #geom_text_repel(data = data_labels, 
      #                #aes(label = household_flood_behaviour), 
      #                nudge_x = 0.5, direction = "y", show.legend = FALSE) +
      labs(y = "% Flood measures") +
      x_axis +
      theme_classic() +
      scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
      theme(legend.position = "none",           
            #axis.title.x = element_blank(),
            #axis.text.x = element_blank(),       
            #axis.ticks.x = element_blank()
      )
    
    }
  grid.arrange(plot_flood_return_period,
               plot_household_flood_per_year, 
               plot_household_flood_accumulated, 
               #plot_household_measures,
               plot_household_measures_flood_proof,
               ncol = 1) 
  
  
  # save plot
  png(paste0("plot/", prefix, "_general_answer.png"), width = 2400, height = 3400, res = 300)
  grid.arrange(plot_flood_return_period,
               plot_household_flood_per_year, 
               plot_household_flood_accumulated, 
               #plot_household_measures,
               plot_household_measures_flood_proof,
               ncol = 1,  widths = c(2)) 
  dev.off()
  
  # plot in viewer
  simplify_flood <- function(flood_layer) {
    flood_layer %>%
      select(max_depth) %>%
      filter(max_depth > 0) %>%
      na.omit() %>%
      vect() %>%
      terra::aggregate() %>%
      st_as_sf() # convertir de terra a sf
  }
  
  # flood_2023 y flood_2024
  flood_2003 <- simplify_flood(list_floodseries$'2003')
  flood_2024 <- simplify_flood(list_floodseries$'2024')
  
  map <- map1 + map2 + map3 + 
   mapview(flood_2003, col.regions = 'skyblue') + 
   mapview(flood_2024, col.regions = 'skyblue')
  map
  #return(df_observed_simulated)
}
compare_observed_simulated(df_household_calibration, list_model_calibration, "calibration")

# validation ----
list_model_validation <- model_integrated(list_development_area, n_switch, n, df_household_validation, 
                                          pfb, household_compliance_threshold, ua_2, f, learning, r)
compare_observed_simulated(df_household_validation, list_model_validation, "validation")


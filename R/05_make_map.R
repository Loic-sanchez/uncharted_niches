make_map = function(x) {
 
  source(here::here("R", "02_prepare_rasters.R"))
  source(here::here("R", "00_join_obitools.R"))
  source(here::here("R", "01_get_occus.R"))
  
  ever_sampled = terra::rast(here::here("outputs", "ever_sampled.asc"))
  focal_rasts = focal_rasters()
  
  r = terra::rast(extent = terra::ext(ever_sampled), 
                   resolution = 0.5)
  ever_sampled_0.5 = terra::aggregate(ever_sampled, fact = 5, fun = sum)
  # terra::plot(log10(ever_sampled_0.5))
  # 
  get_gbif_obis()
  load(here::here("outputs", "df_dna_meta.Rdata"))
  
  # Join obis/gbif, and get unique occurrences
  
  OBIS5 = OBIS4 |> 
    dplyr::select(species, latitude, longitude) 
  GBIF1 = unique_gbif_occus |> 
    dplyr::select(species, latitude, longitude)
  
  all_no_DNA = rbind(OBIS5,
                     GBIF1)
  
  unique_no_DNA = unique(all_no_DNA |> 
                           dplyr::select(longitude, latitude))
  
  # Filter DNA, unique etc
  
  eDNA4b = df_dna_meta  |>
    dplyr::rename(species = "species_name_corrected",
                  longitude = "longitude_start",
                  latitude = "latitude_start") |> 
    dplyr::select(species, latitude, longitude)
  
  eDNA5 = unique(eDNA4b |> 
                   dplyr::select(species, latitude, longitude) |> 
                   tidyr::drop_na(longitude)) |> 
    tidyr::drop_na(latitude)
  
  vec_sp = sort(unique(eDNA5$species))
  
  unique_DNA = unique(eDNA5 |> dplyr::select(longitude, latitude))
  
  rast_noDNA = terra::rasterize(as.matrix(unique_no_DNA), 
                                r, 
                                fun = "length")
  rast_noDNA[is.na(rast_noDNA[])] <- 0
  terra::plot(log10(rast_noDNA))
  
  ever_samp_vect = terra::as.points(ever_sampled_0.5)
  ever_samp_df = terra::as.data.frame(ever_sampled_0.5, xy = T)
  no_dna_df_qt = ever_samp_df |> 
    dplyr::mutate(Sampling_magnitude = ifelse(values_length %in% 1, 1,
                                       ifelse(values_length %in% 2:10, 2,
                                              ifelse(values_length %in% 11:100, 3, 
                                                     ifelse(values_length%in% 101:1000, 4,
                                                            ifelse(values_length %in% 1001:10000, 5, 0))))))
  library(ggplot2)
  
  ggplot() +  
    geom_tile(data = no_dna_df_qt,
              aes(x = x,
                  y = y,
                  fill = Sampling_magnitude),
              alpha = 1)  +
    # geom_point(data = dna_df,
    #            aes(x = x,
    #                y = y),
    #            color = "red",
    #            size = 1.2) +
    scale_fill_viridis_b(option = "D", direction = 1) +
    theme_minimal() +
    theme(legend.position = "bottom")

  }

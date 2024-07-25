join_obi = function() {
  
  list_dna = list.files(
    here::here("data", "DNA_data"),
    pattern = "*.csv",
    all.files = TRUE,
    full.names = TRUE
  )
  
  all_dna_csv = lapply(list_dna, function(i) {
    read.csv(i)
  })
  
  clean_all_dna = lapply(1:length(all_dna_csv), function(k) {
    get_csv = all_dna_csv[[k]]
    current_names <- unique(get_csv$species_name_corrected)
    valid_names <- rfishbase::validate_names(current_names)
    
    eDNAa = data.frame()
    
    for (j in 1:length(current_names)) {
      
      if (is.na(current_names[j])) {
        next  # Skip to the next iteration if item is NA
      }
      
      df <- get_csv %>%
        subset(species_name_corrected == current_names[j])
      df$species_name_corrected <- valid_names[j]
      eDNAa <- rbind(eDNAa, df)
      
    }
    
    eDNAb = eDNAa %>%
      filter(species_name_corrected != "Genus Species") %>%
      drop_na(species_name_corrected)
    
    cols_rep = colnames(eDNAb %>%
                          dplyr::select(all_of(contains(
                            "Nb_replicates"
                          ))))
    if (sum(stringr::str_detect(names(eDNAb), "sum_replicates")) >= 1) {
      eDNAc = eDNAb %>%
        dplyr::select(-contains("Nb_reads")) %>%
        dplyr::select(-c("sum_replicates", "sequence")) %>%
        pivot_longer(cols_rep) %>%
        filter(value > 0)
    }
    else{
      eDNAb %>%
        dplyr::select(-contains("Nb_reads")) %>%
        dplyr::select(-c("sequence")) %>%
        pivot_longer(cols_rep) %>%
        filter(value > 0)
    }
    
    eDNAc = eDNAb %>%
      dplyr::select(-contains("Nb_reads")) %>%
      dplyr::select(-c("sequence")) %>%
      pivot_longer(cols_rep) %>%
      filter(value > 0)
    
    
    eDNAc$value = 1
    
    eDNAd = unique(eDNAc) # drop double detections of a species
    
    eDNAd$code_spygen = eDNAd$name
    eDNAd$code_spygen = gsub(".Nb_replicates", "", eDNAd$code_spygen)
    
    return(eDNAd)
    
  })
  
  df_all_dna = data.table::rbindlist(clean_all_dna, fill = T) %>%
    select(species_name_corrected, code_spygen, value)
  
  meta_b = read.csv(here::here("data", "meta", "metadata_globales_190623.csv"), sep = ";")[, -1]
  meta_nisk = meta_b |> dplyr::filter(sample_method == "niskin" | sample_type == "niskin")
  spy_nisk = unique(meta_nisk$code_spygen)
  
  meta_b = meta_b |> dplyr::filter(!code_spygen %in% spy_nisk)
  
  meta2 = read.csv(here::here("data", "meta", "eDNA_Patagonia.csv"), sep = ",")
  
  no_autho_filt = c(
    "SPY222104",
    "SPY222102",
    "SPY222090",
    "SPY222096",
    "SPY222103",
    "SPY222105",
    "SPY222092",
    "SPY222091"
  ) # filters without legal authorization from the country
  
  meta3 = rbind(meta_b |>     
                  dplyr::filter(!code_spygen %in% no_autho_filt) |> 
                  drop_na(longitude_start) |>
                  filter(country != "NewCaledonia") |>
                  dplyr::select(code_spygen, 
                                longitude_start, 
                                latitude_start, 
                                project), 
                meta2 |> 
                  dplyr::rename(
                    code_spygen = "Code",
                    longitude_start = "Longitude_start",
                    latitude_start = "Latitude_start") |> 
                  dplyr::select(code_spygen,
                                longitude_start,
                                latitude_start, 
                                project))
  
  df_dna_meta = df_all_dna %>%
    left_join(meta3,
              by = "code_spygen")
  
  df_dna_meta$longitude_start = as.numeric(gsub(",", ".", df_dna_meta$longitude_start))
  df_dna_meta$latitude_start = as.numeric(gsub(",", ".", df_dna_meta$latitude_start))
  
  eco_fish = rfishbase::species(unique(df_dna_meta$species_name_corrected))
  salt_only = eco_fish %>% filter(Saltwater > 0)
  df_dna_meta = df_dna_meta %>% filter(species_name_corrected %in% salt_only$Species)

  save(df_dna_meta, file = here::here("outputs", "df_dna_meta.Rdata"))
  return(df_dna_meta)

}
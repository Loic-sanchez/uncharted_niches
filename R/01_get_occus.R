load(here::here("data", "Data.Rdata"))

download_gbif = function(){

vec_sp = sort(unique(df_dna_meta_conf$species)) # All sp. that were detected with DNA)

gbif_occurrences = lapply(1:length(vec_sp), function(i){
  
  tryCatch({sp_gbif = geodata::sp_occurrence(genus = stringr::word(vec_sp[i], 1),
                                             species = stringr::word(vec_sp[i], 2),
                                             removeZeros = T,
                                             geo = T,
                                             nrecs = 10^10,
                                             args = c("year = 2006, 2023"))}, error = function(e) sp_gbif <<- NULL)
  
  if(!is.null(sp_gbif)) {sp_gbif = sp_gbif %>%
    dplyr::select(species, lat, lon, year, occurrenceStatus, datasetKey) %>%
    filter(occurrenceStatus == "PRESENT") %>%
    dplyr::select(species, lat, lon, year, datasetKey) %>%
    rename(longitude = "lon", latitude = "lat")}
  
  cat(i, "\n")
  return(sp_gbif)
  
})

rbound_list = data.table::rbindlist(gbif_occurrences)
unique_gbif_occus = unique(rbound_list)

eco_fish = rfishbase::species(unique(unique_gbif_occus$species))
salt_only = eco_fish %>% filter(Saltwater > 0)
unique_gbif_occus = unique_gbif_occus %>% filter(species %in% salt_only$Species)

save(unique_gbif_occus, file = here::here("data", "OBIS_GBIF", "unique_gbif_occu_07.Rdata"))

}

download_obis = function(){

  library(worrms)
  aphias = worrms::wm_name2id_(vec_all1)
  aphias_df = cbind(vec_all, do.call(rbind.data.frame, aphias))
  colnames(aphias_df)[2] = "ID"
  
  #### Manual corrections for 11 sp
  
  aphias_df$ID[aphias_df$vec_sp == "Aphareus rutilans"] = 218468  
  aphias_df$ID[aphias_df$vec_sp == "Balistes capriscus"] = 154721
  # aphias_df$ID[aphias_df$vec_sp == "Chaetodon lunulatus"] = 398549
  aphias_df$ID[aphias_df$vec_sp == "Harengula humeralis"] = 367197
  # aphias_df$ID[aphias_df$vec_sp == "Gerres longirostris"] = 367242
  aphias_df$ID[aphias_df$vec_sp == "Gobius niger"] = 126892  
  aphias_df$ID[aphias_df$vec_sp == "Malacocephalus laevis"] = 272392  
  # aphias_df$ID[aphias_df$vec_sp == "Mobula tarapacana"] = 105859
  aphias_df$ID[aphias_df$vec_sp == "Naso brevirostris"] = 219671    
  # aphias_df$ID[aphias_df$vec_sp == "Plotosus lineatus"] = 217659  
  aphias_df$ID[aphias_df$vec_sp == "Pterois volitans"] = 159559  
  aphias_df$ID[aphias_df$vec_sp == "Raja brachyura"] = 367297  
  aphias_df$ID[aphias_df$vec_sp == "Sphyraena barracuda"] = 345843  
  aphias_df$ID[aphias_df$vec_sp == "Symphodus ocellatus"] = 273572  
  aphias_df$ID[aphias_df$vec_sp == "Trachinotus falcatus"] = 367285  
  
  np = !is.null(curl::nslookup("r-project.org", error = FALSE))
  assign("has_internet_via_proxy", np, environment(curl::has_internet))
  
  results <- purrr::imap(aphias_df$ID, function(x, i) {
    
    message(i, "/", nrow(aphias_df))
    
    robis::occurrence(taxonid = x)
    
  })
  
  OBIS1 <- dplyr::bind_rows(results)
  
  OBIS2 <- OBIS1 |>  
    dplyr::select(aphiaID, species, decimalLongitude, decimalLatitude, eventDate) |> 
    dplyr::left_join(aphias_df, by = c("aphiaID" = "ID"))
  
  OBIS3 = unique(OBIS2)
  OBIS3b <- OBIS3 |> 
    dplyr::mutate(year = as.numeric(substr(eventDate, 1, 4))) |> 
    tidyr::drop_na(year) |> 
    dplyr::filter(year < 2024 & year > 1742)
  
  current_names <- unique(OBIS3$species)
  valid_names <- rfishbase::validate_names(current_names)
  
  OBIS3b = data.frame()
  
  for (j in 1:length(current_names)) {
    
    df <- OBIS3 |> 
      subset(species == current_names[j])
    df$species <- valid_names[j]
    
    OBIS3b <- rbind(OBIS3b, df)
    
  }
  
  # Get unique species names
  unique_species <- unique(OBIS3$species)
  
  # Validate all unique species names
  validated_species <- rfishbase::validate_names(unique_species)
  
  # Create a mapping data frame
  name_mapping <- data.frame(species = unique_species, valid_species = validated_species)
  
  # Merge the original data frame with the mapping data frame
  
  OBIS3b <- OBIS3 |> 
    dplyr::inner_join(name_mapping, by = "species") |> 
    dplyr::select(-species) |> 
    dplyr::rename(species = valid_species)
  
  OBIS4 = OBIS3b |>  
    dplyr::rename(longitude = "decimalLongitude", latitude = "decimalLatitude")
  
  eco_fish = rfishbase::species(unique(OBIS4$species))
  salt_only = eco_fish |> dplyr::filter(Saltwater > 0)
  OBIS4 = OBIS4 |> dplyr::filter(species %in% salt_only$Species)
  save(OBIS4, file = here::here("data", "OBIS_GBIF", "obis_extract.RData"))
  
}

get_gbif_obis = function(){
  
  list_files = list.files(here::here("data", "OBIS_GBIF"), 
                          pattern = "*RData",
                          full.names = T,
                          all.files = T) 
  lapply(list_files, 
         load, 
         .GlobalEnv)
  
}

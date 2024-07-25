library(rfishbase)

load(here::here("data", "GAPeDNA", "data_for_GAPeDNA.Rdata"))
list_sp_cefe_eth = read.delim(here::here("data", 
                                         "morgane_seq", 
                                         "db_CEFE_ETH_species_list_23.txt"),
                              header = F)

sp_cefe_eth = list_sp_cefe_eth$V1
val12S = c(list_ecopcr_df$Valentini_12S$species_name, sp_cefe_eth)
val12S = gsub("_", " ", val12S)

current_names <- unique(val12S)

get_valid_names <- function(current_names) {
  valid_names <- rfishbase::validate_names(current_names)
  return(valid_names)
}

valid_names_list <- list()

for (current_name in current_names) {
  valid_names_list[[current_name]] <- get_valid_names(current_name)
}

clean_list <- valid_names_list[!sapply(valid_names_list, function(x) any(is.na(x)))]
clean_names <- as.data.frame(do.call(rbind, clean_list)) |> dplyr::rename(valid_names = "V2")

# region check list 

all_occurences_list$Species_name = as.character(all_occurences_list$Species_name)
all_occurences_list$Species_name = gsub("_", " ", all_occurences_list$Species_name)

current_names <- unique(all_occurences_list$Species_name)
valid_names_list <- list()

for (current_name in current_names) {
  valid_names_list[[current_name]] <- get_valid_names(current_name)
}

clean_checklist = valid_names_list[!sapply(valid_names_list, function(x) any(is.na(x)))]
clean_names_checklist <- as.data.frame(do.call(rbind, clean_checklist)) |> 
  dplyr::rename(valid_names = "V2")
clean_names_checklist$V1 = rownames(clean_names_checklist) 
# save(clean_names_checklist, file = here::here("output", "clean_names_checklist.Rdata"))

all_occurrences_cor = all_occurences_list |> 
  dplyr::left_join(clean_names_checklist, by = c("Species_name" = "V1"), multiple = "first")

check_sp = df_dna_meta |> 
  tidyr::drop_na(longitude_start) |> 
  dplyr::filter(species_name_corrected %in% clean_names$valid_names)
unique(check_sp$species_name_corrected)

not_base = unique(df_dna_meta |> 
                    tidyr::drop_na(longitude_start) |> 
                    dplyr::anti_join(check_sp, by = "species_name_corrected"))
not_base_sp <- c(unique(not_base$species_name_corrected), "Lethotremus muticus")

check_12S_sequence <- function(species_name) {
  # Define search term
  search_term <- paste(species_name, "[organism]", "12S", sep=" ")
  database <- "nucleotide"
  
  # Search NCBI using Entrez
  search_results <- entrez_search(db=database, term=search_term)
  
  # Check if any results were returned
  if (search_results$count > 0) {
    return(TRUE)  # Sequence is available
  } else {
    return(FALSE) # Sequence is not available
  }
}

check_12S_genbank <- function(species_name) {
  # Define search term
  search_term <- paste(species_name, "[organism]", "12S", sep=" ")
  database <- "nucleotide"
  
  # Search NCBI using Entrez
  search_results <- reutils::esearch(db=database, term=search_term)
  # Check if any results were returned
  if (is.na(search_results[1])) {
    return(F)  # Sequence is available
  } else {
    return(T) # Sequence is not available
  }
}

# Check availability of 12S sequences for each species

availability <- sapply(not_base_sp, check_12S_genbank)

# Print results

NCBI <- data.frame(Species = not_base_sp, Available = availability)

# Provinces list 

provs = read.csv(here::here("data", "GAPeDNA", "provinces.csv"), sep = ";")
provs$Project = as.factor(provs$Project)
meds = c("Gombessa 6+", "CERCA", "DeepHeart", "eRef", "corse_OFB", "corse_OEC", 
         "ANGE", "Corse_OEC", "PISCIS", "med", "gombessa")

df_dna_meta = df_dna_meta |> 
  tidyr::drop_na(longitude_start) |> 
  dplyr::mutate(Project = ifelse(project %in% meds, "med", project))

unique(df_dna_meta$Project)
df_dna_meta$Project = as.factor(df_dna_meta$Project)

# Get project names 

vec_project = unique(df_dna_meta$Project)

# Loop over 'em

loop_project = lapply(1:length(vec_project), function(i){

df_proj_i = df_dna_meta |>
  dplyr::select(species_name_corrected, Project) |> 
  dplyr::filter(Project == vec_project[i]) # Select project i

vec_species = unique(df_proj_i$species_name_corrected) # What species in this project ?

loop_sp = lapply(1:length(vec_species), function(j) { # Loop over these species

genus = stringr::word(vec_species[j], 1) # Genus of species j

# Filter this species, and join the provinces of interest for the project

df_sp_j = unique(df_proj_i |> 
  dplyr::filter(species_name_corrected == vec_species[j]) |> 
  dplyr::left_join(provs, by = dplyr::join_by(Project), relationship = "many-to-many"))

# Regional checklist filtered by these provinces

occu_province = all_occurrences_cor |> 
  dplyr::filter(spatial_resolution == "Provinces") |> 
  dplyr::filter(RegionName %in% df_sp_j$Main.provinces)

# Obtain all species in this genus

occu_genus = occu_province |> 
  dplyr::filter(stringr::str_detect(valid_names, genus))

occu_genus = all_occurrences_cor |>
  dplyr::filter(stringr::str_detect(valid_names, genus))

if(nrow(occu_genus) == 0){return(NULL)}

# list 'em

genus_in_prov = unique(occu_genus$valid_names)

# Which sp. are sequenced in this genus ?

genus_seq_BG = sapply(genus_in_prov, check_12S_genbank)
genus_sequenced_BG = genus_in_prov[genus_seq_BG]

genus_seq_custom = unique(clean_names |> 
                            dplyr::filter(valid_names %in% genus_in_prov))

genus_all_bases = unique(c(genus_seq_custom$valid_names, genus_sequenced_BG))

# nb_tot_genus_province = sum(stringr::str_count(occu_province$valid_names, genus))
nb_tot_genus_province = length(unique(occu_genus$valid_names))
nb_tot_genus_db = length(genus_all_bases)

completeness = nb_tot_genus_db/nb_tot_genus_province

completeness_df = data.frame(Species = vec_species[j], 
                             Project = vec_project[i], 
                             Completeness = completeness)

cat(j, "\n")
return(completeness_df)

})

completeness_df = data.table::rbindlist(loop_sp)
cat(i, "\n")
return(completeness_df)

})

completeness_projs = data.table::rbindlist(loop_project)
save(completeness_projs, file = here::here("outputs", "completeness_projs.Rdata"))

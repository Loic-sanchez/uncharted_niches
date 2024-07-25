source(here::here("R", "01_get_occus.R"))

load(here::here("outputs", "SPY_extract.Rdata"))
load(here::here("outputs", "df_dna_meta.Rdata"))
load(here::here("outputs", "df_dna_meta_conf.Rdata"))
load(here::here("outputs", "completeness_projs.Rdata"))

get_gbif_obis()

grid_spy_dists = function(){
  
  rast_noDNA_spy = terra::rasterize(as.matrix(all_no_DNA |> dplyr::select(longitude, latitude)), 
                                    r, 
                                    fun = "length", 
                                    background = 0)  
  
  SPY_coords = unique(df_dna_meta %>% 
                        dplyr::select(code_spygen, longitude_start, latitude_start))
  
  SPY_points = terra::vect(SPY_coords, geom = c("longitude_start", "latitude_start"))
  
  SPY_extract_dists = as.data.frame(terra::extract(rast_noDNA_spy,
                                                   SPY_points, 
                                                   cells = T, 
                                                   bind = T))
  
  save(SPY_extract_dists, file = here::here("outputs", "SPY_extract_dists.Rdata"))
  
}

r_focal = terra::rast(here::here("data", "treated_rasters", "focal_rasts.tif"))
canals = terra::vect(here::here("data", "poly_canal", "detroits.shp"))

OBIS5 = OBIS4 |> 
  dplyr::select(species, latitude, longitude)
GBIF1 = unique_gbif_occus |> 
  dplyr::select(species, latitude, longitude)

all_no_DNA = rbind(OBIS5,
                   GBIF1)

# Filter DNA, unique etc

meds = c("Gombessa 6+", "CERCA", "DeepHeart", "eRef", "corse_OFB", "corse_OEC", 
         "ANGE", "Corse_OEC", "PISCIS", "med", "gombessa")

df_dna_meta = df_dna_meta |> 
  tidyr::drop_na(longitude_start) |> 
  dplyr::mutate(Project = ifelse(project %in% meds, "med", project))

safe_sharkies = c("Isurus oxyrinchus", 
                  "Lamna nasus", 
                  "Carcharinus obscurus",
                  "Alopias pelagicus", 
                  "Alopias vulpinus", 
                  "Alopias superciliosus",
                  "Odontaspis ferox",
                  "Mobula mobular", 
                  "Mobula thurstoni", 
                  "Mobula birostris",
                  "Mobula tarapacana")

df_dna_meta_isu = df_dna_meta |> dplyr::filter(species_name_corrected %in% safe_sharkies)
df_dna_meta_isu$Completeness = 1

df_dna_meta_conf = df_dna_meta |> 
  dplyr::left_join(completeness_projs,
                   by = c("species_name_corrected" = "Species", 
                          "Project")) |> 
  dplyr::filter(Completeness >= 0.95)

df_dna_meta_conf = rbind(df_dna_meta_conf, df_dna_meta_isu)

eDNA4b = df_dna_meta_conf |> 
  dplyr::filter(Completeness >= 0.95) |> 
  dplyr::rename(species = "species_name_corrected",
                longitude = "longitude_start",
                latitude = "latitude_start") |> 
  dplyr::select(species, latitude, longitude)

eDNA5 = unique(eDNA4b |> 
                 dplyr::select(species, latitude, longitude) |>  
                 tidyr::drop_na(longitude)) |> 
  tidyr::drop_na(latitude)

vec_sp = sort(unique(eDNA5$species))

r = terra::rast(nrows = 1800, ncols = 3600, crs = terra::crs(r_focal))

r_focal = terra::project(r_focal, r)

longlat = "+proj=longlat +datum=WGS84"
terra::crs(r) = longlat
terra::crs(r_focal) = longlat

all_sp_dists = lapply(1:length(vec_sp), function(i){

  # Filter by species i
  
  eDNA6 = eDNA5 |>  
    dplyr::filter(species == vec_sp[i])
  
  all_no_DNA1 = all_no_DNA |>  
    dplyr::filter(species == vec_sp[i])
  
  unique_no_DNA = unique(all_no_DNA1 |> dplyr::select(longitude, latitude))
  unique_DNA = unique(eDNA6 |> dplyr::select(longitude, latitude))
  
  if(nrow(unique_DNA) < 1) {return(NULL)}
  if(nrow(unique_no_DNA) < 1) {return(NULL)}
  
  # Create a raster of presence/absence 
  
  rast_noDNA = terra::rasterize(as.matrix(unique_no_DNA), 
                                r, 
                                fun = "length", 
                                background = 0)

  rast_DNA = terra::rasterize(as.matrix(unique_DNA),
                              r, 
                              fun = 'length')
  
  rast_noDNA_mask_pres = rast_noDNA
  rast_noDNA_mask_pres[rast_noDNA_mask_pres > 0] = 1
  mask_focal = terra::focal(r_focal$NPP_median, w = 3, fun = mean, na.rm = T)
  
  rast_noDNA_mask_pres = terra::mask(rast_noDNA_mask_pres,
                                     mask_focal,
                                     inverse = F)
  
  canal_r = terra::rasterize(canals, rast_noDNA_mask_pres)
  rast_noDNA_mask_pres[canal_r == 1] = 0 # Panama 
  rast_noDNA_mask_pres[canal_r == 2] = 0 # Suez
  
  dists = terra::gridDist(x = rast_noDNA_mask_pres,
                          target = 1,
                          scale = 1000)
  
  logdists = log10(dists+1)
  # terra::writeRaster(logdists, file = here::here("outputs", "marie", "p_macrolepis", "rast_p_macrolepis.asc"),
  #                    overwrite = T)
  # 
  # terra::plot((dists))

  dist_to_pres = raster::extract(dists,
                                 unique_DNA,
                                 df = T,
                                 cells = T)
 
  max_dist = max(dist_to_pres$values_length, na.rm = T)
  
  Species = vec_sp[i]
  
  df_dist = data.frame(Species, max_dist)
  cat(i, "\n")
  
  return(df_dist)
  
})

library(ggplot2) 
library(tidyterra)

ggplot() +
  geom_spatraster(data = (logdists), aes(fill = values_length)) +
  geom_point(data = unique_DNA,
             aes(x = longitude,
                 y = latitude),
             colour = "white") +
  viridis::scale_fill_viridis(option = "B") +
  theme(legend.position = "bottom") +
  theme(legend.key.width = unit(4, "cm"))

df_dists = data.table::rbindlist(all_sp_dists)
false_positives = c("Anampses chrysocephalus", #Labridae 
                    "Merluccius productus", #Consumption
                    "Sardinella longiceps", #Consumption
                    "Merluccius merluccius", #Consumption
                    "Melanogrammus aeglefinus", #Consumption
                    "Chelon ramada", # 12 reads
                    "Encrasicholina punctifer", #Consumption + 58 reads
                    "Sparisoma viride", #Scaridae
                    "Sparisoma aurofrenatum", #Scaridae
                    "Sparisoma chrysopterum", #Scaridae
                    "Chlorurus strongylocephalus", #Scaridae
                    "Enchelyopus cimbrius") #Error because NA

df_dists_filt = df_dists |>
  dplyr::filter(!Species %in% false_positives)

# save(df_dists_filt, file = here::here("outputs", "df_dists_filt.RData"))

median(df_dists_filt$max_dist)
hist(log10(df_dists_filt$max_dist+1))
hist((df_dists_filt$max_dist))

df_dists_filt$fact = "fact"

d = ggplot(df_dists_filt |> dplyr::filter(max_dist > 0)) +
      geom_violin(aes(x = fact, y = max_dist),
                  fill = "#FFB06B",
                  alpha = 0.6,
                  width = 0.5) +
      geom_jitter(aes(x = fact, y = max_dist),
                  size = 3, 
                  alpha = 0.8) +
      geom_boxplot(aes(x = fact, y = max_dist), 
                   fill = "#FFB06B",
                   alpha = 0.6,
                   width = 0.2) +
      scale_y_continuous(trans = "log10") +
      ylab("Maximum distance to the nearest GBIF/OBIS occurrence") +
      theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank())

ggsave(filename = here::here("outputs", "suppl_fig_y.jpg"),
       plot = d,
       dpi = 300,
       width = 1800,
       height = 2548,
       units = "px")


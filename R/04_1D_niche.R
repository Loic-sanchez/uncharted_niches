grid_spy = function(){
  
  ever_sampled = terra::rast(here::here("outputs", "ever_sampled.asc"))
  
  SPY_coords = unique(df_dna_meta %>% 
                        dplyr::select(code_spygen, longitude_start, latitude_start))
  
  SPY_points = terra::vect(SPY_coords, geom = c("longitude_start", "latitude_start"))
  
  SPY_extract = as.data.frame(terra::extract(ever_sampled,
                                             SPY_points, 
                                             cells = T, bind = T))
  
  save(SPY_extract, file = here::here("outputs", "SPY_extract.Rdata"))
  
}

compute_niche = function(){
  
  source(here::here("R", "01_get_occus.R"))

  load(here::here("outputs", "SPY_extract.Rdata"))
  r_focal = terra::rast(here::here("data", "treated_rasters", "focal_rasts.tif"))
  ever_sampled = terra::rast(here::here("outputs", "ever_sampled.asc"))
  names(ever_sampled) = "ever_sampled"
  names(r_focal)[6] = "Cumulative_impact"
 
  get_gbif_obis()
  load(here::here("outputs", "df_dna_meta.Rdata"))
  load(here::here("outputs", "df_dna_meta_conf.Rdata"))
  load(here::here("outputs", "completeness_projs.Rdata"))
  
  meds = c("Gombessa 6+", "CERCA", "DeepHeart", "eRef", "corse_OFB", "corse_OEC", 
           "ANGE", "Corse_OEC", "PISCIS", "med", "gombessa")
  
  df_dna_meta = df_dna_meta |> 
    tidyr::drop_na(longitude_start) |> 
    dplyr::mutate(Project = ifelse(project %in% meds, "med", project))
  
  df_dna_meta_conf = df_dna_meta |> 
    dplyr::left_join(completeness_projs,
                     by = c("species_name_corrected" = "Species", 
                            "Project")) |> 
    dplyr::filter(Completeness >= 0.95)
  
  safe_sharkies = c("Isurus oxyrinchus", 
                  "Lamna nasus", 
                  "Carcharhinus obscurus",
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
  
  df_dna_meta_conf = rbind(df_dna_meta_conf, df_dna_meta_isu)
  # save(df_dna_meta_conf, here::here("outputs", "df_dna_meta_conf.Rdata"))
  
  # Join obis/gbif, and get unique occurrences
  
  OBIS5 = OBIS4 |> 
    dplyr::select(species, latitude, longitude)
  GBIF1 = unique_gbif_occus |> 
    dplyr::select(species, latitude, longitude)
  
  all_no_DNA = rbind(OBIS5,
                     GBIF1)
  
  # Filter DNA, unique etc
  
  eDNA4b = df_dna_meta_conf |> 
    dplyr::rename(species = "species_name_corrected",
                  longitude = "longitude_start",
                  latitude = "latitude_start") |> 
    dplyr::select(species, latitude, longitude)
  
  eDNA5 = unique(eDNA4b |> 
                   dplyr::select(species, latitude, longitude) |>  
                   tidyr::drop_na(longitude)) |> 
    tidyr::drop_na(latitude)
  
  vec_sp = sort(unique(eDNA5$species))
  
  # Create an empty raster
  
  r = terra::rast(nrows = 1800, ncols = 3599)
  r = terra::project(r, r_focal)

  # Compute 1D niche increases 
  
  load(here::here("outputs", "SPY_extract.Rdata"))
  
  all_sp_incr = lapply(1:length(vec_sp), function(i){
    
    # Filter by species i
    
    eDNA6 = eDNA5 |> 
      dplyr::filter(species == vec_sp[i])
    
    all_no_DNA1 = all_no_DNA |> 
      dplyr::filter(species == vec_sp[i])
    
    unique_no_DNA = unique(all_no_DNA1 |> dplyr::select(longitude, latitude))
    unique_DNA = unique(eDNA6 |> dplyr::select(longitude, latitude))
    
    SPY_with_sp = df_dna_meta_conf$code_spygen[df_dna_meta_conf$species_name_corrected == vec_sp[i]]
    SPY_extract_sp = SPY_extract |> dplyr::filter(code_spygen %in% SPY_with_sp)
    
    if(nrow(unique_DNA) < 1) {return(NULL)}
    if(nrow(unique_no_DNA) < 5) {return(NULL)}
    
    # Create a raster of presence/absence 
    
    rast_noDNA = terra::rasterize(as.matrix(unique_no_DNA), 
                                  r, 
                                  fun = "length")
    # terra::plot(log(rast_noDNA))
    
    rast_DNA = terra::rasterize(as.matrix(unique_DNA),
                                r, 
                                fun = 'length')
    
    # terra::plot(rast_DNA)
    
    rast_noDNA[is.na(rast_noDNA[])] <- 0
    raster_noDNA_masked = terra::mask(rast_noDNA, r_focal$NPP_median, inverse = F)
    
    # Back to Sp
    
    dna_df = terra::as.points(rast_DNA)
    no_dna_df = terra::as.points(raster_noDNA_masked)
    
    # terra::plot(no_dna_df)
    
    # Extract PC for GBIF presences and DNA outside of these
    
    # Add presences to the environmental and fishing
    
    pres_env = c(r_focal$SST_median, 
                 r_focal$NPP_median,
                 r_focal$SBT_median, 
                 r_focal$SBDO_median,
                 r_focal$Cumulative_impact, 
                 raster_noDNA_masked, 
                 ever_sampled)
    
    # Extract PC for GBIF presences and DNA outside of these
    
    dna_env = terra::extract(pres_env,
                             dna_df,
                             cells = T, 
                             xy = T)
    
    no_dna_env = terra::extract(pres_env,
                                no_dna_df)
    
    dna_out_env = dna_env |> dplyr::filter(ever_sampled < 1)
    dna_in_env = dna_env |>  dplyr::filter(ever_sampled > 0)
    no_dna_env_pres = no_dna_env |> dplyr::filter(values_length > 0)
    
    if(nrow(no_dna_env_pres) < 5){return(NULL)}
    if(nrow(dna_env) < 1){return(NULL)}
    
    # 1D increases
    
    SST_new_range = range(c(no_dna_env_pres$SST, dna_env$SST), na.rm = T)
    SST_old_range = range(no_dna_env_pres$SST, na.rm = T)
    SST_Lower_incr = ((SST_old_range[1]-SST_new_range[1])/(SST_old_range[2]-SST_old_range[1]))
    SST_Upper_incr = ((SST_new_range[2]-SST_old_range[2])/(SST_old_range[2]-SST_old_range[1]))

    NPP_new_range = range(c(no_dna_env_pres$NPP, dna_env$NPP), na.rm = T)
    NPP_old_range = range(no_dna_env_pres$NPP, na.rm = T)
    NPP_Lower_incr = ((NPP_old_range[1]-NPP_new_range[1])/(NPP_old_range[2]-NPP_old_range[1]))
    NPP_Upper_incr = ((NPP_new_range[2]-NPP_old_range[2])/(NPP_old_range[2]-NPP_old_range[1]))
   
    SBT_new_range = range(c(no_dna_env_pres$SBT, dna_env$SBT), na.rm = T)
    SBT_old_range = range(no_dna_env_pres$SBT, na.rm = T)
    SBT_Lower_incr = ((SBT_old_range[1]-SBT_new_range[1])/(SBT_old_range[2]-SBT_old_range[1]))
    SBT_Upper_incr = ((SBT_new_range[2]-SBT_old_range[2])/(SBT_old_range[2]-SBT_old_range[1]))
    
    SBDO_new_range = range(c(no_dna_env_pres$SBDO, dna_env$SBDO), na.rm = T)
    SBDO_old_range = range(no_dna_env_pres$SBDO, na.rm = T)
    SBDO_Lower_incr = ((SBDO_old_range[1]-SBDO_new_range[1])/(SBDO_old_range[2]-SBDO_old_range[1]))
    SBDO_Upper_incr = ((SBDO_new_range[2]-SBDO_old_range[2])/(SBDO_old_range[2]-SBDO_old_range[1]))
    
    Impact_new_range = range(c(no_dna_env_pres$Cumulative_impact, dna_env$Cumulative_impact), na.rm = T)
    Impact_old_range = range(no_dna_env_pres$Cumulative_impact, na.rm = T)
    Impact_Upper_incr = ((Impact_new_range[2]-Impact_old_range[2])/(Impact_old_range[2]-Impact_old_range[1]))
    Impact_Lower_incr = ((Impact_old_range[1]-Impact_new_range[1])/(Impact_old_range[2]-Impact_old_range[1]))
    
    # 1D increases with DNA 
    
    DNA_SST_new_range = range(c(no_dna_env_pres$SST, dna_in_env$SST), na.rm = T)
    DNA_SST_old_range = range(no_dna_env_pres$SST, na.rm = T)
    DNA_SST_Lower_incr = ((DNA_SST_old_range[1]-DNA_SST_new_range[1])/(DNA_SST_old_range[2]-DNA_SST_old_range[1]))
    DNA_SST_Upper_incr = ((DNA_SST_new_range[2]-DNA_SST_old_range[2])/(DNA_SST_old_range[2]-DNA_SST_old_range[1]))
    
    DNA_NPP_new_range = range(c(no_dna_env_pres$NPP, dna_in_env$NPP), na.rm = T)
    DNA_NPP_old_range = range(no_dna_env_pres$NPP, na.rm = T)
    DNA_NPP_Lower_incr = ((DNA_NPP_old_range[1]-DNA_NPP_new_range[1])/(DNA_NPP_old_range[2]-DNA_NPP_old_range[1]))
    DNA_NPP_Upper_incr = ((DNA_NPP_new_range[2]-DNA_NPP_old_range[2])/(DNA_NPP_old_range[2]-DNA_NPP_old_range[1]))
    
    DNA_SBDO_new_range = range(c(no_dna_env_pres$SBDO, dna_in_env$SBDO), na.rm = T)
    DNA_SBDO_old_range = range(no_dna_env_pres$SBDO, na.rm = T)
    DNA_SBDO_Lower_incr = ((DNA_SBDO_old_range[1]-DNA_SBDO_new_range[1])/(DNA_SBDO_old_range[2]-DNA_SBDO_old_range[1]))
    DNA_SBDO_Upper_incr = ((DNA_SBDO_new_range[2]-DNA_SBDO_old_range[2])/(DNA_SBDO_old_range[2]-DNA_SBDO_old_range[1]))
    
    DNA_SBT_new_range = range(c(no_dna_env_pres$SBT, dna_in_env$SBT), na.rm = T)
    DNA_SBT_old_range = range(no_dna_env_pres$SBT, na.rm = T)
    DNA_SBT_Lower_incr = ((DNA_SBT_old_range[1]-DNA_SBT_new_range[1])/(DNA_SBT_old_range[2]-DNA_SBT_old_range[1]))
    DNA_SBT_Upper_incr = ((DNA_SBT_new_range[2]-DNA_SBT_old_range[2])/(DNA_SBT_old_range[2]-DNA_SBT_old_range[1]))
    
    DNA_Impact_new_range = range(c(no_dna_env_pres$Cumulative_impact, dna_in_env$Cumulative_impact), na.rm = T)
    DNA_Impact_old_range = range(no_dna_env_pres$Cumulative_impact, na.rm = T)
    DNA_Impact_Upper_incr = ((DNA_Impact_new_range[2]-DNA_Impact_old_range[2])/(DNA_Impact_old_range[2]-DNA_Impact_old_range[1]))
    DNA_Impact_Lower_incr = ((DNA_Impact_old_range[1]-DNA_Impact_new_range[1])/(Impact_old_range[2]-Impact_old_range[1]))
    
    # 1D increases with outside 
    
    out_SST_new_range = range(c(no_dna_env_pres$SST, dna_out_env$SST), na.rm = T)
    out_SST_old_range = range(no_dna_env_pres$SST, na.rm = T)
    out_SST_Lower_incr = ((out_SST_old_range[1]-out_SST_new_range[1])/(out_SST_old_range[2]-out_SST_old_range[1]))
    out_SST_Upper_incr = ((out_SST_new_range[2]-out_SST_old_range[2])/(out_SST_old_range[2]-out_SST_old_range[1]))
    
    out_NPP_new_range = range(c(no_dna_env_pres$NPP, dna_out_env$NPP), na.rm = T)
    out_NPP_old_range = range(no_dna_env_pres$NPP, na.rm = T)
    out_NPP_Lower_incr = ((out_NPP_old_range[1]-out_NPP_new_range[1])/(out_NPP_old_range[2]-out_NPP_old_range[1]))
    out_NPP_Upper_incr = ((out_NPP_new_range[2]-out_NPP_old_range[2])/(out_NPP_old_range[2]-out_NPP_old_range[1]))
    
    out_SBDO_new_range = range(c(no_dna_env_pres$SBDO, dna_out_env$SBDO), na.rm = T)
    out_SBDO_old_range = range(no_dna_env_pres$SBDO, na.rm = T)
    out_SBDO_Lower_incr = ((out_SBDO_old_range[1]-out_SBDO_new_range[1])/(out_SBDO_old_range[2]-out_SBDO_old_range[1]))
    out_SBDO_Upper_incr = ((out_SBDO_new_range[2]-out_SBDO_old_range[2])/(out_SBDO_old_range[2]-out_SBDO_old_range[1]))
    
    out_SBT_new_range = range(c(no_dna_env_pres$SBT, dna_out_env$SBT), na.rm = T)
    out_SBT_old_range = range(no_dna_env_pres$SBT, na.rm = T)
    out_SBT_Lower_incr = ((out_SBT_old_range[1]-out_SBT_new_range[1])/(out_SBT_old_range[2]-out_SBT_old_range[1]))
    out_SBT_Upper_incr = ((out_SBT_new_range[2]-out_SBT_old_range[2])/(out_SBT_old_range[2]-out_SBT_old_range[1]))
    
    out_Impact_new_range = range(c(no_dna_env_pres$Cumulative_impact, dna_out_env$Cumulative_impact), na.rm = T)
    out_Impact_old_range = range(no_dna_env_pres$Cumulative_impact, na.rm = T)    
    out_Impact_Upper_incr = ((out_Impact_new_range[2]-out_Impact_old_range[2])/(out_Impact_old_range[2]-out_Impact_old_range[1]))
    out_Impact_Lower_incr = ((out_Impact_old_range[1]-out_Impact_new_range[1])/(out_Impact_old_range[2]-out_Impact_old_range[1]))
    
    # Final item 
    
    increases = data.frame(species = vec_sp[i],
                           NPP_Upper_incr,
                           NPP_Lower_incr,
                           SST_Upper_incr,
                           SST_Lower_incr,
                           SBDO_Upper_incr,
                           SBDO_Lower_incr,
                           SBT_Upper_incr,
                           SBT_Lower_incr,
                           Impact_Upper_incr,
                           Impact_Lower_incr,
                           DNA_NPP_Lower_incr,
                           DNA_NPP_Upper_incr,
                           DNA_SST_Lower_incr,
                           DNA_SST_Upper_incr,
                           DNA_SBT_Lower_incr,
                           DNA_SBT_Upper_incr,
                           DNA_SBDO_Lower_incr,
                           DNA_SBDO_Upper_incr,
                           DNA_Impact_Lower_incr,
                           DNA_Impact_Upper_incr,
                           out_NPP_Lower_incr,
                           out_NPP_Upper_incr,
                           out_SST_Lower_incr,
                           out_SST_Upper_incr,
                           out_SBT_Lower_incr,
                           out_SBT_Upper_incr,
                           out_SBDO_Lower_incr,
                           out_SBDO_Upper_incr,
                           out_Impact_Lower_incr,
                           out_Impact_Upper_incr
                           )
    
    cat(i, "\n")
    return(increases)
    
  })
  
  all_incr = data.table::rbindlist(all_sp_incr)
  # save(all_incr, file = here::here("outputs", "all_incr.Rdata"))
  
}

vec_sp_useful = unique(all_incr_long$species) # In case you should re-run only with the species that display increases

# For confidence ind

false_positives = c("Anampses chrysocephalus", #Labridae 
                    "Merluccius productus", #Consumption
                    "Merluccius merluccius", #Consumption
                    "Melanogrammus aeglefinus", #Consumption
                    "Sardinella longiceps", #Consumption
                    "Chelon ramada", # 12 reads
                    "Encrasicholina punctifer", #Consumption + 58 reads
                    "Sparisoma viride", #Scaridae
                    "Sparisoma aurofrenatum", #Scaridae
                    "Sparisoma chrysopterum", #Scaridae
                    "Chlorurus strongylocephalus") #Scaridae

all_incr_long <- all_incr[, c(1:11)] |>
  tidyr::gather(key, value, -species) |>
  tidyr::separate(key, into = c("variable", "direction"), sep = "_") |> 
  dplyr::filter(value > 0) |>
  dplyr::mutate(percent_incr = (value*100)) |> 
  dplyr::filter(!species %in% false_positives)

library(ggplot2)
library(tidyverse)

all_incr_long$direction <- factor(all_incr_long$direction)

# Preprocess data
transformed_data <- all_incr_long %>%
  mutate(log_percent_incr = log10(percent_incr + 1),
         log_percent_incr = ifelse(direction == "Lower", -log_percent_incr, log_percent_incr),
         direction = ifelse(direction == "Lower", "Lower (log10)", "Upper (log10)"),
         natural_scale = ifelse(direction == "Lower", -percent_incr, percent_incr))

# Calculate the range for the y-axis
y_range <- range(transformed_data$natural_scale)

# Determine the highest power of 10 needed to display the data
max_power <- ceiling(log10(max(abs(y_range)) + 1))

# Define breaks and labels
upper_breaks <- c(0, 1, 2, 3)
lower_breaks <- c(-3, -2, -1, 0)
labels_upper <- c(0, 10^1, 10^2, 10^3)
labels_lower <- rev(labels_upper)

# Plot

g = ggplot(transformed_data) +
  geom_jitter(aes(x = variable, 
                  y = log_percent_incr, 
                  colour = direction),
              colour = "black",
              position = position_jitter(width = 0.2), 
              size = 3, 
              alpha = 0.8,
              show.legend = F) +
  geom_boxplot(data = transformed_data |> dplyr::filter(direction == "Upper (log10)"),
               aes(x = variable, 
                   y = log_percent_incr),
               fill = "#FF4242",
               alpha = 0.75, 
               width = 0.25,
               show.legend = FALSE) +
  geom_boxplot(data = transformed_data |> dplyr::filter(direction == "Lower (log10)"),
               aes(x = variable, 
                   y = log_percent_incr),
               fill = "#0C6B37",
               alpha = 0.75, 
               width = 0.25,
               show.legend = FALSE) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14)) +
  ylab("Niche range expansions (%)") +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(breaks = c(lower_breaks, upper_breaks),
                     labels = c(labels_lower, labels_upper),
                     limits = c(-max_power, max_power)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

ggsave(filename = here::here("outputs", "fig_3.jpg"),
       plot = g,
       dpi = 300,
       width = 1800,
       height = 2548,
       units = "px")


# DNA vs OUT

all_incr_long_type <- all_incr[, c(1, 12:31)] |>
  tidyr::gather(key, value, -species) |>
  tidyr::separate(key, into = c("type", "variable", "direction"), sep = "_") |> 
  dplyr::filter(value > 0) |> 
  dplyr::mutate(percent_incr = (value*100)) 

all_incr_long_DNA = all_incr_long_type |> 
  dplyr::filter(type == "DNA") |> 
  dplyr::rename(DNA = "value") # length 54 sp
all_incr_long_out = all_incr_long_type |> 
  dplyr::filter(type == "out") |> 
  dplyr::rename(out = "value") # length 243 sp

all_incr_full_join = dplyr::full_join(all_incr_long_DNA |> dplyr::select(-type), 
                                      all_incr_long_out |> dplyr::select(-type), 
                                      by = c("species", "variable", "direction"))
all_incr_full_join[is.na(all_incr_full_join)] = 0

traits = rfishbase::species(unique(all_incr$species))

all_incr_traits = all_incr_full_join |> 
  dplyr::left_join(traits |> 
                     dplyr::rename(species = "Species") |> 
                     dplyr::select(species, Length, Subfamily) |> 
                     dplyr::mutate(log_length = log10(Length)),
                   by = "species")

all_incr_trait = all_incr_full_join |> 
  dplyr::left_join(traits |> 
                     dplyr::rename(species = "Species") |> 
                     dplyr::select(species, Length, Subfamily) |> 
                     dplyr::mutate(log_length = log10(Length)),
                   by = "species")

ggplot(all_incr_traits, 
       aes(x = out, 
           y = DNA,
           colour = decile_rank)) +
  geom_point() +
  theme_minimal()

DNA = ggplot(all_incr_long_DNA |> dplyr::left_join(traits|> 
                                                     dplyr::rename(species = "Species") |> 
                                                     dplyr::select(species, Length) |> 
                                                     dplyr::mutate(log_length = log10(Length)),
                                                   by = "species")) +
  geom_histogram(aes(log_length), colour = "black", fill = "grey", binwidth = 0.2) +
  theme_minimal() +
  xlab("") +
  ylab("Counts (DNA)") +
  xlim(0,3)

out = ggplot(all_incr_long_out |> dplyr::left_join(traits|> 
                                                     dplyr::rename(species = "Species") |> 
                                                     dplyr::select(species, Length) |> 
                                                     dplyr::mutate(log_length = log10(Length)),
                                                   by = "species")) +
  geom_histogram(aes(log_length), colour = "black", fill = "grey", binwidth = 0.2) +
  theme_minimal() +
  xlim(0,3) +
  xlab("Log10(Length)") +
  ylab("Counts (out)")

gridExtra::grid.arrange(DNA, out)

all_incr_long_traits = all_incr_trait |>
  dplyr::select(species, variable, direction, DNA, out, log_length) |>
  tidyr::pivot_longer(cols = c("DNA", "out")) |>
  dplyr::filter(value > 0) |>
  tidyr::drop_na(log_length)
# 
# ggplot(all_incr_long_traits, aes(x = log_length, fill = name)) +
#   geom_histogram(position = "identity", 
#                  alpha = 0.75, 
#                  bins = 15,
#                  colour = "black") +
#   scale_fill_manual(values = c("orange", "darkblue")) +
#   theme_minimal() +
#   labs(x = "Value", y = "Frequency", title = "Histogram Comparison") +
#   theme(legend.title = element_blank())

supp_fig_x = ggplot(all_incr_long_traits, aes(x = name, y = log_length, fill = name)) +
   geom_violin(alpha = 0.8) +
   geom_jitter(aes(x = name, y = log_length), 
              colour = "black",
              position = position_jitter(width = 0.2), 
              size = 3, 
              alpha = 0.8,
              show.legend = F) +
  geom_boxplot(width = 0.1, 
               alpha = 0.7,
               outlier.shape = NA,
               show.legend = F) +
  scale_fill_manual(values = c("#FFE5B4", "#FF7A21"), labels = c("GBIF/OBIS + eDNA", 
                                                                 "eDNA only")) +
  scale_y_continuous(breaks = c(0, 1, 2),
                     labels = c(0, 10, 100),
                     limits = c(0, 2.7)) +
  theme_minimal() +
  labs(fill = "Cell sampling state") +
  ylab("Maximum length (cm)") +
  # xlab("Cell sampling state") +
  theme(legend.position = c(0.73, 0.15),
        legend.title = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        legend.background = element_rect(linewidth = 0.6,
                                         colour = "black"),
        legend.text = element_text(size = 12),
        # legend.key.spacing.y = unit(0.25, "lines"),
        legend.spacing.y = unit(1, "lines"))

ggsave(filename = here::here("outputs", "fig_4.jpg"),
       plot = supp_fig_x,
       dpi = 300,
       width = 1800,
       height = 2548,
       units = "px")

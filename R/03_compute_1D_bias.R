library(ggplot2)

uni_bias = function(){
  
  source(here::here("R", "02_prepare_rasters.R"))
  source(here::here("R", "00_join_obitools.R"))
  source(here::here("R", "01_get_occus.R"))
  
  # get all rasters 
  
  # focal_rasts = focal_rasters()
  r_focal = terra::rast(here::here("data", "treated_rasters", "focal_rasts.tif"))
  
  # get all occus
  
  get_gbif_obis()
  load(here::here("outputs", "df_dna_meta_conf.Rdata"))
  
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
  
  unique_DNA = unique(eDNA5 |> dplyr::select(longitude, latitude))
  
  # Rasterize the occurrences 
  
  # Rasterize 
  
  r = terra::rast(nrows = 1800, ncols = 3599)
  r = terra::project(r, r_focal)
  
  rast_noDNA = terra::rasterize(as.matrix(unique_no_DNA), 
                                r, 
                                fun = "length", background = 0)
  rast_noDNA[is.na(rast_noDNA[])] <- 0
  ever_sampled = log10(rast_noDNA+1)
  terra::writeRaster(ever_sampled,
                     filename = here::here("outputs",
                                           "ever_sampled.asc"),
                     datatype = 'FLT4S',
                     overwrite=TRUE)

  rast_DNA = terra::rasterize(as.matrix(unique(unique_DNA)),
                              r, 
                              fun = 'length')
  
  # count in out
   
  conventional_values <- terra::values(rast_noDNA)
  edna_values <- terra::values(rast_DNA)
  
  # Find cells with positive values in both rasters
  
  positive_in_both <- conventional_values > 0 & edna_values > 0
  
  # Count the number of cells with positive values in both rasters
  
  positive_in_both_count <- sum(positive_in_both, na.rm = TRUE)
  
  # Combine rasters with occurrence raster : 
  
  focal_rasts = r_focal
  pres_env = c(focal_rasts$SST_median, 
               focal_rasts$NPP_median, 
               focal_rasts$SBDO_median,
               focal_rasts$SBT_median,
               focal_rasts$global_cumul_impact_2013_all_layers, 
               rast_noDNA)
  names(pres_env) = c("SST", "NPP", "SBDO", "SBT", "Human Impact", "values_length")
  
  # Extract env. and fishing values for each point : 
  
  dna_df = terra::as.points(rast_DNA)
  no_dna_df = terra::as.points(rast_noDNA)
  
  dna_PC = terra::extract(pres_env,
                          dna_df, 
                          xy = T)
  no_dna_PC = terra::extract(pres_env,
                             no_dna_df)

  # Separate DNA points that fall where OBIS/GBIF already sampled or not
  
  dna_out_PC = dna_PC |> dplyr::filter(values_length < 1)
  dna_in_PC = dna_PC |> dplyr::filter(values_length > 0)
  
  # Get values for GBIF/OBIS
  
  no_dna_PC_pres = no_dna_PC |> dplyr::filter(values_length > 0) 
  
  dna_in_PC$type = "eDNA (in)"
  dna_out_PC$type = "eDNA (new)"
  no_dna_PC_pres$type = "GBIF/OBIS"
  
  # Join and melt into a single DF
  
  all_PC = rbind(dna_out_PC, no_dna_PC_pres, dna_in_PC)

  all_melt = reshape2::melt(all_PC, 
                            measure.vars = names(pres_env)[-6], 
                            na.rm = T) 
  
  calc_boxplot_stat <- function(x) {
    coef <- 1.5
    n <- sum(!is.na(x))
    
    # calculate quantiles
    
    stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
    names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
    iqr <- diff(stats[c(2, 4)])
    
    # set whiskers
    
    outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
    if (any(outliers)) {
      stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
    }
    
    return(stats)
    
  }
  
  library(ggplot2)
  
  cols = c("#357ba3", "#ff7801")
  all_melt$type <- forcats::fct_rev(all_melt$type)
  
  ggplot(all_melt |> dplyr::filter(!type == "eDNA (in)"), aes(x = type, y = value, fill = type)) +
    stat_summary(fun.data = calc_boxplot_stat, geom = "boxplot", show.legend = T) + 
    labs(fill = "Sampling type") +
    scale_fill_manual(values = cols) +
    xlab("Variable") +
    ylab("Value") +
    facet_wrap(~variable, scales = "free_y", nrow = 1) +
    theme_bw() + 
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(size = 14, 
                                    face = "bold"),
          axis.title = element_text(size = 14, 
                                    face = "bold"),
          axis.title.y = element_text(vjust = 0.25),
          legend.title = element_text(size = 12, face = "bold", margin = margin(b = 5)),
          legend.position = c(0.922, 0.85),
          legend.background = element_rect(linewidth = 0.2,
                                           colour = "black"),
          legend.text = element_text(size = 10),
          legend.spacing.y = unit(5, "mm"),
          panel.spacing = unit(1.25, "lines")) +
    guides(fill = guide_legend(keywidth = unit(1, "cm"),
                               keyheight = unit(1, "cm"), 
                               nrow = 2))
  
  ggsave(filename = here::here("outputs", "fig1b.png"), width = 3507, height = 1650, units = "px", dpi = 300)


all_melt_out = all_melt |> dplyr::filter(!type == "eDNA (in)")

t.test(all_melt_out$value[all_melt_out$type == "GBIF/OBIS" & all_melt_out$variable == "SST"],
       all_melt_out$value[all_melt_out$type == "eDNA (new)" & all_melt_out$variable == "SST"])
t.test(all_melt_out$value[all_melt_out$type == "GBIF/OBIS" & all_melt_out$variable == "NPP"],
       all_melt_out$value[all_melt_out$type == "eDNA (new)" & all_melt_out$variable == "NPP"])
t.test(all_melt_out$value[all_melt_out$type == "GBIF/OBIS" & all_melt_out$variable == "SBDO"],
       all_melt_out$value[all_melt_out$type == "eDNA (new)" & all_melt_out$variable == "SBDO"])
t.test(all_melt_out$value[all_melt_out$type == "GBIF/OBIS" & all_melt_out$variable == "SBT"],
       all_melt_out$value[all_melt_out$type == "eDNA (new)" & all_melt_out$variable == "SBT"])
t.test(all_melt_out$value[all_melt_out$type == "GBIF/OBIS" & all_melt_out$variable == "Human Impact"],
       all_melt_out$value[all_melt_out$type == "eDNA (new)" & all_melt_out$variable == "Human Impact"])

}

values_nodna = terra::values(ever_sampled)

rast_DNA[is.na(rast_DNA[])] <- 0
values_dna = terra::values(rast_DNA)

values_nodna = values_nodna[values_nodna > 0]
values_dna = values_dna[values_dna > 0]


dna_nodna = c(rep("GBIF/OBIS", 223956), rep("DNA", 275))
df = data.frame(Sampling_effort = c(values_nodna, values_dna),
                Type = dna_nodna)

library(ggplot2)
ggplot(df) +
  geom_boxplot(aes(x = Type, y = Sampling_effort)) +
  scale_y_continuous(transform = "log10") +
  theme_minimal()

# Median all rasts

median_rasters = function(){

  NPP_files = list.files(here::here("data",
                                  "monthly0.1", 
                                  "NPP"), 
                       pattern = "*.asc",
                       all.files = TRUE, 
                       full.names = TRUE)
  
  SST_files = list.files(here::here("data",
                                    "monthly0.1", 
                                    "SST"), 
                         pattern = "*.asc",
                         all.files = TRUE, 
                         full.names = TRUE)
  
  SSS_files = list.files(here::here("data",
                                    "monthly0.1", 
                                    "SSS"), 
                         pattern = "*.asc",
                         all.files = TRUE, 
                         full.names = TRUE)
  
  SBT_files = list.files(here::here("data",
                                    "monthly0.1", 
                                    "SBT"), 
                         pattern = "*.asc",
                         all.files = TRUE, 
                         full.names = TRUE)
  
  SBDO_files = list.files(here::here("data",
                                    "monthly0.1", 
                                    "SBDO"), 
                         pattern = "*.asc",
                         all.files = TRUE, 
                         full.names = TRUE)

  NPP_stack = terra::rast(NPP_files)
  NPP_median = terra::app(NPP_stack, median)
  SST_stack = terra::rast(SST_files)
  SST_median = terra::app(SST_stack, median)
  SSS_stack = terra::rast(SSS_files)
  SSS_median = terra::app(SSS_stack, median)
  SBT_stack = terra::rast(SBT_files)
  SBT_median = terra::app(SBT_stack, median)
  SBDO_stack = terra::rast(SBDO_files)
  SBDO_median = terra::app(SBDO_stack, median)

  terra::setGDALconfig("GDAL_PAM_ENABLED", "FALSE")
  
  terra::writeRaster(NPP_median, filename = here::here("data", "0.1summarized", "NPP_median.asc"), datatype='INT4S', overwrite=TRUE)
  terra::writeRaster(SST_median, filename = here::here("data", "0.1summarized", "SST_median.asc"), datatype='INT4S', overwrite=TRUE)
  terra::writeRaster(SSS_median, filename = here::here("data", "0.1summarized", "SSS_median.asc"), datatype='INT4S', overwrite=TRUE)
  terra::writeRaster(SBT_median, filename = here::here("data", "0.1summarized", "SBT_median.asc"), datatype='INT4S', overwrite=TRUE)
  terra::writeRaster(SBDO_median, filename = here::here("data", "0.1summarized", "SBDO_median.asc"), datatype='INT4S', overwrite=TRUE)
  
}

focal_rasters = function(){
  
  rastrz = list.files(here::here("data", "0.1summarized"), 
                      pattern = "*.asc",
                      all.files = TRUE, 
                      full.names = TRUE)
  
  impact = terra::rast(here::here("data", "0.1summarized", "global_cumul_impact_2013_all_layers.tif"))
  rasts = terra::rast(rastrz)
  impact_proj = terra::project(impact, rasts)
  all_rasts = c(rasts, impact_proj)
  r_focal = terra::focal(all_rasts, w = matrix(1/9, nc = 3, nr = 3), fun = mean, na.rm = F, na.policy = "only")
  terra::writeRaster(r_focal, 
                     file = here::here("data", "treated_rasters", "focal_rasts.tif"), 
                     overwrite = T, 
                     filetype = "GTiff")

  return(r_focal)
  
}

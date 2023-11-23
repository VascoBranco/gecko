#####GECKO - Geographical Ecology and Conservation Knowledge Online
#####Version 1.0.0 (2022-05-28)
#####By Vasco Branco, Pedro Cardoso, Luís Correia
#####Maintainer: vasco.branco@helsinki.fi
#####Changed from v0.1.2:
#####Full release. Added new SPECTRE functionalities.
#####Renamed gecko.examples to gecko.data

################################################################################
##################################MAIN FUNCTIONS################################
################################################################################

#' Make a raster layer SPECTRE compatible
#' @description Transform a given raster object to the resolution, datum, 
#' projection and extent used in SPECTRE.
#' @param layers SpatRaster. A raster object that you would like to be SPECTRE compatible.
#' @param continuous logical. Whether the data present in \code{layers} is continuous. 
#' If \code{TRUE} bilinear interpolation will be used in the case of resampling and reprojection. 
#' if \code{FALSE} nearest neighbour will be used instead. 
#' See \code{\link[terra:resample]{terra::resample()}} for more information on interpolation methods.
#' @param filepath character. Optional file path to where the final raster layer 
#' should be saved, in the format "folder/file.tif". If \code{filepath} is \code{NULL} 
#' your layer will be saved to your current working directory.
#' @return SpatRaster.
#' @examples \dontrun{
#' # For the sake of demonstration we will transform our raster layer "range".
#' distribution = gecko.data("range")
#' standard_dist = spectrify(distribution)
#' terra::plot(standard_dist)
#' }
#' @export
spectrify = function(layers, continuous = TRUE, filepath = NULL){
  # Check input validity -------------------------------------------------------
  if (!is.null(filepath)){
    if (!file.exists(dirname(filepath))){
      warning("Directory listed in 'filepath' does not exist.")
      return(NULL)
    }
  }
  if (!(is(continuous, "logical") && length(continuous) == 1)){
    warning("'continuous' must be a single logical.")
    return(NULL)
  }
  if(length(dir(paste0(gecko.getDir(), "spectretemplate"), full.names = T)) == 0){
    warning("One or more methods you have selected require downloading large data.
    This is a large dataset that is prone to fail by timeout if downloaded
    through R. You can do this process now or later by running gecko.setDir() 
    and placing the files needed inside it after download them through your browser.
    See detailed instructions for this in the documention of gecko.template().
    Would you like to continue?\n")
    answer <- readline("Proceed [y/n]?")
    if (tolower(answer) == "y"){
      if (is.null(gecko.getDir())){
        gecko.setDir()
      }
      spectre.template() # download the data
      print(paste0("Data was saved to: '", gecko.getDir(), "spectretemplate'!" ))
    } else {
      return(NULL)
    }
  } 
  
  # paste0(gecko.getDir(), "spectretemplate")
  template = terra::rast(paste0(gecko.getDir(), "spectretemplate/spectre.template.tif"))

  if(!file.exists(paste0(gecko.getDir(), "temp_data"))){
    dir.create(paste0(gecko.getDir(), "temp_data"))
  }

  if (continuous){
    res_type = "bilinear"
  } else {
    res_type = "near"
  }
  
  file_args = list(datatype = "FLT4S", filetype = "GTiff", gdal = c("COMPRESS=LZW"),
                   overwrite = TRUE, NAflag = -3.4e+38)
  
  proj_args = list(x = layers, y = template, method = res_type,
                   filename = paste0(gecko.getDir(), "temp_data/temp_pa_raster.tif"),
                   datatype = "FLT4S", filetype = "GTiff", overwrite=TRUE, NAflag = -3.4e+38)
  
  # Cropping ----------------------------------------------------------------
  layer_ext = terra::ext(layers)
  if (!(layer_ext[1] > -180 && layer_ext[2] < 180 && layer_ext[3] > -60 && layer_ext[4] < 90)){
    layers = terra::crop(layers, c(-180, 180, -60, 90))
  }

  # Reprojection / Resampling --------------------------------------------------
  warning("Checking dimensions.")
  if (terra::crs(layers, proj = TRUE) != "+proj=longlat +datum=WGS84 +no_defs"){
    layers = do.call(terra::project, proj_args)
  } else if (!terra::compareGeom(layers, template, crs = FALSE)){
    layers = do.call(terra::resample, proj_args)
  }
  
  # Check missing values  ------------------------------------------------------
  warning(paste("Checking for patchable missing cells in:", names(layers)))
  # Following operations saving to disk
  PA_raster <- do.call(terra::logic, c(file_args, list(
    x = layers, oper = "is.na",
    filename = paste0(gecko.getDir(), "temp_data/temp_sum_1.tif")
  )))

  PA_template <- do.call(terra::logic, c(file_args, list(
    x = template, oper = "is.na",
    filename = paste0(gecko.getDir(), "temp_data/temp_sum_2.tif")
  )))

  sum_raster <- do.call(
    terra::app,
    c(list(
      x = c(PA_template, PA_raster), fun = sum, overwrite = TRUE,
      filename = paste0(gecko.getDir(), "temp_data/temp_sum_3.tif"),
      wopt = file_args
    ))
  )

  sum_raster <- do.call(
    terra::app,
    c(list(
      x = template, fun = function(i) i == 1, overwrite = TRUE,
      filename = paste0(gecko.getDir(), "temp_data/temp_sum_3.tif"),
      wopt = file_args
    ))
  )
  
  # now mask with template, points with info that the template does not have will be discarded.
  sum_raster_2 <- do.call(
    terra::mask,
    c(file_args, list(
      x = sum_raster, mask = template,
      filename = paste0(gecko.getDir(), "temp_data/temp_sum_mask.tif")
    ))
  )
  start_missing = terra::freq(sum_raster_2, value = TRUE)[1,3]
  print(paste0("There are ", start_missing, " points found in the template missing from your layer."))
  # Mask with template ------------------------------------------------------
  layers <- do.call(
    terra::mask,
    c(file_args, list(
      x = layers, mask = template,
      filename = paste0(gecko.getDir(), "temp_data/temp_sum_5.tif")
      ))
  )
  # Saving ------------------------------------------------------------------
  if (is.null(filepath)){
    prev_files = grep("spectrified_layer", dir(getwd()))
    filepath = paste0(getwd(), "/spectrified_layer_", length(prev_files) + 1, ".tif")
  } 
  layers = do.call(terra::writeRaster, c(file_args, list(x = layers, filename = filepath)))
  file.remove(dir(paste0(gecko.getDir(), "temp_data"), full.names = TRUE))
  return(layers)
}

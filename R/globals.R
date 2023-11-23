#####GECKO - Geographical Ecology and Conservation Knowledge Online
#####Version 1.0.0 (2022-05-28)
#####By Vasco Branco, Pedro Cardoso, Luís Correia
#####Maintainer: vasco.branco@helsinki.fi
#####Changed from v0.1.2:
#####Full release. Added new SPECTRE functionalities.
#####Renamed gecko.examples to gecko.data

############################################################################
################################## PACKAGES - fn_spatial ###################
############################################################################
#' @importFrom terra global rast minmax ncell as.points as.data.frame writeRaster hist freq varnames ext merge layerCor values spatSample predict mask trim terrain rasterize crds classify distance vect plot mean crs extract logic
#' @importFrom grDevices dev.copy dev.off chull pdf
#' @importFrom sp SpatialPolygons Polygons Polygon
#' @importFrom graphics par text lines points title
#' @importFrom stats as.dist dist prcomp
#' @importFrom utils data download.file txtProgressBar setTxtProgressBar unzip
#' @importFrom geosphere areaPolygon
#' @importFrom methods is
#' @importFrom red map.draw
NULL
#> NULL 

############################################################################
################################## PACKAGES - fn_access ####################
############################################################################
#' @importFrom stats quantile
NULL
#> NULL 

############################################################################
################################## PACKAGES - outlier_fn ###################
############################################################################
#' @importFrom biomod2 BIOMOD_FormatingData predict
#' @importFrom kernlab ksvm
NULL
#> NULL 

############################################################################
##################################DATASETS##################################
############################################################################
#' Example data packaged with gecko
#' @description Load data included in the package. This includes \strong{records},
#' a matrix of longitude and latitude (two columns) occurrence records for
#' Hogna maderiana (Walckenaer, 1837); \strong{range}, a SpatRaster object, as
#' defined by package terra, of the geographic range of Hogna maderiana
#' (Walckenaer, 1837); \strong{layers}, a SpatRaster object with layers 
#' representing the average annual temperature, total annual precipitation,
#' altitude and landcover for Madeira Island
#' (Fick & Hijmans 2017, Tuanmu & Jetz 2014); \strong{threat}, a layer of mean 
#' fire occurence in Madeira between 2006 and 2016; and \strong{worldborders} is a
#' simplified version of the vector of world country borders created by
#'  \href{https://github.com/victorcazalis/RedList_countries}{Victor Cazalis}.
#' @param data character. String of one of the data names mentioned in the description, e.g.: \code{"gecko.records"}.
#' If \code{NULL}, the example files will be listed.
#' @examples
#' \dontrun{
#' gecko.data()
#' gecko.data("range")
#' }
#' @source This function is inspired by \code{\link[palmerpenguins:path_to_file]{palmerpanguins::path_to_file()}}
#' which in turn is based on \code{\link[readxl:readxl_example]{readxl::readxl_example()}}.
#' @export
gecko.data <- function(data = NULL) {
  if (is.null(data)) {
    print(
      c(
        "records", "range", "layers", "worldborders" 
      )
    )
    return(NULL)
  } else {
    if (data == "records"){
      path = system.file(paste0("extdata/gecko.records.csv"), package = "gecko")
      out = utils::read.csv(path)
    } else if (data == "range") {
      path = system.file(paste0("extdata/gecko.range.tif"), package = "gecko")
      out = terra::rast(x = path)
    } else if (data == "layers") {
      path = system.file(paste0("extdata/gecko.layers.", c(1:4), ".tif"), package = "gecko")
      out = terra::rast(x = path)
    } else if (data == "threat") {  
      path = system.file(paste0("extdata/gecko.threat.tif"), package = "gecko")
      out = terra::rast(x = path)
    } else if (data == "worldborders") {
      path = system.file(paste0("extdata/worldborders"), package = "gecko")
      out = terra::vect(x = path)
    } else {
      warning("Invalid data name. Run gecko.data() for a full list of options.")
      return(NULL)
    }
    return(out)
  }
  
}

############################################################################
###################################SETUP####################################
############################################################################

#' Setup GIS directory.
#' @description Setup directory where GIS files are stored.
#' @param gisPath Path to the directory where the gis files are stored.
#' @details Writes a txt file in the red directory allowing the package to always access the world GIS files directory.
#' @export
gecko.setDir <- function(gisPath = NULL){
  if (is.null(gisPath)) {
    # why not save to find.package("gecko")?
    gisPath = paste0(find.package("gecko"), "/downloaded_data")
    dir.create(gisPath)
  } else {
    gisPath <- readline("Input directory for storing world gis layers:")
  }
  gisPath <- paste0(gisPath, "/")
  globalsFile <- paste0(find.package("gecko"), "/globals.txt")
  dput(gisPath, globalsFile)
}

#' Read GIS directory.
#' @description Read directory where GIS files are stored.
#' @details Reads a txt file pointing to where the world GIS files are stored.
#' @export
gecko.getDir <- function(){
  gecko_file <- paste(find.package("gecko"), "/globals.txt", sep = "")
  if (file.exists(gecko_file)) { # if there is already a file read from it
    dir <- dget(gecko_file)
  } else {
    warning(paste(gecko_file, "not found, please run gecko.setDir()"))
    return(NULL)
  }
  return(dir)
}

#' Download worldclim files.
#' @description Download the latest version of worldclim to your gecko work directory. 
#' If you have not yet setup a work directory, it will be be setup as if running 
#' \code{\link[gecko:gecko.setDir]{gecko::gecko.setDir()}} with \code{gisPath = NULL}.
#' This is a large dataset that is prone to fail by timeout if downloaded 
#' through R. Instead of using this function you can run gecko.setDir() (if you 
#' haven't yet) and download the files at 
#' https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_bio.zip or 
#' https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_bio.zip. 
#' Unzip their contents correspondingly to the folders "./worldclim/1 km" or 
#' "./worldclim/10 km" inside the folder returned by gecko.getDir().
#' @param res character. Specifies the resolution of environmental data used.
#' @details Reads a txt file pointing to where the world GIS files are stored.
#' @examples
#' \dontrun{
#' gecko.worldclim("10 km")
#' }
#' @export
gecko.worldclim = function(res){
  worldclim_refs = list(url_1km = "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_bio.zip",
                        url_10km = "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_bio.zip",
                        zip_name = "bioclim2.zip",
                        layer_prefix = "wc2.1_30s_bio_",
                        altitude_url = "http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/alt_30s_bil.zip",
                        altitude_zip = "alt_30s_bil.zip"
  )
  
  if (res == "1 km") {
    url = worldclim_refs$url_1km
  } else if (res == "10 km") {
    url = worldclim_refs$url_10km
  } else {
    warning("Invalid resolution.")
    return(NULL)
  }
  
  if (is.null(gecko.getDir())){
    gecko.setDir()
  }
  
  cat(paste0("WorldClim data will start downloading.",
                 " Due to its size it is prone to failure.",
                 " You can try to circunvent this by downloading the zip file",
                 " directly from: ", worldclim_refs$url,
                 " and then placing it in: ", gecko.getDir(), "worldclim \n"))
  
  ##download and process bioclim
  status = download.file(url, paste0(gecko.getDir(), worldclim_refs$zip_name))
  
  if (status != 0){
    warning("Download failed.")
    return(NULL)
  }
  
  if (!file.exists(paste0(gecko.getDir(), "worldclim"))){
    dir.create(paste0(gecko.getDir(), "worldclim"))
  }
  
  if (!file.exists(paste0(gecko.getDir(), "worldclim/", res))){
    dir.create(paste0(gecko.getDir(), "worldclim/", res))
  }
  
  
  utils::unzip(zipfile = paste0(gecko.getDir(), worldclim_refs$zip_name),
               exdir = paste0(gecko.getDir(), "worldclim/", res) )
  
  file.remove(paste0(gecko.getDir(), worldclim_refs$zip_name))
}

#' Download the SPECTRE template.
#' @description Download the raster template for SPECTRE layers to your gecko work directory. 
#' If you have not yet setup a work directory, it will be be setup as if running 
#' \code{\link[gecko:gecko.setDir]{gecko::gecko.setDir()}} with \code{gisPath = NULL}.
#' This is a large dataset that is prone to fail by timeout if downloaded 
#' through R. Instead of using this function you can run gecko.setDir() (if you 
#' haven't yet) and download the file at 
#' https://github.com/VascoBranco/spectre.content/raw/main/spectre.template.zip.
#' Unzip its contents to a folder "./spectretemplate" inside the folder returned by gecko.getDir().
#' @details Reads a txt file pointing to where the world GIS files are stored.
#' @examples
#' \dontrun{
#' spectre.template()
#' }
#' @export
spectre.template = function(){
  template_refs <- list(
    url = "https://github.com/VascoBranco/spectre.content/raw/main/spectre.template.zip",
    zip_name = "spectre.template.zip"
  )

  if (is.null(gecko.getDir())) {
    gecko.setDir()
  }

  if (!file.exists(paste0(gecko.getDir(), "spectretemplate"))) {
    dir.create(paste0(gecko.getDir(), "spectretemplate"))
  }

  options(timeout = 120)
  cat(paste0(
    "The SPECTRE template data will start downloading.",
    " Due to its size it is prone to failure, depending several factors such as ",
    "internet speed.",
    " You can try to circunvent this by downloading the zip file",
    " directly from: ", template_refs$url,
    " and then placing it in: ", gecko.getDir(), "spectretemplate \n"
  ))
  download.file(template_refs$url, paste0(gecko.getDir(), "spectretemplate/", template_refs$zip_name))
  options(timeout = 60)

  utils::unzip(
    zipfile = paste0(gecko.getDir(), "spectretemplate/", template_refs$zip_name),
    exdir = paste0(gecko.getDir(), "spectretemplate")
  )
  file.remove(paste0(gecko.getDir(), "spectretemplate/", template_refs$zip_name))
}

###############################################################################
###################################DEV ONLY####################################
###############################################################################

#' Metadata for available SPECTRE layers
#' @description For dev use. Load the SPECTRE metadata included in the package.
#' @noRd
spectre.metadata <- function() {
  path = system.file(paste0("extdata/spectre.metadata.csv"), package = "gecko")
  return(utils::read.csv(path))
}


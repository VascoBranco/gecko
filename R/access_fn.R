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

#' Get in text citations for SPECTRE layers
#' @description Generate in-text citations for a selection of SPECTRE layers. 
#' @param index numeric. A vector of integers specifying the layers. Refer to the Details section.
#' @details
#' The current layers in SPECTRE are:
#' \enumerate{
#'  \item \strong{MINING_AREA}. Mining density based on the number of known mining properties (pre-operational, operational, and closed) in a 50-cell radius (1x1 km cells).
#'  \item \strong{HAZARD_POTENTIAL}. Number of significant hazards (earthquakes, volcanoes, landslides, floods, drought, cyclones) potentially affecting cells based on hazard frequency data.
#'  \item \strong{HUMAN_DENSITY} Continuous metric of population density.
#'  \item \strong{BUILT_AREA} Percentage metric indicating the built-up presence.
#'  \item \strong{ROAD_DENSITY}. Continuous metric of road density.
#'  \item \strong{FOOTPRINT_PERC}. Percentage metric indicating anthropogenic impacts on the environment.
#'  \item \strong{IMPACT_AREA}. Classification of land into very low impact areas (1), low impact areas (2) and non-low impact areas (3).
#'  \item \strong{MODIF_AREA}. Continuous 0-1 metric that reflects the proportion of a landscape that has been modified.
#'  \item \strong{HUMAN_BIOMES}. Classification of land cover into different anthropogenic biomes of differing pressure such as dense settlements, villages and cropland.
#'  \item \strong{FIRE_OCCUR}. Continuous metric of mean fire occurrence during the years of 2006 and 2016.
#'  \item \strong{CROP_PERC_UNI}. Percentage metric indicating the proportion of cropland in each cell.
#'  \item \strong{CROP_PERC_IIASA}. Percentage metric indicating the proportion of cropland in each cell.
#'  \item \strong{LIVESTOCK_MASS}. Estimated total amount of livestock wet biomass based on global livestock head counts.
#'  \item \strong{FOREST_LOSS_PERC}. Continuous -100 to 100 metric of forest tree cover loss between 2007 and 2017.
#'  \item \strong{FOREST_TREND}. Classification metric of 0 (no loss) or a discrete value from 1 to 17, representing loss (a stand-replacement disturbance or change from a forest to non-forest state) detected primarily in the year 2001-2019, respectively.
#'  \item \strong{NPPCARBON_GRAM}. Quantity of carbon needed to derive food and fiber products (HANPP).
#'  \item \strong{NPPCARBON_PERC}. HANNP as a percentage of local Net Primary Productivity.
#'  \item \strong{LIGHT_MCDM2}. Continuous simulated zenith radiance data.
#'  \item \strong{FERTILIZER_LGHA}. Continuous metric of kilograms of fertilizer used per hectare.
#'  \item \strong{TEMP_TRENDS}. Continuous metric of temperature trends, based on the linear regression coefficients of mean monthly temperature for the years of 1950 to 2019.
#'  \item \strong{TEMP_SIGNIF}. Continuous metric of temperature trend significance, the temperature trends divided by its standard error.
#'  \item \strong{CLIM_EXTREME}. Continuous metric calculated as whatever is the largest of the absolute of the trend coefficients of the months with the lowest or highest mean temperatures.
#'  \item \strong{CLIM_VELOCITY}. Continuous metric of the velocity of climate change, the ratio between TEMP_TRENDS and a local spatial gradient in mean temperature calculated as the slope of a plane fitted to the values of a 3x3 cell neighbourhood centered on each pixel.
#'  \item \strong{ARIDITY_TREND}. Continuous metric of aridity trends, based on the linear regression coefficients of aridity for the years of 1990 to 2019, i.e: MPET/(MPRE+1).
#' }
#' @return list. Contains two elements, both characters: the first a single 
#' character containing the in-text citations, the second a character of 
#' length \code{x} with the bibliographic citations.
#' @examples
#' sources = c(2,3)
#' out = spectre.citations(sources)
#' @export
spectre.citations <- function(index) {
  # Here we get our columns inside unique() to filter out repeated sources.
  metadata <- spectre.metadata()
  work_frame <- unique(metadata[index, c(8, 10, 7, 9, 11)])

  for (i in row.names(work_frame)) {
    year <- work_frame[i, 1]
    if (isFALSE(substr(year, nchar(year), nchar(year)) %in% letters)) {
      # this looks really weird and I'd like a better solution. basically, this
      # checks if there are any rows with the exact same 2 cells. it uses matrix()
      # because if you instead use as.matrix here it changes the datatypes and
      # screws up the check.
      check <- work_frame[, 1:2] == matrix(work_frame[i, 1:2],
        nrow = 1,
        ncol = 2, byrow = TRUE
      )

      # If this is true there are repetitions and we need to intervene.
      if (sum(apply(check, FUN = all, 1)) > 1) {
        # We declare two vars, one for
        to_change <- work_frame[rowSums(check) == 2, 1]
        new_years <- letters[1:sum(rowSums(check) == 2)]

        year_n_letter <- vapply(1:length(to_change), function(x) {
          paste(work_frame[rowSums(check) == 2, 1][index], # our years
            letters[1:sum(rowSums(check) == 2)][index], # corresponding letters
            sep = ""
          )
        }, "string")

        work_frame[rowSums(check) == 2, 1] <- year_n_letter
      }
    }
  }
  work_frame <- work_frame[order(work_frame[, 1], decreasing = TRUE), ]

  correction <- work_frame[order(work_frame[, 1], decreasing = FALSE), 1]

  work_frame[, 1][is.na(as.integer(work_frame[, 1]))] <-
    correction[is.na(as.integer(correction))]

  # At this point we have the correct order.

  # Intext citations
  intext <- sapply(
    c(1:nrow(work_frame)),
    function(x) {
      paste(work_frame[x, c(2, 1)], collapse = " ")
    }
  )
  intext <- paste("(", paste(as.vector(intext), collapse = ", "), ")", sep = "")

  # Bibliographic citations
  bibs <- sapply(
    c(1:nrow(work_frame)),
    function(x) {
      paste(work_frame[x, c(3)], " (", work_frame[x, c(1)], ") ",
        work_frame[x, c(4)],
        sep = "", collapse = " "
      )
    }
  )

  return(list(intext, bibs))
}


#' Normalize raster.
#' @description Normalize a raster file according to one three methods, 'standard', 'range' or 'rank'.
#' @param layer SpatRaster. Object with a single layer as defined by package terra.
#' @param method character. Specifying \code{'standard'}, \code{'range'} or \code{'rank'}. 
#' @param filepath character. Optional, specifies a path to the output file.
#' @details  The three options, "standard" standardizes data to a mean = 0 and sd = 1,
#' "range" standardizes to a range of 0 to 1, and "rank" similarly standardizes to
#' a range of 0 to 1 but does so after ranking all points.
#' @return A raster layer.
#' @examples
#' \dontrun{
#' region = gecko.data("layers")[[1]]
#' ranked_region = normalize(region, method = "rank")
#' }
#' @export
normalize = function(layer, method = "standard", filepath = NULL){
  if ((is(filepath, "character") && length(filepath) != 1)) {
    warning("Filepath must be a single character string.")
    return(invisible())
  }

  if (!is(layer, "SpatRaster")) {
    warning("layer must be a SpatRaster object.")
    return(invisible())
  } else {
    if (dim(layer)[3] != 1) {
      warning("SpatRaster object must have a single layer.")
      return(invisible())
    }
  }

  # mean 0 and sd 1
  if (method == "standard") {
    output_raster <- (layer - as.numeric(terra::global(layer, "mean", na.rm = TRUE))) /
      as.numeric(terra::global(layer, "sd", na.rm = TRUE))

    # range 0-1
  } else if (method == "range") {
    lminmax <- terra::minmax(layer)
    output_raster <- (layer - lminmax[1]) / (lminmax[2] - lminmax[1])

    # ranking of values
  } else if (method == "rank") {
    output_raster <- layer
    output_raster[1:terra::ncell(output_raster)] <- rank(terra::as.data.frame(layer, na.rm = FALSE)[, 1], ties.method = "average", na.last = "keep")
    output_raster <- normalize(output_raster, method = "range", filepath = filepath)
  } else {
    warning("Method not recognized, must be one of 'standard', 'range' or 'rank'.")
    return(invisible())
  }
  
  if (is.null(filepath)){
    filepath = paste0(tempdir(), "\\", "norm_rast", ".tif")
  }

  terra::writeRaster(output_raster,
                     filename = filepath,
                     datatype = "FLT4S", filetype = "GTiff",
                     gdal = c("COMPRESS=LZW"), overwrite = TRUE,
                     NAflag = -3.4e+38)
  
  return(output_raster)
}

#' Get a short summary of a given raster segment.
#' @description Return a set of descriptive statistics of the given layer,
#' either a specific one (minimum, q1, median, q3, maximum,
#' median absolute deviation (mad), mean, standard deviation (sd)) or all of them.
#' @param layer SpatRaster. Raster object, as defined by package terra, with a single layer. 
#' @param plot logical. If TRUE, a histogram of raster values is drawn.
#' @return data.frame. If plot is TRUE, also outputs a histogram of the layer.
#' @examples
#' region = gecko.data("layers")
#' stats(region[[1]])
#' @export
stats = function(layer, plot = FALSE){
  if (plot){
    # if(layer@file@name == ""){
    #   layer_name = "temp"
    # } else {
    #   layer_name = layer@file@name
    # }
    # raster::hist(layer, main = paste("Histogram of", layer_name), xlab = "Value bins")
    terra::hist(layer, main = paste("Histogram of", terra::varnames(layer)), xlab = "Value bins")
  }
  
  nas = terra::freq(layer, value = NA)[,3]
  quant = stats::quantile(as.vector(layer), na.rm = TRUE)
  mad = stats::quantile(
    as.vector(layer - as.numeric(terra::global(layer, "mean", na.rm = TRUE))), na.rm = TRUE)[3]
  out = data.frame(stat = c("Max", "Min", "Q1", "Median", 
                            "Q3", "MAD", "Mean", "SD", "NA's"))
  
  lminmax = terra::minmax(layer)
  
  vals = sapply(c(lminmax[2], lminmax[1], quant[2], quant[3], quant[4],
           mad,
           as.numeric(terra::global(layer, "mean", na.rm = TRUE)),
           as.numeric(terra::global(layer, "sd", na.rm = TRUE)),
           nas), 
         round, digits = 5)
  out["value"] = vals
  return(out)
}

#' Get SPECTRE raster segments.
#' @description Downloads SPECTRE segments according to a bounding box selection.
#' @param index numeric. A vector of integers specifying the layers. Refer to the list.
#' @param ext numeric or SpatExtent. A vector of \code{xmin}, \code{xmax}, \code{ymin}, \code{ymax} or a \code{terra} 
#' spatial extent object (See \code{\link[terra:ext]{terra::ext()}}). If no 
#' input is given, an extent of \code{xmin = -180, xmax = 180, ymin = -60, ymax = 90} is selected.
#' @param normalize character or logical. Either logical on whether data should be normalized 
#' for the given interval or a character specifying a type of normalization. Type
#' default to "standard". Check \code{\link[gecko:normalize]{gecko::normalize()}}
#' for more info.
#' @param filepath character. An optional user defined path for the final output. If \code{NULL}, requested files are left in the current temp directory.
#' @return SpatRaster.
#' @examples
#' \dontrun{
#' regional_threats = spectre.area(3, terra::ext(-17.3,-16.6,32.6,32.9), normalize = FALSE)
#' terra::plot(regional_threats[[1]], main = "Human Density")
#' }
#' @export
spectre.area = function(index, ext = c(-180, 180, -60, 90), normalize = FALSE, filepath = NULL){
  if(!(all(index < 24) && all(index > 0) )){
    warning("Incorrect layer indexes chosen.")
    return(invisible())
  }
  # Check for incorrect extents
  if (isFALSE(ext[1] >= -180 && ext[2] <= 180 &&
              ext[3] >= -60  && ext[4] <= 90)){
    warning("Invalid coordinate selection.")
    return(invisible())
  } 
  if(length(normalize) > 1){
    warning("Normalize must be a single element.")
    return(invisible())
  } else {
    if(!(is(normalize, "logical") || c("standard", "rank", "range")%in%c(normalize)) ){
      warning("Normalize must be logical or a character representing normalization type.
              Check gecko::normalization for more info.")
      return(invisible())
    }
  }
  metadata = spectre.metadata()
  
  # max area: 54000, the cutoff is about at half so let's set it to that
  if (get.area(ext) >= 13500){
    continue = TRUE
    frac = 2
    while(continue){
      intervals = breakage(ext[1]:ext[2], frac)
      
      new_ext = lapply(1:length(intervals),
                       function(i){ c(intervals[[i]], ext[3:4]) })
      
      size_check = sapply( 1:length(new_ext), 
                           function(i){ get.area(new_ext[[i]]) >= 13500 }  )
      if (any(size_check)){
        continue = TRUE
        frac = frac + 1 
      } else {
        continue = FALSE
      }
    }
    original_ext = ext
    ext = new_ext
  }
  
  lout = list()
  for (i in index){
    # Adjust depending on whether or not we are using sub segments.
    if(exists("original_ext")){
      temp_path = paste(tempdir(), "\\", metadata$layer_name[i], "_",
                        paste(as.vector(original_ext), collapse = "_"), ".tif", 
                        sep = "")
      # Download our sub segments and save the locations of our files.
      subseg_paths = list()
      for (i in 1:length(ext)){
        subseg_paths[[i]] = area.download(terra::ext(ext[[i]]), i)
      }
      # Merge our sub segments back.
      merge_layer = terra::rast(subseg_paths[[1]])
      z = 1
      while(z < length(subseg_paths)){
        merge_layer = terra::merge(merge_layer, terra::rast(subseg_paths[[z+1]]), 
                                    filename = temp_path, overwrite = TRUE)
        z = z + 1
      }
      # Clear space
      sapply(open, FUN = unlink)
    } else {
      temp_path = area.download(ext, i)
    }
    
    if (!isFALSE(normalize)){
      if (isTRUE(normalize)){
        type = "standard"
      } else {
        type = normalize
      }
      normalize(terra::rast(temp_path), method = type, filepath = temp_path)
    }
    if(!is.null(filepath)){
      file.copy(from = temp_path, to = paste0(filepath, "\\", basename(temp_path)))
    }
    lout[[metadata[i,2]]] = terra::rast(temp_path)
  }
  return(lout)
}


#' Get SPECTRE data from points.
#' @description Downloads SPECTRE layer data according to a selection of points.
#' @param index numeric. A vector of integers specifying the layers. Refer to the documentation of 
#' \code{\link[gecko:spectre.citations]{gecko::spectre.citations()}} for a list
#' of available layers.
#' @param points data.frame or matrix. Containing point data coordinates, organized in longitude, latitude (longlat).
#' @return data.frame or matrix. Contains both the points given as well as 
#' their respective values for each layer specified.
#' @examples
#' \dontrun{
#' localities = gecko.data("records")
#' local_threats = spectre.points(c(2,3), localities)
#' }
#' @export
spectre.points = function(index, points){
  # Check for incorrect input types.
  if (!(is(points, "data.frame") || is(points, "matrix"))) {
    message("Error: Points must a data.frame or matrix object.")
    return(invisible())
  }
  # Check for incorrect input dimensions.
  if (ncol(points) != 2 || nrow(points) == 0) {
    message("Error: Incorrect input dimensions.")
    return(invisible())
  }
  # Check for out of bounds points.
  if (isFALSE(all(points[, 1] >= -180) && all(points[, 1] <= 180) &&
    all(points[, 2] >= -60) && all(points[, 2] <= 90))) {
    message("Error: One or more points are out of bounds.")
    return(invisible())
  }
  
  metadata = spectre.metadata()
  
  # Create an output matrix, assign the first two columns to our input points.
  output_mat <- matrix(
    nrow = nrow(points), ncol = 2 + length(index),
    dimnames = list(NULL, c("x", "y", metadata[index, 2]))
  )
  output_mat[, 1:2] <- as.matrix(points)

  # PB
  pb = utils::txtProgressBar(min = 1, max = length(index),
                      style = 3, width = 50, char = "=") 
  
  # Start working layer by layer.
  for (l in 1:length(index)) {
    utils::setTxtProgressBar(pb, l)
    
    layer <- substr(
      metadata[index[l], 3], 10,
      nchar(metadata[index[l], 3])
    )

    urls <- apply(points, 1, function(i) {
      paste("https://paituli.csc.fi/geoserver/wms?SERVICE=WMS&VERSION=1.1.1",
        "&REQUEST=GetFeatureInfo&FORMAT=image%2Fpng&TRANSPARENT=true",
        "&QUERY_LAYERS=paituli%3A", layer, "&LAYERS=paituli%3A", layer,
        "&INFO_FORMAT=text%2Fplain", "&X=", 0, "&Y=", 0,
        "&WIDTH=", 50, "&HEIGHT=", 50,
        "&SRS=EPSG%3A4326&STYLES=&BBOX=", i[1], "%2C", i[2] - 0.1,
        "%2C", i[1] + 0.1, "%2C", i[2],
        sep = ""
      )
    })

    vals <- tryCatch(
      {
        # account for single problems in requests

        vapply(urls, function(i) {
          
          if (Sys.info()["sysname"] == "Linux"){
            down = suppressMessages(download.file(i, file.path(tempdir(), "gecko_temp.txt"),
                                                  quiet = TRUE))
          } else {
            # Run bash to bypass certificates from geoserver and save.
            down = suppressMessages(system(paste('curl', i, '--ssl-no-revoke --output', 
                                                 file.path(tempdir(), "gecko_temp.txt")),
                                           show.output.on.console = FALSE))
          }
          # Return what you read from the file.
          return(as.double(scan(file.path(tempdir(), "gecko_temp.txt"),
            what = "double", quiet = TRUE
          )[9]))
        }, 0.1)
      },
      warning = function(cond) {
        message(paste0("Error: Could't find data for ", metadata[l,2], "."))
        message("Check your internet connection.")
        return(invisible())
      }
    )
    if (is.null(vals)) {
      return(invisible())
    }

    output_mat[, 2 + l] <- vals
  }
  # Switch the NA value of -3.4e+38 for a proper NA value
  output_mat <- apply(output_mat, 2, function(i) {
    i[i < -3.3e+38] <- NA
    return(i)
  })
  # this line does the same thing but assumes that columns will have NA on the same
  # rows.
  # TO DELETE:
  #output_mat[output_mat[, 3] < -3.3e+38, 3:ncol(output_mat)] <- rep(NA, ncol(output_mat) - 2)
  
  # Delete temp and return output.
  unlink(file.path(tempdir(), "gecko_temp.txt"))
  return(output_mat)
}
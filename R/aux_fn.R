#####GECKO - Geographical Ecology and Conservation Knowledge Online
#####Version 1.0.0 (2022-05-28)
#####By Vasco Branco, Pedro Cardoso, Luís Correia
#####Maintainer: vasco.branco@helsinki.fi
#####Changed from v0.1.2:
#####Full release. Added new SPECTRE functionalities.
#####Renamed gecko.examples to gecko.data

################################################################################
##############################AUX FUNCTIONS#####################################
################################################################################

#' Extent area.
#' @description Get the rough size of a certain raster extent.
#' @param x SpatExtent. An extent object specifying the range of a layer.
#' @noRd
get.area = function(x){
  return(abs(x[1] - x[2]) * abs(x[3] - x[4]))
}

#' Break vector
#' @description Split a vector of integers neatly according to a chosen fraction.
#' @param x An extent object specifying the range of a layer.
#' @param fraction An extent object specifying the range of a layer
#' @noRd
breakage = function(x, fraction){
  # If I had a task on a forloop going through a sequence of ints and
  # and I had to divide them neatly, e.g: to allocate the task to another PC. How
  # could I split them as evenly as possible?
  # i = 10:50
  # breakage(i, 3)
  ls = list()
  out = c(x[1])
  for (i in 1:(fraction-1)){
    out = c(out, out[length(out)] + length(x) %/% fraction)
  }
  out = c(out, x[length(x)])
  for (i in 1:fraction){
    ls[[i]] = c(out[i], out[i+1])  
  }
  return(ls)
  
}


#' Area download.
#' @description Creates maps ready to print in pdf or other formats.
#' @param ext An extent object specifying the range of a layer to download.
#' @param index An integer specifying a layer from spectre. See spectre.metadata.
#' @noRd
area.download = function(ext, index){
  metadata = spectre.metadata() 
  
  raster_url = paste("https://paituli.csc.fi/geoserver/wcs?version=2.0.1",
                     "&request=GetCoverage&coverageId=", metadata[index,3],
                     "&&subset=Long(", ext[1], ",", ext[2],
                     ")&subset=Lat(", ext[3], ",", ext[4],
                     ")&format=image/tiff", sep = "")
  
  # It is crucial that we download using "curl". I don't know why yet but if
  # we don't our data will be distorted. Additonally, sometimes downloading
  # fails because it can't get the certificate from geoserver. So far the 
  # resolution is to stop using download.file() entirely and download using
  # the "curl" bash command with --ssl-no-revoke
  file_name = paste(tempdir(), "\\", metadata[index,2], "_",
                    paste(as.vector(ext), collapse = "_"), ".tif", 
                    sep = "")
  
  if (Sys.info()["sysname"] == "Linux"){
    down = download.file(raster_url, file_name)
  } else {
    down = system(paste('curl', raster_url, '--ssl-no-revoke --output', 
                        file_name), show.output.on.console = FALSE)
  }
  
  if (down == 6){
    warning(paste("Could not resolve host. Service might be temporarily", 
                  "unavailable. Check your internet connection."))
    return(invisible())
  }
  
  tryCatch({terra::rast(file_name)},
           error = function(cond){
             stop("Raster layer is invalid. It might be corrupted or the download might have failed.")
             return(invisible())
           }
  )
  
  # Check here to see if there's a file with that name on the folder?
  
  return(file_name)
}

#' Euclidean distance between two vectors.
#' @description Extract geographic coordinates from strings containing location names
#' through usage of a gazzeteer for Python.
#' @param x A dataframe with the same formatting as 'predicted', containg
#'  some sort of classification data. 
#' @param y A dataframe with the same formatting as 'trained', containg
#'  the predicted classifications of a model trained over the data in 'trained'.
#' @noRd
dist.lite = function(x, y){
  # Euclidean distance ---------------------------------------------------------
  v_inter <- c(0)
  for (i in 1:length(x)) {
    v_inter <- c(v_inter + (as.numeric(x[, i]) - as.numeric(y[i]))^2)
  }
  out <- sqrt(v_inter)
  return(out)
  # Gower distance -------------------------------------------------------------
  # To be implemented in a future update
}



#' Mean euclidean distance between points
#' @description This function calculates the mean euclidean distance between a
#' point and all others in a given set, for every point.
#' @param x data.frame. With two columns containing latitude and longitude, describing
#' the predicted locations of a species, which may contain outliers.
#' @details Euclidean distance calculated is given by:
#' \eqn{d(p,q) =  \sqrt{\sum_{i=1}^{n}(q_i - p_i)²}},
#' with \eqn{p}, \eqn{q} being two points in Euclidean space, creating vectors 
#' \eqn{q_i}, and \eqn{p_i} and \eqn{n} the natural numbers set.
#' @examples
#' occurence = data.frame(X = runif(60), Y = runif(60))
#' gecko:::dist.euclid(occurence)
#' @return numeric.
#' @noRd
dist.euclid = function(x, threshold = 0.05){
  na_rows <- c(1:nrow(x))[rowSums(is.na(x[, 1:ncol(x)])) > 0]
  x <- apply(x, 2, FUN = as.double)
  x <- as.data.frame(x)
  
  if (length(na_rows) == 0) {
    order <- c(1:nrow(x))
  } else {
    order <- c(1:nrow(x))[-na_rows]
  }
  
  means <- rep(NA, nrow(x))
  
  for (r in order) {
    rest <- x[order, ]
    rest <- rest[row.names(rest) != as.character(r), ]
    out <- c()
    for (k in row.names(rest)) {
      out <- c(out, dist.lite(rest[k, ], x[as.character(r), ])) ###
    }
    means[r] <- mean(out)
  }
  
  qvalue <- quantile(means, 0.95, na.rm = TRUE)
  return(means >= qvalue)
}

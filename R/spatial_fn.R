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

#' Uniformize raster layers.
#' @description Crop raster layers to minimum size possible and uniformize \code{NA} values across layers.
#' @param layers SpatRaster. As defined in package terra, see \code{\link[terra:rast]{terra::rast()}}.
#' @details Excludes all marginal rows and columns with only \code{NA} values and change values to \code{NA} if they are \code{NA} in any of the layers.
#' @return SpatRaster. Same class as layers.
#' @examples region = gecko.data("layers")
#' terra::plot(clean(region))
#' @export
clean <- function(layers){
  ##apply mask to have NAs everywhere where any layer has NAs
  maskLayer <- sum(layers)
  maskLayer[!is.na(maskLayer)] <- 1
  layers <- terra::mask(layers, maskLayer)
  
  ##crop by excluding external rows and columns with NAs only
  layers <- terra::trim(layers)
  return(layers)
}


#' Create eastness layer.
#' @description Create a layer depicting eastness based on an elevation layer.
#' @param layers SpatRaster. A layer of elevation (a digital elevation model - DEM). 
#' As defined in package terra, see \code{\link[terra:rast]{terra::rast()}}.
#' @details Using elevation, aspect can be calculated. Yet, it is a circular variable (0 = 360) and has to be converted to northness and eastness to be useful for modelling.
#' @return SpatRaster.
#' @examples region = gecko.data("layers")
#' terra::plot(create.east(region[[3]]))
#' @export
create.east <- function(layers){
  asp <- terra::terrain(layers, v = "aspect")
  return(sin(asp))
}


#' Create latitude layer.
#' @description Create a layer depicting latitude based on any other.
#' @param layers SpatRaster. As defined in package terra, see \code{\link[terra:rast]{terra::rast()}}.
#' @details Using latitude (and longitude) in models may help limiting the extrapolation of the predicted area much beyond known areas.
#' @return SpatRaster.
#' @examples region = gecko.data("layers")
#' terra::plot(create.lat(region[[1]]))
#' @export
create.lat <- function(layers){
  if(dim(layers)[3] > 1){
    layers <- layers[[1]]
  }
  points <- terra::as.points(layers)[,1:2]
  lat <- terra::rasterize(terra::crds(points), y = layers, values = terra::crds(points)[,1] )
  lat <- terra::mask(lat, layers)
  names(lat) <- "latitude"
  return(lat)
}


#' Create longitude layer.
#' @description Create a layer depicting longitude based on any other.
#' @param layers SpatRaster. As defined in package terra, see \code{\link[terra:rast]{terra::rast()}}.
#' @details Using longitude (and latitude) in models may help limiting the extrapolation of the predicted area much beyond known areas.
#' @return SpatRaster.
#' @examples region = gecko.data("layers")
#' terra::plot(create.long(region))
#' @export
create.long <- function(layers){
  if(dim(layers)[3] > 1) {
    layers <- layers[[1]]
  }
  points <- terra::as.points(layers)[,1:2]
  long <- terra::rasterize(terra::crds(points), y = layers, values = terra::crds(points)[,2] )
  long <- terra::mask(long, layers)
  names(long) <- "longitude"
  return(long)
}


#' Create northness layer.
#' @description Create a layer depicting northness based on an elevation layer.
#' @param layers SpatRaster. A layer of elevation (a digital elevation model - DEM). 
#' As defined in package terra, see \code{\link[terra:rast]{terra::rast()}}.
#' @details Using elevation, aspect can be calculated. Yet, it is a circular variable (0 = 360) and has to be converted to northness and eastness to be useful for modelling.
#' @return SpatRaster.
#' @examples region = gecko.data("layers")
#' terra::plot(create.north(region[[3]]))
#' @export
create.north <- function(layers){
  asp <- terra::terrain(layers, v = "aspect")
  return(cos(asp))
}


#' Create distance layer.
#' @description Creates a layer depicting distances to records using the minimum, average, distance to the minimum convex polygon or distance taking into account a cost surface.
#' @param longlat matrix. Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param layers SpatRaster. As defined in package terra, see \code{\link[terra:rast]{terra::rast()}}. To serve as model to create distance layer.
#' @param type character. text string indicating whether the output should be the "minimum", "average" or "mcp" distance to all records. "mcp" means the distance to the minimum convex polygon encompassing all records.
#' @details Using distance to records in models may help limiting the extrapolation of the predicted area much beyond known areas.
#' @return SpatRaster.
#' @examples userpar <- par(no.readonly = TRUE) 
#' region = gecko.data("layers")
#' alt = region[[3]]
#' localities = gecko.data("records")
#' par(mfrow=c(3,2))
#' terra::plot(alt)
#' points(localities)
#' terra::plot(distance(localities, alt))
#' terra::plot(distance(localities, alt, type = "average"))
#' par(userpar)
#' @export
distance <- function(longlat, layers, type = "minimum"){
  if(dim(layers)[3] > 1){
    layers <- layers[[1]]
  }
  layers_template = terra::classify(!is.na(layers), c(TRUE, 0))
  
  if(type == "average"){
    for (d in 1:nrow(longlat)) {
      # layers <- layers + raster::distanceFromPoints(layers, x[d,]) #  he distance unit is in meters if the coordinate reference system (crs) of the Raster* object is (+proj=x) or assumed to be if the crs is NA. In all other cases it is in the units defined by the crs (which typically is meters).
      layers <- c(
        layers,
        terra::distance(
          layers_template,
          terra::vect(longlat[d, ],
                      geom = colnames(longlat),
                      crs = terra::crs(layers)
          )
        )
      )
    }
    layers = layers[[2:dim(layers)[3]]]
    # No mean?
    layers = terra::mean(layers)
    
    # layers <- layers/nrow(x)
    names(layers) <- "average distance"
  } else {
    #layers <- raster::mask(raster::distanceFromPoints(layers, x), layers)
    layers <- terra::mask(
      terra::distance(
        layers,
        terra::vect(longlat,
                    geom = colnames(longlat),
                    crs = terra::crs(layers)
        )
      ),
      layers
    )
    names(layers) <- "minimum distance"
  }
  return(layers)
}


#' Spatial thinning of occurrence records.
#' @description Thinning of records with minimum distances either absolute or relative to the species range.
#' @param longlat matrix. Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param distance numeric. Distance either in relative terms (proportion of maximum distance between any two records) or in raster units.
#' @param relative logical. If \code{TRUE}, represents the proportion of maximum distance between any two records. If \code{FALSE}, is in raster units.
#' @param runs numeric. Number of runs
#' @details Clumped distribution records due to ease of accessibility of sites, emphasis of sampling on certain areas in the past, etc. may bias species distribution models.
#' The algorithm used here eliminates records closer than a given distance to any other record. The choice of records to eliminate is random, so a number of runs are made and the one keeping more of the original records is chosen.
#' @return A matrix of species occurrence records separated by at least the given distance.
#' @examples userpar <- par(no.readonly = TRUE)
#' occ_points <- matrix(sample(100), ncol = 2)
#' par(mfrow=c(1,2))
#' graphics::plot(occ_points)
#' occ_points <- thin(occ_points, 0.1)
#' graphics::plot(occ_points)
#' par(userpar)
#' @export
thin <- function(longlat, distance = 0.01, relative = TRUE, runs = 100){
  longlat = longlat[!duplicated(longlat),]                #first, remove duplicate rows
  nSites = nrow(longlat)
  if(nSites < 4)
    return(longlat)
  
  ##if relative, calculate maxDist between any two points
  if(relative){
    if(nSites < 40){ #if limited number of sites use all data
      maxDist = 0
      for(i in 1:(nSites-1)){
        for(j in (i+1):nSites){
          maxDist = max(maxDist,((longlat[i,1]-longlat[j,1])^2+
                                   (longlat[i,2]-longlat[i,2])^2)^.5)
        }
      }
    } else { #if many sites use hypothenusa of square encompassing all of them
      horiDist = max(longlat[,1]) - min(longlat[,1])
      vertDist = max(longlat[,2]) - min(longlat[,2])
      maxDist = (horiDist^2 + vertDist^2)^0.5
    }
    distance = maxDist*distance
  }
  
  listSites = matrix(longlat[1,], ncol=2, byrow = TRUE)
  for (r in 1:runs){
    longlat = longlat[sample(nSites),]       ##shuffle rows (sites)
    rndSites = longlat[1,]                   ##start with first random site
    for(newSite in 2:nSites){
      for(oldSite in 1:(newSite-1)){
        addSite = TRUE
        dist = ((longlat[newSite,1]-longlat[oldSite,1])^2+
                  (longlat[newSite,2]-longlat[oldSite,2])^2)^.5
        if(dist < distance){
          addSite = FALSE
          break
        }
      }
      if(addSite)
        rndSites = rbind(rndSites, longlat[newSite,])
    }
    if(nrow(rndSites) > nrow(listSites))
      listSites = rndSites
  }
  return(as.matrix(listSites))
}


#' Move records to closest non-NA cell.
#' @description Identifies and moves presence records to cells with environmental values.
#' @param longlat matrix. Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param layers SpatRaster. As defined in package terra, see \code{\link[terra:rast]{terra::rast()}}.
#' @param buffer numeric. Maximum distance in map units that a record will move. If 0 all \code{NA} records will be changed.
#' @details Often records are in coastal or other areas for which no environmental data is available. This function moves such records to the closest cells with data so that no information is lost during modelling.
#' @return A matrix with new coordinate values.
#' @examples region <- terra::rast(matrix(c(rep(NA,100), rep(1,100), rep(NA,100)), ncol = 15))
#' presences <- cbind(runif(100, 0, 0.55), runif(100, 0, 1))
#' terra::plot(region)
#' points(presences)
#' presences <- move(presences, region)
#' terra::plot(region)
#' points(presences)
#' @export
move <- function(longlat, layers, buffer = 0){
  if(dim(layers)[3] > 1){
    layers <- layers[[1]]
  }
  
  if(is(longlat, "matrix")){
    longlat = as.data.frame(longlat)
  }
  
  if(terra::crs(layers) == ""){
    terra::crs(layers) = "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
  }
  
  # layers <- layers[[1]]
  values <- terra::extract(layers, longlat)   #get values of each record
  suppressWarnings(
    for(i in which(is.na(values))){    #if a value is NA, move it
      # Distance does not work when the crs is "". needs an exception
      distRaster <- terra::distance(
        layers,
        terra::vect(as.data.frame(longlat)[i, ], # remove if enforced at start
                    geom = colnames(longlat),
                    crs = terra::crs(layers)
        )
      )
      distRaster <- terra::mask(distRaster, layers)
      vmin <- terra::where.min(distRaster)
      
      if(buffer <= 0 || buffer > vmin){
        # vmin = terra::as.points(distRaster, function(i) i == vmin)
        longlat[i,] = terra::xyFromCell(distRaster, vmin[2]) # vmin[1,1:2]
      }
      
      
    }
  )
  return(longlat)
}

#' Reduce dimensionality of raster layers.
#' @description Reduce the number of layers by either performing a PCA on them or by eliminating highly correlated ones.
#' @param layers SpatRaster. As defined in package terra, see \code{\link[terra:rast]{terra::rast()}}.
#' @param method character. Either Principal Components Analysis ("pca", default) or Pearson's correlation ("cor").
#' @param n numeric. Number of layers to reduce to.
#' @param thres numeric. Value for pairwise Pearson's correlation above which one of the layers (randomly selected) is eliminated.
#' @details Using a large number of explanatory variables in models with few records may lead to overfitting. This function allows to avoid it as much as possible.
#' If both n and thres are given, n has priority. If method is not recognized and layers come from read function, only landcover is reduced by using only the dominating landuse of each cell.
#' @return SpatRaster.
#' @export
reduce <- function(layers, method = "pca", n = NULL, thres = NULL){
  ##method = "pca, cor", if unrecognized method only reduce landcover but not climate
  
  if(dim(layers)[3] == 33){          ##check if layers are obtained with read
    out <- c(layers[[33]])
    layers = layers[[1:19]]
  }
  # HANDLE COR
  if(method == "cor"){                       ##if correlation
    if(is.null(n)){
      if(is.null(thres))
        thres = 0.7
      for(i in 1:dim(layers)[3]){ ##delete layers until none are correlated above threshold
        # cor = as.matrix(as.dist(layerStats(layers, 'pearson', na.rm = TRUE)[[1]]))
        cor = as.matrix(
          as.dist(terra::layerCor(layers, 'pearson', na.rm = TRUE)[[1]])
        )
        
        if(max(cor) < thres)
          break
        corLayer = sample(which(cor == max(cor), arr.ind = TRUE)[,1],1)
        layers = layers[[-corLayer]]
      }
    } else {
      while (dim(layers)[3] > n){                   ##delete layers until reaching n layers
        cor = abs(as.matrix(as.dist(
          terra::layerCor(layers, 'pearson', na.rm = TRUE)[[1]])))
        
        corLayer = sample(which(cor == max(cor), arr.ind = TRUE)[,1],1)
        layers = layers[[-corLayer]]
      }
    }
  } else if(method == "pca"){                                  ##if pca
    if(is.null(n))
      n = 3
    if(sum(!is.na(terra::values(layers[[1]], mat = FALSE))) > 2000)
      sr <- terra::spatSample(layers, 1000)
    else
      sr <- terra::spatSample(layers, as.integer(sum(!is.na(terra::values(layers[[1]], mat = FALSE)))/2), na.rm = TRUE) # added na.rm
    pca <- prcomp(sr)
    layers <- terra::predict(layers, pca, index = 1:n)
    for(i in 1:n){
      names(layers[[i]]) <- paste("pca",i)
    }
  }
  
  if(dim(layers)[3] == 33){
    out <- c(layers, out)
  } else {
    out <- layers
  }
  return(out)
}


#' Visual detection of outliers.
#' @description Draws plots of sites in geographical (longlat) and environmental (2-axis PCA) space.
#' @param longlat matrix. Matrix of longitude and latitude or eastness and northness (two columns in this order) of species occurrence records.
#' @param layers SpatRaster. As defined in package terra, see \code{\link[terra:rast]{terra::rast()}}. It can be any set of environmental layers thought to allow the identification of environmental outliers.
#' @details Erroneous data sources or errors in transcriptions may introduce outliers that can be easily detected by looking at simple graphs of geographical or environmental space.
#' @return data.frame. Contains coordinate values and distance to centroid in pca. Two plots are drawn for visual inspection. The environmental plot includes row numbers for easy identification of possible outliers.
#' @examples localities = gecko.data("records")
#' region = gecko.data("layers")
#' outliers.visualize(localities, region[[1:3]])
#' @export
outliers.visualize <- function(longlat, layers){
  userpar <- par(no.readonly = TRUE) 
  on.exit(par(userpar))
  if(dim(layers)[3] == 33){
    pca <- reduce(layers[[1:19]], n = 2) #if layers come from read
  } else {
    pca <- reduce(layers, n = 2)
  }
  
  
  ##extract pca values from x
  pca <- as.data.frame(terra::extract(pca, longlat))
  goodRows <-  which(!is.na(pca[,1]))
  pca <- pca[goodRows,]
  longlat <- longlat[goodRows,]
  par(mfrow = c(1,2))
  red::map.draw(longlat, layers[[1]], spName = "Geographical")
  # plot(pca, main = "Environmental", type = "n")
  plot(pca, main = "Environmental")
  centroid = colMeans(pca)
  text(centroid[1], centroid[2], label = "X")
  for (i in 1:nrow(pca)){
    text(pca[i,1], pca[i,2], label = row.names(longlat)[i])
  }
  
  ##build new matrix ordered by distance to centroid
  dist2centroid = apply(pca, 1, function(i) dist(rbind(i, centroid)))
  out = as.data.frame(cbind(longlat, dist2centroid))
  out = out[order(-dist2centroid),]
  return(out)
}
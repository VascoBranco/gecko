###############################################################################
##############################AUX FUNCTIONS####################################
###############################################################################

#' Map creation.
#' @description Creates maps ready to print in pdf or other formats.
#' @param longlat Matrix of longitude and latitude or eastness and northness (two columns in this order) of each occurrence record.
#' @param layer RasterLayer object representing the presence/absence map for the species.
#' @param spName String of species name.
#' @param borders If TRUE country borders are drawn.
#' @param scale If TRUE a distance scale in km is drawn.
#' @param legend If TRUE the legend for the map is drawn.
#' @param sites If TRUE the record locations are drawn.
#' @param mcp If TRUE the minimum convex polygon representing the Extent of Occurrence is drawn.
#' @param print If TRUE a pdf is saved instead of the output to the console.
map.draw <- function(longlat = NULL, layer, spName,  borders = FALSE, scale = TRUE, legend = FALSE, sites = TRUE, mcp = FALSE, print = FALSE){
  worldborders <- NULL
  data(worldborders, envir = environment())
  if (borders){
    layer[layer == 0] <- NA
    raster::plot(layer, main = spName, legend = legend, xlab = "longitude", ylab = "latitude", col = "forestgreen")
    lines(worldborders)
  } else {
    raster::plot(layer, main = spName, legend = legend, colNA = "lightblue", xlab = "longitude", ylab = "latitude")
  }
  if (scale){
    width = (raster::xmax(layer) - raster::xmin(layer))
    d = round(width/10^(nchar(width)-1))*10^(nchar(width)-2)
    raster::scalebar(d = d, type="bar", divs = 2)
  }
  if (sites && !is.null(longlat))
    points(longlat, pch = 19)
  if (mcp){
    e <- rasterToPoints(layer, fun = function(dat){dat == 1})   ##convert raster to points
    vertices <- chull(e[,1], e[,2])
    vertices <- c(vertices, vertices[1])
    vertices <- e[vertices,c(1,2)]
    poly <- SpatialPolygons(list(Polygons(list(Polygon(vertices)),1)))
    raster::plot(poly, add = TRUE)
  }
  if(print){
    dev.copy(device = pdf, file = paste(toString(spName), ".pdf", sep=""))
    dev.off()
  }
}

#' Extent of Occurrence (EOO).
#' @description Calculates the Extent of Occurrence of a species based on either records or predicted distribution.
#' @param spData spData One of three options: 1) matrix of longitude and latitude (two columns) of each occurrence record; 2) matrix of easting and northing (two columns, e.g. UTM) of each occurrence record in meters;  3) RasterLayer object of predicted distribution (either 0/1 or probabilistic values).
#' @details EOO is calculated as the minimum convex polygon covering all known or predicted sites for the species.
#' @return A single value in km2 or a vector with lower confidence limit, consensus and upper confidence limit (probabilities 0.975, 0.5 and 0.025 respectively).
eoo <- function(spData){
  if(is(spData, "RasterLayer")){
    if(!all(raster::as.matrix(spData) == floor(raster::as.matrix(spData)), na.rm = TRUE)){ #if probabilistic map
      upMap <- raster::reclassify(spData, matrix(c(0,0.025,0,0.025,1,1), ncol = 3, byrow = TRUE))
      consensusMap <- raster::reclassify(spData, matrix(c(0,0.499,0,0.499,1,1), ncol = 3, byrow = TRUE))
      downMap <- raster::reclassify(spData, matrix(c(0,0.975,0,0.975,1,1), ncol = 3, byrow = TRUE))
      area <- c(eoo(downMap), eoo(consensusMap), eoo(upMap))
    } else {
      if (raster::xmax(spData) <= 180) {  #if longlat data
        e <- rasterToPoints(spData, fun = function(dat){dat == 1})   ##convert raster to points
        vertices <- chull(e[,1], e[,2])
        if(length(vertices) < 3) return(0)
        vertices <- c(vertices, vertices[1])
        vertices <- e[vertices,c(1,2)]
        area = geosphere::areaPolygon(vertices)/1000000
      } else {
        spData[spData < 1] <- NA
        spData <- rasterToPoints(spData)
        vertices <- chull(spData)
        if(length(vertices) < 3) return(0)
        vertices <- c(vertices, vertices[1])
        vertices <- spData[vertices,]
        area = 0
        for(i in 1:(nrow(vertices)-1))
          area = area + (as.numeric(vertices[i,1])*as.numeric(vertices[(i+1),2]) - as.numeric(vertices[i,2])*as.numeric(vertices[(i+1),1]))
        area = abs(area/2000000)
      }
    }
  } else if (ncol(spData) == 2){
    vertices <- chull(spData)
    if(length(vertices) < 3) return(0)
    vertices <- c(vertices, vertices[1])
    vertices <- spData[vertices,]
    if(max(spData) <= 180) {  #if longlat data
      area = geosphere::areaPolygon(vertices)/1000000
    } else { #if square data in meters
      area = 0
      for(i in 1:(nrow(vertices)-1))
        area = area + (as.numeric(vertices[i,1])*as.numeric(vertices[(i+1),2]) - as.numeric(vertices[i,2])*as.numeric(vertices[(i+1),1]))
      area = abs(area/2000000)
    }
  } else {
    return(warning("Data format not recognized"))
  }
  return(round(area))
}
# -------------------- HOMERANGE FUNCTIONS -------------------------------------
# General functions for creation, analysis, and simulation of home ranges
# ------------------------------------------------------------------------------

# CreateHomeRangeKernels Function ----------------------------------------------

###  Creates RasterLayers of homerange kernels based on home centroids
###  Usage: CreateHomeRangeKernels(df_all, df_home, base, max_r, 
###   home_inflection, home_scale, avoid_inflection, avoid_scale, output_dir, 
###   write_dist, write_homerange)
###  Arguments: df_all = dataframe of all homerange centroids 
###             df_home = dataframe of homerange centroids to calculate
###               homerange kernels (these must be a subset of the df_all 
###               dataframe), Default is to use 'df_all' dataframe
###             base = base Raster that sets the projection, extent, and 
###               dimensions of the study area 
###             max_r = maximum radius to calculate the homerange raster from 
###               each df_home centroid
###             home_inflection = inflection point of the Logistic function that
###               governs the home kernel
###             home_scale = scale parameter of the Logistic function that 
###               governs the scale parameter
###             avoid_inflection = inflection point of the Logistic function
###               that governs the conspecific avoidance kernel
###             avoid_scale = scale parameter of the Logistic function that 
###               governs the conspecific avoidance kernel  
###  Returns: A list containing homerange kernel Rasters for all the df_home 
###    centroids
###  Notes: 
###  Blake Massey and Javan Bauder
###  2015.02.06 

CreateHomeRangeKernels <- function(df_all,
                                   df_home = df_all,
                                   base,
                                   max_r,
                                   home_inflection, 
                                   home_scale,
                                   avoid_inflection, 
                                   avoid_scale,
                                   output_dir,
                                   write_distance = FALSE,
                                   write_homerange = FALSE) {
  source('C:/Work/R/Functions/sim/move.R')
  cellsize <- res(base)[1]
  max_r_cells <- ceiling(max_r/cellsize)
  xmin <- xmin(base)
  ymin <- ymin(base)
  df_sp <- SpatialPointsDataFrame(df_all[,c("x","y")], df_all, 
    proj4string=crs(base))
  writeLines(noquote(paste("Calculating global distance")))
  if (write_distance == TRUE) {
    ifelse(exists("i"), i <- i , i <- 1)
    filename <- paste0(output_dir,"/global_dist_", sprintf("%03d", i), ".tif")
    global_dist <- distanceFromPoints(base, df_sp, filename=filename, 
        overwrite=TRUE)
    writeLines(noquote(paste("Writing:", filename)))
  } else {
    global_dist <- distanceFromPoints(base, df_sp)
  }
  homerange_ids <- df_home$nest_id
  total <- length(homerange_ids)
  homerange_kernels <- as.list(setNames(rep(NA,length(homerange_ids)), 
    homerange_ids), homerange_ids)   
  for (j in 1:nrow(df_home)) {
    writeLines(noquote(paste("Calculating homerange", j, "of", total)))
    home <- df_home[j,]
    home_sp <- SpatialPointsDataFrame(home[,c("x","y")], home, 
      proj4string=crs(base))
    xy <- CenterXYInCell(home_sp@coords[,"x"], home_sp@coords[,"y"], 
      xmin, ymin, cellsize)
    cell_extent <- extent(xy[1]-(cellsize/2), xy[1]+(cellsize/2), xy[2]-
      (cellsize/2), xy[2]+(cellsize/2))
    cell <- setValues(raster(cell_extent, crs=projection(base), res=cellsize),j)
    home_ext <- extend(cell, c(max_r_cells, max_r_cells), value=NA)
    home_dist <- distanceFromPoints(home_ext, home[,c("x","y")])
    home_kern <- calc(home_dist, fun = function(x){(1/(exp((-(x - 
      home_inflection)) / home_scale) + 1))})
    global_dist_crop <- crop(global_dist, home_dist)
    cent_dist <- overlay(home_dist, global_dist_crop, fun=function (x,y)
      {ifelse(x != y, NA, x)})
    cent_bounds <- boundaries(cent_dist)
    cent_bounds <- subs(cent_bounds, data.frame(from=c(0,1), to=c(NA,1)))
    edge_dist_abs <- distance(cent_bounds)
    edge_dist <- overlay(cent_dist, edge_dist_abs, fun=function(x,y)
      {ifelse(!is.na(x), y*-1, y*1)})    
    avoid_kern <- calc(edge_dist, fun = function(x){(1/(exp((-(x - 
      avoid_inflection)) / avoid_scale) + 1))})
    homerange_kern <- overlay(avoid_kern, home_kern, fun=function(x,y){
      p <- y*x #    y+x-1
    #       ifelse(p<=0, 0, p)
      })
    if (write_distance == TRUE) {
      k <- home[,"id"]
      filename <- paste0(output_dir,"/homerange_", k, "_",sprintf("%03d",
        i), ".tif")  
      writeLines(noquote(paste0("Writing: ", filename)))
      writeRaster(homerange_kern, filename=filename, overwrite=TRUE)
    }
    homerange_kernels[[j]] <- homerange_kern
    names(homerange_kernels[[j]]) <- df_home[j,"nest_id"]
  }
  return(homerange_kernels)
}



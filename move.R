# ----------------------------- MOVE FUNCTIONS ---------------------------------
# Movement process functions for individual-based model simulation
# ------------------------------------------------------------------------------

# CenterXYInCell Function ------------------------------------------------------

###  Centers x and y values into the center of a raster cell based on the xmin, 
###    ymin, and cellsize parameters
###  Usage: CenterXYInCell(x, y, xmin, ymin, cellsize)
###  Arguments: x = x value
###             y = y value
###             xmin = minimum value of x that establishes the grid arrangement
###             ymin = minimum value of y that establishes the grid arrangement
###             cellsize = cell size in units of the x and y values
###  Returns: vector of centered x and y (i.e., (x,y))   
###  Notes:
###  Blake Massey
###  2015.01.07

CenterXYInCell <- function(x, 
                           y, 
                           xmin, 
                           ymin, 
                           cellsize) { 
    x <- xmin + (floor(((x-xmin)/cellsize))*cellsize) + (cellsize/2)
    y <- ymin + (floor(((y-ymin)/cellsize))*cellsize) + (cellsize/2)
    output <- c(x, y)
    return(output)
}

# ConvertAngle Function --------------------------------------------------------

###  Converts a radian angle that is outside of the Unit Circle range (i.e., 
###    [0, 2pi]) to the equivalent value that within the Unit Circle range. 
###  Usage: ConvertAngle(x)
###  Arguments: x = radian value
###  Returns: vector of one value 
###  Notes: Used in several simulations so that the angles never go outside of 
###    the Unit Circle range when angles are being added or subtracted from each 
###    other over and over during a simulation. 
###  Blake Massey
###  2015.01.07

ConvertAngle <- function(x) {
  if (x > 2*pi) x <- x-(2*pi)
  if (x < 0) x <- x+(2*pi)
  return(x)
}

# CreateParetoKernel Function --------------------------------------------------

###  Creates a probability matrix based on a Pareto distribution.
###  Usage: CreateParetoKernel(scale, shape, max_r, cellsize)
###  Arguments: scale = scale parameter of Pareto distribution
###             shape = shape parameter of Pareto distribution
###             max_r = maximum radius of kernel in meters, default = 100
###             cellsize = cell size in meters, default = 1
###  Returns: RasterLayer  
###  Notes:
###  Javan Bauder and Blake Massey
###  2015.01.03

CreateParetoKernel <- function(scale, 
                               shape, 
                               max_r = 100, 
                               cellsize = 1) {
  require("VGAM")
  max_r_cells <- max_r/cellsize
  size = ceiling(max_r_cells) * 2 + 1
  center = ceiling(max_r_cells) + 1
  kernel <- new("matrix", 0, size, size)
  for (i in 1:size) for (j in 1:size) {
    r = sqrt((i - center)^2 + (j - center)^2) * cellsize  
    if (r <= max_r) 
      kernel[i, j] <- VGAM::dgpd(r,scale=scale,shape=shape,log=FALSE)
  }
  kernel[center, center] <- 1/scale  
  kernel <- kernel / sum(kernel)
  # This last part deletes the cells at the edge if they are all zero
  if (all(kernel[1, ] == 0, kernel[, 1] == 0, 
          kernel[nrow(kernel),] == 0, kernel[, ncol(kernel)] == 0)) {
    kernel <- kernel[2:(nrow(kernel) - 1), 2:(ncol(kernel) - 1)]
  }
  return(kernel)
}

# CreateRedistKernel Function -------------------------------------------------

###  Create a redistribution kernel matrix based on a wrapped Cauchy 
###    distribution for direction and a Pareto distribution for distance.
###  Usage: CreateRedistKernel(max_r, cellsize, mu, rho, shape, scale, 
###    ignore_cauchy, ignore_pareto)
###  Arguments: max_r = maximum radius of kernel in meters, default = 300
###             cellsize = cell size in meters, default = 30
###             mu = mu parameter of wrapped Cauchy distribution, 0 radians is 
###               due east because everything is based on the Unit Circle
###             rho = rho parameter of wrapped Cauchy distribution
###             shape = shape parameter of Pareto distribution
###             scale = scale parameter of Pareto distribution
###             ignore_cauchy = logical, removes cauchy kernel's contribution to 
###               output raster. Default is FALSE.
###             ignore_pareto = logical, removes pareto kernel's contribution to
###               output raster. Default is FALSE.
###  Returns: matrix
###  Notes:
###  Javan Bauder and Blake Massey
###  2015.01.03

CreateRedistKernel <- function(max_r = 300,
                               cellsize = 30,
                               mu,
                               rho,
                               shape,
                               scale, 
                               ignore_cauchy = FALSE,
                               ignore_pareto = FALSE) {
  suppressPackageStartupMessages(require(circular))
  suppressPackageStartupMessages(require(raster))
  suppressPackageStartupMessages(require(texmex))
  AngleToPoint <- function(origin_x, 
                         origin_y, 
                         target_x, 
                         target_y){
    dx <- c(target_x - origin_x)
    dy <- c(target_y - origin_y)
    abs_angle <- atan2(dy, dx)
    abs_angle <- ifelse(abs_angle < 0, (2*pi) + abs_angle, abs_angle)
  } 
  # Create the empty kernel objects
  max_r_cells <- ceiling(max_r/cellsize)
  size <- max_r_cells * 2 + 1
  center <- max_r_cells + 1
  wrpc_kernel <- new("matrix", 0, size, size)
  gpd_kernel <- new("matrix", 0, size, size)
  for (i in 1:size) {
    for (j in 1:size) {
      r = sqrt((i - center)^2 + (j - center)^2) * cellsize  
      b = AngleToPoint(center, center, j, i)
      if(r <= max_r){
        wrpc_kernel[i, j] <- round(suppressWarnings(dwrappedcauchy(b, mu=mu, 
          rho=rho)), 5)
        gpd_kernel[i, j] <- dgpd(r, sigma=scale, xi=shape, log=FALSE)
      }
    }
  } 
  wrpc_kernel <- apply(wrpc_kernel, 2, rev)
  gpd_kernel[center, center] <- 1/scale
  # This last part deletes the cells at the edge if they are all zero
  if (all(wrpc_kernel[1, ] == 0, wrpc_kernel[, 1] == 0, 
    wrpc_kernel[nrow(wrpc_kernel),] == 0, wrpc_kernel[, ncol(wrpc_kernel)] ==0))
    wrpc_kernel <- wrpc_kernel[2:(nrow(wrpc_kernel) - 1), 2:(ncol(wrpc_kernel) 
      - 1)]
  if (all(gpd_kernel[1, ] == 0, gpd_kernel[, 1] == 0, 
    gpd_kernel[nrow(gpd_kernel),] == 0, gpd_kernel[, ncol(gpd_kernel)] == 0))
    gpd_kernel <- gpd_kernel[2:(nrow(gpd_kernel) - 1), 2:(ncol(gpd_kernel) - 1)]
  # Multiply the two kernels together and re-normalize
  if (ignore_cauchy) wrpc_kernel <- 1
  if (ignore_pareto) gpd_kernel <- 1
  redist_kernel <- gpd_kernel*wrpc_kernel
  redist_kernel <- redist_kernel/sum(redist_kernel)
  return(redist_kernel)
}

# MovementSubModel Function ----------------------------------------------------

###  Movement submodel that predicts movements based on current location, 
###    homerange_kernel, nest_return probability, and others.
###  Usage: MovementSubModel()
###  Arguments: agent_states = 'agent_states' list object
###             step_data = step_data dataframe
###             step = step interval
###  Returns: 'step_data' object 
###  Notes: 
###  Blake Massey
###  2015.03.28

MovementSubModel <- function(agent_states = agent_states,
                             step_data = step_data,
                             step = step) {
  suppressPackageStartupMessages(require(circular))
  suppressPackageStartupMessages(require(ggplot2))
  suppressPackageStartupMessages(require(grid))
  suppressPackageStartupMessages(require(raster))
  suppressPackageStartupMessages(require(rasterVis))
  suppressPackageStartupMessages(require(sampling))
  suppressPackageStartupMessages(require(VGAM))
  source('C:/Work/R/Functions/gen.R')
  source('C:/Work/R/Functions/gps.R')  
  source('C:/Work/R/Functions/pars.R')
  source('C:/Work/R/Functions/sim/home.R') 
  cellsize <- res(sim$spatial$base)[1]
  step_start <- int_start(step)
  step_end <- int_end(step)
  sex <- agent_states$sex
  season <- FindSeasonFromDatetime(step_start, sim$pars$global$sim_seasons)
  nest_return <- sim$pars$classes[[sex]]$julian[yday(int_start(step)), 
    "nest_return"]
  max_r <- sim$pars$classes[[sex]]$season[[season]]$step_max_r
  step_cauchy_mu <- sim$pars$classes[[sex]]$constant$fixed$step_cauchy_mu
  step_cauchy_rho <- sim$pars$classes[[sex]]$constant$fixed$step_cauchy_rho
  step_pareto_shape <-sim$pars$classes[[sex]]$season[[season]]$step_pareto_shape
  step_pareto_scale <-sim$pars$classes[[sex]]$season[[season]]$step_pareto_scale
  homerange_kernel <- sim$spatial$homerange_kernel[[agent_states$nest_id]]
  if (nrow(step_data) == 1) {
    i <- 1
    step_data[i, "datetime"] <- int_start(step)
    step_data[i+1, "datetime"] <- int_end(step)
    x <- CenterXYInCell(step_data[i, "x"], step_data[i, "y"], 
      xmin(homerange_kernel), ymin(homerange_kernel), cellsize)[1]
    y <- CenterXYInCell(step_data[i, "x"], step_data[i, "y"], 
      xmin(homerange_kernel), ymin(homerange_kernel), cellsize)[2]
    step_data[i, "exp_angle"] <- sample(x=seq(from=0, to=(2*pi), by=(2*pi/360)), 
      size=1)     
    redist <- CreateRedistKernel(max_r=max_r, cellsize=cellsize, mu=
      step_data[i, "exp_angle"], rho=step_cauchy_rho, shape=step_pareto_shape, 
      scale=step_pareto_scale)
    r <- (cellsize*((nrow(redist)-1)/2))+(cellsize/2)
    redist_raster <- raster(redist, xmn=-r, xmx=r, ymn=-r, ymx=r)
    redist_shift <- shift(redist_raster, x=x, y=y)
### PLACE TO ADD IN OTHER PROBABILITY LAYERS
                
    homerange_crop <- crop(homerange_kernel, redist_shift)
    prob_raster <- overlay(redist_shift, homerange_crop, fun=function(x,y) 
      {return(x*y)}, recycle=FALSE) 
              
### END OF OTHER PROBABILITY LAYERS
    destination_cell <- suppressWarnings(strata(data=data.frame(cell=
      1:ncell(prob_raster)), stratanames=NULL, size=1, method="systematic", 
      pik=prob_raster@data@values))
    destination_xy <- as.vector(xyFromCell(prob_raster,destination_cell[1,1]))
    step_data[i+1, "x"] <- destination_xy[1]
    step_data[i+1, "y"] <- destination_xy[2]
    step_data$abs_angle[i] <- CalculateAngleToPoint2(step_data$x[i], 
      step_data$y[i], step_data$x[i+1], step_data$y[i+1])
  } else {
    i <- nrow(step_data)
    step_data[i+1, "datetime"] <- int_end(step)
    go_nest <- rbinom(1, 1, nest_return) 
    if (go_nest == TRUE) {
      nests_df <- data.frame(sim$spatial$nests)
      nests_xy <- nests_df[which(nests_df$nest_id == agent_states$nest_id),
        c("x", "y")]
      nest_xy <- CenterXYInCell(nests_xy[1], nests_xy[2], xmin(homerange_kernel), 
        ymin(homerange_kernel), cellsize)
      step_data$x[i+1] <- nest_xy[[1]]
      step_data$y[i+1] <- nest_xy[[2]]        
      step_data$abs_angle[i] <- CalculateAngleToPoint2(step_data$x[i], 
        step_data$y[i], step_data$x[i+1], step_data$y[i+1])    
    } else {
      step_data$exp_angle[i] <- step_data$abs_angle[i-1]
      redist <- CreateRedistKernel(max_r=max_r, cellsize=cellsize, 
        mu=step_data$exp_angle[i], rho=step_cauchy_rho, shape=step_pareto_shape, 
        scale=step_pareto_scale)
      r <- (cellsize*((nrow(redist)-1)/2))+(cellsize/2)
      redist_raster <- raster(redist, xmn=-r, xmx=r, ymn=-r, ymx=r)
      redist_shift <- shift(redist_raster, x=step_data$x[i], y=step_data$y[i])
        
### PLACE TO ADD IN OTHER PROBABILITY LAYERS
                
      homerange_crop <- crop(homerange_kernel, redist_shift, snap="out")
      prob_raster <- overlay(redist_shift, homerange_crop, fun=function(x,y) 
        {return(x*y)}, recycle=FALSE) 
                
### END OF OTHER PROBABILITY LAYERS
      destination_cell <- suppressWarnings(strata(data=data.frame(cell=
        1:ncell(prob_raster)), stratanames=NULL, size=1, method="systematic", 
        pik=prob_raster@data@values))
      destination_xy <- xyFromCell(prob_raster, destination_cell[1,1])
      step_data[i+1, "x"] <- destination_xy[1]
      step_data[i+1, "y"] <- destination_xy[2]
      step_data$abs_angle[i] <- CalculateAngleToPoint2(step_data$x[i], 
        step_data$y[i], step_data$x[i+1], step_data$y[i+1])
    }
  } 
  step_data[i+1, "id"] <- step_data[i, "id"]
  return(step_data)
}
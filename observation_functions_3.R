# edited 27th Feb 2021!
# functions to determine which observation locations yield presence and absence

library(sp)




set_observation_squares <- function(r, square_size, selection_type, N_observation_locs, habitat_cells, N_observation_cells, treshold_distance, habitat, systematic){
   
    # OccupiedLocations: binary list of observed presences and absences (this should be compared with the data)
    # Females: how many females have their core areas in the location
    # Males: how many males have their core areas there. 
    # FemalesMissed: how many Females are left totally out of the study area 
    # MalesMissed: (same)
  
  if (is.empty(treshold_distance)){
    treshold_distance <- 3
  }
    
    ObsSquares <- list()
    
    if (FALSE){
    N_observation_cells <- 100
    }
  
    rr <- r
    if (habitat== TRUE){
      rr[rr==0] <- NA
    }
    W <-  c(0,100,0,100)
    
    if (systematic==FALSE){
      r_sample <- as.data.frame(sampleRandom(x=rr, size = 20*N_observation_cells, na.rm = TRUE, xy = TRUE))
      new_sample <- discrete.sample(cbind(r_sample$x, r_sample$y), N_observation_cells, treshold_distance, k = 0)
    }else{
      r_sample <- as.data.frame(sampleRegular(x=rr, size = N_observation_cells, na.rm = TRUE, xy = TRUE))
      #new_sample <- discrete.sample(cbind(r_sample$x, r_sample$y), N_observation_cells, treshold_distance, k = 0) 
    }
    
    if (FALSE){
      plot(r)
      points(r_sample$x, r_sample$y)
    }
    
    for (i in seq(1, N_observation_cells)){
      
      #buf <- gBuffer(ObservationLocations[i], width = ObservationRadius)
      radius <- square_size/2
      if (systematic==TRUE){
        centroidy <- r_sample$y[i]
        centroidx <- r_sample$x[i]}else{
          centroidy <- new_sample[i, 2]
          centroidx <- new_sample[i, 1]  
        } 
      
      yPlus <- centroidy +radius
      xPlus <- centroidx +radius
      yMinus <- centroidy-radius
      xMinus <- centroidx-radius
      
      square=rbind(cbind(xMinus,yPlus),  # NW corner
                   cbind(xPlus, yPlus),  # NE corner
                   cbind(xPlus,yMinus),  # SE corner
                   cbind(xMinus,yMinus), # SW corner
                   cbind(xMinus,yPlus))  # NW corner again - close ploygon
      
      s <- Polygon(square)
      buf <- SpatialPolygons(list(Polygons(list(s), i)))
      ObsSquares[i] <- buf
    }
    
    joined = SpatialPolygons(lapply(ObsSquares, function(x){x@polygons[[1]]}))
    return(joined)
}

set_observation_transects <- function(r, square_size, selection_type, N_observation_locs, habitat_cells, N_observation_cells, treshold_distance, habitat, systematic){
  
  # OccupiedLocations: binary list of observed presences and absences (this should be compared with the data)
  # Females: how many females have their core areas in the location
  # Males: how many males have their core areas there. 
  # FemalesMissed: how many Females are left totally out of the study area 
  # MalesMissed: (same)
  
  buffer <- 20
  xmin <- r@extent@xmin+buffer
  xmax <- r@extent@xmax-buffer
  
  
  if (is.empty(treshold_distance)){
    treshold_distance <- 3
  }
  
  ObsTransects <- list()
  
  if (FALSE){
    N_observation_cells <- 100
  }
  
  rr <- r
  if (habitat== TRUE){
  rr[rr==0] <- NA
  }
  
  if (systematic==TRUE){
  locs <- seq(from= xmin, to =xmax, length.out = N_observation_locs)
  }else{
  locs <- seq(from= xmin, to =xmax, length.out = N_observation_locs)
  }

  for (i in seq(1, N_observation_locs)){

    xPlus <- locs[i]
    xMinus <- locs[i]
    yPlus <- xmax+buffer
    yMinus <- 0
    
    transect=rbind(cbind(xMinus, yMinus),  # NW corner
                 cbind(xPlus, yPlus))
    
    
    s <- Line(transect)
    buf <- SpatialLines(list(Lines(list(s), i)))
    ObsTransects[i] <- buf
  }
  
  joined = SpatialLines(lapply(ObsTransects, function(x){x@lines[[1]]}))
  return(joined)
}

set_individual_locations <- function(r, home_ranges, Nhabitat, Nindividuals, always_in_habitat, N_locations, habitat_dependency){
  
  # samples random spatial points at which the individual is located
  # always_in_habitat is either TRUE or FALSE and it defines whether the locations at which 
  # individuals are observable are always in habitat
  # N_locations can be >=1. If one, we interpret this as the individual being seen, if larger could be interpreted as droppings,
  # resting sites.. etc
  # habitat_dependency determines if the N_locations depends on the size of the home range, or if it is a constant
  
  dx <- diff(c(xmin(r), xmax(r))) / ncol(r) / 2 # Half of horizontal width
  dy <- diff(c(ymin(r), ymax(r))) / nrow(r) / 2 # Half of vertical width
  n <- 1
  
  track_individual_id <- numeric(0)
  individual_coords_x <- numeric(0)
  individual_coords_y <- numeric(0)
  cc <- 1
  Nindividuals <- length(home_ranges)
  #print(N_locations)
  if (habitat_dependency ==FALSE){
    for (i in seq(1,Nindividuals)){
      #print(home_ranges)
      #print(home_ranges)
      #print(home_ranges[[i]])
      cells_of_obs <- sample(home_ranges[[i]], N_locations)
      xy <- xyFromCell(r, cells_of_obs)
      xy <- xy + c(runif(n, -dx, dx), runif(n, -dy, dy)) 
      for (j in seq(1, N_locations)){
        track_individual_id[cc] <- i
        individual_coords_x[cc] <- as.numeric(xy[j,1])
        individual_coords_y[cc] <- as.numeric(xy[j,2])
        cc <- cc+1
      }
    }
  }
  
  if (habitat_dependency ==TRUE){
    N_locations_v <- numeric(0)
    for (i in seq(1,Nindividuals)){
      N_locations_v[i] <- ceiling(length(home_ranges[[i]])*N_locations) 
      cells_of_obs <- sample(home_ranges[[i]], N_locations_v[i])
      xy <- xyFromCell(r, cells_of_obs)
      xy <- xy + c(runif(n, -dx, dx), runif(n, -dy, dy)) 
      for (j in seq(1, N_locations_v[i])){
        track_individual_id[cc] <- i
        #print(xy)
        #print(cc)
        individual_coords_x[cc] <- as.numeric(xy[j,1])
        individual_coords_y[cc] <- as.numeric(xy[j,2])
        cc <- cc+1
      }
    }
    
  }
  
  res <- list('individual_locations_xy' = cbind(individual_coords_x, individual_coords_y), 'individual_ids' = track_individual_id)
  return(res)
  
}

allocate_locations_within_squares <- function(individual_coords, ObsSquares){
 N_coords<- length(individual_coords) 
 Nsquares <-length(ObsSquares)
 locations_per_square <- numeric(0)
 locations_per_square[seq(1, Nsquares)] <- 0
 
 individual_coords <- SpatialPoints(individual_coords)
 squares <- numeric(0) # for each coordinate, gives id of the observation square it falls within
 query_result <- as.numeric(over(individual_coords, ObsSquares, returnList = FALSE))
 for (i in seq(1, Nsquares)){
     locations_per_square[i] <- sum(1*(query_result==i),na.rm = T)
 }
 res <- list('squares_for_locations' = query_result, 'locations_per_square' = locations_per_square)
 return(res)
}

compute_distances_to_transects <- function(individual_coords, ObsTransects){
  
  require(rgeos)
  individual_coords <- SpatialPoints(individual_coords)
  Ncoords<- length(individual_coords) 
  Ntransects <- length(ObsTransects)

  distance_matrix <- matrix(nrow = Ncoords, ncol= Ntransects)
  
  for (i in seq(1, Ncoords)){# for each coordinate, gives id of the observation square it falls within
  for (j in seq(1, Ntransects)){
  distance_matrix[i,] <- gDistance(individual_coords[i], ObsTransects[j]) #  get a matrix of distances from each point to each line
  }
  }
  
  #distance_matrix[i] <- query_result
  #res <- list('distance_matrix' = distance_matrix)
  return(distance_matrix)
}

observe_individuals_within_squares <- function(locations_per_square, observation_probability, Nsquares){
 observed_tracks <-  rbinom(rep(1, Nsquares),locations_per_square,  observation_probability)
 #print(observed_tracks)
 res <- list('presence_observed' = 1*(observed_tracks>0), 'number_of_tracks_observed' = observed_tracks)
 return(res)
}
  
half_normal <- function(distance){
  prob <- exp(-1*(distance)^2/(2*(0.5)^2))
  return(prob)
}

observe_individuals_along_transects_2 <- function(Transects, individual_xys){
  # rows = individuals, 
  # cols = transects
  # observed_tracks counts how many presences observed per each transect:
  Ntransects <- length(Transects)
  Nindividuals <- length(individual_xys[,1])
  observed_individuals <- rep(0,Ntransects)
  observed_sightings <- rep(0,Ntransects)
  for (i in seq(1, Ntransects)){
    inds_observed <- rep(0, Nindividuals)
    transect <- Transects[i]
    x <- extent(transect)@xmin
    ymin <- extent(transect)@ymin
    ymax <- extent(transect)@ymax
    ypoints <- seq(from= ymin, to =ymax, by = 1)
    for (yy in ypoints){
      distances_x <- abs(individual_xys[,1]-x)
      distances_y <- abs(individual_xys[,2]-yy)
      angles <- atan(distances_x/distances_y)
      dist_to_transect <- distances_y/cos(angles)
      probs <- half_normal(dist_to_transect)
      sightings <- rbinom(Nindividuals,1,  p= probs)
      inds_observed <- inds_observed +sightings
    }
    observed_individuals[i] <- sum(inds_observed>0)
    observed_sightings[i] <- sum(inds_observed) 
  }
  res <- list('presence_observed' = observed_individuals, 'number_of_tracks_observed' = observed_sightings) 
}


# observe_individuals_along_transects <- function(distance_matrix, Ntransects){
#   # rows = individuals, 
#   # cols = transects
#   # observed_tracks counts how many presences observed per each transect:
#   Ntransects <- length(distance_matrix[1,]) 
#   Nlocations<- length(distance_matrix[,1]) 
#   observed_tracks <- numeric(0)
#   alpha <- 1 # decay in observation rate (w.r.t. distance to transect)
#   for (i in seq(1, Ntransects)){
#     distances_to_transect <- distance_matrix[, i]
#     obs_probs <- exp(-1*alpha*(distances_to_transect))
#     observed_tracks[i] <- sum(rbinom(1, Nlocations, p= obs_probs))
#   }
#   res <- list('presence_observed' = 1*(observed_tracks>0), 'number_of_tracks_observed' = observed_tracks) 
# }
# 
# 
# observe_individuals_along_transects <- function(distance_matrix, Ntransects){
#   # rows = individuals, 
#   # cols = transects
#   # observed_tracks counts how many presences observed per each transect:
#   Ntransects <- length(distance_matrix[1,]) 
#   Nlocations<- length(distance_matrix[,1]) 
#   observed_tracks <- numeric(0)
#   alpha <- 1 # decay in observation rate (w.r.t. distance to transect)
#   for (i in seq(1, Ntransects)){
#   distances_to_transect <- distance_matrix[, i]
#   obs_probs <- exp(-1*alpha*(distances_to_transect))
#   observed_tracks[i] <- sum(rbinom(1, Nlocations, p= obs_probs))
#   }
#   res <- list('presence_observed' = 1*(observed_tracks>0), 'number_of_tracks_observed' = observed_tracks) 
# }
#   
observe_individuals_intersecting_homerange_transects <- function(r, home_ranges, ObsTransects){
  Nindividuals <- length(home_ranges)
  Ntransects <- length(ObsTransects)
  habitat_intersects_with_transect <- numeric(0)
  home_range_cells <- unlist(home_ranges) 
  all_home_range_cells <- unique(home_range_cells)
  for (i in seq(1, Ntransects)){
    habitat_intersects_with_transect[i] <-0
    raster_cells <-  cellFromLine(r, ObsTransects[i])
    habitat_intersects_with_transect[i] <-  length(intersect(all_home_range_cells, raster_cells[[1]]))
  }
  observed_tracks <- rbinom(1, sum(habitat_intersects_with_transect), 1/3)
  res <- list('presence_observed' = observed_tracks)
  return(res)
}

observe_individuals_intersecting_homerange_squares <- function(r, home_ranges, ObsSquares){
  #xx <- rasterize(ObsTransects, r, fun='count')
  habitat_intersects_with_transect <- numeric(0)
  habitat_usage_within_cell <- numeric(0)
  Ntransects <- length(ObsSquares)
  Nindividuals <- length(home_ranges)
  individual_intersects_with_transect <- matrix(0, nrow = Ntransects, ncol = Nindividuals)
  number_of_individuals_observed  <- numeric(0)
  presence_observed  <- numeric(0)
  all_cells_in_home_ranges <- unique(unlist(home_ranges))
  
  for (i in seq(1, Ntransects)){
    habitat_intersects_with_transect[i] <-0
    habitat_usage_within_cell[i] <-0
    raster_cells <- cellFromPolygon(r, ObsSquares[i])
    #print(raster_cells)
    habitat_usage_within_cell[i] <- length(intersect(all_cells_in_home_ranges, raster_cells[[1]]))
    individuals_usage_of_cell <-  vapply(home_ranges, function(x) sum(is.element(x,raster_cells[[1]])), 0)
    presence_observed[i] <- rbinom(1,habitat_usage_within_cell[i], prob=1/3)>0
    number_of_individuals_observed[i] <- sum(rbinom(rep(1,Nindividuals),individuals_usage_of_cell, prob=1/3)>0)  
    } 
  
  res <- list('presence_observed' = presence_observed, 'number_of_individuals_observed' = number_of_individuals_observed, 
              'landscape_exploitation' =rbinom(Ntransects,habitat_usage_within_cell, prob=1/3)) 
  return(res)
}

plot_observations_intersecting_squares <- function(r, home_ranges, ObsSquares, ObsResult){
  rr <- r
  values(rr) <- 0
  rr[unlist(home_ranges)] <- 1
  sienna <-  brocolors("crayons")["Dandelion"]
  plot(rr,axes=FALSE, box=FALSE,  legend=FALSE, col = c('white', sienna))
  col1 <-  brocolors("crayons")["Jazzberry Jam"]
  col2 <-  brocolors("crayons")["Blue"]
  palette <- rep('white', length(ObsResult))
  palette[which(ObsResult==1)] <- col1
  palette[which(ObsResult==0)] <- col2
  print(ObsResult)
  print(palette)
  for (i in seq(1, length(ObsResult))){
  lines(ObsSquares[i], col = palette[i], add = T)}
}

survey_FS_population_polygons <- function(ObservationLocations, Population, border, Parameters){
  
  # OccupiedLocations: binary list of observed presences and absences (this should be compared with the data)
  # Females: how many females have their core areas in the location
  # Males: how many males have their core areas there. 
  # FemalesMissed: how many Females are left totally out of the study area 
  # MalesMissed: (same)
  
  NLoc <- length(ObservationLocations)
  OccupiedLocations <- rep(0, NLoc)
  Females <- rep(0, NLoc)
  Males <- rep(0, NLoc)
  
  allfemales <- unique(Population$FemaleID)
  allmales <- unique(Population$MaleID)
  
  femaleextinction <- Population$FemaleExtinction
  maleextinction <- Population$MaleExtinction
  
  Nfemales <- (femaleextinction!=1)*length(allfemales)
  Nmales <- (maleextinction!=1)*length(allmales)
  
  # if (Nmales>0){MaleAll <- SpatialPoints(Population$MaleAll)} else{MaleAll <- NULL}
  #  if (Nfemales>0){FemaleAll <- SpatialPoints(Population$FemaleAll)}else{FemaleAll <- NULL}
  
  FemaleAll <- Population$FemalePolygons
  MaleAll <- Population$MalePolygons
  
  observedfemales <- NULL
  observedmales <- NULL
  
  for (i in seq(1, NLoc)){
    
    #buf <- gBuffer(ObservationLocations[i], width = ObservationRadius)
    radius <- 150
    centroidy <- ObservationLocations@coords[i, 2]
    centroidx <- ObservationLocations@coords[i, 1] 
    
    yPlus <- centroidy +radius
    xPlus <- centroidx +radius
    yMinus <- centroidy-radius
    xMinus <- centroidx-radius
    
    square=rbind(cbind(xMinus,yPlus),  # NW corner
                 cbind(xPlus, yPlus),  # NE corner
                 cbind(xPlus,yMinus),  # SE corner
                 cbind(xMinus,yMinus), # SW corner
                 cbind(xMinus,yPlus))  # NW corner again - close ploygon
    
    s <- Polygon(square)
    buf <- SpatialPolygons(list(Polygons(list(s), 1)))
    
    # the indices for which the core area falls within the polygon
    if (Nmales>0){maleinds <-   !is.na(over(MaleAll, buf, returnList = F))}else{maleinds <- NULL}
    if (Nfemales>0){femaleinds <- !is.na(over(FemaleAll,buf,  returnList = F))}else{femaleinds <- NULL}
    
    OccupiedLocations[i] <- 1*((1*(sum(maleinds)>0) + 1*(sum(femaleinds)>0))>0) # is the location recorded in the data as occupied? 
    Females[i] <- sum(femaleinds)
    Males[i] <- sum(maleinds)
  }
  
  SurveyResults <-list( 'OccupiedLocations' = OccupiedLocations, 'Females' = Females, 'Males' = Males) 
  return(SurveyResults)
}

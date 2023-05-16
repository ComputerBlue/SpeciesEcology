# Each home range function has parameter max_diameter that defines how flexibly and how large the home ranmge may be

# write home-range functions that utilize the spread function, investigate the execution time. 

# fires <- spread(r, loci = locs, 
#                 0.5*r, 0, NULL, 1e8, 4, iterations = 2, id = TRUE, returnIndices = TRUE)
# fires
# 
# split_fires <- split(fires[, c('id', 'indices')], by = 'id', keep.by = FALSE, flatten = TRUE)
# ll <- lapply(lapply(split_fires, unlist), as.numeric)

# use mask when modeling fires <- spread(r, loci = locs, 
# 0.9*mask, 0, NULL, 1e8, 4, iterations = 20, id = TRUE)

# CHECK THE FUNCTIONS AND SPREAD PROBABILITIES! (and how the habitat influences them!)
# -> THE STOP RULES (ITERATIONS OR DIAMETER)

get_home_range_overlaps <- function(home_ranges){
  l <- home_ranges
  Nindividuals <- length(l)
  shared_territory_cells <- sapply(seq_len(length(l)), function(x) 
    sapply(seq_len(length(l)), function(y) length(intersect(unlist(l[x]), unlist(l[y])))))
  shared_cells <- numeric(0)
  for (i in seq(1, Nindividuals)){
    shared_cells[i] <- 1+sum(shared_territory_cells[i,])-shared_territory_cells[i,i]
  }
  return(shared_cells)
}

get_home_range_sizes <- function(home_ranges){
  home_range_sizes <- unlist(lapply(home_ranges, 'length'))
  return(home_range_sizes)
}



make_convex_polygon <- function(home_range, r){
  xys <- xyFromCell(r, home_range)
  dat <- chull(xys)
  coords <- xys[c(dat, dat[1]), ]
  a <- Polygons(list(Polygon(coords)), ID=1)
  return(a)
}

spread_random_overlapping <-function(locus, raster, spread_prob2, spread_prob){
  result <- spread(raster, loci = locus, spreadProb = spread_prob,  persistence = 0, directions= 4, iterations = 4, maxSize = 18, id = TRUE, returnIndices = TRUE, allowOverlap = 0, lowMemory = TRUE)
  return(result)
}


set_home_ranges_constant_size <- function(r, individual_cells, Nhabitat, Nindividuals,  max_diameter){
  
  home_ranges <- list()
  
  # to get overlapping, we use lapply, spread prob is evenly 0.75 throughout the landscape
  # stop rule is 3 iterations (max). 
  
  fires <- lapply(individual_cells, spread_random_overlapping,  r,  0.75, 1)
  ff <- sqrt(2)/2
  for (h in seq(1, length(individual_cells))){
    # sample Nhabitat cells cells and make a convex set out of them. 
    xys <- xyFromCell(r, fires[[h]]$indices)
    dat <- chull(xys)
    xys <- rbind(xys+c(ff, ff), xys+c(ff, -ff), xys+c(-ff, -ff), xys+c(-ff, ff))
    coords <- xys[c(dat, dat[1]), ]
    sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)))
    ee <- raster::extract(r, sp_poly,cellnumbers = TRUE)
    home_ranges[[h]] <- ee[[1]][,1]
  } 
  
  return(home_ranges) 
}




set_home_ranges_by_habitat <- function(r, individual_cells, Nhabitat, Nindividuals,  max_diameter){
  
  home_ranges <- list()
  #r[r==0] <- 0.5
  spread_prob_raster <- r
  spread_prob_raster[spread_prob_raster ==0] <- 0.65
  spread_prob_raster[spread_prob_raster ==1] <- 0.95
  # to get overlap, apply lapply. Spread rule is 0.75 to non-habitat and 1 to habitat
  # stop rule is stoprule (Nhabitat) or maxsize 30!
  
  fires <- lapply(individual_cells, spread_habitat_overlapping,  r,  spread_prob_raster, Nhabitat)
  ff <- sqrt(2)/2
  for (h in seq(1, length(individual_cells))){
    # sample Nhabitat cells cells and make a convex set out of them. 
    xys <- xyFromCell(r, fires[[h]]$indices)
    xys <- rbind(xys+c(ff, ff), xys+c(ff, -ff), xys+c(-ff, -ff), xys+c(-ff, ff))
    dat <- chull(xys)
    coords <- xys[c(dat, dat[1]), ]
    sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)))
    ee <- raster::extract(r, sp_poly,cellnumbers = TRUE)
    home_ranges[[h]] <- ee[[1]][,1]
  } 
  
  return(home_ranges)
}

spread_habitat_overlapping <-function(locus, raster, spread_prob_raster, Nhabitat){
  stopRule2 <- function(landscape) sum(landscape) > Nhabitat
  result <- spread(landscape = raster, loci = locus, spreadProb= spread_prob_raster, persistence=0, mask= NULL, directions=4, maxSize = 18, id = TRUE, returnIndices = TRUE, stopRule= stopRule2,stopRuleBehavior = "excludePixel", allowOverlap = 0,  lowMemory = TRUE)
  return(result)
}




permute_to_unoccupied_cell <- function(r, individual_cells){
  
  dup_ids <- which(duplicated(individual_cells))
  if (length(dup_ids)>0){
  dup_cells <- individual_cells[dup_ids]
  
  for (j in seq(1, length(dup_ids))){
  neighbours <- adjacent(r, dup_cells[j], directions=16)
  free_neighbours <- setdiff(neighbours[,2], individual_cells)
  if (length(free_neighbours)>0){
    individual_cells[dup_ids[j]] <-sample(free_neighbours,1) 
  }
  }
  }
  return(individual_cells)
}


set_home_ranges_no_overlap_even <- function(r, individual_cells, Nhabitat, Nindividuals,  max_diameter){
  
  # for some reason this might lose one individual?
  # need to check why someone disappears. Could be the split or lapply? 
  stopRule2 <- function(landscape) sum(landscape) > Nhabitat
  r2 <- r 
  r2[r2==0] <- 0.55
  r2[r2==1] <- 0.85
  r3 <- r 
  values(r3) <- 1
  individual_cells <- permute_to_unoccupied_cell(r, individual_cells)
  fires <- spread(landscape = r, loci = individual_cells, spreadProb= r2, persistence = 0, mask= NULL, maxSize =16, directions=4,  id = TRUE, stopRule= stopRule2, returnIndices = TRUE,stopRuleBehavior = "excludePixel", allowOverlap = 0, relativeSpreadProb = FALSE, lowMemory = TRUE)
  split_fires <- split(fires[, c('id', 'indices')], by = 'id', keep.by = FALSE, flatten = TRUE)
  home_ranges <- lapply(lapply(split_fires, unlist), as.numeric)
  poly_list <- list()
  ff <- sqrt(2)/2
  for (h in seq(1, length(home_ranges))){
  #sample Nhabitat cells cells and make a convex set out of them.
  xys <- xyFromCell(r, home_ranges[[h]])
  xys <- rbind(xys+c(ff, ff), xys+c(ff, -ff), xys+c(-ff, -ff), xys+c(-ff, ff))
  dat <- chull(xys)
  coords <- xys[c(dat, dat[1]), ]
  poly_list[[h]] <- Polygons(list(Polygon(coords)), ID=h)
  }
  # write a check that if someone is without a home-range, it is allocated with its focal cell (and immediate neighbours). 
  hh <- h
  dup_ids <- which(duplicated(individual_cells))
  if (length(dup_ids)>0){
    for (h in seq(1, length(dup_ids))){
    xys <- xyFromCell(r, individual_cells[dup_ids[h]])
    xys <- rbind(xys+c(ff, ff), xys+c(ff, -ff), xys+c(-ff, -ff), xys+c(-ff, ff))
    dat <- chull(xys)
    coords <- xys[c(dat, dat[1]), ]
    poly_list[[hh+h]] <- Polygons(list(Polygon(coords)), ID=hh+h)
    }
  }
  
  sp_poly <- SpatialPolygons(poly_list)
  ee <- raster::extract(r, sp_poly,cellnumbers = TRUE )

  home_ranges <- lapply(ee, '[', , 1)
  return(home_ranges)
}

spread_habitat_overlapping_core <-function(locus, raster, spread_prob_raster, Nhabitat, spread_prob_later){
  stopRule2 <- function(landscape) sum(landscape) > Nhabitat
  values(spread_prob_raster) <- 1
  result <- spread(landscape = raster, loci = locus, spreadProb= spread_prob_raster, spreadProbLater= spread_prob_later, persistence=0, mask= NULL, directions=4, maxSize = 25, id = TRUE, returnIndices = TRUE, stopRule= stopRule2, stopRuleBehavior = "excludePixel",allowOverlap = 0,lowMemory = TRUE)
  return(result)
}


set_home_ranges_by_non_overlapping_core_areas <- function(r, individual_cells, Nhabitat, Nindividuals,  max_diameter){
  
  home_ranges <- list()

  # take the first adjacent cells to each considered cell, take raster as output (this ensures that each individual will have the core area in its territory)
  
  individual_cells <- permute_to_unoccupied_cell(r, individual_cells)
  initial_spread <- spread(landscape = r, loci = individual_cells, spreadProb= 1, iterations = 1, persistence=0, mask= NULL, directions=4, maxSize = 18, id = TRUE, returnIndices =FALSE,  allowOverlap = 0,  lowMemory = TRUE)
  initial_spread[initial_spread>0] <- 1 

  # make the subsequent spread raster with habitat preference: 
  spread_prob_raster <- r
  spread_prob_raster[spread_prob_raster==0] <- 0.65
  spread_prob_raster[spread_prob_raster==1] <- 0.95
  
  # but omit the spread to someone else's core areas: 
  spread_prob_raster[initial_spread>0] <- 0 
  
  #values(initial_spread) <-1
  
  # before simulating the fires, permute to unoccupied cells, if possible (see how this function works): 

  fires <- lapply(individual_cells, spread_habitat_overlapping_core,  r,  initial_spread, Nhabitat, spread_prob_raster)
  ff <- sqrt(2)/2
  poly_list <- list()
  for (h in seq(1, length(individual_cells))){
    # sample Nhabitat cells cells and make a convex set out of them. 
    xys <- xyFromCell(r, fires[[h]]$indices)
    dat <- chull(xys)
    coords <- xys[c(dat, dat[1]), ]
    xys <- rbind(xys+c(ff, ff), xys+c(ff, -ff), xys+c(-ff, -ff), xys+c(-ff, ff))
    #sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)))
    poly_list[[h]] <- Polygons(list(Polygon(coords)), ID=h)
    # ee <- raster::extract(r, sp_poly,cellnumbers = TRUE)
    # home_ranges[[h]] <- ee[[1]][,1]
  } 

  sp_poly <- SpatialPolygons(poly_list)
  ee <- raster::extract(r, sp_poly,cellnumbers = TRUE )
  home_ranges <- lapply(ee, '[', , 1)
  
  hh <- h

  # # fires were not simulated for 
  # dup_ids <- which(duplicated(individual_cells))
  # if (length(dup_ids)>0){
  #   for (h in seq(1, length(dup_ids))){
  #     xys <- xyFromCell(r, individual_cells[dup_ids[h]])
  #     xys <- rbind(xys+c(ff, ff), xys+c(ff, -ff), xys+c(-ff, -ff), xys+c(-ff, ff))
  #     dat <- chull(xys)
  #     coords <- xys[c(dat, dat[1]), ]
  #     poly_list[[hh+h]] <- Polygons(list(Polygon(coords)), ID=hh+h)
  #   }
  # }
  # 
  return(home_ranges)
}




#########################

plot_home_ranges <- function(home_ranges, raster, Nindividuals, fitness){
  if (FALSE){
    raster <- initial_habitat
  }

    rbPal <- colorRampPalette(c('blue','red'))
    col_scale <- rbPal(10)[as.numeric(cut(fitness,breaks = 10))]
  
  maroon <- brocolors("crayons")["Asparagus"]
  image(raster,axes=FALSE, box=FALSE,  legend=FALSE, col = c('white', maroon))
  # lines(polygons)
  
  for (i in seq(1,Nindividuals)){
    hr <- home_ranges[[i]]
    rr <- raster
    values(rr) <- 0
    rr[hr] <- 1
    polygons <- gUnaryUnion((rasterToPolygons(rr, fun=function(x){x>0})))
    
    plot(polygons, add=TRUE)
  }

}

plot_home_ranges_nocolor <- function(home_ranges, raster, Nindividuals, hr_type, individual_cells){
  
  if (FALSE){
    raster <- initial_habitat
  }
  #rbPal <- colorRampPalette(c('blue','red'))
  #col_scale <- rbPal(10)[as.numeric(cut(fitness,breaks = 10))]
  
  pal <- c('white', 'olivedrab2')
  plot(as.factor(raster),axes=FALSE, box=FALSE,  legend=FALSE, col = pal)
  # lines(polygons)
  
  for (i in seq(1,Nindividuals)){
    hr <- home_ranges[[i]]
    rr <- raster
    values(rr) <- 0
    rr[hr] <- 1
    #poly <- rasterToPolygons(rr, fun=function(x){x>0})
    poly <- gUnaryUnion((rasterToPolygons(rr, fun=function(x){x>0})))
   # if (hr_type!='nonoverlapping_uneven'){
    poly <- gConvexHull(poly)
   # }
    plot(poly, add=TRUE, border = 'slateblue4')
    individual_cells_xy <- xyFromCell(raster, individual_cells)
    points(individual_cells_xy, col = 'slateblue4', pch = 19, cex = 0.8)
  }
  
}


plot_home_ranges_nocolor2 <- function(home_ranges, raster, Nindividuals){
  
  if (FALSE){
    raster <- initial_habitat
  }
  #rbPal <- colorRampPalette(c('blue','red'))
  #col_scale <- rbPal(10)[as.numeric(cut(fitness,breaks = 10))]
  seagreen <- brocolors("crayons")["Spring Green"]
  
  pal <- c('white', seagreen)
  plot(as.factor(raster),axes=FALSE, box=FALSE,  legend=FALSE, col = pal)
  # lines(polygons)
  
  for (i in seq(1,Nindividuals)){
    hr <- home_ranges[[i]]
    rr <- raster
    values(rr) <- 0
    rr[hr] <- 1
    #poly <- rasterToPolygons(rr, fun=function(x){x>0})
    poly <- gUnaryUnion((rasterToPolygons(rr, fun=function(x){x>0})))
    # if (hr_type!='nonoverlapping_uneven'){
    poly <- gConvexHull(poly)
    # }
    col1 <-  brocolors("crayons")["Midnight Blue"]
    plot(poly, add=TRUE, border = col1)
   
  }
  
}


plot_home_ranges_nocolor3 <- function(home_ranges, raster, Nindividuals, individual_cells){
  
  if (FALSE){
    raster <- initial_habitat
  }
  #rbPal <- colorRampPalette(c('blue','red'))
  #col_scale <- rbPal(10)[as.numeric(cut(fitness,breaks = 10))]
  seagreen <- brocolors("crayons")["Spring Green"]
  
  Nindividuals <- min(length(home_ranges), length(individual_cells))
  pal <- c('white', seagreen)
  plot(as.factor(raster),axes=FALSE, box=FALSE,  legend=FALSE, col = pal)
  # lines(polygons)
  
  for (i in seq(1,Nindividuals)){
    hr <- home_ranges[[i]]
    rr <- raster
    values(rr) <- 0
    rr[hr] <- 1
    #poly <- rasterToPolygons(rr, fun=function(x){x>0})
    poly <- gUnaryUnion((rasterToPolygons(rr, fun=function(x){x>0})))
    # if (hr_type!='nonoverlapping_uneven'){
    poly <- gConvexHull(poly)
    # }
    col1 <-  brocolors("crayons")["Midnight Blue"]
    plot(poly, add=TRUE, border = col1)
    xys <- xyFromCell(rr, individual_cells)
    points(xys[seq(1,Nindividuals),], add=TRUE, col = col1, pch = 19, cex = 0.5)
    
  }
  
}


plot_home_ranges_nocolor4 <- function(home_ranges, raster, Nindividuals, individual_cells){
  
  if (FALSE){
    raster <- initial_habitat
  }
  #rbPal <- colorRampPalette(c('blue','red'))
  #col_scale <- rbPal(10)[as.numeric(cut(fitness,breaks = 10))]
  seagreen <- brocolors("crayons")["Spring Green"]
  
  Nindividuals <- min(length(home_ranges), length(individual_cells))
  pal <- c('white', 'white')
  plot(as.factor(raster),axes=FALSE, box=FALSE,  legend=FALSE, col = pal)
  # lines(polygons)
  
  for (i in seq(1,Nindividuals)){
    hr <- home_ranges[[i]]
    rr <- raster
    values(rr) <- 0
    rr[hr] <- 1
    #poly <- rasterToPolygons(rr, fun=function(x){x>0})
    poly <- gUnaryUnion((rasterToPolygons(rr, fun=function(x){x>0})))
    # if (hr_type!='nonoverlapping_uneven'){
    poly <- gConvexHull(poly)
    # }
    col1 <-  brocolors("crayons")["Midnight Blue"]
    plot(poly, add=TRUE, border = col1)
    xys <- xyFromCell(rr, individual_cells)
    points(xys[seq(1,Nindividuals),], add=TRUE, col = col1, pch = 19, cex = 0.5)
    
  }
  
}


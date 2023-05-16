

simulate_passive_dispersal <- function(NB,dsts, individual_cells, habitat_cells){
  
  Nindividuals <- length(individual_cells)
  new_location_cells <- rep(0, Nindividuals)
  for (i in seq(1, Nindividuals)){
    sample_neighbour <- sample(NB[[individual_cells[i]]],1, prob = exp(-1*dsts[[individual_cells[i]]]/5))
    new_location_cells[i] <-  sample_neighbour
  }
  return(new_location_cells)
}

simulate_conspecific_dispersal <- function(NB, dsts, individual_cells, raster, NB_perception, dsts_perception, sum_habitat){
  
  Nindividuals <- length(individual_cells)
  new_location_cells <- rep(0, Nindividuals)
  values(raster) <- 0
  
  for (i in seq(1, Nindividuals)){
    neighbouring_cells <- NB[[individual_cells[i]]]
    sample_neighbour <- sample(neighbouring_cells,1, prob = exp(-1*dsts[[individual_cells[i]]]/5)*(population_density_induced_fitness(raster[neighbouring_cells])+0.01)*(sum_habitat[neighbouring_cells]+0.01))
    new_location_cells[i] <-  sample_neighbour
    
    raster[sample_neighbour] <- raster[sample_neighbour] +1
    raster[NB_perception[[sample_neighbour]]] <-raster[NB_perception[[sample_neighbour]]] +1/dsts_perception[[sample_neighbour]]
  }
  return(new_location_cells)
}

simulate_mid_density_driven_dispersal <- function(NB, dsts,individual_cells, raster, NB_perception, dsts_perception, sum_habitat){
  
  Nindividuals <- length(individual_cells)
  new_location_cells <- rep(0, Nindividuals)
  values(raster) <- 0
  
  for (i in seq(1, Nindividuals)){
    neighbouring_cells <- NB[[individual_cells[i]]]
    sample_neighbour <- sample(neighbouring_cells,1, prob = exp(-1*dsts[[individual_cells[i]]]/5)*(population_density_induced_fitness_gamme(raster[neighbouring_cells])+0.01)*(sum_habitat[neighbouring_cells]+0.01))
    new_location_cells[i] <-  sample_neighbour
    raster[sample_neighbour] <- raster[sample_neighbour] +1
    raster[NB_perception[[sample_neighbour]]] <-raster[NB_perception[[sample_neighbour]]] +1/dsts_perception[[sample_neighbour]]
  }
  return(new_location_cells)
}



simulate_avoidance_driven_dispersal <- function(NB, dsts,individual_cells, raster, NB_perception, dsts_perception, sum_habitat){
  
  Nindividuals <- length(individual_cells)
  new_location_cells <- rep(0, Nindividuals)
  values(raster) <- 0
  
  for (i in seq(1, Nindividuals)){
    neighbouring_cells <- NB[[individual_cells[i]]]
    sample_neighbour <- sample(neighbouring_cells,1, prob = exp(-1*dsts[[individual_cells[i]]]/5)*(inverse_density_induced_fitness(raster[neighbouring_cells])+0.01)*(sum_habitat[neighbouring_cells]+0.01))
    new_location_cells[i] <-  sample_neighbour
    raster[sample_neighbour] <- raster[sample_neighbour] +1
    raster[NB_perception[[sample_neighbour]]] <-raster[NB_perception[[sample_neighbour]]] +1/dsts_perception[[sample_neighbour]]
  }
  return(new_location_cells)
}


simulate_habitat_driven_dispersal <- function(NB, dsts,individual_cells, habitat_cells, sum_habitat){
  Nindividuals <- length(individual_cells)
  new_location_cells <- rep(0,  Nindividuals)
  for (i in seq(1, Nindividuals)){
    sample_neighbour <- sample(NB[[individual_cells[i]]],1, prob =exp(-1*dsts[[individual_cells[i]]]/5)* (sum_habitat[NB[[individual_cells[i]]]]+0.01))
    in_habitat <- (sum_habitat[sample_neighbour]>0)
    new_location_cells[i] <-  sample_neighbour
  }
  return(new_location_cells)
}






initialize_dispersal <- function(dispersal_type, initial_habitat, dispersal_distance, NB,xy_individuals, cell_ids_individuals, habitat_cells, sum_habitat, sum_individuals){
  
  xy_individuals<- sample_individuals_to_landscape(initial_habitat, Nindividuals)

  
  #compute the population densities to use in the dispersal
  NSimYears <- 10
  W <- matrix(1,3,3)
  
  for (y in seq(1, NSimYears)){
    if (y>1){
      xy_individuals <- xyFromCell(initial_habitat, new_locs)
    }  
    tab <- table(cellFromXY(initial_habitat, xy_individuals))
    pd <- raster(ncol=100, nrow=100, xmn=1, xmx=100, ymn=1, ymx=100)
    cell_ids_individuals <- cellFromXY(initial_habitat, xy_individuals)
    habitat_cells <-  cellFromXY(initial_habitat, cbind(initial_habitat_df$x[initial_habitat_df$habitat==1], initial_habitat_df$y[initial_habitat_df$habitat==1]))
    
    # raster giving counts of individuals in each cells
    values(pd) <- 0
    pd[as.numeric(names(tab))] <- tab
    sum_individuals <- rasterLocalSums(pd,W)
    sum_habitat <- rasterLocalSums(initial_habitat,W)
    habitat_cells <- which(values(initial_habitat)>0)
    if (y==1){
      new_locs <- simulate_species_dispersal('passive', initial_habitat, dispersal_distance, NB,xy_individuals, cell_ids_individuals, habitat_cells, sum_habitat, sum_individuals) 
    }
  }
  xy_individuals <- xyFromCell(initial_habitat, new_locs)
  res <- list('xy_individuals' = xy_individuals, 'cells_individuals' = new_locs )
  return(res)
}




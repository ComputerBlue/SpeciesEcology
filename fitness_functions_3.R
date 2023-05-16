
# functions that relate the population density, habitat & home ranges to reproductive/survival fitness weights.

# fitness is a >=0 number that can be related to reproductive output. 

# dev.off()
# par(mfrow = c(1,4))
# library(broman)
# col1 <- brocolors("crayons")["Pine Green"]
# col2 <- brocolors("crayons")["Mango Tango"]

# 
# png('FitnessComponents.png', height=  2.5, width = 9, units = 'in', res= 300)
# par(mfrow = c(1,4))
# x <- seq(0, 20, by= 0.2)
# y <- logistic_function_habitat(x)
# plot(x, y, 'l', bty="n", col = col1, lwd= 1.5, xlab= 'Habitat within home range', ylab = 'Fitness effect')
# 
# x <- seq(0, 3*18, by= 0.2)
# y <- inverse_density_induced_fitness(x)
# plot(x, y, 'l', bty="n",  col = col2, lwd= 1.5, xlab= 'Home range sharing', ylab = 'Fitness effect')
# 
# x <- seq(0, 4*18, by= 0.2)
# y <- population_density_induced_fitness(x)
# plot(x, y, 'l', bty="n",  col = col2, lwd= 1.5, xlab= 'Home range sharing', ylab = 'Fitness effect')
# 
# x <- seq(0, 4*18, by= 0.2)
# y <- population_density_induced_fitness_gamma(x)
# plot(x, y, 'l', bty="n", col = col2,lwd= 1.5, xlab= 'Home range sharing', ylab = 'Fitness effect')
# dev.off()


population_density_induced_fitness_gamma <- function(overlaps){
  # for each individual, number of times a cell is shared with another individual
  #degree_of_shared_territory <- numeric(0)
  #l <- home_ranges
  #shared_territory_cells <- sapply(seq_len(length(l)), function(x) 
  #sapply(seq_len(length(l)), function(y) length(intersect(unlist(l[x]), unlist(l[y])))))
  L <- 1
  k <- 0.05
  x_0 <- 2
  
  # shared_cells <- numeric(0)
  # Nindividuals <- length(overlaps)
  # 
  # for (i in seq(1, Nindividuals)){
  #   shared_cells[i] <-logistic_function_inv(overlaps[i], L, k, x_0)
  # }
  #individual_fitness <- sapply(overlaps, logistic_function, L, k, x_0, simplify = TRUE)
  #individual_fitness <- logistic_function(overlaps, L, k, x_0)
  individual_fitness <-  dgamma(overlaps, shape=x_0, rate = k,  log = FALSE)
  
  return(individual_fitness)
}

population_density_induced_fitness <- function(overlaps){
  # for each individual, number of times a cell is shared with another individual
  #degree_of_shared_territory <- numeric(0)
  #l <- home_ranges
  #shared_territory_cells <- sapply(seq_len(length(l)), function(x) 
  #sapply(seq_len(length(l)), function(y) length(intersect(unlist(l[x]), unlist(l[y])))))
  L <- 1
  k <- 0.2
  x_0 <- 18
  
  # shared_cells <- numeric(0)
  # Nindividuals <- length(overlaps)
  # 
  # for (i in seq(1, Nindividuals)){
  #   shared_cells[i] <-logistic_function_inv(overlaps[i], L, k, x_0)
  # }
  #individual_fitness <- sapply(overlaps, logistic_function, L, k, x_0, simplify = TRUE)
  individual_fitness <- logistic_function(overlaps, L, k, x_0)
  #individual_fitness <-  dgamma(overlaps, shape=x_0, rate = k,  log = FALSE)
  
  return(individual_fitness)
}

inverse_density_induced_fitness <- function(overlaps){
  # for each individual, number of times a cell is shared with another individual
  shared_cells <- numeric(0)
  L <- 1
  k <- 0.2
  x_0 <- 18
  Nindividuals <- length(overlaps)
  
  # for (i in seq(1, Nindividuals)){
  #   shared_cells[i] <-logistic_function_inv(overlaps[i], L, k, x_0)
  # }
  #individual_fitness <- unlist(lapply(overlaps, logistic_function_inv, L, k, x_0))
  #print(shared_cells)
  
  individual_fitness <- logistic_function_inv(overlaps, L, k, x_0)
  return(individual_fitness)
}

logistic_function_habitat <- function(value){
  L <- 1
  k <- 0.5
  x_0 <- 3
  output <- L/(1+ exp(-k*(value-x_0)))
  return(output)
}

calculate_cell_usage <- function(home_ranges, r){
  #gives how many individuals are using that cell: 
  values(r) <- 0
  for (h in seq(1, length(home_ranges))){
    r[home_ranges[[h]]] <- r[home_ranges[[h]]]+1
  }
  return(r)
}

neutral_fitness <- function(Nindividuals){
  individual_fitness <- rep(1, Nindividuals)
  return(individual_fitness) 
}

habitat_induced_fitness <- function(r, home_ranges, Nindividuals){
 individual_fitness <- numeric(0)
 for(i in seq(1, Nindividuals)){
 individual_fitness[i] <- sum(r[home_ranges[[i]]], na.rm=T)  # need to check why there is NA's
 }
 individual_fitness <- logistic_function_habitat(individual_fitness)
 return(individual_fitness) 
}

hab_unit_area <- function(home_range, r){
individual_fitness<- (sum(r[home_range], na.rm=T))/(1+length(home_range))  # need to check why there is NA's
return(individual_fitness)  
}


habitat_per_unit_area_fitness <- function(r, home_ranges, Nindividuals){
  individual_fitness <- sapply(home_ranges, hab_unit_area, r, simplify = TRUE)
  return(individual_fitness) 
}


hab_consumer <- function(home_range, habitat, cell_usage){
  individual_fitness<- sum(habitat[home_range] %/% cell_usage[home_range])  # need to check why there is NA's
  return(individual_fitness)  
}


habitat_per_consumer<- function(r, home_ranges, Nindividuals){
  cell_usage <- calculate_cell_usage(home_ranges, r)
  individual_fitness <- sapply(home_ranges, hab_consumer, r, cell_usage, simplify = TRUE)
  return(individual_fitness) 
}

habitat_per_consumer_fitness <- function(r, home_ranges, Nindividuals){
  habitat_per_consumer_v <- habitat_per_consumer(r, home_ranges, Nindividuals)
  individual_fitness <- logistic_function_habitat(habitat_per_consumer_v)
  return(individual_fitness) 
}


logistic_function <- function(value, L, k, x_0){
  output <- L/(1+ exp(-k*(value-x_0)))
  return(output)
}

logistic_function_inv <- function(value, L, k, x_0){
  output <- L-L/(1+ exp(-k*(value-x_0)))
  return(output)
}



fitness_induced_reproduction <- function(individual_cells, individual_fitness, Nindividuals){
  Nindividuals <- length(individual_cells)
  modifier <- runif(1, min =0.8, max = 1.2) # the effect of the year
  #print(length(individual_fitness))
  #print(Nindividuals)
  probs <- individual_fitness*runif(Nindividuals, min =0.95, max = 1.05) # the effect of the year
  NewPopsize <- min(floor(modifier*Nindividuals), 1.2*Nindividuals)
  NewPopsize <-  max(NewPopsize, 0.8*Nindividuals)
  offspring_cells <- sample(individual_cells, NewPopsize, prob = probs, replace=TRUE)
  #print(length(offspring_cells))
  return(offspring_cells)
}

fitness_induced_reproduction_reference <- function(individual_cells, individual_fitness, Nindividuals, reference_fitness){

  individual_fitness[individual_fitness>(1.5*reference_fitness)] <- 1.5*reference_fitness
  individual_fitness[individual_fitness<(0.5*reference_fitness)] <- 0.5*reference_fitness
  modifier <-  runif(1, min =0.95, max = 1.05)
  expected_offspring <- modifier*floor(sum(individual_fitness/reference_fitness))
  offspring_cells <- sample(individual_cells, expected_offspring, prob = individual_fitness, replace=TRUE)
  return(offspring_cells)
  
}

plot_home_ranges_fitness <- function(home_ranges, raster, Nindividuals, fitness){
  
  if (FALSE){
    raster <- initial_habitat
  }
  
  
  
  plot(raster,axes=FALSE, box=FALSE)
  # lines(polygons)
  
  for (i in seq(1,Nindividuals)){
    hr <- home_ranges[[i]]
    rr <- raster
    values(rr) <- 0
    rr[hr] <- 1
    polygons <- gUnaryUnion((rasterToPolygons(rr, fun=function(x){x>0})))
    plot(polygons, A)
  }
  
}



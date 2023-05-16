#######################

# individual functions:

simulate_observations_environmental_change <- function(initial_habitat_raster, habitat_change_df, individual_cells, Nyears, plotting, hr, fitness, dispersal, NB,dsts, max_diameter, Nhabitat, Wr, obs_list, fitness_assumption,NB_perception, dsts_perception){
  

  N_locations <-1 # (these refer to individuals' behavior)
  habitat_dependency  <- FALSE  #check the exact definition of this parameter. 
  

  if (fitness_assumption == 'changing'){
    result <- simulate_time_series_environmental_change(initial_habitat_raster, habitat_change_df, individual_cells, hr, fitness, dispersal, NB, dsts, max_diameter, Nhabitat, Wr)
  } else {
    result <- simulate_time_series_environmental_change_stable_fitness(initial_habitat_raster, habitat_change_df, individual_cells, hr, fitness, dispersal, NB, dsts, max_diameter, Nhabitat, Wr, NB_perception, dsts_perception) 
  }
  
  
  raster_end <- generate_raster_from_cells(initial_habitat_raster, habitat_change_df, 5)
  Nindividuals_end <- length(result$individual_cells)
  loc_xy <- set_individual_locations(raster_end, result$home_ranges, Nhabitat, Nindividuals_end, always_in_habitat, N_locations, habitat_dependency)
  
  # compute summary statistics for the spatial distribution of the individuals
  
  ppp <- as.ppp(loc_xy$individual_locations_xy, c(0,100,0,100))
  nn1 <- nndist(ppp, k=1)
  observation_probability <- 1
  
  # o1:
  
  obs_within_squares_o1 <- allocate_locations_within_squares(loc_xy$individual_locations_xy, obs_list$obs_cells5_o1 )
  res_squares_intersecting_o1 <- observe_individuals_intersecting_homerange_squares(raster_end, result$home_ranges, obs_list$obs_cells5_o1)
  Nsquares  <- length(obs_list$obs_cells5_o1)
  Ntransects  <- length(obs_list$obs_transects5_o1)
  res_o1 <- observe_individuals_within_squares(obs_within_squares_o1$locations_per_square, observation_probability, Nsquares)
  res_transects_o1 <- observe_individuals_along_transects_2(obs_list$obs_transects5_o1, loc_xy$individual_locations_xy)
  res_transects_intersecting_o1 <- observe_individuals_intersecting_homerange_transects(raster_end, result$home_ranges, obs_list$obs_transects5_o1)
  
  
  # o2:
  
  obs_within_squares_o2 <- allocate_locations_within_squares(loc_xy$individual_locations_xy, obs_list$obs_cells5_o2)
  res_squares_intersecting_o2 <- observe_individuals_intersecting_homerange_squares(raster_end, result$home_ranges, obs_list$obs_cells5_o2)
  Nsquares  <- length(obs_list$obs_cells5_o2)
  Ntransects  <- length(obs_list$obs_transects5_o2)
  res_o2 <- observe_individuals_within_squares(obs_within_squares_o2$locations_per_square, observation_probability, Nsquares)

  

  
  # o4:
  
  obs_within_squares_o4 <- allocate_locations_within_squares(loc_xy$individual_locations_xy, obs_list$obs_cells5_o4 )
  res_squares_intersecting_o4 <- observe_individuals_intersecting_homerange_squares(raster_end, result$home_ranges, obs_list$obs_cells5_o4)
  Nsquares  <- length(obs_list$obs_cells5_o4)
  Ntransects  <- length(obs_list$obs_transects5_o4)
  res_o4 <- observe_individuals_within_squares(obs_within_squares_o4$locations_per_square, observation_probability, Nsquares)

 
  # observations at time point 10: 
  
  raster_end10 <- generate_raster_from_cells(initial_habitat_raster, habitat_change_df, 10)
  Nindividuals_end10 <- length(result$individual_cells10)
  loc_xy10 <- set_individual_locations(raster_end10, result$home_ranges10, Nhabitat, Nindividuals_end10, always_in_habitat, N_locations, habitat_dependency)
  ppp <- as.ppp(loc_xy10$individual_locations_xy, c(0,100,0,100))

  nn110 <- nndist(ppp, k=1)

  # o1:
  
  obs_within_squares_o1_10 <- allocate_locations_within_squares(loc_xy10$individual_locations_xy, obs_list$obs_cells10_o1 )
  res_squares_intersecting_o1_10 <- observe_individuals_intersecting_homerange_squares(raster_end10, result$home_ranges10, obs_list$obs_cells10_o1)
  Nsquares  <- length(obs_list$obs_cells10_o1)
  Ntransects  <- length(obs_list$obs_transects10_o1)
  res_o1_10 <- observe_individuals_within_squares(obs_within_squares_o1_10$locations_per_square, observation_probability, Nsquares)
  res_transects_o1_10 <- observe_individuals_along_transects_2(obs_list$obs_transects10_o1, loc_xy10$individual_locations_xy)
  res_transects_intersecting_o1_10 <- observe_individuals_intersecting_homerange_transects(raster_end10, result$home_ranges10, obs_list$obs_transects10_o1)

  # o2:
  
  obs_within_squares_o2_10 <- allocate_locations_within_squares(loc_xy10$individual_locations_xy, obs_list$obs_cells10_o2 )
  res_squares_intersecting_o2_10 <- observe_individuals_intersecting_homerange_squares(raster_end10, result$home_ranges10, obs_list$obs_cells10_o2)
  Nsquares  <- length(obs_list$obs_cells10_o2)
  Ntransects  <- length(obs_list$obs_transects10_o2)
  res_o2_10 <- observe_individuals_within_squares(obs_within_squares_o2_10$locations_per_square, observation_probability, Nsquares)

  

  
  # o4:
  obs_within_squares_o4_10 <- allocate_locations_within_squares(loc_xy10$individual_locations_xy, obs_list$obs_cells10_o4 )
  res_squares_intersecting_o4_10 <- observe_individuals_intersecting_homerange_squares(raster_end10, result$home_ranges10, obs_list$obs_cells10_o4)
  Nsquares  <- length(obs_list$obs_cells10_o4)
  Ntransects  <- length(obs_list$obs_transects10_o4)
  res_o4_10 <- observe_individuals_within_squares(obs_within_squares_o4_10$locations_per_square, observation_probability, Nsquares)

  

  
  res <-list('pop_size' = length(loc_xy$individual_locations_xy[,1]), 
             'pop_size10' = length(loc_xy10$individual_locations_xy[,1]),
             'home_ranges' = result$home_ranges, 
             'Nsquares' = Nsquares, 
             'home_ranges10' = result$home_ranges10,
             'fitness10'= result$fitness10,
             
             'presence_observed_o1' = res_o1$presence_observed, 
             'number_of_tracks_observed_o1' = res_o1$number_of_tracks_observed, 
             'presence_observed_transects_o1' = res_transects_o1$presence_observed, 
             'number_of_tracks_observed_transects_o1' = res_transects_o1$number_of_tracks_observed, 
             'presence_observed_intersecting_squares_o1' = res_squares_intersecting_o1$presence_observed,
             'number_individuals_observed_intersecting_squares_o1' = res_squares_intersecting_o1$number_of_individuals_observed,
             'landscape_exploitation_o1' = res_squares_intersecting_o1$landscape_exploitation,
             'presence_observed_intersecting_transects_o1' = res_transects_intersecting_o1$presence_observed, 
             
             
             'presence_observed_o2' = res_o2$presence_observed, 
             'number_of_tracks_observed_o2' = res_o2$number_of_tracks_observed, 
             'presence_observed_intersecting_squares_o2' = res_squares_intersecting_o2$presence_observed,
             'number_individuals_observed_intersecting_squares_o2' = res_squares_intersecting_o2$number_of_individuals_observed,
             'landscape_exploitation_o2' = res_squares_intersecting_o2$landscape_exploitation,
             

             'presence_observed_o4' = res_o4$presence_observed, 
             'number_of_tracks_observed_o4' = res_o4$number_of_tracks_observed, 
             'presence_observed_intersecting_squares_o4' = res_squares_intersecting_o4$presence_observed,
             'number_individuals_observed_intersecting_squares_o4' = res_squares_intersecting_o4$number_of_individuals_observed,
             'landscape_exploitation_o4' = res_squares_intersecting_o4$landscape_exploitation,
  
             'presence_observed_o1_10' = res_o1_10$presence_observed, 
             'number_of_tracks_observed_o1_10' = res_o1_10$number_of_tracks_observed, 
             'presence_observed_transects_o1_10' = res_transects_o1_10$presence_observed, 
             'number_of_tracks_observed_transects_o1_10' = res_transects_o1_10$number_of_tracks_observed, 
             'presence_observed_intersecting_squares_o1_10' = res_squares_intersecting_o1_10$presence_observed,
             'number_individuals_observed_intersecting_squares_o1_10' = res_squares_intersecting_o1_10$number_of_individuals_observed,
             'landscape_exploitation_o1_10' = res_squares_intersecting_o1_10$landscape_exploitation,
             'presence_observed_intersecting_transects_o1_10' = (res_transects_intersecting_o1_10$presence_observed), 
             
             
             'presence_observed_o2_10' = res_o2_10$presence_observed, 
             'number_of_tracks_observed_o2_10' = res_o2_10$number_of_tracks_observed, 
             'presence_observed_intersecting_squares_o2_10' = res_squares_intersecting_o2_10$presence_observed,
             'number_individuals_observed_intersecting_squares_o2_10' = res_squares_intersecting_o2_10$number_of_individuals_observed,
             'landscape_exploitation_o2_10' = res_squares_intersecting_o2_10$landscape_exploitation,
             

             'presence_observed_o4_10' = res_o4_10$presence_observed, 
             'number_of_tracks_observed_o4_10' = res_o4_10$number_of_tracks_observed, 
             'presence_observed_intersecting_squares_o4_10' = res_squares_intersecting_o4_10$presence_observed,
             'number_individuals_observed_intersecting_squares_o4_10' = res_squares_intersecting_o4_10$number_of_individuals_observed,
             'landscape_exploitation_o4_10' = res_squares_intersecting_o4_10$landscape_exploitation,
             
             'nn1'= nn1,
             'nn110'= nn110

  )
  return(res)
}

simulate_time_series_environmental_change_stable_fitness <- function(initial_habitat_raster, habitat_change_df, individual_cells, hr, fitness, dispersal, NB, dsts, max_diameter, Nhabitat, Wr, NB_perception, dsts_perception){
  
  
  
  HabitatPerceptionW <- matrix(c( 0, 0,  0.25, 0.25, 0.25,  0, 0,  
                                 0,  0.5,  0.5, 0.5,  0.5,  0.5, 0,  
                                 0.25,  0.5,  1, 1,  1,  0.5, 0.25, 
                                 0.25,  0.5,  1, 1,  1,  0.5, 0.25,  
                                 0.25,  0.5,  1, 1,  1,  0.5, 0.25, 
                                 0,  0.5,  0.5, 0.5,  0.5,  0.25, 0,  
                                 0, 0, 0.25 , 0.25,  0.25,  0, 0), ncol=7, byrow=TRUE)
  
 
  W <- Wr
  Nyears <- 10
  Nindividuals <- length(individual_cells)

  for (y in seq(1, 10)){
    
    if (y==5){
      individual_cells_5 <- individual_cells
    }
    
    if (y==10){
      individual_cells_10 <- individual_cells
    }
    
    Nimmigration <- sample(seq(1, 10), 1)
    habitat_raster <- generate_raster_from_cells(initial_habitat_raster, habitat_change_df, y)
    
   

# 1 .home ranges -------------------------------------------------------------

    Nindividuals <- length(individual_cells)
    
    if (hr == 'constant_size'){ 
      home_ranges1 <- set_home_ranges_constant_size(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    if (hr == 'habitat'){ 
      home_ranges1 <- set_home_ranges_by_habitat(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    if (hr == 'nonoverlapping_habitat'){ 
      home_ranges1 <- set_home_ranges_by_habitat_no_overlap(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    if (hr == 'nonoverlapping_even'){ 
      home_ranges1 <- set_home_ranges_no_overlap_even(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    if (hr == 'nonoverlapping_core'){ 
      home_ranges1 <- set_home_ranges_by_non_overlapping_core_areas(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    if (hr == 'nonoverlapping_core'){ 
      home_ranges1 <- set_home_ranges_by_non_overlapping_core_areas(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    

    #######################
    

# 2 .fitness -----------------------------------------------------------------

    
    if (fitness=='neutral'){
      fitness1 <- neutral_fitness(length(home_ranges1))
      }  
      
    if (fitness=='habitat_density'){
      fitness0 <- habitat_induced_fitness(habitat_raster, home_ranges1, length(home_ranges1))
      overlaps <- get_home_range_overlaps(home_ranges1)
      fitness2 <- population_density_induced_fitness(overlaps)
      fitness1 <- fitness0*fitness2
      rm(overlaps)
    }
    
    if (fitness=='habitat_mid_density'){
      fitness0 <- habitat_induced_fitness(habitat_raster, home_ranges1, length(home_ranges1))
      overlaps <- get_home_range_overlaps(home_ranges1)
      fitness2 <- population_density_induced_fitness_gamma(overlaps)
      fitness1 <- fitness0*fitness2
      rm(overlaps)
    }
    
    if (fitness=='habitat_inverse_density'){
      fitness0 <- habitat_induced_fitness(habitat_raster, home_ranges1, length(home_ranges1))
      overlaps <- get_home_range_overlaps(home_ranges1)
      fitness2 <- inverse_density_induced_fitness(overlaps)
      fitness1 <- fitness0*fitness2
      rm(overlaps)
    }
    
    if (fitness=='habitat_per_consumer_density'){
      fitness0 <- habitat_per_consumer_fitness(habitat_raster, home_ranges1, length(home_ranges1))
      overlaps <- get_home_range_overlaps(home_ranges1)
      fitness2 <- population_density_induced_fitness(overlaps)
      fitness1 <- fitness0*fitness2
      rm(overlaps)
    }
    
    if (fitness=='habitat_per_consumer_mid_density'){
      fitness0 <- habitat_per_consumer_fitness(habitat_raster, home_ranges1, length(home_ranges1))
      overlaps <- get_home_range_overlaps(home_ranges1)
      fitness2 <- population_density_induced_fitness_gamma(overlaps)
      fitness1 <- fitness0*fitness2
      rm(overlaps)
    }
    
    if (fitness=='habitat_per_consumer_inverse_density'){
      fitness0 <- habitat_per_consumer_fitness(habitat_raster, home_ranges1, length(home_ranges1))
      overlaps <- get_home_range_overlaps(home_ranges1)
      fitness2 <- inverse_density_induced_fitness(overlaps)
      fitness1 <- fitness0*fitness2
      rm(overlaps)
    }
    
    
   ########################
    
    # ( reproduction ) 
    
    fitness1[is.na(fitness1)] <- 0
    individual_cells <- fitness_induced_reproduction(individual_cells, fitness1, Nindividuals)
    
    # ( immigration ) 
    
    xy_immigration<- sample_individuals_to_landscape(initial_habitat, Nimmigration)
    individual_cells_immigrants <-cellFromXY(initial_habitat, xy_immigration)
    individual_cells <- c(individual_cells, individual_cells_immigrants)
    
    habitat_df <- habitat_change_df[which(habitat_change_df$year==y), ]
    
    if (dispersal!='passive'){
    sum_habitat <- rasterLocalSums(habitat_raster, HabitatPerceptionW)
    sum_habitat <- logistic_function_habitat(sum_habitat)
    habitat_cells <-  cellFromXY(initial_habitat, cbind(habitat_df$x[habitat_df$habitat==1], habitat_df$y[habitat_df$habitat==1]))
    }
    
# 4. simulate dispersal!-------------------------------------------------------------------------

    
    if (y < 10){ # why? 
    if (dispersal == 'passive'){ 
      individual_cells <- simulate_passive_dispersal(NB, dsts, individual_cells, habitat_cells)
    }
    
    if (dispersal == 'conspecific'){ 
      individual_cells <- simulate_conspecific_dispersal(NB, dsts, individual_cells, habitat_raster, NB_perception, dsts_perception, sum_habitat)
    }
    if (dispersal == 'mid_density'){ 
        individual_cells <- simulate_conspecific_dispersal(NB, dsts, individual_cells, habitat_raster, NB_perception, dsts_perception, sum_habitat)
    }
    
    if (dispersal == 'habitat_driven'){ 
      individual_cells <- simulate_habitat_driven_dispersal(NB, dsts, individual_cells, habitat_cells, sum_habitat)
      rm(sum_habitat)
      }

    if (dispersal == 'conspecific_avoidance_driven'){ 
      individual_cells <- simulate_avoidance_driven_dispersal(NB, dsts,individual_cells, sum_habitat, NB_perception, dsts_perception, sum_habitat)
      rm(sum_habitat)
      }
    }
   

    
    if (y ==5){
      individual_cells5 <- individual_cells
      home_ranges5 <-home_ranges1 
    }
    
    #print('Dispersal occurred')
    #print(length(individual_cells))
    
    
  }
  
  res <- list('individual_cells' = individual_cells_5, 'home_ranges' = home_ranges5,
              'individual_cells10' = individual_cells_10, 'home_ranges10' = home_ranges1, 'fitness10' = fitness1)
  return(res)
}

simulate_time_series_environmental_change <- function(initial_habitat_raster, habitat_change_df, individual_cells, hr, fitness, dispersal, NB,dsts, max_diameter, Nhabitat, Wr){
  
  W <- Wr
  Nyears <- 5
  #Nhabitat <- 6 
  
  Nindividuals <- length(individual_cells)
  plotting_home_ranges <- TRUE
  
  for (y in seq(1, 5)){
    
    Nimmigration <- sample(seq(1, 10), 1)
    habitat_raster <- generate_raster_from_cells(initial_habitat_raster, habitat_change_df, y)
   
    if (y==5){
      individual_cells_5 <- individual_cells
    }
    
    # 1 set the home ranges
    
    if (hr == 'constant_size'){ 
      home_ranges1 <- set_home_ranges_constant_size(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    if (hr == 'habitat'){ 
    home_ranges1 <- set_home_ranges_by_habitat(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    if (hr == 'nonoverlapping_habitat'){ 
    home_ranges1 <- set_home_ranges_by_habitat_no_overlap(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    if (hr == 'nonoverlapping_uneven'){ 
    home_ranges1 <- set_home_ranges_no_overlap_uneven(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    if (hr == 'nonoverlapping_core'){ 
      
      home_ranges1 <- set_home_ranges_by_non_overlapping_core_areas(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    if (hr == 'nonoverlapping_even'){ 
      home_ranges1 <- set_home_ranges_no_overlap_even(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    # 2. set the fitness
    
    if (fitness=='habitat_induced'){
    fitness1 <- habitat_induced_fitness(habitat_raster, home_ranges1, length(home_ranges1))
    }
    
    if (fitness=='density_induced'){
      overlaps <- get_home_range_overlaps(home_ranges1)
      #sum_individuals <- get_individual_densities(habitat_raster, xy_individuals, Wr)
      fitness1 <- population_density_induced_fitness(overlaps)
    }
    
    if (fitness=='inverse_density'){
      overlaps <- get_home_range_overlaps(home_ranges1)
      #sum_individuals <- get_individual_densities(habitat_raster, xy_individuals, Wr)
      fitness1 <- inverse_density_induced_fitness(overlaps)
    }
    
    if (fitness=='neutral'){
   # print(home_ranges1)
    fitness1 <- neutral_fitness(length(home_ranges1))}  
    
   # hist(fitness1)
   # plot(sum_individuals, main = sum(values(sum_individuals)))
   #  points(xy_individuals, cex = fitness1)
    if (y==5){
      # this is used for letting the population to grow or decay according to realized fitnesses
      reference_fitness <- mean(fitness1)
    }
    
    #print('Fitness defined')
    
    # 3. let individuals reproduce

    reproduction <- fitness_induced_reproduction(individual_cells, fitness1, Nindividuals)
    #print('Reproduction occurred')
    
    individual_cells <- reproduction
    xy_immigration<- sample_individuals_to_landscape(initial_habitat, Nimmigration)
    individual_cells_immigrants <-cellFromXY(initial_habitat, xy_immigration)
    #print("Immigration occurred!")
    individual_cells <- c(individual_cells, individual_cells_immigrants)
    
    habitat_df <- habitat_change_df[which(habitat_change_df$year==y), ]
    sum_habitat <- rasterLocalSums(habitat_raster, Wr)
    habitat_cells <-  cellFromXY(initial_habitat, cbind(habitat_df$x[habitat_df$habitat==1], habitat_df$y[habitat_df$habitat==1]))
    
    #print("Something occurred!")
    # 4. simulate dispersal!
    # carefully define whether only newborns, or all individuals disperse (and how this affects the different dispersal rules)
    if (dispersal == 'passive'){ 
      individual_cells <- simulate_passive_dispersal(NB, dsts,individual_cells, habitat_cells)
    }
    
    if (dispersal == 'conspecific'){ 
      #xy_individuals <- xyFromCell(habitat_raster, individual_cells)
      #sum_individuals <- get_individual_densities(habitat_raster, xy_individuals, Wr)
      individual_cells <- simulate_conspecific_dispersal(NB, dsts, individual_cells,  habitat_raster)
    }
    
    if (dispersal == 'habitat_driven'){ 
      individual_cells <- simulate_habitat_driven_dispersal(NB, dsts, individual_cells, habitat_cells, sum_habitat)
    }
    
    # if (dispersal == 'habitat_and_conspecific_driven'){  
    #   #xy_individuals <- xyFromCell(habitat_raster, individual_cells)
    #   #sum_individuals <- get_individual_densities(habitat_raster, xy_individuals, Wr)
    #   individual_cells <- simulate_habitat_conspecific_dispersal(NB, dsts, individual_cells, habitat_cells, sum_individuals, sum_habitat)
    # }
    # 
    # if (dispersal == 'habitat_conspecific_avoidance_driven'){ 
    #   xy_individuals <- xyFromCell(habitat_raster, individual_cells)
    #   sum_individuals <- get_individual_densities(habitat_raster, xy_individuals, Wr)
    #   individual_cells <- simulate_habitat_avoidance_driven_dispersal(NB, dsts,individual_cells, habitat_cells, sum_individuals, sum_habitat)
    # }
    
    if (dispersal == 'conspecific_avoidance_driven'){ 
      #print(individual_cells)
      #print(habitat_raster)
      #xy_individuals <- xyFromCell(habitat_raster, individual_cells)
      #print(xy_individuals)
      #sum_individuals <- get_individual_densities(habitat_raster, xy_individuals, Wr)
      individual_cells <- simulate_avoidance_driven_dispersal(NB, dsts,individual_cells, sum_habitat)
    }
    
    Nindividuals <- length(individual_cells)
    #print(Nindividuals)
  }
 
  ##print("Dispersal occurred")
  # 
  # # (next step is for storing the intermediate home ranges)
  # 
  # if (hr == 'constant_size'){ 
  #   home_ranges1 <- set_home_ranges_constant_size(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
  # }
  # 
  # if (hr == 'habitat'){ 
  #   home_ranges1 <- set_home_ranges_by_habitat(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
  # }
  # if (hr == 'nonoverlapping_habitat'){ 
  #   home_ranges1 <- set_home_ranges_by_habitat_no_overlap(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
  #   }
  # if (hr == 'nonoverlapping_uneven'){ 
  #   home_ranges1 <- set_home_ranges_no_overlap_uneven(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
  # }
  # if (hr == 'nonoverlapping_even'){ 
  #   home_ranges1 <- set_home_ranges_no_overlap_even(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
  # }

  individual_cells5 <- individual_cells
  home_ranges5 <-home_ranges1 
  
  for (y in seq(6, 10)){
    
    habitat_raster <- generate_raster_from_cells(initial_habitat_raster, habitat_change_df, y)
    
    if (y==10){
      individual_cells_10 <- individual_cells
    }
    
    
    if (hr == 'constant_size'){ 
      home_ranges1 <- set_home_ranges_constant_size(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
      # some NA's appear here. Is it due to no habitat there? then perhaps the individuals die. 
    }
    if (hr == 'habitat'){ 
      home_ranges1 <- set_home_ranges_by_habitat(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
      # some NA's appear here. Is it due to no habitat there? then perhaps the individuals die. 
    }
    if (hr == 'nonoverlapping_habitat'){ 
      home_ranges1 <- set_home_ranges_by_habitat_no_overlap(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)}
    
    if (hr == 'nonoverlapping_uneven'){ 
      home_ranges1 <- set_home_ranges_no_overlap_uneven(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    if (hr == 'nonoverlapping_core'){ 
      home_ranges1 <- set_home_ranges_by_non_overlapping_core_areas(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    if (hr == 'nonoverlapping_even'){ 
      home_ranges1 <- set_home_ranges_no_overlap_even(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    ##print('Home ranges done')
    
    if (fitness=='habitat_induced'){
      fitness1 <- habitat_induced_fitness(habitat_raster, home_ranges1, length(home_ranges1))
    }
    
    if (fitness=='density_induced'){
      overlaps <- get_home_range_overlaps(home_ranges1)
      #sum_individuals <- get_individual_densities(habitat_raster, xy_individuals, Wr)
      fitness1 <- population_density_induced_fitness(overlaps)
    }
    
    if (fitness=='inverse_density'){
      overlaps <- get_home_range_overlaps(home_ranges1)
      #sum_individuals <- get_individual_densities(habitat_raster, xy_individuals, Wr)
      fitness1 <- inverse_density_induced_fitness(overlaps)
    }
  
    if (fitness=='neutral'){
      fitness1 <- neutral_fitness(length(home_ranges1))}  
    
    
    reproduction <- fitness_induced_reproduction_reference(individual_cells, fitness1, Nindividuals, reference_fitness)
    
    individual_cells <- reproduction
    xy_immigration<- sample_individuals_to_landscape(initial_habitat, Nimmigration)
    individual_cells_immigrants <-cellFromXY(initial_habitat, xy_immigration)
    individual_cells <- c(individual_cells, individual_cells_immigrants)
    
    habitat_df <- habitat_change_df[which(habitat_change_df$year==y), ]
    sum_habitat <- rasterLocalSums(habitat_raster, W)
    habitat_cells <-  cellFromXY(initial_habitat, cbind(habitat_df$x[habitat_df$habitat==1], habitat_df$y[habitat_df$habitat==1]))
    
    # carefully define whether only newborns, or all individuals disperse (and how this affects the different dispersal rules)
    if (y < 10){
    if (dispersal == 'passive'){ 
      individual_cells <- simulate_passive_dispersal(NB, dsts, individual_cells, habitat_cells)
    }
    if (dispersal == 'conspecific'){ 
      #xy_individuals <- xyFromCell(habitat_raster, individual_cells)
      #sum_individuals <- get_individual_densities(habitat_raster, xy_individuals, Wr)
      individual_cells <- simulate_conspecific_dispersal(NB, dsts,individual_cells, habitat_raster)
    }
    if (dispersal == 'habitat_driven'){ 
      
      individual_cells <- simulate_habitat_driven_dispersal(NB,dsts, individual_cells, habitat_cells, sum_habitat)
    }
    # if (dispersal == 'habitat_and_conspecific_driven'){  
    #   xy_individuals <- xyFromCell(habitat_raster, individual_cells)
    #   sum_individuals <- get_individual_densities(habitat_raster, xy_individuals, Wr)
    #   individual_cells <- simulate_habitat_conspecific_dispersal(NB,dsts, individual_cells, habitat_cells, sum_individuals, sum_habitat)
    # }
    # if (dispersal == 'habitat_conspecific_avoidance_driven'){ 
    #   xy_individuals <- xyFromCell(habitat_raster, individual_cells)
    #   sum_individuals <- get_individual_densities(habitat_raster, xy_individuals, Wr)
    #   individual_cells <- simulate_habitat_avoidance_driven_dispersal(NB,dsts, individual_cells, habitat_cells, sum_individuals, sum_habitat)
    # }
    if (dispersal == 'conspecific_avoidance_driven'){ 
     # xy_individuals <- xyFromCell(habitat_raster, individual_cells)
     # sum_individuals <- get_individual_densities(habitat_raster, xy_individuals, Wr)
      individual_cells <- simulate_avoidance_driven_dispersal(NB,dsts, individual_cells, sum_habitat)
    }
    }
    Nindividuals <- length(individual_cells)
    #print(Nindividuals)
  }
    #print('Dispersal done')
  # 
  # if (hr == 'constant_size'){ 
  #   home_ranges1 <- set_home_ranges_constant_size(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
  # }
  # 
  # if (hr == 'habitat'){ 
  #   home_ranges1 <- set_home_ranges_by_habitat(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
  # }
  # if (hr == 'nonoverlapping_habitat'){ 
  #   home_ranges1 <- set_home_ranges_by_habitat_no_overlap(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)}
  # if (hr == 'nonoverlapping_uneven'){ 
  #   home_ranges1 <- set_home_ranges_no_overlap_uneven(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
  # }
  #   if (hr == 'nonoverlapping_even'){ 
  #     home_ranges1 <- set_home_ranges_no_overlap_even(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
  #   }
  # 
  
  # if (plotting_home_ranges ==TRUE){
  #   if (fitness=='habitat_induced'){
  #     fitness1 <- habitat_induced_fitness(habitat_raster, home_ranges1, length(home_ranges1))
  #   }
  #   if (fitness=='density_induced'){
  #     fitness1 <- population_density_induced_fitness(home_ranges1, length(home_ranges1))}
  #   home_range_name <- paste('hrtype', hr, '.png')
  #   png(home_range_name, height = 7, width = 7, units = 'in', res = 300)
  #   plot_home_ranges(home_ranges1,  habitat_raster,length(home_ranges1), fitness1)
  #   dev.off()
  # }
  #
  res <- list('individual_cells' = individual_cells_5, 'home_ranges' = home_ranges5,
              'individual_cells10' = individual_cells_10, 'home_ranges10' = home_ranges1, 'fitness10'= fitness1)

  return(res)
}

simulate_time_series <- function(initial_habitat, individual_cells){
  
Nyears <- 10
Nhabitat <- 6 

for (y in seq(1, Nyears)){
  
  home_ranges1 <- set_home_ranges_by_habitat(initial_habitat, individual_cells, Nhabitat, Nindividuals)
  #home_ranges2 <- set_home_ranges_by_habitat_no_overlap(initial_habitat, individual_cells, Nhabitat, Nindividuals)
  #home_ranges3 <- set_home_ranges_no_overlap(initial_habitat, individual_cells, Nhabitat, Nindividuals)
  
  fitness1 <- habitat_induced_fitness(initial_habitat, home_ranges1, length(home_ranges1))
  #fitness2 <- population_density_induced_fFitness(home_ranges1, length(home_ranges1))
  reproduction <- fitness_induced_reproduction(individual_cells, fitness1, Nindividuals)
  individual_cells <- reproduction
        
}
home_ranges1 <- set_home_ranges_by_habitat(initial_habitat, individual_cells, Nhabitat, Nindividuals)
res <- list('individual_cells' = individual_cells, 'home_ranges' = home_ranges1)
return(res)
}



simulate_time_series_environmental_change_stable_fitness_plot5 <- function(initial_habitat_raster, habitat_change_df, individual_cells, hr, fitness, dispersal, NB, dsts, max_diameter, Nhabitat, Wr){
  
  W <- Wr
  Nyears <- 10
  Nindividuals <- length(individual_cells)
  
  
  
  for (y in seq(1, 10)){
    
    if (y==5){
      individual_cells_5 <- individual_cells
    }
    
    if (y==10){
      individual_cells_10 <- individual_cells
    }
    
    Nimmigration <- sample(seq(1, 10), 1)
    habitat_raster <- generate_raster_from_cells(initial_habitat_raster, habitat_change_df, y)
    
    
    
    # 1 .home ranges -------------------------------------------------------------
    
    Nindividuals <- length(individual_cells)
    
    if (hr == 'constant_size'){ 
      home_ranges1 <- set_home_ranges_constant_size(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    if (hr == 'habitat'){ 
      home_ranges1 <- set_home_ranges_by_habitat(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    if (hr == 'nonoverlapping_habitat'){ 
      home_ranges1 <- set_home_ranges_by_habitat_no_overlap(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    if (hr == 'nonoverlapping_uneven'){ 
      home_ranges1 <- set_home_ranges_no_overlap_uneven(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    if (hr == 'nonoverlapping_even'){ 
      home_ranges1 <- set_home_ranges_no_overlap_even(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    if (hr == 'nonoverlapping_core'){ 
      
      home_ranges1 <- set_home_ranges_by_non_overlapping_core_areas(habitat_raster, individual_cells, Nhabitat, Nindividuals, max_diameter)
    }
    
    # print('Home ranges occurred')
    # print(length(individual_cells))
    # print(length(home_ranges1))
    # 
    
    # 2 .fitness -----------------------------------------------------------------
    
    
    if (fitness=='habitat_induced'){
      fitness1 <- habitat_induced_fitness(habitat_raster, home_ranges1, length(home_ranges1))
    }
    
    if (fitness=='density_induced'){
      overlaps <- get_home_range_overlaps(home_ranges1)
      #sum_individuals <- get_individual_densities(habitat_raster, xy_individuals, Wr)
      fitness1 <- population_density_induced_fitness(overlaps)
    }
    
    if (fitness=='inverse_density'){
      overlaps <- get_home_range_overlaps(home_ranges1)
      #sum_individuals <- get_individual_densities(habitat_raster, xy_individuals, Wr)
      fitness1 <- inverse_density_induced_fitness(overlaps)
    }
    
    if (fitness=='neutral'){
      # print(home_ranges1)
      fitness1 <- neutral_fitness(length(home_ranges1))}  
    
    # hist(fitness1)
    # plot(sum_individuals, main = sum(values(sum_individuals)))
    #  points(xy_individuals, cex = fitness1)
    
    #print('Fitness defined')
    
    # 3. let individuals reproduce
    #print(fitness1)
    #print(length(individual_cells))
    #print(length(fitness1))
    #print(y)
    
    reproduction <- fitness_induced_reproduction(individual_cells, fitness1, Nindividuals)
    
    # print('Reproduction occurred')
    # individual_cells <- reproduction
    # print(length(individual_cells))
    
    xy_immigration<- sample_individuals_to_landscape(initial_habitat, Nimmigration)
    individual_cells_immigrants <-cellFromXY(initial_habitat, xy_immigration)
    #print("Immigration occurred!")
    #print(length(individual_cells))
    
    individual_cells <- c(individual_cells, individual_cells_immigrants)
    
    habitat_df <- habitat_change_df[which(habitat_change_df$year==y), ]
    sum_habitat <- rasterLocalSums(habitat_raster, Wr)
    habitat_cells <-  cellFromXY(initial_habitat, cbind(habitat_df$x[habitat_df$habitat==1], habitat_df$y[habitat_df$habitat==1]))
    
    #print("Something occurred!")
    #print(length(individual_cells))
    
    
    # 4. simulate dispersal!-------------------------------------------------------------------------
    
    if (y < 10){
      if (dispersal == 'passive'){ 
        individual_cells <- simulate_passive_dispersal(NB, dsts, individual_cells, habitat_cells)
      }
      
      if (dispersal == 'conspecific'){ 
        # xy_individuals <- xyFromCell(habitat_raster, individual_cells)
        #  sum_individuals <- get_individual_densities(habitat_raster, xy_individuals, Wr)
        individual_cells <- simulate_conspecific_dispersal(NB, dsts, individual_cells, habitat_raster)
      }
      
      if (dispersal == 'habitat_driven'){ 
        individual_cells <- simulate_habitat_driven_dispersal(NB, dsts, individual_cells, habitat_cells, sum_habitat)
      }
      
      if (dispersal == 'conspecific_avoidance_driven'){ 
        individual_cells <- simulate_avoidance_driven_dispersal(NB, dsts,individual_cells,  sum_habitat)
      }
    }
    
    if (y ==5){
      individual_cells5 <- individual_cells
      home_ranges5 <-home_ranges1 
    }
    #print('Dispersal occurred')
    #print(length(individual_cells))
    
    if (y > 5){
      hab <- generate_raster_from_cells(initial_habitat_raster, habitat_change_df, y)
      #plot_home_ranges_nocolor2(home_ranges1, hab, length(home_ranges1))
      if (y == 6){
        home_ranges6 <- home_ranges1
        individual_cells6 <- individual_cells
      }
      if (y == 7){
        home_ranges7 <- home_ranges1
        individual_cells7 <- individual_cells
      }
      if (y == 8){
        home_ranges8 <- home_ranges1
        individual_cells8 <- individual_cells
      }
      if (y == 9){
        home_ranges9 <- home_ranges1
        individual_cells9 <- individual_cells
      }
    } 
    
  }
  
  
  res <- list('individual_cells' = individual_cells_5, 'home_ranges' = home_ranges5,
              'home_ranges6' = home_ranges6, 'home_ranges7' = home_ranges7,
              'home_ranges8' = home_ranges8,'home_ranges9' = home_ranges9,
              'individual_cells6' = individual_cells6,'individual_cells7' = individual_cells7,
              'individual_cells8' = individual_cells8,'individual_cells9' = individual_cells9,
              'individual_cells10' = individual_cells, 'home_ranges10' = home_ranges1, 'fitness10' = fitness1)
  
  return(res)
}

get_individual_densities <- function(r, xy_individuals, Wr){
  tab <- table(cellFromXY(r, xy_individuals))
  pd <- r
  #pd <- raster(ncol=100, nrow=100, xmn=1, xmx=100, ymn=1, ymx=100)
  # raster giving counts of individuals in each cells
  values(pd) <- 0
  pd[as.numeric(names(tab))] <- tab
  #plot(pd)
  
  # compute the moving average habitat-density-matrix: 
  #W <- matrix(1,3,3)
  #W <- focalWeight(r, 2, "Gauss")
  #?rasterLocalSums
  par(mfrow=c(1,2))
  #plot(pd)
  sum_individuals <- rasterLocalSums(pd,Wr)
  #plot(sum_individuals)
  return(sum_individuals)
}

sample_individuals_to_landscape <- function(r, Nindividuals){
  # give habitat as a raster
  ones = xyFromCell(r,1:prod(dim(r)))[getValues(r)==1,]
  #head(ones)
  loc_ids = sample(nrow(ones), Nindividuals, replace=TRUE)
  return(ones[loc_ids,])
}


simulate_observations_environmental_change_plot <- function(initial_habitat_raster, habitat_change_df, individual_cells, Nyears, plotting, hr, fitness, dispersal, NB, dsts, max_diameter, Nhabitat, Wr, obs_list, fitness_assumption){
  
  #Nhabitat <- 6 # check the exact definition of this parameter. 
  N_locations <-1 # (these refer to individuals' behavior)
  habitat_dependency  <- FALSE  #check the exact definition of this parameter. 
  
  #result <- simulate_time_series(initial_habitat, individual_cells) 
  
  # here we can compute spatial summary statistics for the spatial clustering of individuals
  
  
  if (fitness_assumption == 'changing'){
    result <- simulate_time_series_environmental_change(initial_habitat_raster, habitat_change_df, individual_cells, hr, fitness, dispersal, NB,dsts, max_diameter, Nhabitat, Wr)
  } else {
    result <- simulate_time_series_environmental_change_stable_fitness(initial_habitat_raster, habitat_change_df, individual_cells, hr, fitness, dispersal, NB,dsts, max_diameter, Nhabitat, Wr) 
  }
  
  
  raster_end <- generate_raster_from_cells(initial_habitat_raster, habitat_change_df, 5)
  Nindividuals_end <- length(result$individual_cells)
  loc_xy <- set_individual_locations(raster_end, result$home_ranges, Nhabitat, Nindividuals_end, always_in_habitat, N_locations, habitat_dependency)
  
  # compute summary statistics for the spatial distribution of the individuals
  
  #observations at time point 5: 
  
  ppp <- as.ppp(loc_xy$individual_locations_xy, c(0,100,0,100))
  # habitat_Wr <- focal(raster_end, Wr)
  # habitat_image <- as.im(as.matrix(habitat_Wr))
  # rhohat <- rhohat(ppp, habitat_image)
  
  nn1 <- nndist(ppp, k=1)
  #Kest <- Kest(ppp)
  observation_probability <- 1
  
  # o1:
  
  obs_within_squares_o1 <- allocate_locations_within_squares(loc_xy$individual_locations_xy, obs_list$obs_cells5_o1 )
  res_squares_intersecting_o1 <- observe_individuals_intersecting_homerange_squares(raster_end, result$home_ranges, obs_list$obs_cells5_o1)
  Nsquares  <- length(obs_list$obs_cells5_o1)
  Ntransects  <- length(obs_list$obs_transects5_o1)
  res_o1 <- observe_individuals_within_squares(obs_within_squares_o1$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy$individual_locations_xy, obs_list$obs_transects5_o1)
  res_transects_o1 <- observe_individuals_along_transects(distance_matrix, Ntransects)
  res_transects_intersecting_o1 <- observe_individuals_intersecting_homerange_transects(raster_end, result$home_ranges, obs_list$obs_transects5_o1,  Nindividuals_end, Ntransects)
  
  #print('o1')
  # o2:
  
  obs_within_squares_o2 <- allocate_locations_within_squares(loc_xy$individual_locations_xy, obs_list$obs_cells5_o2 )
  res_squares_intersecting_o2 <- observe_individuals_intersecting_homerange_squares(raster_end, result$home_ranges, obs_list$obs_cells5_o2)
  Nsquares  <- length(obs_list$obs_cells5_o2)
  Ntransects  <- length(obs_list$obs_transects5_o2)
  res_o2 <- observe_individuals_within_squares(obs_within_squares_o2$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy$individual_locations_xy, obs_list$obs_transects5_o2)
  res_transects_o2 <- observe_individuals_along_transects(distance_matrix, Ntransects)
  res_transects_intersecting_o2 <- observe_individuals_intersecting_homerange_transects(raster_end, result$home_ranges, obs_list$obs_transects5_o2,  Nindividuals_end, Ntransects)
  
  #print('o2') 
  # o3:
  
  obs_within_squares_o3 <- allocate_locations_within_squares(loc_xy$individual_locations_xy, obs_list$obs_cells5_o3 )
  res_squares_intersecting_o3 <- observe_individuals_intersecting_homerange_squares(raster_end, result$home_ranges, obs_list$obs_cells5_o3)
  Nsquares  <- length(obs_list$obs_cells5_o3)
  Ntransects  <- length(obs_list$obs_transects5_o3)
  res_o3 <- observe_individuals_within_squares(obs_within_squares_o3$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy$individual_locations_xy, obs_list$obs_transects5_o3)
  res_transects_o3 <- observe_individuals_along_transects(distance_matrix, Ntransects)
  res_transects_intersecting_o3 <- observe_individuals_intersecting_homerange_transects(raster_end, result$home_ranges, obs_list$obs_transects5_o3,  Nindividuals_end, Ntransects)
  
  #print('o3')
  
  # o4:
  
  obs_within_squares_o4 <- allocate_locations_within_squares(loc_xy$individual_locations_xy, obs_list$obs_cells5_o4 )
  res_squares_intersecting_o4 <- observe_individuals_intersecting_homerange_squares(raster_end, result$home_ranges, obs_list$obs_cells5_o4)
  Nsquares  <- length(obs_list$obs_cells5_o4)
  Ntransects  <- length(obs_list$obs_transects5_o4)
  res_o4 <- observe_individuals_within_squares(obs_within_squares_o4$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy$individual_locations_xy, obs_list$obs_transects5_o4)
  res_transects_o4 <- observe_individuals_along_transects(distance_matrix, Ntransects)
  res_transects_intersecting_o4 <- observe_individuals_intersecting_homerange_transects(raster_end, result$home_ranges, obs_list$obs_transects5_o4,  Nindividuals_end, Ntransects)
  
  #print('o4')
  
  # o5:
  
  obs_within_squares_o5 <- allocate_locations_within_squares(loc_xy$individual_locations_xy, obs_list$obs_cells5_o5 )
  res_squares_intersecting_o5 <- observe_individuals_intersecting_homerange_squares(raster_end, result$home_ranges, obs_list$obs_cells5_o5)
  Nsquares  <- length(obs_list$obs_cells5_o5)
  Ntransects  <- length(obs_list$obs_transects5_o5)
  res_o5 <- observe_individuals_within_squares(obs_within_squares_o5$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy$individual_locations_xy, obs_list$obs_transects5_o5)
  res_transects_o5 <- observe_individuals_along_transects(distance_matrix, Ntransects)
  res_transects_intersecting_o5 <- observe_individuals_intersecting_homerange_transects(raster_end, result$home_ranges, obs_list$obs_transects5_o5,  Nindividuals_end, Ntransects)
  
  #print('o5')
  
  
  
  
  # observations at time point 10: 
  
  raster_end10 <- generate_raster_from_cells(initial_habitat_raster, habitat_change_df, 10)
  Nindividuals_end10 <- length(result$individual_cells10)
  loc_xy10 <- set_individual_locations(raster_end10, result$home_ranges10, Nhabitat, Nindividuals_end10, always_in_habitat, N_locations, habitat_dependency)
  
  
  #print('now here!')
  ppp <- as.ppp(loc_xy10$individual_locations_xy, c(0,100,0,100))
  #print('now here! 2')
  #habitat_Wr <- focal(raster_end10, Wr)
  #print('now here! 3')
  #habitat_image <- as.im(as.matrix(habitat_Wr))
  #print('now here! 4')
  #rhohat10 <- rhohat(ppp, habitat_image)
  #print('now here! 5')
  nn110 <- nndist(ppp, k=1)
  #print('now here! 6')
  #Kest10 <- Kest(ppp)
  # if (TRUE){
  #   plot(habitat_Wr)
  #   points(loc_xy10$individual_locations_xy)  
  # } 
  
  # o1:
  
  obs_within_squares_o1_10 <- allocate_locations_within_squares(loc_xy10$individual_locations_xy, obs_list$obs_cells10_o1 )
  res_squares_intersecting_o1_10 <- observe_individuals_intersecting_homerange_squares(raster_end10, result$home_ranges10, obs_list$obs_cells10_o1)
  Nsquares  <- length(obs_list$obs_cells10_o1)
  Ntransects  <- length(obs_list$obs_transects10_o1)
  res_o1_10 <- observe_individuals_within_squares(obs_within_squares_o1_10$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy10$individual_locations_xy, obs_list$obs_transects10_o1)
  res_transects_o1_10 <- observe_individuals_along_transects(distance_matrix, Ntransects)
  res_transects_intersecting_o1_10 <- observe_individuals_intersecting_homerange_transects(raster_end10, result$home_ranges10, obs_list$obs_transects10_o1,  Nindividuals_end, Ntransects)
  #print('o1')
  # o2:
  
  obs_within_squares_o2_10 <- allocate_locations_within_squares(loc_xy10$individual_locations_xy, obs_list$obs_cells10_o2 )
  #print('now here! 68')
  res_squares_intersecting_o2_10 <- observe_individuals_intersecting_homerange_squares(raster_end, result$home_ranges10, obs_list$obs_cells10_o2)
  Nsquares  <- length(obs_list$obs_cells10_o2)
  Ntransects  <- length(obs_list$obs_transects10_o2)
  #print('now here! 69')
  res_o2_10 <- observe_individuals_within_squares(obs_within_squares_o2_10$locations_per_square, observation_probability, Nsquares)
  #print('now here! 70')
  distance_matrix <-compute_distances_to_transects(loc_xy10$individual_locations_xy, obs_list$obs_transects10_o2)
  #print('now here! 71')
  res_transects_o2_10 <- observe_individuals_along_transects(distance_matrix, Ntransects)
  #print('now here! 72')
  res_transects_intersecting_o2_10 <- observe_individuals_intersecting_homerange_transects(raster_end, result$home_ranges10, obs_list$obs_transects10_o2,  Nindividuals_end, Ntransects)
  #print('o2')
  #print('now here! 73')
  
  # o3:
  
  obs_within_squares_o3_10 <- allocate_locations_within_squares(loc_xy10$individual_locations_xy, obs_list$obs_cells10_o3 )
  res_squares_intersecting_o3_10 <- observe_individuals_intersecting_homerange_squares(raster_end10, result$home_ranges10, obs_list$obs_cells10_o3)
  Nsquares  <- length(obs_list$obs_cells10_o3)
  Ntransects  <- length(obs_list$obs_transects10_o3)
  res_o3_10 <- observe_individuals_within_squares(obs_within_squares_o3_10$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy10$individual_locations_xy, obs_list$obs_transects10_o3)
  res_transects_o3_10 <- observe_individuals_along_transects(distance_matrix, Ntransects)
  res_transects_intersecting_o3_10 <- observe_individuals_intersecting_homerange_transects(raster_end10, result$home_ranges10, obs_list$obs_transects10_o3,  Nindividuals_end, Ntransects)
  #print('o3')
  
  # o4:
  
  obs_within_squares_o4_10 <- allocate_locations_within_squares(loc_xy10$individual_locations_xy, obs_list$obs_cells10_o4 )
  res_squares_intersecting_o4_10 <- observe_individuals_intersecting_homerange_squares(raster_end10, result$home_ranges10, obs_list$obs_cells10_o4)
  Nsquares  <- length(obs_list$obs_cells10_o4)
  Ntransects  <- length(obs_list$obs_transects10_o4)
  res_o4_10 <- observe_individuals_within_squares(obs_within_squares_o4_10$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy10$individual_locations_xy, obs_list$obs_transects10_o4)
  res_transects_o4_10 <- observe_individuals_along_transects(distance_matrix, Ntransects)
  res_transects_intersecting_o4_10 <- observe_individuals_intersecting_homerange_transects(raster_end10, result$home_ranges10, obs_list$obs_transects10_o4,  Nindividuals_end, Ntransects)
  #print('o4')
  
  # o5:
  
  obs_within_squares_o5_10 <- allocate_locations_within_squares(loc_xy10$individual_locations_xy, obs_list$obs_cells10_o5 )
  res_squares_intersecting_o5_10 <- observe_individuals_intersecting_homerange_squares(raster_end10, result$home_ranges10, obs_list$obs_cells10_o5)
  Nsquares  <- length(obs_list$obs_cells10_o5)
  Ntransects  <- length(obs_list$obs_transects10_o5)
  res_o5_10 <- observe_individuals_within_squares(obs_within_squares_o5_10$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy10$individual_locations_xy, obs_list$obs_transects10_o5)
  res_transects_o5_10 <- observe_individuals_along_transects(distance_matrix, Ntransects)
  res_transects_intersecting_o5_10 <- observe_individuals_intersecting_homerange_transects(raster_end10, result$home_ranges10, obs_list$obs_transects10_o5,  Nindividuals_end, Ntransects)
  #print('o5')
  
  Nindividuals <- length(loc_xy10$individual_locations[,1])
  plot_home_ranges_nocolor(result$home_ranges10, raster_end10, Nindividuals, hr, result$individual_cells10)
  lines(obs_cells10_o2)
  
  res <-list('pop_size' = length(loc_xy$individual_locations_xy[,1]), 
             'pop_size10' = length(loc_xy10$individual_locations_xy[,1]),
             'home_ranges' = result$home_ranges, 
             'Nsquares' = Nsquares, 
             'home_ranges10' = result$home_ranges10,
             
             'presence_observed_o1' = res_o1$presence_observed, 
             'number_of_tracks_observed_o1' = res_o1$number_of_tracks_observed, 
             'presence_observed_transects_o1' = res_transects_o1$presence_observed, 
             'number_of_tracks_observed_transects_o1' = res_transects_o1$number_of_tracks_observed, 
             'presence_observed_intersecting_squares_o1' = res_squares_intersecting_o1$presence_observed,
             'presence_observed_intersecting_transects_o1' = 1*(res_transects_intersecting_o1$presence_observed>0), 
             
             
             'presence_observed_o2' = res_o2$presence_observed, 
             'number_of_tracks_observed_o2' = res_o2$number_of_tracks_observed, 
             'presence_observed_transects_o2' = res_transects_o2$presence_observed, 
             'number_of_tracks_observed_transects_o2' = res_transects_o2$number_of_tracks_observed, 
             'presence_observed_intersecting_squares_o2' = res_squares_intersecting_o2$presence_observed,
             'presence_observed_intersecting_transects_o2' = 1*(res_transects_intersecting_o2$presence_observed>0), 
             
             'presence_observed_o3' = res_o3$presence_observed, 
             'number_of_tracks_observed_o3' = res_o3$number_of_tracks_observed, 
             'presence_observed_transects_o3' = res_transects_o3$presence_observed, 
             'number_of_tracks_observed_transects_o3' = res_transects_o3$number_of_tracks_observed, 
             'presence_observed_intersecting_squares_o3' = res_squares_intersecting_o3$presence_observed,
             'presence_observed_intersecting_transects_o3' = 1*(res_transects_intersecting_o3$presence_observed>0), 
             
             'presence_observed_o4' = res_o4$presence_observed, 
             'number_of_tracks_observed_o4' = res_o4$number_of_tracks_observed, 
             'presence_observed_transects_o4' = res_transects_o4$presence_observed, 
             'number_of_tracks_observed_transects_o4' = res_transects_o4$number_of_tracks_observed, 
             'presence_observed_intersecting_squares_o4' = res_squares_intersecting_o4$presence_observed,
             'presence_observed_intersecting_transects_o4' = 1*(res_transects_intersecting_o4$presence_observed>0), 
             
             'presence_observed_o5' = res_o5$presence_observed, 
             'number_of_tracks_observed_o5' = res_o5$number_of_tracks_observed, 
             'presence_observed_transects_o5' = res_transects_o5$presence_observed, 
             'number_of_tracks_observed_transects_o5' = res_transects_o5$number_of_tracks_observed, 
             'presence_observed_intersecting_squares_o5' = res_squares_intersecting_o5$presence_observed,
             'presence_observed_intersecting_transects_o5' = 1*(res_transects_intersecting_o5$presence_observed>0), 
             
             'presence_observed_o1_10' = res_o1_10$presence_observed, 
             'number_of_tracks_observed_o1_10' = res_o1_10$number_of_tracks_observed, 
             'presence_observed_transects_o1_10' = res_transects_o1_10$presence_observed, 
             'number_of_tracks_observed_transects_o1_10' = res_transects_o1_10$number_of_tracks_observed, 
             'presence_observed_intersecting_squares_o1_10' = res_squares_intersecting_o1_10$presence_observed,
             'presence_observed_intersecting_transects_o1_10' = 1*(res_transects_intersecting_o1_10$presence_observed>0), 
             
             
             'presence_observed_o2_10' = res_o2_10$presence_observed, 
             'number_of_tracks_observed_o2_10' = res_o2_10$number_of_tracks_observed, 
             'presence_observed_transects_o2_10' = res_transects_o2_10$presence_observed, 
             'number_of_tracks_observed_transects_o2_10' = res_transects_o2_10$number_of_tracks_observed, 
             'presence_observed_intersecting_squares_o2_10' = res_squares_intersecting_o2_10$presence_observed,
             'presence_observed_intersecting_transects_o2_10' = 1*(res_transects_intersecting_o2_10$presence_observed>0), 
             
             'presence_observed_o3_10' = res_o3_10$presence_observed, 
             'number_of_tracks_observed_o3_10' = res_o3_10$number_of_tracks_observed, 
             'presence_observed_transects_o3_10' = res_transects_o3_10$presence_observed, 
             'number_of_tracks_observed_transects_o3_10' = res_transects_o3_10$number_of_tracks_observed, 
             'presence_observed_intersecting_squares_o3_10' = res_squares_intersecting_o3_10$presence_observed,
             'presence_observed_intersecting_transects_o3_10' = 1*(res_transects_intersecting_o3_10$presence_observed>0), 
             
             'presence_observed_o4_10' = res_o4_10$presence_observed, 
             'number_of_tracks_observed_o4_10' = res_o4_10$number_of_tracks_observed, 
             'presence_observed_transects_o4_10' = res_transects_o4_10$presence_observed, 
             'number_of_tracks_observed_transects_o4_10' = res_transects_o4_10$number_of_tracks_observed, 
             'presence_observed_intersecting_squares_o4_10' = res_squares_intersecting_o4_10$presence_observed,
             'presence_observed_intersecting_transects_o4_10' = 1*(res_transects_intersecting_o4_10$presence_observed>0), 
             
             'presence_observed_o5_10' = res_o5_10$presence_observed, 
             'number_of_tracks_observed_o5_10' = res_o5_10$number_of_tracks_observed, 
             'presence_observed_transects_o5_10' = res_transects_o5_10$presence_observed, 
             'number_of_tracks_observed_transects_o5_10' = res_transects_o5_10$number_of_tracks_observed, 
             'presence_observed_intersecting_squares_o5_10' = res_squares_intersecting_o5_10$presence_observed,
             'presence_observed_intersecting_transects_o5_10' = 1*(res_transects_intersecting_o5_10$presence_observed>0), 
             
             
             
             # 'presence_observed10' = res_o10$presence_observed, 
             # 'number_of_tracks_observed10' = res_o10$number_of_tracks_observed, 
             # 'presence_observed_transects10' = res_transects10$presence_observed, 
             # 'number_of_tracks_observed_transects10' = res_transects10$number_of_tracks_observed, 
             # 'presence_observed_intersecting_squares10' = res_squares_intersecting10$presence_observed,
             # 'presence_observed_intersecting_transects10' = 1*(res_transects_intersecting10$presence_observed>0), 
             
             # 'rhohat'= rhohat,
             'nn1'= nn1,
             #'Kest' = Kest,
             #'rhohat10'= rhohat10,
             'nn110'= nn110
             #'Kest10' = Kest10
  )
  return(res)
}


simulate_observations_environmental_change_plot5 <- function(initial_habitat_raster, habitat_change_df, individual_cells, Nyears, plotting, hr, fitness, dispersal, NB, dsts, max_diameter, Nhabitat, Wr, obs_list, fitness_assumption, label){
  
  #Nhabitat <- 6 # check the exact definition of this parameter. 
  N_locations <-1 # (these refer to individuals' behavior)
  habitat_dependency  <- FALSE  #check the exact definition of this parameter. 
  
  
  
  
  if (fitness_assumption == 'changing'){
    result <- simulate_time_series_environmental_change(initial_habitat_raster, habitat_change_df, individual_cells, hr, fitness, dispersal, NB,dsts, max_diameter, Nhabitat, Wr)
  } else {
    result <- simulate_time_series_environmental_change_stable_fitness_plot5(initial_habitat_raster, habitat_change_df, individual_cells, hr, fitness, dispersal, NB,dsts, max_diameter, Nhabitat, Wr) 
  }
  
  
  raster_end5 <- generate_raster_from_cells(initial_habitat_raster, habitat_change_df, 5)
  raster_end6 <- generate_raster_from_cells(initial_habitat_raster, habitat_change_df, 6)
  raster_end7 <- generate_raster_from_cells(initial_habitat_raster, habitat_change_df, 7)
  raster_end8 <- generate_raster_from_cells(initial_habitat_raster, habitat_change_df, 8)
  raster_end9 <- generate_raster_from_cells(initial_habitat_raster, habitat_change_df, 9)
  raster_end10 <- generate_raster_from_cells(initial_habitat_raster, habitat_change_df, 10)
  raster_end <- raster_end10
  
  Nindividuals_end <- length(result$individual_cells)
  loc_xy <- set_individual_locations(raster_end, result$home_ranges, Nhabitat, Nindividuals_end, always_in_habitat, N_locations, habitat_dependency)
  
  # compute summary statistics for the spatial distribution of the individuals
  
  #observations at time point 5: 
  
  ppp <- as.ppp(loc_xy$individual_locations_xy, c(0,100,0,100))
  # habitat_Wr <- focal(raster_end, Wr)
  # habitat_image <- as.im(as.matrix(habitat_Wr))
  # rhohat <- rhohat(ppp, habitat_image)
  
  nn1 <- nndist(ppp, k=1)
  #Kest <- Kest(ppp)
  observation_probability <- 1
  
  # o1:
  
  obs_within_squares_o1 <- allocate_locations_within_squares(loc_xy$individual_locations_xy, obs_list$obs_cells5_o1 )
  res_squares_intersecting_o1 <- observe_individuals_intersecting_homerange_squares(raster_end, result$home_ranges, obs_list$obs_cells5_o1)
  Nsquares  <- length(obs_list$obs_cells5_o1)
  Ntransects  <- length(obs_list$obs_transects5_o1)
  res_o1 <- observe_individuals_within_squares(obs_within_squares_o1$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy$individual_locations_xy, obs_list$obs_transects5_o1)
  res_transects_o1 <- observe_individuals_along_transects_2(distance_matrix, Ntransects)
  res_transects_intersecting_o1 <- observe_individuals_intersecting_homerange_transects(raster_end, result$home_ranges, obs_list$obs_transects5_o1,  Nindividuals_end, Ntransects)
  
  #print('o1')
  # o2:
  
  obs_within_squares_o2 <- allocate_locations_within_squares(loc_xy$individual_locations_xy, obs_list$obs_cells5_o2 )
  res_squares_intersecting_o2 <- observe_individuals_intersecting_homerange_squares(raster_end, result$home_ranges, obs_list$obs_cells5_o2)
  Nsquares  <- length(obs_list$obs_cells5_o2)
  Ntransects  <- length(obs_list$obs_transects5_o2)
  res_o2 <- observe_individuals_within_squares(obs_within_squares_o2$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy$individual_locations_xy, obs_list$obs_transects5_o2)
  res_transects_o2 <- observe_individuals_along_transects_2(distance_matrix, Ntransects)
  res_transects_intersecting_o2 <- observe_individuals_intersecting_homerange_transects(raster_end, result$home_ranges, obs_list$obs_transects5_o2,  Nindividuals_end, Ntransects)
  
  #print('o2') 
  # o3:
  
  obs_within_squares_o3 <- allocate_locations_within_squares(loc_xy$individual_locations_xy, obs_list$obs_cells5_o3 )
  res_squares_intersecting_o3 <- observe_individuals_intersecting_homerange_squares(raster_end, result$home_ranges, obs_list$obs_cells5_o3)
  Nsquares  <- length(obs_list$obs_cells5_o3)
  Ntransects  <- length(obs_list$obs_transects5_o3)
  res_o3 <- observe_individuals_within_squares(obs_within_squares_o3$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy$individual_locations_xy, obs_list$obs_transects5_o3)
  res_transects_o3 <- observe_individuals_along_transects_2(distance_matrix, Ntransects)
  res_transects_intersecting_o3 <- observe_individuals_intersecting_homerange_transects(raster_end, result$home_ranges, obs_list$obs_transects5_o3,  Nindividuals_end, Ntransects)
  
  #print('o3')
  
  # o4:
  
  obs_within_squares_o4 <- allocate_locations_within_squares(loc_xy$individual_locations_xy, obs_list$obs_cells5_o4 )
  res_squares_intersecting_o4 <- observe_individuals_intersecting_homerange_squares(raster_end, result$home_ranges, obs_list$obs_cells5_o4)
  Nsquares  <- length(obs_list$obs_cells5_o4)
  Ntransects  <- length(obs_list$obs_transects5_o4)
  res_o4 <- observe_individuals_within_squares(obs_within_squares_o4$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy$individual_locations_xy, obs_list$obs_transects5_o4)
  res_transects_o4 <- observe_individuals_along_transects_2(distance_matrix, Ntransects)
  res_transects_intersecting_o4 <- observe_individuals_intersecting_homerange_transects(raster_end, result$home_ranges, obs_list$obs_transects5_o4,  Nindividuals_end, Ntransects)
  
  #print('o4')
  
  # o5:
  
  obs_within_squares_o5 <- allocate_locations_within_squares(loc_xy$individual_locations_xy, obs_list$obs_cells5_o5 )
  res_squares_intersecting_o5 <- observe_individuals_intersecting_homerange_squares(raster_end, result$home_ranges, obs_list$obs_cells5_o5)
  Nsquares  <- length(obs_list$obs_cells5_o5)
  Ntransects  <- length(obs_list$obs_transects5_o5)
  res_o5 <- observe_individuals_within_squares(obs_within_squares_o5$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy$individual_locations_xy, obs_list$obs_transects5_o5)
  res_transects_o5 <- observe_individuals_along_transects_2(distance_matrix, Ntransects)
  res_transects_intersecting_o5 <- observe_individuals_intersecting_homerange_transects(raster_end, result$home_ranges, obs_list$obs_transects5_o5,  Nindividuals_end, Ntransects)
  
  #print('o5')
  
  
  
  
  # observations at time point 10: 
  
  raster_end10 <- generate_raster_from_cells(initial_habitat_raster, habitat_change_df, 10)
  Nindividuals_end10 <- length(result$individual_cells10)
  loc_xy10 <- set_individual_locations(raster_end10, result$home_ranges10, Nhabitat, Nindividuals_end10, always_in_habitat, N_locations, habitat_dependency)
  
  
  #print('now here!')
  ppp <- as.ppp(loc_xy10$individual_locations_xy, c(0,100,0,100))
  #print('now here! 2')
  #habitat_Wr <- focal(raster_end10, Wr)
  #print('now here! 3')
  #habitat_image <- as.im(as.matrix(habitat_Wr))
  #print('now here! 4')
  #rhohat10 <- rhohat(ppp, habitat_image)
  #print('now here! 5')
  nn110 <- nndist(ppp, k=1)
  #print('now here! 6')
  #Kest10 <- Kest(ppp)
  # if (TRUE){
  #   plot(habitat_Wr)
  #   points(loc_xy10$individual_locations_xy)  
  # } 
  
  # o1:
  
  obs_within_squares_o1_10 <- allocate_locations_within_squares(loc_xy10$individual_locations_xy, obs_list$obs_cells10_o1 )
  res_squares_intersecting_o1_10 <- observe_individuals_intersecting_homerange_squares(raster_end10, result$home_ranges10, obs_list$obs_cells10_o1)
  Nsquares  <- length(obs_list$obs_cells10_o1)
  Ntransects  <- length(obs_list$obs_transects10_o1)
  res_o1_10 <- observe_individuals_within_squares(obs_within_squares_o1_10$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy10$individual_locations_xy, obs_list$obs_transects10_o1)
  res_transects_o1_10 <- observe_individuals_along_transects_2(distance_matrix, Ntransects)
  res_transects_intersecting_o1_10 <- observe_individuals_intersecting_homerange_transects(raster_end10, result$home_ranges10, obs_list$obs_transects10_o1,  Nindividuals_end, Ntransects)
  #print('o1')
  # o2:
  
  obs_within_squares_o2_10 <- allocate_locations_within_squares(loc_xy10$individual_locations_xy, obs_list$obs_cells10_o2 )
  
  #print('now here! 68')
  res_squares_intersecting_o2_10 <- observe_individuals_intersecting_homerange_squares(raster_end, result$home_ranges10, obs_list$obs_cells10_o2)
  Nsquares  <- length(obs_list$obs_cells10_o2)
  Ntransects  <- length(obs_list$obs_transects10_o2)
  #print('now here! 69')
  res_o2_10 <- observe_individuals_within_squares(obs_within_squares_o2_10$locations_per_square, observation_probability, Nsquares)
  #print('now here! 70')
  distance_matrix <-compute_distances_to_transects(loc_xy10$individual_locations_xy, obs_list$obs_transects10_o2)
  #print('now here! 71')
  res_transects_o2_10 <- observe_individuals_along_transects_2(distance_matrix, Ntransects)
  #print('now here! 72')
  res_transects_intersecting_o2_10 <- observe_individuals_intersecting_homerange_transects(raster_end, result$home_ranges10, obs_list$obs_transects10_o2,  Nindividuals_end, Ntransects)
  #print('o2')
  #print('now here! 73')
  
  # o3:
  
  obs_within_squares_o3_10 <- allocate_locations_within_squares(loc_xy10$individual_locations_xy, obs_list$obs_cells10_o3 )
  res_squares_intersecting_o3_10 <- observe_individuals_intersecting_homerange_squares(raster_end10, result$home_ranges10, obs_list$obs_cells10_o3)
  Nsquares  <- length(obs_list$obs_cells10_o3)
  Ntransects  <- length(obs_list$obs_transects10_o3)
  res_o3_10 <- observe_individuals_within_squares(obs_within_squares_o3_10$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy10$individual_locations_xy, obs_list$obs_transects10_o3)
  res_transects_o3_10 <- observe_individuals_along_transects_2(distance_matrix, Ntransects)
  res_transects_intersecting_o3_10 <- observe_individuals_intersecting_homerange_transects(raster_end10, result$home_ranges10, obs_list$obs_transects10_o3,  Nindividuals_end, Ntransects)
  #print('o3')
  
  # o4:
  
  obs_within_squares_o4_10 <- allocate_locations_within_squares(loc_xy10$individual_locations_xy, obs_list$obs_cells10_o4 )
  res_squares_intersecting_o4_10 <- observe_individuals_intersecting_homerange_squares(raster_end10, result$home_ranges10, obs_list$obs_cells10_o4)
  Nsquares  <- length(obs_list$obs_cells10_o4)
  Ntransects  <- length(obs_list$obs_transects10_o4)
  res_o4_10 <- observe_individuals_within_squares(obs_within_squares_o4_10$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy10$individual_locations_xy, obs_list$obs_transects10_o4)
  res_transects_o4_10 <- observe_individuals_along_transects_2(distance_matrix, Ntransects)
  res_transects_intersecting_o4_10 <- observe_individuals_intersecting_homerange_transects(raster_end10, result$home_ranges10, obs_list$obs_transects10_o4,  Nindividuals_end, Ntransects)
  #print('o4')
  
  # o5:
  
  obs_within_squares_o5_10 <- allocate_locations_within_squares(loc_xy10$individual_locations_xy, obs_list$obs_cells10_o5 )
  res_squares_intersecting_o5_10 <- observe_individuals_intersecting_homerange_squares(raster_end10, result$home_ranges10, obs_list$obs_cells10_o5)
  Nsquares  <- length(obs_list$obs_cells10_o5)
  Ntransects  <- length(obs_list$obs_transects10_o5)
  res_o5_10 <- observe_individuals_within_squares(obs_within_squares_o5_10$locations_per_square, observation_probability, Nsquares)
  distance_matrix <-compute_distances_to_transects(loc_xy10$individual_locations_xy, obs_list$obs_transects10_o5)
  res_transects_o5_10 <- observe_individuals_along_transects_2(distance_matrix, Ntransects)
  res_transects_intersecting_o5_10 <- observe_individuals_intersecting_homerange_transects(raster_end10, result$home_ranges10, obs_list$obs_transects10_o5,  Nindividuals_end, Ntransects)
  #print('o5')
  
  Nindividuals <- length(loc_xy10$individual_locations[,1])
  
  par(mfrow =c(2, 3))
  
  library(broman)
  springgreen <- brocolors("crayons")["White"]
  maroon <- brocolors("crayons")["Asparagus"]
  sienna <- brocolors("crayons")["Caribbean Green"]
  asparagus <- brocolors("crayons")["Pine Green"]
  sienna <- brocolors("crayons")["Mountain Meadow"]
  
  png(paste(label, 'l6.png'), height = 6, width = 6, units = 'in', res=300)
  plot(as.factor(raster_end6),axes=FALSE, box=FALSE,  legend=FALSE, col = c('white', sienna))
  dev.off()
  
  png(paste(label, 'l7.png'), height = 6, width = 6, units = 'in', res=300)
  plot(as.factor(raster_end7),axes=FALSE, box=FALSE,  legend=FALSE, col = c('white', sienna))
  dev.off()
  
  png(paste(label, 'l8.png'), height = 6, width = 6, units = 'in', res=300)
  plot(as.factor(raster_end8),axes=FALSE, box=FALSE,  legend=FALSE, col = c('white', sienna))
  dev.off()
  
  png(paste(label, 'l9.png'), height = 6, width = 6, units = 'in', res=300)
  plot(as.factor(raster_end9),axes=FALSE, box=FALSE,  legend=FALSE, col = c('white', sienna))
  dev.off()
  
  png(paste(label, 'l10.png'), height = 6, width = 6, units = 'in', res=300)
  plot(as.factor(raster_end10),axes=FALSE, box=FALSE,  legend=FALSE, col = c('white', sienna))
  dev.off()
  
  png(paste(label, 'i6.png'), height = 6, width = 6, units = 'in', res=300)
  plot_home_ranges_nocolor2(result$home_ranges6, raster_end6, length(result$home_ranges6))
  dev.off()
  
  png(paste(label, 'i7.png'), height = 6, width = 6, units = 'in', res=300)
  plot_home_ranges_nocolor2(result$home_ranges7, raster_end7, length(result$home_ranges7))
  dev.off()
  
  png(paste(label, 'i8.png'), height = 6, width = 6, units = 'in', res=300)
  plot_home_ranges_nocolor2(result$home_ranges8, raster_end8, length(result$home_ranges8))
  dev.off()
  
  png(paste(label, 'i9.png'), height = 6, width = 6, units = 'in', res=300)
  plot_home_ranges_nocolor2(result$home_ranges10, raster_end9, length(result$home_ranges10))
  #points(result$individual_cells9 )
  dev.off()
  
  png(paste(label, '_only_territories_i9.png'), height = 6, width = 6, units = 'in', res=300)
  plot_home_ranges_nocolor4(result$home_ranges10, raster_end9, length(result$home_ranges10), result$individual_cells9)
  #points(result$individual_cells9 )
  dev.off()
  
  png(paste(label, 'i9.png'), height = 6, width = 6, units = 'in', res=300)
  plot_home_ranges_nocolor2(result$home_ranges10, raster_end9, length(result$home_ranges10))
  #points(result$individual_cells9 )
  dev.off()
  
  png(paste(label, 'i10.png'), height = 6, width = 6, units = 'in', res=300)
  plot_home_ranges_nocolor2(result$home_ranges10, raster_end10, length(result$home_ranges10))
  dev.off()
  
  
  loc_xy6 <- set_individual_locations(raster_end6, result$home_ranges6, Nhabitat, length(result$home_ranges6), always_in_habitat, N_locations, habitat_dependency)
  loc_xy7 <- set_individual_locations(raster_end7, result$home_ranges7, Nhabitat, length(result$home_ranges7), always_in_habitat, N_locations, habitat_dependency)
  loc_xy8 <- set_individual_locations(raster_end8, result$home_ranges8, Nhabitat, length(result$home_ranges8), always_in_habitat, N_locations, habitat_dependency)
  loc_xy9 <- set_individual_locations(raster_end9, result$home_ranges9, Nhabitat, length(result$home_ranges9), always_in_habitat, N_locations, habitat_dependency)
  loc_xy10 <- set_individual_locations(raster_end10, result$home_ranges10, Nhabitat, length(result$home_ranges10), always_in_habitat, N_locations, habitat_dependency)
  
  res_6 <- observe_individuals_intersecting_homerange_squares(raster_end6, result$home_ranges6, obs_list$obs_cells10_o2)
  res_8 <- observe_individuals_intersecting_homerange_squares(raster_end10, result$home_ranges10, obs_list$obs_cells10_o2)
  res_10 <- observe_individuals_intersecting_homerange_squares(raster_end10, result$home_ranges10, obs_list$obs_cells10_o2)
  res_9 <- observe_individuals_intersecting_homerange_squares(raster_end9, result$home_ranges9, obs_list$obs_cells10_o2)
  res_7 <- observe_individuals_intersecting_homerange_squares(raster_end7, result$home_ranges7, obs_list$obs_cells10_o2)
  
  dt <- data.frame('popsize' = NULL, 'census' = NULL)
  dt[1,c('popsize', 'census')] <- c(length(result$home_ranges6), sum(res_6$presence_observed))
  dt[2,c('popsize', 'census')] <- c(length(result$home_ranges7), sum(res_7$presence_observed))
  dt[3,c('popsize', 'census')] <- c(length(result$home_ranges8), sum(res_8$presence_observed))
  dt[4,c('popsize', 'census')] <- c(length(result$home_ranges9), sum(res_9$presence_observed))
  dt[5,c('popsize', 'census')] <- c(length(result$home_ranges10), sum(res_10$presence_observed))
  
  res9 <- dt
  
  populationColor <- brocolors("crayons")["Purple Heart"]
  censusColor <- brocolors("crayons")["Mango Tango"]
  res9$prop_census <- res9$census/50
  head(res9)
  res9$year <- seq(1,5)
  coeff <- 150
  
  png(paste(label, 'PopCensusLines.png'), height = 3, width = 3.5, units = 'in', res =300)
  ggplot(res9, aes(x=year)) +
    
    geom_line( aes(y=prop_census), color=populationColor) + 
    geom_line( aes(y= popsize / coeff),  color=censusColor) + # Divide by 10 to get the same range than the temperature
    
    scale_y_continuous(
      name = "Population census",
      sec.axis = sec_axis(~.*coeff, name="Population size")
    )     + theme_classic()  +  theme(
      axis.title.y = element_text(color = populationColor, size=13),
      axis.title.y.right = element_text(color = censusColor, size=13)
    )   
  dev.off()
  
  
  dev.off()
  
  png(paste(label, 'o6.png'), height = 6, width = 6, units = 'in', res=300)
  plot_observations_intersecting_squares(raster_end10, result$home_ranges6, obs_list$obs_cells10_o2, res_6$presence_observed)
  dev.off()
  
  png(paste(label, 'o7.png'), height = 6, width = 6, units = 'in', res=300)
  plot_observations_intersecting_squares(raster_end7, result$home_ranges7, obs_list$obs_cells10_o2, res_7$presence_observed)
  dev.off()
  
  png(paste(label, 'o8.png'), height = 6, width = 6, units = 'in', res=300)
  plot_observations_intersecting_squares(raster_end10, result$home_ranges8, obs_list$obs_cells10_o2, res_8$presence_observed)
  dev.off()
  
  png(paste(label, 'o9.png'), height = 6, width = 6, units = 'in', res=300)
  plot_observations_intersecting_squares(raster_end9, result$home_ranges9, obs_list$obs_cells10_o2, res_9$presence_observed)
  dev.off()
  
  png(paste(label, 'o10.png'), height = 6, width = 6, units = 'in', res=300)
  plot_observations_intersecting_squares(raster_end10, result$home_ranges10, obs_list$obs_cells10_o2, res_10$presence_observed)
  dev.off()
  
  ppp <- as.ppp(loc_xy$individual_locations_xy, c(0,100,0,100))
  
  #lines(obs_cells10_o2)
  
  return(dt)
}



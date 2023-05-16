# testing how the functions work together. 
rm(list = ls())
setwd("~/packs")
arg <- commandArgs(trailingOnly = TRUE)
set.seed(arg)
###############
# packages
.libPaths("~/packs")

#setwd("/home/local/sroto/oa/code/final")
library(PrevMap)
library(MASS)
library(spatstat)
library(spdep)
library(ggplot2)
library(Rfast)
library(fields)
library(raster)
library(gdistance)
library(rasterKernelEstimates)
library(spdep)
library(rgdal)
library(rgdal)
library(rgeos)
library(dplyr)
library(spdep)
library(R.utils)
library(SDMTools)
library(spatstat)
library(SpaDES)
library(ff)
library(ffbase)

setwd("~/oa/oa/oa2")
#setwd("/home/local/sroto/oa")
source("habitat_functions_3.R")
source("ecology_functions_3.R")
source("home_range_functions_3.R")
source("fitness_functions_3.R")
source("dispersal_functions_3.R")
source('observation_functions_3.R')

# 1. create initial habitat:

n_x <- 150
n_y <- 150
Nyears <- 10

xy <- set_coordinates(n_x, n_y)
N_xy <- n_x*n_y
W <- matrix(1,3,3)

#DM <- Dist(xy, method = "euclidean", square = TRUE, p = 0,vector = FALSE)
#AM <- generate_adjacency_matrix(DM, treshold= 1) # needed to model the habitat change

# super big matrices, should we tune this? 


# 2.  Combine habitat change with species dynamics 
# first: set initial observation locations based on where the habitat is located

occupancy_level <- numeric(0)
population_size <- numeric(0)

# **** set observation schemes *****

# Q: check if the number od these matters (and if they should be in habitat throughout the simulation)

N_observation_locs <- 100
N_observation_cells <- N_observation_locs
Ntransects <- N_observation_locs

Nyears <- 10
Nhabitat <- 6 # check the exact definition of this parameter. 
N_locations <-1 # (these refer to individuals' behavior)
habitat_dependency  <- FALSE  #check the exact definition of this parameter. 

habitat_scenario_id <- 1

###########

df_sim <-data.frame() 

# ECOLOGICAL PROFILES:

hr_v <-    c('constant_size', 'habitat',  'nonoverlapping_core', 'nonoverlapping_even')
fitness_v <- c('neutral', 'habitat_density', 'habitat_inverse_density', 'habitat_mid_density', 'habitat_per_consumer_density','habitat_per_consumer_inverse_density', 'habitat_per_consumer_mid_density')
dispersal_v <-  c('passive', 'conspecific', 'mid_density', 'habitat_driven', 'conspecific_avoidance_driven') 
fitness_assumption_v <- c('stable')


initial_habitat_v <- c('scattered', 'metapopulation', 'solid') 
habitat_change_v <- c('scattered_decay', 'scattered_growth', 'clustered_decay', 'clustered_growth', 'no_change') 


if (FALSE){
hr_v <-    c('constant_size')
fitness_v <- c('neutral')
dispersal_v <-  c('passive') 
fitness_assumption_v <- c('stable')
}


habitat_v <- c(9) # (habitat requirements)
k_values <- c(2) # (dispersal) (go through these values, see how affects )
max_diameter_values <- 3

Nsims <- 8  # how many different landscape change scenarios
Nreps <- 2  # how many repetitions (of all ecological profiles per scenario)

counter <- 1
hroverlap5_list <- list()
hroverlap10_list <- list()
hrsize5_list <- list()
hrsize10_list <- list()
nn1_list <- list()
nn1_list10 <- list()

time_d <- numeric(0)
sigma2 <- 1
hab1 <- generate_initial_habitat_gp(n_x, n_y, sigma2, 0.5)
hab2 <- generate_initial_habitat_gp(n_x, n_y, sigma2, 2)
hab3 <- generate_initial_habitat_gp(n_x, n_y, sigma2, 15)



print('Habitats created')
#setwd("/home/local/sroto/oa")
#load('ExampleHabitats.Rdata')


initial_habitat <- hab1
k_values <- 2
k <- k_values 
max_diameter <- max_diameter_values 

max_dispersal_distance  <- 10
xy_cells <- xyFromCell(initial_habitat, cell= seq(1, n_x*n_y))


HabitatPerception <- matrix(c( NA, NA,  1, 1, 1,  NA, NA,  
                               NA,  1,  1, 1,  1,  1, NA,  
                               1,  1,  1, 1,  1,  1, 1, 
                               1,  1,  1, 0,  1,  1, 1,  
                               1,  1,  1, 1,  1,  1, 1, 
                               NA,  1,  1, 1,  1,  1, NA,  
                               NA, NA, 1 , 1,  1,  NA, NA), ncol=7, byrow=TRUE)
all_cells <- seq(1, n_x*n_y)

# Create neighbour lists used for dispersal
# NB is for the possible dispersal
# NB_perception is the radius for evaluating the dispersal destination. 

NB <- knn2nb(knearneigh(xy_cells, k = 4*(max_dispersal_distance)^2)) # see if this is used in the end
NB_perception <- dnearneigh(xy_cells, d1 = 0, d2 = 3)

dsts <-nbdists(NB, xy_cells)   
dsts_perception <-nbdists(NB_perception , xy_cells)  
square_size <- 4 # the lenght of one side of the square

print('distances computed!')

setwd("~/oa/oa/simresults")
  
for (reps in seq(1, Nreps)){
  
 for (h in initial_habitat_v){
      if (h== 'scattered'){
        hab <- hab1}
      if (h== 'metapopulation'){
        hab <- hab2}
      if (h== 'solid'){
        hab <- hab3}
      initial_habitat <- hab
      initial_habitat_df <- generate_cells_from_rasters(hab, xy)
      habitat_df <- generate_habitat_df(xy[,1], xy[,2], initial_habitat_df$habitat, Nyears)
      
      sum_habitat <- rasterLocalSums(initial_habitat, W)
      
      cluster_mean_size <- sample(seq(1, 30), 1) 
      lambda <- runif(1, min = 0.0001, max = 0.005)    
      cluster_mean_sizeNew <- sample(seq(1, 30), 1) 
      lambdaNew <- runif(1, min = 0.0001, max = 0.005) 
      
  for (hh in habitat_change_v){
  
  if (hh =='no_change'){
  habitat_df_p9 <- habitat_df 
  }else{
   # set parameters based on change type:
    
  ###############
  if (hh =='scattered_decay'){
  if (h=='scattered'){
    lambda <- 0.02
    lambdaNew <- 0
    p <- 0.8
    pNew <- 0.5
    d <- 1
    dNew <- 10
    }
    if (h=='metapopulation'){
      lambda <- 0.02
      lambdaNew <- 0
      p <- 0.8
      pNew <- 0.5
      d <- 1
      dNew <- 10
    }
    if (h=='solid'){
      lambda <- 0.02
      lambdaNew <- 0
      p <- 0.8
      pNew <- 0.5
      d <- 1
      dNew <- 10
    }
  }
  if (hh=='scattered_growth'){
    if (h=='scattered'){
      lambda <- 0
      lambdaNew <- 0.02
      p <- 0.8
      pNew <- 0.8
      d <- 1
      dNew <- 1
    }
    if (h=='metapopulation'){
      lambda <- 0
      lambdaNew <- 0.02
      p <- 0.8
      pNew <- 0.8
      d <- 1
      dNew <- 1
    }
    if (h=='solid'){
      lambda <- 0
      lambdaNew <- 0.02
      p <- 0.8
      pNew <- 0.8
      d <- 1
      dNew <- 1
    }
  }
  if (hh== 'clustered_decay'){
    if (h=='scattered'){
      lambda <- 0.005
      lambdaNew <- 0
      p <- 0.8
      pNew <- 0.5
      d <- 8
      dNew <- 10
    }
    if (h=='metapopulation'){
      lambda <- 0.005
      lambdaNew <- 0
      p <- 0.7
      pNew <- 0.5
      d <- 8
      dNew <- 10
    }
    if (h=='solid'){
      lambda <- 0.005
      lambdaNew <- 0
      p <- 0.6
      pNew <- 0.5
      d <- 8
      dNew <- 10
    }
    }
  if (hh=='clustered_growth'){
    if (h=='scattered'){
      lambda <- 0
      lambdaNew <- 0.005
      p <- 0.8
      pNew <- 0.8
      d <- 8
      dNew <- 8
    }
    if (h=='metapopulation'){
      lambda <- 0
      lambdaNew <- 0.005
      p <- 0.7
      pNew <- 0.7
      d <- 8
      dNew <- 8
    }
    if (h=='solid'){
      lambda <- 0
      lambdaNew <- 0.005
      p <- 0.6
      pNew <- 0.5
      d <- 8
      dNew <- 10
    }
  }
  ###############  
  habitat_df_p9 <- generate_habitat_change_fires(habitat_df, hab, 
                                      Nyears, N_xy, 
                                      lambda, lambdaNew, 
                                      p, pNew,
                                      d, dNew)
  }
    
  
    
  Nindividuals <- sample(seq(25, 800), 1)
  #Nindividuals <- 100
  print('simulation')
  print(reps)
  xy_individuals<- sample_individuals_to_landscape(initial_habitat, Nindividuals)
  individual_cells <-cellFromXY(initial_habitat, xy_individuals)
  
  for (max_diameter in max_diameter_values){
  
  # set the habitat perception matrix based on the home range sizes. 
  Wr <- focalWeight(initial_habitat, max_diameter, type = "Gauss")
  pp  <- get_individual_densities(initial_habitat, xy_individuals, Wr)
  
  
  raster5 <- generate_raster_from_cells(initial_habitat, habitat_df_p9, 5)
  stats5 <- PatchStat(raster5)
  raster10 <- generate_raster_from_cells(initial_habitat, habitat_df_p9, 10)
  stats10 <- PatchStat(raster10)

  perim.area.ratio5 <- stats5$perim.area.ratio[stats5$patchID==1]
  perim.area.ratio10 <- stats10$perim.area.ratio[stats10$patchID==1]
  perim.core.area5 <- stats5$core.area.index[stats5$patchID==1]
  perim.core.area10 <- stats10$core.area.index[stats10$patchID==1]
  habitatcells5 <- stats5$n.cell[stats5$patchID==1]
  habitatcells10 <- stats10$n.cell[stats10$patchID==1]
  print(habitatcells10 - habitatcells5)

  # systematic sampling, permanent plots (o1)
  
  # obs_cells5_o1 <- set_observation_squares(raster5, square_size, selection_type, N_observation_locs, habitat_cells, N_observation_cells, 5, FALSE, TRUE)
  # obs_cells10_o1 <-obs_cells5_o1
  
  N_observation_locs_1 <- 50
  
  obs_cells5_o1 <- set_observation_squares(raster5, square_size, selection_type, N_observation_locs_1, habitat_cells, N_observation_locs_1, 10, TRUE, FALSE)
  obs_cells10_o1 <- set_observation_squares(raster10, square_size, selection_type, N_observation_locs_1, habitat_cells, N_observation_locs_1, 10, TRUE, FALSE)
  
  obs_transects5_o1 <- set_observation_transects(raster5, square_size, selection_type, ceiling(N_observation_locs/16), habitat_cells, N_observation_cells, 5, FALSE, TRUE)
  obs_transects10_o1 <-obs_transects5_o1 
  
  #'permanent_plots_random_sampling_small_d' (o2)
  # obs_cells5_o2 <- set_observation_squares(raster5, square_size, selection_type, N_observation_locs, habitat_cells, N_observation_cells, 10, FALSE, FALSE)
  # obs_cells10_o2 <-obs_cells5_o2
  
  N_observation_locs_2 <- 75  
  obs_cells5_o2 <- set_observation_squares(raster5, square_size, selection_type, N_observation_locs_2, habitat_cells, N_observation_locs_2, 10, TRUE, FALSE)
  obs_cells10_o2 <- set_observation_squares(raster10, square_size, selection_type, N_observation_locs_2, habitat_cells, N_observation_locs_2, 10, TRUE, FALSE)
  
  # 'dynamic_plots_habitat_random_small_d' (o4)
  N_observation_locs_4 <- 100 
  obs_cells5_o4 <- set_observation_squares(raster5, square_size, selection_type, N_observation_locs_4, habitat_cells, N_observation_locs_4, 9, TRUE, FALSE)
  obs_cells10_o4 <- set_observation_squares(raster10, square_size, selection_type, N_observation_locs_4, habitat_cells, N_observation_locs_4, 9, TRUE, FALSE)

  obs_list <- list('obs_cells5_o1' = obs_cells5_o1, 'obs_cells10_o1' = obs_cells10_o1,
                   'obs_transects5_o1' = obs_transects5_o1, 'obs_transects10_o1' = obs_transects10_o1,
                   'obs_cells5_o2' = obs_cells5_o2, 'obs_cells10_o2' = obs_cells10_o2,
                   'obs_cells5_o4' = obs_cells5_o4, 'obs_cells10_o4' = obs_cells10_o4)
  
  

    for (k in k_values){ 
      for(nh in habitat_v){
        for (fitness in fitness_v){
          for (dispersal in dispersal_v){
            for (hr in hr_v){
              for (fitness_assumption in fitness_assumption_v){
                
                print(paste(dispersal, hr, fitness))
                
                if (((hr == 'nonoverlapping_uneven') | (hr =='nonoverlapping_even')) & ((fitness =='density_induced') | (fitness =='inverse_density'))){
                  print('Dont simulate this!')
                }else{

                res9 <- try(simulate_observations_environmental_change(initial_habitat, habitat_df_p9, individual_cells, 10, TRUE, hr, fitness, dispersal, NB, dsts, max_diameter, nh, Wr, obs_list, fitness_assumption, NB_perception, dsts_perception))

               if (class(res9) != "try-error"){
                 
                 df_sim[counter,'cluster_mean_size'] <- cluster_mean_size
                 df_sim[counter,'lambda'] <-lambda
                 
                 # population size and observations at year 5:
                 df_sim[counter,'popsize'] <-res9$pop_size 
                 df_sim[counter,'occupancy_o1'] <-sum(res9$presence_observed_o1)/N_observation_locs_1
                 df_sim[counter,'occupancy_o2'] <-sum(res9$presence_observed_o2)/N_observation_locs_2
                 df_sim[counter,'occupancy_o4'] <-sum(res9$presence_observed_o4)/N_observation_locs_4
                 
                 df_sim[counter,'abundance_o1'] <-sum(res9$number_of_tracks_observed_o1)
                 df_sim[counter,'abundance_o2'] <-sum(res9$number_of_tracks_observed_o2)
                 df_sim[counter,'abundance_o4'] <-sum(res9$number_of_tracks_observed_o4)
                 
                 df_sim[counter,'individuals_o1'] <-sum(res9$number_individuals_observed_intersecting_squares_o1)
                 df_sim[counter,'individuals_o2'] <-sum(res9$number_individuals_observed_intersecting_squares_o2)
                 df_sim[counter,'individuals_o4'] <-sum(res9$number_individuals_observed_intersecting_squares_o4)
                 
                 df_sim[counter,'landscape_exploitation_o1'] <-sum(res9$landscape_exploitation_o1)
                 df_sim[counter,'landscape_exploitation_o2'] <-sum(res9$landscape_exploitation_o2)
                 df_sim[counter,'landscape_exploitation_o4'] <-sum(res9$landscape_exploitation_o4)
                 
                 df_sim[counter,'abundance_transects_o1'] <-sum(res9$number_of_tracks_observed_transects_o1)
                 df_sim[counter,'abundance_transects_o1_10'] <-sum(res9$number_of_tracks_observed_transects_o1_10)
                 
                 df_sim[counter,'occupancy_transects_o1'] <-sum(res9$presence_observed_transects_o1)
                 ## the following line was wrong:
                 df_sim[counter,'occupancy_transects_o1_10'] <-sum(res9$presence_observed_transects_o1_10) 
                 
                 df_sim[counter,'occupancy_intersecting_o1'] <-sum(res9$presence_observed_intersecting_squares_o1)/N_observation_locs_1
                 df_sim[counter,'occupancy_intersecting_o2'] <-sum(res9$presence_observed_intersecting_squares_o2)/N_observation_locs_2
                 df_sim[counter,'occupancy_intersecting_o4'] <-sum(res9$presence_observed_intersecting_squares_o4)/N_observation_locs_4

                 df_sim[counter,'occupancy_transects_intersecting_o1'] <-sum(res9$presence_observed_intersecting_transects_o1)
                 df_sim[counter,'occupancy_transects_intersecting_o1_10'] <-sum(res9$presence_observed_intersecting_transects_o1_10)
                 
                 
                 
                 
                # # population size and observations at year 10:
                 df_sim[counter,'popsize10'] <-res9$pop_size10 
                 
                 df_sim[counter,'occupancy_o1_10'] <-sum(res9$presence_observed_o1_10)/N_observation_locs_1
                 df_sim[counter,'occupancy_o2_10'] <-sum(res9$presence_observed_o2_10)/N_observation_locs_2
                 df_sim[counter,'occupancy_o4_10'] <-sum(res9$presence_observed_o4_10)/N_observation_locs_4


                 
                 df_sim[counter,'occupancy_intersecting_o1_10'] <-sum(res9$presence_observed_intersecting_squares_o1_10)/N_observation_locs_1
                 df_sim[counter,'occupancy_intersecting_o2_10'] <-sum(res9$presence_observed_intersecting_squares_o2_10)/N_observation_locs_2
                 df_sim[counter,'occupancy_intersecting_o4_10'] <-sum(res9$presence_observed_intersecting_squares_o4_10)/N_observation_locs_4

                 

                 df_sim[counter,'abundance_o1_10'] <-sum(res9$number_of_tracks_observed_o1_10)
                 df_sim[counter,'abundance_o2_10'] <-sum(res9$number_of_tracks_observed_o2_10)
                 df_sim[counter,'abundance_o4_10'] <-sum(res9$number_of_tracks_observed_o4_10)

                 df_sim[counter,'individuals_o1_10'] <-sum(res9$number_individuals_observed_intersecting_squares_o1_10)
                 df_sim[counter,'individuals_o2_10'] <-sum(res9$number_individuals_observed_intersecting_squares_o2_10)
                 df_sim[counter,'individuals_o4_10'] <-sum(res9$number_individuals_observed_intersecting_squares_o4_10)


                 df_sim[counter,'landscape_exploitation_o1_10'] <-sum(res9$landscape_exploitation_o1_10)
                 df_sim[counter,'landscape_exploitation_o2_10'] <-sum(res9$landscape_exploitation_o2_10)
                 df_sim[counter,'landscape_exploitation_o4_10'] <-sum(res9$landscape_exploitation_o4_10)

                 df_sim[counter,'hr'] <-hr
                 df_sim[counter,'fitness'] <-fitness
                 df_sim[counter,'fitness_assumption'] <-fitness_assumption
                 
                 fitnesses <- scale(res9$fitness10, center = FALSE)
                 home_range_sizes <- unlist(lapply(res9$home_ranges10, length))
                 fitnesses5 <- scale(res9$fitness, center = FALSE)
                 home_range_sizes5 <- unlist(lapply(res9$home_ranges, length))
                 
                 df_sim[counter,'mean_fitness'] <- mean(fitnesses, na.rm=T)
                 df_sim[counter,'var_fitness'] <- var(fitnesses, na.rm=T)
                 df_sim[counter,'fitness_q1'] <- quantile(fitnesses, 0.05)
                 df_sim[counter,'fitness_q2'] <- quantile(fitnesses, 0.95)
                 df_sim[counter,'mean_hrs'] <- mean(home_range_sizes, na.rm=T)
                 df_sim[counter,'var_hrs'] <- var(home_range_sizes, na.rm=T)
                 df_sim[counter,'hrs_q1'] <- quantile(home_range_sizes, 0.05)
                 df_sim[counter,'hrs_q2'] <- quantile(home_range_sizes, 0.95)
                 
                 df_sim[counter,'mean_fitness5'] <- mean(fitnesses5, na.rm=T)
                 df_sim[counter,'var_fitness5'] <- var(fitnesses5, na.rm=T)
                 df_sim[counter,'fitness_q15'] <- quantile(fitnesses5, 0.05)
                 df_sim[counter,'fitness_q25'] <- quantile(fitnesses5, 0.95)
                 df_sim[counter,'mean_hrs5'] <- mean(home_range_sizes5, na.rm=T)
                 df_sim[counter,'var_hrs5'] <- var(home_range_sizes5, na.rm=T)
                 df_sim[counter,'hrs_q15'] <- quantile(home_range_sizes5, 0.05)
                 df_sim[counter,'hrs_q25'] <- quantile(home_range_sizes5, 0.95)

                 df_sim[counter,'dispersal'] <-dispersal
                 df_sim[counter,'sim_id'] <-counter
                 df_sim[counter,'lambda'] <- lambda
                 df_sim[counter,'cluster_mean_size'] <-cluster_mean_size
                 df_sim[counter,'perim.area.ratio5'] <-perim.area.ratio5
                 df_sim[counter,'perim.area.ratio10'] <-perim.area.ratio10
                 df_sim[counter,'perim.core.area5'] <-perim.core.area5
                 df_sim[counter,'perim.core.area10']<-perim.core.area10
                 df_sim[counter,'habitatcells5'] <-habitatcells5
                 df_sim[counter,'habitatcells10'] <- habitatcells10
                 df_sim[counter,'coveredarea5'] <- length(unique(unlist(res9$home_ranges)))
                 df_sim[counter,'coveredarea10'] <-length(unique(unlist(res9$home_ranges10)))
                 df_sim[counter,'max_diameter'] <- max_diameter
                 df_sim[counter,'rep'] <- reps
                 df_sim[counter,'initial_habitat'] <- h
                 df_sim[counter,'habitat_change'] <- hh
                 df_sim[counter,'habitat_scenario_id'] <- habitat_scenario_id 
                 hrsize5_list[[counter]] <- unlist(lapply(res9$home_ranges, length))
                 hrsize10_list[[counter]] <- unlist(lapply(res9$home_ranges10, length))
                 hroverlap5_list[[counter]] <- get_home_range_overlaps(res9$home_ranges)
                 hroverlap10_list[[counter]] <-get_home_range_overlaps(res9$home_ranges10)
                 nn1_list[[counter]] <- res9$nn1
                 nn1_list10[[counter]] <- res9$nn110
		             counter <- counter+1
              }else{
                
              }
             # counter <- counter+1
  #            print(paste('separate simulations', counter))
              }
            }
          }
        }
      }
    }
  }
  }
  }
      habitat_scenario_id <- habitat_scenario_id + 1
  }
 # print(paste('independent simulations', nn, 'parallel', arg))
   filename <- paste('simdata_new', reps, arg, '.Rdata', collapse = "") 
   save(file = filename, df_sim, 
        hroverlap5_list,
        hroverlap10_list, 
        hrsize5_list,
        hrsize10_list,
        nn1_list,
        nn1_list10) 
}



setwd("~/oa/oa/simresults")
filename <- paste('simdata_new', arg, '.Rdata', collapse = "") 

save(file = filename, df_sim, 
     hroverlap5_list,
     hroverlap10_list, 
     hrsize5_list,
     hrsize10_list,
     nn1_list,
     nn1_list10)  

names(df_sim)


ggplot(df_sim, aes(occupancy_o1_10, occupancy_intersecting_o1_10))+geom_point()
ggplot(df_sim, aes(occupancy_o1_10, individuals_o1_10))+geom_point()
ggplot(df_sim, aes(occupancy_o1_10, landscape_exploitation_o1))+geom_point()
ggplot(df_sim, aes(abundance_transects_o1_10, occupancy_transects_o1_10))+geom_point()
ggplot(df_sim, aes(abundance_transects_o1_10, occupancy_transects_intersecting_o1_10))+geom_point()


"occupancy_o1"                          
[5] "occupancy_o2"                           "occupancy_o4"                          
[7] "abundance_o1"                           "abundance_o2"                          
[9] "abundance_o4"                           "individuals_o1"                        
[11] "individuals_o2"                         "individuals_o4"                        
[13] "landscape_exploitation_o1"              "landscape_exploitation_o2"             
[15] "landscape_exploitation_o4"              "abundance_transects_o1"                
[17] "abundance_transects_o1_10"              "occupancy_transects_o1"                
[19] "occupancy_transects_o1_10"              "occupancy_intersecting_o1"             
[21] "occupancy_intersecting_o2"              "occupancy_intersecting_o4"             
[23] "occupancy_transects_intersecting_o1"    "occupancy_transects_intersecting_o1_10"
[25] "popsize10"                              "occupancy_o1_10"                       
[27] "occupancy_o2_10"                        "occupancy_o4_10"                       
[29] "occupancy_intersecting_o1_10"           "occupancy_intersecting_o2_10"          
[31] "occupancy_intersecting_o4_10"           "abundance_o1_10"                       
[33] "abundance_o2_10"                        "abundance_o4_10"                       
[35] "individuals_o1_10"                      "individuals_o2_10"                     
[37] "individuals_o4_10"                      "landscape_exploitation_o1_10"          
[39] "landscape_exploitation_o2_10"           "landscape_exploitation_o4_10" 

# ggplot(df_sim, aes(occupancy_o1_10, popsize10))+geom_point()
# ggplot(df_sim, aes(occupancy_o2_10, popsize10))+geom_point()
# ggplot(df_sim, aes(occupancy_intersecting_o2_10, popsize10))+geom_point()
# ggplot(df_sim, aes(abundance_o2_10, popsize10))+geom_point()
# 
# ggplot(df_sim, aes(occupancy_transects_o1_10, popsize10))+geom_point()
# ggplot(df_sim, aes(occupancy_transects_intersecting_o1_10, popsize10))+geom_point()
# ggplot(df_sim, aes(abundance_transects_o1_10, popsize10))+geom_point()
# ggplot(df_sim, aes(occupancy_transects_o1_10, popsize10))+geom_point()
# ggplot(df_sim, aes(individuals_o2_10, popsize10))+geom_point()
# ggplot(df_sim, aes(landscape_exploitation_o2_10, popsize10))+geom_point()
# ggplot(df_sim, aes(popsize, popsize10))+geom_point()
# 
# df_sim$ecology <- interaction(interaction(df_sim$hr, df_sim$fitness), df_sim$dispersal)
# ggplot(df_sim[which(df_sim$fitness!='neutral'), ], aes(ecology, mean_fitness))+geom_point()+coord_flip()
# 
# 
# dd <- data.frame('hr_size' = NULL, 'ecology'= NULL, 'hr' = NULL)
# for (i in seq(1, length(df_sim))){
#   dd <- rbind(dd, data.frame('hr_size' =hrsize10_list[[i]], 'ecology'= rep(df_sim$ecology[i], length(hrsize10_list[[i]])),  'hr'= rep(df_sim$hr[i], length(hrsize10_list[[i]])))) 
# }
# ggplot(dd, aes(x=hr, y = hr_size)) +geom_boxplot()+coord_flip()

#ggplot(df_sim, aes(occupancy_o2_10, abundance_o5_10))+geom_point()
# ggplot(df_sim, aes(occupancy_o1_10, occupancy_o5_10))+geom_point()
# ggplot(df_sim, aes(occupancy_o3_10, occupancy_o5_10))+geom_point()
# ggplot(df_sim, aes(occupancy_o4_10, occupancy_o5_10))+geom_point()



#df_sim$landscape_change <- df_sim$habitatcells10-df_sim$habitatcells5
#df_sim$perim.core.area_change <- df_sim$perim.core.area10-df_sim$perim.core.area5    
#df_sim$perim.area.ratio_change <- df_sim$perim.area.ratio10-df_sim$perim.area.ratio5    

# to do: think about growth and decay as functions of fitness. 
# especially w.r.t. density dependent fitness. 
# make a training run - think what should be monitored.         

setwd("~/oa/oa/simresults")




###################################

hist(hroverlap10_list[[1]])
head(df_sim)


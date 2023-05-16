
library(MASS)
library(ggplot2)
library(Rfast)
library(fields)
library(raster)
library(gdistance)
library(rasterKernelEstimates)

plot_habitat_time_series <- function(habitat_df_in, group.colors, lambda, cluster_mean_size){
  ggplot(habitat_df_in, aes(x, y, col = as.factor(habitat)))+theme_classic()+geom_point()+facet_wrap(.~year, ncol = 5)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+  theme(axis.title.y = element_blank(),
                                                axis.text.y =element_blank(),
                                                axis.ticks.y =element_blank()) +
    theme(strip.background = element_blank(), strip.text = element_blank())+ scale_color_manual(values=group.colors)+ guides(color=FALSE)+ggtitle(paste(lambda, cluster_mean_size))
}

set_coordinates <- function(n_x, n_y){
  # generates n_x* n_y points in a rectangular grid
  xx <- seq(1, n_x, by = 1)
  yy <- seq(1, n_y, by = 1)
  xy <- expand.grid(x = xx, y = yy)
  return(xy)
}

generate_distance_matrix <- function(x1, x2){
  
  # returns squared distance
  Ndim <- length(x1[,1])
  DM <- matrix(nrow = Ndim, ncol = Ndim)
  
  for (i in 1:(Ndim-1)) {
    for (j in (i+1):Ndim) {
      DM[i, j] <- sum((x1[i,]-x2[j,])^2)
      DM[j, i] <- DM[i, j]
    }
  }
  diag(DM) <- 0
  return(DM)
  
}

generate_adjacency_matrix <- function(DM, treshold){
  
  AM <- matrix(0, nrow = length(DM[,1]), ncol = length(DM[,1]))
  AM[sqrt(DM) <= treshold] <- 1
  diag(AM) <- 1
  return(AM)
  
}

Matern32Covariance <- function(x1,x2,l,sigma2, DM) {
  # K = Matern32Covariance(x1,x2,l,sigma2)
  # 
  # A function that takes matrices x1 and x2 of input locations, lengthscale
  # parameter l and magnitude parameter sigma2 and constructs a covariance
  # matrix between f(x1) and f(x2) corresponding to a Matern with degrees of 
  # freedom nu=3/2
  
  K=matrix(1,nrow(x1),nrow(x2));
  
  if (nrow(x1)==nrow(x2) && all(x1==x2) ) {
    # If symmetric matrix do fast way
    for (i in 1:(nrow(x1)-1)) {
      for (j in (i+1):nrow(x1)) {
        h = sum((DM[i,j])^2)/l^2
        K[i,j] = (1 + sqrt( 3*h ) )*exp(-sqrt( 3*h) )
        K[j,i] = K[i, j]
      }
    }
  }
  else {
    # If non-symmetric matrix
    for (j1 in 1:nrow(x1)){
      for (j2 in 1:nrow(x2)){
        h = sum((x1[j1,]-x2[j2,])^2)/l^2
        K[j1,j2] = (1 + sqrt( 3*h ) )*exp(-sqrt( 3*h) )
      }
    }
  }
  K = sigma2*K
  return(K)
}

set_habitat <- function(n, Kt){
#mu <- rep(0, n)
#cc  <- 1/(1+exp(-1*mvrnorm(n = 1, mu, Kt2, tol = 1e-6, empirical = FALSE)))

u <- rnorm(n, mean = 0, sd = 1)
mu <- rep(0, n)+Kt%*%u
cc  <- 1/(1+exp(-1*mu))

return(cc)
}

generate_initial_habitat <- function(n_x, n_y, sigma2, l){
  
  if (FALSE){
  n_x <- 100
  n_y <- 100
  sigma2 <- 1
  Nindividuals <- 200
  l <- 25
  }
  
  xy <- set_coordinates(n_x, n_y)
  N_xy <- n_x*n_y
  DM <- Dist(xy, method = "euclidean", square = TRUE, p = 0,vector = FALSE)
  AM <- generate_adjacency_matrix(DM, treshold= 5)
  Kt <- Matern32Covariance(xy, xy, l, sigma2, DM)
  habitat <- set_habitat(N_xy, Kt)
  habitat_treshold <- quantile(habitat, 0.5)
  habitat <- 1*(habitat>habitat_treshold)
  r <- raster(ncol=n_x, nrow=n_y, xmn=1, xmx=n_x, ymn=1, ymx=n_y)
  values(r) <- habitat
  initial_habitat <- r
  xy_cells <- xyFromCell(initial_habitat, cell= seq(1, n_x*n_y)) 
  NB <- knn2nb(knearneigh(xy_cells, k = 22)) # see if this is used in the end
  
  initial_habitat_df <- generate_cells_from_rasters(initial_habitat, xy)
  habitat_df <- generate_habitat_df(xy[,1], xy[,2], initial_habitat_df$habitat, Nyears)
  sum_habitat <- rasterLocalSums(initial_habitat, W)
  
  change_direction <- 'decay'
  change_type <- 'localized'
  
  cluster_mean_size <- sample(seq(1, 30), 1) 
  lambda <- runif(1, min = 0.0001, max = 0.005)    
  cluster_mean_sizeNew <- sample(seq(1, 30), 1) 
  lambdaNew <- runif(1, min = 0.0001, max = 0.005)  
  
  habitat_df_p9 <- generate_habitat_change_poissonpoint_process(habitat_df, change_type, Nyears, change_direction, N_xy, AM, cluster_mean_size, lambda, cluster_mean_sizeNew, lambdaNew)
  
  initial_habitat<- generate_raster_from_cells(initial_habitat, habitat_df_p9, 10)
  
  return(initial_habitat)
}

generate_initial_habitat_gp <- function(n_x, n_y, sigma2, l){
  
  if (FALSE){
    n_x <- 100
    n_y <- 100
    sigma2 <- 1
    Nindividuals <- 200
    l <- 25
  }
  
  xy <- set_coordinates(n_x, n_y)
  N_xy <- n_x*n_y
  DM <- Dist(xy, method = "euclidean", square = TRUE, p = 0,vector = FALSE)
  #AM <- generate_adjacency_matrix(DM, treshold= 5)
  Kt <- Matern32Covariance(xy, xy, l, sigma2, DM)
  habitat <- set_habitat(N_xy, Kt)
  habitat_treshold <- quantile(habitat, 0.5)
  habitat <- 1*(habitat>habitat_treshold)
  r <- raster(ncol=n_x, nrow=n_y, xmn=1, xmx=n_x, ymn=1, ymx=n_y)
  values(r) <- habitat
  initial_habitat <- r
  xy_cells <- xyFromCell(initial_habitat, cell= seq(1, n_x*n_y)) 
  #NB <- knn2nb(knearneigh(xy_cells, k = 22)) # see if this is used in the end
  return(r)
}

generate_raster_from_cells <- function(r, habitat_df, year){
  rr <- r 
  values(rr) <- 0
  inds_year <- which((habitat_df$year==year) & (habitat_df$habitat==1))
  habitat_cells <- cellFromXY(rr, cbind(habitat_df[inds_year, 'x'],habitat_df[inds_year, 'y']) )
  rr[habitat_cells] <- 1
  return(rr)
}

generate_cells_from_rasters <- function(r, xy){
  Nxy <- length(xy[,1])
  habitat_df <- data.frame('x' = seq(1,Nxy), 'y' =seq(1,Nxy), 'habitat' = seq(1,Nxy), 'cell_number' = seq(1,Nxy))
  cell_numbers <- cellFromXY(r, xy)
  habitat_df$x <-xy[,1]
  habitat_df$y <-xy[,2]
  habitat_df$habitat <-values(r) 
  habitat_df$cell_number <-cell_numbers
  return(habitat_df)
}

generate_habitat_df <- function(x, y, habitat, Nyears){
  # generates the habitat 10-year time-series, assuming no change in it. 
  df1<- data.frame('x' = xy[,1], 'y'= xy[,2], 'habitat' = habitat, 'year' = 1, 'loc_id' = seq(1, length(xy[,1])))
  for (y in seq(2, Nyears)){
    dft <- data.frame('x' = xy[,1], 'y'= xy[,2], 'habitat' = habitat, 'year' = y, 'loc_id' = seq(1, length(xy[,1])))
    df1 <- rbind(df1, dft)
  }
  return(df1)
}


generate_habitat_change_fires <- function(habitat_df, hab, 
                                          Nyears, N_xy, 
                                          lambda, lambdaNew, 
                                          p, pNew,
                                          d, dNew){
  # lambdas give the rate of decay or emergence of habitat, 
  # hab is the raster
  # p is the probability of fire spread
  # d is the max diameter of a fire. 
  mirror_hab <- 1*(hab==0)
  
  if (FALSE){
    habitat_df <- generate_habitat_df(xy[,1], xy[,2], initial_habitat_df$habitat, Nyears)
    lambda <- 0.01
    cluster_size <- 15
  }
  
  # generates the habitat 10-year time-series, assuming no change in it.
  suitable_habitat_cells <- which((habitat_df$habitat==1) & (habitat_df$year==1))
  initial_habitat <- sum(habitat_df[suitable_habitat_cells, 'habitat'])
 ########################### 
  # if (FALSE){
  #   ggplot(habitat_df, aes(x, y, col = as.factor(habitat)))+theme_classic()+geom_point()+facet_wrap(.~year, ncol = 5)+
  #     theme(axis.title.x=element_blank(),
  #           axis.text.x=element_blank(),
  #           axis.ticks.x=element_blank())+  theme(axis.title.y = element_blank(),
  #                                                 axis.text.y =element_blank(),
  #                                                 axis.ticks.y =element_blank()) +
  #     theme(strip.background = element_blank(), strip.text = element_blank())+ scale_color_manual(values=group.colors)+ guides(color=FALSE)
  # }
 ###########################
  
  # first habitat disappearance:
  suitable_habitat_cells_left <- suitable_habitat_cells 
  all_affected_habitat_cells <- NULL
  all_affected_nonhabitat_cells <- NULL
  
  for (i in seq(6, Nyears)){
    rr <- runif(1, min = 0, max = 1)
    if (rr>0.5){
      Ncellsleft <- length(suitable_habitat_cells_left)
      Naffected  <- rpois(1, lambda*Ncellsleft)
      Naffected <- min(Ncellsleft, Naffected)
      affected_habitat_cells <- sample(suitable_habitat_cells_left, Naffected, replace=FALSE)
      if (length(affected_habitat_cells)>0){
      fires <-  spread(hab, affected_habitat_cells, 
                       p*hab, 0, NULL, 1e8, 4, iterations = d, id = TRUE, returnIndices = TRUE)
      affected_habitat_cells <- union(affected_habitat_cells, fires$indices)
      inds_affected <- which((habitat_df$loc_id %in% affected_habitat_cells)&(habitat_df$year >=i)) 
      habitat_df$habitat[inds_affected] <- 0
      hab[inds_affected] <- 0
      mirror_hab[inds_affected] <- 1
      suitable_habitat_cells_left <- setdiff(suitable_habitat_cells_left, affected_habitat_cells)
      }
      
      suitable_nonhabitat_cells <- habitat_df$loc_id[which((habitat_df$habitat==0) & (habitat_df$year==i))]
      suitable_nonhabitat_cells_left <- suitable_nonhabitat_cells
      Ncellsleft <- length(suitable_nonhabitat_cells_left)
      Naffected  <- rpois(1, lambdaNew*Ncellsleft)
      Naffected <- min(Ncellsleft, Naffected)
      affected_nonhabitat_cells <- sample(suitable_nonhabitat_cells_left, Naffected, replace=FALSE)
      if (length(affected_nonhabitat_cells)>0){
      fires <-  spread(mirror_hab, affected_nonhabitat_cells, 
                       pNew*mirror_hab, 0, NULL, 1e8, 4, iterations = dNew, id = TRUE, returnIndices = TRUE)
      affected_nonhabitat_cells <- union(affected_nonhabitat_cells, fires$indices)
      inds_affected <- which((habitat_df$loc_id %in% affected_nonhabitat_cells)&(habitat_df$year >=i)) 
      habitat_df$habitat[inds_affected] <- 1
      mirror_hab[inds_affected] <- 0
      hab[inds_affected] <- 1
      suitable_nonhabitat_cells_left <- setdiff(suitable_nonhabitat_cells_left, affected_nonhabitat_cells)
      }
    }
    if(rr <0.5){
      suitable_nonhabitat_cells <- habitat_df$loc_id[which((habitat_df$habitat==0) & (habitat_df$year==i))]
      suitable_nonhabitat_cells_left <- suitable_nonhabitat_cells
      Ncellsleft <- length(suitable_nonhabitat_cells_left)
      Naffected  <- rpois(1, lambdaNew*Ncellsleft)
      Naffected <- min(Ncellsleft, Naffected)
      affected_nonhabitat_cells <- sample(suitable_nonhabitat_cells_left, Naffected, replace=FALSE)
      if (length(affected_nonhabitat_cells)>0){
      fires <-  spread(mirror_hab, affected_nonhabitat_cells, 
                       pNew*mirror_hab, 0, NULL, 1e8, 4, iterations = dNew, id = TRUE, returnIndices = TRUE)
      affected_nonhabitat_cells <- union(affected_nonhabitat_cells, fires$indices)
      inds_affected <- which((habitat_df$loc_id %in% affected_nonhabitat_cells)&(habitat_df$year >=i)) 
      habitat_df$habitat[inds_affected] <- 1
      mirror_hab[inds_affected] <- 0
      hab[inds_affected] <- 1
      inds_affected <- which((habitat_df$loc_id %in% affected_nonhabitat_cells)&(habitat_df$year >=i)) 
      habitat_df$habitat[inds_affected] <- 1
      suitable_nonhabitat_cells_left <- setdiff(suitable_nonhabitat_cells_left, affected_nonhabitat_cells)
      }
      
      Ncellsleft <- length(suitable_habitat_cells_left)
      Naffected  <- rpois(1, lambda*Ncellsleft)
      Naffected <- min(Ncellsleft, Naffected)
      affected_habitat_cells <- sample(suitable_habitat_cells_left, Naffected, replace=FALSE)
      if (length(affected_habitat_cells)>0){
      fires <-  spread(hab, affected_habitat_cells, 
                       p*hab, 0, NULL, 1e8, 4, iterations = d, id = TRUE, returnIndices = TRUE)

      affected_habitat_cells <- union(affected_habitat_cells, fires$indices)
      inds_affected <- which((habitat_df$loc_id %in% affected_habitat_cells)&(habitat_df$year >=i)) 
      habitat_df$habitat[inds_affected] <- 0
      hab[inds_affected] <- 0
      mirror_hab[inds_affected] <- 1
      suitable_habitat_cells_left <- setdiff(suitable_habitat_cells_left, affected_habitat_cells)
      }
    }
    }

  if (FALSE){
    
    ggplot(habitat_df, aes(x, y, col = as.factor(habitat)))+theme_classic()+geom_point()+facet_wrap(.~year, ncol = 5)+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+  theme(axis.title.y = element_blank(),
                                                  axis.text.y =element_blank(),
                                                  axis.ticks.y =element_blank()) +
      theme(strip.background = element_blank(), strip.text = element_blank())+ scale_color_manual(values=group.colors)+ guides(color=FALSE)
    
  }
  return(habitat_df)
}


generate_habitat_change_poissonpoint_process <- function(habitat_df, change_type, Nyears, change_direction, N_xy, AM, cluster_mean_size, lambda, cluster_mean_sizeNew, lambdaNew){
  
  if (FALSE){
    habitat_df <- generate_habitat_df(xy[,1], xy[,2], initial_habitat_df$habitat, Nyears)
  }
  
  # generates the habitat 10-year time-series, assuming no change in it.
  suitable_habitat_cells <- which((habitat_df$habitat==1) & (habitat_df$year==1))
  initial_habitat <- sum(habitat_df[suitable_habitat_cells, 'habitat'])
  
  # if (FALSE){
  #   ggplot(habitat_df, aes(x, y, col = as.factor(habitat)))+theme_classic()+geom_point()+facet_wrap(.~year, ncol = 5)+
  #     theme(axis.title.x=element_blank(),
  #           axis.text.x=element_blank(),
  #           axis.ticks.x=element_blank())+  theme(axis.title.y = element_blank(),
  #                                                 axis.text.y =element_blank(),
  #                                                 axis.ticks.y =element_blank()) +
  #     theme(strip.background = element_blank(), strip.text = element_blank())+ scale_color_manual(values=group.colors)+ guides(color=FALSE)
  # }
  
  # two ideas combined:
  # 1. Poisson process with marks, that determine the size of the destruction.
  # 2. Poisson process that describes the intensity of the habitat change.  
  
  # cluster size := 
  # max_cluster_radius := 
  
  if (FALSE){
    lambda <- 0.01
    cluster_size <- 15
  }
  
  max_cluster_radius <- cluster_mean_size
  max_cluster_radiusNew <- cluster_mean_sizeNew

  # first habitat disappearance:
  suitable_habitat_cells_left <- suitable_habitat_cells 
  all_affected_habitat_cells <- NULL
  all_affected_nonhabitat_cells <- NULL

  for (i in seq(2, Nyears)){
    rr <- runif(1, min = 0, max = 1)
    if (rr>0.5){
    Ncellsleft <- length(suitable_habitat_cells_left)
    Naffected  <- rpois(1, lambda*Ncellsleft)
    Naffected <- min(Ncellsleft, Naffected)
    affected_habitat_cells <- sample(suitable_habitat_cells_left, Naffected, replace=FALSE)
    for (c in affected_habitat_cells){
      neighbours_to_cell <- which(AM[c, ]==1) # go through the affected habitat cells 
      # to construct a neighbourhood of size = cluster_size
      total_neighbourhood <- union(neighbours_to_cell, c)
      habitatcellswithinneighbours <- intersect(suitable_habitat_cells_left, total_neighbourhood) 
      Nhabitatcellswithinneighbours <- length(habitatcellswithinneighbours)
      cluster_size <- rpois(1, cluster_mean_size)
      Nattempts <- 1
      while (Nattempts <= max_cluster_radius){
        if (Nhabitatcellswithinneighbours<=cluster_size){
          if(Nhabitatcellswithinneighbours>1){
            neighbours_to_cell <- which(AM[total_neighbourhood, ]==1, arr.ind = T)[,2]}else{
              neighbours_to_cell <- which(AM[total_neighbourhood, ]==1, arr.ind = T)[2]} 
          # neighbours to all the cells so far selected to the neighbourhood. 
          total_neighbourhood <- union(neighbours_to_cell, total_neighbourhood)
          habitatcellswithinneighbours <- intersect(suitable_habitat_cells_left, total_neighbourhood) 
          Nhabitatcellswithinneighbours <- length(habitatcellswithinneighbours)
          Nattempts <- Nattempts+1
        }else{break}
      }
      Nselectedcells<- min(cluster_size, Nhabitatcellswithinneighbours) 
      affected_habitat_cells <- sample(habitatcellswithinneighbours, Nselectedcells) 
      inds_affected <- which((habitat_df$loc_id %in% affected_habitat_cells)&(habitat_df$year >=i)) 
      habitat_df$habitat[inds_affected] <- 0
      suitable_habitat_cells_left <- setdiff(suitable_habitat_cells_left, affected_habitat_cells)
    }
        
    suitable_nonhabitat_cells <- habitat_df$loc_id[which((habitat_df$habitat==0) & (habitat_df$year==i))]
    suitable_nonhabitat_cells_left <- suitable_nonhabitat_cells
    Ncellsleft <- length(suitable_nonhabitat_cells_left)
    Naffected  <- rpois(1, lambdaNew*Ncellsleft)
    Naffected <- min(Ncellsleft, Naffected)
    affected_nonhabitat_cells <- sample(suitable_nonhabitat_cells_left, Naffected, replace=FALSE)
    for (c in affected_nonhabitat_cells){
      neighbours_to_cell <- which(AM[c, ]==1) # go through the affected habitat cells 
      # to construct a neighbourhood of size = cluster_size
      total_neighbourhood <- union(neighbours_to_cell, c)
      nonhabitatcellswithinneighbours <- intersect(suitable_nonhabitat_cells_left, total_neighbourhood) 
      Nnonhabitatcellswithinneighbours <- length(nonhabitatcellswithinneighbours)
      cluster_size <- rpois(1, cluster_mean_sizeNew)
      Nattempts <- 1
      while (Nattempts <= max_cluster_radiusNew){
        if (Nnonhabitatcellswithinneighbours<=cluster_size){
          if(Nnonhabitatcellswithinneighbours>1){
            neighbours_to_cell <- which(AM[total_neighbourhood, ]==1, arr.ind = T)[,2]}else{
              neighbours_to_cell <- which(AM[total_neighbourhood, ]==1, arr.ind = T)[2]} 
          # neighbours to all the cells so far selected to the neighbourhood. 
          total_neighbourhood <- union(neighbours_to_cell, total_neighbourhood)
          nonhabitatcellswithinneighbours <- intersect(suitable_nonhabitat_cells_left, total_neighbourhood) 
          Nnonhabitatcellswithinneighbours <- length(nonhabitatcellswithinneighbours)
          Nattempts <- Nattempts+1
        }else{break}
      }
      Nselectedcells<- min(cluster_size, Nnonhabitatcellswithinneighbours) 
      affected_nonhabitat_cells <- sample(nonhabitatcellswithinneighbours, Nselectedcells) 
      inds_affected <- which((habitat_df$loc_id %in% affected_nonhabitat_cells)&(habitat_df$year >=i)) 
      habitat_df$habitat[inds_affected] <- 1
      suitable_nonhabitat_cells_left <- setdiff(suitable_nonhabitat_cells_left, affected_nonhabitat_cells)
    }
    }
    # IN REVERSE ORDER:

    if(rr <0.5){
    suitable_nonhabitat_cells <- habitat_df$loc_id[which((habitat_df$habitat==0) & (habitat_df$year==i))]
    suitable_nonhabitat_cells_left <- suitable_nonhabitat_cells
    Ncellsleft <- length(suitable_nonhabitat_cells_left)
    Naffected  <- rpois(1, lambdaNew*Ncellsleft)
    Naffected <- min(Ncellsleft, Naffected)
    affected_nonhabitat_cells <- sample(suitable_nonhabitat_cells_left, Naffected, replace=FALSE)
    for (c in affected_nonhabitat_cells){
      neighbours_to_cell <- which(AM[c, ]==1) # go through the affected habitat cells 
      # to construct a neighbourhood of size = cluster_size
      total_neighbourhood <- union(neighbours_to_cell, c)
      nonhabitatcellswithinneighbours <- intersect(suitable_nonhabitat_cells_left, total_neighbourhood) 
      Nnonhabitatcellswithinneighbours <- length(nonhabitatcellswithinneighbours)
      cluster_size <- rpois(1, cluster_mean_sizeNew)
      Nattempts <- 1
      while (Nattempts <= max_cluster_radiusNew){
        if (Nnonhabitatcellswithinneighbours<=cluster_size){
          if(Nnonhabitatcellswithinneighbours>1){
            neighbours_to_cell <- which(AM[total_neighbourhood, ]==1, arr.ind = T)[,2]}else{
              neighbours_to_cell <- which(AM[total_neighbourhood, ]==1, arr.ind = T)[2]} 
          # neighbours to all the cells so far selected to the neighbourhood. 
          total_neighbourhood <- union(neighbours_to_cell, total_neighbourhood)
          nonhabitatcellswithinneighbours <- intersect(suitable_nonhabitat_cells_left, total_neighbourhood) 
          Nnonhabitatcellswithinneighbours <- length(nonhabitatcellswithinneighbours)
          Nattempts <- Nattempts+1
        }else{break}
      }
      Nselectedcells<- min(cluster_size, Nnonhabitatcellswithinneighbours) 
      affected_nonhabitat_cells <- sample(nonhabitatcellswithinneighbours, Nselectedcells) 
      inds_affected <- which((habitat_df$loc_id %in% affected_nonhabitat_cells)&(habitat_df$year >=i)) 
      habitat_df$habitat[inds_affected] <- 1
      suitable_nonhabitat_cells_left <- setdiff(suitable_nonhabitat_cells_left, affected_nonhabitat_cells)
    }
    
    Ncellsleft <- length(suitable_habitat_cells_left)
    Naffected  <- rpois(1, lambda*Ncellsleft)
    Naffected <- min(Ncellsleft, Naffected)
    affected_habitat_cells <- sample(suitable_habitat_cells_left, Naffected, replace=FALSE)
    for (c in affected_habitat_cells){
      neighbours_to_cell <- which(AM[c, ]==1) # go through the affected habitat cells 
      # to construct a neighbourhood of size = cluster_size
      total_neighbourhood <- union(neighbours_to_cell, c)
      habitatcellswithinneighbours <- intersect(suitable_habitat_cells_left, total_neighbourhood) 
      Nhabitatcellswithinneighbours <- length(habitatcellswithinneighbours)
      cluster_size <- rpois(1, cluster_mean_size)
      Nattempts <- 1
      while (Nattempts <= max_cluster_radius){
        if (Nhabitatcellswithinneighbours<=cluster_size){
          if(Nhabitatcellswithinneighbours>1){
            neighbours_to_cell <- which(AM[total_neighbourhood, ]==1, arr.ind = T)[,2]}else{
              neighbours_to_cell <- which(AM[total_neighbourhood, ]==1, arr.ind = T)[2]} 
          # neighbours to all the cells so far selected to the neighbourhood. 
          total_neighbourhood <- union(neighbours_to_cell, total_neighbourhood)
          habitatcellswithinneighbours <- intersect(suitable_habitat_cells_left, total_neighbourhood) 
          Nhabitatcellswithinneighbours <- length(habitatcellswithinneighbours)
          Nattempts <- Nattempts+1
        }else{break}
      }
      Nselectedcells<- min(cluster_size, Nhabitatcellswithinneighbours) 
      affected_habitat_cells <- sample(habitatcellswithinneighbours, Nselectedcells) 
      inds_affected <- which((habitat_df$loc_id %in% affected_habitat_cells)&(habitat_df$year >=i)) 
      habitat_df$habitat[inds_affected] <- 0
      suitable_habitat_cells_left <- setdiff(suitable_habitat_cells_left, affected_habitat_cells)
    }
  }
    
    
  } 
  
  if (FALSE){
    
    ggplot(habitat_df, aes(x, y, col = as.factor(habitat)))+theme_classic()+geom_point()+facet_wrap(.~year, ncol = 5)+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+  theme(axis.title.y = element_blank(),
                                                  axis.text.y =element_blank(),
                                                  axis.ticks.y =element_blank()) +
      theme(strip.background = element_blank(), strip.text = element_blank())+ scale_color_manual(values=group.colors)+ guides(color=FALSE)
    
  }
  return(habitat_df)
}



generate_habitat_change <- function(habitat_df, change_type, Nyears, change_direction, N_xy, AM){
  # generates the habitat 10-year time-series, assuming no change in it.
  suitable_habitat_cells <- which((habitat_df$habitat==1) & (habitat_df$year==1))
  initial_habitat <- sum(habitat_df[suitable_habitat_cells, 'habitat'])
  
  if (change_direction =='decay'){
    if (change_type == 'random'){
      suitable_habitat_cells_left <- suitable_habitat_cells 
     for (i in seq(2, Nyears)){
      Ncellsleft <- length(suitable_habitat_cells_left)
      Naffected <- sample(seq(1,floor(Ncellsleft/3)), 1)
      affected_habitat_cells <- sample(suitable_habitat_cells_left, Naffected, replace=FALSE)
      inds_affected <- which((habitat_df$loc_id %in% affected_habitat_cells)&(habitat_df$year >=i)) 
      habitat_df$habitat[inds_affected] <- 0
      suitable_habitat_cells_left <- setdiff(suitable_habitat_cells_left, affected_habitat_cells)
     } 
    }
    if (change_type=='localized'){
      suitable_habitat_cells_left <- suitable_habitat_cells 
      all_affected_habitat_cells <- NULL
      for (i in seq(2, Nyears)){
        Ncellsleft <- length(suitable_habitat_cells_left)
        Naffected <- sample(seq(1,floor(Ncellsleft/3)), 1)
        if (i==2){
        affected_habitat_cells <- sample(suitable_habitat_cells_left, Naffected, replace=FALSE)
        }
        all_affected_habitat_cells <- union(all_affected_habitat_cells, affected_habitat_cells)
        neighbours_to_affected_cells <- union(suitable_habitat_cells, setdiff(unique(AM[all_affected_habitat_cells, ]==1), all_affected_habitat_cells))
        Nneighbours <- length(neighbours_to_affected_cells)
        Naffectedneighbours <- sample(seq(1, floor(Nneighbours/3)))
        affected_habitat_cells <- sample(neighbours_to_affected_cells, Naffectedneighbours)
        inds_affected <- which((habitat_df$loc_id %in% affected_habitat_cells)&(habitat_df$year >=i)) 
        habitat_df$habitat[inds_affected] <- 0
        suitable_habitat_cells_left <- setdiff(suitable_habitat_cells_left, affected_habitat_cells)
      } 
    }
  }
  return(habitat_df)
}


generate_habitat_dynamics <- function(habitat_style, habitat_distribution,  l, habitat_change,  habitat_trend, n_x, n_y, Kt, DM, AM){
  
  N_xy <- n_x*n_y
  plotting <- FALSE
  
  sigma2 <- 1
  xy <- set_coordinates(n_x, n_y)
  
  DM <- Dist(xy, method = "euclidean", square = TRUE, p = 0,vector = FALSE)
  AM <- generate_adjacency_matrix(DM, treshold= 5)
  Kt <- Matern32Covariance(xy, xy, l, sigma2, DM)
  
  
  if (FALSE){
    habitat_style = 'binary'
    habitat_distribution = 'low_range'
  }
  
  if (habitat_style=='binary'){ #is habitat described by binary suitability, or continuously distributed values describing its quality?
   # binary habitat quality: 
  # if (habitat_distribution=='low_range'){
   #   l <- 0.5
   # }else{
   #   l <- 5
   # }
   #  sigma2 <- 1
   #  Kt <- Matern32Covariance(xy, xy, l, sigma2)
    habitat <- set_habitat(N_xy, Kt)
    habitat_treshold <- quantile(habitat, 0.5)
    habitat <- 1*(habitat>habitat_treshold)
    
  }else{
    # continuous habitat quality
    if (habitat_distribution=='low_range'){
      l <- 0.5
    }else{
      l <- 5
    }
    #sigma2 <- 1
    #Kt <- Matern32Covariance(xy, xy, l, sigma2)
    habitat <- set_habitat(n_x*n_y, Kt)
  }
  
  # if (FALSE){
  #   df <- data.frame('x' = xy[,1], 'y' = xy[,2], 'habitat' = habitat)
  #   ggplot(df, aes('x' = x, 'y' =y, 'col'  = habitat))+geom_point()
  # }
  # 
  initial_habitat <- habitat
  #next generate 10 year habitat time-series
  Nyears <- 10 
 
  if (habitat_trend =='decay'){
    if (habitat_change== 'localized'){
      # simulate localized decay in the habitat:
      if (habitat_style=='binary'){
        habitat_df <- generate_habitat_df(xy[,1], xy[,2], initial_habitat, Nyears)
        habitat_df <- generate_habitat_change(habitat_df, 'localized', 10, 'decay', n_x*n_y, AM)
        if (plotting==TRUE){ggplot(habitat_df, aes(x =x, y = y, col = habitat))+geom_point()+facet_wrap(.~year)}
      }else{ # (continuous habitat)
        
      }
    }else{
      # simulate randomly distributed decay in the habitat:
      habitat_df <- generate_habitat_df(xy[,1], xy[,2], initial_habitat, Nyears)
      habitat_df <- generate_habitat_change(habitat_df, 'random', 10, 'decay', n_x*n_y, AM)
      if (plotting==TRUE){ggplot(habitat_df, aes(x =x, y = y, col = habitat))+geom_point()+facet_wrap(.~year)}
    }
  }
  return(habitat_df)
}



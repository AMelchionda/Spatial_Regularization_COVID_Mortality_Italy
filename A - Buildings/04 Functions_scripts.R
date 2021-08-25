###########################################################
######################## FUNCTIONS ########################
###########################################################


# main_path <- 'C:/Users/aless/Documents/NAPDE project'
# main_path <- '/Users/Cheesecake/Desktop/NAPDE/Project'
# main_path <- 'C:/Users/Tommaso/Desktop/NAPDE/Project'
main_path <- ...

setwd(paste(main_path, '/Resources', sep = ''))

# COMPUTE_Z_MATRIX --------------------------------------------------------

compute_z_matrix <- function(mainland, island, nx=150, ny=150){
  # EXPORT MODEL RESULTS FOR PYTHON PLOTTING
  # Called by 'plot_results'
  #
  # INPUTS:
  # - mainland = output of smooth.FEM for the mainland
  # - island   = output of smooth.FEM for Sardegna
  # - nx, ny   = refinement in x and y dimension
  #
  # OUTPUT:
  # list used by 'create_single_plot' or 'plot_results' function
  # to show the results
  
  library(fdaPDE)
  
  # Create the refinement in x and y direction
  ymin <- 36.7
  ymax <- 47.1
  xmin <- 6.5
  xmax <- 18.5
  x_coord <- ((1:nx)/nx)*(xmax - xmin) + xmin
  y_coord <- ((1:ny)/ny)*(ymax - ymin) + ymin
  
  # Evaluate the solution in the meshgrid
  Z = matrix(nrow = nx, ncol = ny)
  for(i in 1:nx){
    x = x_coord[i]
    coords = cbind(x, y_coord)
    Z[i,] <- eval.FEM(mainland$fit.FEM, coords)
  }
  
  ids_island_x <- which(x_coord >= 8 & x_coord <= 10)
  ids_island_y <- which(y_coord >= 39 & y_coord <= 42)
  for(i in ids_island_x){
    x = x_coord[i]
    coords = cbind(x, y_coord[ids_island_y])
    Z[i, ids_island_y] <- eval.FEM(island$fit.FEM, coords)
  }
  
  return(list(x = x_coord, y = y_coord, z = Z))
}


# COMPUTE_Z_MATRIX_ALT --------------------------------------------------------

compute_z_matrix_alt <- function(solution, locations, nx=150, ny=150){
  # EXPORT MODEL RESULTS FOR PYTHON PLOTTING
  # Called by 'plot_results'
  #
  # INPUTS:
  # - mainland = output of smooth.FEM for the mainland
  # - island   = output of smooth.FEM for Sardegna
  # - nx, ny   = refinement in x and y dimension
  #
  # OUTPUT:
  # list used by 'create_single_plot' or 'plot_results' function
  # to show the results
  
  library(fdaPDE)
  
  # Create the refinement in x and y direction
  ymin <- min(locations[,2])
  ymax <- max(locations[,2])
  xmin <- min(locations[,1])
  xmax <- max(locations[,1])
  x_coord <- ((1:nx)/nx)*(xmax - xmin) + xmin
  y_coord <- ((1:ny)/ny)*(ymax - ymin) + ymin
  
  # Evaluate the solution in the meshgrid
  Z = matrix(nrow = nx, ncol = ny)
  for(i in 1:nx){
    x = x_coord[i]
    coords = cbind(x, y_coord)
    Z[i,] <- eval.FEM(solution$fit.FEM, coords)
  }
  
  return(list(x = x_coord, y = y_coord, z = Z))
}


# CREATE_SINGLE_PLOT ------------------------------------------------------

create_single_plot <- function(solution, lambda, title, limits=NULL,
                               avoid_bar = FALSE, bar_height = 15,
                               col_values = NULL, difference = FALSE,
                               show_lambdas = TRUE, flag1 = FALSE,
                               flag2 = FALSE){
  
  # CREATE THE PLOT OF A PARTICULAR RESULT
  # called by the 'plot_results' function
  #
  # INUTS:
  # - solution   = output of the 'compute_z_matrix' function
  # - title      = title of the plot
  # - limits     = to make it comparable with other plots
  # - col_low    = to control the colorbar
  # - col_high   = to control the colorbar
  # - avoid_bar  = TRUE if you don't want the colorbar
  # - bar_height = controls the height of the colorbar
  #
  # OUTPUT:
  # ggplot object, you can save it to modify it or just call 
  # the function to plot the results
  
  library(ggplot2)
  library(scales)
  
  # Convert data in percentage and prepare them for plotting
  x <- solution[[1]]
  y <- solution[[2]]
  z <- solution[[3]]
  
  z <- as.vector(z)*100
  
  df <- expand.grid(x, y)
  
  if(is.null(limits)){
    zmax <- max(na.omit(z))
    zmin <- min(na.omit(z))
    plot_labels <- round(c(0, 0.25, 0.5, 0.75, 1)*(zmax - zmin) + zmin,
                         digits = 2)
  }else{
    zmax <- max(limits)
    zmin <- min(limits)
    plot_labels <- round(c(0, 0.25, 0.5, 0.75, 1)*(zmax - zmin) + zmin,
                         digits = 2)
    
  }
  
  lambda <- scientific(lambda, digits = 2)
  
  # Initialize the plot
  p1 <- ggplot(df, aes(x = Var1, y = Var2, z = z, fill = z)) +
    geom_raster() + #interpolate for success 
    theme(plot.caption = element_text(size = 10, vjust = 98.5),
          panel.background = element_rect(fill = "gray98"),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "gray", fill=NA, size=0.5),
          aspect.ratio = 1,
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) 
  
  if(flag1){
    p1 <- p1 + coord_flip()
  }
  
  # Add the title
  if(show_lambdas){
    p1 <- p1 + labs(title = title,
                    subtitle = bquote(lambda*''['Mainland']*'='*.(lambda[1])*
                                        '       '*
                                        lambda*''['Sardinia']*'='*.(lambda[2])))
  }else{
    
    p1 <- p1 + labs(title = title)
  }
  
  # Color it
  if(difference){
    if(!flag2){
      if(!is.null(col_values)){
        inf_pt <- (col_values[1]-min(na.omit(z)))/(max(na.omit(z))-min(na.omit(z)))
        sup_pt <- (col_values[2]-min(na.omit(z)))/(max(na.omit(z))-min(na.omit(z)))
      }else{
        sup_pt <- 1
        inf_pt <- 0
      }
      zero_pt <- (0.001-min(na.omit(z)))/(max(na.omit(z))-min(na.omit(z)))
      inf_int_pt <- (inf_pt + zero_pt)/2
      sup_int_pt <- (sup_pt + zero_pt)/2
    }else{
      inf_pt <- 0
      sup_pt <- 1
      zero_pt <- (0.001-col_values[1])/(col_values[2]-col_values[1])
      inf_int_pt <- (inf_pt + zero_pt)/2
      sup_int_pt <- (sup_pt + zero_pt)/2
    }
    
    colors <- hcl.colors(5, palette = 'Lisbon', rev = F)
    if(flag2){
      colors[3:5] <- hcl.colors(3, palette = 'Inferno', rev = F)
      colors[4] <- 'orangered3'
    }
      
    
    p1 <- p1 + scale_fill_gradientn(colors = colors,
                                    values = c(inf_pt, 
                                               inf_int_pt, 
                                               zero_pt, 
                                               sup_int_pt,
                                               sup_pt),
                                    na.value = 'transparent',
                                    breaks = plot_labels,
                                    labels = paste(plot_labels, '%'),
                                    limits = plot_labels[c(1, 5)])
  }else{
    colors <- hcl.colors(3, palette = 'Inferno', rev = F)
    colors[2] <- 'orangered3'
    p1 <- p1 + scale_fill_gradientn(colors = colors,
                                    values = c(0, 0.5, 1),
                                    na.value = 'transparent',
                                    breaks = plot_labels,
                                    labels = paste(plot_labels, '%'),
                                    limits = plot_labels[c(1, 5)])
  }
  
  if(avoid_bar){
    p1 <- p1 + guides(fill = FALSE)
  }else{
    p1 <- p1 + guides(fill = guide_colourbar(barwidth = 1, 
                                             barheight = bar_height,
                                             title = 'Mortality \nincidence\n'))
  }
  
  return(p1)
}



# PLOT_RESULTS ------------------------------------------------------------

plot_results <- function(solution_list, reverse = F, nx = 200, ny = 200, 
                         type=NULL, other_info = '', case = NULL,
                         show_lambda = TRUE, file_name = NULL){
  # PLOT RESULTS
  # 
  # INPUTS:
  # - solution_list = list of outputs of 'smooth.fem'. The list should be
  #                   of the form (mailand1, island1, mainland2, island2, ...)
  # - reverse       = set to TRUE if you want to reverse the colors
  # - type          = optional in general, mandatory if plotting a single
  #                    result: it should be a vector of the same length of
  #                    'plot_data_list', used in the title to say if the data
  #                    are pre covid, from 2020 or the difference
  # - other_info    = additional info to be displayed in the title
  #                    
  # OUTPUT:
  # ggplot of the results, with maximum three plots in the same window
  
  library(ggplot2)
  library(gridExtra)
  
  # Extract number of required plots
  if(length(solution_list) %% 2 != 0){
    stop("'solution_list' should be of te form
    (mainland1, island1, mainland2, island2, ...),
         and should contain an EVEN number of objects")
  }
  
  n_plots <- length(solution_list) / 2
  if(n_plots != 3){
    stop("Length of 'solution' list should be exactly 3")
  }
  
  # Compute z matrix
  print('Constructing matrices')
  plot_data_list <- list()
  for(i in 1:n_plots){
    print(paste(i, '/', n_plots, sep = ''))
    j <- 2*i
    sol_i <- compute_z_matrix(solution_list[[j-1]], solution_list[[j]], 
                              nx = nx, ny = ny)
    plot_data_list <- append(plot_data_list, list(sol_i))
  }
  print('Done')
  
  # Start plotting
  if(is.null(type)){
    type <- c('2011 - 2019', '2020', 'Difference')
  }
  if(!is.null(case)){
    other_info <- extract_title(paste(case, 'tot', sep = '_'))
  }
  x1 <- plot_data_list[[1]][[1]]
  y1 <- plot_data_list[[1]][[2]]
  z1 <- plot_data_list[[1]][[3]]
  x2 <- plot_data_list[[2]][[1]]
  y2 <- plot_data_list[[2]][[2]]
  z2 <- plot_data_list[[2]][[3]]
  x3 <- plot_data_list[[3]][[1]]
  y3 <- plot_data_list[[3]][[2]]
  z3 <- plot_data_list[[3]][[3]]
  limits <- c(min(na.omit(as.vector(z1)), 
                  na.omit(as.vector(z2))),
              max(na.omit(as.vector(z1)), 
                  na.omit(as.vector(z2))))
  limits <- limits*100
  col_values <- c(0, max(na.omit(as.vector(z1)), 
                         na.omit(as.vector(z2))))
  col_values <- col_values*100
  lambda_mainland <- solution_list[[1]]$optimization$lambda_solution
  lambda_sardegna <- solution_list[[2]]$optimization$lambda_solution
  lambda <- c(lambda_mainland, lambda_sardegna)
  plot1 <- create_single_plot(list(x1, y1, z1), 
                              paste(type[1], '  |  ', other_info), 
                              limits = limits, lambda = lambda, 
                              avoid_bar = F, bar_height = 10,
                              difference = F, col_values = col_values,
                              show_lambda = show_lambda)
  lambda_mainland <- solution_list[[3]]$optimization$lambda_solution
  lambda_sardegna <- solution_list[[4]]$optimization$lambda_solution
  lambda <- c(lambda_mainland, lambda_sardegna)
  plot2 <- create_single_plot(list(x2, y2, z2),  
                              paste(type[2], '  |  ', other_info), 
                              limits = limits, lambda = lambda, 
                              bar_height = 10, difference = F,
                              col_values = col_values,
                              show_lambda = show_lambda)
  lambda_mainland <- solution_list[[5]]$optimization$lambda_solution
  lambda_sardegna <- solution_list[[6]]$optimization$lambda_solution
  lambda <- c(lambda_mainland, lambda_sardegna)
  col_values <- c(-max(abs(na.omit(as.vector(z3)))),
                   max(abs(na.omit(as.vector(z3)))))
  col_values <- col_values*100
  plot3 <- create_single_plot(list(x3, y3, z3), 
                              paste(type[3], '  |  ', other_info), 
                              lambda = lambda, bar_height = 10,
                              difference = T, col_values = col_values,
                              show_lambda = show_lambda)
  plot1 <- plot1 + theme(plot.margin=unit(c(0,0,0,0.5),"cm"))
  plot2 <- plot2 + theme(plot.margin=unit(c(0,0.5,0,0),"cm"))
  plot3 <- plot3 + theme(plot.margin=unit(c(0,0,0,0),"cm"))
  
  plot_combined <- grid.arrange(plot1, plot2, plot3, nrow = 1, 
                                heights=unit(0.4, "npc"))
  
  if(!is.null(file_name)){
    ggsave(file_name, plot_combined, width=14, height=8.3) 
  }
  
  return(plot_combined)
    
  
}

# plot_results(solutions, reverse = T, other_info = 'First wave')


# CREATE Z MATRIX TIME ----------------------------------------------------

compute_z_matrix_time <- function(mainland, island, nx=150, ny=150, nt = NULL,
                                  case = NULL, tmax = NULL){
  # EXPORT MODEL RESULTS FOR PYTHON PLOTTING
  # Called by 'plot_results'
  #
  # INPUTS:
  # - mainland  = output of smooth.FEM for the mainland
  # - island    = output of smooth.FEM for Sardegna
  # - nx, ny    = refinement in x and y dimension
  # - nt        = refinement in time
  # - time_unit = '1st', 'biweek', 'week' or 'day
  #
  # OUTPUT:
  # list used by 'create_single_plot' or 'plot_results' function
  # to show the results
  
  library(fdaPDE)
  
  # Adjust time
  if(!is.null(case)){
    if(case == '1st'){
      tmax <- 15
    }
    if(case == '2nd'){
      tmax <- 13
    }
    if(case == 'year'){
      tmax <- 12
    }
    
    nt <- tmax
  }
  
  if(is.null(nt)){
    nt <- tmax
  }
  
  timestamps <- 1 + (0:(nt-1)) * (tmax-1)/(nt-1)
  timestamps <- round(timestamps, digits = 3)
  
  # Create the refinement in x and y direction
  ymin <- 36.7
  ymax <- 47.1
  xmin <- 6.5
  xmax <- 18.5
  x_coord <- ((1:nx)/nx)*(xmax - xmin) + xmin
  y_coord <- ((1:ny)/ny)*(ymax - ymin) + ymin
  
  coords <- cbind(x_coord[1], y_coord)
  for(i in 2:nx){
    coords <- rbind(coords, cbind(x_coord[i], y_coord))
  }
  idx_main <- inside_indexes(coords, mainland$fit.FEM.time$FEMbasis)
  idx_island <- inside_indexes(coords, island$fit.FEM.time$FEMbasis)
  coords_main <- coords[idx_main,]
  coords_island <- coords[idx_island,]
  
  # Evaluate the solution in the meshgrid
  Z = matrix(nrow = nx*ny, ncol = nt)
  
  for(t in 1:nt){
    instant = timestamps[t]
    print(paste('t =', instant, '/', tmax))
    Z[idx_main, t] <- eval.FEM.time(mainland$fit.FEM.time, 
                                    locations = coords_main,
                                    time.instants = instant, 
                                    lambdaS = mainland$bestlambda[1], 
                                    lambdaT = mainland$bestlambda[2])
    
    Z[idx_island,t] <- eval.FEM.time(island$fit.FEM.time, 
                                     locations = coords_island,
                                     time.instants = instant, 
                                     lambdaS = island$bestlambda[1], 
                                     lambdaT = island$bestlambda[2])
  }
  
  Z <- as.data.frame(Z)
  colnames(Z) <- timestamps
  return(list(x = x_coord, y = y_coord, z = Z))
}


# CREATE SINGLE ANIM ------------------------------------------------------

create_single_anim <- function(solution, title, limits=NULL,
                               col_low='yellow', col_high='orangered3',
                               avoid_bar = FALSE, bar_height = 15,
                               reverse = F, threecolors = F,
                               col_values = NULL, case = NULL){
  library(ggplot2)
  library(scales)
  library(gganimate)
  
  x <- solution[[1]]
  y <- solution[[2]]
  z <- solution[[3]]
  
  pts <- expand.grid(y, x)
  df <- as.data.frame(z[,1])
  colnames(df) <- 'z'
  df$x <- pts$Var2
  df$y <- pts$Var1
  df$tt <- as.numeric(colnames(z)[1])
  for(i in 2:dim(z)[2]){
    temp <- as.data.frame(z[,i])
    colnames(temp) <- 'z'
    temp$x <- pts$Var2
    temp$y <- pts$Var1
    temp$tt <- as.numeric(colnames(z)[i])
    # temp$tt <- factor(temp$tt, levels = extract_title(temp$tt, '1st'))
    df <- rbind(df, temp)
  }
  df$tt <- as.factor(df$tt)
  if(!is.null(case)){
    levels(df$tt) <- extract_title(case)
  }
  
  if(is.null(limits)){
    zmax <- max(na.omit(z))
    zmin <- min(na.omit(z))
    plot_labels <- round(c(0, 0.25, 0.5, 0.75, 1)*(zmax - zmin) + zmin,
                         digits = 2)
  }else{
    zmax <- max(limits)
    zmin <- min(limits)
    plot_labels <- round(c(0, 0.25, 0.5, 0.75, 1)*(zmax - zmin) + zmin,
                         digits = 2)
    col_values = c(-max(abs(c(zmin, zmax))), 
                   max(abs(c(zmin, zmax))))
    
  }
  
  # Reverse colors if asked
  if(reverse){
    temp <- col_low
    col_low <-col_high
    col_high <- temp
  }
  
  # Initialize the plot
  p1 <- ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_raster() + 
    transition_states(tt) +
    theme(plot.caption = element_text(size = 10, vjust = 98.5),
          panel.background = element_rect(fill = "gray98"),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "gray", fill=NA, size=0.5),
          aspect.ratio = 1,
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  # Add the title
  p1 <- p1 + labs(title = paste(title, '|', '{closest_state}'))
  
  # Color it
  if(threecolors){
    if(!is.null(col_values)){
      inf_pt <- (col_values[1]-zmin)/(zmax-zmin)
      sup_pt <- (col_values[2]-zmin)/(zmax-zmin)
    }else{
      sup_pt <- 1
      inf_pt <- 0
    }
    zero_pt <- (0.001-zmin)/(zmax-zmin)
    
    p1 <- p1 + scale_fill_gradientn(colors = c('darkblue',
                                               'orangered3',
                                               'yellow1'),
                                    values = c(inf_pt, zero_pt, sup_pt),
                                    na.value = 'transparent',
                                    breaks = plot_labels,
                                    labels = paste(plot_labels, '%'),
                                    limits = plot_labels[c(1, 5)])
  }else{
    p1 <- p1 + scale_fill_gradient(low = col_low, high = col_high, 
                                   na.value = 'transparent',
                                   breaks = plot_labels,
                                   labels = paste(plot_labels, '%'),
                                   limits = plot_labels[c(1, 5)])
  }
  if(avoid_bar){
    p1 <- p1 + guides(fill = FALSE)
  }else{
    p1 <- p1 + guides(fill = guide_colourbar(barwidth = 1,
                                             barheight = bar_height,
                                             title = 'Mortality \nincidence\n'))
  }
  
  return(p1)
}


# EXTRACT_TITLE -----------------------------------------------------------



extract_title <- function(case){
  
  if(case == '1st_tot'){
    return('Jan 29 - May 12')
  }
  
  if(case == '2nd_tot'){
    return('Sep 30 - Dec 29')
  }
  
  if(case == 'year_tot'){
    return('Jan 01 - Dec 30')
  }
  
  if(case == '1st'){
    return(c('Jan 29 - Feb  4', 
             'Feb  5 - Feb 11',
             'Feb 12 - Feb 18',
             'Feb 19 - Feb 25',
             'Feb 26 - Mar  3',
             'Mar  4 - Mar 10',
             'Mar 11 - Mar 17',
             'Mar 18 - Mar 24',
             'Mar 25 - Mar 31',
             'Apr  1 - Apr  7',
             'Apr  8 - Apr 14',
             'Apr 15 - Apr 21',
             'Apr 22 - Apr 28',
             'Apr 29 - May  5',
             'May  6 - May 12'))
  }
  
  if(case == '2nd'){
    return(c('Sep 30 - Oct  6', 
             'Oct  7 - Oct 13',
             'Oct 14 - Oct 20',
             'Oct 21 - Oct 27',
             'Oct 28 - Nov  3',
             'Nov  4 - Nov 10',
             'Nov 11 - Nov 17',
             'Nov 18 - Nov 24',
             'Nov 25 - Dec  1',
             'Dec  2 - Dec  8',
             'Dec  9 - Dec 15',
             'Dec 16 - Dec 22',
             'Dec 23 - Dec 29'))
  }
  
  if(case == 'year'){
    return(c('January',
             'February',
             'March',
             'April',
             'May',
             'June',
             'July',
             'August',
             'September',
             'October',
             'November',
             'December'))
  }
  
}


# CREATE ANIMATIONS -------------------------------------------------------

create_animations <- function(solution_list, reverse = F, 
                              nx = 150, ny = 150, nt = 100,
                              case = NULL, tmax = NULL,
                              type=NULL, other_info = ''){
  # PLOT RESULTS TIME
  # 
  # INPUTS:
  # - solution_list = list of outputs of 'smooth.fem'. The list should be
  #                   of the form (mailand1, island1, mainland2, island2, ...)
  # - reverse       = set to TRUE if you want to reverse the colors
  # - type          = optional in general, mandatory if plotting a single
  #                    result: it should be a vector of the same length of
  #                    'plot_data_list', used in the title to say if the data
  #                    are pre covid, from 2020 or the difference
  # - other_info    = additional info to be displayed in the title
  #                    
  # OUTPUT:
  # ggplot of the results, with maximum three plots in the same window
  
  library(ggplot2)
  library(gganimate)
  
  # Extract number of required plots
  if(length(solution_list) %% 2 != 0){
    stop("'solution_list' should be of te form
    (mainland1, island1, mainland2, island2, ...),
         and should contain an EVEN number of objects")
  }
  
  n_plots <- length(solution_list) / 2
  if(n_plots > 3){
    stop("Length of 'solution' list should be at most 3")
  }
  
  # Compute z matrix
  print('Constructing matrices')
  plot_data_list <- list()
  for(i in 1:n_plots){
    print(paste(i, '/', n_plots, sep = ''))
    j <- 2*i
    sol_i <- compute_z_matrix_time(solution_list[[j-1]], 
                                   solution_list[[j]], 
                                   case = case,
                                   tmax = tmax,
                                   nx = nx, ny = ny, nt = nt)
    plot_data_list <- append(plot_data_list, list(sol_i))
  }
  print('Done')
  
  if(n_plots == 3){
    if(is.null(type)){
      type <- c('2011-2019', '2020', 'Difference')
    }
    x1 <- plot_data_list[[1]][[1]]
    y1 <- plot_data_list[[1]][[2]]
    z1 <- plot_data_list[[1]][[3]]*100
    x2 <- plot_data_list[[2]][[1]]
    y2 <- plot_data_list[[2]][[2]]
    z2 <- plot_data_list[[2]][[3]]*100
    x3 <- plot_data_list[[3]][[1]]
    y3 <- plot_data_list[[3]][[2]]
    z3 <- plot_data_list[[3]][[3]]*100
    limits <- c(min(na.omit(as.vector(z1)),
                           na.omit(as.vector(z2))),
                max(na.omit(as.vector(z1)), 
                    na.omit(as.vector(z2))))
    col_values <- c(-max(abs(limits)), max(abs(limits)))
    col_values <- col_values
    plot1 <- create_single_anim(list(x1, y1, z1), 
                                # paste(type[1], 'incidence', '-', other_info), 
                                type[1],
                                limits = limits, bar_height = 10, 
                                reverse = reverse, threecolors = T,
                                col_values = col_values, case = case)
    plot2 <- create_single_anim(list(x2, y2, z2), 
                                # paste(type[2], 'incidence', '-', other_info), 
                                type[2],
                                limits = limits, bar_height = 10, 
                                reverse = reverse, threecolors = T,
                                col_values = col_values, case = case)
    col_values <- c(-max(abs(na.omit(as.vector(z3)))),
                    max(abs(na.omit(as.vector(z3)))))
    col_values <- col_values
    plot3 <- create_single_anim(list(x3, y3, z3), 
                                # paste(type[3], '', '-', other_info), 
                                type[3],
                                bar_height = 10,
                                reverse = reverse, threecolors = T,
                                col_values = col_values, case = case)
    
    return(list(plot1, plot2, plot3))
    
  }
  
}


# SHOW ANIMATIONS ---------------------------------------------------------

show_animations <- function(animations, file = '', n_weeks = 15){
  
  library(magick)
  duration <- n_weeks * (20/15)
  
  print('Animating 1st plot')
  animation1 <- animate(animations[[1]], 
                        width = 300, height = 300, 
                        duration = 20)
  print('Animating 2nd plot')
  animation2 <- animate(animations[[2]], 
                        width = 300, height = 300, 
                        duration = 20)
  print('Animating 3rd plot')
  animation3 <- animate(animations[[3]], 
                        width = 300, height = 300, 
                        duration = 20)
  
  print('Wrapping up')
  a_mgif <- image_read(animation1)
  b_mgif <- image_read(animation2)
  c_mgif <- image_read(animation3)
  
  new_gif <- image_append(c(a_mgif[1], b_mgif[1], c_mgif[1]))
  for(i in 2:length(a_mgif)){
    combined <- image_append(c(a_mgif[i], b_mgif[i], c_mgif[i]))
    new_gif <- c(new_gif, combined)
  }
  
  filename1 <- paste('results/average', file, sep = '_')
  filename2 <- paste('results/2020', file, sep = '_')
  filename3 <- paste('results/difference', file, sep = '_')
  filename4 <- paste('results/complete', file, sep = '_')
  
  save_animation(animation1, file = paste(filename1, '.gif', sep=''))
  save_animation(animation2, file = paste(filename2, '.gif', sep=''))
  save_animation(animation3, file = paste(filename3, '.gif', sep=''))
  image_write_gif(new_gif, delay = duration/length(new_gif),
                  path = paste(filename4, '.gif', sep=''))
  
  return(new_gif)
  
}




# PLOT TIMESTAMPS ---------------------------------------------------------

plot_timestamps <- function(solution_list, nx = 200, ny = 200,
                            instants = NULL, type=NULL, other_info = '', 
                            case = NULL, show_lambda = TRUE, 
                            file_name = NULL){
  # PLOT RESULTS
  # 
  # INPUTS:
  # - solution_list = list of outputs of 'smooth.fem'. The list should be
  #                   of the form (mailand1, island1, mainland2, island2, ...)
  # - reverse       = set to TRUE if you want to reverse the colors
  # - type          = optional in general, mandatory if plotting a single
  #                    result: it should be a vector of the same length of
  #                    'plot_data_list', used in the title to say if the data
  #                    are pre covid, from 2020 or the difference
  # - other_info    = additional info to be displayed in the title
  #                    
  # OUTPUT:
  # ggplot of the results, with maximum three plots in the same window
  
  library(ggplot2)
  library(gridExtra)
  
  # Extract number of required plots
  if(length(solution_list) %% 2 != 0){
    stop("'solution_list' should be of te form
    (mainland1, island1, mainland2, island2, ...),
         and should contain an EVEN number of objects")
  }
  
  n_plots <- length(solution_list) / 2
  if(n_plots != 2){
    stop("Length of 'solution' list for time series should be exactly 2")
  }
  
  # Compute z matrix
  print('Constructing matrices')
  plot_data_list <- list()
  for(i in 1:n_plots){
    print(paste(i, '/', n_plots, sep = ''))
    j <- 2*i
    sol_i <- compute_z_matrix_time(solution_list[[j-1]], solution_list[[j]], 
                                   nx = nx, ny = ny, case = case)
    plot_data_list <- append(plot_data_list, list(sol_i))
  }
  print('Done')
  
  # Start plotting
  if(is.null(type)){
    type <- c('2011 - 2019', '2020')
  }
  if(!is.null(case)){
    other_info <- extract_title(case)
  }
  x1 <- plot_data_list[[1]][[1]]
  y1 <- plot_data_list[[1]][[2]]
  z1 <- plot_data_list[[1]][[3]]/ifelse(case == 'year', 100, 1)
  x2 <- plot_data_list[[2]][[1]]
  y2 <- plot_data_list[[2]][[2]]
  z2 <- plot_data_list[[2]][[3]]/ifelse(case == 'year', 100, 1)
  limits <- c(min(na.omit(as.vector(z1)), 
                  na.omit(as.vector(z2))),
              max(na.omit(as.vector(z1)), 
                  na.omit(as.vector(z2))))
  limits <- limits*100
  col_values <- c(min(na.omit(as.vector(z1)), 
                      na.omit(as.vector(z2))),
                  max(na.omit(as.vector(z1)), 
                      na.omit(as.vector(z2))))
  col_values <- col_values*100
  
  plots <- list()
  
  if(case == '1st'){
    tmax <- 15
  }
  if(case == '2nd'){
    tmax <- 13
  }
  if(case == 'year'){
    tmax <- 12
  }
  
  i = 1
  np = 1
  flag <- ifelse(case == 'year', T, F)
  while(i < 2*tmax){
    plot1 <- create_single_plot(list(x1, y1, z1[,np]), 
                                paste(type[1], '  |  ', other_info[np]), 
                                limits = limits, lambda = NULL, 
                                avoid_bar = F, bar_height = 9,
                                difference = flag, col_values = col_values,
                                show_lambda = F, flag1 = T,
                                flag2 = flag)
    plot2 <- create_single_plot(list(x2, y2, z2[,np]),  
                                paste(type[2], '  |  ', other_info[np]), 
                                limits = limits, lambda = NULL, 
                                bar_height = 9, difference = flag,
                                col_values = col_values,
                                show_lambda = F, flag1 = T,
                                flag2 = flag)
    plot1 <- plot1 + theme(plot.margin=unit(c(0,0,0,0.5),"cm"))
    plot2 <- plot2 + theme(plot.margin=unit(c(0,0.5,0,0),"cm"))
    
    plots[[i]] <- plot1
    i = i+1
    plots[[i]] <- plot2
    i = i+1
    
    np = np+1
    
  }
  
  return(plots)
  
}
# OUTSIDE_INDEXES ---------------------------------------------------------

outside_indexes <- function(locations, fem_basis){
  # OUTSIDE_INDEXES
  #
  # Returns the indexes of the elements outside the mesh
  #
  # INPUT:
  # - locations = matrix with the coordinates
  # - fem_basis = 'FEM.basis' object generated by the mesh
  #
  # OUTPUT:
  # vector containing the indexes to remove by setting e.g.
  # indexes <- outside_indexes(locations, fem_basis)

  coeff=rep(0,fem_basis$nbasis)
  fun=FEM(coeff,fem_basis)
  
  return(which(is.na(eval.FEM(fun, locations))))
}


# INSIDE_INDEXES ---------------------------------------------------------

inside_indexes <- function(locations, fem_basis){
  # INSIDEE_INDEXES
  #
  # Returns the indexes of the elements outside the mesh
  #
  # INPUT:
  # - locations = matrix with the coordinates
  # - fem_basis = 'FEM.basis' object generated by the mesh
  #
  # OUTPUT:
  # vector containing the indexes to remove by setting e.g.
  # indexes <- outside_indexes(locations, fem_basis)

  coeff=rep(0,fem_basis$nbasis)
  fun=FEM(coeff,fem_basis)
  
  return(which(!is.na(eval.FEM(fun, locations))))
}

#POSITION ----
position<-function(v,x,y){
  
  # INPUT:
  # - v = set of points (n x 2 matrix)
  # - x = x-coordinate of the point
  # - y = y-coordinate of the point
  
  # OUTPUT:
  # - p = index corresponding to the row of the matrix v which contains (x,y)
  #        (if it does not exist p is set to 0)
  
  n=dim(v)[1]
  w=1:n
  p=w[v==c(x,y)][1]
  if( is.na(p)){
    p=0
  }
  
  return(p)
}

#MERGE ------------------------------------
merge<-function(nodes,boundary){
  # INPUTS:
  # - nodes    = nodes (nn x 2 matrix)
  # - boundary = boundary nodes (nb x 2 matrix)
  
  # OUTPUTS: (list)
  # - output[[1]] = new set of nodes made merging nodes and boundary (avoiding duplicates)
  # - output[[2]] = list of indexes corresponding to the boundary nodes in output[[1]]
  # - output[[3]] = the number of boundary points in output[[1]] that do not belong to the nodes
  
  nn=dim(nodes)[1]
  nb=dim(boundary)[1]
  
  j=1;
  indexes=rep(0,nb)
  new_bd_indexes=rep(F,nb)
  
  for(i in 1:nb){
    p=position(nodes,boundary[i,1],boundary[i,2])
    if(p==0){
      indexes[i]=nn+j
      new_bd_indexes[i]=T;
      j=j+1
    }
    else{
      indexes[i]=p
    }
  }
  
  new_boundary=boundary[new_bd_indexes,]
  new_nodes=rbind(nodes,new_boundary)
  
  output=list(new_nodes,indexes,j-1)
  
  return(output)
}

#CREATE_PLOT_MESH --------------------------------

create_plot_mesh<-function(boundary_plot,boundary,nodes){
  
  # Create a mesh merging boundary_plot and  nodes (avoiding duplicates)
  #
  # INPUT:
  # - boundary_plot = the boundary of the new mesh
  # - boundary      = the former boundary
  # - nodes         = the nodes of the new mesh
  # 
  # OUTPUT: (list)
  # - out[[1]] = number of nodes of the new mesh
  # - out[[2]] = the number of points formerly belonging to nodes
  # - out[[3]] = the number of boundary points
  # - out[[4]] = segments structure
  # - out[[5]] = the mesh object
  # - out[[6]] = the fem basis object
  
  m=dim(boundary)[1]
  segments=cbind(1:m,c(2:m,1))
  
  mesh=create.mesh.2D(nodes = boundary,segments = segments)
  fem_basis=create.FEM.basis(mesh)
  indexes_out=outside_indexes(nodes,fem_basis)
  
  
  locations=nodes[-indexes_out,]
  
  out=merge(locations,boundary_plot)
  
  new_nodes=out[[1]]
  indexes_bd=out[[2]]
  
  n=dim(new_nodes)[1]
  
  m=dim(locations)[1]
  
  mb=length(indexes_bd)
  
  new_segments=cbind(indexes_bd[1:mb],indexes_bd[c(2:mb,1)])
  
  mesh=create.mesh.2D(nodes=new_nodes,segments = new_segments)
  
  fem_basis=create.FEM.basis(mesh)
  
  out=list(n,m,mb,new_segments,mesh,fem_basis)
  
  return(out)
}

## create_plot_mesh(island_boundary_plot,island_boundary_new, nodes)
## create_plot_mesh(mainland_boundary_new,mainland_boundary_new, nodes)



# CREATE_MATRIX_Z_COVARIATES ---------------------------

compute_z_matrix_covariates<-function(xmin,xmax,ymin,ymax,fun_main,fun_is,nx=100,ny=100){
  
  # INPUTS:
  # - locations = set of points
  # - fun_main  = function to plot (FEM object)
  # - fun_is    = function to plot (FEM object)
  # - nx, ny    = number of points along x and y directions
  
  # OUTPUT: 
  # - output[[1]] = vector of x coordinates
  # - output[[2]] = vector of y coordinates
  # - output[[3]] = matrix z that contains the evaluation of the function
  
  x=seq(xmin,xmax,(-xmin+xmax)/nx)
  y=seq(ymin,ymax,(-ymin+ymax)/ny)
  z=matrix(nrow=length(x),ncol = length(y))
  
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      z[i,j]=eval.FEM(fun_main,as.matrix(t(c(x[i],y[j]))))
      if(is.na(z[i,j])){
        z[i,j]=eval.FEM(fun_is,as.matrix(t(c(x[i],y[j]))))
      }
    }
  }
  
  output<-list(x,y,z)
}


# CREATE_SINGLE_PLOT_COVARIATES ------------------------------------------------------

create_single_plot_covariates <- function(xyz, lambda, beta, title, limits=NULL,
                                          avoid_bar = FALSE, bar_height = 15,
                                          col_values = NULL, difference = FALSE,
                                          show_lambdas = TRUE, show_betas = FALSE, flag1 = FALSE,
                                          flag2 = FALSE){
  
  # INPUTS:
  # - solution         = output of the compute_z_matrix_covariates function
  # - lambda           = vector of lambdas
  # - beta             = vector of betas
  # - title            = title of the plot
  # - limits           = to make it comparable with other plots
  # - col_values[[1]]  = to control the colorbar (col_low)
  # - col_values[[2]]  = to control the colorbar (col_high)
  # - avoid_bar        = TRUE if you don't want the colorbar
  # - bar_height       = controls the height of the colorbar
  # - show_lambdas     = show lambdas if TRUE
  # - show_beta        = show betas if TRUE
  #
  # OUTPUT:
  # ggplot object, you can save it to modify it or just call 
  # the function to plot the results
  
  library(ggplot2)
  library(scales)
  
  # Convert data in percentage and prepare them for plotting
  x <- xyz[[1]]
  y <- xyz[[2]]
  z <- xyz[[3]]
  
  z <- as.vector(z)*100
  
  df <- expand.grid(x, y)
  
  if(is.null(limits)){
    zmax <- max(na.omit(z))
    zmin <- min(na.omit(z))
    plot_labels <- round(c(0, 0.25, 0.5, 0.75, 1)*(zmax - zmin) + zmin,
                         digits = 2)
  }else{
    zmax <- max(limits)
    zmin <- min(limits)
    plot_labels <- round(c(0, 0.25, 0.5, 0.75, 1)*(zmax - zmin) + zmin,
                         digits = 2)
    
  }
  
  lambda <- scientific(lambda, digits = 2)
  beta <- scientific(beta, digits = 2)
  
  # Initialize the plot
  p1 <- ggplot(df, aes(x = Var1, y = Var2, z = z, fill = z)) +
    geom_raster() + #interpolate for success 
    theme(plot.caption = element_text(size = 10, vjust = 98.5),
          panel.background = element_rect(fill = "gray98"),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = "gray", fill=NA, size=0.5),
          aspect.ratio = 1,
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) 
  
  if(flag1){
    p1 <- p1 + coord_flip()
  }
  
  # Add the title
  if(show_lambdas && show_betas){
    p1 <- p1 + labs(title = title,
                    subtitle = bquote(lambda*''['Mainland']*'='*.(lambda[1])*
                                        '       '*
                                        beta*''['Mainland']*'='*.(beta[1])*
                                        '       '*
                                        lambda*''['Sardinia']*'='*.(lambda[2])*
                                        '       '*
                                        beta*''['Sardinia']*'='*.(beta[2])))
  }else{
    if(show_lambdas){
      p1 <- p1 + labs(title = title,
                      subtitle = bquote(lambda*''['Mainland']*'='*.(lambda[1])*
                                          '       '*
                                          lambda*''['Sardinia']*'='*.(lambda[2])))
    }else{
      if(show_betas){
        p1 <- p1 + labs(title = title, subtitle = bquote(beta*''['Mainland']*'='*.(beta[1])*
                                                           '       '*
                                                           beta*''['Sardinia']*'='*.(beta[2])))
      }else{
        p1 <- p1 + labs(title = title)
      }
    }
  }
  
  
  
  # Color it
  if(difference){
    if(!flag2){
      if(!is.null(col_values)){
        inf_pt <- (col_values[1]-min(na.omit(z)))/(max(na.omit(z))-min(na.omit(z)))
        sup_pt <- (col_values[2]-min(na.omit(z)))/(max(na.omit(z))-min(na.omit(z)))
      }else{
        sup_pt <- 1
        inf_pt <- 0
      }
      zero_pt <- (0.001-min(na.omit(z)))/(max(na.omit(z))-min(na.omit(z)))
      inf_int_pt <- (inf_pt + zero_pt)/2
      sup_int_pt <- (sup_pt + zero_pt)/2
    }else{
      inf_pt <- 0
      sup_pt <- 1
      zero_pt <- (0.001-col_values[1])/(col_values[2]-col_values[1])
      inf_int_pt <- (inf_pt + zero_pt)/2
      sup_int_pt <- (sup_pt + zero_pt)/2
    }
    
    colors <- hcl.colors(5, palette = 'Lisbon', rev = F)
    if(flag2){
      colors[3:5] <- hcl.colors(3, palette = 'Inferno', rev = F)
      colors[4] <- 'orangered3'
    }
    
    
    p1 <- p1 + scale_fill_gradientn(colors = colors,
                                    values = c(inf_pt, 
                                               inf_int_pt, 
                                               zero_pt, 
                                               sup_int_pt,
                                               sup_pt),
                                    na.value = 'transparent',
                                    breaks = plot_labels,
                                    labels = paste(plot_labels, '%'),
                                    limits = plot_labels[c(1, 5)])
  }else{
    colors <- hcl.colors(3, palette = 'Inferno', rev = F)
    colors[2] <- 'orangered3'
    p1 <- p1 + scale_fill_gradientn(colors = colors,
                                    values = c(0, 0.5, 1),
                                    na.value = 'transparent',
                                    breaks = plot_labels,
                                    labels = paste(plot_labels, '%'),
                                    limits = plot_labels[c(1, 5)])
  }
  
  if(avoid_bar){
    p1 <- p1 + guides(fill = FALSE)
  }else{
    p1 <- p1 + guides(fill = guide_colourbar(barwidth = 1, 
                                             barheight = bar_height,
                                             title = 'Mortality \nincidence\n'))
  }
  
  return(p1)
}



# PLOT_RESULTS_COVARIATES -------

plot_results_covariates <- function(solution, island_boundary_plot, island_boundary,
                                    mainland_boundary_plot, mainland_boundary, nodes,
                                    xmin=NA, xmax=NA, ymin=NA, ymax=NA, nx=150, ny=150, 
                                    title, limits=NULL,
                                    avoid_bar = FALSE, bar_height = 15,
                                    col_values = NULL, difference = FALSE,
                                    show_lambdas = TRUE, show_betas = FALSE, flag1 = FALSE,
                                    flag2 = FALSE){
  
  # INPUT:
  # - solution                       = pair of solutions (mainland, island)
  # - xmin, xmax, ymin, ymax         = coordinates defining the region to plot
  # - nx, ny                         = number of points along x and y directions
  # - boundary_plot(mainland/island) = the boundary of the new mesh
  # - boundary(mainland/island)      = the boundary of the model
  # - nodes                          = the locations of the new mesh
  # - title                          = title of the plot
  # - limits                         = to make it comparable with other plots
  # - col_low                        = to control the colorbar (col_low)
  # - col_high                       = to control the colorbar (col_high)
  # - avoid_bar                      = TRUE if you don't want the colorbar
  # - bar_height                     = controls the height of the colorbar
  # 
  # ggplot object, you can save it to modify it or just call 
  # the function to plot the results
  
  if(is.na(xmin)||is.na(ymin)||is.na(xmax)||is.na(ymax)){
    xmin=min(nodes[,1])
    xmax=max(nodes[,1])
    ymin=min(nodes[,2])
    ymax=max(nodes[,2])
  }
  
  #mainland-create mesh for plot
  out=create_plot_mesh(mainland_boundary_plot,mainland_boundary, nodes)
  
  n_main=out[[1]]
  plot_mesh_main=out[[5]]
  plot_femb_main=out[[6]]
  
  #mainland-defines the functions and the parameters
  solution_main=solution[[1]]
  beta_main=solution_main$solution$beta
  lambda_main=solution_main$optimization$lambda_solution
  z_hat_main=solution_main$solution$z_hat
  coeff_main=c(z_hat_main,rep(0,n_main-length(z_hat_main)))
  fun_main=FEM(coeff_main,plot_femb_main)
  
  #island-create mesh for plot
  out=create_plot_mesh(island_boundary_plot,island_boundary, nodes)
  
  n_is=out[[1]]
  plot_mesh_is=out[[5]]
  plot_femb_is=out[[6]]
  
  #island-defines the functions and the parameters
  solution_is=solution[[2]]
  beta_is=solution_is$solution$beta
  lambda_is=solution_is$optimization$lambda_solution
  z_hat_is=solution_is$solution$z_hat
  coeff_is=c(z_hat_is,rep(0,n_is-length(z_hat_is)))
  fun_is=FEM(coeff_is,plot_femb_is)
  
  #compute the matrix and set the parameters
  out<-compute_z_matrix_covariates(xmin,xmax,ymin,ymax,fun_main,fun_is,nx,ny)
  beta=c(beta_main,beta_is)
  lambda=c(lambda_main,lambda_is)
  
  #actual plotting
  p1 <- create_single_plot_covariates(out, lambda, beta, title, limits=limits,
                                      avoid_bar = avoid_bar, bar_height = bar_height,
                                      col_values = col_values, difference = difference,
                                      show_lambdas = show_lambdas, show_betas = show_betas, flag1 = flag1,
                                      flag2 = flag2)
  return(p1)
}

# plot_results_covariates(solution, island_boundary_plot, island_boundary_new,
#                         mainland_boundary_new, mainland_boundary_new, nodes,
#                         xmin=NA, xmax=NA, ymin=NA, ymax=NA, nx=150, ny=150, 
#                         title = 'insert title', limits=NULL,
#                         avoid_bar = FALSE, bar_height = 15,
#                         col_values = NULL, difference = FALSE,
#                         show_lambdas = TRUE, show_beta = FALSE, flag1 = FALSE,
#                         flag2 = FALSE)



# Save functions in .RData file -------------------------------------------

setwd('')
save(list = ls(), file = 'Functions.RData')


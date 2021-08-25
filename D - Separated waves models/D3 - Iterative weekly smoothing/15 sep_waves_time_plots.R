

# main_path <- 'C:/Users/aless/Documents/NAPDE project'
# main_path <- '/Users/Cheesecake/Desktop/NAPDE/Project'
# main_path <- 'C:/Users/Tommaso/Desktop/NAPDE/Project'
main_path <- ...

setwd(paste(main_path, '/Resources', sep = ''))
load('Functions.RData')
load('AllBoundaries.RData')
Data <- read.csv('weekly_dataset.csv', encoding = 'UTF-8')

setwd('../Results/First wave')
load('solutions_wave1.RData')
setwd('../Second wave')
load('solutions_wave2.RData')



#### Plots of the spatial distribution (1st wave) --------------------------
setwd('../First wave')

solutions <- list(mp1_solution, ip1_solution,
                  mc1_solution, ic1_solution)

# Produce GIF
produce_animation(solutions, case = '1st', file_name='1st_wave_animation')

# Produce timestamps
# timestamps <- plot_timestamps(solutions, case = '1st',
#                               file_name = '1st_wave_animation')
# plot_combined <- grid.arrange(timestamps[[1]], timestamps[[2]],
#                               timestamps[[5]], timestamps[[6]],
#                               timestamps[[9]], timestamps[[10]],
#                               timestamps[[13]], timestamps[[14]],
#                               timestamps[[17]], timestamps[[19]],
#                               timestamps[[21]], timestamps[[22]],
#                               timestamps[[25]], timestamps[[26]],
#                               timestamps[[29]], timestamps[[30]],
#                               ncol = 4)
# ggsave('Results/1st_wave_timestamp.pdf', plot_combined, width=16, height=13.5) 



#### Plots of the spatial distribution (2nd wave) --------------------------
setwd('../Second wave')

solutions <- list(mp2_solution, ip2_solution,
                  mc2_solution, ic2_solution)

# Produce GIF
produce_animation(solutions, case = '2nd', file_name='2nd_wave_animation')

# Produce timestamps
# timestamps <- plot_timestamps(solutions, case = '2nd',
#                               file_name = '2nd_wave_animation')
# plot_combined <- grid.arrange(timestamps[[1]], timestamps[[2]],
#                               timestamps[[5]], timestamps[[6]],
#                               timestamps[[9]], timestamps[[10]],
#                               timestamps[[13]], timestamps[[14]],
#                               timestamps[[17]], timestamps[[19]],
#                               timestamps[[21]], timestamps[[22]],
#                               timestamps[[25]], timestamps[[26]],
#                               timestamps[[29]], timestamps[[30]],
#                               ncol = 4)
# ggsave('Results/2nd_wave_timestamp.pdf', plot_combined, width=16, height=13.5) 


#### Mainland + Sicilia mesh -------------------------------------------

main_mesh <- create.mesh.2D(italy_boundary, nodesattributes = NA, 
                            segments = italy_segments, 
                            holes = NA,
                            triangles = NA, 
                            order = 1, 
                            verbosity = 0)
main_refinement <- refine.mesh.2D(main_mesh, minimum_angle = 25,
                                  maximum_area = 0.01, delaunay = T,
                                  verbosity = 0)

#### Sardegna mesh -------------------------------------------

island_mesh <- create.mesh.2D(sardegna_boundary, nodesattributes = NA, 
                              segments = sardegna_segments, 
                              holes = NA,
                              triangles = NA, 
                              order = 1, 
                              verbosity = 0)
island_refinement <- refine.mesh.2D(island_mesh, minimum_angle = 25,
                                    maximum_area = 0.01, delaunay = T,
                                    verbosity = 0)

#### Smoothing settings -------------------------------------------

# FEM basis
mainFEMbasis <- create.FEM.basis(main_refinement)
islandFEMbasis <- create.FEM.basis(island_refinement)
# Reduce datasets
mainland <- Data[-outside_indexes(cbind(Data$x, Data$y),
                                  mainFEMbasis), ]
island <- Data[-outside_indexes(cbind(Data$x, Data$y),
                                islandFEMbasis), ]
# Locations
mainloc <- mainland[c('x', 'y')]
islandloc <- island[c('x', 'y')]
# Only the observations change below



# Plot municipalities behaviour (1st wave) -------------------------------------------



mainland_pre <- mainland[, 12:26]
island_pre   <-   island[, 12:26]
mainland_cov <- mainland[, 64:78]
island_cov   <-   island[, 64:78]

ids <- NULL
mun <- c('Milano', 'Pesaro', 'Roma', 'Napoli', 'Cagliari',
         'Codogno', 'Mondolfo', 'Castelfidardo', 'Baronissi', 'Tortolì')
for(i in 1:10){
  if(i == 5 | i == 10){
    ids[i] <- which(island$name == mun[i])
  }else{
    ids[i] <- which(mainland$name == mun[i])
  }
}
# ids[5:8] <- sample(1:dim(mainloc)[1], size = 4)

time_instants <- 1+(0:100)*(12/100)
time_instants[length(time_instants)] <- 13
pred1 <- pred2 <- matrix(nrow = 10, ncol = 101)
for(i in 1:length(ids)){
  id <- ids[i]
  if(i == 5 | i == 10){
    pred1[i,] <- eval.FEM.time(ip1_solution$fit.FEM.time,
                               locations = islandloc[id,],
                               time.instants = time_instants, 
                               lambdaS = ip1_time_solution$bestlambda[1], 
                               lambdaT = ip1_time_solution$bestlambda[2])
    
    pred2[i,] <- eval.FEM.time(ic1_solution$fit.FEM.time,
                               locations = islandloc[id,],
                               time.instants = time_instants, 
                               lambdaS = ic1_time_solution$bestlambda[1], 
                               lambdaT = ic1_time_solution$bestlambda[2])
  }
  else{
    pred1[i,] <- eval.FEM.time(mp1_solution$fit.FEM.time,
                               locations = mainloc[id,],
                               time.instants = time_instants, 
                               lambdaS = mp1_time_solution$bestlambda[1], 
                               lambdaT = mp1_time_solution$bestlambda[2])
    
    pred2[i,] <- eval.FEM.time(mc1_solution$fit.FEM.time,
                               locations = mainloc[id,],
                               time.instants = time_instants, 
                               lambdaS = mc1_time_solution$bestlambda[1], 
                               lambdaT = mc1_time_solution$bestlambda[2])
  }
  
}

limits <- c(min(min(pred1, pred2), 
                min(mainland_pre[ids[c(1:4, 6:9)], ]),
                min(mainland_cov[ids[c(1:4, 6:9)], ]),
                min(island_pre[ids[c(5, 10)], ]),
                min(island_cov[ids[c(5, 10)], ])),
            max(max(pred1, pred2),
                max(mainland_pre[ids[c(1:4, 6:9)], ]),
                max(mainland_cov[ids[c(1:4, 6:9)], ]),
                max(island_pre[ids[c(5, 10)], ]),
                max(island_cov[ids[c(5, 10)], ]))
)*100
limits = c(0, 0.22)

par(mfrow = c(2,5))
for(i in 1:length(ids)){
  id <- ids[i]
  
  if(i == 5 | i == 10){
    plot(time_instants, pred1[i,]*100, type = 'l', col = 'orange',
         main = island$name[id], xaxt = 'n',
         xlab = '', ylab = '', ylim = limits, lwd = 2)
    points(1:15, island_pre[id, ]*100, pch = 16, col = 'orange')
    points(time_instants, pred2[i,]*100, type = 'l', col = 'orangered3',
           main = island$name[id], 
           xlab = '', ylab = '', ylim = limits, lwd = 2)
    points(1:15, island_cov[id, ]*100, pch = 16, col = 'orangered3')
  }else{
    plot(time_instants, pred1[i,]*100, type = 'l', col = 'orange',
         main = mainland$name[id], xaxt = 'n',
         xlab = '', ylab = '', ylim = limits, lwd = 2)
    points(1:15, mainland_pre[id, ]*100, pch = 16, col = 'orange')
    points(time_instants, pred2[i,]*100, type = 'l', col = 'orangered3',
           main = mainland$name[id], 
           xlab = '', ylab = '', ylim = limits, lwd = 2)
    points(1:15, mainland_cov[id, ]*100, pch = 16, col = 'orangered3')
  }
}
# 6.5 x 13.5




# Errors 1st wave ---------------------------------------------------------

epm <- ecm <- numeric(dim(mainloc)[1])
epi <- eci <- numeric(dim(islandloc)[1])
for(i in 1:dim(mainloc)[1]){
  print(paste(i, '/', dim(mainloc)[1]), sep = '')
  pred <- eval.FEM.time(mp2_solution$fit.FEM.time,
                        locations = mainloc[i,],
                        time.instants = 1:13,
                        lambdaS = mp2_solution$bestlambda[1], 
                        lambdaT = mp2_solution$bestlambda[2])
  epm[i] <- sum((pred - mainland_pre[i, ])^2)/13
  pred <- eval.FEM.time(mc2_solution$fit.FEM.time,
                        locations = mainloc[i,],
                        time.instants = 1:13,
                        lambdaS = mc2_solution$bestlambda[1], 
                        lambdaT = mc2_solution$bestlambda[2])
  ecm[i] <- sum((pred - mainland_cov[i, ])^2)/13
}
for(i in 1:dim(islandloc)[1]){
  print(paste(i, '/', dim(islandloc)[1]), sep = '')
  pred <- eval.FEM.time(ip2_solution$fit.FEM.time,
                        locations = islandloc[i,],
                        time.instants = 1:13,
                        lambdaS = ip2_solution$bestlambda[1], 
                        lambdaT = ip2_solution$bestlambda[2])
  epi[i] <- sum((pred - island_pre[i, ])^2)/13
  pred <- eval.FEM.time(ic2_solution$fit.FEM.time,
                        locations = islandloc[i,],
                        time.instants = 1:13,
                        lambdaS = ic2_solution$bestlambda[1], 
                        lambdaT = ic2_solution$bestlambda[2])
  eci[i] <- sum((pred - island_cov[i, ])^2)/13
}

scientific(sum(epm))
scientific(sum(epi))
scientific(sum(ecm))
scientific(sum(eci))
scientific(median(epm))
scientific(median(epi))
scientific(median(ecm))
scientific(median(eci))


# Plot municipalities behaviour (2nd wave) -------------------------------------------



mainland_pre <- mainland[, 45:57]
island_pre   <-   island[, 45:57]
mainland_cov <- mainland[, 97:109]
island_cov   <-   island[, 97:109]

ids <- NULL
mun <- c('Milano', 'Pesaro', 'Roma', 'Napoli', 'Cagliari',
         'Codogno', 'Mondolfo', 'Castelfidardo', 'Baronissi', 'Tortolì')
for(i in 1:10){
  if(i == 5 | i == 10){
    ids[i] <- which(island$name == mun[i])
  }else{
    ids[i] <- which(mainland$name == mun[i])
  }
}
# ids[5:8] <- sample(1:dim(mainloc)[1], size = 4)

time_instants <- 1+(0:100)*(12/100)
time_instants[length(time_instants)] <- 13
pred1 <- pred2 <- matrix(nrow = 10, ncol = 101)
for(i in 1:length(ids)){
  id <- ids[i]
  if(i == 5 | i == 10){
    pred1[i,] <- eval.FEM.time(ip2_solution$fit.FEM.time,
                               locations = islandloc[id,],
                               time.instants = time_instants, 
                               lambdaS = ip2_time_solution$bestlambda[1], 
                               lambdaT = ip2_time_solution$bestlambda[2])
    
    pred2[i,] <- eval.FEM.time(ic2_solution$fit.FEM.time,
                               locations = islandloc[id,],
                               time.instants = time_instants, 
                               lambdaS = ic2_time_solution$bestlambda[1], 
                               lambdaT = ic2_time_solution$bestlambda[2])
  }
  else{
    pred1[i,] <- eval.FEM.time(mp2_solution$fit.FEM.time,
                               locations = mainloc[id,],
                               time.instants = time_instants, 
                               lambdaS = mp2_time_solution$bestlambda[1], 
                               lambdaT = mp2_time_solution$bestlambda[2])
    
    pred2[i,] <- eval.FEM.time(mc2_solution$fit.FEM.time,
                               locations = mainloc[id,],
                               time.instants = time_instants, 
                               lambdaS = mc2_time_solution$bestlambda[1], 
                               lambdaT = mc2_time_solution$bestlambda[2])
  }
  
}

limits <- c(min(min(pred1, pred2), 
                min(mainland_pre[ids[c(1:4, 6:9)], ]),
                min(mainland_cov[ids[c(1:4, 6:9)], ]),
                min(island_pre[ids[c(5, 10)], ]),
                min(island_cov[ids[c(5, 10)], ])),
            max(max(pred1, pred2),
                max(mainland_pre[ids[c(1:4, 6:9)], ]),
                max(mainland_cov[ids[c(1:4, 6:9)], ]),
                max(island_pre[ids[c(5, 10)], ]),
                max(island_cov[ids[c(5, 10)], ]))
)*100
limits = c(0, 0.22)

par(mfrow = c(2,5))
for(i in 1:length(ids)){
  id <- ids[i]
  
  if(i == 5 | i == 10){
    plot(time_instants, pred1[i,]*100, type = 'l', col = 'orange',
         main = island$name[id], xaxt = 'n',
         xlab = '', ylab = '', ylim = limits, lwd = 2)
    points(1:13, island_pre[id, ]*100, pch = 16, col = 'orange')
    points(time_instants, pred2[i,]*100, type = 'l', col = 'orangered3',
           main = island$name[id], 
           xlab = '', ylab = '', ylim = limits, lwd = 2)
    points(1:13, island_cov[id, ]*100, pch = 16, col = 'orangered3')
  }else{
    plot(time_instants, pred1[i,]*100, type = 'l', col = 'orange',
         main = mainland$name[id], xaxt = 'n',
         xlab = '', ylab = '', ylim = limits, lwd = 2)
    points(1:13, mainland_pre[id, ]*100, pch = 16, col = 'orange')
    points(time_instants, pred2[i,]*100, type = 'l', col = 'orangered3',
           main = mainland$name[id], 
           xlab = '', ylab = '', ylim = limits, lwd = 2)
    points(1:13, mainland_cov[id, ]*100, pch = 16, col = 'orangered3')
  }
}
# 6.5 x 13.5


# Errors 2nd wave ---------------------------------------------------------

epm <- ecm <- numeric(dim(mainloc)[1])
epi <- eci <- numeric(dim(islandloc)[1])
for(i in 1:dim(mainloc)[1]){
  print(paste(i, '/', dim(mainloc)[1]), sep = '')
  pred <- eval.FEM.time(mp2_solution$fit.FEM.time,
                        locations = mainloc[i,],
                        time.instants = 1:13,
                        lambdaS = mp2_solution$bestlambda[1], 
                        lambdaT = mp2_solution$bestlambda[2])
  epm[i] <- sum((pred - mainland_pre[i, ])^2)/13
  pred <- eval.FEM.time(mc2_solution$fit.FEM.time,
                        locations = mainloc[i,],
                        time.instants = 1:13,
                        lambdaS = mc2_solution$bestlambda[1], 
                        lambdaT = mc2_solution$bestlambda[2])
  ecm[i] <- sum((pred - mainland_cov[i, ])^2)/13
}
for(i in 1:dim(islandloc)[1]){
  print(paste(i, '/', dim(islandloc)[1]), sep = '')
  pred <- eval.FEM.time(ip2_solution$fit.FEM.time,
                        locations = islandloc[i,],
                        time.instants = 1:13,
                        lambdaS = ip2_solution$bestlambda[1], 
                        lambdaT = ip2_solution$bestlambda[2])
  epi[i] <- sum((pred - island_pre[i, ])^2)/13
  pred <- eval.FEM.time(ic2_solution$fit.FEM.time,
                        locations = islandloc[i,],
                        time.instants = 1:13,
                        lambdaS = ic2_solution$bestlambda[1], 
                        lambdaT = ic2_solution$bestlambda[2])
  eci[i] <- sum((pred - island_cov[i, ])^2)/13
}

scientific(sum(epm))
scientific(sum(epi))
scientific(sum(ecm))
scientific(sum(eci))
scientific(median(epm))
scientific(median(epi))
scientific(median(ecm))
scientific(median(eci))


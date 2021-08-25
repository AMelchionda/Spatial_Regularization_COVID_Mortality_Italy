###########################################################
############## WEEKLY TIME SERIES SMOOTHING ###############
###########################################################

#### READ ME ####
# The idea will be to iterate in coarse grid, then fine grid
# For the coarse we'll only save vectors and indices
# This to not be saving entire solutions on R since we have a lot to run
# For the fine in time we'll CHECK convergence but only save vectors and indices
# For the fine in space we'll CHECK convergence

#### Load data + package ----------------------------------

# main_path <- 'C:/Users/aless/Documents/NAPDE project'
# main_path <- '/Users/Cheesecake/Desktop/NAPDE/Project'
# main_path <- 'C:/Users/Tommaso/Desktop/NAPDE/Project'
main_path <- ...

setwd(paste(main_path, '/Resources', sep = ''))
load('AllBoundaries.RData')
load('Functions.RData')
Data <- read.csv('weekly_dataset.csv')

setwd('../Results')

# Update package with GitHub: Download and install from local file
#install.packages("fdaPDE-master", repo=NULL, type="source")

library(fdaPDE)

#### Mainland + Sicilia mesh -------------------------------------------

main_mesh <- create.mesh.2D(italy_boundary, nodesattributes = NA, 
                            segments = italy_segments, 
                            holes = NA,
                            triangles = NA, 
                            order = 1, 
                            verbosity = 0)

plot(main_mesh)

main_refinement <- refine.mesh.2D(main_mesh, minimum_angle = 25,
                                  maximum_area = 0.01, delaunay = T,
                                  verbosity = 0)

plot(main_refinement)

#### Sardegna mesh -------------------------------------------

island_mesh <- create.mesh.2D(sardegna_boundary, nodesattributes = NA, 
                              segments = sardegna_segments, 
                              holes = NA,
                              triangles = NA, 
                              order = 1, 
                              verbosity = 0)

plot(island_mesh)

island_refinement <- refine.mesh.2D(island_mesh, minimum_angle = 25,
                                    maximum_area = 0.01, delaunay = T,
                                    verbosity = 0)

plot(island_refinement)

#### Smoothing settings + split mainland and Sardegna -------------------------------------------

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

# Observations
obs_mp1 <- mainland[,12:26]
obs_ip1 <- island[,12:26]
obs_mc1 <- mainland[,64:78]
obs_ic1 <- island[,64:78]
obs_mp2 <- mainland[,45:57]
obs_ip2 <- island[,45:57]
obs_mc2 <- mainland[,97:109]
obs_ic2 <- island[,97:109]

# Clear memory
rm(Data, mainland, island,
   island_mesh, island_refinement, main_mesh, main_refinement,
   italy_boundary, italy_segments, sardegna_boundary, sardegna_segments)

#### Iteration 1: Coarse grid in time ----

l_time <- seq(from = 10^-3, to = 10^2, length.out = 30)

# 1st wave ----

rm(obs_mp2, obs_ip2, obs_mc2, obs_ic2)

l_space_mp1 <- 3e-4
l_space_ip1 <- 3.1e-4
l_space_mc1 <- 7.9e-3
l_space_ic1 <- 8.7e-3

start <- Sys.time()
mp1_time_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mp1, mainFEMbasis, 
                                     time_locations = 1:15, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation = 'stochastic',
                                     lambdaS = l_space_mp1, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
end <- Sys.time()
print(paste("Smoothing complete with 30 lambda options", end - start))
mp1_time_solution <- mp1_time_solution[["bestlambda"]][2]

start <- Sys.time()
ip1_time_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ip1, islandFEMbasis, 
                                     time_locations = 1:15, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_ip1, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
end <- Sys.time()
ip1_time_solution <- ip1_time_solution[["bestlambda"]][2]

start <- Sys.time()
mc1_time_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mc1, mainFEMbasis, 
                                     time_locations = 1:15, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_mc1, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
end <- Sys.time()
mc1_time_solution <- mc1_time_solution[["bestlambda"]][2]

ic1_time_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ic1, islandFEMbasis, 
                                     time_locations = 1:15, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_ic1, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
ic1_time_solution <- ic1_time_solution[["bestlambda"]][2]

setwd('./First wave')
save(mp1_time_solution, mc1_time_solution,
     ip1_time_solution, ic1_time_solution,
     l_space_mp1, l_space_mc1, l_space_ip1, l_space_ic1, l_time,
     file = 'solutions_wave1_time.RData')
setwd('../')

# 2nd wave ----

rm(obs_mp1, obs_ip1, obs_mc1, obs_ic1)

l_space_mp2 <- 0.0003063419
l_space_ip2 <- 0.006756503
l_space_mc2 <- 0.008786872
l_space_ic2 <- 0.008545056

mp2_time_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mp2, mainFEMbasis, 
                                     time_locations = 1:13, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_mp2, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
mp2_time_solution <- mp2_time_solution[["bestlambda"]][2]

ip2_time_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ip2, islandFEMbasis, 
                                     time_locations = 1:13, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_ip2, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
ip2_time_solution <- ip2_time_solution[["bestlambda"]][2]

mc2_time_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mc2, mainFEMbasis, 
                                     time_locations = 1:13, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_mc2, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
mc2_time_solution <- mc2_time_solution[["bestlambda"]][2]

ic2_time_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ic2, islandFEMbasis, 
                                     time_locations = 1:13, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_ic2, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
ic2_time_solution <- ic2_time_solution[["bestlambda"]][2]

setwd('./Second wave')
save(mp2_time_solution, mc2_time_solution,
     ip2_time_solution, ic2_time_solution,
     l_space_mp2, l_space_mc2, l_space_ip2, l_space_ic2, l_time,
     file = 'solutions_wave2_time.RData')
setwd('../')

#### Iteration 2: Coarse grid in space ----

l_space <- seq(from = 10^-5, to = 1, length.out = 40)

# 1st wave ----

setwd('./First wave')
rm(obs_mp2, obs_ip2, obs_mc2, obs_ic2)
load('solutions_wave1_time.RData')
l_time_mp1 <- l_time[mp1_time_solution]
l_time_ip1 <- l_time[ip1_time_solution]
l_time_mc1 <- l_time[mc1_time_solution]
l_time_ic1 <- l_time[ic1_time_solution]

start <- Sys.time()
mp1_space_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mp1, mainFEMbasis, 
                                     time_locations = 1:15, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation = 'stochastic',
                                     lambdaS = l_space, lambdaT = l_time_mp1,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
end <- Sys.time()
mp1_space_solution <- mp1_space_solution[["bestlambda"]][1]

ip1_space_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ip1, islandFEMbasis, 
                                     time_locations = 1:15, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time_ip1,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
ip1_space_solution <- ip1_space_solution[["bestlambda"]][1]

mc1_space_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mc1, mainFEMbasis, 
                                     time_locations = 1:15, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time_mc1,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
mc1_space_solution <- mc1_space_solution[["bestlambda"]][1]

ic1_space_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ic1, islandFEMbasis, 
                                     time_locations = 1:15, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time_ic1,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
ic1_space_solution <- ic1_space_solution[["bestlambda"]][1]

save(mp1_space_solution, mc1_space_solution,
     ip1_space_solution, ic1_space_solution,
     l_time_mp1, l_time_mc1, l_time_ip1, l_time_ic1, l_space,
     file = 'solutions_wave1_space.RData')
setwd('../')

# 2nd wave ----

setwd('./Second wave')
rm(obs_mp1, obs_ip1, obs_mc1, obs_ic1)
load('solutions_wave2_time.RData')
l_time_mp2 <- l_time[mp2_time_solution]
l_time_ip2 <- l_time[ip2_time_solution]
l_time_mc2 <- l_time[mc2_time_solution]
l_time_ic2 <- l_time[ic2_time_solution]

mp2_space_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mp2, mainFEMbasis, 
                                     time_locations = 1:13, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time_mp2,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
mp2_space_solution <- mp2_space_solution[["bestlambda"]][1]

ip2_space_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ip2, islandFEMbasis, 
                                     time_locations = 1:13, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time_ip2,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
ip2_space_solution <- ip2_space_solution[["bestlambda"]][1]

mc2_space_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mc2, mainFEMbasis, 
                                     time_locations = 1:13, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time_mc2,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
mc2_space_solution <- mc2_space_solution[["bestlambda"]][1]

ic2_space_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ic2, islandFEMbasis, 
                                     time_locations = 1:13, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time_ic2,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
ic2_space_solution <- ic2_space_solution[["bestlambda"]][1]

save(mp2_space_solution, mc2_space_solution,
     ip2_space_solution, ic2_space_solution,
     l_time_mp2, l_time_mc2, l_time_ip2, l_time_ic2, l_space,
     file = 'solutions_wave2_space.RData')
setwd('../')

#### Iteration 3: Fine grid in time ----

l_time <- seq(from = 1, to = 500, length.out = 40)

# 1st wave ----

setwd('./First wave')
rm(obs_mp2, obs_ip2, obs_mc2, obs_ic2)
load('solutions_wave1_space.RData')
l_space_mp1 <- l_space[mp1_space_solution]
l_space_ip1 <- l_space[ip1_space_solution]
l_space_mc1 <- l_space[mc1_space_solution]
l_space_ic1 <- l_space[ic1_space_solution]

start <- Sys.time()
mp1_time_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mp1, mainFEMbasis, 
                                     time_locations = 1:15, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation = 'stochastic',
                                     lambdaS = l_space_mp1, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
end <- Sys.time()
mp1_time_solution <- mp1_time_solution[["bestlambda"]][2]

ip1_time_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ip1, islandFEMbasis, 
                                     time_locations = 1:15, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_ip1, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
ip1_time_solution <- ip1_time_solution[["bestlambda"]][2]

ic1_time_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ic1, islandFEMbasis, 
                                     time_locations = 1:15, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_ic1, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
ic1_time_solution <- ic1_time_solution[["bestlambda"]][2]

mc1_time_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mc1, mainFEMbasis, 
                                     time_locations = 1:15, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_mc1, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
mc1_time_solution[["GCV"]]
mc1_time_solution <- mc1_time_solution[["bestlambda"]][2]

save(mp1_time_solution, mc1_time_solution,
     ip1_time_solution, ic1_time_solution,
     l_space_mp1, l_space_mc1, l_space_ip1, l_space_ic1, l_time,
     file = 'solutions_wave1_time_fine.RData')
setwd('../')

# 2nd wave ----

setwd('./Second wave')
rm(obs_mp1, obs_ip1, obs_mc1, obs_ic1)
load('solutions_wave2_space.RData')
l_space_mp2 <- l_space[mp2_space_solution]
l_space_ip2 <- l_space[ip2_space_solution]
l_space_mc2 <- l_space[mc2_space_solution]
l_space_ic2 <- l_space[ic2_space_solution]

mp2_time_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mp2, mainFEMbasis, 
                                     time_locations = 1:13, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_mp2, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
mp2_time_solution <- mp2_time_solution[["bestlambda"]][2]

ip2_time_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ip2, islandFEMbasis, 
                                     time_locations = 1:13, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_ip2, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
ip2_time_solution <- ip2_time_solution[["bestlambda"]][2]

mc2_time_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mc2, mainFEMbasis, 
                                     time_locations = 1:13, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_mc2, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
mc2_time_solution <- mc2_time_solution[["bestlambda"]][2]

ic2_time_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ic2, islandFEMbasis, 
                                     time_locations = 1:13, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_ic2, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
ic2_time_solution <- ic2_time_solution[["bestlambda"]][2]

save(mp2_time_solution, mc2_time_solution,
     ip2_time_solution, ic2_time_solution,
     l_space_mp2, l_space_mc2, l_space_ip2, l_space_ic2, l_time,
     file = 'solutions_wave2_time_fine.RData')
setwd('../')

#### Iteration 4: Fine grid in space ----

l_space <- seq(from = 10^-2, to = 1, length.out = 40)
l_space_mp <- seq(from = 10^-7, to = 10^-5, length.out = 40)  # only different ones

# 1st wave ----

setwd('./First wave')
rm(obs_mp2, obs_ip2, obs_mc2, obs_ic2)
load('solutions_wave1_time_fine.RData')
l_time_mp1 <- l_time[mp1_time_solution]
l_time_ip1 <- l_time[ip1_time_solution]
l_time_mc1 <- l_time[mc1_time_solution]
l_time_ic1 <- l_time[ic1_time_solution]

mp1_space_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mp1, mainFEMbasis, 
                                      time_locations = 1:15, covariates = NULL,
                                      PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                      areal.data.avg = F, DOF.evaluation = 'stochastic',
                                      lambdaS = l_space_mp, lambdaT = l_time_mp1,
                                      lambda.selection.criterion='grid',
                                      lambda.selection.lossfunction = 'GCV')
mp1_space_solution <- mp1_space_solution[["bestlambda"]][1]

ip1_space_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ip1, islandFEMbasis, 
                                      time_locations = 1:15, covariates = NULL,
                                      PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                      areal.data.avg = F, DOF.evaluation ='stochastic',
                                      lambdaS = l_space, lambdaT = l_time_ip1,
                                      lambda.selection.criterion='grid',
                                      lambda.selection.lossfunction = 'GCV')
ip1_space_solution <- ip1_space_solution[["bestlambda"]][1]

start <- Sys.time()
mc1_space_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mc1, mainFEMbasis, 
                                      time_locations = 1:15, covariates = NULL,
                                      PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                      areal.data.avg = F, DOF.evaluation ='stochastic',
                                      lambdaS = l_space, lambdaT = l_time_mc1,
                                      lambda.selection.criterion='grid',
                                      lambda.selection.lossfunction = 'GCV')
end <- Sys.time()
mc1_space_solution <- mc1_space_solution[["bestlambda"]][1]

ic1_space_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ic1, islandFEMbasis, 
                                      time_locations = 1:15, covariates = NULL,
                                      PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                      areal.data.avg = F, DOF.evaluation ='stochastic',
                                      lambdaS = l_space, lambdaT = l_time_ic1,
                                      lambda.selection.criterion='grid',
                                      lambda.selection.lossfunction = 'GCV')
ic1_space_solution <- ic1_space_solution[["bestlambda"]][1]

save(mp1_space_solution, mc1_space_solution,
     ip1_space_solution, ic1_space_solution,
     l_time_mp1, l_time_mc1, l_time_ip1, l_time_ic1, l_space,
     file = 'solutions_wave1_space_fine.RData')
setwd('../')

# 2nd wave ----

setwd('./Second wave')
rm(obs_mp1, obs_ip1, obs_mc1, obs_ic1)
load('solutions_wave2_time_fine.RData')
l_time_mp2 <- l_time[mp2_time_solution]
l_time_ip2 <- l_time[ip2_time_solution]
l_time_mc2 <- l_time[mc2_time_solution]
l_time_ic2 <- l_time[ic2_time_solution]

mp2_space_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mp2, mainFEMbasis, 
                                      time_locations = 1:13, covariates = NULL,
                                      PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                      areal.data.avg = F, DOF.evaluation ='stochastic',
                                      lambdaS = l_space_mp, lambdaT = l_time_mp2,
                                      lambda.selection.criterion='grid',
                                      lambda.selection.lossfunction = 'GCV')
mp2_space_solution <- mp2_space_solution[["bestlambda"]][1]

ip2_space_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ip2, islandFEMbasis, 
                                      time_locations = 1:13, covariates = NULL,
                                      PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                      areal.data.avg = F, DOF.evaluation ='stochastic',
                                      lambdaS = l_space, lambdaT = l_time_ip2,
                                      lambda.selection.criterion='grid',
                                      lambda.selection.lossfunction = 'GCV')
ip2_space_solution <- ip2_space_solution[["bestlambda"]][1]

ic2_space_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ic2, islandFEMbasis, 
                                      time_locations = 1:13, covariates = NULL,
                                      PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                      areal.data.avg = F, DOF.evaluation ='stochastic',
                                      lambdaS = l_space, lambdaT = l_time_ic2,
                                      lambda.selection.criterion='grid',
                                      lambda.selection.lossfunction = 'GCV')
ic2_space_solution <- ic2_space_solution[["bestlambda"]][1]

mc2_space_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mc2, mainFEMbasis, 
                                      time_locations = 1:13, covariates = NULL,
                                      PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                      areal.data.avg = F, DOF.evaluation ='stochastic',
                                      lambdaS = l_space, lambdaT = l_time_mc2,
                                      lambda.selection.criterion='grid',
                                      lambda.selection.lossfunction = 'GCV')
mc2_space_solution <- mc2_space_solution[["bestlambda"]][1]

save(mp2_space_solution, mc2_space_solution,
     ip2_space_solution, ic2_space_solution,
     l_time_mp2, l_time_mc2, l_time_ip2, l_time_ic2, l_space,
     file = 'solutions_wave2_space_fine.RData')
setwd('../')

#### What we got so far ----

# 1st wave ----

setwd('./First wave')
load('solutions_wave1_time_fine.RData')
l_time_mp1 <- l_time[mp1_time_solution]
l_time_ip1 <- l_time[ip1_time_solution]
l_time_mc1 <- l_time[mc1_time_solution]
l_time_ic1 <- l_time[ic1_time_solution]
load('solutions_wave1_space_fine.RData')
l_space_mp1 <- l_space[mp1_space_solution]
l_space_ip1 <- l_space[ip1_space_solution]
l_space_mc1 <- l_space[mc1_space_solution]
l_space_ic1 <- l_space[ic1_space_solution]

mp1_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mp1, mainFEMbasis, 
                                      time_locations = 1:15, covariates = NULL,
                                      PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                      areal.data.avg = F, DOF.evaluation = 'stochastic',
                                      lambdaS = l_space_mp1, lambdaT = l_time_mp1,
                                      lambda.selection.criterion='grid',
                                      lambda.selection.lossfunction = 'GCV')

ip1_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ip1, islandFEMbasis, 
                                      time_locations = 1:15, covariates = NULL,
                                      PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                      areal.data.avg = F, DOF.evaluation ='stochastic',
                                      lambdaS = l_space_ip1, lambdaT = l_time_ip1,
                                      lambda.selection.criterion='grid',
                                      lambda.selection.lossfunction = 'GCV')

mc1_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mc1, mainFEMbasis, 
                                      time_locations = 1:15, covariates = NULL,
                                      PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                      areal.data.avg = F, DOF.evaluation ='stochastic',
                                      lambdaS = l_space_mc1, lambdaT = l_time_mc1,
                                      lambda.selection.criterion='grid',
                                      lambda.selection.lossfunction = 'GCV')

ic1_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ic1, islandFEMbasis, 
                                      time_locations = 1:15, covariates = NULL,
                                      PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                      areal.data.avg = F, DOF.evaluation ='stochastic',
                                      lambdaS = l_space_ic1, lambdaT = l_time_ic1,
                                      lambda.selection.criterion='grid',
                                      lambda.selection.lossfunction = 'GCV')

save(mp1_solution, mc1_solution, ip1_solution, ic1_solution,
     l_time_mp1, l_time_mc1, l_time_ip1, l_time_ic1,
     l_space_mp1, l_space_mc1, l_space_ip1, l_space_ic1,
     file = 'solutions_wave1.RData')
setwd('../')

# 2nd wave ----
setwd('./Second wave')

load('solutions_wave2_time_fine.RData')
l_time_mp2 <- l_time[mp2_time_solution]
l_time_ip2 <- l_time[ip2_time_solution]
l_time_mc2 <- l_time[mc2_time_solution]
l_time_ic2 <- l_time[ic2_time_solution]
load('solutions_wave2_space_fine.RData')
l_space_mp2 <- l_space[mp2_space_solution]
l_space_ip2 <- l_space[ip2_space_solution]
l_space_mc2 <- l_space[mc2_space_solution]
l_space_ic2 <- l_space[ic2_space_solution]

mp2_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mp2, mainFEMbasis, 
                                      time_locations = 1:13, covariates = NULL,
                                      PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                      areal.data.avg = F, DOF.evaluation ='stochastic',
                                      lambdaS = l_space_mp2, lambdaT = l_time_mp2,
                                      lambda.selection.criterion='grid',
                                      lambda.selection.lossfunction = 'GCV')

ip2_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ip2, islandFEMbasis, 
                                      time_locations = 1:13, covariates = NULL,
                                      PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                      areal.data.avg = F, DOF.evaluation ='stochastic',
                                      lambdaS = l_space_ip2, lambdaT = l_time_ip2,
                                      lambda.selection.criterion='grid',
                                      lambda.selection.lossfunction = 'GCV')

ic2_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ic2, islandFEMbasis, 
                                      time_locations = 1:13, covariates = NULL,
                                      PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                      areal.data.avg = F, DOF.evaluation ='stochastic',
                                      lambdaS = l_space_ic2, lambdaT = l_time_ic2,
                                      lambda.selection.criterion='grid',
                                      lambda.selection.lossfunction = 'GCV')

mc2_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mc2, mainFEMbasis, 
                                      time_locations = 1:13, covariates = NULL,
                                      PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                      areal.data.avg = F, DOF.evaluation ='stochastic',
                                      lambdaS = l_space_mc2, lambdaT = l_time_mc2,
                                      lambda.selection.criterion='grid',
                                      lambda.selection.lossfunction = 'GCV')

save(mp2_solution, mc2_solution, ip2_solution, ic2_solution,
     l_time_mp2, l_time_mc2, l_time_ip2, l_time_ic2,
     l_space_mp2, l_space_mc2, l_space_ip2, l_space_ic2,
     file = 'solutions_wave2.RData')
setwd('../')

#### Iteration 5: Another grid in time NOT DONE ----

l_time <- seq(from = 1, to = 500, length.out = 40)

# 1st wave ----
setwd('./First wave')
rm(obs_mp2, obs_ip2, obs_mc2, obs_ic2)
load('solutions_wave1_space_fine.RData')
l_space_mp1 <- l_space[mp1_space_solution]
l_space_ip1 <- l_space[ip1_space_solution]
l_space_mc1 <- l_space[mc1_space_solution]
l_space_ic1 <- l_space[ic1_space_solution]

start <- Sys.time()
mp1_time_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mp1, mainFEMbasis, 
                                     time_locations = 1:15, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation = 'stochastic',
                                     lambdaS = l_space_mp1, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
end <- Sys.time()
mp1_time_solution <- mp1_time_solution[["bestlambda"]][2]

ip1_time_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ip1, islandFEMbasis, 
                                     time_locations = 1:15, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_ip1, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
ip1_time_solution <- ip1_time_solution[["bestlambda"]][2]

ic1_time_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ic1, islandFEMbasis, 
                                     time_locations = 1:15, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_ic1, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
ic1_time_solution <- ic1_time_solution[["bestlambda"]][2]

mc1_time_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mc1, mainFEMbasis, 
                                     time_locations = 1:15, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_mc1, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
mc1_time_solution[["GCV"]]
mc1_time_solution <- mc1_time_solution[["bestlambda"]][2]

save(mp1_time_solution, mc1_time_solution,
     ip1_time_solution, ic1_time_solution,
     l_space_mp1, l_space_mc1, l_space_ip1, l_space_ic1, l_time,
     file = 'solutions_wave1_time_finer.RData')
setwd('../')

# 2nd wave ----

setwd('./Second wave')
rm(obs_mp1, obs_ip1, obs_mc1, obs_ic1)
load('solutions_wave2_space_fine.RData')
l_space_mp2 <- l_space[mp2_space_solution]
l_space_ip2 <- l_space[ip2_space_solution]
l_space_mc2 <- l_space[mc2_space_solution]
l_space_ic2 <- l_space[ic2_space_solution]

mp2_time_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mp2, mainFEMbasis, 
                                     time_locations = 1:13, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_mp2, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
mp2_time_solution <- mp2_time_solution[["bestlambda"]][2]

ip2_time_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ip2, islandFEMbasis, 
                                     time_locations = 1:13, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_ip2, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
ip2_time_solution <- ip2_time_solution[["bestlambda"]][2]

mc2_time_solution <- smooth.FEM.time(locations = mainloc, observations = obs_mc2, mainFEMbasis, 
                                     time_locations = 1:13, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_mc2, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
mc2_time_solution <- mc2_time_solution[["bestlambda"]][2]

ic2_time_solution <- smooth.FEM.time(locations = islandloc, observations = obs_ic2, islandFEMbasis, 
                                     time_locations = 1:13, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space_ic2, lambdaT = l_time,
                                     lambda.selection.criterion='grid',
                                     lambda.selection.lossfunction = 'GCV')
ic2_time_solution <- ic2_time_solution[["bestlambda"]][2]

save(mp2_time_solution, mc2_time_solution,
     ip2_time_solution, ic2_time_solution,
     l_space_mp2, l_space_mc2, l_space_ip2, l_space_ic2, l_time,
     file = 'solutions_wave2_time_finer.RData')
setwd('../')

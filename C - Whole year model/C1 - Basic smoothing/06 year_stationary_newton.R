###########################################################
##################### full year spatial ###################
###########################################################

#### Load data + package + set lambda ----------------------------------

# main_path <- 'C:/Users/aless/Documents/NAPDE project'
# main_path <- '/Users/Cheesecake/Desktop/NAPDE/Project'
# main_path <- 'C:/Users/Tommaso/Desktop/NAPDE/Project'
main_path <- ...

setwd(paste(main_path, '/Resources', sep = ''))
Data <- read.csv('static_dataset.csv')
load('AllBoundaries.RData')
load('Functions.RData')

setwd('../')

# Update package with GitHub
# install.packages("fdaPDE-master", repo=NULL, type="source")

library(fdaPDE)


# Lambda can be chosen by fdaPDE

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

lambda_search <- 'newton'

#### Smoothing pre covid -----------
# Notation:
# 1st letter is i/m island/mainland
# 2nd letter is p pre-covid


observations <- mainland$mean_norm_deaths
mp1_solution <- smooth.FEM(locations = mainloc, observations[], mainFEMbasis, 
                          covariates = NULL, PDE_parameters = NULL, BC = NULL,
                          incidence_matrix = NULL, areal.data.avg = F,
                          DOF.evaluation ='stochastic', 
                          lambda.selection.criterion = lambda_search, 
                          lambda.selection.lossfunction = 'GCV')

# image(mp_solution$fit.FEM)

observations <- island$mean_norm_deaths
ip1_solution <- smooth.FEM(locations = islandloc, observations[], islandFEMbasis, 
                          covariates = NULL, PDE_parameters = NULL, BC = NULL,
                          incidence_matrix = NULL, areal.data.avg = F,
                          DOF.evaluation ='stochastic', 
                          lambda.selection.criterion = lambda_search, 
                          lambda.selection.lossfunction = 'GCV')

# image(ip_solution$fit.FEM)

#### Smoothing with covid -----------
# Notation:
# 1st letter is i/m island/mainland
# 2nd letter is c covid

observations <- mainland$mean_norm_deaths_covid
mc1_solution <- smooth.FEM(locations = mainloc, observations[], mainFEMbasis, 
                          covariates = NULL, PDE_parameters = NULL, BC = NULL,
                          incidence_matrix = NULL, areal.data.avg = F,
                          DOF.evaluation ='stochastic', 
                          lambda.selection.criterion = lambda_search, 
                          lambda.selection.lossfunction = 'GCV')

# image(mc_solution$fit.FEM)

observations <- island$mean_norm_deaths_covid
ic1_solution <- smooth.FEM(locations = islandloc, observations[], islandFEMbasis, 
                          covariates = NULL, PDE_parameters = NULL, BC = NULL,
                          incidence_matrix = NULL, areal.data.avg = F,
                          DOF.evaluation ='stochastic', 
                          lambda.selection.criterion = lambda_search, 
                          lambda.selection.lossfunction = 'GCV')

# image(ic_solution$fit.FEM)

#### Smoothing covid increment -----------
# Notation:
# 1st letter is i/m island/mainland
# 2nd letter is d difference

observations <- mainland$mean_norm_deaths_covid - mainland$mean_norm_deaths
md1_solution <- smooth.FEM(locations = mainloc, observations[], mainFEMbasis,
                          covariates = NULL, PDE_parameters = NULL, BC = NULL,
                          incidence_matrix = NULL, areal.data.avg = F,
                          DOF.evaluation ='stochastic', 
                          lambda.selection.criterion = lambda_search, 
                          lambda.selection.lossfunction = 'GCV')

# image(md_solution$fit.FEM)

observations <- island$mean_norm_deaths_covid - island$mean_norm_deaths
id1_solution <- smooth.FEM(locations = islandloc, observations[], islandFEMbasis,
                          covariates = NULL, PDE_parameters = NULL, BC = NULL,
                          incidence_matrix = NULL, areal.data.avg = F,
                          DOF.evaluation ='stochastic', 
                          lambda.selection.criterion = lambda_search, 
                          lambda.selection.lossfunction = 'GCV')

# image(id_solution$fit.FEM)

#### Save results in a .RData file -------------------------------------------

setwd('./Results/Whole year')
save(mp1_solution, ip1_solution, 
     mc1_solution, ic1_solution, 
     md1_solution, id1_solution, 
     file = 'solutions_year_newton.RData')

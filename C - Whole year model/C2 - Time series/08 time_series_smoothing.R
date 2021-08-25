###########################################################
################## TIME SERIES SMOOTHING ##################
###########################################################

#### Load data + package + set lambda ----------------------------------

# main_path <- 'C:/Users/aless/Documents/NAPDE project'
# main_path <- '/Users/Cheesecake/Desktop/NAPDE/Project'
# main_path <- 'C:/Users/Tommaso/Desktop/NAPDE/Project'
main_path <- ...

setwd(paste(main_path, '/Resources', sep = ''))
load('AllBoundaries.RData')
load('Functions.RData')


# Update package with GitHub
#install.packages("fdaPDE-master", repo=NULL, type="source")

library(fdaPDE)

# Lambda can be chosen by fdaPDE
# For time smoothing, a vector of lambdas must be provided though
l_space <- 10^-5
l_time <- 10^-5
#l_space <- 10^-(2:5)
#l_time <- 10^-(2:5)

#### Choosing the data + split mainland and Sardegna ---------------

# Choose the data for its corresponding smoothing stage
freq <- c('monthly', 'biweekly', 'weekly')
# Change number below
choice <- paste(freq[1], 'dataset.csv', sep = '_')

Data <- read.csv(choice)
# Every time you change the data, run Smoothing settings again

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

#### Monthly smoothing -----------
# Notation:
# 1st letter is i/m island/mainland
# 2nd letter is p/c/d pre-covid/covid/difference
# 3rd letter is m/b/w monthly/biweekly/weekly

# First all pre covid
mpm_time_solution <- smooth.FEM.time(locations = mainloc, observations = 100 * mainland[,6:17], mainFEMbasis, 
                                     time_locations = 1:12, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')
# Patience
#plot(mpm_time_solution$fit.FEM, locations = mainloc[mainland$name=="Milano",])

# Errors: Obs - Smoothing
mpm_time_e <- 100 * mainland[,6:17] - eval.FEM.time(mpm_time_solution$fit.FEM.time, locations = mainloc, time.instants = 1:12)
# hist(log(as.matrix(mpm_time_e ^2)))

ipm_time_solution <- smooth.FEM.time(locations = islandloc, observations = 100 * island[,6:17], islandFEMbasis, 
                                     time_locations = 1:12, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')
#plot(ipm_time_solution$fit.FEM)
# Works

# Errors: Obs - Smoothing
ipm_time_e <- 100 * island[,6:17] - eval.FEM.time(ipm_time_solution$fit.FEM.time, locations = islandloc, time.instants = 1:12)


# Second all covid
mcm_time_solution <- smooth.FEM.time(locations = mainloc, observations = 100 * mainland[,18:29], mainFEMbasis, 
                                     time_locations = 1:12, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')
#plot(mcm_time_solution$fit.FEM)

# Errors: Obs - Smoothing
mcm_time_e <- 100 * mainland[,18:29] - eval.FEM.time(mcm_time_solution$fit.FEM.time, locations = mainloc, time.instants = 1:12)

icm_time_solution <- smooth.FEM.time(locations = islandloc, observations = 100 * island[,18:29], islandFEMbasis, 
                                     time_locations = 1:12, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')

# Errors: Obs - Smoothing
icm_time_e <- 100 * island[,18:29] - eval.FEM.time(icm_time_solution$fit.FEM.time, locations = islandloc, time.instants = 1:12)

# Third all differences
mdm_time_solution <- smooth.FEM.time(locations = mainloc, time_locations = 1:12,
                                     observations = 100 * (- mainland[,6:17] + mainland[,18:29]),
                                     mainFEMbasis, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')
#par(mar=c(5,5,5,5))
#plot(mdm_time_solution$fit.FEM)

# Errors: Obs - Smoothing
mdm_time_e <- 100 * (- mainland[,6:17] + mainland[,18:29]) - eval.FEM.time(mdm_time_solution$fit.FEM.time, locations = mainloc, time.instants = 1:12)

idm_time_solution <- smooth.FEM.time(locations = islandloc, time_locations = 1:12,
                                     observations = 100 * (- island[,6:17] + island[,18:29]),
                                     islandFEMbasis, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')

# Errors: Obs - Smoothing
idm_time_e <- 100 * (- island[,6:17] + island[,18:29]) - eval.FEM.time(idm_time_solution$fit.FEM.time, locations = islandloc, time.instants = 1:12)

# Boxplot of errors
par(mar=c(7,5,1,1))
boxplot(as.vector(unlist(mpm_time_e)), as.vector(unlist(ipm_time_e)),
        as.vector(unlist(mcm_time_e)), as.vector(unlist(icm_time_e)),
        as.vector(unlist(mdm_time_e)), as.vector(unlist(idm_time_e)),
        col = c('blue', 'blue', 'green', 'green', 'cyan', 'cyan'), las = 2, ylab = 'Fitting error',
        names = c('Ita 2011-19', 'Sard 2011-19',
                  'Ita 2020', 'Sard 2020',
                  'Ita covid', 'Sard covid'))

#### Monthly smoothing with covariates: WIP -----------

# Fails
raw_complete <- read.csv("raw_complete.csv", header=T)
population <- raw_complete[, c('id', 'pop', 'cod_reg')] %>% unique()
mainpop <- population[population$cod_reg!=20, ]
mainpop2 <- matrix(rep(mainpop[,'pop'], times = 12), nrow=length(mainpop[,'pop']), ncol=12, byrow = FALSE)

# First all pre covid
mpm_cov_solution <- smooth.FEM.time(locations = mainloc, observations = 100 * mainland[,6:17], mainFEMbasis, 
                                     time_locations = 1:12, covariates = mainpop2,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')

ipm_time_solution <- smooth.FEM.time(locations = islandloc, observations = 100 * island[,6:17], islandFEMbasis, 
                                     time_locations = 1:12, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')
#plot(ipm_time_solution$fit.FEM)
# Works

# Second all covid
mcm_time_solution <- smooth.FEM.time(locations = mainloc, observations = 100 * mainland[,18:29], mainFEMbasis, 
                                     time_locations = 1:12, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')
#plot(mcm_time_solution$fit.FEM)

icm_time_solution <- smooth.FEM.time(locations = islandloc, observations = 100 * island[,18:29], islandFEMbasis, 
                                     time_locations = 1:12, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')

# Third all differences
mdm_time_solution <- smooth.FEM.time(locations = mainloc, time_locations = 1:12,
                                     observations = 100 * (- mainland[,6:17] + mainland[,18:29]),
                                     mainFEMbasis, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')
#plot(mdm_time_solution$fit.FEM, locations = mainloc[mainland$name=="Milano",])

idm_time_solution <- smooth.FEM.time(locations = islandloc, time_locations = 1:12,
                                     observations = 100 * (- island[,6:17] + island[,18:29]),
                                     islandFEMbasis, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')

#Errors
beta <- solution$solution$beta
e <- observations - beta[1,1]*covariates - eval.FEM(fun_sol,locations)
hist(e)

#### Biweekly smoothing: Stand by -----------
# Notation:
# 1st letter is i/m island/mainland
# 2nd letter is p/c/d pre-covid/covid/difference
# 3rd letter is m/b/w monthly/biweekly/weekly

# First all pre covid
mpb_time_solution <- smooth.FEM.time(locations = mainloc, observations = 100 * mainland[,6:31], mainFEMbasis, 
                                     time_locations = 1:26, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')

ipb_time_solution <- smooth.FEM.time(locations = islandloc, observations = 100 * island[,6:31], islandFEMbasis, 
                                     time_locations = 1:26, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')
# Works but really slow
plot(ipb_time_solution$fit.FEM.time)

# Second all covid
mcb_time_solution <- smooth.FEM.time(locations = mainloc, observations = 100 * mainland[,32:57], mainFEMbasis, 
                                     time_locations = 1:26, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')

icb_time_solution <- smooth.FEM.time(locations = islandloc, observations = 100 * island[,32:57], islandFEMbasis, 
                                     time_locations = 1:26, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')

# Third all differences
mdb_time_solution <- smooth.FEM.time(locations = mainloc, time_locations = 1:26,
                                     observations = 100 * (- mainland[,6:31] + mainland[,32:57]),
                                     mainFEMbasis, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')

idb_time_solution <- smooth.FEM.time(locations = islandloc, time_locations = 1:26,
                                     observations = 100 * (- island[,6:31] + island[,32:57]),
                                     islandFEMbasis, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')

#### Weekly smoothing: Stand by -----------
# Notation:
# 1st letter is i/m island/mainland
# 2nd letter is p/c/d pre-covid/covid/difference
# 3rd letter is m/b/w monthly/biweekly/weekly

# First all pre covid
mpw_time_solution <- smooth.FEM.time(locations = mainloc, observations = mainland[,6:57], mainFEMbasis, 
                                     time_locations = 1:52, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')

ipw_time_solution <- smooth.FEM.time(locations = islandloc, observations = island[,6:57], islandFEMbasis, 
                                     time_locations = 1:52, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')

# Second all covid
mcw_time_solution <- smooth.FEM.time(locations = mainloc, observations = mainland[,58:109], mainFEMbasis, 
                                     time_locations = 1:52, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')

icw_time_solution <- smooth.FEM.time(locations = islandloc, observations = island[,58:109], islandFEMbasis, 
                                     time_locations = 1:52, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')

# Third all differences
mdw_time_solution <- smooth.FEM.time(locations = mainloc, time_locations = 1:52,
                                     observations = mainland[,6:57] - mainland[,58:109],
                                     mainFEMbasis, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')

idw_time_solution <- smooth.FEM.time(locations = islandloc, time_locations = 1:52,
                                     observations = island[,6:57] - island[,58:109],
                                     islandFEMbasis, covariates = NULL,
                                     PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL,
                                     areal.data.avg = F, DOF.evaluation ='stochastic',
                                     lambdaS = l_space, lambdaT = l_time,
                                     lambda.selection.criterion='grid', lambda.selection.lossfunction = 'GCV')

#### Save results in a .RData file -------------------------------------------

setwd('../Results/Whole Year')
save(mpm_time_solution, mcm_time_solution, mdm_time_solution,
     ipm_time_solution, icm_time_solution, idm_time_solution,
     file = 'solutions_monthly.RData')


###########################################################
##################### RESULTS ANALYSIS ####################
###########################################################

#### Load data + package -----------------------------------------------------


main_path <- 'C:/Users/aless/Documents/NAPDE project'
# main_path <- '/Users/Cheesecake/Desktop/NAPDE/Project'
# main_path <- 'C:/Users/Tommaso/Desktop/NAPDE/Project'
# main_path <- ...

setwd(paste(main_path, '/Resources', sep = ''))
load('Functions.RData')
load('AllBoundaries.RData')
Data <- read.csv('weekly_dataset.csv', encoding = 'UTF-8')

setwd('../Results/First wave')
load('solutions_wave1_newton.RData')


#### Plots of the spatial distribution (NORMALIZED) --------------------------

solutions <- list(mp1_solution, ip1_solution,
                  mc1_solution, ic1_solution,
                  md1_solution, id1_solution)

plot_results(solutions, nx = 500, ny = 500, case = '1st',
             file_name = '1st_wave_static.pdf')

plot_results(solutions, nx = 500, ny = 500, case = '1st', show_lambda = F,
             file_name = '1st_wave_static_NOlambda.pdf')



#### Mainland + Sicilia mesh -------------------------------------------

main_mesh <- create.mesh.2D(italy_boundary, nodesattributes = NA, 
                            segments = italy_segments, 
                            holes = NA,
                            triangles = NA, 
                            order = 1, 
                            verbosity = 0)

# plot(main_mesh)

main_refinement <- refine.mesh.2D(main_mesh, minimum_angle = 25,
                                  maximum_area = 0.01, delaunay = T,
                                  verbosity = 0)

# plot(main_refinement)

#### Sardegna mesh -------------------------------------------

island_mesh <- create.mesh.2D(sardegna_boundary, nodesattributes = NA, 
                              segments = sardegna_segments, 
                              holes = NA,
                              triangles = NA, 
                              order = 1, 
                              verbosity = 0)

# plot(island_mesh)

island_refinement <- refine.mesh.2D(island_mesh, minimum_angle = 25,
                                    maximum_area = 0.01, delaunay = T,
                                    verbosity = 0)

# plot(island_refinement)

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




#### Visualize accuracy ------------------------------------------------------

set.seed(123)
ids_m <- sample(dim(mainloc)[1], size = 40)
ids_i <- sample(dim(islandloc)[1], size = 10)

observations_pre <- c(rowSums(mainland[ids_m, 12:26]),
                      rowSums(island[ids_i, 12:26]))
observations_cov <- c(rowSums(mainland[ids_m, 64:78]),
                      rowSums(island[ids_i, 64:78]))

predictions_pre_m <- eval.FEM(mp1_solution$fit.FEM, 
                              mainloc[ids_m, ])
predictions_pre_i <- eval.FEM(ip1_solution$fit.FEM,
                              islandloc[ids_i, ])
predictions_cov_m <- eval.FEM(mc1_solution$fit.FEM, 
                              mainloc[ids_m, ])
predictions_cov_i <- eval.FEM(ic1_solution$fit.FEM,
                              islandloc[ids_i, ])
predictions_pre <- c(predictions_pre_m, predictions_pre_i)
predictions_cov <- c(predictions_cov_m, predictions_cov_i)

limits_pre <- c(min(min(observations_pre),
                    min(predictions_pre)),
                max(max(observations_pre),
                    max(predictions_pre))
                + 0.001)*100
limits_cov <- c(min(min(observations_cov),
                    min(predictions_cov)),
                max(max(observations_cov),
                    max(predictions_cov))
                + 0.001)*100
labels = c(mainland$name[ids_m], island$name[ids_i])
colors_pre <- NULL
for(i in 1:length(labels)){
  colors_pre[i] <- ifelse(predictions_pre[i] > observations_pre[i],
                          'seagreen3', 'steelblue3')
}
colors_cov <- NULL
for(i in 1:length(labels)){
  colors_cov[i] <- ifelse(predictions_cov[i] > observations_cov[i],
                          'seagreen3', 'steelblue3')
}


par(mar = c(12,5,4,2))
plot(observations_pre*100, pch = 16, col = colors_pre, 
     xlim = c(1, 50), xlab = '', ylab = 'Mortality incidence (%)',
     main = 'Pre COVID prediction errors - First wave',
     xaxt = 'n', ylim = limits_pre)
axis(1, at = 1:length(labels), labels = labels, las = 3)
points(predictions_pre*100, col = colors_pre)
for(i in 1:length(observations_pre)){
  points(c(i, i), c(observations_pre[i], predictions_pre[i])*100,
         type = 'l', col = colors_pre[i])
}
legend('topleft', fill = c('steelblue3', 'seagreen3'),
       legend = c('Underestimated', 'Overestimated'),
       bty = "n")



plot(observations_cov*100, pch = 16, col = colors_cov, 
     xlim = c(1, 50), xlab = '', ylab = 'Mortality incidence (%)',
     main = '2020 Prediction errors - First wave',
     xaxt = 'n', ylim = limits_cov)
axis(1, at = 1:length(labels), labels = labels, las = 3)
points(predictions_cov*100, col = colors_cov)
for(i in 1:length(observations_cov)){
  points(c(i, i), c(observations_cov[i], predictions_cov[i])*100,
         type = 'l', col = colors_cov[i])
}
legend('topleft', fill = c('steelblue3', 'seagreen3'),
       legend = c('Underestimated', 'Overestimated'),
       bty = "n")


# par(mar = c(12,5,4,2))
# plot(-0.1+(1:50), observations_pre*100, pch = 16, col = 'orange', 
#      xlim = c(1, 50), xlab = '', ylab = 'Mortality incidence',
#      main = 'Pre COVID prediction errors - First wave',
#      xaxt = 'n', ylim = limits)
# axis(1, at = 1:length(labels), labels = labels, las = 3)
# points(-0.1+(1:50), predictions_pre*100, col = 'orange')
# for(i in 1:length(observations_pre)){
#   points(-0.1+c(i, i), c(observations_pre[i], predictions_pre[i])*100,
#          type = 'l', col = 'orange')
# }
# 
# points(0.1+(1:50), observations_cov*100, pch = 16, col = 'orangered3')
# axis(1, at = 1:length(labels), labels = labels, las = 3)
# points(0.1+(1:50), predictions_cov*100, col = 'orangered3')
# for(i in 1:length(observations_cov)){
#   points(0.1+c(i, i), c(observations_cov[i], predictions_cov[i])*100,
#          type = 'l', col = 'orangered3')
# }



# Information regarding the error -----------------------------------------

observations_pre <- c(rowSums(mainland[, 12:26]),
                      rowSums(island[, 12:26]))*100
observations_cov <- c(rowSums(mainland[, 64:78]),
                      rowSums(island[, 64:78]))*100

predictions_pre_m <- eval.FEM(mp1_solution$fit.FEM, 
                              mainloc[, ])
predictions_pre_i <- eval.FEM(ip1_solution$fit.FEM,
                              islandloc[, ])
predictions_cov_m <- eval.FEM(mc1_solution$fit.FEM, 
                              mainloc[, ])
predictions_cov_i <- eval.FEM(ic1_solution$fit.FEM,
                              islandloc[, ])
predictions_pre <- c(predictions_pre_m, predictions_pre_i)*100
predictions_cov <- c(predictions_cov_m, predictions_cov_i)*100


# Define errors
errors_pre <- predictions_pre - observations_pre
errors_cov <- predictions_cov - observations_cov

mean(errors_pre)
mean(errors_cov)

mean(errors_pre^2)
mean(errors_cov^2)

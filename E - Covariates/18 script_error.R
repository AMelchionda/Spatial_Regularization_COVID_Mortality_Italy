###########################################################
##################### RESULTS ANALYSIS ####################
###########################################################

# main_path <- 'C:/Users/aless/Documents/NAPDE project'
# main_path <- '/Users/Cheesecake/Desktop/NAPDE/Project'
# main_path <- 'C:/Users/Tommaso/Desktop/NAPDE/Project'
main_path <- ...

setwd(paste(main_path, '/Resources', sep = ''))

#### Load data -----------------------------------------------------

dataset_covariates <- read.csv("dataset_covariates.csv")

#Functions
load("Functions.RData")

#boundary
load("AllBoundaries.RData")

#solution
setwd('../Results/Covariates')
load("solutions_year_covariates_newton.RData")

#nodes
nodes=dataset_covariates[c(4,5)]

####---- mainland mesh ----
mesh=create.mesh.2D(nodes=italy_boundary,segments = italy_segments)
mesh=refine.mesh.2D(mesh,minimum_angle = 25,maximum_area = 0.01)
#fem basis
fem_basis=create.FEM.basis(mesh)

#select the nodes inside the boundary
indexes_out=outside_indexes(nodes,fem_basis)
locations=nodes[-indexes_out,]


#### Visualize accuracy (mainland, 2020, population) ------------------------------------------------------

set.seed(123)

observations=dataset_covariates$mean_norm_deaths_covid[-indexes_out]

ids_m <- sample(length(observations), size = 40)

observations_m <- observations[ids_m]

predictions_m_population <- solution_mp$solution$z_hat[ids_m]

predictions_m_density <- solution_md$solution$z_hat[ids_m]


#population-mainland

limits_cov <- c(min(min(observations_m),
                    min(predictions_m_population)),
                max(max(observations_m),
                    max(predictions_m_population))
                + 0.001)*100

colors_cov <- NULL
for(i in 1:length(labels)){
  colors_cov[i] <- ifelse(predictions_m_population[i] > observations_m[i],
                          'seagreen3', 'steelblue3')
}

labels = c(dataset_covariates$name[ids_m])

plot(observations_m*100, pch = 16, col = colors_cov, 
     xlim = c(1, 40), xlab = '', ylab = 'Mortality incidence (%)',
     main = 'Prediction errors - model with covariates (population)',
     xaxt = 'n', ylim = limits_cov)
axis(1, at = 1:length(labels), labels = labels, las = 3)
points(predictions_m_population*100, col = colors_cov)
for(i in 1:length(observations_m)){
  points(c(i, i), c(observations_m[i], predictions_m_population[i])*100,
         type = 'l', col = colors_cov[i])
}
legend('topleft', fill = c('steelblue3', 'seagreen3'),
       legend = c('Underestimated', 'Overestimated'),
       bty = "n")

# Information regarding the error -----------------------------------------

#
predictions_m_p=solution_mp$solution$z_hat

# Define errors
errors_m_population <- predictions_m_p - observations

mean(errors_m_population)

mean(errors_m_population^2)

#### Visualize accuracy (mainland, 2020, density) ------------------------------------------------------

set.seed(123)

observations=dataset_covariates$mean_norm_deaths_covid[-indexes_out]

ids_m <- sample(length(observations), size = 40)

observations_m <- observations[ids_m]

predictions_m_density <- solution_md$solution$z_hat[ids_m]


#population-mainland

limits_cov <- c(min(min(observations_m),
                    min(predictions_m_density)),
                max(max(observations_m),
                    max(predictions_m_density))
                + 0.001)*100

colors_cov <- NULL
for(i in 1:length(labels)){
  colors_cov[i] <- ifelse(predictions_m_density[i] > observations_m[i],
                          'seagreen3', 'steelblue3')
}

labels = c(dataset_covariates$name[ids_m])

plot(observations_m*100, pch = 16, col = colors_cov, 
     xlim = c(1, 40), xlab = '', ylab = 'Mortality incidence (%)',
     main = 'Prediction errors - model with covariates (density)',
     xaxt = 'n', ylim = limits_cov)
axis(1, at = 1:length(labels), labels = labels, las = 3)
points(predictions_m_density*100, col = colors_cov)
for(i in 1:length(observations_m)){
  points(c(i, i), c(observations_m[i], predictions_m_density[i])*100,
         type = 'l', col = colors_cov[i])
}
legend('topleft', fill = c('steelblue3', 'seagreen3'),
       legend = c('Underestimated', 'Overestimated'),
       bty = "n")

# Information regarding the error -----------------------------------------

#
predictions_m_d=solution_md$solution$z_hat

# Define errors
errors_m_density <- predictions_m_d - observations

mean(errors_m_density)

mean(errors_m_density^2)


####---- island mesh ----
mesh=create.mesh.2D(nodes=sardegna_boundary,segments = sardegna_segments)
mesh=refine.mesh.2D(mesh,minimum_angle = 25,maximum_area = 0.01)
#fem basis
fem_basis=create.FEM.basis(mesh)

#select the nodes inside the boundary
indexes_out=outside_indexes(nodes,fem_basis)
locations=nodes[-indexes_out,]


####

ids_i <- sample(dim(dataset_covariates)[1], size = 10)


observations_cov_i <- dataset_covariates$mean_norm_deaths_covid[ids_i]


predictions_cov_i_population <- solution_ip$solution$z_hat[ids_i]

predictions_cov_i_density <- solution_id$solution$z_hat[ids_i]

#### Visualize accuracy (island, 2020, population) ------------------------------------------------------

set.seed(123)

observations=dataset_covariates$mean_norm_deaths_covid[-indexes_out]

ids_i <- sample(length(observations), size = 10)

observations_i <- observations[ids_i]

predictions_i_population <- solution_ip$solution$z_hat[ids_i]

#population-mainland

limits_cov <- c(min(min(observations_i),
                    min(predictions_i_population)),
                max(max(observations_i),
                    max(predictions_i_population))
                + 0.001)*100

colors_cov <- NULL
for(i in 1:length(labels)){
  colors_cov[i] <- ifelse(predictions_i_population[i] > observations_i[i],
                          'seagreen3', 'steelblue3')
}

labels = c(dataset_covariates$name[ids_i])

plot(observations_i*100, pch = 16, col = colors_cov, 
     xlim = c(1, 10), xlab = '', ylab = 'Mortality incidence (%)',
     main = 'Prediction errors - model with covariates (population)',
     xaxt = 'n', ylim = limits_cov)
axis(1, at = 1:length(labels), labels = labels, las = 3)
points(predictions_i_population*100, col = colors_cov)
for(i in 1:length(observations_i)){
  points(c(i, i), c(observations_i[i], predictions_i_population[i])*100,
         type = 'l', col = colors_cov[i])
}
legend('topleft', fill = c('steelblue3', 'seagreen3'),
       legend = c('Underestimated', 'Overestimated'),
       bty = "n")

# Information regarding the error -----------------------------------------

#
predictions_i_p=solution_ip$solution$z_hat

# Define errors
errors_i_population <- predictions_i_p - observations

mean(errors_i_population)

mean(errors_i_population^2)

#### Visualize accuracy (mainland, 2020, density) ------------------------------------------------------

set.seed(123)

observations=dataset_covariates$mean_norm_deaths_covid[-indexes_out]

ids_i <- sample(length(observations), size = 10)

observations_i <- observations[ids_i]

predictions_i_density <- solution_id$solution$z_hat[ids_i]


#population-mainland

limits_cov <- c(min(min(observations_i),
                    min(predictions_i_density)),
                max(max(observations_i),
                    max(predictions_i_density))
                + 0.001)*100

colors_cov <- NULL
for(i in 1:length(labels)){
  colors_cov[i] <- ifelse(predictions_i_density[i] > observations_i[i],
                          'seagreen3', 'steelblue3')
}

labels = c(dataset_covariates$name[ids_i])

plot(observations_i*100, pch = 16, col = colors_cov, 
     xlim = c(1, 10), xlab = '', ylab = 'Mortality incidence (%)',
     main = 'Prediction errors - model with covariates (density)',
     xaxt = 'n', ylim = limits_cov)
axis(1, at = 1:length(labels), labels = labels, las = 3)
points(predictions_i_density*100, col = colors_cov)
for(i in 1:length(observations_i)){
  points(c(i, i), c(observations_i[i], predictions_i_density[i])*100,
         type = 'l', col = colors_cov[i])
}
legend('topleft', fill = c('steelblue3', 'seagreen3'),
       legend = c('Underestimated', 'Overestimated'),
       bty = "n")

# Information regarding the error -----------------------------------------

#
predictions_i_d=solution_id$solution$z_hat

# Define errors
errors_i_density <- predictions_i_d - observations

mean(errors_i_density)

mean(errors_i_density^2)


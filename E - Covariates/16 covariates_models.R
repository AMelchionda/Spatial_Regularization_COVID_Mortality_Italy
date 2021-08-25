###################################
######## COVARIATES MODEL #########
###################################

# Pre setting + data loading ----

# main_path <- 'C:/Users/aless/Documents/NAPDE project'
# main_path <- '/Users/Cheesecake/Desktop/NAPDE/Project'
# main_path <- 'C:/Users/Tommaso/Desktop/NAPDE/Project'
main_path <- ...

setwd(paste(main_path, '/Resources', sep = ''))

library(fdaPDE)
library(readr)

load("Functions.RData")
load("AllBoundaries.RData")
dataset_covariates <- read_csv("dataset_covariates.csv",
                               col_types = cols(id = col_integer(), 
                                                name = col_character(),
                                                cod_reg = col_integer(),
                                                population = col_double()))

#nodes
nodes=dataset_covariates[c(4,5)]

setwd('../')



####---- mainland mesh ----
mesh=create.mesh.2D(nodes=italy_boundary,segments = italy_segments)
mesh=refine.mesh.2D(mesh,minimum_angle = 25,maximum_area = 0.01)
#fem basis
fem_basis=create.FEM.basis(mesh)

#select the nodes inside the boundary
indexes_out=outside_indexes(nodes,fem_basis)
locations=nodes[-indexes_out,]


####---- smoothing with covariates (population, 2020, mainland)----

#set the variables
observations=as.matrix(dataset_covariates$mean_norm_deaths_covid[-indexes_out])#(<-norm_death_covid)

covariates=as.matrix(dataset_covariates$population[-indexes_out])#(<-population) 

#smoothing
solution_mp=smooth.FEM(locations = locations,
                       observations = observations,
                       FEMbasis = fem_basis,
                       covariates = covariates,
                       lambda.selection.criterion='newton')

####---- smoothing with covariates (density, 2020, mainland) ----

#set the variables
observations=as.matrix(dataset_covariates$mean_norm_deaths_covid[-indexes_out])#(<-norm_death_covid)

covariates=as.matrix(dataset_covariates$density[-indexes_out])#(<-density) 

#smoothing
solution_md=smooth.FEM(locations = locations,
                       observations = observations,
                       FEMbasis = fem_basis,
                       covariates = covariates,
                       lambda.selection.criterion='newton')


####---- island mesh ----
mesh=create.mesh.2D(nodes=sardegna_boundary,segments = sardegna_segments)
mesh=refine.mesh.2D(mesh,minimum_angle = 25,maximum_area = 0.01)
#fem basis
fem_basis=create.FEM.basis(mesh)

#select the nodes inside the boundary
indexes_out=outside_indexes(nodes,fem_basis)
locations=nodes[-indexes_out,]


####---- smoothing with covariates (population, 2020, island)----

#set the variables
observations=as.matrix(dataset_covariates$mean_norm_deaths_covid[-indexes_out])#(<-norm_death_covid)

covariates=as.matrix(dataset_covariates$population[-indexes_out])#(<-population) 

#smoothing
solution_ip=smooth.FEM(locations = locations,
                       observations = observations,
                       FEMbasis = fem_basis,
                       covariates = covariates,
                       lambda.selection.criterion = "newton")


####---- smoothing with covariates (density, 2020, island)----

#set the variables
observations=as.matrix(dataset_covariates$mean_norm_deaths_covid[-indexes_out])#(<-norm_death_covid)

covariates=as.matrix(dataset_covariates$density[-indexes_out])#(<-density) 

#smoothing
solution_id=smooth.FEM(locations = locations,
                       observations = observations,
                       FEMbasis = fem_basis,
                       covariates = covariates,
                       lambda.selection.criterion = "newton")


#### Save the solutions ----
setwd('./Results/Covariates')
save(solution_mp, solution_md, 
     solution_ip, solution_id, 
     file = 'solutions_year_covariates_newton.RData')


# the solutions name follows the the following code:
# m/i <-mainland/island
# p/d <-population/density


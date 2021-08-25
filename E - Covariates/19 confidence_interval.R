###########################################################
######## CONFIDENCE INTERVAL FOR COVARIATES MODEL #########
###########################################################

# Pre setting ----

# main_path <- 'C:/Users/aless/Documents/NAPDE project'
# main_path <- '/Users/Cheesecake/Desktop/NAPDE/Project'
# main_path <- 'C:/Users/Tommaso/Desktop/NAPDE/Project'
main_path <- ...

setwd(paste(main_path, '/Resources', sep = ''))

load('Functions.RData')

#install.packages("geometry")
#install.packages("rgl")
#install.packages("plot3D")
#install.packages("plot3Drgl")
#install.packages("shiny")
#install.packages("RcppEigen")
# Update package with GitHub: Download and install from local file
#install.packages("fdaPDE-dev", repo=NULL, type="source")

library(fdaPDE)

#### Load the data ----

load("AllBoundaries.RData")
dataset_covariates <- read.csv("dataset_covariates.csv")  # load the dataset with the covariates
nodes <- dataset_covariates[c(4,5)]  # nodes

#### Mainland mesh ----------
mesh = create.mesh.2D(nodes=italy_boundary,segments = italy_segments)  # space mesh
mesh = refine.mesh.2D(mesh,minimum_angle = 25,maximum_area = 0.01)
fem_basis = create.FEM.basis(mesh)  # fem basis

indexes_out = outside_indexes(nodes,fem_basis)  # select the nodes inside the boundary
locations = nodes[-indexes_out,]

#### Smoothing with population ----

observations = as.matrix(dataset_covariates$mean_norm_deaths_covid[-indexes_out])  # norm_death_covid
covariates = as.matrix(dataset_covariates$population[-indexes_out])  # population

#?inferenceDataObjectBuilder
infer_obj <- inferenceDataObjectBuilder(test = "simultaneous",
                                        interval = "bonferroni",
                                        type = "wald", exact = "True",  # exact=T always
                                        dim = 1)

output_CPP <- smooth.FEM(locations = locations, observations = observations,
                         covariates = covariates, FEMbasis = fem_basis,
                         lambda.selection.lossfunction = "GCV",
                         R_Inference_Data_Object = infer_obj)
main_pop_pval <- output_CPP$inference$p_values$wald[[1]]
main_pop_ci <- output_CPP$inference$CI$wald[[1]]

#### Smoothing with density ----

observations = as.matrix(dataset_covariates$mean_norm_deaths_covid[-indexes_out])  # norm_death_covid
covariates = as.matrix(dataset_covariates$density[-indexes_out])  # density

infer_obj <- inferenceDataObjectBuilder(test = "simultaneous",
                                        interval = "bonferroni",
                                        type = "wald",
                                        exact = "True",  # exact=T always
                                        dim = 1)

output_CPP <- smooth.FEM(locations = locations, observations = observations,
                         covariates = covariates, FEMbasis = fem_basis,
                         lambda.selection.lossfunction = "GCV",
                         R_Inference_Data_Object = infer_obj)
main_den_pval <- output_CPP$inference$p_values$wald[[1]]
main_den_ci <- output_CPP$inference$CI$wald[[1]]

#### Island mesh ----------


mesh = create.mesh.2D(nodes=sardegna_boundary,segments = sardegna_segments)  # space mesh
mesh = refine.mesh.2D(mesh,minimum_angle = 25,maximum_area = 0.01)
fem_basis=create.FEM.basis(mesh)  # fem basis

indexes_out = outside_indexes(nodes,fem_basis)  # nodes inside the boundary
locations = nodes[-indexes_out,]

#### Smoothing with population ----

observations = as.matrix(dataset_covariates$mean_norm_deaths_covid[-indexes_out])
covariates = as.matrix(dataset_covariates$population[-indexes_out])

infer_obj <- inferenceDataObjectBuilder(test = "simultaneous",
                                        interval = "bonferroni",
                                        type = "wald", exact = "True",  # exact=T always
                                        dim = 1)

output_CPP <- smooth.FEM(locations = locations, observations = observations,
                         covariates = covariates, FEMbasis = fem_basis,
                         lambda.selection.lossfunction = "GCV",
                         R_Inference_Data_Object = infer_obj)

is_pop_pval <- output_CPP$inference$p_values$wald[[1]]
is_pop_ci <- output_CPP$inference$CI$wald[[1]]

#### Smoothing with density ----

observations = as.matrix(dataset_covariates$mean_norm_deaths_covid[-indexes_out])
covariates = as.matrix(dataset_covariates$density[-indexes_out])

infer_obj <- inferenceDataObjectBuilder(test = "simultaneous",
                                        interval = "bonferroni",
                                        type = "wald", exact = "True",  # exact=T always
                                        dim = 1)

output_CPP <- smooth.FEM(locations = locations, observations = observations,
                         covariates = covariates, FEMbasis = fem_basis,
                         lambda.selection.lossfunction = "GCV",
                         R_Inference_Data_Object = infer_obj)

is_den_pval <- output_CPP$inference$p_values$wald[[1]]
is_den_ci <- output_CPP$inference$CI$wald[[1]]

# Save everything ---- 

setwd('../Results/Covariates')
save(main_pop_pval, main_pop_ci, main_den_pval, main_den_ci,
     is_pop_pval, is_pop_ci, is_den_pval, is_den_ci,
     file = 'covars_inference_testing.RData')

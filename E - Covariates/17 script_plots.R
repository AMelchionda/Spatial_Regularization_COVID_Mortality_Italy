
# Pre setting ----

# main_path <- 'C:/Users/aless/Documents/NAPDE project'
# main_path <- '/Users/Cheesecake/Desktop/NAPDE/Project'
# main_path <- 'C:/Users/Tommaso/Desktop/NAPDE/Project'
main_path <- ...

setwd(paste(main_path, '/Resources', sep = ''))


library(fdaPDE)
library(readr)
library(ggplot2)
library(gridExtra)
library(scales)

#### Load the data ----

#boundary
load("AllBoundaries.RData")
load("AllBoundariesPlot.RData")

#load the dataset
dataset_covariates <- read_csv("dataset_covariates.csv",
                               col_types = cols(id = col_integer(), name = col_character(),
                               cod_reg = col_integer(), population = col_double()))
#nodes
nodes=dataset_covariates[c(4,5)]

#functions
load("Functions.RData")

#solutions
setwd('../Results/Covariates')
load("solutions_year_covariates_newton.RData")


####################################### plot ###################################

#density
solution=list(solution_md,solution_id)

pdf('density.pdf',paper = 'letter')
plot_results_covariates(solution, sardegna_boundary_plot, sardegna_boundary,
                        italy_boundary_plot, italy_boundary, nodes,
                        xmin=NA, xmax=NA, ymin=NA, ymax=NA, nx=500, ny=500, 
                        title = '', limits=NULL,
                        avoid_bar = FALSE, bar_height = 20,
                        col_values = NULL, difference = FALSE,
                        show_lambdas = TRUE, show_betas = TRUE, flag1 = FALSE,
                        flag2 = FALSE)
dev.off()

#population
solution=list(solution_mp,solution_ip)

pdf('population.pdf',paper = 'letter')
plot_results_covariates(solution, sardegna_boundary_plot, sardegna_boundary,
                        italy_boundary_plot, italy_boundary, nodes,
                        xmin=NA, xmax=NA, ymin=NA, ymax=NA, nx=500, ny=500, 
                        title = '', limits=NULL,
                        avoid_bar = FALSE, bar_height = 20,
                        col_values = NULL, difference = FALSE,
                        show_lambdas = TRUE, show_betas = TRUE, flag1 = FALSE,
                        flag2 = FALSE)
dev.off()
 

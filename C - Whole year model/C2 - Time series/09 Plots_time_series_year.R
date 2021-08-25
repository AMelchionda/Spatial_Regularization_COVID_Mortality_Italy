


# main_path <- 'C:/Users/aless/Documents/NAPDE project'
# main_path <- '/Users/Cheesecake/Desktop/NAPDE/Project'
# main_path <- 'C:/Users/Tommaso/Desktop/NAPDE/Project'
main_path <- ...

setwd(paste(main_path, '/Resources', sep = ''))
load('AllBoundaries.RData')
load('Functions.RData')
Data <- read.csv('weekly_dataset.csv', encoding = 'UTF-8')

setwd('../Results/Whole year')
load('solutions_monthly.RData')

library(fdaPDE)

l_space <- 10^-5
l_time <- 10^-5


# Plotting the results ----------------------------------------------------

solutions <- list(mpm_time_solution, ipm_time_solution,
                  mcm_time_solution, icm_time_solution)

# Produce GIF
produce_animation(solutions, case = 'year', file_name='whole_year_animation')

# Produce timestamps
# timestamps <- plot_timestamps(solutions, case = 'year',
#                               file_name = 'whole_year_animation')
# plot_combined <- grid.arrange(timestamps[[1]], timestamps[[2]],
#                               timestamps[[5]], timestamps[[6]],
#                               timestamps[[9]], timestamps[[10]],
#                               timestamps[[13]], timestamps[[14]],
#                               timestamps[[17]], timestamps[[19]],
#                               timestamps[[21]], timestamps[[22]],
#                               timestamps[[25]], timestamps[[26]],
#                               timestamps[[29]], timestamps[[30]],
#                               ncol = 4)
# ggsave('Results/whole_year_timestamp.pdf', plot_combined, width=16, height=13.5) 


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



# Plot municipalities behaviour -------------------------------------------



set.seed(222)
mun <- c('Milano', 'Pesaro', 'Roma', 'Napoli')
for(i in 1:4){
  ids[i] <- which(mainland$name == mun[i])
}
ids[5:12] <- sample(1:dim(mainloc)[1], size = 8)

par(mfrow = c(3,4))
for(id in ids){
  pred1 <- eval.FEM.time(mpm_time_solution$fit.FEM.time,
                        locations = mainloc[id,],
                        time.instants = 1+(0:100)*(11/100))
  
  pred2 <- eval.FEM.time(mcm_time_solution$fit.FEM.time,
                        locations = mainloc[id,],
                        time.instants = 1+(0:100)*(11/100))
  
  limits <- c(min(min(pred1, pred2), min(mainland[id, 6:29]*100)),
              max(max(pred1, pred2), max(mainland[id, 6:29]*100)))
  
  plot(1+(0:100)*(11/100), pred1, type = 'l', col = 'orange',
       main = mainland$name[id], 
       xlab = '', ylab = '', ylim = limits, lwd = 2)
  points(1:12, mainland[id, 6:17]*100, pch = 16, col = 'orange')
  points(1+(0:100)*(11/100), pred2, type = 'l', col = 'orangered3',
       main = mainland$name[id], 
       xlab = '', ylab = '', ylim = limits, lwd = 2)
  points(1:12, mainland[id, 18:29]*100, pch = 16, col = 'orangered3')
}


# Plot errors -------------------------------------------------------------

graphics.off()
epm <- ecm <- edm <- numeric(dim(mainloc)[1])
epi <- eci <- edi <- numeric(dim(islandloc)[1])
for(i in 1:dim(mainloc)[1]){
  print(paste(i, '/', dim(mainloc)[1]), sep = '')
  pred <- eval.FEM.time(mpm_time_solution$fit.FEM.time,
                        locations = mainloc[i,],
                        time.instants = 1:12)
  epm[i] <- sum((pred - mainland[i, 6:17])^2)/12
  pred <- eval.FEM.time(mcm_time_solution$fit.FEM.time,
                        locations = mainloc[i,],
                        time.instants = 1:12)
  ecm[i] <- sum((pred - mainland[i, 18:29])^2)/12
  pred <- eval.FEM.time(mdm_time_solution$fit.FEM.time,
                        locations = mainloc[i,],
                        time.instants = 1:12)
  edm[i] <- sum((pred - mainland[i, 18:29] + mainland[i, 6:17])^2)
}
for(i in 1:dim(islandloc)[1]){
  print(paste(i, '/', dim(islandloc)[1]), sep = '')
  pred <- eval.FEM.time(ipm_time_solution$fit.FEM.time,
                        locations = islandloc[i,],
                        time.instants = 1:12)
  epi[i] <- sum((pred - island[i, 6:17])^2)/12
  pred <- eval.FEM.time(icm_time_solution$fit.FEM.time,
                        locations = islandloc[i,],
                        time.instants = 1:12)
  eci[i] <- sum((pred - island[i, 18:29])^2)/12
  pred <- eval.FEM.time(idm_time_solution$fit.FEM.time,
                        locations = islandloc[i,],
                        time.instants = 1:12)
  edi[i] <- sum((pred - island[i, 18:29] + island[i, 6:17])^2)
}

layout(cbind(c(2,1), c(3,1)))
boxplot(log10(epm), log10(epi),
        log10(ecm), log10(eci),
        log10(edm), log10(edi), 
        col = c('orange', 'orange',
                'orangered3', 'orangered3',
                'grey', 'grey'),
        names = c('Mainland\naverage','Sardegna\naverage', 
                  'Mainland\n2020', 'Sardegna\n2020', 
                  'Mainland\ndifference', 'Sardegna\ndifference'),
        main = '',
        ylab = 'log10(err^2)')

transp_col <- c(rgb(255/255, 165/255, 000/255, 0.55),
                rgb(205/255, 055/255, 000/255, 0.55),
                rgb(190/255, 190/255, 190/255, 0.55))

hist(log10(ecm), col = transp_col[2], breaks = 20, xlim = c(-3, 0),
     ylab = '', xlab = '', main = 'Mainland')
hist(log10(epm), col = transp_col[1], breaks = 20, add = T)
hist(log10(edm), col = transp_col[3], breaks = 30, add = T)

hist(log10(eci), col = transp_col[2], breaks = 10, xlim = c(-3, 0),
     ylab = '', xlab = '', main = 'Sardegna')
hist(log10(epi), col = transp_col[1], breaks = 10, add = T)
hist(log10(edi), col = transp_col[3], breaks = 30, add = T)

mtext("Errors distribution in logarithmic scale", 
      side = 3, line = -1.5, outer = TRUE)
    
##############################################################################
############################ EXPLORATORY ANALYSIS ############################
##############################################################################


#### Load the data -----------------------------------------------------------

# main_path <- 'C:/Users/aless/Documents/NAPDE project'
# main_path <- '/Users/Cheesecake/Desktop/NAPDE/Project'
# main_path <- 'C:/Users/Tommaso/Desktop/NAPDE/Project'
main_path <- ...

setwd(paste(main_path, '/Resources', sep = ''))
load('Functions.RData')
load('AllBoundaries.RData')
CompleteData <- read.csv('raw_complete.csv', encoding = 'UTF-8')

setwd('../')

#### ... ####
#### ... ####
#### ... ####
# Plot boundary -----------------------------------------------------------
library(fdaPDE)

italy_boundary <- as.matrix(italy_boundary)
sardegna_boundary <- as.matrix(sardegna_boundary)

x_it = c(italy_boundary[,1], italy_boundary[1,1])
y_it = c(italy_boundary[,2], italy_boundary[1,2])
x_is = c(sardegna_boundary[,1], sardegna_boundary[1,1])
y_is = c(sardegna_boundary[,2], sardegna_boundary[1,2])

locations <- cbind(CompleteData$x, CompleteData$y)
locations <- unique(locations)

main_mesh <- create.mesh.2D(italy_boundary,
                            segments = italy_segments)
island_mesh <- create.mesh.2D(sardegna_boundary,
                              segments = sardegna_segments)
mainFEMbasis <- create.FEM.basis(main_mesh)
islandFEMbasis <- create.FEM.basis(island_mesh)

idx_in <- NULL
idx_in <- inside_indexes(locations, mainFEMbasis)
idx_in <- c(idx_in, inside_indexes(locations, islandFEMbasis))

plot(locations[idx_in, 1], locations[idx_in, 2], pch = 19, asp = 1, cex = 0.3, 
     col = 'gray18', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
points(locations[-idx_in, 1], locations[-idx_in, 2], pch = 19, col = 'orangered3')
points(x_it, y_it, type = 'l', asp = 1, col = 'dodgerblue', lwd = 2)
points(x_is, y_is, type = 'l', col = 'dodgerblue', lwd = 2)
legend('bottomleft', fill = c('dodgerblue', 'black', 'orangered3'), 
       legend = c('Boundary', 'Considered municipalities', 'Discarded municipalities'),
       bty = 'n')
# save as 7 x 9 in

#### ... ####
#### ... ####
#### ... ####
#### Plot daily data for Milano, Morimondo, Nosate (colored pre COVID) -------------------------------------------

hcl.pals(type = 'qualitative')
colors <- hcl.colors(10, palette = 'Dark 3', rev = T)
colors[10] <- 'red3'
legend_labels <- as.character(2011:2020)

months_length <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
corresp_days <- cumsum(months_length)
days_loc <- c(0, corresp_days[c(2, 4, 8, 10, 12)])

graphics.off()
par(mfrow = c(1, 3))
# par(oma = c(4,1,1,1), mfrow = c(1, 3), mar = c(1, 1, 1, 1))

mun <- 'Milano'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

plot(1:366, Data[1,], type = 'l', col = colors[1], 
     ylim = c(0, max(Data)), xlim = c(0,366),
     main = mun, xlab = '', ylab ='Mortality',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:366, Data[i,], type = 'l', col = colors[i])
}
points(1:366, Data[10,], type = 'l', col = colors[10], lwd = 2)


mun <- 'Morimondo'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

plot(1:366, Data[1,], type = 'l', col = colors[1], 
     ylim = c(0, max(Data)), xlim = c(0,366),
     main = mun, xlab = '', ylab ='',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:366, Data[i,], type = 'l', col = colors[i])
}
points(1:366, Data[10,], type = 'l', col = colors[10], lwd = 2)


mun <- 'Nosate'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

plot(1:366, Data[1,], type = 'l', col = colors[1], 
     ylim = c(0, max(Data)), xlim = c(0,366),
     main = mun, xlab = '', ylab ='',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:366, Data[i,], type = 'l', col = colors[i])
}
points(1:366, Data[10,], type = 'l', col = colors[10], lwd = 2)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend_labels, fill = colors, xpd = TRUE, bty = 'n', ncol = 5)




#### Plot daily data for Milano, Morimondo, Nosate (gray pre COVID) -------------------------------------------

hcl.pals(type = 'qualitative')
colors <- c('gray', 'red3')
legend_labels <- c('2011 - 2019', '2020')

months_length <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
corresp_days <- cumsum(months_length)
days_loc <- c(0, corresp_days[c(2, 4, 8, 10, 12)])

graphics.off()
par(mfrow = c(1, 3))
# par(oma = c(4,1,1,1), mfrow = c(1, 3), mar = c(1, 1, 1, 1))

mun <- 'Milano'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

plot(1:366, Data[1,], type = 'l', col = colors[1], 
     ylim = c(0, max(Data)), xlim = c(0,366),
     main = mun, xlab = '', ylab ='Mortality',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:366, Data[i,], type = 'l', col = colors[1])
}
points(1:366, Data[10,], type = 'l', col = colors[2], lwd = 2)



mun <- 'Morimondo'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

plot(1:366, Data[1,], type = 'l', col = colors[1], 
     ylim = c(0, max(Data)), xlim = c(0,366),
     main = mun, xlab = '', ylab ='',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:366, Data[i,], type = 'l', col = colors[1])
}
points(1:366, Data[10,], type = 'l', col = colors[2], lwd = 2)


mun <- 'Nosate'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

plot(1:366, Data[1,], type = 'l', col = colors[1], 
     ylim = c(0, max(Data)), xlim = c(0,366),
     main = mun, xlab = '', ylab ='',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:366, Data[i,], type = 'l', col = colors[1])
}
points(1:366, Data[10,], type = 'l', col = colors[2], lwd = 2)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend_labels, fill = colors, xpd = TRUE,
       bty = 'n', ncol = 2)


#### Plot moving average for Milano, Morimondo, Nosate (colored pre COVID) ----------------------------

library(zoo)
hcl.pals(type = 'qualitative')
colors <- hcl.colors(10, palette = 'Dark 3', rev = T)
colors[10] <- 'red3'
legend_labels <- as.character(2011:2020)

months_length <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
corresp_days <- cumsum(months_length)
days_loc <- c(0, corresp_days[c(2, 4, 8, 10, 12)] - 7)

graphics.off()
par(mfrow = c(1, 3))

smoothing_interval <- 14

mun <- 'Milano'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

Data_movmean <- t(rollmean(t(Data), smoothing_interval))
k <- dim(Data_movmean)[2]

plot(1:k, Data_movmean[1,], type = 'l', col = colors[1],
     ylim = c(0, max(Data_movmean)),
     main = mun, xlab = '', ylab ='Mortality',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:k, Data_movmean[i,], type = 'l', col = colors[i])
}
points(1:k, Data_movmean[10,], type = 'l', col = colors[10], lwd = 3)


mun <- 'Morimondo'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

Data_movmean <- t(rollmean(t(Data), smoothing_interval))
k <- dim(Data_movmean)[2]

plot(1:k, Data_movmean[1,], type = 'l', col = colors[1],
     ylim = c(0, max(Data_movmean)),
     main = mun, xlab = '', ylab ='Mortality',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:k, Data_movmean[i,], type = 'l', col = colors[i])
}
points(1:k, Data_movmean[10,], type = 'l', col = colors[10], lwd = 3)



mun <- 'Nosate'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

Data_movmean <- t(rollmean(t(Data), smoothing_interval))
k <- dim(Data_movmean)[2]

plot(1:k, Data_movmean[1,], type = 'l', col = colors[1],
     ylim = c(0, max(Data_movmean)),
     main = mun, xlab = '', ylab ='Mortality',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:k, Data_movmean[i,], type = 'l', col = colors[i])
}
points(1:k, Data_movmean[10,], type = 'l', col = colors[10], lwd = 3)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend_labels, fill = colors, xpd = TRUE, bty = 'n', ncol = 5)



#### Plot moving average for Milano, Morimondo, Nosate (colored pre COVID) ----------------------------

library(zoo)
hcl.pals(type = 'qualitative')
colors <- c('gray', 'red3')
legend_labels <- as.character(2011:2020)

months_length <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
corresp_days <- cumsum(months_length)
days_loc <- c(0, corresp_days[c(2, 4, 8, 10, 12)] - 7)

graphics.off()
par(mfrow = c(1, 3))

smoothing_interval <- 14

mun <- 'Milano'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

Data_movmean <- t(rollmean(t(Data), smoothing_interval))
k <- dim(Data_movmean)[2]

plot(1:k, Data_movmean[1,], type = 'l', col = colors[1],
     ylim = c(0, max(Data_movmean)),
     main = mun, xlab = '', ylab ='Mortality',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:k, Data_movmean[i,], type = 'l', col = colors[1])
}
points(1:k, Data_movmean[10,], type = 'l', col = colors[2], lwd = 2)


mun <- 'Morimondo'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

Data_movmean <- t(rollmean(t(Data), smoothing_interval))
k <- dim(Data_movmean)[2]

plot(1:k, Data_movmean[1,], type = 'l', col = colors[1],
     ylim = c(0, max(Data_movmean)),
     main = mun, xlab = '', ylab ='Mortality',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:k, Data_movmean[i,], type = 'l', col = colors[1])
}
points(1:k, Data_movmean[10,], type = 'l', col = colors[2], lwd = 2)



mun <- 'Nosate'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

Data_movmean <- t(rollmean(t(Data), smoothing_interval))
k <- dim(Data_movmean)[2]

plot(1:k, Data_movmean[1,], type = 'l', col = colors[1],
     ylim = c(0, max(Data_movmean)),
     main = mun, xlab = '', ylab ='Mortality',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:k, Data_movmean[i,], type = 'l', col = colors[1])
}
points(1:k, Data_movmean[10,], type = 'l', col = colors[2], lwd = 2)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend_labels, fill = colors, xpd = TRUE, bty = 'n', ncol = 5)



#### ... ####
#### ... ####
#### ... ####
#### Couple plots for Milano (colored pre COVID) -------------------------------------------

hcl.pals(type = 'qualitative')
colors <- hcl.colors(10, palette = 'Dark 3', rev = T)
colors[10] <- 'red3'
legend_labels <- as.character(2011:2020)

months_length <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
corresp_days <- cumsum(months_length)
days_loc <- c(0, corresp_days[c(2, 4, 8, 10, 12)])

graphics.off()
par(mfrow = c(1, 2))
# par(oma = c(4,1,1,1), mfrow = c(1, 3), mar = c(1, 1, 1, 1))

mun <- 'Milano'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

plot(1:366, Data[1,], type = 'l', col = colors[1], 
     ylim = c(0, max(Data)), xlim = c(0,366),
     xlab = '', ylab ='Mortality',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:366, Data[i,], type = 'l', col = colors[i])
}
points(1:366, Data[10,], type = 'l', col = colors[10], lwd = 2)


library(zoo)
months_length <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
corresp_days <- cumsum(months_length)
days_loc <- c(0, corresp_days[c(2, 4, 8, 10, 12)] - 7)

smoothing_interval <- 14

Data_movmean <- t(rollmean(t(Data), smoothing_interval))
k <- dim(Data_movmean)[2]

plot(1:k, Data_movmean[1,], type = 'l', col = colors[1],
     ylim = c(0, max(Data_movmean)),
     xlab = '', ylab ='',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:k, Data_movmean[i,], type = 'l', col = colors[i])
}
points(1:k, Data_movmean[10,], type = 'l', col = colors[10], lwd = 3)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend_labels, fill = colors, xpd = TRUE, bty = 'n', ncol = 5)




#### Couple plots for Milano (gray pre COVID) -------------------------------------------

hcl.pals(type = 'qualitative')
colors <- c('gray', 'red3')
legend_labels <- c('2011 - 2019', '2020')

months_length <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
corresp_days <- cumsum(months_length)
days_loc <- c(0, corresp_days[c(2, 4, 8, 10, 12)])

graphics.off()
par(mfrow = c(1, 2))
# par(oma = c(4,1,1,1), mfrow = c(1, 3), mar = c(1, 1, 1, 1))

mun <- 'Milano'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

plot(1:366, Data[1,], type = 'l', col = colors[1], 
     ylim = c(0, max(Data)), xlim = c(0,366),
     xlab = '', ylab ='Mortality',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:366, Data[i,], type = 'l', col = colors[1])
}
points(1:366, Data[10,], type = 'l', col = colors[2], lwd = 2)



library(zoo)
months_length <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
corresp_days <- cumsum(months_length)
days_loc <- c(0, corresp_days[c(2, 4, 8, 10, 12)] - 7)

smoothing_interval <- 14

Data_movmean <- t(rollmean(t(Data), smoothing_interval))
k <- dim(Data_movmean)[2]

plot(1:k, Data_movmean[1,], type = 'l', col = colors[1],
     ylim = c(0, max(Data_movmean)),
     xlab = '', ylab ='',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:k, Data_movmean[i,], type = 'l', col = colors[1])
}
points(1:k, Data_movmean[10,], type = 'l', col = colors[2], lwd = 2)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend_labels, fill = colors, xpd = TRUE,
       bty = 'n', ncol = 2)


#### Couple plots for Morimondo (colored pre COVID) -------------------------------------------

hcl.pals(type = 'qualitative')
colors <- hcl.colors(10, palette = 'Dark 3', rev = T)
colors[10] <- 'red3'
legend_labels <- as.character(2011:2020)

months_length <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
corresp_days <- cumsum(months_length)
days_loc <- c(0, corresp_days[c(2, 4, 8, 10, 12)])

graphics.off()
par(mfrow = c(1, 2))
# par(oma = c(4,1,1,1), mfrow = c(1, 3), mar = c(1, 1, 1, 1))

mun <- 'Morimondo'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

plot(1:366, Data[1,], type = 'l', col = colors[1], 
     ylim = c(0, max(Data)), xlim = c(0,366),
     xlab = '', ylab ='Mortality',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:366, Data[i,], type = 'l', col = colors[i])
}
points(1:366, Data[10,], type = 'l', col = colors[10], lwd = 2)


library(zoo)
months_length <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
corresp_days <- cumsum(months_length)
days_loc <- c(0, corresp_days[c(2, 4, 8, 10, 12)] - 7)

smoothing_interval <- 14

Data_movmean <- t(rollmean(t(Data), smoothing_interval))
k <- dim(Data_movmean)[2]

plot(1:k, Data_movmean[1,], type = 'l', col = colors[1],
     ylim = c(0, max(Data_movmean)),
     xlab = '', ylab ='',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:k, Data_movmean[i,], type = 'l', col = colors[i])
}
points(1:k, Data_movmean[10,], type = 'l', col = colors[10], lwd = 3)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend_labels, fill = colors, xpd = TRUE, bty = 'n', ncol = 5)




#### Couple plots for Morimondo (gray pre COVID) -------------------------------------------

hcl.pals(type = 'qualitative')
colors <- c('gray', 'red3')
legend_labels <- c('2011 - 2019', '2020')

months_length <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
corresp_days <- cumsum(months_length)
days_loc <- c(0, corresp_days[c(2, 4, 8, 10, 12)])

graphics.off()
par(mfrow = c(1, 2))
# par(oma = c(4,1,1,1), mfrow = c(1, 3), mar = c(1, 1, 1, 1))

mun <- 'Morimondo'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

plot(1:366, Data[1,], type = 'l', col = colors[1], 
     ylim = c(0, max(Data)), xlim = c(0,366),
     xlab = '', ylab ='Mortality',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:366, Data[i,], type = 'l', col = colors[1])
}
points(1:366, Data[10,], type = 'l', col = colors[2], lwd = 2)



library(zoo)
months_length <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
corresp_days <- cumsum(months_length)
days_loc <- c(0, corresp_days[c(2, 4, 8, 10, 12)] - 7)

smoothing_interval <- 14

Data_movmean <- t(rollmean(t(Data), smoothing_interval))
k <- dim(Data_movmean)[2]

plot(1:k, Data_movmean[1,], type = 'l', col = colors[1],
     ylim = c(0, max(Data_movmean)),
     xlab = '', ylab ='',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:k, Data_movmean[i,], type = 'l', col = colors[1])
}
points(1:k, Data_movmean[10,], type = 'l', col = colors[2], lwd = 2)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend_labels, fill = colors, xpd = TRUE,
       bty = 'n', ncol = 2)


#### Couple plots for Nosate (colored pre COVID) -------------------------------------------

hcl.pals(type = 'qualitative')
colors <- hcl.colors(10, palette = 'Dark 3', rev = T)
colors[10] <- 'red3'
legend_labels <- as.character(2011:2020)

months_length <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
corresp_days <- cumsum(months_length)
days_loc <- c(0, corresp_days[c(2, 4, 8, 10, 12)])

graphics.off()
par(mfrow = c(1, 2))
# par(oma = c(4,1,1,1), mfrow = c(1, 3), mar = c(1, 1, 1, 1))

mun <- 'Nosate'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

plot(1:366, Data[1,], type = 'l', col = colors[1], 
     ylim = c(0, max(Data)), xlim = c(0,366),
     xlab = '', ylab ='Mortality',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:366, Data[i,], type = 'l', col = colors[i])
}
points(1:366, Data[10,], type = 'l', col = colors[10], lwd = 2)


library(zoo)
months_length <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
corresp_days <- cumsum(months_length)
days_loc <- c(0, corresp_days[c(2, 4, 8, 10, 12)] - 7)

smoothing_interval <- 14

Data_movmean <- t(rollmean(t(Data), smoothing_interval))
k <- dim(Data_movmean)[2]

plot(1:k, Data_movmean[1,], type = 'l', col = colors[1],
     ylim = c(0, max(Data_movmean)),
     xlab = '', ylab ='',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:k, Data_movmean[i,], type = 'l', col = colors[i])
}
points(1:k, Data_movmean[10,], type = 'l', col = colors[10], lwd = 3)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend_labels, fill = colors, xpd = TRUE, bty = 'n', ncol = 5)




#### Couple plots for Nosate (gray pre COVID) -------------------------------------------

hcl.pals(type = 'qualitative')
colors <- c('gray', 'red3')
legend_labels <- c('2011 - 2019', '2020')

months_length <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
corresp_days <- cumsum(months_length)
days_loc <- c(0, corresp_days[c(2, 4, 8, 10, 12)])

graphics.off()
par(mfrow = c(1, 2))
# par(oma = c(4,1,1,1), mfrow = c(1, 3), mar = c(1, 1, 1, 1))

mun <- 'Nosate'
Data <- CompleteData[CompleteData$name == mun, 8:dim(CompleteData)[2]]

plot(1:366, Data[1,], type = 'l', col = colors[1], 
     ylim = c(0, max(Data)), xlim = c(0,366),
     xlab = '', ylab ='Mortality',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:366, Data[i,], type = 'l', col = colors[1])
}
points(1:366, Data[10,], type = 'l', col = colors[2], lwd = 2)



library(zoo)
months_length <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
corresp_days <- cumsum(months_length)
days_loc <- c(0, corresp_days[c(2, 4, 8, 10, 12)] - 7)

smoothing_interval <- 14

Data_movmean <- t(rollmean(t(Data), smoothing_interval))
k <- dim(Data_movmean)[2]

plot(1:k, Data_movmean[1,], type = 'l', col = colors[1],
     ylim = c(0, max(Data_movmean)),
     xlab = '', ylab ='',
     xaxt = 'n')
axis(1, at = c(0, corresp_days), labels = 0:12)
for(i in 2:9){
  points(1:k, Data_movmean[i,], type = 'l', col = colors[1])
}
points(1:k, Data_movmean[10,], type = 'l', col = colors[2], lwd = 2)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend_labels, fill = colors, xpd = TRUE,
       bty = 'n', ncol = 2)


#### ... ####
#### ... ####
#### ... ####
#### Plot total deaths per year -------------------------------------------

DeathsPerYear <- CompleteData[, c(1, 7)]
DeathsPerYear$deaths <- rowSums(CompleteData[, 8:dim(CompleteData)[2]])

library(dplyr)
average_deaths_per_year <- DeathsPerYear %>% group_by(year) %>% 
                           summarise(deaths = mean(deaths))

head(DeathsPerYear)



n_municipalities <- length(unique(DeathsPerYear$id))  # losing 6 municipalitis?!
samp <- 1:n_municipalities

graphics.off()
par(mfrow = c(1, 2))

plot(NA, NA, xlim = c(2011, 2020), ylim = c(0, 1500),
     xlab = '', ylab = 'Mortality',
     main = 'Total deaths over the years')
for(municipality in DeathsPerYear$id[samp]){
  points(2011:2020, 
         DeathsPerYear$deaths[DeathsPerYear$id == municipality],
         type = 'l', col = 'lightgray')
}
points(2011:2020, average_deaths_per_year$deaths,
       type = 'l', lwd = 3, col = 'red')

ids <- unique(DeathsPerYear$id)
for(id in ids){
  DeathsPerYear$pop[DeathsPerYear$id == CompleteData$id] <-
    CompleteData$pop[DeathsPerYear$id == CompleteData$id]
}
DeathsPerYear <- unique(DeathsPerYear)
DeathsPerYear$norm_deaths <- DeathsPerYear$deaths/DeathsPerYear$pop

average_deaths_per_year_normalized <- DeathsPerYear %>% group_by(year) %>% 
                                      summarise(deaths = mean(norm_deaths))

head(DeathsPerYear)



plot(NA, NA, xlim = c(2011, 2020), ylim = c(0, 0.06),
     xlab = '', ylab = 'Mortality incidence',
     main = 'Fraction of deaths over the years')
for(municipality in DeathsPerYear$id[samp]){
  points(2011:2020, 
         DeathsPerYear$norm_deaths[DeathsPerYear$id == municipality],
         type = 'l', col = 'lightgray')
}
points(2011:2020, average_deaths_per_year_normalized$deaths,
       type = 'l', lwd = 3, col = 'red')


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom', c('Municipalities', 'Average'), 
       fill = c('gray', 'red3'), xpd = TRUE,
       bty = 'n', ncol = 2)



#### ... ####
#### ... ####
#### ... ####
#### Define the t_col function -----------------------------------------------
# source: https://www.dataanalytics.org.uk/make-transparent-colors-in-r/

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}



#### Empirical comparison the distribution of the annual deaths --------------
library(dplyr)

average_deaths_mun_2011_1019 <- DeathsPerYear[DeathsPerYear$year != 2020,] %>% 
  group_by(id) %>% 
  summarise(norm_deaths = mean(norm_deaths))

hcl.pals(type = 'qualitative')
colors <- hcl.colors(10, palette = 'Dark 3', rev = T)
colors[10] <- 'red3'
legend_labels <- as.character(2011:2020)

par(layout(cbind(c(1,6, 11),
                 c(2,7, 11),
                 c(3,8, 11),
                 c(4,9, 11),
                 c(5,10,11))), 
    mai = c(0.4, 0.1, 0.4, 0.1))
for(i in 1:9){
  hist(DeathsPerYear$norm_deaths[DeathsPerYear$year==2010+i], 
       col = colors[i], breaks = 40, 
       ylim = c(0, 1600), xlim = c(0, 0.05),
       main = paste(2010+i),
       yaxt = 'n', xlab = '', ylab = '')
}
hist(DeathsPerYear$norm_deaths[DeathsPerYear$year==2020], 
     col = colors[10], breaks = 80, 
     ylim = c(0, 1600), xlim = c(0, 0.05),
     main = paste(2020),
     yaxt = 'n', xlab = '', ylab = '')
boxplot(DeathsPerYear$norm_deaths[DeathsPerYear$year==2011],
        DeathsPerYear$norm_deaths[DeathsPerYear$year==2012],
        DeathsPerYear$norm_deaths[DeathsPerYear$year==2013],
        DeathsPerYear$norm_deaths[DeathsPerYear$year==2014],
        DeathsPerYear$norm_deaths[DeathsPerYear$year==2015],
        DeathsPerYear$norm_deaths[DeathsPerYear$year==2016],
        DeathsPerYear$norm_deaths[DeathsPerYear$year==2017],
        DeathsPerYear$norm_deaths[DeathsPerYear$year==2018],
        DeathsPerYear$norm_deaths[DeathsPerYear$year==2019],
        average_deaths_mun_2011_1019$norm_deaths,
        DeathsPerYear$norm_deaths[DeathsPerYear$year==2020],
        col = c(colors[1:9], 'lightgray', colors[10]), 
        ylim = c(0, 0.03), bty = 'n', yaxt = 'n',
        names = c(2011:2019, 'average', 2020))
# 6 x 9


colors_comparison <- c(t_col(colors[10], perc = 30), t_col('lightgray', perc = 50))
layout(rbind(c(1,2,2,2), c(1,2,2,2)))
boxplot(DeathsPerYear$norm_deaths[DeathsPerYear$year==2020],
        average_deaths_mun_2011_1019$norm_deaths,
        col = colors_comparison, ylim = c(0, 0.05),
        names = c(2020, 'average'))
hist(DeathsPerYear$norm_deaths[DeathsPerYear$year==2020], 
     col = colors_comparison[1], breaks = 80, 
     ylim = c(0, 1600), xlim = c(0, 0.05),
     main = '', xlab = '', ylab = '', yaxt = 'n')
hist(average_deaths_mun_2011_1019$norm_deaths, add = T,
     col = colors_comparison[2], breaks = 40)
legend('topright', c('2020', '2011-2019 average'),
       fill = colors_comparison, bty = 'n', cex = 1.3)




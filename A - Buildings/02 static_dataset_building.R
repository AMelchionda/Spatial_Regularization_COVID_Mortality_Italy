###########################################################
##################### STATIC DATASET #####################
###########################################################

#### Load the data -----------------------------------------------------------

# main_path <- 'C:/Users/aless/Documents/NAPDE project'
# main_path <- '/Users/Cheesecake/Desktop/NAPDE/Project'
# main_path <- 'C:/Users/Tommaso/Desktop/NAPDE/Project'
main_path <- ...

setwd(paste(main_path, '/Resources', sep = ''))
Data <- read.csv("raw_complete.csv", encoding = "UTF-8")

setwd('../')

#### Compute averages ---------------------------

library(dplyr)

DeathsPerYear <- Data[, c('id', 'pop', 'year')]
DeathsPerYear$deaths <- rowSums(Data[, 8:dim(Data)[2]])  # Add the total
DeathsPerYear$norm_deaths <- DeathsPerYear$deaths/DeathsPerYear$pop

averages <- DeathsPerYear[DeathsPerYear$year != 2020, c('id', 'deaths', 'norm_deaths')] %>%
  group_by(id) %>%
  summarize(mean_deaths = mean(deaths), mean_norm_deaths = mean(norm_deaths))

covid <- DeathsPerYear[DeathsPerYear$year == 2020, c('id', 'deaths', 'norm_deaths')] %>%
  group_by(id) %>%
  summarize(mean_deaths_covid = mean(deaths), mean_norm_deaths_covid = mean(norm_deaths))

#### Merge final information ---------------------------

info <- unique(Data[c('id', 'name', 'x', 'y', 'cod_reg')])
reg <- merge(averages, covid)
reg <- merge(info, reg)

#### Export dataset ---------------------------

setwd(paste(main_path, '/Resources', sep = ''))
write.csv(reg, "static_dataset.csv", row.names = FALSE)

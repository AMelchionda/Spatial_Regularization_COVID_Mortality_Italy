###########################################################
##################### TIME SERIES DATASET #################
###########################################################

#### Load the data -----------------------------------------------------------

# main_path <- 'C:/Users/aless/Documents/NAPDE project'
# main_path <- '/Users/Cheesecake/Desktop/NAPDE/Project'
# main_path <- 'C:/Users/Tommaso/Desktop/NAPDE/Project'
main_path <- ...

setwd(paste(main_path, '/Resources', sep = ''))
Data <- read.csv("raw_complete.csv", encoding = "UTF-8")

setwd('../')

library(dplyr)

DeathsPerYear <- Data[, c('id', 'pop', 'year')]

#### Compute averages by month ---------------------------

start <- 8  # First 7 columns are not death counts
end <- 39  # January initialization

for (i in 1:12) {  # Months
  col_n <- paste('month', i, sep = '')
  if(i==1){  # January
    DeathsPerYear[col_n] <- rowSums(Data[, start:end])  # Add the accumulated deaths
    start <- start + 31
    end <- end + 29
  }
  else if(i==2){  # February
    DeathsPerYear[col_n] <- rowSums(Data[, start:end])  # Add the accumulated deaths
    start <- start + 29
    end <- end + 31
  }
  else if(i%%2==1){  # Odd: 31 days
    DeathsPerYear[col_n] <- rowSums(Data[, start:end])  # Add the accumulated deaths
    start <- start + 31
    end <- end + 30
  }
  else if(i%%2==0){  # Even: 30 days
    DeathsPerYear[col_n] <- rowSums(Data[, start:end])  # Add the accumulated deaths
    start <- start + 30
    end <- end + 31
  }
}

DeathsPerYear[,4:15] <- DeathsPerYear[,4:15] / DeathsPerYear$pop

averages <- DeathsPerYear[DeathsPerYear$year != 2020, c(1, 4:15)] %>%
  group_by(id) %>%
  summarise_all(mean)

covid <- DeathsPerYear[DeathsPerYear$year == 2020, c(1, 4:15)] %>%
  group_by(id) %>%
  summarise_all(mean)

colnames(covid)[2:13] <- paste(colnames(covid)[2:13], 'covid', sep = "_")

# Merge final information 
info <- unique(Data[c('id', 'name', 'x', 'y', 'cod_reg')])
reg <- merge(averages, covid)
reg <- merge(info, reg)

# Export dataset 
setwd(paste(main_path, '/Resources', sep = ''))
write.csv(reg, "monthly_dataset.csv", row.names = FALSE)

#### Compute averages by fortnight ---------------------------
# NOT USED IN THE PROJECT

start <- 8  # First 7 columns are not death counts
end <- 8 + 14

while (end < 373) {  # By 2 weeks
  col_n <- paste('week', start - 7, end - 7, sep = '_')
  DeathsPerYear[col_n] <- rowSums(Data[, start:end])  # Add the accumulated deaths
  start <- start + 14
  end <- end + 14
}
# We lose the last day but it makes little difference

DeathsPerYear[,4:29] <- DeathsPerYear[,4:29] / DeathsPerYear$pop

averages <- DeathsPerYear[DeathsPerYear$year != 2020, c(1, 4:29)] %>%
  group_by(id) %>%
  summarise_all(mean)

covid <- DeathsPerYear[DeathsPerYear$year == 2020, c(1, 4:29)] %>%
  group_by(id) %>%
  summarise_all(mean)

colnames(covid)[2:27] <- paste(colnames(covid)[2:27], 'covid', sep = "_")

# Merge final information 
info <- unique(Data[c('id', 'name', 'x', 'y', 'cod_reg')])
reg <- merge(averages, covid)
reg <- merge(info, reg)

# Export dataset 
setwd(paste(main_path, '/Resources', sep = ''))
write.csv(reg, "biweekly_dataset.csv", row.names = FALSE)


#### Compute averages by week ---------------------------

start <- 8  # First 7 columns are not death counts
end <- 8 + 7

while (end < 373) {  # By 2 weeks
  col_n <- paste('week', start - 7, end - 7, sep = '_')
  DeathsPerYear[col_n] <- rowSums(Data[, start:end])  # Add the accumulated deaths
  start <- start + 7
  end <- end + 7
}
# We lose the last day but it makes little difference

DeathsPerYear[,4:55] <- DeathsPerYear[,4:55] / DeathsPerYear$pop

averages <- DeathsPerYear[DeathsPerYear$year != 2020, c(1, 4:55)] %>%
  group_by(id) %>%
  summarise_all(mean)

covid <- DeathsPerYear[DeathsPerYear$year == 2020, c(1, 4:55)] %>%
  group_by(id) %>%
  summarise_all(mean)

colnames(covid)[2:53] <- paste(colnames(covid)[2:53], 'covid', sep = "_")

# Merge final information 
info <- unique(Data[c('id', 'name', 'x', 'y', 'cod_reg')])
reg <- merge(averages, covid)
reg <- merge(info, reg)

# Export dataset 
setwd(paste(main_path, '/Resources', sep = ''))
write.csv(reg, "weekly_dataset.csv", row.names = FALSE)



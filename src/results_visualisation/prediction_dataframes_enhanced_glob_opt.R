
## Making prediction dataframes

# load libraries 

library(readxl)
library(LarsChill)
library(chillR)
library(dplyr)
library(tidyverse)
library(patchwork)
install.packages("BiocManager")
BiocManager::install("MEIGOR")
source("code/custom_ess.R")
source("code/load_save_fitting_results.R")
source("code/helper_functions.R")

### Location-specific ####

###### Spain Enhanced global opt #####
#read param
param_sum_es_new <- read.csv('new_fitter/ES_fit/param_sum_es.csv')
#read temperature data
temp <- read_excel('MaxMin1978and2020.xlsx')

#make day and month numeric
temp$Month <- as.numeric(temp$Month)
temp$Day <- as.numeric(temp$Day)

#check if the temperature data is complete
check_temperature_record(temp)

#get hourly data from daily data
hourtemps <- stack_hourly_temps(temp, latitude=43.48)

#only take dataframe with data
hourtemps <- hourtemps[[1]]


# Initialize a list to store prediction dataframes for different cultivars
prediction_df_es_new <- list()

# Iterate over all cultivars
cultivar_list <- c("Blanquina", "Clara", "Collaos", "Coloradona", "DelaRiega", 
                   "Perezosa", "Perico","Raxao", "Teorica", "Verdialona", "Xuanina")

for (cultivar in cultivar_list) {
  # Initialize a list to store prediction dataframes for different rounds
  prediction_df_cultivar <- list()
  
  for (i in 1:10) {
    # Load data for the current round
    cal_data <- read.csv(paste0("new_fitter/ES_fit/data_split/cal_", i, "_", cultivar, ".csv"))
    val_data <- read.csv(paste0("new_fitter/ES_fit/data_split/val_", i, "_", cultivar, ".csv"))
    SeasonList <- genSeasonList(hourtemps, mrange = c(8, 6), years = cal_data$Year)
    eval_SeasonList <- genSeasonList(hourtemps, mrange = c(8, 6), years = val_data$Year)
    
    fit_res_formatted <- sprintf("Fit_res_%02d", i)
    
    pheno_cal_pred <- return_predicted_days(
      as.numeric(param_sum_es_new$Value[param_sum_es_new$Fit_res == fit_res_formatted & param_sum_es_new$Cultivar == cultivar]), 
      modelfn = custom_PhenoFlex_GDHwrapper, 
      SeasonList = SeasonList
    )
    
    pheno_val_pred <- return_predicted_days(
      as.numeric(param_sum_es_new$Value[param_sum_es_new$Fit_res == fit_res_formatted & param_sum_es_new$Cultivar == cultivar]), 
      modelfn = custom_PhenoFlex_GDHwrapper, 
      SeasonList = eval_SeasonList
    )
    
    # Create prediction dataframe for the current round and cultivar
    cal_data_df <- data.frame(
      location = "Spain",
      data_type = "Calibration",
      Species = "Apple",
      Cultivar = cultivar,
      round = i,
      Year = cal_data$Year,
      Observed = cal_data$pheno,
      Predicted = pheno_cal_pred
    )
    
    val_data_df <- data.frame(
      location = "Spain",
      data_type = "Validation",
      Species = "Apple",
      Cultivar = cultivar,
      round = i,
      Year = val_data$Year,
      Observed = val_data$pheno,
      Predicted = pheno_val_pred
    )
    
    # Store the dataframes for the current round and cultivar in the list
    prediction_df_cultivar[[i]] <- dplyr::bind_rows(cal_data_df, val_data_df)
  }
  
  # Combine all dataframes from different rounds and cultivars into a single dataframe
  prediction_df_cultivar_combined <- dplyr::bind_rows(prediction_df_cultivar)
  prediction_df_es_new[[cultivar]] <- prediction_df_cultivar_combined
}

# Combine all dataframes from different cultivars into a single dataframe
prediction_df_es_new <- dplyr::bind_rows(prediction_df_es_new)

## save 
write.csv(prediction_df_es_new,"enhanced_global_opt/ES_fit/prediction_df_es.csv", row.names = FALSE)

#  run with the lowest calibration error
performance_df_es_best_calibration <- prediction_df_es_new %>%
  filter(data_type == "Calibration") %>%
  group_by(round) %>%
  summarise(RMSE = RMSEP(Predicted,Observed))

#round 1 had the lowest error

###### Germany Enhanced global opt #####
#read parameters
param_sum_de_new <- read.csv('new_fitter/DE_fit/param_sum_de.csv')

#homogenize with data split name files 
param_sum_de_new <- param_sum_de_new %>%
  mutate(Cultivar = recode( Cultivar,
                            Berlepsch = 'berlepsch',
                            `Cox Orange` = 'cox_orange',
                            `Golden Delicious` = 'golden_delicious',
                            `James Grieve` = 'james_grieve',
                            `Roter Boskoop` = 'roter_boskoop'))

#read temperature data
temp_de <- read.csv('Edu/weather_data_final.csv')

#get hourly data from daily data
hourtemps_de <- stack_hourly_temps(temp_de,latitude = 50.62)

#only take dataframe with data
hourtemps_de <- hourtemps_de[[1]]


# Initialize a list to store prediction dataframes for different cultivars
prediction_df_de_new <- list()

# Iterate over all cultivars
cultivar_list_de <- c('berlepsch', 'cox_orange', 'golden_delicious', 'james_grieve', 'roter_boskoop')

for (cultivar in cultivar_list_de) {
  # Initialize a list to store prediction dataframes for different rounds
  prediction_df_cultivar_de <- list()
  
  for (i in 1:10) {
    # Load data for the current round
    cal_data <- read.csv(paste0("new_fitter/DE_fit/data_split/cal_", i, "_", cultivar, ".csv"))
    val_data <- read.csv(paste0("new_fitter/DE_fit/data_split/val_", i, "_", cultivar, ".csv"))
    SeasonList <- genSeasonList(hourtemps_de, mrange = c(8, 6), years = cal_data$Year)
    #error in generating the seasonlists. somehow only some years have data!
    eval_SeasonList <- genSeasonList(hourtemps_de, mrange = c(8, 6), years = val_data$Year)
    #same error here
    #also mixing up years between calibration and validation
    #
    
    fit_res_formatted <- sprintf("Fit_res_%02d", i)
    
    pheno_cal_pred <- return_predicted_days(
      as.numeric(param_sum_de_new$Value[param_sum_de_new$Fit_res == fit_res_formatted & param_sum_de_new$Cultivar == cultivar]), 
      modelfn = custom_PhenoFlex_GDHwrapper, 
      SeasonList = SeasonList
    )
    
    pheno_val_pred <- return_predicted_days(
      as.numeric(param_sum_de_new$Value[param_sum_de_new$Fit_res == fit_res_formatted & param_sum_de_new$Cultivar == cultivar]), 
      modelfn = custom_PhenoFlex_GDHwrapper, 
      SeasonList = eval_SeasonList
    )
    
    # Create prediction dataframe for the current round and cultivar
    cal_data_df <- data.frame(
      location = "Germany",
      data_type = "Calibration",
      Species = "Apple",
      Cultivar = cultivar,
      round = i,
      Year = cal_data$Year,
      Observed = cal_data$Bloom_full_JDay,
      Predicted = pheno_cal_pred
    )
    
    val_data_df <- data.frame(
      location = "Germany",
      data_type = "Validation",
      Species = "Apple",
      Cultivar = cultivar,
      round = i,
      Year = val_data$Year,
      Observed = val_data$Bloom_full_JDay,
      Predicted = pheno_val_pred
    )
    
    # Store the dataframes for the current round and cultivar in the list
    prediction_df_cultivar_de[[i]] <- dplyr::bind_rows(cal_data_df, val_data_df)
  }
  
  # Combine all dataframes from different rounds and cultivars into a single dataframe
  prediction_df_cultivar_combined_de <- dplyr::bind_rows(prediction_df_cultivar_de)
  prediction_df_de_new[[cultivar]] <- prediction_df_cultivar_combined_de
}

# Combine all dataframes from different cultivars into a single dataframe
prediction_df_de_new <- dplyr::bind_rows(prediction_df_de_new)

prediction_df_de_new <- prediction_df_de_new %>%
  mutate(Cultivar = recode( Cultivar,
                            berlepsch = 'Berlepsch',
                            cox_orange = 'Cox Orange',
                            golden_delicious = 'Golden Delicious',
                            james_grieve = 'James Grieve',
                            roter_boskoop = 'Roter Boskoop'))
## save 
write.csv(prediction_df_de_new,"enhanced_global_opt/DE_fit/prediction_df_de.csv", row.names = FALSE)


#run with the lowest calibration error
performance_df_de_best_calibration <- prediction_df_de_new %>%
  filter(data_type == "Calibration") %>%
  group_by(round) %>%
  summarise(RMSE = RMSEP(Predicted,Observed))


# round 1 has the lowest error 


### Species-specific ####


#read in parameters 
param_sum_apple_new <- read.csv('new_fitter/de_es_fit/param_sum_de_es.csv')


###### Spain  #####
# Initialize a list to store prediction dataframes for different cultivars
prediction_df_apple_new_es <- list()

# Iterate over all cultivars
cultivar_list_apple_es <- c("Blanquina","Clara", "Collaos", "Coloradona", "DelaRiega", "Perezosa",
                            "Perico","Raxao", "Teorica", "Verdialona", "Xuanina")

for (cultivar in cultivar_list_apple_es) {
  # Initialize a list to store prediction dataframes for different rounds
  prediction_df_cultivar <- list()
  
  for (i in 1:10) {
    # add a step for a dataframe
    
    # Load data for the current round
    cal_data <- read.csv(paste0("new_fitter/de_es_fit/data_split/cal_", i, "_", cultivar, ".csv"))
    val_data <- read.csv(paste0("new_fitter/de_es_fit/data_split/val_", i, "_", cultivar, ".csv"))
    
    
    SeasonList <- genSeasonList(hourtemps, mrange = c(8, 6), years = cal_data$Year)
    eval_SeasonList <- genSeasonList(hourtemps, mrange = c(8, 6), years = val_data$Year)
    
    fit_res_formatted <- sprintf("Fit_res_%02d", i)
    
    pheno_cal_pred <- return_predicted_days(
      as.numeric(param_sum_apple_new$Value[param_sum_apple_new$Fit_res == fit_res_formatted & param_sum_apple_new$Cultivar == cultivar]), 
      modelfn = custom_PhenoFlex_GDHwrapper, 
      SeasonList = SeasonList
    )
    
    pheno_val_pred <- return_predicted_days(
      as.numeric(param_sum_apple_new$Value[param_sum_apple_new$Fit_res == fit_res_formatted & param_sum_apple_new$Cultivar == cultivar]), 
      modelfn = custom_PhenoFlex_GDHwrapper, 
      SeasonList = eval_SeasonList
    )
    
    # Create prediction dataframe for the current round and cultivar
    cal_data_df <- data.frame(
      location = "Spain",
      data_type = "Calibration",
      Species = "Apple",
      Cultivar = cultivar,
      round = i,
      Year = cal_data$Year,
      Observed = cal_data$pheno,
      Predicted = pheno_cal_pred
    )
    
    val_data_df <- data.frame(
      location = "Spain",
      data_type = "Validation",
      Species = "Apple",
      Cultivar = cultivar,
      round = i,
      Year = val_data$Year,
      Observed = val_data$pheno,
      Predicted = pheno_val_pred
    )
    
    # Store the dataframes for the current round and cultivar in the list
    prediction_df_cultivar[[i]] <- dplyr::bind_rows(cal_data_df, val_data_df)
  }
  
  # Combine all dataframes from different rounds and cultivars into a single dataframe
  prediction_df_cultivar_combined <- dplyr::bind_rows(prediction_df_cultivar)
  prediction_df_apple_new_es[[cultivar]] <- prediction_df_cultivar_combined
}

# Combine all dataframes from different cultivars into a single dataframe
prediction_df_apple_new_es <- dplyr::bind_rows(prediction_df_apple_new_es)

###### Germany #####
#homogenize with data split name files 
param_sum_apple_new <- param_sum_apple_new %>%
  mutate(Cultivar = recode( Cultivar,
                            Berlepsch = 'berlepsch',
                            `Cox Orange` = 'cox_orange',
                            `Golden Delicious` = 'golden_delicious',
                            `James Grieve` = 'james_grieve',
                            `Roter Boskoop` = 'roter_boskoop'))


# Initialize a list to store prediction dataframes for different cultivars
prediction_df_apple_new_de <- list()

cultivar_list_apple_de <- c("berlepsch","cox_orange","golden_delicious","james_grieve","roter_boskoop")

for (cultivar in cultivar_list_apple_de) {
  # Initialize a list to store prediction dataframes for different rounds
  prediction_df_cultivar <- list()
  
  for (i in 1:10) {
    
    # Load data for the current round
    cal_data <- read.csv(paste0("new_fitter/de_es_fit/data_split/cal_", i, "_", cultivar, ".csv"))
    val_data <- read.csv(paste0("new_fitter/de_es_fit/data_split/val_", i, "_", cultivar, ".csv"))
    
    
    SeasonList <- genSeasonList(hourtemps_de, mrange = c(8, 6), years = cal_data$Year)
    eval_SeasonList <- genSeasonList(hourtemps_de, mrange = c(8, 6), years = val_data$Year)
    
    fit_res_formatted <- sprintf("Fit_res_%02d", i)
    
    pheno_cal_pred <- return_predicted_days(
      as.numeric(param_sum_apple_new$Value[param_sum_apple_new$Fit_res == fit_res_formatted & param_sum_apple_new$Cultivar == cultivar]), 
      modelfn = custom_PhenoFlex_GDHwrapper, 
      SeasonList = SeasonList
    )
    
    pheno_val_pred <- return_predicted_days(
      as.numeric(param_sum_apple_new$Value[param_sum_apple_new$Fit_res == fit_res_formatted & param_sum_apple_new$Cultivar == cultivar]), 
      modelfn = custom_PhenoFlex_GDHwrapper, 
      SeasonList = eval_SeasonList
    )
    
    # Create prediction dataframe for the current round and cultivar
    cal_data_df <- data.frame(
      location = "Germany",
      data_type = "Calibration",
      Species = "Apple",
      Cultivar = cultivar,
      round = i,
      Year = cal_data$Year,
      Observed = cal_data$Bloom_full_JDay,
      Predicted = pheno_cal_pred
    )
    
    val_data_df <- data.frame(
      location = "Germany",
      data_type = "Validation",
      Species = "Apple",
      Cultivar = cultivar,
      round = i,
      Year = val_data$Year,
      Observed = val_data$Bloom_full_JDay,
      Predicted = pheno_val_pred
    )
    
    # Store the dataframes for the current round and cultivar in the list
    prediction_df_cultivar[[i]] <- dplyr::bind_rows(cal_data_df, val_data_df)
  }
  
  # Combine all dataframes from different rounds and cultivars into a single dataframe
  prediction_df_cultivar_combined <- dplyr::bind_rows(prediction_df_cultivar)
  prediction_df_apple_new_de[[cultivar]] <- prediction_df_cultivar_combined
}

# Combine all dataframes from different cultivars into a single dataframe
prediction_df_apple_new_de <- dplyr::bind_rows(prediction_df_apple_new_de)

###### Spain and Gemrany together#####
#combine all together 
prediction_df_apple_new <- dplyr::bind_rows(prediction_df_apple_new_es,prediction_df_apple_new_de)

#change names for German cultivars
prediction_df_apple_new <- prediction_df_apple_new %>%
  mutate(Cultivar = recode( Cultivar,
                            berlepsch = 'Berlepsch',
                            cox_orange = 'Cox Orange',
                            golden_delicious = 'Golden Delicious',
                            james_grieve = 'James Grieve',
                            roter_boskoop = 'Roter Boskoop'))
## save 
write.csv(prediction_df_apple_new,"enhanced_global_opt/de_es_fit/prediction_df_apple.csv", row.names = FALSE)

# run with the lowest calibration error
performance_df_apple_new_best_calibration <- prediction_df_apple_new %>%
  filter(data_type == "Calibration") %>%
  group_by(round) %>%
  summarise(RMSE = RMSEP(Predicted,Observed))

# round 3  had the lowest error 
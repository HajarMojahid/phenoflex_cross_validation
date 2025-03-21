# Import the required libraries
library(chillR)
library(readxl)
library(ggplot2)
library(tidyverse)

# Source the helpers, including the multiple PhenoFlex fitter
source('code/99_helper_functions.R')

# Read temperature data
temps_ES <- as.data.frame(read_excel("data/Hajars_data/weather/MaxMin1978and2020.xlsx"))

# Make day and month numeric
temps_ES$Month <- as.numeric(temps_ES$Month)
temps_ES$Day <- as.numeric(temps_ES$Day)

# Check if the temperature data is complete
check_temperature_record(temps_ES)

# Get hourly data from daily data
hourtemps_ES <- stack_hourly_temps(temps_ES, latitude = 43.48)

# Only take dataframe with data
hourtemps_ES <- hourtemps_ES[[1]]

# Read phenology data and prepare it for fitting
cal_data_ES <- load_temperature_scenarios("data/Hajars_data/calibration/", prefix = "")

# Validation data
val_data_ES <- load_temperature_scenarios("data/Hajars_data/validation/", prefix = "")

# Merge all pheno
pheno_all <- bind_rows(bind_rows(cal_data_ES, .id = "Cultivar"), bind_rows(val_data_ES, .id = "Cultivar"))

# Select only relevant columns
pheno_all <- select(pheno_all, Cultivar, Year, pheno)

# order the DF
pheno_all <- pheno_all[order(pheno_all$Cultivar, pheno_all$Year), ]

# Generate the list of calibration and validation seasons
# First define the share of calibration and validation splits
percentage <- 0.8

# Define the name of the cultivars (includes the data thing)
cvs <- unique(pheno_all$Cultivar)

# Placeholders to save the results
calibration_list <- list()
validation_list <- list()

# Sample the order of cultivars
index_vector <- sample(1 : length(unique(pheno_all$Cultivar)))

# Loop to extract specific subsets of DFs
for (i in 1 : length(index_vector)){
  
  # Generate a temporary df according to cultivar
  temp_df <- pheno_all[pheno_all$Cultivar == cvs[index_vector[[i]]], ]
  
  # Sample the index for calibration
  index <- sample(nrow(temp_df), nrow(temp_df) * percentage)
  
  # Subset the rows for calibration and validation
  calibration_list[[i]] <- temp_df[sort(index), ]
  validation_list[[i]] <- temp_df[-sort(index), ]
  
  # Name the elements of the list according to cultivar
  names(calibration_list)[i] <- unique(temp_df$Cultivar)
  names(validation_list)[i] <- unique(temp_df$Cultivar)
}


# Split hourly temperature data to the respective seasons (start in august, end in july of following year)
SeasonList <- list()
eval_SeasonList <- list()

for (i in 1 : length(calibration_list)) {
  
  # Extract the years for calibration and validation
  years_cal <- calibration_list[[i]]$Year
  years_val <- validation_list[[i]]$Year
  
  # Generate the seasons according to the years
  SeasonList[[i]] <- genSeasonList(hourtemps_ES, mrange = c(8, 6), years = years_cal)
  eval_SeasonList[[i]] <- genSeasonList(hourtemps_ES, mrange = c(8, 6), years = years_val)
}

# Fix the names in both lists
names(SeasonList) <- names(calibration_list)
names(eval_SeasonList) <- names(validation_list)



###### Run 1 #######
# You need to specify how many cultivars you want to fit jointly
n_cult <- length(calibration_list)

# Combine phenology observations used for model training
bloomJDays_list <- lapply(calibration_list, function (x) x[["pheno"]])

# Fix the name
names(bloomJDays_list) <- names(calibration_list)

# The number of parameters fitted increase, but overall they decrease because only yc, zc and s1 are cultivar specific
lower_01 <- c(rep(10, n_cult),  rep(30, n_cult), rep(0.05, n_cult),  0, 2300.0,  8900.0, 5300.0,      4.9e13,  0,  0,  0,  0.05)
par_01 <-   c(rep(45, n_cult), rep(220, n_cult),  rep(0.5, n_cult), 25, 3372.8,  9900.3, 6319.5, 5.939917e13,  4, 36,  4,  1.60)
upper_01 <- c(rep(90, n_cult), rep(500, n_cult),  rep(1.0, n_cult), 35, 4300.0, 10900.0, 7300.0,      6.9e13, 15, 45, 15, 50.00)


# Please note that you need to specify the number of cultivars in the function call, otherwise nothing changed
Fit_res_01 <- phenologyFitter_combined_fitting(par.guess = par_01, 
                                               modelfn = PhenoFlex_GDHwrapper,
                                               bloomJDays = bloomJDays_list,
                                               SeasonList = SeasonList,
                                               n_cult = n_cult,
                                               lower = lower_01,
                                               upper = upper_01,
                                               seed = 001,
                                               control = list(smooth = FALSE, verbose = FALSE, maxit = 500,
                                                              nb.stop.improvement = 10))

# The next function returns the bloom dates as a list, one element for each cultivar.
# Each element is a vector with the bloomdates
pbloomJDays_01 <- get_bloom_days(Fit_res_01)
Fit_res_01$par

# Check the calibration error for each cultivar
# First define the lists that will contain the estimations
rmse_list_01 <- list()
rpiq_list_01 <- list()
mean_list_01 <- list()
mean_abs_list_01 <- list ()

# Implement the loop
for (i in 1 : n_cult) {
  
  # Compute the root mean square error
  rmse_list_01[[i]] <- RMSEP(pbloomJDays_01[[i]], bloomJDays_list[[i]], na.rm = TRUE)
  
  # Compute the ratio of performance to interquartile distance
  rpiq_list_01[[i]] <- RPIQ(pbloomJDays_01[[i]], bloomJDays_list[[i]], na.rm = TRUE)
  
  # Compute the mean absolute error (no negative values)
  mean_abs_list_01[[i]] <- mean(abs(bloomJDays_list[[i]]- pbloomJDays_01[[i]]), na.rm = TRUE)
  
  # Compute the mean error (with negative values)
  mean_list_01[[i]] <- mean(bloomJDays_list[[i]]- pbloomJDays_01[[i]], na.rm = TRUE)
  
}

# Fix the names of the lists to know the cultivars
names(rmse_list_01) <- names(SeasonList)
names(rpiq_list_01) <- names(SeasonList)
names(mean_abs_list_01) <- names(SeasonList)
names(mean_list_01) <- names(SeasonList)

# Check the temp response curve 
temp_values <- seq(-5, 55, 0.1)

# Implement the plot. Note that all cultivars are forced to show the same response to chill
# and heat accumulation
temp_plot_01 <- get_temp_response_plot(Fit_res_01$par[c(1, n_cult + 1, (n_cult * 2) + 1,
                                                        ((n_cult * 3) + 1) : length(par_01))],
                                       temp_values)
temp_plot_01

# Check initial RMSE for calibration
unlist(rmse_list_01)

# Median and mean RMSE for calibration
median(unlist(rmse_list_01))
mean(unlist(rmse_list_01))



###### Run 2 #######
# Define the new sets of parameters for the next fitting
# Need to take a look at the previously fitted parameters to define new lower and upper bounds
# yc
c(min(Fit_res_01$par[1 : n_cult]), max(Fit_res_01$par[1 : n_cult]))
(yc_01 <- Fit_res_01$par[1 : n_cult])
  
# zc
c(min(Fit_res_01$par[(n_cult + 1) : (n_cult * 2)]), max(Fit_res_01$par[(n_cult + 1) : (n_cult * 2)]))
(zc_01 <- Fit_res_01$par[(n_cult + 1) : (n_cult * 2)])

# s1
c(min(Fit_res_01$par[((n_cult * 2) + 1) : (n_cult * 3)]), max(Fit_res_01$par[((n_cult * 2) + 1) : (n_cult * 3)]))
(s1_01 <- Fit_res_01$par[((n_cult * 2) + 1) : (n_cult * 3)])

# Remaining pars
Fit_res_01$par[((n_cult * 3) + 1) : length(par_01)]

# Adjust the values
lower_02 <- c(rep(11, n_cult), rep(110, n_cult), rep(0.01, n_cult), 20, 2300.0,  8900.0, 5300.0,      4.9e13, -1, 30,  0,  0.05)
par_02 <-   Fit_res_01$par
upper_02 <- c(rep(74, n_cult), rep(580, n_cult), rep(2.50, n_cult), 36, 4300.0, 10900.0, 7300.0,      6.9e13, 15, 46, 15, 62.00)


# Please note that you need to specify the number of cultivars in the function call, otherwise nothing changed
Fit_res_02 <- phenologyFitter_combined_fitting(par.guess = par_02, 
                                               modelfn = PhenoFlex_GDHwrapper,
                                               bloomJDays = bloomJDays_list,
                                               SeasonList = SeasonList,
                                               n_cult = n_cult,
                                               lower = lower_02,
                                               upper = upper_02,
                                               seed = 001,
                                               control = list(smooth = FALSE, verbose = FALSE, maxit = 500,
                                                              nb.stop.improvement = 10))

# The next function returns the bloom dates as a list, one element for each cultivar.
# Each element is a vector with the bloomdates
pbloomJDays_02 <- get_bloom_days(Fit_res_02)
Fit_res_02$par

# Check the calibration error for each cultivar
# First define the lists that will contain the estimations
rmse_list_02 <- list()
rpiq_list_02 <- list()
mean_list_02 <- list()
mean_abs_list_02 <- list ()

# Implement the loop
for (i in 1 : n_cult) {
  
  # Compute the root mean square error
  rmse_list_02[[i]] <- RMSEP(pbloomJDays_02[[i]], bloomJDays_list[[i]], na.rm = TRUE)
  
  # Compute the ratio of performance to interquartile distance
  rpiq_list_02[[i]] <- RPIQ(pbloomJDays_02[[i]], bloomJDays_list[[i]], na.rm = TRUE)
  
  # Compute the mean absolute error (no negative values)
  mean_abs_list_02[[i]] <- mean(abs(bloomJDays_list[[i]]- pbloomJDays_02[[i]]), na.rm = TRUE)
  
  # Compute the mean error (with negative values)
  mean_list_02[[i]] <- mean(bloomJDays_list[[i]]- pbloomJDays_02[[i]], na.rm = TRUE)
  
}

# Fix the names of the lists to know the cultivars
names(rmse_list_02) <- names(SeasonList)
names(rpiq_list_02) <- names(SeasonList)
names(mean_abs_list_02) <- names(SeasonList)
names(mean_list_02) <- names(SeasonList)

# Implement the plot. Note that all cultivars are forced to show the same response to chill
# and heat accumulation
temp_plot_02 <- get_temp_response_plot(Fit_res_02$par[c(1, n_cult + 1, (n_cult * 2) + 1,
                                                        ((n_cult * 3) + 1) : length(par_02))],
                                       temp_values)
temp_plot_02

# Check if RMSE for calibration increased or decreased
unlist(rmse_list_02) - unlist(rmse_list_01)

# Median and mean RMSE for calibration
median(unlist(rmse_list_02))
mean(unlist(rmse_list_02))



###### Run 3 #######
# Define the new sets of parameters for the next fitting
# Need to take a look at the previously fitted parameters to define new lower and upper bounds
# yc
c(min(Fit_res_02$par[1 : n_cult]), max(Fit_res_02$par[1 : n_cult]))
(yc_02 <- Fit_res_02$par[1 : n_cult])

# zc
c(min(Fit_res_02$par[(n_cult + 1) : (n_cult * 2)]), max(Fit_res_02$par[(n_cult + 1) : (n_cult * 2)]))
(zc_02 <- Fit_res_02$par[(n_cult + 1) : (n_cult * 2)])

# s1
c(min(Fit_res_02$par[((n_cult * 2) + 1) : (n_cult * 3)]), max(Fit_res_02$par[((n_cult * 2) + 1) : (n_cult * 3)]))
(s1_02 <- Fit_res_02$par[((n_cult * 2) + 1) : (n_cult * 3)])

# Remaining pars
Fit_res_02$par[((n_cult * 3) + 1) : length(par_02)]

# Adjust the values
lower_03 <- c(rep(11, n_cult), rep(110, n_cult), rep(0.01, n_cult), 20, 2300.0,  8900.0, 5300.0,      4.9e13, 1, 25,  0,  0.05)
par_03 <-   Fit_res_02$par
upper_03 <- c(rep(74, n_cult), rep(580, n_cult), rep(3.5, n_cult), 33, 4300.0, 10900.0, 7300.0,      6.9e13, 12, 40, 15, 61.00)

# Please note that you need to specify the number of cultivars in the function call, otherwise nothing changed
Fit_res_03 <- phenologyFitter_combined_fitting(par.guess = par_03, 
                                               modelfn = PhenoFlex_GDHwrapper,
                                               bloomJDays = bloomJDays_list,
                                               SeasonList = SeasonList,
                                               n_cult = n_cult,
                                               lower = lower_03,
                                               upper = upper_03,
                                               seed = 001,
                                               control = list(smooth = FALSE, verbose = FALSE, maxit = 1000,
                                                              nb.stop.improvement = 25))

# The next function returns the bloom dates as a list, one element for each cultivar.
# Each element is a vector with the bloomdates
pbloomJDays_03 <- get_bloom_days(Fit_res_03)
Fit_res_03$par

# Check the calibration error for each cultivar
# First define the lists that will contain the estimations
rmse_list_03 <- list()
rpiq_list_03 <- list()
mean_list_03 <- list()
mean_abs_list_03 <- list ()

# Implement the loop
for (i in 1 : n_cult) {
  
  # Compute the root mean square error
  rmse_list_03[[i]] <- RMSEP(pbloomJDays_03[[i]], bloomJDays_list[[i]], na.rm = TRUE)
  
  # Compute the ratio of performance to interquartile distance
  rpiq_list_03[[i]] <- RPIQ(pbloomJDays_03[[i]], bloomJDays_list[[i]], na.rm = TRUE)
  
  # Compute the mean absolute error (no negative values)
  mean_abs_list_03[[i]] <- mean(abs(bloomJDays_list[[i]]- pbloomJDays_03[[i]]), na.rm = TRUE)
  
  # Compute the mean error (with negative values)
  mean_list_03[[i]] <- mean(bloomJDays_list[[i]]- pbloomJDays_03[[i]], na.rm = TRUE)
  
}

# Fix the names of the lists to know the cultivars
names(rmse_list_03) <- names(SeasonList)
names(rpiq_list_03) <- names(SeasonList)
names(mean_abs_list_03) <- names(SeasonList)
names(mean_list_03) <- names(SeasonList)

# Implement the plot. Note that all cultivars are forced to show the same response to chill
# and heat accumulation
temp_plot_03 <- get_temp_response_plot(Fit_res_03$par[c(1, n_cult + 1, (n_cult * 2) + 1,
                                                        ((n_cult * 3) + 1) : length(par_03))],
                                       temp_values)
temp_plot_03

# Check if RMSE for calibration increased or decreased
unlist(rmse_list_03) - unlist(rmse_list_02)

# Median and mean RMSE for calibration
median(unlist(rmse_list_03))
mean(unlist(rmse_list_03))



###### Run 4 #######
# Define the new sets of parameters for the next fitting
# Need to take a look at the previously fitted parameters to define new lower and upper bounds
# yc
c(min(Fit_res_03$par[1 : n_cult]), max(Fit_res_03$par[1 : n_cult]))
(yc_03 <- Fit_res_03$par[1 : n_cult])

# zc
c(min(Fit_res_03$par[(n_cult + 1) : (n_cult * 2)]), max(Fit_res_03$par[(n_cult + 1) : (n_cult * 2)]))
(zc_03 <- Fit_res_03$par[(n_cult + 1) : (n_cult * 2)])

# s1
c(min(Fit_res_03$par[((n_cult * 2) + 1) : (n_cult * 3)]), max(Fit_res_03$par[((n_cult * 2) + 1) : (n_cult * 3)]))
(s1_03 <- Fit_res_03$par[((n_cult * 2) + 1) : (n_cult * 3)])

# Remaining pars
Fit_res_03$par[((n_cult * 3) + 1) : length(par_03)]

# Adjust the values
lower_04 <- c(rep(10, n_cult), rep(50, n_cult), rep(0.001, n_cult), 20, 2700.0,  9100.0, 5500.0, 5.1e13,  0, 21,  0,  0.05)
par_04 <-   Fit_res_03$par
upper_04 <- c(rep(75, n_cult), rep(650, n_cult), rep(4.500, n_cult), 30, 4100.0, 10700.0, 7100.0, 6.7e13, 13, 32, 12, 51.00)


# Please note that you need to specify the number of cultivars in the function call, otherwise nothing changed
Fit_res_04 <- phenologyFitter_combined_fitting(par.guess = par_04, 
                                               modelfn = PhenoFlex_GDHwrapper,
                                               bloomJDays = bloomJDays_list,
                                               SeasonList = SeasonList,
                                               n_cult = n_cult,
                                               lower = lower_04,
                                               upper = upper_04,
                                               seed = 001,
                                               control = list(smooth = FALSE, verbose = FALSE, maxit = 1500,
                                                              nb.stop.improvement = 30))

# The next function returns the bloom dates as a list, one element for each cultivar.
# Each element is a vector with the bloomdates
pbloomJDays_04 <- get_bloom_days(Fit_res_04)
Fit_res_04$par

# Check the calibration error for each cultivar
# First define the lists that will contain the estimations
rmse_list_04 <- list()
rpiq_list_04 <- list()
mean_list_04 <- list()
mean_abs_list_04 <- list ()

# Implement the loop
for (i in 1 : n_cult) {
  
  # Compute the root mean square error
  rmse_list_04[[i]] <- RMSEP(pbloomJDays_04[[i]], bloomJDays_list[[i]], na.rm = TRUE)
  
  # Compute the ratio of performance to interquartile distance
  rpiq_list_04[[i]] <- RPIQ(pbloomJDays_04[[i]], bloomJDays_list[[i]], na.rm = TRUE)
  
  # Compute the mean absolute error (no negative values)
  mean_abs_list_04[[i]] <- mean(abs(bloomJDays_list[[i]]- pbloomJDays_04[[i]]), na.rm = TRUE)
  
  # Compute the mean error (with negative values)
  mean_list_04[[i]] <- mean(bloomJDays_list[[i]]- pbloomJDays_04[[i]], na.rm = TRUE)
  
}

# Fix the names of the lists to know the cultivars
names(rmse_list_04) <- names(SeasonList)
names(rpiq_list_04) <- names(SeasonList)
names(mean_abs_list_04) <- names(SeasonList)
names(mean_list_04) <- names(SeasonList)

# Implement the plot. Note that all cultivars are forced to show the same response to chill
# and heat accumulation
temp_plot_04 <- get_temp_response_plot(Fit_res_04$par[c(1, n_cult + 1, (n_cult * 2) + 1,
                                                        ((n_cult * 3) + 1) : length(par_04))],
                                       temp_values)
temp_plot_04

# Check if RMSE for calibration increased or decreased
unlist(rmse_list_04) - unlist(rmse_list_03)

# Median and mean RMSE for calibration
median(unlist(rmse_list_04))
mean(unlist(rmse_list_04))



###### Run 5 #######
# Define the new sets of parameters for the next fitting
# Need to take a look at the previously fitted parameters to define new lower and upper bounds
# yc
c(min(Fit_res_04$par[1 : n_cult]), max(Fit_res_04$par[1 : n_cult]))
(yc_04 <- Fit_res_04$par[1 : n_cult])

# zc
c(min(Fit_res_04$par[(n_cult + 1) : (n_cult * 2)]), max(Fit_res_04$par[(n_cult + 1) : (n_cult * 2)]))
(zc_04 <- Fit_res_04$par[(n_cult + 1) : (n_cult * 2)])

# s1
c(min(Fit_res_04$par[((n_cult * 2) + 1) : (n_cult * 3)]), max(Fit_res_04$par[((n_cult * 2) + 1) : (n_cult * 3)]))
(s1_04 <- Fit_res_04$par[((n_cult * 2) + 1) : (n_cult * 3)])

# Remaining pars
Fit_res_04$par[((n_cult * 3) + 1) : length(par_04)]

# Adjust the values
lower_05 <- c(rep(35, n_cult), rep(250, n_cult), rep(0.001, n_cult), 22, 3000.0,  9600.0, 6000.0, 5.6e13,  0, 22,  0,  0.05)
par_05 <-   Fit_res_04$par
upper_05 <- c(rep(46, n_cult), rep(460, n_cult), rep(6.500, n_cult), 30, 3600.0, 10200.0, 6600.0, 6.2e13, 10, 32, 10, 41.00)

# Please note that you need to specify the number of cultivars in the function call, otherwise nothing changed
Fit_res_05 <- phenologyFitter_combined_fitting(par.guess = par_05, 
                                               modelfn = PhenoFlex_GDHwrapper,
                                               bloomJDays = bloomJDays_list,
                                               SeasonList = SeasonList,
                                               n_cult = n_cult,
                                               lower = lower_05,
                                               upper = upper_05,
                                               seed = 001,
                                               control = list(smooth = FALSE, verbose = FALSE, maxit = 2000,
                                                              nb.stop.improvement = 40))

# The next function returns the bloom dates as a list, one element for each cultivar.
# Each element is a vector with the bloomdates
pbloomJDays_05 <- get_bloom_days(Fit_res_05)
Fit_res_05$par

# Check the calibration error for each cultivar
# First define the lists that will contain the estimations
rmse_list_05 <- list()
rpiq_list_05 <- list()
mean_list_05 <- list()
mean_abs_list_05 <- list ()

# Implement the loop
for (i in 1 : n_cult) {
  
  # Compute the root mean square error
  rmse_list_05[[i]] <- RMSEP(pbloomJDays_05[[i]], bloomJDays_list[[i]], na.rm = TRUE)
  
  # Compute the ratio of performance to interquartile distance
  rpiq_list_05[[i]] <- RPIQ(pbloomJDays_05[[i]], bloomJDays_list[[i]], na.rm = TRUE)
  
  # Compute the mean absolute error (no negative values)
  mean_abs_list_05[[i]] <- mean(abs(bloomJDays_list[[i]]- pbloomJDays_05[[i]]), na.rm = TRUE)
  
  # Compute the mean error (with negative values)
  mean_list_05[[i]] <- mean(bloomJDays_list[[i]]- pbloomJDays_05[[i]], na.rm = TRUE)
  
}

# Fix the names of the lists to know the cultivars
names(rmse_list_05) <- names(SeasonList)
names(rpiq_list_05) <- names(SeasonList)
names(mean_abs_list_05) <- names(SeasonList)
names(mean_list_05) <- names(SeasonList)

# Implement the plot. Note that all cultivars are forced to show the same response to chill
# and heat accumulation
temp_plot_05 <- get_temp_response_plot(Fit_res_05$par[c(1, n_cult + 1, (n_cult * 2) + 1,
                                                        ((n_cult * 3) + 1) : length(par_05))],
                                       temp_values)
temp_plot_05

# Check if RMSE for calibration increased or decreased
unlist(rmse_list_05) - unlist(rmse_list_04)

# Median and mean RMSE for calibration
median(unlist(rmse_list_05))
mean(unlist(rmse_list_05))


###### Extract the summary of calibration ######
# Use bind_rows to merge the calibration_list
summary_calibration <- bind_rows(calibration_list)

# Add the predicted JDay for bloom using the fifth run
summary_calibration[["Predicted"]] <- unlist(pbloomJDays_05)

# Compute the error
summary_calibration[["Error"]] <- summary_calibration$pheno - summary_calibration$Predicted



###### Validation ######
# Use the estimated parameters to compute bloom dates in the validation datasets. Validation datasets are located in
# validation_list. In addition, eval_SeasonList contains all weather data for the seasons used in the validation
# datasets.

# Define a list to allocate the outputs of the validation
bloom_list_val <- list()
  
# Implement a loop to estimate the bloom dates in the validation datasets
for(i in 1 : length(Fit_res_05$SeasonList)){
  
  # First define the parameters for each cultivar
  par <- Fit_res_05$par[c(i, n_cult + i, (n_cult * 2) + i, ((n_cult * 3) + 1) : length(Fit_res_05$par))]
    
  # Calculate the predicted flower dates
  bloom_list_val[[i]] <- unlist(lapply(eval_SeasonList[[i]], PhenoFlex_GDHwrapper, par = par))
}
  
# Merge the lists of observed and predicted values for the validation seasons
summary_validation <- data.frame(bind_rows(validation_list),
                                 Predicted = unlist(bloom_list_val))

# Compute the error for further metrics
summary_validation[["Error"]] <- summary_validation$pheno - summary_validation$Predicted

# Summarize the data frame to compute RMSEP, RPIQ and other metrics
summary_metrics <- summary_validation %>% group_by(Cultivar) %>% 
  summarize(RMSEP_val = RMSEP(Predicted, pheno),
            RPIQ_val = RPIQ(Predicted, pheno),
            Error_mean_val = mean(Error),
            Error_median_val = median(Error),
            Error_abs_val = mean(abs(Error)),
            IQR_val = IQR(Error))

# Add the calibration metrics
summary_metrics[["RMSEP_cal"]] <- unlist(rmse_list_05)

# Add the calibration metrics
summary_metrics[["RPIQ_cal"]] <- unlist(rpiq_list_05)



###### Collect the parameters from all runs #######
# Create one intermediary list to save the outputs from run "j"
pars_sum <- list()

# Create a list to save the results of the double loop
par_sums <- list()

# Implement the loop over Fit_res and then over cultvivars
for (j in 1 : 5){
  for (i in 1 : n_cult){
    
    # Name the parameters
    par <- data.frame(Par = c("yc", "zc", "s1", "Tu", "E0", "E1", "A0", "A1", "Tf", "Tc", "Tb", "s"))
    
    # First define the parameters for each cultivar
    par[["Value"]] <- eval(sym(paste0("Fit_res_0", j)))$par[c(i, i + n_cult, (n_cult * 2) + i, ((n_cult * 3) + 1) :
                                                                     length(eval(sym(paste0("Fit_res_0", j)))$par))]
    
    # Calculate the predicted flower dates
    pars_sum[[i]] <- par
  }
  
  # Name the elements of the intermediary list "j"
  names(pars_sum) <- names(validation_list)
  
  # Allocate the results of the loop "j" into the final list
  par_sums[[j]] <- pars_sum}

# Name the elements of the final list according to the Fit_res element
names(par_sums) <- paste0("Fit_res_0", 1 : 5)

# Summarize the evolution of parameters
summary_params <- bind_rows(lapply(par_sums, bind_rows, .id = "Cultivar"), .id = "Fit_res")


# Save all outputs to folder as csv
write.csv(summary_calibration, "results/02_summary_cal_ES.csv", row.names = FALSE)
write.csv(summary_validation, "results/05_summary_val_ES.csv", row.names = FALSE)
write.csv(summary_params, "results/08_summary_pars_ES.csv", row.names = FALSE)
write.csv(summary_metrics, "results/11_summary_metrics_ES.csv", row.names = FALSE)



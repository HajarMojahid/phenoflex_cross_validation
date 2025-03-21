# Import the required libraries
library(chillR)
library(readxl)
library(ggplot2)
library(tidyverse)

# Source the helpers, including the multiple PhenoFlex fitter
source('code/99_helper_functions.R')

# Read temperature data
temps_ES <- as.data.frame(read_excel("data/Hajars_data/weather/MaxMin1978and2020.xlsx"))
temps_DE <- read.csv('data/weather_data_final.csv')

# Make day and month numeric
temps_ES$Month <- as.numeric(temps_ES$Month)
temps_ES$Day <- as.numeric(temps_ES$Day)

# Check if the temperature data is complete
check_temperature_record(temps_ES)

# Get hourly data from daily data
hourtemps_ES <- stack_hourly_temps(temps_ES, latitude = 43.48)
hourtemps_DE <- stack_hourly_temps(temps_DE, latitude = 50.62)

# Only take dataframe with data
hourtemps_ES <- hourtemps_ES[[1]]
hourtemps_DE <- hourtemps_DE[[1]]

# Read phenology data and prepare it for fitting
cal_data_ES <- load_temperature_scenarios("data/Hajars_data/calibration/", prefix = "")
cal_data_DE <- load_temperature_scenarios("data/outputs/pheno/", prefix = "")
cal_data_DE <- cal_data_DE[str_detect(names(cal_data_DE), "calibration_FB")]
names(cal_data_DE) <- paste0("data_", c("Berlepsch", "CoxOrange", "GoldenDelicious",
                                        "JamesGrieve", "RoterBoskoop"))

# Validation data
val_data_ES <- load_temperature_scenarios("data/Hajars_data/validation/", prefix = "")
val_data_DE <- load_temperature_scenarios("data/outputs/pheno/", prefix = "")
val_data_DE <- val_data_DE[str_detect(names(val_data_DE), "validation_FB")]
names(val_data_DE) <- names(cal_data_DE)

# Merge data from Spain and Germany into one list
cal_data <- c(cal_data_ES, cal_data_DE)
val_data <- c(val_data_ES, val_data_DE)

# Split hourly temperature data to the respective seasons (start in august, end in july of following year)
SeasonList <- list()
eval_SeasonList <- list()

for (i in 1 : length(cal_data)) {
  
  if (i < 12) weather <- hourtemps_ES else weather <- hourtemps_DE
  
  years_cal <- cal_data[[i]]$Year
  years_val <- val_data[[i]]$Year
  
  SeasonList[[i]] <- genSeasonList(weather, mrange = c(8, 6), years = years_cal)
  eval_SeasonList[[i]] <- genSeasonList(weather, mrange = c(8, 6), years = years_val)
}

# Add the names of the cultivars to the lists generated
names(SeasonList) <- names(cal_data)

# The same for the validation list
names(eval_SeasonList) <- names(val_data)


###### Run 1 #######
# You need to specify how many cultivars you want to fit jointly
n_cult <- length(cal_data)

# Homogenize name of column for pheno
for (i in 1 : length(cal_data)) names(cal_data[[i]])[which(names(cal_data[[i]]) == "Bloom_full_JDay")] <- "pheno"

# Combine phenology observations used for model training
bloomJDays_list <- lapply(cal_data, function (x) x[["pheno"]])

# The number of parameters fitted increase, but overall they decrease because only yc, zc and s1 are cultivar specific
lower_01 <- c(rep( 10, n_cult), rep( 20, n_cult), rep(0.05, n_cult),  0, 2300.0,  8900.0, 5300.0,      4.9e13,  0,  0,  0,  0.05)
par_01 <-   c(rep( 50, n_cult), rep(230, n_cult), rep(0.50, n_cult), 25, 3372.8,  9900.3, 6319.5, 5.939917e13,  4, 36,  4,  1.60)
upper_01 <- c(rep(100, n_cult), rep(550, n_cult), rep(1.00, n_cult), 40, 4300.0, 10900.0, 7300.0,      6.9e13, 20, 50, 20, 50.00)


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
temp_values <- seq(-5, 40, 0.1)

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
yc_01 <- Fit_res_01$par[1 : n_cult]

# zc
c(min(Fit_res_01$par[(n_cult + 1) : (n_cult * 2)]), max(Fit_res_01$par[(n_cult + 1) : (n_cult * 2)]))
zc_01 <- Fit_res_01$par[(n_cult + 1) : (n_cult * 2)]

# s1
c(min(Fit_res_01$par[((n_cult * 2) + 1) : (n_cult * 3)]), max(Fit_res_01$par[((n_cult * 2) + 1) : (n_cult * 3)]))
s1_01 <- Fit_res_01$par[((n_cult * 2) + 1) : (n_cult * 3)]

# Remaining pars
Fit_res_01$par[((n_cult * 3) + 1) : length(par_01)]

# Adjust the values
lower_02 <- c(rep(10, n_cult),  rep(50, n_cult), rep(0.01, n_cult), 17, 2372.0,  8900.0, 5292.0, 4.9e13,  0, 17,  0,  0.05)
par_02 <-   Fit_res_01$par
upper_02 <- c(rep(80, n_cult), rep(550, n_cult), rep(1.50, n_cult), 37, 4372.0, 10900.0, 7292.0, 6.9e13, 16, 37, 12, 40.00)


# Please note that you need to specify the number of cultivars in the function call, otherwise nothing changed
Fit_res_02 <- phenologyFitter_combined_fitting(par.guess = par_02, 
                                               modelfn = PhenoFlex_GDHwrapper,
                                               bloomJDays = bloomJDays_list,
                                               SeasonList = SeasonList,
                                               n_cult = n_cult,
                                               lower = lower_02,
                                               upper = upper_02,
                                               seed = 001,
                                               control = list(smooth = FALSE, verbose = FALSE, maxit = 1000,
                                                              nb.stop.improvement = 20))

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
yc_02 <- Fit_res_02$par[1 : n_cult]

# zc
c(min(Fit_res_02$par[(n_cult + 1) : (n_cult * 2)]), max(Fit_res_02$par[(n_cult + 1) : (n_cult * 2)]))
zc_02 <- Fit_res_02$par[(n_cult + 1) : (n_cult * 2)]

# s1
c(min(Fit_res_02$par[((n_cult * 2) + 1) : (n_cult * 3)]), max(Fit_res_02$par[((n_cult * 2) + 1) : (n_cult * 3)]))
s1_02 <- Fit_res_02$par[((n_cult * 2) + 1) : (n_cult * 3)]

# Remaining pars
Fit_res_02$par[((n_cult * 3) + 1) : length(par_02)]

# Adjust the values
lower_03 <- c(rep(15, n_cult), rep(150, n_cult), rep(0.001, n_cult), 20, 2300.0,  9900.0, 5200.0, 4.9e13, -1, 20,  -1,  0.05)
par_03 <-   Fit_res_02$par
upper_03 <- c(rep(70, n_cult), rep(490, n_cult), rep(1.650, n_cult), 34, 4300.0, 10900.0, 7200.0, 6.9e13, 12, 35, 12, 25.00)

# Please note that you need to specify the number of cultivars in the function call, otherwise nothing changed
Fit_res_03 <- phenologyFitter_combined_fitting(par.guess = par_03, 
                                               modelfn = PhenoFlex_GDHwrapper,
                                               bloomJDays = bloomJDays_list,
                                               SeasonList = SeasonList,
                                               n_cult = n_cult,
                                               lower = lower_03,
                                               upper = upper_03,
                                               seed = 001,
                                               control = list(smooth = FALSE, verbose = FALSE, maxit = 2200,
                                                              nb.stop.improvement = 40))

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
yc_03 <- Fit_res_03$par[1 : n_cult]

# zc
c(min(Fit_res_03$par[(n_cult + 1) : (n_cult * 2)]), max(Fit_res_03$par[(n_cult + 1) : (n_cult * 2)]))
zc_03 <- Fit_res_03$par[(n_cult + 1) : (n_cult * 2)]

# s1
c(min(Fit_res_03$par[((n_cult * 2) + 1) : (n_cult * 3)]), max(Fit_res_03$par[((n_cult * 2) + 1) : (n_cult * 3)]))
s1_03 <- Fit_res_03$par[((n_cult * 2) + 1) : (n_cult * 3)]

# Remaining pars
Fit_res_03$par[((n_cult * 3) + 1) : length(par_03)]

# Adjust the values
lower_04 <- c(rep(19, n_cult), rep(200, n_cult), rep(0.001, n_cult), 22, 2600.0,  9300.0, 5500.0, 5.4e13,   0, 22,  0,  0.05)
par_04 <-   Fit_res_03$par
upper_04 <- c(rep(58, n_cult), rep(430, n_cult), rep(2.000, n_cult), 32, 4000.0, 10500.0, 6900.0, 6.4e13,  10, 32, 10, 15.00)


# Please note that you need to specify the number of cultivars in the function call, otherwise nothing changed
Fit_res_04 <- phenologyFitter_combined_fitting(par.guess = par_04, 
                                               modelfn = PhenoFlex_GDHwrapper,
                                               bloomJDays = bloomJDays_list,
                                               SeasonList = SeasonList,
                                               n_cult = n_cult,
                                               lower = lower_04,
                                               upper = upper_04,
                                               seed = 001,
                                               control = list(smooth = FALSE, verbose = FALSE, maxit = 2300,
                                                              nb.stop.improvement = 50))

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

# zc
c(min(Fit_res_04$par[(n_cult + 1) : (n_cult * 2)]), max(Fit_res_04$par[(n_cult + 1) : (n_cult * 2)]))

# s1
c(min(Fit_res_04$par[((n_cult * 2) + 1) : (n_cult * 3)]), max(Fit_res_04$par[((n_cult * 2) + 1) : (n_cult * 3)]))

# Remaining pars
Fit_res_04$par[((n_cult * 3) + 1) : length(par_04)]

# Adjust the values
lower_05 <- c(rep(24, n_cult), rep(195, n_cult), rep(0.001, n_cult), 22, 2600.0,  9200.0, 5500.0, 5.4e13,  1, 22,  1,  0.05)
par_05 <-   Fit_res_04$par
upper_05 <- c(rep(53, n_cult), rep(424, n_cult), rep(2.750, n_cult), 32, 4000.0, 10600.0, 6900.0, 6.4e13,  9, 32,  9, 15.00)

# Please note that you need to specify the number of cultivars in the function call, otherwise nothing changed
Fit_res_05 <- phenologyFitter_combined_fitting(par.guess = par_05, 
                                               modelfn = PhenoFlex_GDHwrapper,
                                               bloomJDays = bloomJDays_list,
                                               SeasonList = SeasonList,
                                               n_cult = n_cult,
                                               lower = lower_05,
                                               upper = upper_05,
                                               seed = 001,
                                               control = list(smooth = FALSE, verbose = FALSE, maxit = 3500,
                                                              nb.stop.improvement = 120))

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
# Combining data from Hajar's repo to mine is tricky. Will need to organize the column names.
cal_data[1 : 11] <- lapply(cal_data[1 : 11], function (x) select(x, -X))
cal_data[12 : 16] <- lapply(cal_data[12 : 16], function (x) select(x, -Predicted_r10, -Error_r10))

# Add the column species and cultivar to the elements of the list
for (i in 1 : 11){
  
  # Substract the name of the cultivar and then add it to the existing dataframes
  cal_data[[i]] <- data.frame(Species = "Apple",
                              Cultivar = substr(names(cal_data), 6, nchar(names(cal_data)))[i],
                              cal_data[[i]])}


# Use bind_rows to merge the calibration_list
summary_calibration <- bind_rows(cal_data)

# Add the predicted JDay for bloom using the fifth run
summary_calibration[["Predicted"]] <- unlist(pbloomJDays_05)

# Compute the error
summary_calibration[["Error"]] <- summary_calibration$pheno - summary_calibration$Predicted



###### Validation ######
# Use the estimated parameters to compute bloom dates in the validation datasets. Validation datasets are located in
# val_data. In addition, eval_SeasonList contains all weather data for the seasons used in the validation
# datasets.
# Fix some columns to facilitate the further use of the list.
val_data[1 : 11] <- lapply(val_data[1 : 11], function (x) select(x, -X))
val_data[12 : 16] <- lapply(val_data[12 : 16], function (x) select(x, Species, Cultivar, Year, Bloom_full_JDay))

# Add the column species and cultivar to the elements of the list
for (i in 1 : 11){
  
  # Substract the name of the cultivar and then add it to the existing dataframes
  val_data[[i]] <- data.frame(Species = "Apple",
                              Cultivar = substr(names(val_data), 6, nchar(names(val_data)))[i],
                              val_data[[i]])}

# Fix the name of the pheno column in German datasets
for (i in 12 : 16) names(val_data[[i]]) <- c("Species", "Cultivar", "Year", "pheno")


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
summary_validation <- data.frame(bind_rows(val_data),
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
  names(pars_sum) <- names(val_data)
  
  # Allocate the results of the loop "j" into the final list
  par_sums[[j]] <- pars_sum}

# Name the elements of the final list according to the Fit_res element
names(par_sums) <- paste0("Fit_res_0", 1 : 5)

# Summarize the evolution of parameters
summary_params <- bind_rows(lapply(par_sums, bind_rows, .id = "Cultivar"), .id = "Fit_res")


# Save all outputs to folder as csv
write.csv(summary_calibration, "results/01_summary_cal_both.csv", row.names = FALSE)
write.csv(summary_validation, "results/04_summary_val_both.csv", row.names = FALSE)
write.csv(summary_params, "results/07_summary_pars_both.csv", row.names = FALSE)
write.csv(summary_metrics, "results/10_summary_metrics_both.csv", row.names = FALSE)






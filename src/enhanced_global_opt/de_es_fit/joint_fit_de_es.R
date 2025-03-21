

#' Function used in the optimization process to generate bloom data

#' 

#' This is the evaluation function for the fitting process. It takes the "new" parameters

#' (theta_star, theta_c, tau, pie_c) instead of the old ones (E0, E1, A0, A1). This function

#' allows combined fitting of cultivars. It will only estimate cultivar-specific parameters of yc, zc and s1.

#' 

#' This function will be best used in combination with the \link[MEIGOR]{MEIGO} function 

#' to estimate PhenoFlex parameter values. The order of bloomJDays and SeasonList should be the same. 

#' The order will correspond to estimated cultivar-specific parameters 

#' 

#' 

#' @param x traditional model parameters of PhenoFlex in the order yc, zc, s1, Tu, theta_star, theta_c, tau, pie_c, Tf, Tc, Tb, slope. The parameters

#' yc, zc and s1 need to be replicated as often as cultivars are fitted. 

#' @param modelfn function used within the evaluation function to calculate the actual bloomday, often we use

#' the 'custom_GDH_wrapper' function for that

#' @param bloomJDays list of numeric vectors containing the observed bloom dates in day of the year format

#' @param SeasonList list of lists, containing the hourly temperatures for the individual phenological seasons. One element per cultivar.

#' The list of each cultivar should contain data.frames for eachs season

#' with the columns "Temp" (for the hourly temperature) and "JDay" for the corresponding Julian day. Is usually

#' generated using \link[chillR]{genSeasonList}

#' @param n_cult Numeric, specifying the numbers of cultivars fitted

#' @param na_penalty numeric, value which is used when the model fails to generate a prediction

#' for the bloom date. By default 365

#' @param return_bloom_days boolean, by default set FALSE. Controls if the evaluation score or the predicted bloom dates get returned by the function

#' @return list with two elements. First is called 'f' and contains the residual sum of squares of the model. The 

#' second is 'g' which is the values of the additional model constraints defined in the function.

#' 

#' @author Lars Caspersen

#' @keywords utility

#' @importFrom nleqslv nleqslv

#'  

#' @export evaluation_function_meigo_nonlinear_combined



#bloomJDays is now a dataframe with the column pheno

#seasonlist elements must be in the same order as bloomJDays elements (also within the elements)

evaluation_function_meigo_nonlinear_combined <- function(x, 
                                                         
                                                         modelfn,
                                                         
                                                         bloomJDays,
                                                         
                                                         SeasonList,
                                                         
                                                         n_cult,
                                                         
                                                         na_penalty = 365,
                                                         
                                                         return_bloom_days = FALSE){
  
  
  
  #innput:
  
  #         x is the parameters in meigo
  
  #         modelfn is the function used to calculate the bloomdays
  
  #         SeasonList contains the weather data for the different years
  
  #         na_penalty is the value used if the model failed to predict any day with the current set of parameters for a particular year
  
  
  
  #output: inequality constraints g
  
  #        model performance value F
  
  
  
  
  
  
  
  #instead of A0, A1, E0 and E1 we have now theta*, theta_c, Tau(thetha*) and pie_c
  
  #--> we need to solve now a non-linear system to calculate A0, A1, E0 and E1 based on the parameters
  
  #    ('yc', 'zc', 's1', 'Tu', 'theta*', 'theta_c', 'Tau(thetha*)', 'pie_c', 'Tf', 'Tc', 'Tb',  'slope')
  
  #x<- c(40,   190,   0.5,  25,   279,      287,       28,             26,       4,   36,    4,    1.60)
  
  
  
  params<-numeric(4)
  
  
  
  params[1] <- x[5+((n_cult -1) * 3)]   #theta*
  
  params[2] <- x[6+((n_cult -1) * 3)]    #theta_c
  
  params[3] <- x[7+((n_cult -1) * 3)]    #Tau(thetha*)
  
  params[4] <- x[8+((n_cult -1) * 3)]     #pi_c
  
  
  
  
  
  output<-nleqslv::nleqslv(c(500, 15000), solve_nle, jac=NULL, params, xscalm="auto", method="Newton",
                           
                           control=list(trace=0,allowSingular=TRUE))
  
  
  
  
  
  #This is a numerical method which can produce non-convergence. Check this
  
  if (output$termcd >= 3){
    
    #if the nle algorithm has stalled just discard this solution
    
    E0<-NA; E1<-NA; A0<-NA; A1<-NA
    
    return(list(F=10^6, g=rep(10^6,5)))
    
    
    
    #You would add here a flag to let your optimization procedure know
    
    #That this solution should be ignored by lack of convergence
    
    
    
  } else {
    
    
    
    E0 <- output$x[1]
    
    E1 <- output$x[2]
    
    
    
    #A1 and A0 can be calculated through Equations 36 and 37
    
    
    
    q=1/params[1]-1/params[2]
    
    
    
    A1 <- -exp(E1/params[1])/params[3]*log(1-exp((E0-E1)*q))
    
    A0 <- A1*exp((E0-E1)/params[2])
    
  }
  
  
  
  
  
  #change the name of the parameters, dont know if necessary
  
  par <- x
  
  par[(5:8)+((n_cult -1) * 3)] <- c(E0, E1, A0, A1)
  
  
  
  #split the parameters and reconstruct them for the different cultivar
  
  
  
  pred_bloom <- NULL
  
  
  
  #loop over the cultivars, calculate predicted days for each cultivar
  
  rss <- purrr::map_dbl(1:length(SeasonList), function(i){
    
    
    
    #extract the parameters
    
    par_cult <- par[c(i, 1+((n_cult -1))+i, 2+((n_cult -1)*2)+i, (length(par)-8):length(par))]
    
    #predict the bloom
    
    pred_bloom <- unlist(lapply(X = SeasonList[[i]], FUN = modelfn, par = par_cult))
    
    pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, no = pred_bloom)
    
    #calculate for each cultivar the rss
    
    rss <- sum((pred_bloom - bloomJDays[[i]])^2)
    
    return(rss)
    
    
    
  })
  
  
  
  
  
  #calculate the model performance value
  
  F <- sum(rss)
  
  
  
  
  
  
  
  #####
  
  #inequality constraints
  
  #####
  
  
  
  #this is the vector containing the values for the inequality constraints
  
  #at first initialize the vector
  
  g <- rep(0,5)
  
  
  
  
  
  #equality constraints should be always stated before inequality constraints, according to meigo vignette
  
  #(but we dont have any in this case)
  
  
  
  
  
  #inequality constraints are reformulated as differences
  
  
  
  #Tu >= Tb
  
  g[1] <- par[length(par)-8] - par[length(par)-1]
  
  #Tx >= Tb
  
  g[2] <- par[length(par)-2] - par[length(par)-1]
  
  #Tc >= Tu
  
  g[3] <- par[length(par)-2] - par[length(par)-8]
  
  
  
  
  
  #the q10 value should be between 1.5 and 3.5
  
  #q10 formation = exp((10*E0)/(T2-T1))
  
  #q10 destruction = exp((10*E1)/(T2-T1))
  
  #T1 = 297
  
  #T2 = 279
  
  #parameters T1 and T2 chosen after Egea 2020
  
  #ranges for q10 are provided in c_L and c_U
  
  g[4] <- exp((10 * par[length(par)-7]) / (297 * 279))
  
  
  
  g[5] <- exp((10 * par[length(par)-6]) / (297 * 279))
  
  
  
  
  
  if(return_bloom_days == FALSE){
    
    #output
    
    return(list(F=F, g=g))
    
  } else{
    
    return(pred_bloom)
    
  }
  
}



library(LarsChill)
library(Metrics)
library(readxl)
library(tidyverse)
install.packages("BiocManager")
BiocManager::install("MEIGOR")
source("code/custom_ess.R")
source("code/load_save_fitting_results.R")
source("code/helper_functions.R")

#read temperature data
temp <- read_excel('MaxMin1978and2020.xlsx')
temp_de <- read.csv('Edu/weather_data_final.csv')

#make day and month numeric
temp$Month <- as.numeric(temp$Month)
temp$Day <- as.numeric(temp$Day)

#check if the temperature data is complete
check_temperature_record(temp)

#get hourly data from daily data
hourtemps <- stack_hourly_temps(temp, latitude = 43.48)
hourtemps_de <- stack_hourly_temps(temp_de, latitude = 50.62)

#only take dataframe with data
hourtemps_de <- hourtemps_de[[1]]
hourtemps <- hourtemps[[1]]

#read phenology data 
pheno_berlepsch <- read.csv("pheno_data_de/berlepsch.csv")
pheno_cox_orange <- read.csv("pheno_data_de/cox-orange.csv")
pheno_golden_delicious <- read.csv("pheno_data_de/golden-delicious.csv")
pheno_james_grieve <- read.csv("pheno_data_de/james-grieve.csv")
pheno_roter_boskoop <- read.csv("pheno_data_de/roter-boskoop.csv")

pheno_Blanquina <- read_excel("pheno_data/Blanquina.xlsx")
pheno_Clara <- read_excel("pheno_data/Clara.xlsx")
pheno_Collaos <- read_excel("pheno_data/Collaos.xlsx")
pheno_Coloradona <- read_excel("pheno_data/Coloradona.xlsx")
pheno_DelaRiega <- read_excel("pheno_data/DelaRiega.xlsx")
pheno_Perezosa <- read_excel("pheno_data/Perezosa.xlsx")
pheno_Perico <- read_excel("pheno_data/Perico.xlsx")
pheno_Raxao <- read_excel("pheno_data/Raxao.xlsx")
pheno_Teorica <- read_excel("pheno_data/Teorica.xlsx")
pheno_Verdialona <- read_excel("pheno_data/Verdialona.xlsx")
pheno_Xuanina <- read_excel("pheno_data/Xuanina.xlsx")

## #drop rows without data
pheno_data_clean_Blanquina <- na.omit(pheno_Blanquina)
pheno_data_clean_Clara <- na.omit(pheno_Clara)
pheno_data_clean_Collaos <- na.omit(pheno_Collaos)
pheno_data_clean_Coloradona <- na.omit(pheno_Coloradona)
pheno_data_clean_DelaRiega <- na.omit(pheno_DelaRiega)
pheno_data_clean_Perezosa <- na.omit(pheno_Perezosa)
pheno_data_clean_Perico <- na.omit(pheno_Perico)
pheno_data_clean_Raxao <- na.omit(pheno_Raxao)
pheno_data_clean_Teorica <- na.omit(pheno_Teorica)
pheno_data_clean_Verdialona <- na.omit(pheno_Verdialona)
pheno_data_clean_Xuanina <- na.omit(pheno_Xuanina)

#
# #change years to numeric
pheno_data_clean_Blanquina$Year <- as.numeric(pheno_data_clean_Blanquina$Year)
pheno_data_clean_Clara$Year <- as.numeric(pheno_data_clean_Clara$Year)
pheno_data_clean_Collaos$Year <- as.numeric(pheno_data_clean_Collaos$Year)
pheno_data_clean_Coloradona$Year <- as.numeric(pheno_data_clean_Coloradona$Year)
pheno_data_clean_DelaRiega$Year <- as.numeric(pheno_data_clean_DelaRiega$Year)
pheno_data_clean_Perezosa$Year <- as.numeric(pheno_data_clean_Perezosa$Year)
pheno_data_clean_Perico$Year <- as.numeric(pheno_data_clean_Perico$Year)
pheno_data_clean_Raxao$Year <- as.numeric(pheno_data_clean_Raxao$Year)
pheno_data_clean_Teorica$Year <- as.numeric(pheno_data_clean_Teorica$Year)
pheno_data_clean_Verdialona$Year <- as.numeric(pheno_data_clean_Verdialona$Year)
pheno_data_clean_Xuanina$Year <- as.numeric(pheno_data_clean_Xuanina$Year)

#set seed for reproducibility 
set.seed(51)
cal_data_list <- list()
val_data_list <- list()

for(i in 1:10){
  
  #i <- 1  
  
  # Create validation data by excluding the calibration data
  # calibration
  cal_berlepsch <- pheno_berlepsch[sample(nrow(pheno_berlepsch), 0.8 * nrow(pheno_berlepsch)), ]
  cal_cox_orange <- pheno_cox_orange[sample(nrow(pheno_cox_orange), 0.8 * nrow(pheno_cox_orange)), ]
  cal_golden_delicious <- pheno_golden_delicious[sample(nrow(pheno_golden_delicious), 0.8 * nrow(pheno_golden_delicious)), ]
  cal_james_grieve <- pheno_james_grieve[sample(nrow(pheno_james_grieve), 0.8 * nrow(pheno_james_grieve)), ]
  cal_roter_boskoop <- pheno_roter_boskoop[sample(nrow(pheno_roter_boskoop), 0.8 * nrow(pheno_roter_boskoop)), ]
  
  val_berlepsch <- pheno_berlepsch[!(rownames(pheno_berlepsch) %in% rownames(cal_berlepsch)), ]
  val_cox_orange <- pheno_cox_orange[!(rownames(pheno_cox_orange) %in% rownames(cal_cox_orange)), ]
  val_golden_delicious <- pheno_golden_delicious[!(rownames(pheno_golden_delicious) %in% rownames(cal_golden_delicious)), ]
  val_james_grieve <- pheno_james_grieve[!(rownames(pheno_james_grieve) %in% rownames(cal_james_grieve)), ]
  val_roter_boskoop <- pheno_roter_boskoop[!(rownames(pheno_roter_boskoop) %in% rownames(cal_roter_boskoop)), ]
  
  val_years_Blanquina <- sample(pheno_data_clean_Blanquina$Year, 7)
  val_Blanquina <- pheno_data_clean_Blanquina[pheno_data_clean_Blanquina$Year %in% val_years_Blanquina,]
  
  val_years_Clara <- sample(pheno_data_clean_Clara$Year, 7)
  val_Clara <- pheno_data_clean_Clara[pheno_data_clean_Clara$Year %in% val_years_Clara,]
  
  val_years_Collaos <- sample(pheno_data_clean_Collaos$Year, 7)
  val_Collaos <- pheno_data_clean_Collaos[pheno_data_clean_Collaos$Year %in% val_years_Collaos,]
  
  val_years_Coloradona <- sample(pheno_data_clean_Coloradona$Year, 7)
  val_Coloradona <- pheno_data_clean_Coloradona[pheno_data_clean_Coloradona$Year %in% val_years_Coloradona,]
  
  val_years_DelaRiega <- sample(pheno_data_clean_DelaRiega$Year, 7)
  val_DelaRiega <- pheno_data_clean_DelaRiega[pheno_data_clean_DelaRiega$Year %in% val_years_DelaRiega,]
  
  val_years_Perezosa <- sample(pheno_data_clean_Perezosa$Year, 7)
  val_Perezosa <- pheno_data_clean_Perezosa[pheno_data_clean_Perezosa$Year %in% val_years_Perezosa,]
  
  val_years_Perico <- sample(pheno_data_clean_Perico$Year, 7)
  val_Perico <- pheno_data_clean_Perico[pheno_data_clean_Perico$Year %in% val_years_Perico,]
  
  val_years_Raxao <- sample(pheno_data_clean_Raxao$Year, 7)
  val_Raxao <- pheno_data_clean_Raxao[pheno_data_clean_Raxao$Year %in% val_years_Raxao,]
  
  val_years_Teorica <- sample(pheno_data_clean_Teorica$Year, 7)
  val_Teorica <- pheno_data_clean_Teorica[pheno_data_clean_Teorica$Year %in% val_years_Teorica,]
  
  val_years_Verdialona <- sample(pheno_data_clean_Verdialona$Year, 7)
  val_Verdialona <- pheno_data_clean_Verdialona[pheno_data_clean_Verdialona$Year %in% val_years_Verdialona,]
  
  val_years_Xuanina <- sample(pheno_data_clean_Xuanina$Year, 7)
  val_Xuanina <- pheno_data_clean_Xuanina[pheno_data_clean_Xuanina$Year %in% val_years_Xuanina,]
  
  val_data_list[[i]] <- list(val_berlepsch,val_cox_orange, val_golden_delicious,val_james_grieve, val_roter_boskoop,
                             val_Blanquina,val_Clara,val_Collaos,val_Coloradona,val_DelaRiega,val_Perezosa,val_Perico,val_Raxao,
                             val_Teorica,val_Verdialona,val_Xuanina)
  
  # calibration
  cal_years_Blanquina <- pheno_data_clean_Blanquina$Year[!(pheno_data_clean_Blanquina$Year %in% val_years_Blanquina)]
  cal_Blanquina <- pheno_data_clean_Blanquina[pheno_data_clean_Blanquina$Year %in% cal_years_Blanquina,]
  
  cal_years_Clara <- pheno_data_clean_Clara$Year[!(pheno_data_clean_Clara$Year %in% val_years_Clara)]
  cal_Clara <- pheno_data_clean_Clara[pheno_data_clean_Clara$Year %in% cal_years_Clara,]
  
  cal_years_Collaos <- pheno_data_clean_Collaos$Year[!(pheno_data_clean_Collaos$Year %in% val_years_Collaos)]
  cal_Collaos <- pheno_data_clean_Collaos[pheno_data_clean_Collaos$Year %in% cal_years_Collaos,]
  
  cal_years_Coloradona <- pheno_data_clean_Coloradona$Year[!(pheno_data_clean_Coloradona$Year %in% val_years_Coloradona)]
  cal_Coloradona <- pheno_data_clean_Coloradona[pheno_data_clean_Coloradona$Year %in% cal_years_Coloradona,]
  
  cal_years_DelaRiega <- pheno_data_clean_DelaRiega$Year[!(pheno_data_clean_DelaRiega$Year %in% val_years_DelaRiega)]
  cal_DelaRiega <- pheno_data_clean_DelaRiega[pheno_data_clean_DelaRiega$Year %in% cal_years_DelaRiega,]
  
  cal_years_Perezosa <- pheno_data_clean_Perezosa$Year[!(pheno_data_clean_Perezosa$Year %in% val_years_Perezosa)]
  cal_Perezosa <- pheno_data_clean_Perezosa[pheno_data_clean_Perezosa$Year %in% cal_years_Perezosa,]
  
  cal_years_Perico <- pheno_data_clean_Perico$Year[!(pheno_data_clean_Perico$Year %in% val_years_Perico)]
  cal_Perico <- pheno_data_clean_Perico[pheno_data_clean_Perico$Year %in% cal_years_Perico,]
  
  cal_years_Raxao <- pheno_data_clean_Raxao$Year[!(pheno_data_clean_Raxao$Year %in% val_years_Raxao)]
  cal_Raxao <- pheno_data_clean_Raxao[pheno_data_clean_Raxao$Year %in% cal_years_Raxao,]
  
  cal_years_Teorica <- pheno_data_clean_Teorica$Year[!(pheno_data_clean_Teorica$Year %in% val_years_Teorica)]
  cal_Teorica <- pheno_data_clean_Teorica[pheno_data_clean_Teorica$Year %in% cal_years_Teorica,]
  
  cal_years_Verdialona <- pheno_data_clean_Verdialona$Year[!(pheno_data_clean_Verdialona$Year %in% val_years_Verdialona)]
  cal_Verdialona <- pheno_data_clean_Verdialona[pheno_data_clean_Verdialona$Year %in% cal_years_Verdialona,]
  
  cal_years_Xuanina <- pheno_data_clean_Xuanina$Year[!(pheno_data_clean_Xuanina$Year %in% val_years_Xuanina)]
  cal_Xuanina <- pheno_data_clean_Xuanina[pheno_data_clean_Xuanina$Year %in% cal_years_Xuanina,]
  
  cal_data_list[[i]] <- list(cal_berlepsch, cal_cox_orange, cal_golden_delicious, cal_james_grieve, cal_roter_boskoop,
                             cal_Blanquina,cal_Clara,cal_Collaos,cal_Coloradona,cal_DelaRiega,cal_Perezosa,cal_Perico,cal_Raxao,
                             cal_Teorica,cal_Verdialona,cal_Xuanina)
  
  # #save validation and calibration data
  write.csv(cal_berlepsch, paste0("enhanced_global_opt/de_es_fit/data_split/cal_", i,"_berlepsch.csv")) 
  write.csv(val_berlepsch, paste0("enhanced_global_opt/de_es_fit/data_split/val_", i,"_berlepsch.csv"))
  
  write.csv(cal_cox_orange, paste0("enhanced_global_opt/de_es_fit/data_split/cal_", i, "_cox_orange.csv")) 
  write.csv(val_cox_orange, paste0("enhanced_global_opt/de_es_fit/data_split/val_", i, "_cox_orange.csv"))
  
  write.csv(cal_golden_delicious,paste0("enhanced_global_opt/de_es_fit/data_split/cal_", i, "_golden_delicious.csv") )
  write.csv(val_golden_delicious,paste0("enhanced_global_opt/de_es_fit/data_split/val_", i, "_golden_delicious.csv"))
  
  write.csv(cal_james_grieve,paste0("enhanced_global_opt/de_es_fit/data_split/cal_", i, "_james_grieve.csv") )
  write.csv(val_james_grieve,paste0("enhanced_global_opt/de_es_fit/data_split/val_", i, "_james_grieve.csv"))
  
  write.csv(cal_roter_boskoop,paste0("enhanced_global_opt/de_es_fit/data_split/cal_", i, "_roter_boskoop.csv")) 
  write.csv(val_roter_boskoop,paste0("enhanced_global_opt/de_es_fit/data_split/val_", i, "_roter_boskoop.csv"))
  
  write.csv(cal_Blanquina, paste0("enhanced_global_opt/de_es_fit/data_split/cal_", i,"_Blanquina.csv")) 
  write.csv(val_Blanquina, paste0("enhanced_global_opt/de_es_fit/data_split/val_", i,"_Blanquina.csv"))
  
  write.csv(cal_Clara, paste0("enhanced_global_opt/de_es_fit/data_split/cal_", i,"_Clara.csv")) 
  write.csv(val_Clara, paste0("enhanced_global_opt/de_es_fit/data_split/val_", i,"_Clara.csv"))
  
  write.csv(cal_Collaos, paste0("enhanced_global_opt/de_es_fit/data_split/cal_", i,"_Collaos.csv")) 
  write.csv(val_Collaos, paste0("enhanced_global_opt/de_es_fit/data_split/val_", i,"_Collaos.csv"))
  
  write.csv(cal_Coloradona, paste0("enhanced_global_opt/de_es_fit/data_split/cal_", i,"_Coloradona.csv")) 
  write.csv(val_Coloradona, paste0("enhanced_global_opt/de_es_fit/data_split/val_", i,"_Coloradona.csv"))
  
  write.csv(cal_DelaRiega, paste0("enhanced_global_opt/de_es_fit/data_split/cal_", i,"_DelaRiega.csv")) 
  write.csv(val_DelaRiega, paste0("enhanced_global_opt/de_es_fit/data_split/val_", i,"_DelaRiega.csv"))
  
  write.csv(cal_Perezosa, paste0("enhanced_global_opt/de_es_fit/data_split/cal_", i,"_Perezosa.csv")) 
  write.csv(val_Perezosa, paste0("enhanced_global_opt/de_es_fit/data_split/val_", i,"_Perezosa.csv"))
  
  write.csv(cal_Perico, paste0("enhanced_global_opt/de_es_fit/data_split/cal_", i,"_Perico.csv")) 
  write.csv(val_Perico, paste0("enhanced_global_opt/de_es_fit/data_split/val_", i,"_Perico.csv"))
  
  write.csv(cal_Raxao, paste0("enhanced_global_opt/de_es_fit/data_split/cal_", i,"_Raxao.csv")) 
  write.csv(val_Raxao, paste0("enhanced_global_opt/de_es_fit/data_split/val_", i,"_Raxao.csv"))
  
  write.csv(cal_Teorica, paste0("enhanced_global_opt/de_es_fit/data_split/cal_", i,"_Teorica.csv")) 
  write.csv(val_Teorica, paste0("enhanced_global_opt/de_es_fit/data_split/val_", i,"_Teorica.csv"))
  
  write.csv(cal_Verdialona, paste0("enhanced_global_opt/de_es_fit/data_split/cal_", i,"_Verdialona.csv")) 
  write.csv(val_Verdialona, paste0("enhanced_global_opt/de_es_fit/data_split/val_", i,"_Verdialona.csv"))
  
  write.csv(cal_Xuanina, paste0("enhanced_global_opt/de_es_fit/data_split/cal_", i,"_Xuanina.csv")) 
  write.csv(val_Xuanina, paste0("enhanced_global_opt/de_es_fit/data_split/val_", i,"_Xuanina.csv"))
  
  
  
}


#contain the repetitions
res_loop <- list()
SeasonList_list <- list()
eval_SeasonList_list <- list()
bloomJDays_eval_list_squared <- list()
bloomJDays_list_squared <- list()

for(i in 1:10){
  
  #i <-1
  
  #split hourly temperature data to the respective seasons (start in august, end in july of following year)
  SeasonList_berlepsch <- genSeasonList(hourtemps_de, mrange = c(8, 6), years = cal_data_list[[i]][[1]]$Year)
  eval_SeasonList_berlepsch <- genSeasonList(hourtemps_de, mrange = c(8, 6), years = val_data_list[[i]][[1]]$Year)
  
  SeasonList_cox_orange <- genSeasonList(hourtemps_de, mrange = c(8, 6), years = cal_data_list[[i]][[2]]$Year)
  eval_SeasonList_cox_orange <- genSeasonList(hourtemps_de, mrange = c(8, 6), years = val_data_list[[i]][[2]]$Year)
  
  SeasonList_golden_delicious <- genSeasonList(hourtemps_de, mrange = c(8, 6), years = cal_data_list[[i]][[3]]$Year)
  eval_SeasonList_golden_delicious <- genSeasonList(hourtemps_de, mrange = c(8, 6), years = val_data_list[[i]][[3]]$Year)
  
  SeasonList_james_grieve <- genSeasonList(hourtemps_de, mrange = c(8, 6), years = cal_data_list[[i]][[4]]$Year)
  eval_SeasonList_james_grieve <- genSeasonList(hourtemps_de, mrange = c(8, 6), years = val_data_list[[i]][[4]]$Year)
  
  SeasonList_roter_boskoop <- genSeasonList(hourtemps_de, mrange = c(8, 6), years = cal_data_list[[i]][[5]]$Year)
  eval_SeasonList_roter_boskoop <- genSeasonList(hourtemps_de, mrange = c(8, 6), years = val_data_list[[i]][[5]]$Year)
  
  SeasonList_Blanquina <- genSeasonList(hourtemps, mrange = c(8, 6), years = cal_data_list[[i]][[6]]$Year)
  eval_SeasonList_Blanquina <- genSeasonList(hourtemps, mrange = c(8, 6), years = val_data_list[[i]][[6]]$Year)
  
  SeasonList_Clara <- genSeasonList(hourtemps, mrange = c(8, 6), years = cal_data_list[[i]][[7]]$Year)
  eval_SeasonList_Clara <- genSeasonList(hourtemps, mrange = c(8, 6), years = val_data_list[[i]][[7]]$Year)
  
  SeasonList_Collaos <- genSeasonList(hourtemps, mrange = c(8, 6), years = cal_data_list[[i]][[8]]$Year)
  eval_SeasonList_Collaos <- genSeasonList(hourtemps, mrange = c(8, 6), years = val_data_list[[i]][[8]]$Year)
  
  SeasonList_Coloradona <- genSeasonList(hourtemps, mrange = c(8, 6), years = cal_data_list[[i]][[9]]$Year)
  eval_SeasonList_Coloradona <- genSeasonList(hourtemps, mrange = c(8, 6), years = val_data_list[[i]][[9]]$Year)
  
  SeasonList_DelaRiega <- genSeasonList(hourtemps, mrange = c(8, 6), years = cal_data_list[[i]][[10]]$Year)
  eval_SeasonList_DelaRiega <- genSeasonList(hourtemps, mrange = c(8, 6), years = val_data_list[[i]][[10]]$Year)
  
  SeasonList_Perezosa <- genSeasonList(hourtemps, mrange = c(8, 6), years = cal_data_list[[i]][[11]]$Year)
  eval_SeasonList_Perezosa <- genSeasonList(hourtemps, mrange = c(8, 6), years = val_data_list[[i]][[11]]$Year)
  
  SeasonList_Perico <- genSeasonList(hourtemps, mrange = c(8, 6), years = cal_data_list[[i]][[12]]$Year)
  eval_SeasonList_Perico <- genSeasonList(hourtemps, mrange = c(8, 6), years = val_data_list[[i]][[12]]$Year)
  
  SeasonList_Raxao <- genSeasonList(hourtemps, mrange = c(8, 6), years = cal_data_list[[i]][[13]]$Year)
  eval_SeasonList_Raxao <- genSeasonList(hourtemps, mrange = c(8, 6), years = val_data_list[[i]][[13]]$Year)
  
  SeasonList_Teorica <- genSeasonList(hourtemps, mrange = c(8, 6), years = cal_data_list[[i]][[14]]$Year)
  eval_SeasonList_Teorica <- genSeasonList(hourtemps, mrange = c(8, 6), years = val_data_list[[i]][[14]]$Year)
  
  SeasonList_Verdialona <- genSeasonList(hourtemps, mrange = c(8, 6), years = cal_data_list[[i]][[15]]$Year)
  eval_SeasonList_Verdialona <- genSeasonList(hourtemps, mrange = c(8, 6), years = val_data_list[[i]][[15]]$Year)
  
  SeasonList_Xuanina <- genSeasonList(hourtemps, mrange = c(8, 6), years = cal_data_list[[i]][[16]]$Year)
  eval_SeasonList_Xuanina <- genSeasonList(hourtemps, mrange = c(8, 6), years = val_data_list[[i]][[16]]$Year)
  
  
  #specify number of cultivars
  n_cult <- 16
  cultivar_names <- c('Berlepsch', 'Cox Orange', 'Golden Delicious', 'James Grieve', 'Roter Boskoop',
                      'Blanquina', 'Clara', 'Collaos', 'Coloradona', 'DelaRiega','Perezosa','Perico',
                      'Raxao','Teorica','Verdialona', 'Xuanina')
  
  #season list made for each cultivar need to be combined in a giant list
  SeasonList <- list(SeasonList_berlepsch,SeasonList_cox_orange,SeasonList_golden_delicious, 
                     SeasonList_james_grieve, SeasonList_roter_boskoop, SeasonList_Blanquina, SeasonList_Clara, 
                     SeasonList_Collaos,SeasonList_Coloradona,SeasonList_DelaRiega,
                     SeasonList_Perezosa, SeasonList_Perico, SeasonList_Raxao,SeasonList_Teorica,SeasonList_Verdialona,
                     SeasonList_Xuanina)
  
  eval_SeasonList <- list(eval_SeasonList_berlepsch,eval_SeasonList_cox_orange,eval_SeasonList_golden_delicious, 
                          eval_SeasonList_james_grieve, eval_SeasonList_roter_boskoop, eval_SeasonList_Blanquina, 
                          eval_SeasonList_Clara, eval_SeasonList_Collaos,eval_SeasonList_Coloradona,
                          eval_SeasonList_DelaRiega,eval_SeasonList_Perezosa, eval_SeasonList_Perico, eval_SeasonList_Raxao,
                          eval_SeasonList_Teorica,eval_SeasonList_Verdialona,eval_SeasonList_Xuanina)
  
  #same for phenology observations used for model training
  bloomJDays_list <- list( cal_data_list[[i]][[1]]$Bloom_full_JDay,cal_data_list[[i]][[2]]$Bloom_full_JDay,
                           cal_data_list[[i]][[3]]$Bloom_full_JDay,cal_data_list[[i]][[4]]$Bloom_full_JDay,
                           cal_data_list[[i]][[5]]$Bloom_full_JDay, cal_data_list[[i]][[6]]$pheno,
                           cal_data_list[[i]][[7]]$pheno,
                           cal_data_list[[i]][[8]]$pheno,cal_data_list[[i]][[9]]$pheno,
                           cal_data_list[[i]][[10]]$pheno,cal_data_list[[i]][[11]]$pheno,
                           cal_data_list[[i]][[12]]$pheno,cal_data_list[[i]][[13]]$pheno,
                           cal_data_list[[i]][[14]]$pheno,cal_data_list[[i]][[15]]$pheno,
                           cal_data_list[[i]][[16]]$pheno)
  
  bloomJDays_eval_list <- list( val_data_list[[i]][[1]]$Bloom_full_JDay,val_data_list[[i]][[2]]$Bloom_full_JDay,
                                val_data_list[[i]][[3]]$Bloom_full_JDay,val_data_list[[i]][[4]]$Bloom_full_JDay,
                                val_data_list[[i]][[5]]$Bloom_full_JDay, val_data_list[[i]][[6]]$pheno,
                                val_data_list[[i]][[7]]$pheno,
                                val_data_list[[i]][[8]]$pheno,val_data_list[[i]][[9]]$pheno,
                                val_data_list[[i]][[10]]$pheno,val_data_list[[i]][[11]]$pheno,
                                val_data_list[[i]][[12]]$pheno,val_data_list[[i]][[13]]$pheno,
                                val_data_list[[i]][[14]]$pheno,val_data_list[[i]][[15]]$pheno,
                                val_data_list[[i]][[16]]$pheno)
  
  SeasonList_list[[i]] <- SeasonList
  eval_SeasonList_list[[i]] <- eval_SeasonList
  bloomJDays_list_squared[[i]] <- bloomJDays_list
  bloomJDays_eval_list_squared[[i]] <- bloomJDays_eval_list
  
  #values for plotting temperature response
  temp_values = seq(-5, 40, 0.1)
  
  #      ('yc', 'zc', 's1', 'Tu', 'theta*', 'theta_c', 'Tau(thetha*)', 'pie_c',   'Tf', 'Tc', 'Tb',  'slope')
  x_0 <- c(rep(40, n_cult), rep(190, n_cult), rep(0.5, n_cult),  25,   278.1,    285.1,   47.7,           30.6,        4,   36,     4,    1.60)
  x_U <- c(rep(80, n_cult), rep(500, n_cult), rep(1.0, n_cult),  30,   281,      287,       48,             50,       10,   55,    10,    5.00)
  x_L <- c(rep(20, n_cult), rep(100, n_cult), rep(0.1, n_cult),   0,   276,      284,       16,             24,        2,    0,     2,    1)
  
  #limits for the inequality constraints
  #these are used to make sure you don't get completely nonsense parameter combinations
  #like upper temperature threshold lower than the lower temperature threshold and so on
  #         #gdh parameters   #q10 for E0 and E1
  c_L <- c(  0,   0,   0,     1.5, 1.5)
  c_U <- c(Inf, Inf, Inf,     3.5, 3.5)
  
  #this list is needed by the optimizer, the function f
  #is the function the optimizer inserts the parameters into
  #it is defined in the LarsChill package. you can also custimize it (like I did)
  #or you can run it for the x_0 parameters if you are curious and want to see what it returns
  
  problem<-list(f="evaluation_function_meigo_nonlinear_combined",
                x_0 = x_0,
                x_L = x_L,
                x_U = x_U,
                c_L = c_L, 
                c_U = c_U)
  opts<-list(maxeval = 30000,
             #maxtime = 10 * 1, 
             local_solver = 'DHC', #
             local_bestx = 1,
             inter_save = 0)
  
  #####round 1######
  #set.seed(123456789)
  #this is the goal of the whole procedure, here we start the actual fitting of PhenoFlex
  res_loop[[i]]<-custom_essr(problem, #this was defined earlier, contains the function, the initial guess and the search space
                             opts, #different hyperparameters for the global optimization
                             
                             #the other arguments are not needed by meigo, but our evaluation function, which we specified in the problems list
                             modelfn = custom_PhenoFlex_GDHwrapper, #this is the function to calculate the bloom date of an individual year
                             bloomJDays = bloomJDays_list, #vector containing the bloom date in day of the year format
                             SeasonList = SeasonList,
                             n_cult = n_cult) #list of hourly temperatures, generated by genseasonlist and stackhourlytemperatures functions from chillr
  #                                                                                    needs to be of the same order as bloomdays
  
}





par_list <- list()
for(rep in 1:length(res_loop)){
  par_list[[rep]] <- list()
  for(i in 1:n_cult){
    par_list[[rep]][[i]] <- res_loop[[rep]]$xbest[c(i, n_cult + i, (n_cult*2)+i,
                                                    (length(res_loop[[rep]]$xbest)-8):length(res_loop[[rep]]$xbest))]
  }
}


prediction_df <- data.frame()


for(rep in 1:length(res_loop)){
  for(i in 1:n_cult){
    
    #calibration
    pheno_cal_pred <- return_predicted_days(convert_parameters(par_list[[rep]][[i]]), 
                                            modelfn = custom_PhenoFlex_GDHwrapper, 
                                            SeasonList = SeasonList_list[[rep]][[i]])
    
    cal_obs <- bloomJDays_list_squared[[rep]][[i]]
    
    
    prediction_df <- rbind(prediction_df,
                           data.frame(predicted = pheno_cal_pred,
                                      observed = cal_obs,
                                      cultivar = cultivar_names[i], 
                                      split = 'Calibration',
                                      round = rep)
    )
    
    
    
    pheno_val_pred <- return_predicted_days(convert_parameters(par_list[[rep]][[i]]), 
                                            modelfn = custom_PhenoFlex_GDHwrapper, 
                                            SeasonList = eval_SeasonList_list[[rep]][[i]])
    
    
    
    eval_obs <- bloomJDays_eval_list_squared[[rep]][[i]]
    
    prediction_df <- rbind(prediction_df,
                           data.frame(predicted = pheno_val_pred,
                                      observed = eval_obs,
                                      cultivar = cultivar_names[i], 
                                      split = 'Validation',
                                      round = rep)
    )
    
  }
}


library(tidyverse)

iqr_df <- prediction_df %>% 
  group_by(cultivar) %>% 
  summarise(iqr = IQR(observed))

iqr_all <- IQR(prediction_df$observed)

#calculate performance on cultivar level
performance_df <- prediction_df %>%
  group_by(cultivar, split) %>% 
  summarise(rmse = RMSEP(predicted, observed)) %>% 
  merge(iqr_df, by = c('cultivar')) %>% 
  mutate(rpiq = iqr / rmse)


#calculate performance across cultivars
performance_df_all <- prediction_df %>%
  group_by(split) %>% 
  summarise(rmse = RMSEP(predicted, observed)) %>% 
  mutate(rpiq = iqr_all / rmse)

#save_fitting_list(list("Germany_Spain" = res_loop),path ="enhanced_global_opt/de_es_fit/", prefix = "fit_list" )

names(res_loop) <- paste0('rep', 1:10)
save_fitting_list( res_loop,path ="enhanced_global_opt/de_es_fit/", prefix = "test" )



#get temperature response
for (i in 1:10) {
  p <- get_temp_response_plot(convert_parameters(par_list[[i]][[1]]), temp_values)
  ggsave(filename = paste0("enhanced_global_opt/de_es_fit/temp_response_run_",i, ".jpeg"),
         plot = p,
         device = "jpeg",
         width = 17,
         height = 12,
         units = "cm",
         path = getwd())
}

p <- get_temp_response_plot(convert_parameters(par_list[[2]][[1]]), temp_values)

ggsave(filename = "enhanced_global_opt/de_es_fit/temp_response_run_2.jpeg",
       plot = p,
       device = "jpeg",
       width = 17,
       height = 12,
       units = "cm",
       path = getwd())
















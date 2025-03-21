#this function takes the weather data (x) and the model parameters (par) and returns the expected day of bloom 
custom_PhenoFlex_GDHwrapper <- function (x, par, log_par = FALSE){
  #x is one of the elements in season List
  #par are the parameters of the model
  
  if(log_par){
    #transfrom log to normla value of parameters
    par[7]<-exp(par[7])
    par[8]<-exp(par[8])
  }

  
  #in case the parameters do not make sense, return NA
  # if (par[4] <= par[11]){
  #   return(NA)
  # } else if(par[10] <= par[4]){
  #   return(NA)
  # } else if(exp((10 * par[5]) / (297 * 279)) < 1.5 |  exp((10 * par[5]) / (297 * 279)) > 3.5 ){
  #   #also when the q10 criterion is outside the limits of 1.5 and 3.5
  #   return(NA)
  # } else if(exp((10 * par[6]) / (297 * 279)) < 1.5 |  exp((10 * par[6]) / (297 * 279)) > 3.5 ){
  #   return(NA)
  # } else{
    #calculat the bloom day
    bloomindex <- PhenoFlex(temp = x$Temp, times = seq_along(x$Temp), 
                            yc = par[1], zc = par[2], s1 = par[3], Tu = par[4], E0 = par[5], 
                            E1 = par[6], A0 = par[7], A1 = par[8], Tf = par[9], Tc = par[10], 
                            Tb = par[11], slope = par[12], Imodel = 0L, basic_output = TRUE)$bloomindex
    
    #return values
    if (bloomindex == 0){
      return(NA)
    } 
    
    JDay <- x$JDay[bloomindex]
    JDaylist <- which(x$JDay == JDay)
    n <- length(JDaylist)
    if (n == 1){
      return(JDay)
    } 
    return(JDay + which(JDaylist == bloomindex)/n - 1/(n/ceiling(n/2)))
  #}
  
}


evaluation_function_meigo <- function(x, 
                                      modelfn,
                                      bloomJDays,
                                      SeasonList,
                                      na_penalty = 365,
                                      return_bloom_days = FALSE){
  
  #innput:
  #         x is the parameters in meigo
  #         modelfn is the function used to calculate the bloomdays
  #         SeasonList contains the weather data for the different years
  #         na_penalty is the value used if the model failed to predict any day with the current set of parameters for a particular year
  
  #output: inequality constraints g
  #        model performance value F
  
  #change the name of the parameters, dont know if necessary
  par <- x
  
  #calculate the predicted flower dates
  pred_bloom <- unlist(lapply(X = SeasonList, FUN = modelfn, par = par))
  
  #if the model returns no bloom day, then give penalty term instead
  pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, no = pred_bloom)
  
  #calculate the model performance value
  F <- sum((pred_bloom - bloomJDays)^2)
  
  
  
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
  g[1] <- x[4] - x[11]
  #Tx >= Tb
  g[2] <- x[10] - x[11]
  #Tc >= Tu
  g[3] <- x[10] - x[4]
  
  
  #the q10 value should be between 1.5 and 3.5
  #q10 formation = exp((10*E0)/(T2-T1))
  #q10 destruction = exp((10*E1)/(T2-T1))
  #T1 = 297
  #T2 = 279
  #parameters T1 and T2 chosen after Egea 2020
  #ranges for q10 are provided in c_L and c_U
  g[4] <- exp((10 * x[5]) / (297 * 279))
  
  g[5] <- exp((10 * x[6]) / (297 * 279))
  
  
  if(return_bloom_days == FALSE){
    #output
    return(list(F=F, g=g))
  } else{
    return(pred_bloom)
  }
  
}

evaluation_function_meigo_vns <- function(x, 
                                          modelfn,
                                          bloomJDays,
                                          SeasonList,
                                          na_penalty = 365,
                                          return_bloom_days = FALSE){
  
  #innput:
  #         x is the parameters in meigo
  #         modelfn is the function used to calculate the bloomdays
  #         SeasonList contains the weather data for the different years
  #         na_penalty is the value used if the model failed to predict any day with the current set of parameters for a particular year
  
  #output: inequality constraints g
  #        model performance value F
  
  #change the name of the parameters, dont know if necessary
  par <- x
  
  #calculate the predicted flower dates
  pred_bloom <- unlist(lapply(X = SeasonList, FUN = modelfn, par = par))
  
  #if the model returns no bloom day, then give penalty term instead
  pred_bloom <- ifelse(is.na(pred_bloom), yes = na_penalty, no = pred_bloom)
  
  #calculate the model performance value
  F <- sum((pred_bloom - bloomJDays)^2)
  
  
  
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
  g[1] <- x[4] - x[11]
  #Tx >= Tb
  g[2] <- x[10] - x[11]
  #Tc >= Tu
  g[3] <- x[10] - x[4]
  
  
  #the q10 value should be between 1.5 and 3.5
  #q10 formation = exp((10*E0)/(T2-T1))
  #q10 destruction = exp((10*E1)/(T2-T1))
  #T1 = 297
  #T2 = 279
  #parameters T1 and T2 chosen after Egea 2020
  #ranges for q10 are provided in c_L and c_U
  g[4] <- exp((10 * x[5]) / (297 * 279))
  
  g[5] <- exp((10 * x[6]) / (297 * 279))
  
  
  if(return_bloom_days == FALSE){
    #output
    return(F)
  } else{
    return(pred_bloom)
  }
  
}
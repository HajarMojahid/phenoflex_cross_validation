#contains three functions: (1) ensemble_predictions_with_failure: makes predictions from parameters, returns weighted mean
#                          (2) weighted_mean_with_fail:           takes predictions of repetitions, returns weighted mean
#                          (3) return_weights:                    helper function for the functions above
#
#biggest difference between 1) and 2): (1) makes predictions from scratch, (2) works with already made predictions (in case you saved them somewhere and don't want to redo them)
#                                      both yield the same output, they just have different starting points

# functions (1) and (2) make very specific assumptions on the input, you may wanna adjust that to your own convenience, 
# the input of parameters and confidence are a bit different between (1) and (2)



#input:
#     par_list:   list of the parameters 
#                   entries contain one read of the read with the function LarsChill::load_fitting_result(), each entry is a repetition of the same cultivar
#                   parameters are extracted from the 'xbest' entry within the individual list members
#                   I assumed ten model parameters with theta_star and Tc being fixed, but you can change that in the code after the line par <- x$xbest
#                   order of entries is identical with order in confidence (see next input)
#     confidence: gives the weights to the individual predictions
#                   numeric vector of same length as the par_list
#                   assumes that bigger is better
#     temp:       seasonlist containing hourly temperature data for the predictions
#                   generated with chillR::genSeasonList
#     theta_star: fixed value for theta_star
#     Tc:         fixed value for Tc
#     return_se:  logical, decides if standard deviation of the predictions around the weighted mean should be returned as well
#     n_fail:     numeric, decides the cut-off number of failure predictions of the weighted mean members so that the weighted mean also returns a failure
#                  in case n_fail = 5: if we have 4 or less failure predictions --> get ignored and the weighted mean is calculated based on remaining results
#                                      if we have 5 or more failure predictions --> weighted mean is failure as well
#
#output:
#     when return_se = TRUE --> list with weighted mean, sd and individual model predictions
#     when return_se = FALSE --> vector with weighted means


ensemble_prediction_with_failure <- function(par_list, confidence, temp, theta_star = 279, Tc = 36, return_se = TRUE, n_fail = 5,
                                             max_weight = NULL){
  
  
  predicted <- purrr::map(par_list, function(x){
    par <- x$xbest
    par <- c(par[1:4], theta_star, par[5:8], Tc, par[9:10])
    
    return_predicted_days(convert_parameters(par), 
                          modelfn = custom_PhenoFlex_GDHwrapper, 
                          SeasonList =temp,
                          na_penalty = NA)
  }) %>% 
    stats::setNames(1:length(par_list)) %>% 
    dplyr::bind_cols() %>% 
    as.matrix()
  
  #if as many as n_fail model runs indicate failure, then mark as NA, otherwise discard them and calculate mean of the remaining ones
  pred_na <- apply(predicted, MARGIN = 1, FUN = function(x) sum(is.na(x)) >= n_fail) 
  
  
  #create a weight data.frame with the same dimensions as predicted
  weights <- return_weights(predicted = predicted,
                            confidence = confidence, 
                            n_fail = n_fail, 
                            max_weight = max_weight)
  
  #calculate sd, excluding NA values
  sd_pred <- predicted %>% 
    t() %>% 
    as.data.frame() %>% 
    magrittr::set_colnames( 1:nrow(predicted)) %>% 
    reshape2::melt(id.vars = NULL) %>% 
    dplyr::group_by(.data$variable) %>% 
    dplyr::summarise(sd = stats::sd(.data$value, na.rm = TRUE)) %>% 
    dplyr::pull(.data$sd)
  
  #replace predicted NA with 0s
  predicted_tmp <- predicted %>% replace(is.na(.), 0)
  
  #calculate weighted individual pred, then get the sum
  weighted_pred <- predicted_tmp * weights 
  pred_out <- rowSums(weighted_pred)
  
  #in case too many models indicated NA, then the whole prediction becomes NA
  pred_out[pred_na] <- NA
  
  
  
  if(return_se){
    return(list(predicted = pred_out, sd = sd_pred,
                individual_pred = predicted))
  } else{
    return(weighted_pred)
  }
  
}


#input: 
#1) predicted: data.frame with individual predictions, need to contain column called id
#               role of id column is to match predictions with the entries of the confidence input
#               so each cultivar would get their own id, can be the name, or something else
#              assume wide format --> so one column for each prediction
#               that means for one particular year we have ten columns with the ten predictions for that year (assuming you have ten repetitions)
#               function assumes that the names of the columns containing the predictions will be also present in the confidence input
#               so if the predictions are stored in columns R1 to R10, these column names need to be present in confidence, too. 
#              
#2) confidence: data.frame, needs to contain same columns as predicted (through predicted can have many more additional columns)
#               assume wide format, so 10 columns when ten repitions, plus id
#               assumes that larger values mean more confidence
#               id entries need to match with id entries of the first input ('predicted')
#3) theta_star: fixed value for theta_star
#4) Tc:         fixed value for Tc
#5) return_se:  logical, decides if standard deviation of the predictions around the weighted mean should be returned as well
#6) n_fail:     numeric, decides the cut-off number of failure predictions of the weighted mean members so that the weighted mean also returns a failure
#                  in case n_fail = 5: if we have 4 or less failure predictions --> get ignored and the weighted mean is calculated based on remaining results
#                                      if we have 5 or more failure predictions --> weighted mean is failure as well
#output:
# same as 'predicted' but with additional column containing weighted mean and standard deviation
#


weighted_mean_with_fail <- function(predicted, confidence, return_se = TRUE, n_fail = 5, max_weight = NULL, .progress = FALSE){
  
  #check if the required colnames are present and that they match
  stopifnot('id' %in% colnames(confidence))
  stopifnot('id' %in% colnames(predicted))
  
  other_colnames <- colnames(confidence)[colnames(confidence) != 'id']
  stopifnot(all(other_colnames %in% colnames(predicted)))
  
  #save the old predicted
  predicted_old <- predicted
  confidence_old <- confidence
  
  #iterate other the different ids
  pred_out <- purrr::map(unique(predicted_old$id), function(i){
    
  #i <- unique(predicted_old$id)[2]
    
    #subset to the columns we need
    predicted <- predicted_old %>% 
      filter(id == i) %>% 
      select(all_of(other_colnames)) %>% 
      as.matrix()
    
    confidence <- confidence_old %>% 
      filter(id == i) %>% 
      select(all_of(other_colnames)) %>% 
      unlist() %>% 
      unname()
    
    #if as many as n_fail model runs indicate failure, then mark as NA, otherwise discard them and calculate mean of the remaining ones
    pred_na <- apply(predicted, MARGIN = 1, FUN = function(x) sum(is.na(x)) >= n_fail) 
    
    
    #create a weight data.frame with the same dimensions as predicted
    weights <- return_weights(predicted = predicted,
                              confidence = confidence, 
                              n_fail = n_fail, 
                              max_weight = max_weight)
    
    #calculate sd, excluding NA values
    sd_pred <- predicted %>% 
      t() %>% 
      as.data.frame() %>% 
      magrittr::set_colnames( 1:nrow(predicted)) %>% 
      reshape2::melt(id.vars = NULL) %>% 
      dplyr::group_by(.data$variable) %>% 
      dplyr::summarise(sd = stats::sd(.data$value, na.rm = TRUE)) %>% 
      dplyr::pull(.data$sd)
    
    #replace predicted NA with 0s
    predicted_tmp <- predicted %>% replace(is.na(.), 0)
    
    #calculate weighted individual pred, then get the sum
    weighted_pred <- predicted_tmp * weights 
    pred_out <- rowSums(weighted_pred)
    
    #in case too many models indicated NA, then the whole prediction becomes NA
    pred_out[pred_na] <- NA
    
    predicted_old %>% 
      filter(id == i) %>% 
      mutate(weighted_pred = pred_out,
             sd = sd_pred) %>% 
      return()
  }, .progress = .progress) %>% 
    bind_rows()
  
  if(return_se == FALSE){
    pred_out <- pred_out %>% 
      select(-sd)
  }
  
  return(pred_out)
}




return_weights <- function(predicted, confidence, n_fail = 5, max_weight = NULL){
  
  apply(predicted, MARGIN = 1, function(x){
    
    #create intermediate object which I can manipulate
    conf <- confidence
    
    if(sum(is.na(x)) != 0){
      #position of NA
      pos.na <- which(is.na(x))
      
      #calculate replace the confidence at position of NA with zero
      conf[pos.na] <- 0
    }
    
    
    conf_out <- rep(0, length(conf))
    #calculate the relative weight of confidence
    if(sum(conf) !=0) conf_out <- conf / sum(conf)
    
    #in case we restrict the maximum weight
    if(is.null(max_weight) == FALSE & all(conf == 0) == FALSE){
      
      #position where confidence exceeds max value
      pos_toohigh <- conf_out > max_weight
      
      #check if the remaining entries have any confidences (are not NA)
      #if so, skip the routine and keep weights as they are
      if(sum(conf_out[pos_toohigh == FALSE]) != 0){
        
        #calculate how much confidence to allocate to remaining entries
        exceed_conf <- ifelse(conf_out > max_weight, yes = conf_out - max_weight, no = 0) %>% sum()
        
        
        #distribute exceed confidence
        weight_rest_relative <- conf_out / sum(conf_out[pos_toohigh == FALSE])
        weight_rest_relative[pos_toohigh] <- 0
        
        extra_conf <- weight_rest_relative * exceed_conf
        conf_out[pos_toohigh] <- max_weight
        conf_out <- conf_out + extra_conf
      }
      
    }
    #make sure the sum is always 1 (usual case), 0 (when all predictions are NA or max_weight)
    stopifnot(round(sum(conf_out), digits = 2) %in% c(1, 0, max_weight))
    
    return(conf_out)
  }) %>% 
    t() %>% 
    return()
}

PhenoFlex_custom_wrapper <- function(x, par,
                                     A1 = 5.94e+13, #activation energy destroying PDBF
                                     A0 = 6.34e+3,
                                     E0 = 3.37e+3, #optimal temperature PDBF formation
                                     E1 = 9.91e+3, #optimal temperature PDBF destruction
                                     slope = 1.08, #slope of sigmoidal function to permanently store PDBF as CP
                                     Tf = 5.38, #transition parameter of sigmpodial function to permanently store PDBF as CP
                                     Tu = 2.57e+1, #optimal temperature for forcing
                                     Tb = 5.23, #base temperature for forcing
                                     Tc =27.7) { #upper temperature limit for forcing
  ## for GDH Tb > Tu and Tu > Tc
  #if((par[4] <= par[6]) || (par[5] <= par[4])) return(NA)
  bloomindex <- PhenoFlex(temp=x$Temp,
                          times=seq_along(x$Temp),
                          yc=par[1],
                          zc=par[2],
                          s1=par[3],
                          Tu=Tu,
                          E0=E0,
                          E1=E1,
                          A0=A0,
                          A1=A1,
                          Tf=Tf,
                          Tc=Tc,
                          Tb=Tb,
                          slope=slope,
                          Imodel=0L,
                          basic_output=TRUE)$bloomindex
  if(bloomindex == 0) return(NA)
  ## attempt to smooth the function a bit
  JDay <- x$JDay[bloomindex]
  JDaylist <- which(x$JDay == JDay)
  n <- length(JDaylist)
  if(n == 1) return(JDay)
  return(JDay + which(JDaylist == bloomindex)/n - 1./(n/ceiling(n/2)))
}

predictBloomDays <- function(par, SeasonList, modelfn, ...) {
  paravail <- requireNamespace("parallel")
  sres <- 
    if(FALSE) { ## currently, parallel seems to be slower
      mc.cores <- as.numeric(Sys.getenv("OMP_NUM_THREADS"))
      if(is.na(mc.cores)) mc.cores <- getOption("mc.cores", default=1L)
      parallel::mclapply(X=SeasonList, FUN=modelfn, par=par, ...,
                         mc.cores = mc.cores)
    }
  else {
    lapply(X=SeasonList, FUN=modelfn, par=par)
  }
  return(invisible(simplify2array(x=sres)))
}



#custom phenology fitter and custom chifull function
phenologyFitter_custom <- function(par.guess=NULL,
                                   modelfn=PhenoFlex_GDHwrapper,
                                   bloomJDays,
                                   SeasonList,
                                   control=list(smooth=FALSE, verbose=TRUE, maxit=1000,
                                                nb.stop.improvement=250),
                                   lower,
                                   upper,
                                   seed = 1235433,
                                   name_par = NULL,
                                   par.fixed = NULL,
                                   name_par.fixed = NULL,
                                   ...) {
  
  control$seed <- seed
  
  ## some sanity checks
  stopifnot(is.list(SeasonList))
  stopifnot(length(SeasonList) == length(bloomJDays))
  
  #if par.fixed is supplied, then name_par and name_par.fixed must be supplied as well
  if(is.null(par.fixed) == F){
    if(is.null(name_par) | is.null(name_par.fixed)){
      stop('If par.fixed is supplied, then name_par and name_par.fixed must be supplied too')
    }
    
    #the names in name_par and name_par.fixed need to contain only the following
    stopifnot(sum(c(name_par, name_par.fixed) %in% c('yc', 'zc', 's1', 'Tu', 'E0', 'E1', 'A0', 'A1', 'Tf', 'Tc', 'Tb',  'slope')) == 12)
  }
  
  
  
  ## generate empty `phenologyFit` and store input objects
  res <- phenologyFit()
  res$par.guess <- par.guess
  res$modelfn  <- modelfn
  res$bloomJDays <- bloomJDays
  res$SeasonList <- SeasonList
  res$lower <- lower
  res$upper <- upper
  res$control <- control
  res$par.fixed <- par.fixed
  res$name_par.fixed <- name_par.fixed
  res$name_par <- name_par
  ##res$... <- ...
  
  ## cure the 365 to 1 jump at new years in JDays and, possibly, bloomJDays
  for(i in c(1:length(SeasonList))) {
    minJDay <- SeasonList[[i]]$JDay[1]
    maxJDay <- SeasonList[[i]]$JDay[length(SeasonList[[i]]$JDay)]
    ## some more sanity checks
    if(maxJDay > minJDay) {
      stop(paste0("Season ", i, " is overlapping with the previous or following one. Aborting!"))
    }
    if(bloomJDays[i] > maxJDay && bloomJDays[i] < minJDay) {
      stop(paste0("In season ", i, " the bloomJDay is outside the provided JDay vector. Aborting!"))
    }
    ## determine the index of the jump and its value
    dx <- diff(SeasonList[[i]]$JDay)
    mx <- min(dx)
    kmx <- which(dx == mx)
    SeasonList[[i]]$JDay[1:kmx] <- SeasonList[[i]]$JDay[1:kmx] + mx - 1
    if(bloomJDays[i] > minJDay) bloomJDays[i]  <- bloomJDays[i] + mx - 1
    ## store this also in res for later usage
    res$SeasonList[[i]]$JDayunwrapped <- SeasonList[[i]]$JDay
    res$bloomJDaysunwrapped <- bloomJDays
  }
  
  ## now we fit
  res$model_fit <- GenSA::GenSA(par=par.guess, fn=chifull_custom, bloomJDays=bloomJDays,
                                SeasonList=SeasonList, 
                                modelfn=modelfn,
                                control=control,
                                lower=lower, upper=upper,
                                par.fixed = par.fixed, 
                                name_par = name_par, 
                                name_par.fixed = name_par.fixed)
  res$par <- res$model_fit$par
  
  if(is.null(par.fixed) == F & is.null(name_par.fixed) == F & is.null(name_par) == F){
    #append par and par.fixed
    test <- c(res$par, par.fixed)
    test_name <- c(name_par, name_par.fixed)
    
    correct_order <- purrr::map_dbl(c('yc', 'zc', 's1', 'Tu', 'E0', 'E1', 'A0', 'A1', 'Tf', 'Tc', 'Tb',  'slope'), .f =  function(x) which(x == test_name))
    res$par_full <- test[correct_order]
    
    res$pbloomJDays <- predictBloomDays(par=res$par_full,
                                       SeasonList=res$SeasonList,
                                       modelfn=res$modelfn)
  } else{
    res$pbloomJDays <- predictBloomDays(par=res$par,
                                        SeasonList=res$SeasonList,
                                        modelfn=res$modelfn)
  }
  

  return(res)
}


chifull_custom <- function(par,
                           modelfn,
                           bloomJDays,
                           SeasonList,
                           na_penalty=365,
                           name_par,
                           par.fixed,
                           name_par.fixed,
                           ...) {
  # #---------#
  # #temp
  # #---------#
  # par <-   c(40, 190, 0.5, 5.939917e13)
  # name_par <- c('yc', 'zc', 's1', 'A1')
  # par.fixed <- c(6.34e+3,3.37e+3, 9.91e+3, 1.08, 5.38, 2.57, 5.23, 27.7)
  # name_par.fixed <- c('A0', 'E0', 'E1', 'slope', 'Tf', 'Tu', 'Tb', 'Tc')
  # #---------#
  # #temp
  # #---------#
  
  #if fixed parameter are supplied, bring par in correct order
  if(is.null(par.fixed) == F & is.null(name_par.fixed) == F & is.null(name_par) == F){
    #append par and par.fixed
    test <- c(par, par.fixed)
    test_name <- c(name_par, name_par.fixed)
    
    correct_order <- purrr::map_dbl(c('yc', 'zc', 's1', 'Tu', 'E0', 'E1', 'A0', 'A1', 'Tf', 'Tc', 'Tb',  'slope'), .f =  function(x) which(x == test_name))
    par <- test[correct_order]
  }
  
  
  
  
  
  
  sres <- predictBloomDays(par=par, SeasonList=SeasonList,
                           modelfn=modelfn,...)
  s <- (sres-bloomJDays)
  ## we penalise if no bloomday could be found
  nai <- which(is.na(sres))
  s[nai] <- na_penalty
  ## Remaining NAs can only come from the data, not from the predictions
  return(sum(s[!is.na(s)]^2))
}

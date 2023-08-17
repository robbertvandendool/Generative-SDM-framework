###
### Performance benchmark
###

permutevarFinal_maxnet = function(presdata, backdata, maxnetmodel, variable, repeats){
  output = data.frame("id" = seq_len(repeats), 
                      "NLLp" = rep(0, repeats), 
                      "AUC" = rep(0, repeats))
  
  combdata = rbind(presdata, backdata)
  combdata$label = c(rep(1,nrow(presdata)),rep(0,nrow(backdata)))
  
  results = furrr::future_map(seq_len(repeats),.options = furrr::furrr_options(seed = TRUE),function(x){
    combdata[,variable] = sample(combdata[,variable])
    
    prob = maxnet:::predict.maxnet(maxnetmodel, combdata, clamp=T, type=c("exponential"))
    psummed = sum(prob)
    presprob = prob[combdata$label==1] 
    presprobscaled = presprob / psummed
    probscaled = prob/psummed
    
    logloss = -sum(log(presprobscaled))
    auc = pROC::auc(response = combdata$label, predictor = probscaled, quiet=TRUE, na.rm=TRUE) 
    c(logloss, auc)
    
  } )
  
  output[,2:3] = do.call(rbind.data.frame, results)
  return(output)
}

mostimportantvars = function(presdataset, backdataset, maxnetmod){
  VarImp = future_map(names(backdataset),.options = furrr::furrr_options(seed = TRUE),function(x) permutevarFinal_maxnet(presdata = presdataset, backdata = backdataset, maxnetmodel=maxnetmod, variable = x,repeats = 50))
  names(VarImp) = names(backdataset)
  
  VarImp_AUC = lapply(VarImp, function(x) myCI(1 - x$AUC)) #how much is AUC higher with the variable not shuffled?
  VarImp_AUC_min = min(sapply(VarImp_AUC, function(x)x[2]),na.rm = T)
  VarImp_AUC_max = max(sapply(VarImp_AUC, function(x)x[2]),na.rm = T)
  VarImp_AUC_relative = lapply(VarImp_AUC , function(y) (y - VarImp_AUC_min)/ (VarImp_AUC_max-VarImp_AUC_min))
  VarImp_AUC_relativemean = sapply(VarImp_AUC_relative, function(x)x[2])
  names(VarImp_AUC_relativemean) = names(VarImp)
  VarImp_AUC_relativesort = sort(VarImp_AUC_relativemean, decreasing=T)
  
  vars = names(VarImp_AUC_relativesort)[1:4]
  
  return(vars)
}


fitmodels_species_all = function(PO, back, PA, env, species, varnames){
  
  tryCatch({
    
    PO=PO
    back=back
    PA=PA
    env=env
    species=species
    varnames=varnames
    
    padat = c(PO[PO$spid==species,which(names(PO)=="occ")], back[,which(names(back)=="occ")])
    envdat = rbind(PO[PO$spid==species,varnames], back[,varnames])
    tpadat = PA[,species]
    tenvdat = env[,varnames]
    
    mins <- sapply(varnames, function(x){min(envdat[,x])})
    maxs <- sapply(varnames, function(x){max(envdat[,x])})
    
    bdata <- as.matrix(back[,varnames])
    pdata <- as.matrix(PO[PO$spid==species,varnames])
    
    #A) vine copula using all variables
    start_time <- Sys.time()
    vine_back <- vine(bdata, 
                      margins_controls = list(mult = 1, xmin = mins, xmax = maxs, bw = NA, deg = 2), #mult=NULL means  log(1 + d) after rvinecopulib:::expand_factors(). This is only true for discrete variables. 
                      copula_controls = list(family_set = "nonparametric", structure = NA, par_method = "mle",
                                             nonpar_method = "quadratic", mult = 1, selcrit = "mbicv", psi0 = 0.9, presel = TRUE,
                                             trunc_lvl = NA, tree_crit = "tau", threshold = 0, keep_data = FALSE, show_trace =
                                               FALSE, cores = 6),
                      cores=6)
    end_time <- Sys.time() 
    end_time - start_time 
    
    vine_pres <- vine(pdata, 
                      margins_controls = list(mult = 1, xmin = mins, xmax = maxs, bw = NA, deg = 2), #mult=NULL means  log(1 + d) after rvinecopulib:::expand_factors(). This is only true for discrete variables. 
                      copula_controls = list(family_set = "nonparametric", structure = NA, par_method = "mle",
                                             nonpar_method = "quadratic", mult = 1, selcrit = "mbicv", psi0 = 0.9, presel = TRUE,
                                             trunc_lvl = NA, tree_crit = "tau", threshold = 0, keep_data = FALSE, show_trace =
                                               FALSE, cores = 6),
                      cores=6)
    
    
    pdens = dvine(tenvdat, vine=vine_pres, cores=6)
    pdens[pdens==0] = min(pdens[pdens>0],na.rm=T)
    bdens = dvine(tenvdat, vine=vine_back, cores=6)
    bdens[bdens==0] = min(bdens[bdens>0],na.rm=T)
    
    prob = pdens/bdens
    nprob = prob / sum(prob,na.rm=T)
    
    vine_preds <- data.frame(nprob=nprob, pa=tpadat)
    
    vine_AUC = dismo::evaluate(p=as.vector(vine_preds[vine_preds$pa==1,1]), a=as.vector(vine_preds[vine_preds$pa==0,1]), tr=0.2)@auc
    
    #B) Maxnet - Four variables with highest permutation importance
    maxnet_model = maxnet(padat, envdat, addsamplestobackground=T)
    maxnet_preds = predict(maxnet_model, tenvdat, clamp=T, type=c("exponential")) #"link","exponential","cloglog","logistic"
    maxnet_preds = maxnet_preds / sum(maxnet_preds,na.rm=T)
    maxnet_preds = data.frame(nprob=maxnet_preds, pa=tpadat)
    maxnet_AUC = dismo::evaluate(p=as.vector(maxnet_preds[maxnet_preds$pa==1,1]), a=as.vector(maxnet_preds[maxnet_preds$pa==0,1]), tr=0.2)@auc
    
    #C) Four variables with highest permutation importance
    presdat = PO[PO$spid==species,varnames]
    fourvars = mostimportantvars(presdataset=presdat, backdataset=back[,varnames], maxnetmod=maxnet_model)
    
    mins4 <- sapply(fourvars, function(x){min(envdat[,x])})
    maxs4 <- sapply(fourvars, function(x){max(envdat[,x])})
    
    bdata4 <- as.matrix(back[,fourvars])
    pdata4 <- as.matrix(PO[PO$spid==species,fourvars])
    envdat4 = envdat[,fourvars]
    tenvdat4 = tenvdat[,fourvars]
    
    #D) Maxnet with four variables
    maxnet4_model = maxnet(padat, envdat4, addsamplestobackground=T)
    maxnet4_preds = predict(maxnet4_model, tenvdat4, clamp=T, type=c("exponential")) #"link","exponential","cloglog","logistic"
    maxnet4_preds = maxnet4_preds/sum(maxnet4_preds, na.rm=T)
    maxnet4_preds = data.frame(nprob=maxnet4_preds, pa=tpadat)
    
    maxnet4_AUC = dismo::evaluate(p=as.vector(maxnet4_preds[maxnet4_preds$pa==1,1]), a=as.vector(maxnet4_preds[maxnet4_preds$pa==0,1]), tr=0.2)@auc
    
    #E) Vine with four variables
    vine4_back <- vine(bdata4, 
                       margins_controls = list(mult = 1, xmin = mins4, xmax = maxs4, bw = NA, deg = 2), #mult=NULL means  log(1 + d) after rvinecopulib:::expand_factors(). This is only true for discrete variables. 
                       copula_controls = list(family_set = "nonparametric", structure = NA, par_method = "mle",
                                              nonpar_method = "quadratic", mult = 1, selcrit = "mbicv", psi0 = 0.9, presel = TRUE,
                                              trunc_lvl = NA, tree_crit = "tau", threshold = 0, keep_data = FALSE, show_trace =
                                                FALSE, cores = 6),
                       cores=6)
    
    vine4_pres <- vine(pdata4, 
                       margins_controls = list(mult = 1, xmin = mins4, xmax = maxs4, bw = NA, deg = 2), #mult=NULL means  log(1 + d) after rvinecopulib:::expand_factors(). This is only true for discrete variables. 
                       copula_controls = list(family_set = "nonparametric", structure = NA, par_method = "mle",
                                              nonpar_method = "quadratic", mult = 1, selcrit = "mbicv", psi0 = 0.9, presel = TRUE,
                                              trunc_lvl = NA, tree_crit = "tau", threshold = 0, keep_data = FALSE, show_trace =
                                                FALSE, cores = 6),
                       cores=6)
    
    pdens = dvine(tenvdat4, vine=vine4_pres, cores=6)
    pdens[pdens==0] = min(pdens[pdens>0],na.rm=T)
    bdens = dvine(tenvdat4, vine=vine4_back, cores=6)
    bdens[bdens==0] = min(bdens[bdens>0],na.rm=T)
    
    prob = pdens/bdens
    nprob = prob / sum(prob,na.rm=T)
    
    vine4_preds <- data.frame(nprob=nprob, pa=tpadat)
    
    vine4_AUC = dismo::evaluate(p=as.vector(vine4_preds[vine4_preds$pa==1,1]), a=as.vector(vine4_preds[vine4_preds$pa==0,1]), tr=0.2)@auc
    
    
    #F) Fit Multivariate KDE
    mins <- sapply(names(envdat4), function(x){min(envdat4[,x])})
    maxs <- sapply(names(envdat4), function(x){max(envdat4[,x])})
    
    bdata <- as.matrix(envdat4)
    pdata <- as.matrix(PO[PO$spid==species,fourvars])
    
    BHpi3 <- Hpi(x = bdata)*3
    kde_back <- kde(bdata,xmin=mins,xmax=maxs, binned=T, H=BHpi3, bgridsize=rep(10,4))
    kde_back$estimate[kde_back$estimate<0 | kde_back$estimate==0] = min(kde_back$estimate[kde_back$estimate>0])
    
    PHpi3 <- Hpi(x = pdata)*3
    kde_pres <- kde(pdata,xmin=mins,xmax=maxs, binned=T, H=PHpi3, bgridsize=rep(10,4))
    kde_pres$estimate[kde_pres$estimate<0 | kde_pres$estimate==0] = min(kde_pres$estimate[kde_pres$estimate>0])
    
    dkde_back <- dkde(fhat=kde_back,as.matrix(tenvdat4))
    dkde_pres <- dkde(fhat=kde_pres,as.matrix(tenvdat4)) #sum(dkde_pres)
    
    kde_preds <- (dkde_pres)/(dkde_back)
    kde_npreds = kde_preds/sum(kde_preds, na.rm=T)
    kde_preds <- data.frame(nprob=kde_npreds, pa=tpadat)
    
    kde_AUC = dismo::evaluate(p=as.vector(kde_preds[kde_preds$pa==1,1]), a=as.vector(kde_preds[kde_preds$pa==0,1]), tr=0.2)@auc
    
    
    #results
    preds = data.frame("maxnet" = maxnet_preds$nprob, "maxnet4" = maxnet4_preds$nprob, "kde4" = kde_preds$nprob, "vine" = vine_preds$nprob, "vine4" = vine4_preds$nprob, "PA" = tpadat)
    results = c("maxnet" = as.numeric(maxnet_AUC), "maxnet4"= as.numeric(maxnet4_AUC), "kde4" = as.numeric(kde_AUC), "vine"= as.numeric(vine_AUC), "vine4"= as.numeric(vine4_AUC))
    
    #preds = data.frame("maxnet" = maxnet_preds$nprob, "maxnet4" = maxnet4_preds$nprob, "kde4" = kde_preds$nprob, "vine" = vine_preds$nprob, "vine4" = vine4_preds$nprob, "PA" = tpadat)
    results = c("maxnet" = as.numeric(maxnet_AUC), 
                "maxnet4"= as.numeric(maxnet4_AUC), 
                "kde4" = as.numeric(kde_AUC), 
                "vine"= as.numeric(vine_AUC), 
                "vine4"= as.numeric(vine4_AUC)
    )
    
    return(results)
    
  },error=function(e){
    # if(!exists("abc")) { abc <- NA }
    # if(!exists("def")) { def <- 5 }
    # results = c("maxnet" = maxnet_AUC, "maxnet4" = maxnet4_AUC, "kde" = kde_AUC, "species" = species)
    return(c("maxnet" = NA, 
             "maxnet4"= NA, 
             "kde4" = NA, 
             "vine"= NA, 
             "vine4"= NA))
  }
  )}






###
### Case study 1
###
findoptimalMaxent <- function(data, features, regularization){
  tryCatch({
    data <- data
    
    adjustgrid <- expand.grid(features=features, reg=regularization)
    #prior <- length(data@pa[data@pa==1]) / length(data@pa)
    
    
    outputdf <- data.frame(features=rep(0,length(data)), regularization=rep(0,length(data)), logloss=rep(0,length(data)),
                           logsum=rep(0,length(data)))
    
    
    for (i in 1:length(data)){
      svMisc::progress(i)
      #datapartition
      internaldf <- cbind(adjustgrid, data.frame(logloss=rep(0,nrow(adjustgrid)), logsum=rep(0,nrow(adjustgrid))))
      
      #predictdata  <- as.data.frame(valNB[,1])
      #names(predictdata)[1] <- variable
      
      #fit different adjusts
      for (j in 1:nrow(adjustgrid)){
        
        model.maxnet <- train(method = "Maxnet", data = data[[i]][[1]], fc = as.character(adjustgrid[j,1]), reg=adjustgrid[j,2])
        
        
        output <- data.frame(prob=predict(model.maxnet, data = data[[i]][[2]], type = "exponential"), pa=data[[i]][[2]]@pa)
        output$pnorm <- output$prob/sum(output$prob)
        logloss <- -sum(log(output[output$pa==1,3]))
        logsum <- sum(output[output$pa==1,3])
        
        internaldf[j,3:4] <- c(logloss, logsum)
      }
      
      best <- internaldf[which(internaldf[,3]==min(internaldf[,3]))[1],]#if there are multiple, it will take the first one. 
      outputdf[i,1] <- as.character(best[,1])
      outputdf[i,2:4] <- best[,2:4]
      
    }
    return(outputdf)
    
  },error=function(e){})
}


fitoptimalMaxentEach <- function(data, settings){
  tryCatch({
    data <- data
    settings <- settings
    
    outputlist <- vector(mode = "list", length = length(data))
    
    for (i in 1:length(data)){
      svMisc::progress(i)
      
      dataset <- data[[i]]
      
      train <- dataset[[1]]
      test <- dataset[[3]]
      
      model.maxnet <- train(method = "Maxnet", data = train, fc = as.character(settings[i,1]), reg=as.numeric(settings[i,2]))
      outputlist[[i]] <- data.frame(prob=predict(model.maxnet, data = test, type = "exponential"), pa=test@pa)
    }
    return(outputlist)
  },error=function(e){})
}

fitMaxentNoTune <- function(data){
  tryCatch({
    data <- data
    settings <- settings
    
    outputlist <- vector(mode = "list", length = length(data))
    
    for (i in 1:length(data)){
      svMisc::progress(i)
      
      dataset <- data[[i]]
      
      train <- dataset[[1]]
      test <- dataset[[3]]
      
      model.maxnet <- train(method = "Maxnet", data = train)
      outputlist[[i]] <- data.frame(prob=predict(model.maxnet, data = test, type = "exponential"), pa=test@pa)
    }
    return(outputlist)
  },error=function(e){})
}

myCI = function(x, ci=0.95){
  a <- mean(x, na.rm=T)
  s <- sd(x, na.rm=T)
  n <- length(x[!is.na(x)])
  error <- qt(ci + (1 - ci)/2, df = n - 1) * s/sqrt(n)
  return(c(upper = a + error, mean = a, lower = a - error))
}

getmetrics2 = function(data) {
  data$nprob = data$prob/sum(data$prob)
  
  pbest <- data$pa / sum(data$pa)
  pbase <- rep(1/length(data$pa), length(data$pa))
  
  Pmean <- mean(data$nprob[data$pa==1])
  Psd <- sd(data$nprob[data$pa==1])
  Bmean <- mean(data$nprob[data$pa==0])
  Bsd <- sd(data$nprob[data$pa==0])
  nloglike <- -sum(log(data$nprob[data$pa==1]))
  psum <- sum(data$nprob[data$pa==1]) 
  auc <- dismo::evaluate(p=as.vector(data[data$pa==1,3]), a=as.vector(data[data$pa==0,3]), tr=0.2)@auc
  prauc <- PRROC::pr.curve(data$nprob[data$pa==1], data$nprob[data$pa==0])$auc.integral
  me <- Metrics::bias(actual=data$nprob, predicted= pbest)
  mae <- Metrics::mae(actual=data$nprob, predicted= pbest)
  mse <- Metrics::mse(actual=data$nprob, predicted= pbest) #brier score
  mses <- 1-mse/Metrics::mse(actual=pbase, predicted= pbest)
  pme <- Metrics::bias(actual=data$nprob[data$pa==1], predicted= pbest[data$pa==1]) #mean(actual-predicted) # 0 means perfect model. 
  pmae <- Metrics::mae(actual=data$nprob[data$pa==1], predicted= pbest[data$pa==1]) #absolute, so negative values ignored / absolute. 
  pmse <- Metrics::mse(actual=data$nprob[data$pa==1], predicted= pbest[data$pa==1]) #brier score
  pmses <- 1-pmse/Metrics::mse(actual=pbase[data$pa==1], predicted= pbest[data$pa==1])
  
  output <- data.frame(
    Pmean = Pmean,
    Psd = Psd,
    Bmean = Bmean,
    Bsd = Bsd,
    nloglike = nloglike,
    psum = psum, 
    auc = auc,
    prauc = prauc,
    me = me,
    mae = mae,
    mse = mse,
    mses = mses,
    pme = pme,
    pmae = pmae,
    pmse = pmse,
    pmses = pmses
  )
}

calclistmetrics2 <- function(predlist, seeds){
  #outputdf <- data.frame(seed=1:length(predlist), logloss=rep(0,length(predlist)), logsum=rep(0,length(predlist)))
  outputdf <- matrix(0, nrow=length(predlist), ncol=17)
  outputdf[,1] <- seeds
  for (i in 1:length(predlist)){
    
    listdf <- predlist[[i]]
    listdf$pa <- as.numeric(as.character(listdf$pa))
    listdf$nprob <- listdf$prob/sum(listdf$prob)
    
    metrics = getmetrics2(listdf)
    outputdf[i,2:17] <- as.numeric(metrics)
    #outputdf[i,2] <- -sum(log(nprob[which(listdf$pa==1)]))
    #outputdf[i,3] <- sum(nprob[which(listdf$pa==1)])
  }
  colnames(outputdf) <- c("seed",colnames(metrics))
  outputdf <- as.data.frame(outputdf)
  return(outputdf)
}


findoptimalmethod3b_optim2 <- function(data, par, maxit=500){
  tryCatch({
    data <- data
    
    prior <- 0.001
    
    
    #outputdf <- data.frame(seed=seeds,adjusts0=rep(0,length(seeds)), adjusts1=rep(0,length(seeds)), logloss=rep(0,length(seeds)),
    #                       logsum=rep(0,length(seeds)))
    outputdf <-as.data.frame(matrix(0, ncol = (3+2*ncol(data[[1]][[1]]@data)), nrow = length(data)))
    outputdf[,1] <- seq_len(length(data))
    for (i in 1:length(data)){
      svMisc::progress(i)
      #datapartition
      
      datasets <- data[[i]]
      
      train <- datasets[[1]]
      val <- datasets[[2]]
      
      trainNB <- cbind(train@data, pa=as.factor(train@pa))
      valNB <- cbind(val@data, pa=as.factor(val@pa))
      
      #predictdata  <- as.data.frame(valNB[,1])
      #names(predictdata)[1] <- variable
      
      #Fit model with adjust=1
      grid.NB = data.frame(usekernel=TRUE,laplace=T,adjust=1)
      
      model.NB = caret::train(pa ~ .,data=trainNB,method="naive_bayes",
                              trControl=trainControl(method="none"),
                              tuneGrid=grid.NB) #class(model.NB$finalModel)
      
      #fit different adjusts
      NLL.function <- function(par){
        
        parexp <- exp(par)
        for (k in 1:ncol(train@data)){
          
          density0 <-  density(trainNB[trainNB$pa==0,k], adjust=parexp[2*k])
          density1 <-  density(trainNB[trainNB$pa==1,k], adjust=parexp[2*k-1])
          
          model.NB <- updateMBvar(model=model.NB, variable=colnames(train@data)[k], densityobject0=density0, densityobject1=density1)
        }
        
        output <- data.frame(prob=predict.nonparametric_naive_bayes2(object=model.NB$finalModel, newdata=as.matrix(valNB[,-ncol(valNB)]), prior=c(1, prior), type="prob")[,2], pa=valNB$pa)
        
        output$pnorm <- output$prob/sum(output$prob)          
        logloss <- -sum(log(output[output$pa==1,3]))
        
      }
      
      par <- par
      optim.result <- optim(par,NLL.function,method="Nelder-Mead", control=list(maxit=maxit, parscale=exp(par)))  #conv=0, value=946.1976
      outputdf[i,2:(2*ncol(train@data)+2)] <- c(exp(optim.result$par), optim.result$value)
      outputdf[i,(2*ncol(train@data)+3)] <- optim.result$convergence
    }
    return(outputdf)
    
  },error=function(e){})
}

predict.nonparametric_naive_bayes2 <- function (object, newdata = NULL, prior, type = c("class", "prob"), threshold = 0.001, eps = 0, ...) {
  
  if (is.null(newdata))
    newdata <- object$data$x
  if (!is.matrix(newdata))
    stop("predict.nonparametric_naive_bayes(): newdata must be a numeric matrix with at least one row and two columns.", call. = FALSE)
  if (mode(newdata) != "numeric")
    stop("predict.nonparametric_naive_bayes(): newdata must be a numeric matrix.", call. = FALSE)
  
  if (threshold < 0)
    stop("predict.nonparametric_naive_bayes(): threshold must be non-negative.", call. = FALSE)
  if (eps < 0)
    stop("predict.nonparametric_naive_bayes(): eps must be non-negative.", call. = FALSE)
  
  type <- match.arg(type)
  lev <- object$levels
  n_lev <- length(lev)
  n_obs <- dim(newdata)[1L]
  prior <- prior
  col_names <- colnames(newdata)
  tables <- object$table
  features <- colnames(newdata)[colnames(newdata) %in% names(tables)]
  n_tables <- length(tables)
  n_features <- length(features)
  n_features_newdata <- ncol(newdata)
  
  if (n_features == 0) {
    warning(paste0("predict.nonparametric_naive_bayes(): no feature in newdata corresponds to ",
                   "features defined in the object. Classification is based on prior probabilities"), call. = FALSE)
    if (type == "class") {
      return(factor(rep(lev[which.max(prior)], n_obs), levels = lev))
    } else {
      return(matrix(prior, ncol = n_lev, nrow = n_obs, byrow = TRUE, dimnames = list(NULL, lev)))
    }
  }
  if (n_features < n_tables) {
    warning(paste0("predict.nonparametric_naive_bayes(): only ", n_features, " feature(s) in newdata could be matched ",
                   "with ", n_tables, " feature(s) defined in the object."), call. = FALSE)
  }
  if (n_features_newdata > n_features) {
    warning(paste0("predict.nonparametric_naive_bayes(): newdata contains feature(s) that could not be matched ",
                   "with (", n_features, ") feature(s) defined in the object. Only matching features are used for calculation."), call. = FALSE)
    newdata <- newdata[ ,features, drop = FALSE]
  }
  NAx <- anyNA(newdata)
  if (NAx) {
    len_na <- sum(is.na(newdata))
    if (len_na > 0)
      warning(paste0("predict.nonparametric_naive_bayes(): ", len_na, " missing", ifelse(len_na == 1, " value", " values"),
                     " discovered in the newdata. ", ifelse(len_na == 1, "It is", "They are"), " not included into the calculation."), call. = FALSE)
    na <- apply(newdata, 2, anyNA)
  }
  eps <- ifelse(eps == 0, log(.Machine$double.xmin), log(eps))
  threshold <- log(threshold)
  
  post <- matrix(log(prior), nrow = n_obs, ncol = n_lev, byrow = TRUE)
  colnames(post) <- lev
  
  for (var in features) {
    V <- newdata[ ,var]
    tab <- tables[[var]]
    logp <- sapply(lev, function(z) {
      dens <- tab[[z]]
      log(stats::approx(dens$x, dens$y, xout = V, rule = 2, ties = "ordered")$y)
    })
    if (NAx) { if (na[var]) { logp[is.na(logp)] <- 0 }}
    logp[logp <= eps] <- threshold
    post <- post + logp
  }
  
  if (type == "class") {
    if (n_obs == 1) {
      return(factor(lev[which.max(post)], levels = lev))
    } else {
      return(factor(lev[max.col(post, "first")], levels = lev))
    }
  } else {
    if (n_obs == 1) {
      post <- t(as.matrix(apply(post, 2, function(x) { exp(x-post[,1]) })))
      colnames(post) <- lev
      return(post)
    } else {
      return(apply(post, 2, function(x) { exp(x-post[,1]) }))
    }
  }
}

#updateMBvar
updateMBvar <- function(model, variable, densityobject0=NULL, densityobject1=NULL){
  if(is.null(model)) stop("non-existing model")
  if (is.null(densityobject0) && is.null(densityobject1)) stop("Error: both densities can't be null")
  
  model <- model
  
  if(is.null(eval(parse(text=paste("model$finalModel$tables$",variable, sep=""))))) stop("variable not in model")
  
  variable <- variable 
  
  if(!is.null(densityobject0)){
    densityobject0 <- densityobject0
    eval(parse(text=paste("model$finalModel$tables$",variable,"$'0'$x",  " <- densityobject0$x", sep="")))
    eval(parse(text=paste("model$finalModel$tables$",variable,"$'0'$y",  " <- densityobject0$y", sep="")))
  }
  if(!is.null(densityobject1)){
    densityobject1 <- densityobject1
    eval(parse(text=paste("model$finalModel$tables$",variable,"$'1'$x",  " <- densityobject1$x", sep="")))
    eval(parse(text=paste("model$finalModel$tables$",variable,"$'1'$y",  " <- densityobject1$y", sep="")))  
  }
  return(model)
}


fitoptimalmodel.method1 <- function(data, settings){
  tryCatch({
    data <- data
    settings <- settings
    seeds <- seq_len(length(data))
    varnames <- colnames(data[[1]][[1]]@data)
    prior <- 0.01
    
    outputlist <- vector(mode = "list", length = length(data))
    
    for (i in 1:length(seeds)){
      svMisc::progress(i)
      #datapartition
      datasets <- data[[i]] 
      
      train <- datasets[[1]]
      test <- datasets[[3]]
      trainNB <- cbind(train@data, pa=as.factor(train@pa))
      testNB <- cbind(test@data, pa=as.factor(test@pa))
      #eval(parse(text=paste("trainNB <- data.frame(",variable, "=train@data$\"", variable, "\",pa=as.factor(train@pa))", sep=""))) 
      #eval(parse(text=paste("testNB <- data.frame(",variable, "=test@data$\"", variable, "\",pa=as.factor(test@pa))", sep=""))) 
      #predictdata  <- as.data.frame(testNB[,1])
      #names(predictdata)[1] <- variable
      
      #Fit model with adjust=1
      grid.NB = data.frame(usekernel=TRUE,laplace=T,adjust=1)
      
      model.NB = caret::train(pa ~ .,data=trainNB,method="naive_bayes",
                              trControl=trainControl(method="none"),
                              tuneGrid=grid.NB) #class(model.NB$finalModel)
      
      #fit different adjusts
      for (j in 1:length(settings)){
        
        density0 <-  density(trainNB[trainNB$pa==0,j], adjust=settings[[j]]$adjusts0[i])
        density1 <-  density(trainNB[trainNB$pa==1,j], adjust=settings[[j]]$adjusts1[i])
        
        model.NB <- updateMBvar(model=model.NB, variable=varnames[j], densityobject0=density0, densityobject1=density1)
        
      }
      
      outputlist[[i]] <- data.frame(prob=predict.nonparametric_naive_bayes2(object=model.NB$finalModel, newdata=as.matrix(testNB[,1:length(settings)]), prior=c(1, prior), type="prob")[,2], pa=testNB$pa)
      
    }
    return(outputlist)
    
  },error=function(e){})
}

######Gaussian copula

##old functions
findoptimalmethod3bgauscop_optim <- function(data, par, seeds, maxit=500){
  
  data <- data
  maxit <- maxit
  prior <- 1
  ncols = ncol(data[[1]][[1]]@data)
  
  #outputdf <- data.frame(seed=seeds,adjusts0=rep(0,length(seeds)), adjusts1=rep(0,length(seeds)), logloss=rep(0,length(seeds)),
  #                       logsum=rep(0,length(seeds)))
  outputdf <-as.data.frame(matrix(0, ncol = (4+2*ncols), nrow = length(seeds)))
  outputdf[,1] <- seeds
  
  for (i in 1:length(seeds)){
    tryCatch({
      svMisc::progress(i)
      
      #datapartition
      
      train <- data[[i]][[1]]
      val <- data[[i]][[2]]
      
      trainNB <- cbind(train@data, pa=as.factor(train@pa))
      valNB <- cbind(val@data, pa=as.factor(val@pa))
      
      #predictdata  <- as.data.frame(valNB[,1])
      #names(predictdata)[1] <- variable
      
      #Fit model with adjust=1
      grid.NB = data.frame(usekernel=TRUE,laplace=T,adjust=1)
      
      model.NB = caret::train(pa ~ .,data=trainNB,method="naive_bayes",
                              trControl=trainControl(method="none"),
                              tuneGrid=grid.NB) #class(model.NB$finalModel)
      
      copula1_C = cor(train@data[which(train@pa==1),], method = "spearman")
      copula1_rhos = corexport(copula1_C)
      copula1_est_rho = as.numeric(2*sin(copula1_rhos*pi/6))
      
      copula0_C = cor(train@data[which(train@pa==0),], method = "spearman")
      copula0_rhos = corexport(copula0_C)
      copula0_est_rho = as.numeric(2*sin(copula0_rhos*pi/6))
      
      
      #fit different adjusts
      NLL.function <- function(par){
        
        parexp <- exp(par)
        for (k in 1:ncols){
          
          density1 <-  density(trainNB[trainNB$pa==1,k], adjust=parexp[2*k-1])
          density0 <-  density(trainNB[trainNB$pa==0,k], adjust=parexp[2*k])
          
          
          model.NB <- updateMBvar(model=model.NB, variable=colnames(data[[1]][[1]]@data)[k], densityobject0=density0, densityobject1=density1)
        }
        
        mprob <- predict.nonparametric_naive_bayes2(object=model.NB$finalModel, newdata=as.matrix(valNB[,-ncol(valNB)]), prior=c(1, prior), type="prob")[,2]
        #mprob <- mprob/sum(mprob)
        #sum(is.na(mprob))
        #sum(is.infinite(mprob))
        #sum(mprob==0)
        
        pu <- getunifPKDE2(data=val@data, naivebayesobject=model.NB, class=1)
        bu <- getunifPKDE2(data=val@data, naivebayesobject=model.NB, class=0)
        
        # copula1_density_p <- dCopula(as.matrix(pu[val@pa==1,]), copula=normalCopula(param = copula1_est_rho, dim = ncols, dispstr = "un"))
        # copula1_density_b <- dCopula(as.matrix(bu[val@pa==0,]), copula=normalCopula(param = copula1_est_rho, dim = ncols, dispstr = "un"))
        # copula1_density <- c(copula1_density_p, copula1_density_b)
        copula1_density <- dCopula(as.matrix(pu), copula=normalCopula(param = copula1_est_rho, dim = ncols, dispstr = "un"))
        
        # copula0_density_p <- dCopula(as.matrix(pu[val@pa==1,]), copula=normalCopula(param = copula0_est_rho, dim = ncols, dispstr = "un"))
        # copula0_density_b <- dCopula(as.matrix(bu[val@pa==0,]), copula=normalCopula(param = copula0_est_rho, dim = ncols, dispstr = "un"))
        # copula0_density <- c(copula0_density_p, copula0_density_b)
        copula0_density <- dCopula(as.matrix(bu), copula=normalCopula(param = copula0_est_rho, dim = ncols, dispstr = "un"))
        
        cprob <- copula1_density / copula0_density
        #cprob <- cprob/sum(cprob) 
        #sum(is.na(cprob))
        #sum(is.infinite(cprob))
        #sum(cprob==0)
        
        output <- data.frame(prob=mprob*cprob, pa=valNB$pa)
        output$pnorm <- output$prob/sum(output$prob)
        
        logloss <- -sum(log(output[output$pa==1,3]))
        return(logloss)
        #sum(is.na(output$pnorm)); sum(is.infinite(output$pnorm))
      }
      #NLL.function(rep(log(6),14))
      par <- par
      #      optim.result <- optim(par,NLL.function,method="Nelder-Mead", control=list(maxit=maxit, parscale=exp(par)))  #conv=0, value=946.1976
      #      outputdf[i,2:(2*ncol(data@data)+2)] <- c(exp(optim.result$par), optim.result$value)
      #      outputdf[i,(2*ncol(data@data)+3)] <- optim.result$convergence
      
      #optim.result <- optim(par, NLL.function,method="L-BFGS-B", lower=rep(log(0.1),14), upper=rep(log(50), 14), control=list(maxit=maxit, parscale=exp(par)))    #, parscale=exp(par)
      optim.result <- optim(par, NLL.function,method="Nelder-Mead", control=list(maxit=maxit, parscale=exp(par)))   
      
      #optim.result <- optimParallel(par=par, fn=NLL.function,
      #                               method = "L-BFGS-B", lower=rep(log(0.1),14), upper=rep(log(50), 14),control=list(maxit=maxit, parscale=exp(par)))
      
      outputdf[i,2:(2*ncols+3)] <- c(exp(optim.result$par), optim.result$value, optim.result$counts[1])
      outputdf[i,(2*ncols+4)] <- optim.result$convergence
    },error=function(e){})
  }
  
  return(outputdf)
}

findoptimalmethod3bgauscop_optimparallel <- function(data, par, seeds, maxit=500){
  
  data <- data
  maxit <- maxit
  prior <- 1
  ncols = ncol(data[[1]][[1]]@data)
  
  #outputdf <- data.frame(seed=seeds,adjusts0=rep(0,length(seeds)), adjusts1=rep(0,length(seeds)), logloss=rep(0,length(seeds)),
  #                       logsum=rep(0,length(seeds)))
  outputdf <-as.data.frame(matrix(0, ncol = (4+2*ncols), nrow = length(seeds)))
  outputdf[,1] <- seeds
  
  for (i in 1:length(seeds)){
    tryCatch({
      svMisc::progress(i)
      
      #datapartition
      
      train <- data[[i]][[1]]
      val <- data[[i]][[2]]
      
      trainNB <- cbind(train@data, pa=as.factor(train@pa))
      valNB <- cbind(val@data, pa=as.factor(val@pa))
      
      #predictdata  <- as.data.frame(valNB[,1])
      #names(predictdata)[1] <- variable
      
      #Fit model with adjust=1
      grid.NB = data.frame(usekernel=TRUE,laplace=T,adjust=1)
      
      model.NB = caret::train(pa ~ .,data=trainNB,method="naive_bayes",
                              trControl=trainControl(method="none"),
                              tuneGrid=grid.NB) #class(model.NB$finalModel)
      
      copula1_C = cor(train@data[which(train@pa==1),], method = "spearman")
      copula1_rhos = corexport(copula1_C)
      copula1_est_rho = as.numeric(2*sin(copula1_rhos*pi/6))
      
      copula0_C = cor(train@data[which(train@pa==0),], method = "spearman")
      copula0_rhos = corexport(copula0_C)
      copula0_est_rho = as.numeric(2*sin(copula0_rhos*pi/6))
      
      
      #fit different adjusts
      NLL.function <- function(par){
        
        parexp <- exp(par)
        for (k in 1:ncols){
          
          density1 <-  density(trainNB[trainNB$pa==1,k], adjust=parexp[2*k-1])
          density0 <-  density(trainNB[trainNB$pa==0,k], adjust=parexp[2*k])
          
          
          model.NB <- updateMBvar(model=model.NB, variable=colnames(data[[1]][[1]]@data)[k], densityobject0=density0, densityobject1=density1)
        }
        
        mprob <- predict.nonparametric_naive_bayes2(object=model.NB$finalModel, newdata=as.matrix(valNB[,-ncol(valNB)]), prior=c(1, prior), type="prob")[,2]
        #mprob <- mprob/sum(mprob)
        #sum(is.na(mprob))
        #sum(is.infinite(mprob))
        #sum(mprob==0)
        
        pu <- getunifPKDE2(data=val@data, naivebayesobject=model.NB, class=1)
        bu <- getunifPKDE2(data=val@data, naivebayesobject=model.NB, class=0)
        
        # copula1_density_p <- dCopula(as.matrix(pu[val@pa==1,]), copula=normalCopula(param = copula1_est_rho, dim = ncols, dispstr = "un"))
        # copula1_density_b <- dCopula(as.matrix(bu[val@pa==0,]), copula=normalCopula(param = copula1_est_rho, dim = ncols, dispstr = "un"))
        # copula1_density <- c(copula1_density_p, copula1_density_b)
        copula1_density <- dCopula(as.matrix(pu), copula=normalCopula(param = copula1_est_rho, dim = ncols, dispstr = "un"))
        
        # copula0_density_p <- dCopula(as.matrix(pu[val@pa==1,]), copula=normalCopula(param = copula0_est_rho, dim = ncols, dispstr = "un"))
        # copula0_density_b <- dCopula(as.matrix(bu[val@pa==0,]), copula=normalCopula(param = copula0_est_rho, dim = ncols, dispstr = "un"))
        # copula0_density <- c(copula0_density_p, copula0_density_b)
        copula0_density <- dCopula(as.matrix(bu), copula=normalCopula(param = copula0_est_rho, dim = ncols, dispstr = "un"))
        
        cprob <- copula1_density / copula0_density
        #cprob <- cprob/sum(cprob) 
        #sum(is.na(cprob))
        #sum(is.infinite(cprob))
        #sum(cprob==0)
        
        output <- data.frame(prob=mprob*cprob, pa=valNB$pa)
        output$pnorm <- output$prob/sum(output$prob)
        
        logloss <- -sum(log(output[output$pa==1,3]))
        return(logloss)
        #sum(is.na(output$pnorm)); sum(is.infinite(output$pnorm))
      }
      #NLL.function(rep(log(6),14))
      par <- par
      #      optim.result <- optim(par,NLL.function,method="Nelder-Mead", control=list(maxit=maxit, parscale=exp(par)))  #conv=0, value=946.1976
      #      outputdf[i,2:(2*ncol(data@data)+2)] <- c(exp(optim.result$par), optim.result$value)
      #      outputdf[i,(2*ncol(data@data)+3)] <- optim.result$convergence
      
      #optim.result <- optim(par, NLL.function,method="L-BFGS-B", lower=rep(log(0.1),14), upper=rep(log(50), 14), control=list(maxit=maxit, parscale=exp(par)))    #, parscale=exp(par)
      #optim.result <- optim(par, NLL.function,method="Nelder-Mead", control=list(maxit=maxit, parscale=exp(par)))   
      
      optim.result <- optimParallel(par=par, fn=NLL.function,
                                    method = "L-BFGS-B", lower=rep(log(0.1),14), upper=rep(log(50), 14),control=list(maxit=maxit, parscale=exp(par)))
      
      outputdf[i,2:(2*ncols+3)] <- c(exp(optim.result$par), optim.result$value, optim.result$counts[1])
      outputdf[i,(2*ncols+4)] <- optim.result$convergence
    },error=function(e){})
  }
  
  return(outputdf)
}


fitoptimalmodel.mbgauscop <- function(data, settings){
  
  data <- data
  settings <- settings
  seeds <- settings[[1]]$seed
  varnames <- colnames(data[[1]][[1]]@data)
  prior <- 1
  ncols = ncol(data[[1]][[1]]@data)
  
  outputlist <- vector(mode = "list", length = length(seeds))
  
  for (i in 1:length(seeds)){
    tryCatch({
      svMisc::progress(i)
      #datapartition
      train <- data[[i]][[1]]
      test <- data[[i]][[3]]
      trainNB <- cbind(train@data, pa=as.factor(train@pa))
      testNB <- cbind(test@data, pa=as.factor(test@pa))
      
      #copula cors MoM here 
      copula1_C = cor(train@data[which(train@pa==1),], method = "spearman")
      copula1_rhos = corexport(copula1_C)
      copula1_est_rho = as.numeric(2*sin(copula1_rhos*pi/6))
      
      copula0_C = cor(train@data[which(train@pa==0),], method = "spearman")
      copula0_rhos = corexport(copula0_C)
      copula0_est_rho = as.numeric(2*sin(copula0_rhos*pi/6))
      #/copula cors MoM here
      
      #Fit model with adjust=1
      grid.NB = data.frame(usekernel=TRUE,laplace=T,adjust=1)
      
      model.NB = caret::train(pa ~ .,data=trainNB,method="naive_bayes",
                              trControl=trainControl(method="none"),
                              tuneGrid=grid.NB) #class(model.NB$finalModel)
      
      #fit different adjusts
      for (j in 1:length(settings)){
        
        density0 <-  density(trainNB[trainNB$pa==0,j], adjust=settings[[j]]$adjusts0[i])
        density1 <-  density(trainNB[trainNB$pa==1,j], adjust=settings[[j]]$adjusts1[i])
        
        model.NB <- updateMBvar(model=model.NB, variable=varnames[j], densityobject0=density0, densityobject1=density1)
        
      }
      
      mprob <- predict.nonparametric_naive_bayes2(object=model.NB$finalModel, newdata=as.matrix(testNB[,-ncol(testNB)]), prior=c(1, prior), type="prob")[,2]
      
      pu <- getunifPKDE2(data=test@data, naivebayesobject=model.NB, class=1)
      bu <- getunifPKDE2(data=test@data, naivebayesobject=model.NB, class=0)
      
      
      # copula1_density_p <- dCopula(as.matrix(pu[test@pa==1,]), copula=normalCopula(param = copula1_est_rho, dim = ncols, dispstr = "un"))
      # copula1_density_b <- dCopula(as.matrix(bu[test@pa==0,]), copula=normalCopula(param = copula1_est_rho, dim = ncols, dispstr = "un"))
      # copula1_density <- c(copula1_density_p, copula1_density_b)
      copula1_density <- dCopula(as.matrix(pu), copula=normalCopula(param = copula1_est_rho, dim = ncols, dispstr = "un"))
      
      # copula0_density_p <- dCopula(as.matrix(pu[test@pa==1,]), copula=normalCopula(param = copula0_est_rho, dim = ncols, dispstr = "un"))
      # copula0_density_b <- dCopula(as.matrix(bu[test@pa==0,]), copula=normalCopula(param = copula0_est_rho, dim = ncols, dispstr = "un"))
      # copula0_density <- c(copula0_density_p, copula0_density_b)
      copula0_density <- dCopula(as.matrix(bu), copula=normalCopula(param = copula0_est_rho, dim = ncols, dispstr = "un"))
      
      cprob <- copula1_density / copula0_density
      
      
      #outputlist[[i]] <- data.frame(prob=predict.nonparametric_naive_bayes2(object=model.NB$finalModel, newdata=as.matrix(testNB[,1:length(settings)]), prior=c(1, prior), type="prob")[,2], pa=testNB$pa)
      outputlist[[i]] <- data.frame(prob=mprob*cprob, pa=testNB$pa)
    },error=function(e){})
  }
  return(outputlist)
  
}

fitoptimalmodel.mbgauscop_correct <- function(data, settings){
  
  data <- data
  settings <- settings
  seeds <- settings[[1]]$seed
  varnames <- colnames(data[[1]][[1]]@data)
  prior <- 1
  ncols = ncol(data[[1]][[1]]@data)
  
  outputlist <- vector(mode = "list", length = length(seeds))
  
  for (i in 1:length(seeds)){
    tryCatch({
      svMisc::progress(i)
      #datapartition
      train <- data[[i]][[1]]
      test <- data[[i]][[3]]
      trainNB <- cbind(train@data, pa=as.factor(train@pa))
      testNB <- cbind(test@data, pa=as.factor(test@pa))
      
      #copula cors MoM here 
      copula1_C = cor(train@data[which(train@pa==1),], method = "spearman")
      copula1_rhos = corexport(copula1_C)
      copula1_est_rho = as.numeric(2*sin(copula1_rhos*pi/6))
      
      copula0_C = cor(train@data[which(train@pa==0),], method = "spearman")
      copula0_rhos = corexport(copula0_C)
      copula0_est_rho = as.numeric(2*sin(copula0_rhos*pi/6))
      #/copula cors MoM here
      
      #Fit model with adjust=1
      grid.NB = data.frame(usekernel=TRUE,laplace=T,adjust=1)
      
      model.NB = caret::train(pa ~ .,data=trainNB,method="naive_bayes",
                              trControl=trainControl(method="none"),
                              tuneGrid=grid.NB) #class(model.NB$finalModel)
      
      #fit different adjusts
      for (j in 1:length(settings)){
        
        density0 <-  density(trainNB[trainNB$pa==0,j], adjust=settings[[j]]$adjusts0[i])
        density1 <-  density(trainNB[trainNB$pa==1,j], adjust=settings[[j]]$adjusts1[i])
        
        model.NB <- updateMBvar(model=model.NB, variable=varnames[j], densityobject0=density0, densityobject1=density1)
        
      }
      
      mprob <- predict.nonparametric_naive_bayes2(object=model.NB$finalModel, newdata=as.matrix(testNB[,-ncol(testNB)]), prior=c(1, prior), type="prob")[,2]
      
      pu <- getunifPKDE2(data=test@data, naivebayesobject=model.NB, class=1)
      bu <- getunifPKDE2(data=test@data, naivebayesobject=model.NB, class=0)
      
      
      #copula1_density_p <- dCopula(as.matrix(pu[test@pa==1,]), copula=normalCopula(param = copula1_est_rho, dim = ncols, dispstr = "un"))
      #copula1_density_b <- dCopula(as.matrix(bu[test@pa==0,]), copula=normalCopula(param = copula1_est_rho, dim = ncols, dispstr = "un"))
      #copula1_density <- c(copula1_density_p, copula1_density_b)
      copula1_density <- dCopula(as.matrix(pu), copula=normalCopula(param = copula1_est_rho, dim = ncols, dispstr = "un"))
      
      #copula0_density_p <- dCopula(as.matrix(pu[test@pa==1,]), copula=normalCopula(param = copula0_est_rho, dim = ncols, dispstr = "un"))
      #copula0_density_b <- dCopula(as.matrix(bu[test@pa==0,]), copula=normalCopula(param = copula0_est_rho, dim = ncols, dispstr = "un"))
      #copula0_density <- c(copula0_density_p, copula0_density_b)
      copula0_density <- dCopula(as.matrix(bu), copula=normalCopula(param = copula0_est_rho, dim = ncols, dispstr = "un"))
      
      cprob <- copula1_density / copula0_density
      
      
      #outputlist[[i]] <- data.frame(prob=predict.nonparametric_naive_bayes2(object=model.NB$finalModel, newdata=as.matrix(testNB[,1:length(settings)]), prior=c(1, prior), type="prob")[,2], pa=testNB$pa)
      outputlist[[i]] <- data.frame(prob=mprob*cprob, pa=testNB$pa)
    },error=function(e){})
  }
  return(outputlist)
  
}

getunifPKDE2 <- function(data, naivebayesobject, class){
  data <- as.matrix(data)
  output <- matrix(ncol=ncol(data), nrow=nrow(data))
  colnames(output) <- colnames(data)
  
  if(class==1){
    for (i in 1:ncol(data)){
      output[,i] <- approxPKDE2(data[,i], naivebayesobject$finalModel$tables[[i]]$'1')
    }
  }
  
  if(class==0){
    for (i in 1:ncol(data)){
      output[,i] <- approxPKDE2(data[,i], naivebayesobject$finalModel$tables[[i]]$'0')
    }
  }
  output
}

approxPKDE2  <- function(xnew, density){
  x = density$x 
  y = density$y
  dx = diff(x)[1] 
  total = sum(rowMeans(cbind(y[-1], y[-length(y)]))*dx)
  cump = cumsum(rowMeans(cbind(y[-1], y[-length(y)]))*dx)/total
  output <- approx(x, c(0,cump), xnew, rule=2, ties="ordered")$y
  output[output==0] <- 0.0001 #from 0.000001 to 0.001 (0 density issue) #to 0.01 #to 0.0001 and 1 value as well
  output[output==1] <- 0.9999 #1 is also not allowed
  output
}

#
predictMaxBayes.3bgauscop4 <- function(data, settings1, settings0, indep){
  data <- data 
  prior <- length(data@pa[data@pa==1]) / length(data@pa)
  
  outputdf <-as.data.frame(matrix(0, ncol = (3+2*ncol(data@data)), nrow = 1))
  trainNB <- cbind(data@data, pa=as.factor(data@pa))
  
  grid.NB = data.frame(usekernel=TRUE,laplace=T,adjust=1)
  model.NB = caret::train(pa ~ .,data=trainNB,method="naive_bayes",
                          trControl=trainControl(method="none"),
                          tuneGrid=grid.NB) 
  
  
  for (k in 1:ncol(data@data)){
    
    density1 <-  density(trainNB[trainNB$pa==1,k], adjust=settings1[k])
    density0 <-  density(trainNB[trainNB$pa==0,k], adjust=settings0[k])
    
    model.NB <- updateMBvar(model=model.NB, variable=colnames(data@data)[k], densityobject0=density0, densityobject1=density1)
  }
  
  mprob <- predict.nonparametric_naive_bayes2(object=model.NB$finalModel, newdata=as.matrix(trainNB[,-ncol(trainNB)]), prior=c(1, prior), type="prob")[,2]
  
  pu <- getunifPKDE2(data=data@data, naivebayesobject=model.NB, class=1)
  bu <- getunifPKDE2(data=data@data, naivebayesobject=model.NB, class=0)
  
  if(indep==T){
    copula1_density <- dCopula(as.matrix(pu), copula=indepCopula(dim=ncol(data@data)))
    copula0_density <- dCopula(as.matrix(bu), copula=indepCopula(dim=ncol(data@data)))
  }
  
  if(indep==F){
    copula1_C = cor(data@data[which(data@pa==1),], method = "spearman")
    copula1_rhos = corexport(copula1_C)
    copula1_est_rho = as.numeric(2*sin(copula1_rhos*pi/6))
    copula1_density <- dCopula(as.matrix(pu), copula=normalCopula(param = copula1_est_rho, dim = ncol(data@data), dispstr = "un"))

    copula0_C = cor(data@data[which(data@pa==0),], method = "spearman")
    copula0_rhos = corexport(copula0_C)
    copula0_est_rho = as.numeric(2*sin(copula0_rhos*pi/6))
    copula0_density <- dCopula(as.matrix(bu), copula=normalCopula(param = copula0_est_rho, dim = ncol(data@data), dispstr = "un"))
  }
  
  cprob <- copula1_density / copula0_density
  
  output <- data.frame(prob=mprob*cprob, pa=trainNB$pa)
  output$nprob <- output$prob/sum(output$prob)
  
  return(output)
}

getavpar <- function(datas){
  output <- matrix(0,nrow=1, ncol=15)
  for (i in 1:15){
    
    colnames(output) <- c("seed", rep(colnames(data@data),each=2))
    output[1,i] <- mean(datas[,i])
  }
  output <- as.data.frame(output)
  return(output)
}



###
### Case Study 2
###



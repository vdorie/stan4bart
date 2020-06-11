
rstanbartarmglmerglmer_fit <- function(x, z, iter, warmup, samplingArgs, bartControl)
{
  y <- samplingArgs$data$y
  n <- samplingArgs$data$N
  q <- samplingArgs$data$q
  
  bartSampler <- dbarts::dbarts(x, y, resid.prior = fixed(1.0), control = bartControl, sigma = 1.0)
  bartModel <- bartSampler$model
  bartSampler$sampleTreesFromPrior()
  bartSamples <- bartSampler$run()
  
  samplingArgs$data$y <- y - bartSamples$train[,1L]
  y.scale <- diff(range(y))
        
  stanFit <- suppressWarnings(do.call(rstan::sampling, samplingArgs))
  stanSamples <- extract(stanFit)
  
  adaptation_info <- rstan::get_adaptation_info(stanfit) # currently not referenced
  samplingArgs$check_data <- FALSE
  samplingArgs$init <- list(stanSamples)
  
  if (warmup > 0L) for (i in seq_len(warmup)) {
    # use ranef as offset
    bartSampler$setOffset(as.vector(z %*% stanSamples[["b"]][1L,]))
    
    # TODO: figure out something less messy than having to scale by the data when
    #       fixing the resid sd; see todo in dbarts/bartFit.cpp
    
    # fix resid sd to sample provided by Stan
    bartModel@resid.prior@value <- (stanSamples[["aux"]][1L] / y.scale)^2
    bartSampler$setModel(bartModel)
    
    # draw trees
    bartSamples <- bartSampler$run()
    
    # draw ranef using Stan with partial resid as response
    samplingArgs$data$y <- y - bartSamples$train[,1L] 
    stanFit <- suppressWarnings(do.call(rstan::sampling, samplingArgs))
    stanSamples <- extract(stanFit)
    samplingArgs$init <- list(stanSamples)
  }
  
  # the pars for each iter should line-up such that the random effects
  # used when sampling trees are the same; this creates an annoying off-by one situation
  fixef <- matrix(NA_real_, iter, n)
  ranef <- matrix(NA_real_, iter, q, dimnames = list(NULL, b_names(names(stanFit), value = TRUE)))
  sigma <- rep_len(NA_real_, iter)
  
  if (iter > 1L) {
    for (i in seq_len(iter - 1L)) {
      # use ranef as offset
      bartSampler$setOffset(as.vector(z %*% stanSamples[["b"]][1L,]))
      
      # fix resid sd to sample provided by Stan
      bartModel@resid.prior@value <- (stanSamples[["aux"]][1L] / y.scale)^2
      bartSampler$setModel(bartModel)
      
      # draw trees
      bartSamples <- bartSampler$run()
      
      # store current samples such that trees are matched with their ranefs
      fixef[i,] <- bartSamples$train[,1L]
      ranef[i,] <- stanSamples[["b"]][1L,]
      sigma[i]  <- stanSamples[["aux"]][1L]
      
      # draw ranef using Stan with partial resid as response
      samplingArgs$data$y <- y - bartSamples$train[,1L]
      stanFit <- suppressWarnings(do.call(sampling, samplingArgs))
      stanSamples <- extract(stanFit)
      samplingArgs$init <- list(stanSamples)
    }
    
    bartSampler$setOffset(as.vector(z %*% stanSamples[["b"]][1L,]))
    bartModel@resid.prior@value <- (stanSamples[["aux"]][1L] / y.scale)^2
    bartSampler$setModel(bartModel)
    bartSamples <- bartSampler$run()
     
    fixef[iter,] <- bartSamples$train[,1L]
    ranef[iter,] <- stanSamples[["b"]][1L,]
    sigma[iter]  <- stanSamples[["aux"]][1L]
  }
  
  nlist(fixef, ranef, sigma)
}


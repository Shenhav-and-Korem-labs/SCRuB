rescale <- function(m) m / sum(m)
get_rescaled_mat <- function(Q) ( Q %>% apply(MARGIN=1, rescale) %>% t() )

rarefy <- function(x,maxdepth=10000){
  nr <- nrow(x)
  nc <- ncol(x)
  
  for(i in 1:nrow(x)){
    if(sum(x[i,]) > maxdepth){
      prev.warn <- options()$warn
      options(warn=-1)
      s <- sample(nc, size=maxdepth, prob=x[i,], replace=T)
      options(warn=prev.warn)
      x[i,] <- hist(s,breaks=seq(.5,nc+.5,1), plot=FALSE)$counts
    }
  }
  return(x)
}

LSQ_init <- function(sinks, sources, COVERAGE=10000){
  n_samples <- nrow(sinks)
  n_controls <- nrow(sources)
  alpha_inits <- matrix(data = rep(0, ( n_controls + 1)  * n_samples),
                        nrow = n_samples,
                        ncol = (n_controls + 1) )
  S_unk_obs = sinks * 0
  
  for( i in (1:n_samples) ){
    rare_sink <- rescale( rarefy( t( sinks[i,] ), COVERAGE ) )  * COVERAGE
    inits <- lsq_procedure(
      sources,
      rare_sink
    )
    
    alpha_inits[i, ] <- inits$alpha
    S_unk_obs[i, ] <- rescale( inits$unknown )  * COVERAGE
  }
  
  return(list(
    alpha_inits=alpha_inits,
    S_unk_obs = S_unk_obs ) )
}



initialize_contaminant <- function(sinks, control_samples, COVERAGE=10000){
  # goal - initialize contaminant as weighted average of control samples,
  # weight is determined by alphas from lsq_init
  n_controls <- nrow(control_samples)
  lsq_out <- LSQ_init(sinks, control_samples, COVERAGE=10000)
  alpha_inits <- lsq_out$alpha_inits
  contam_init <-  rescale( (alpha_inits %>% apply(MARGIN=2, mean) )[1:n_controls] ) %*% control_samples
  return(contam_init)
}

run_linear_setup <- function(samples, source, lsq_out, COVERAGE){
  n_sources <- nrow(source)
  # this function takes as input as lsq output,
  # and organizes the info into a list
  gamma_known <- apply( source , 1, rescale ) %>% t()
  gamma_known[gamma_known %>% is.na()] = 0
  
  #getting one gamma_unknown for each row of X
  gamma_unknowns = apply(lsq_out$S_unk_obs, MARGIN = 1, function(x) rescale( (x) ) ) %>% t()
  
  #Giving S_unk_obs the normalized coverage
  S_unks <- gamma_unknowns * COVERAGE
  
  #making the gammas into a list, to be compatiable with FEAST's do_EM function
  gamma_known_list <- lapply(seq_len(nrow(gamma_known)), function(i) as.matrix( gamma_known[i,] ) %>% t() )
  
  #making the S_known into a list
  S_known_list <- lapply(seq_len(n_sources), function(i) as.matrix(source[i,]) %>% t() )
  
  return(list(samples=samples,
              torets = lsq_out$alpha_inits,
              sources_known = gamma_known_list,
              sources_unk = gamma_unknowns,
              observed_known = S_known_list,
              observed_unknown = S_unks)
  )
}


set_up_SCRUB <- function(samples,
                         control_samples,
                         use_vi = FALSE,
                         COVERAGE = 10000
                         ){
  

  contam <- initialize_contaminant( samples, control_samples, COVERAGE)
  
  n_samples <- nrow(samples)
  
  sink_lsq <- LSQ_init(samples, contam, COVERAGE)
  
  sink_setup <-  run_linear_setup(samples, contam, sink_lsq, COVERAGE)
  
  return(list(sink_setup = sink_setup,
              contam=contam )
  )
}



update_scheme <- function(X, gb, G, p, L){
  R <- (p*G) / ( p*G + ( (1-p) %*% t(gb) ) )
  R[is.na(R)] <- 0
  new_gamma <- ( X * R ) %>% get_rescaled_mat()
  new_gb <- rescale( ( X * ( 1-R ) ) %>% apply(MARGIN=2, sum) + L %>% apply(MARGIN=2, sum) )
  new_p <- apply( ( X * R ) / apply(X, MARGIN=1, sum), MARGIN=1, sum )
  
  return(list(new_gamma=new_gamma,
              new_gb=new_gb,
              new_p=new_p))
}


log_likelihood_SCRUB <- function(X, Y, G, gb, p){
  t1 <- X * log( p*G + ( (1-p) %*% t(gb) ) )
  t2 <- Y * log( gb )
  return( sum(t1*is.finite(t1) , na.rm=T) + sum(t2*is.finite(t2) , na.rm=T) )
}

#' Decontaminate a set of samples
#' SCRuB removes contamination from an inputted set of samples using microbial source-tracking techniques.
#' @param samples n_samples x n_taxa -- a count matrix representing the observed reads from samples to be decontaminated
#' @param controls n_controls x n_taxa -- a count matrix representing the observed reads from controls whose contents represent the contamination community to be removed from `samples`
#' @param COVERAGE 10000, integer representing the normalizing variable to be used for SCRuB's initialization scheme.
#' @param print_loglikelihood Boolean, TRUE of FALSE. Determines if SCRuB should print the calculated log-likelihood during each iteration
#' @return A list containing:
#' 1) decontaminated_samples - a n_samples x n_taxa count matrix, representing the decontaminated samples
#' 2) p - The fitted p parameter, as described in SCRuB's methods.
#' An n_sample vector representing the estimate proportion of each observe sample that was not contamination
#' A dataset that had no contamination would have a p of 1s, while a dataset of entirely contamination would have a p of 0
#' 3) gamma - the $\gamma$ parameter described in SCRuB's methods. An n_taxa vector representing the estimated relative abundance of the contamination community
#' 4) loglikelihood - float. The log-likelihood of the inputted dataset based on SCRuB's fitted parameters.
#' @export
SCRUB_no_spatial <- function(samples,
                             controls,
                             COVERAGE=10000,
                             print_loglikelihood = F
){
  
  
  setup <- set_up_SCRUB(samples,
                        controls)
  
  
  f_setup <- setup$sink_setup
  contam <- setup$contam
  
  X <- f_setup$samples
  gb <- contam[1, ] %>% rescale
  G <- f_setup$sources_unk
  p <- f_setup$torets[,1]
  L <- controls
  
  eps <- 1e-7
  patience <- 0
  max_iters <- 500
  a <- 0
  while(a<max_iters){
    out <- update_scheme(X, gb, G, p, L)
    gb <- out$new_gb
    diff <- abs( G-out$new_gamma ) %>% mean()
    # early stopping criteria
    if(diff<eps){
      patience <- patience+1
      if(patience>5) a <- max_iters + 1
    }else{
      patience <- 0
      a <- a+1
    }
    
    G <- out$new_gamma
    p <- out$new_p
    
    if(print_loglikelihood)  print( paste('Log-Likelihood of',
                                          log_likelihood_SCRUB(X, L, G, gb, p ),
                                          'after', a, 'iterations!') )
  }
  
  return( list(p = p,
               decontaminated_samples = G  %>%
                 sweep(MARGIN=1, samples %>% apply(MARGIN = 1, sum) %>% as.vector(), `*` ) %>% round(),
               gamma = gb,
               loglikelihood=log_likelihood_SCRUB(X, L, G, gb, p )
  ))
}


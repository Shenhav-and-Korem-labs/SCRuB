
## Spatial SCRUB Functions

spatial_SCRUB_EM_update <- function(X, Y, gb, G, p, alpha, n, m){
  ww <- torch_unsqueeze(p, -1)*G # making this it's own line so it doesn't need to be calculated twice
  R <- ww$divide(ww + torch_unsqueeze(1-p, -1)*torch_unsqueeze(gb, 1) )
  R[torch_isnan(R)] <- 0
  
  Q <- nnf_normalize( torch_unsqueeze(alpha, -1) *
                        torch_unsqueeze( torch_vstack( c(G, torch_unsqueeze(gb, 1) ) ), 1 ) , 
                      p=1, dim=2 )
  
  G <- nnf_normalize( X*R + torch_sum( torch_unsqueeze(Y, 2) * Q[,1:n,], dim=1) , 
                      p=1, dim=2 )
  
  gb <- ( torch_sum( X * (1-R), dim=1) +  torch_sum( Y*Q[,n+1,], dim=1) ) %>%
    nnf_normalize(p=1, dim=1)
  
  p <- torch_sum( X*R / torch_sum(X, dim=2, keepdim = T) , dim=2) 
  
  alpha <- torch_divide(  ( torch_unsqueeze(Y, dim=2) * Q )$sum(dim=3) ,
                          torch_sum(Y, dim=2, keepdim = T) )
  return(list(G = G, 
              gb=gb,
              p=p, 
              alpha=alpha))
}

log_likelihood_spatial_SCRUB <- function(X, Y, G, gb, p, alpha){
  t1 <- as.matrix( X * torch_log( p$unsqueeze( -1)*G + torch_unsqueeze(1-p, -1)*torch_unsqueeze(gb, 1)  ) )
  t2 <- as.matrix(Y) * log( as.matrix(alpha) %*% as.matrix( torch_vstack( c(G, torch_unsqueeze(gb, 1)  ) ) ) )
  return( sum(t1*is.finite(t1) , na.rm=T) + sum(t2*is.finite(t2) , na.rm=T) )
}

initialize_w2w_params <- function(samples, 
                                  controls, 
                                  cont_samp_dists, 
                                  dist_threshold=1){
  n <- nrow(samples)
  m <- ncol(samples)
  l <- nrow(controls)
  alpha <- matrix(0, nrow=l, ncol=n+1 )
  row.names(alpha) <-  row.names(controls)
  colnames(alpha) <- c(row.names(samples), 'background')
  
  cont_tmp = controls * 0 
  for(k in row.names(controls)){
    nearby_samps <- which( cont_samp_dists[k,] <= dist_threshold ) %>% names()
    
    if(is_empty(nearby_samps)){
      alpha[k, 'background'] <- 1
      cont_tmp[k, ] <-  controls[k, ]
    }else{
      well_conts <- samples[nearby_samps, ]
      
      if( length(nearby_samps)==1 ) well_conts %<>% t()
      
      inits <- lsq_procedure(
        well_conts,
        t(controls[k,])
      )
      
      alpha[k, c(nearby_samps, 'background')] <- inits$alpha 
      cont_tmp[k, ] <-  inits$unknown 
    }
  }
  
  return(list(alpha=alpha, 
              cont_tmp=cont_tmp
  ))
}

#' Decontaminate a set of samples
#' SCRuB removes contamination from an inputted set of samples using microbial source-tracking techniques. 
#' @param data ( n_samples + n_controls ) x n_taxa -- a count matrix representing the observed reads from the controls, and the samples to be decontaminated
#' @param is_control boolean vecotr, length n_samples+n_controls -- a booolean vector indicating which rows of the data matrix correspond to negative controls, whose contents represent the contamination community to be removed from the non-control samples
#' @param well_dists numeric matrix. Deontes the spatial distance between each pairwise samples, names rows and columns must correst to the rownames of `data`
#' @param dist_threshold float - Determines the maximum euclidean distance between samples and controls which SCRuB determines as potential sources of well leakage. Default of 1.5 
#' @param print_loglikelihood Boolean, TRUE of FALSE. Determines if SCRuB should print the calculated log-likelihood during each iteration
#' @return A list containing:
#' 1) decontaminated_samples - a n_samples x n_taxa count matrix, representing the decontaminated samples
#' 2) p - The fitted p parameter, as described in SCRuB's methods. 
#' An n_sample vector representing the estimate proportion of each observe sample that was not contamination
#' A dataset that had no contamination would have a p of 1s, while a dataset of entirely contamination would have a p of 0
#' 3) alpha - The fitted alpha parameter, as described in SCRuB's methods. 
#' An n_control x ( n_sample + 1 ) matrix, representing the estimated contribution of the contaminant and each sample to each control, where the (n_sample + 1)th column represents the contribution from the contamination to the control.
#' Each row of alpha sums to 1, with each entry of the (n_sample + 1)th  column being 1 means there is zero estimated well leakge, while entries close to zero would indicate there is a high level of well leakage
#' 4) gamma - the gamma parameter described in SCRuB's methods. An n_taxa vector representing the estimated relative abundance of the contamination community
#' 5) loglikelihood - float. The log-likelihood of the inputted dataset based on SCRuB's fitted parameters.
#' @export
spatial_SCRUB <- function(data, 
                          is_control, 
                          well_dists,
                          dist_threshold=1.5, 
                          print_loglikelihood=F
                          ){
  
  
  samples <- data[is_control == FALSE, ]
  if(sum(is_control==F)==1){ 
      samples <- t(samples)
      row.names(samples) <- c(row.names(data)[is_control==F])
  }
  controls <- data[is_control, ]
  
  if(sum(is_control)==1){ 
    controls <- controls %>% t()
    row.names(controls) <- row.names(data)[which(is_control)]
    cont_samp_dists <- well_dists[row.names(controls),row.names(samples)] %>% t()
    row.names(cont_samp_dists) <- row.names(controls)
  }else{
    cont_samp_dists <- well_dists[row.names(controls),row.names(samples)]
    if(nrow(samples)==1){ 
        cont_samp_dists <- t(cont_samp_dists) %>% t()
        colnames(cont_samp_dists) <-row.names(samples)
        }
  }

  w2w_inits <- initialize_w2w_params(samples,
                                     controls, 
                                     cont_samp_dists, 
                                     dist_threshold)
  alpha <- w2w_inits$alpha
  tmp_controls <- w2w_inits$cont_tmp
  if(sum(tmp_controls) == 0 ) print("Controls' LSQ init removed everything from a control -- will instead set initialization to the observed control")
  if(sum(tmp_controls) == 0 )  tmp_controls <- controls
  
  setup <- set_up_SCRUB(samples, 
                           tmp_controls) 
  
  f_setup <- setup$sink_setup
  contam <- setup$contam
  
  # Initiallize all parameters
  X <- torch_tensor( f_setup$samples, requires_grad = FALSE) 
  gb <- torch_tensor( contam[1, ] %>% rescale, requires_grad = FALSE)
  G <- torch_tensor( f_setup$sources_unk, requires_grad = FALSE)
  p <- torch_tensor( f_setup$torets[,1], requires_grad = FALSE)
  Y <- torch_tensor( controls, requires_grad = FALSE)
  
  n <- nrow(X)
  m <- ncol(X)
  l <- nrow(controls)
  
  alpha <- torch_tensor(alpha, requires_grad = FALSE)
  
  eps <- 1e-7
  patience <- 0
  max_iters <- 500
  a <- 0
  while(a<max_iters){
    out <- spatial_SCRUB_EM_update(X, Y, gb, G, p, alpha, n, m)
    diff <- torch_dist(G, out$G, p=1)$item()
    
    # early stopping criteria
    if(diff<eps){
      patience <- patience+1
      if(patience>5) a <- max_iters + 1
    }else{
      patience <- 0
      a <- a+1
    }
    gb <- out$gb
    G <- out$G
    p <- out$p
    alpha <- out$alpha
    
    if( print_loglikelihood ) print( paste('Log-Likelihood of', 
                                           log_likelihood_spatial_SCRUB(X, Y, G, gb, p, alpha), 
                                           'after', a, 'iterations!') )
    
  }
  return( list( 
              decontaminated_samples = as.matrix(G) %>% 
                  sweep(MARGIN=1, samples %>% apply(MARGIN = 1, sum) %>% as.vector(), `*` ) %>% round(), 
               p = p %>% as.matrix() %>% c(),
               alpha=as.matrix(alpha),
               gamma = gb %>% as.matrix() %>% c(), 
               loglikelihood=log_likelihood_spatial_SCRUB(X, Y, G, gb, p, alpha)
  ))
}

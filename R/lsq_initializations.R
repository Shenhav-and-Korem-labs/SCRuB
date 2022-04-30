
get_rowwise_rmses <- function(mat1, mat2){
  ans <-  c()
  for(i in 1:nrow(mat1)) ans <- c(ans, get_rmse(mat1[i, ], mat2[i, ]))
  return(ans)
}

lsq_glmnet_l1l2 <- function(
                            sink, 
                            sources,
                            l1l2=1, 
                            lambda=1e-6,
                            normalize=T
                            ){
  
  if (normalize) {
    sources <- sources / rowSums(sources)
    sink <- sink / sum(sink)
  }
  
  Amat <- model.matrix(~Sources, list(Sources=t(sources)))
  result <- glmnet(Amat, sink, alpha=l1l2, lambda=lambda, lower.limits=0)
  
  contr <- result$beta[2:length(result$beta)]
  
  return (contr)
}

lsq_procedure <- function(sources, sink, unknown_eps=0.01) {
  weights_l1 <- lsq_glmnet_l1l2(
    sink,
    sources,
    normalize=F,
    l1l2=1)
  
  # normalizing
  if(sum(weights_l1) != 0){
    lsq_init <- weights_l1 / sum(weights_l1)
  } else lsq_init <- weights_l1
  
  # 1. The alpha init
  alpha_init <- c(lsq_init, 0)
  alpha_init[alpha_init < unknown_eps] <- 0.01/length(alpha_init)
  alpha_init <- alpha_init / sum(alpha_init)
  
  # 2. The unknown init
  max_source <- sources[which.max(lsq_init),]
  
  # We scale the max source according to its given weight
  mixed_max <- max_source * max(lsq_init)
  
  
  
  unknown_init <- sink - as.vector( mixed_max )
  
  unknown_init[unknown_init < 0] <- 0   # clip at zero
  
  return (list(
    lsq=weights_l1,
    alpha=alpha_init,
    unknown=unknown_init
  ))
}





SCRuB_wrapper_no_spatial <- function(data, control_idcs, verbose=F){
  
  any_cont_type <- ( control_idcs %>% rowSums() ) > 0
  samples <- data[any_cont_type==F, ]
  n_smps <- nrow(samples)
  
  inner_scrub_iterations <- list()
  for(i in 1:ncol(control_idcs) ){
    print(paste( 'SCRuBbing away contamination in the',  colnames(control_idcs)[i] , 'controls...') )
    if(i==1){
        if(sum(control_idcs[, i])==1){
          cont_tmp <- t( data[control_idcs[,i],] )
          row.names(cont_tmp) <- row.names(data)[control_idcs[,i]]
        }else{
          cont_tmp <- data[control_idcs[, i], ]
        }
      
      
      inner_scrub_iterations[[ colnames(control_idcs)[i] ]] <- SCRUB_no_spatial( samples, 
                                                                                 cont_tmp,
                                                                                 print_loglikelihood = verbose
                                                                                 )
      
      cumulative_p <- inner_scrub_iterations[[ colnames(control_idcs)[i] ]]$p
    }else{
      
      if(sum(control_idcs[, i])==1){
        cont_tmp <- t( data[control_idcs[,i],] )
        row.names(cont_tmp) <- row.names(data)[control_idcs[,i]]
      }else{
        cont_tmp <- data[control_idcs[, i], ]
      }
      
      inner_scrub_iterations[[ colnames(control_idcs)[i] ]] <- SCRUB_no_spatial( 
        inner_scrub_iterations[[ colnames(control_idcs)[i-1] ]]$decontaminated_samples,
        cont_tmp,
        print_loglikelihood = verbose)
      
      cumulative_p <- inner_scrub_iterations[[ colnames(control_idcs)[i] ]]$p * cumulative_p
    }
    
    row.names( inner_scrub_iterations[[ colnames(control_idcs)[i] ]]$decontaminated_samples ) <- row.names(samples)
  }
  
  return(list(decontaminated_samples=inner_scrub_iterations[[colnames(control_idcs)[ncol(control_idcs)] ]]$decontaminated_samples,
              p = cumulative_p, 
              inner_iterations=inner_scrub_iterations
  )
  )
}


SCRuB_wrapper <- function(data, control_idcs, well_dists, dist_threshold=1.5, verbose=F){
  
  any_cont_type <- ( control_idcs %>% rowSums() ) > 0
  samples <- data[any_cont_type==F, ]
  n_smps <- nrow(samples)
  
  inner_scrub_iterations <- list()
  for( i in 1:ncol(control_idcs) ){
    print(paste( 'SCRuBbing away contamination in the',  colnames(control_idcs)[i] , 'controls...') )
    if(i==1){
      if(sum(control_idcs[, i])==1){
          cont_tmp <- t( data[control_idcs[,i],] )
          row.names(cont_tmp) <- row.names(data)[control_idcs[,i]]
      }else{
        cont_tmp <- data[control_idcs[, i], ]
      }


      inner_scrub_iterations[[ colnames(control_idcs)[i] ]] <- spatial_SCRUB( rbind( samples, cont_tmp ),
                                                                              c( rep(F, nrow(samples)), rep(T, sum(control_idcs[,i]) ) ),
                                                                              well_dists,
                                                                              dist_threshold, 
                                                                              print_loglikelihood = verbose 
                                                                              )
      colnames(inner_scrub_iterations[[ colnames(control_idcs)[i] ]]$decontaminated_samples) <- colnames(data)
      cumulative_p <- inner_scrub_iterations[[ colnames(control_idcs)[i] ]]$p
    }else{
          if(sum(control_idcs[, i])==1){
            cont_tmp <- t( data[control_idcs[,i],] )
            row.names(cont_tmp) <- row.names(data)[control_idcs[,i]]
          }else{
            cont_tmp <- data[control_idcs[, i], ]
          }
          
      
      inner_scrub_iterations[[ colnames(control_idcs)[i] ]] <- spatial_SCRUB( 
                                      rbind( inner_scrub_iterations[[ colnames(control_idcs)[i-1] ]]$decontaminated_samples,
                                             cont_tmp ) , 
                                      c( rep(F, nrow(samples)), rep(T, sum(control_idcs[,i]) ) ),
                                      well_dists,
                                      dist_threshold, 
                                      print_loglikelihood = verbose
                                    )
      colnames(inner_scrub_iterations[[ colnames(control_idcs)[i] ]]$decontaminated_samples) <- colnames(data)
      cumulative_p <- inner_scrub_iterations[[ colnames(control_idcs)[i] ]]$p * cumulative_p
    }
    
    row.names( inner_scrub_iterations[[ colnames(control_idcs)[i] ]]$decontaminated_samples ) <- row.names(samples)
  }
  
  return(list(decontaminated_samples=inner_scrub_iterations[[colnames(control_idcs)[ncol(control_idcs)] ]]$decontaminated_samples,
              p = cumulative_p, 
              inner_iterations=inner_scrub_iterations
  )
  )
}


#' Decontaminate a set of samples
#' SCRuB removes contamination from an inputted set of samples using microbial source-tracking techniques. 
#' @param data ( n_samples + n_controls ) x n_taxa -- a count matrix representing the observed reads from the controls, and the samples to be decontaminated
#' @param metadata a metadata matrix of 2 columns (and a third optional column), where each row name is aligned to the `data` parameter row name. 
#' The first metadata column is boolean, denoting `TRUE` if the corresponding `data` row belongs to a sample representing a contamination source, and `FALSE` otherwise. 
#' The second column is a string, identifying the type of each sample, such that SCRuB can identifying which control samples should be grouped together. 
#' The third (and optional, but highly recommended) column is a string entry identifying the well location of the corresponding sample, which allows SCRuB to track well leakage.
#' This must be in a standard *LETTER**NUMBER* format, i.e. A3, B12, D4...
#' @param control_order vector, default NA. If specified, outlines the order of controls to run decontamination in. Input as a vector, of which each element must also be found in the metadata's second column. If not specified, all control types found in `metadata` will be run sequentially based on their order from that table. 
#' @param dist_threshold float - Determines the maximum euclidean distance between samples and controls which SCRuB determines as potential sources of well leakage. This input is only used if the well location metadata is provided Default of 1.5 
#' @param dist_metric string, default `euclidean`. The distance metric to be used when evaluating samples' physical distance. This input is used in the `stats` library's `dist` function; see their documentation for other options.
#' @param verbose boolean - if TRUE, SCRuB prints the log-likelihood of hte dataset thoughout each iteration. 
#' @return A list containing:
#' 1) decontaminated_samples - a n_samples x n_taxa count matrix, representing the decontaminated samples
#' 2) p - The fitted p parameter, as described in SCRuB's methods. 
#' An n_sample vector representing the estimate proportion of each observe sample that was not contamination
#' A dataset that had no contamination would have a p of 1s, while a dataset of entirely contamination would have a p of 0
#' 3) inner_iterations -- results from SCRuB's intermediary steps, see the `Spatial_SCRUB` and `SCRUB_no_spatial` documentation for more information
#' Additional imports
#' @import tidyverse
#' @import glmnet
#' @import torch
#' @export
SCRuB <- function(data, 
                  metadata, 
                  control_order=NA,
                  dist_threshold=1.5, 
                  dist_metric='euclidean',
                  verbose=F
                  ){
  
  if( ( row.names(data) == row.names(metadata) ) %>% mean() < 1 ){
    stop("The row names of the `data` and `metadata` inputs must be equivalent!")
  }
  
  data <- as.matrix(data)
  metadata <- as.data.frame(metadata)
  
  ## checking the specified controls/control order
  # if(is.character(control_order)){
  #   print('here')
  #   control_order <- c(control_order)
  # }

  if(is.vector(control_order)&(max(is.na(control_order))==F)){
    if(length(control_order) !=length(unique(control_order))) stop('All entries within `control_order` must be unique!')
    for(cont_tp in control_order){
      if( F == (cont_tp %in% metadata[metadata[,1]==T, 2 ] %>% unique() ) ){
        stop(paste0( "All entries in `control_order` must be found in the `metadata` rows where `is_contaminant` is TRUE. No such entries found for `", 
                      cont_tp,
                      "`."  ) )
        }   }
  } else if (min(is.na(control_order))){
    ## sleect control order based on the oder in the metadata
    control_order <- metadata[metadata[,1]==T, 2 ] %>% unique() %>% as.character()
  }
  
  control_mat <- control_order %>%
    sapply( function(x) metadata[,2] == x ) %>% as.matrix()
  
  
  if(ncol(metadata) == 3){
    well_dists <- metadata %>%
      mutate(well = metadata[, 3] %>% sapply( function(x) which( LETTERS == substr(x, 1, 1) ) ),
             indices = metadata[, 3] %>% sapply( function(x) substr(x, 2, 3) %>% as.integer)
      ) %>% 
      select(well,indices) %>% 
      dist(method=dist_metric) %>% as.matrix()
    print('Incorporating the well metadata to track well-to-well leakage!' )
    
    return( SCRuB_wrapper(data, control_mat, well_dists, dist_threshold=dist_threshold, verbose = verbose) )
    
  }else{
    print('Did not find well metadata, running SCRuB without the spatial component')
    return(SCRuB_wrapper_no_spatial(data = data, control_mat, verbose = verbose))
  }
}




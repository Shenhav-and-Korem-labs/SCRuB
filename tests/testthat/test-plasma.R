# 
# data <- read.csv(here::here('tutorial/plasma_data.csv'), row.names=1) %>% as.matrix()
# metadata <- read.csv(here::here('tutorial/plasma_metadata.csv'), row.names=1)

context('test-plasma.R')

test_that(desc = 'testing SCRuB variations on plasma dataset',
          {
            n_feats_considered=100 # to avoid memory errors on GHA
            data <- (read.csv( paste0( test_path(), '/plasma_data.csv'), row.names=1) %>% as.matrix() )[,1:n_feats_considered]
            metadata <- read.csv( paste0( test_path(), '/plasma_metadata.csv'), row.names=1)
            
            set.seed(1)

            scr_out_1 <- SCRuB( data[, 1:n_feats_considered], metadata, c("control blank DNA extraction", "control blank library prep") )
            scr_out_2 <- SCRuB( data[, 1:n_feats_considered], metadata[,1:2], "control blank library prep") 
            scr_out_2 <- SCRuB( data[, 1:n_feats_considered], metadata[,1:2],  c("control blank DNA extraction", "control blank library prep") ) 
           
            print('Testing shortened data')
            scr_out_3 <- SCRuB( 100*data[c(1:50, 81, 54), 1:n_feats_considered], metadata[c(1:50, 81, 54),], verbose=T) 
            scr_out_3 <- SCRuB( data[c(1:50, 81, 54), 1:n_feats_considered], metadata[c(1:50, 81, 54),], verbose=T) 
            scr_out_3 <- SCRuB( 100*data[c(1:50, 80), 1:n_feats_considered], metadata[c(1:50, 80),1:2], verbose=T) 
            
            
            expect_type(scr_out_1, 'list')
            expect_type(scr_out_1$inner_iterations,  'list' )
            
            expect_true( nrow(scr_out_1$decontaminated_samples) == sum( F == metadata$is_control ) )
            expect_true( nrow(scr_out_1$decontaminated_samples) == length(scr_out_1$p))
            expect_true( length(scr_out_1$inner_iterations$`control blank DNA extraction`$gamma) == n_feats_considered )
            expect_true( sum(scr_out_1$inner_iterations$`control blank DNA extraction`$gamma) %>% round(5) == 1 )
            expect_true( scr_out_1$inner_iterations$`control blank DNA extraction`$alpha %>% rowSums() %>% mean() %>% round(5) == 1 )
            expect_true( min( colnames(data)[1:n_feats_considered] == colnames(scr_out_1$decontaminated_samples)) == T )
            expect_true( min( colnames(data)[1:n_feats_considered] == colnames(scr_out_2$decontaminated_samples)) == T )
            
            ## error checks
            expect_error( SCRuB( data[1:50, 1:n_feats_considered], metadata[51:100,]) )
            expect_error( SCRuB( data[1:50, 1:n_feats_considered], metadata[1:50,],  rep("control blank library prep", 2)) )
            expect_error( SCRuB( data[1:50, 1:n_feats_considered], metadata[1:50,],  c("control blank library prep", 'fake control')) )
            metadata$sample_well[1] <- 'B16'
            expect_error( SCRuB( data[1:50, 1:n_feats_considered], metadata[1:50,] ) )
            expect_error( SCRuB( data[1:50, 1:n_feats_considered], metadata[1:50,], NA ) )
            expect_error( SCRuB( data[1:50, 1:n_feats_considered], metadata[1:50,], c('control blank library prep', NA ) ) )
          })




data <- read.csv(here::here('tutorial/plasma_data.csv'), row.names=1) %>% as.matrix()
metadata <- read.csv(here::here('tutorial/plasma_metadata.csv'), row.names=1)

set.seed(1)
scr_out_1 <- SCRuB( data[, 1:500], metadata )
scr_out_2 <- SCRuB( data[, 1:500], metadata[,1:2] )
scr_out_3 <- SCRuB( data[1:50, 1:500], metadata[1:50,] )
scr_out_4 <- SCRuB( data[1:50, 1:500], metadata[1:50,1:3] )
scr_out_5 <- SCRuB( data[1:50, 1:500], metadata[1:50,1:2] )


expect_type(scr_out_1, 'list')
expect_type(scr_out_1$inner_iterations,  'list' )

expect_true( nrow(scr_out_1$decontaminated_samples) == sum( F == metadata$is_control ) )
expect_true( nrow(scr_out_1$decontaminated_samples) == length(scr_out_1$p))
expect_true( length(scr_out_1$inner_iterations$`control blank DNA extraction`$gamma) == 500 )
expect_true( sum(scr_out_1$inner_iterations$`control blank DNA extraction`$gamma) %>% round(5) == 1 )
expect_true( scr_out_1$inner_iterations$`control blank DNA extraction`$alpha %>% rowSums() %>% mean() %>% round(5) == 1 )


data <- read.csv(here::here('tutorial/plasma_data.csv'), row.names=1) %>% as.matrix()
metadata <- read.csv(here::here('tutorial/plasma_metadata.csv'), row.names=1)

set.seed(1)
print('Testing full samples, metadata....')
scr_out_1 <- SCRuB( data[, 1:500], metadata, c("control blank DNA extraction", "control blank library prep") )
print('Testing metadata without spatial information, only control blank library prep')
scr_out_2 <- SCRuB( data[, 1:500], metadata[,1:2], "control blank library prep") 
print('Testing shortened data')
scr_out_3 <- SCRuB( data[1:50, 1:500], metadata[1:50,], c("control blank DNA extraction") )


expect_type(scr_out_1, 'list')
expect_type(scr_out_1$inner_iterations,  'list' )

expect_true( nrow(scr_out_1$decontaminated_samples) == sum( F == metadata$is_control ) )
expect_true( nrow(scr_out_1$decontaminated_samples) == length(scr_out_1$p))
expect_true( length(scr_out_1$inner_iterations$`control blank DNA extraction`$gamma) == 500 )
expect_true( sum(scr_out_1$inner_iterations$`control blank DNA extraction`$gamma) %>% round(5) == 1 )
expect_true( scr_out_1$inner_iterations$`control blank DNA extraction`$alpha %>% rowSums() %>% mean() %>% round(5) == 1 )
expect_true( min( colnames(data)[1:500] == colnames(scr_out_1$decontaminated_samples)) == T )
expect_true( min( colnames(data)[1:500] == colnames(scr_out_2$decontaminated_samples)) == T )



# Source-tracking for Contamination Removal in *micro*Biomes (SCRuB)

SCRuB is a tool designed to help researchers address the common issue of contamination in microbial studies. This package provides an easy to use framework to apply SCRuB to your projects. All you need to get started are n samples x m taxa count matrices for both your samples and controls. 


Support
-----------------------

For support using SCRuB, please use our <a href="https://github.com/korem-lab/SCRuB/issues">issues page</a> or email: gia2105@columbia.edu.


Software Requirements and dependencies
-----------------------

*SCRuB* is implemented in R (>= 3.6.3) and requires the following dependencies: **glmnet**, **torch**, **tidyverse**. Please install and load them prior to trying to install *SCRuB*. 

```
install.packages( c('glmnet', 'torch', 'tidyverse') )
```


Installation
---------------------------

*SCRuB* will be available on QIIME 2 very soon. Until then you can you can simply install *SCRuB* using **devtools**: 
```
devtools::install_github("korem-lab/SCRuB")
```


Usage
___________________

### As input, *SCRuB* takes arguments:

 _data_ - An ( n_samples + n_controls ) x n_taxa count matrix
 
 _is_control_ - A  (n_samples+n_controls) length boolean vector, indicating which rows of the data matrix correspond to negative controls, whose contents represent the contamination community to be removed from the non-control samples
 
(optional) _well_dist_ - An ( n_samples + n_controls ) x ( n_samples + n_controls ) distance matrix, summarizing the pairwise distance between each sample. Both the row names and column names of _well_dist_ must correspond to the row names of the _data_ matrix. See our <a href="https://korem-lab.github.io/SCRuB/tutorial.html">tutorial</a> for a demonstration of how to create this input, using Euclidean distance. 

(optional) _dist_threshold_ float - Determines the maximum distance between samples and controls which SCRuB determines as potential sources of well leakage. Default of 1.5 

(optional) a_init_ float \in (0,1) - The prior assumption representing (1 - level of well leakage into each control). Default is 0.99, i.e. 1% of reads in controls are the result of leakge. 

(optional) _print_loglikelihood_ Boolean, TRUE of FALSE. Determines if SCRuB should print the calculated log-likelihood during each iteration

#### As output, *SCRuB* returns a list containing:

 decontaminated_samples - a n_samples x n_taxa count matrix, representing the decontaminated samples
 
 p - The fitted p parameter, as described in SCRuB's methods. An n_sample vector representing the estimate proportion of each observe sample that was not contamination. A dataset that had no contamination would have a p of 1s, while a dataset of entirely contamination would have a p of 0
 
 alpha - The fitted \alpha parameter, as described in SCRuB's methods. An n_control x ( n_sample + 1 ) matrix, representing the estimated contribution of the contaminant and each sample to each control, where the (n_sample + 1)th column represents the contribution from the contamination to the control. Each row of alpha sums to 1, with each entry of the (n_sample + 1)th  column being 1 means there is zero estimated well leakage, while entries close to zero would indicate there is a high level of well leakage
 
gamma - the $\gamma$ parameter described in SCRuB's methods. An n_taxa vector representing the estimated relative abundance of the contamination community
loglikelihood - float. The log-likelihood of the inputted dataset based on SCRuB's fitted parameters.


Demo
-----------------------
We provide a dataset for an example of *SCRuB* usage. Download the demo files <a href="https://github.com/korem-lab/SCRuB/tree/gh-pages_tmp/tutorial_data">here</a>. We provide an rmarkdown notebook <a href="https://github.com/korem-lab/SCRuB/blob/gh-pages_tmp/tutorial.Rmd">here</a> to follow along with the example below

First load the **SCRuB** packages into R:
```
library(SCRuB)
```

Then, load the datasets:
```
data <- read.csv('tutorial_data/hospital.csv', row.names=1) %>% as.matrix()
```

Next, load the well metadata
```
well_metadata <- read.csv('tutorial_data/well_metadata.csv', row.names=2)[row.names(data), ] %>% 
        select(na_plate_location) %>% 
        mutate(well = na_plate_location %>% sapply( function(x) which( LETTERS == substr(x, 1, 1) ) ),
                               indices = na_plate_location %>% sapply( function(x) substr(x, 2, 3) %>% as.integer)
                                 )
```

Run _SCRuB, saving the output with prefix "demo":

```
scr_out <- spatial_SCRUB(data = data, 
                            is_control=control_indices, 
                            well_dists = distance_matrix, 
                            print_loglikelihood = TRUE
                            )
```


Input format
-----------------------
The input to *SCRuB* is composed of one count matrix, one vector identify the matrix control samples, and a distance matrix illustrating the samples:

(1) data - ( n_samples + n_controls ) x n_taxa count matrix, where m is the number samples and n is the number of taxa. Row names are the sample ids ('SampleID'). Column names are the taxa ids. Every consecutive column contains read counts for each sample. Note that this order must be respected.


count matrix (first 4 rows and columns):

| | taxa_1 | taxa_2 | taxa_3 | taxa_4 |
| ------------- | ------------- |------------- |------------- |------------- |
| ERR525698  |  0 | 5 | 0|20 |
| ERR525693  |  15 | 5 | 0|0 |
| ERR525688  |  0 | 13 | 0| 200 |
| ERR525699  |  4 | 5 | 0|0 |

(2) is_control. A vector identifying which rows of _data_ correspond to control samples:
```
c(F, F, T, F)
```

(3) _well_dist_ - An ( n_samples + n_controls ) x ( n_samples + n_controls ) distance matrix, summarizing the pairwise distance between each sample. Both the row names and column names of _well_dist_ must correspond to the row names of the _data_ matrix.

| | ERR525698 |ERR525693 | ERR525688| ERR525699|
| ------------- | ------------- |------------- |------------- |------------- |
| ERR525698  |  0 | 1 | 1 | 1.4 |
| ERR525693  |  1 | 0 | 1.4| 1 |
| ERR525688  |  1 | 1.4 | 0| 1 |
| ERR525699  |  1.4 | 1 | 1| 0 |


 

Output format
-----------------------
SCRuB outputs a list containg the following:

(1) decontaminated_samples

| | taxa_1 | taxa_2 | taxa_3 | taxa_4 |
| ------------- | ------------- |------------- |------------- |------------- |
| ERR525698  |  0 | 4 | 0 | 0 |
| ERR525693  |  15 | 5 | 0 | 0 |
| ERR525699  |  4 | 5 | 0 | 0 |

(2) p - The estimated fraction of each sample that is not contamination
```
c(0.16, 1, 1)
```

(3) alpha - the estimated composition of each control:
| | ERR525698  | ERR525693 | ERR525699 | contaminant  |  
| ------------- | ------------- | ------------- | ------------- | ------------- |
| ERR525688 | 0.001 | 0.001 | 0.001 | 0.997 |



(4) gamma - The estimated relative abundance of the contamination community
```
c(0, 0.05, 0, 0.95)
```










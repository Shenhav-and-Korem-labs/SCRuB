# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(SCRuB)
Sys.setenv(TORCH_INSTALL=1)
# install_torch() Instead using environment variable following the error described here: https://github.com/mlverse/torch/issues/1017

test_check("SCRuB")

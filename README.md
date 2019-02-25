# LARS-EN
Calculate the regression coefficients for the Elastic-Net family

The algorithm will compute the regression coefficients for the whole path for the penalization parameters 
in the Elastic-Net regression considering a variance-covariance matrix **C=X'X**/n for the (centered) predictors and a covariance (between the predictors and the response) vector **cov=X'y**/n.

The code is a modification to the 'lars' function from the R-package LARS by Hastie and Efron.

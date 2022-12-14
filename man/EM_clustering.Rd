% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM_clustering.R
\name{EM_clustering}
\alias{EM_clustering}
\title{Implements the EM algorithm for the unlabeled dataset (unsupervised clustering with EM)}
\usage{
EM_clustering(data_list, initial_params, iter = 1000)
}
\arguments{
\item{data_list}{A list of matrices of dimensions n_i x d.}

\item{initial_params}{list of length 3 with the initial values of the parameters> The names of list elements should be mu, sigma, and pi.}

\item{iter}{Path to the input file}
}
\value{
output is a list of 3 dimensional arrays where each array (d x d x g) represents a set and the third dimension of the array represents the group in that set.

clust_out$class: the class labels

clust_out$z_hard: Z matrix with elements 0 or 1 where 1 indicates membership of that class.

clust_out$z_soft: the estimated Z matrix with elements between 0 and 1.

clust_out$parameters: the estimated parameters in a list with names mu, sigma, and pi.
}
\description{
Implements the EM algorithm for the unlabeled dataset (unsupervised clustering with EM)
}
\examples{
test_mat = matrix(rep(1, 20), 5, nrow = 5, ncol = 4)
input_list = list(test_mat, test_mat)
mu = c(1:4)
sigma = array(NA, c(4, 4, 2))
sigma1 = diag(c(1,2,3,4))
sigma2 = diag(c(1,2,3,4))
sigma[ , , 1] = sigma1
sigma[ , , 2] = sigma2
pi_vec = c(.3, .7)
initial_params_list = list("mu" = mu, "sigma" = sigma, "pi" = pi)
clust_out = EM_clustering(input_list, initial_params_list, iter = 1000)
}

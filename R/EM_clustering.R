
#' Implements the EM algorithm for the unlabeled dataset (unsupervised clustering with EM)
#'
#'
#' @param data_list A list of matrices of dimensions n_i x d.
#' @param initial_params  list of length 3 with the initial values of the parameters> The names of list elements should be mu, sigma, and pi.
#' @param iter Path to the input file
#'
#' @return output is a list of 3 dimensional arrays where each array (d x d x g) represents a set and the third dimension of the array represents the group in that set.
#' @return clust_out$class: the class labels
#' @return clust_out$z_hard: Z matrix with elements 0 or 1 where 1 indicates membership of that class.
#' @return clust_out$z_soft: the estimated Z matrix with elements between 0 and 1.
#' @return clust_out$parameters: the estimated parameters in a list with names mu, sigma, and pi.
#'
#' @examples test_mat = matrix(rep(1, 20), 5, nrow = 5, ncol = 4)
#' @examples input_list = list(test_mat, test_mat)
#' @examples mu = c(1:4)
#' @examples sigma = array(NA, c(4, 4, 2))
#' @examples sigma1 = diag(c(1,2,3,4))
#' @examples sigma2 = diag(c(1,2,3,4))
#' @examples sigma[ , , 1] = sigma1
#' @examples sigma[ , , 2] = sigma2
#' @examples pi_vec = c(.3, .7)
#' @examples initial_params_list = list("mu" = mu, "sigma" = sigma, "pi" = pi)
#' @examples clust_out = EM_clustering(input_list, initial_params_list, iter = 1000)
#' @export
EM_clustering = function(data_list, initial_params, iter = 1000) {

  if (!is.list(data_list) || length(data_list) > 1) {
    stop("data_list should be a list of length 1 where data_list[[1]] is the unlabelled dataset.")
  }

  if (!setequal(names(initial_params), c("mu", "sigma", "pi"))) {
    stop("Names of the parameters in the initial parameter list are incorrect. Names should be mu, sigma, and pi")
  }

  out = EM_fractional(data_list, Z_list = list(), w_list = list(1), initial_params, iter)
  return(out)
}

#' EM for fractional supervised learning
#'
#' @param data_list A list of dataframes or matrices of dimensions n_i x d.
#' @param Z_list A list of matrices with dimension n_i x g where matrix elements are 0 or 1.
#' @param w_list A list of weights for each set where length(list_data) = length(list_Z) = length(list_w)
#' @param initial_params A matrix of d x g.
#' @param iter Path to the input file
#'
#' @return output is a list of 3 dimensional arrays where each array (d x d x g) represents a set and the third dimension of the array represents the group in that set.
#' @return clust_out$class: the class labels
#' @return clust_out$z_hard: Z matrix with elements 0 or 1 where 1 indicates membership of that class.
#' @return clust_out$z_soft: the estimated Z matrix with elements between 0 and 1.
#' @return clust_out$parameters: the estimated parameters in a list with names mu, sigma, and pi.
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

EM_fractional = function(data_list, Z_list, w_list, initial_params, iter = 1000) {

  if (!is.list(data_list)) {
    stop("data_list should be a list.")
  }

  if (!is.list(Z_list)) {
    stop("Z_list should be a list.")
  }

  if (!is.list(w_list)) {
    stop("w_list should be a list.")
  }

  if (length(Z_list) != 0) {
    if ((length(data_list) != length(Z_list)) || (length(Z_list) != length(w_list))) {
      stop("data_list, Z_list, w_list should be of equal length.")
    }
  } else {
    if ((length(data_list) != length(w_list))) {
      stop("data_list and w_list  should be of equal length.")
    }
  }

  if (!setequal(names(initial_params), c("mu", "sigma", "pi"))) {
    stop("Names of the parameters in the initial parameter list are incorrect. Names should be mu, sigma, and pi")
  }

  if (!is.array(initial_params$sigma)) {
    stop('Sigma should be an array of dimension d x d x g')
  }

  if (!is.matrix(initial_params$mu)) {
    stop('mu should be a matrix of dimension d x g')
  }

  if (!is.numeric(initial_params$pi)) {
    stop('pi should be a vector of length g')
  }

  d = nrow(initial_params$sigma[ , , 1])
  g = length(initial_params$pi)

  mu = initial_params$mu
  s = initial_params$sigma
  pi = initial_params$pi
  unlabelled = length(w_list)
  Z_hat = Z_list

  for (k in 2:iter) {
    Z_unlabelled = E_fn(data_list, mu, s, pi)
    #Z_hat = list(Z_list[[1:(unlabelled-1)]], Z_unlabelled)
    Z_hat[[unlabelled]] = Z_unlabelled
    par_updates = M_fn(data_list, Z_hat, w_list)
    mu = par_updates$mu
    s = par_updates$sigma
    pi = par_updates$pi
  }
  all_param = list("mu" = mu, "sigma" = s, "pi" = pi)
  predicted_Z = map_fn(Z_unlabelled)
  class_labels = label_from_Z(predicted_Z)
  out = list("class" = class_labels, "z_hard" = predicted_Z, "z_soft" = Z_unlabelled, 'parameters' = all_param)
  return(out)
}

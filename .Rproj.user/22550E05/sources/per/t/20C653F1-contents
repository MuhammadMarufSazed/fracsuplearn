E_fn<-function(data_list, mu, sigma, pi_vec){
  unlabelled_set = length(data_list) # assuming that the unlabelled set will the last in the list data_list
  Z_hat = matrix(0, nrow = nrow(data_list[[unlabelled_set]]), ncol = length(pi_vec)) # n_{ul} by g matrix where n_{ul} is the number of unlabelled observations

  for (g in 1:length(pi_vec)) {
    Z_hat[ , g] = mvtnorm::dmvnorm(data_list[[unlabelled_set]], mu[ , g], sigma[ , , g])*pi_vec[g]
  }

  Z_hat_sum = apply(Z_hat, 1, sum)
  Z_hat_stand = Z_hat/Z_hat_sum

  return(Z_hat_stand)
}


pi_fn = function(Z, w) {
  Z_new = Z*w
  sum_z = apply(Z_new, 2, sum)
  return(sum_z)
}

pi_fn_all = function(list_Z, list_w) {
  out_all = list()

  for (i in 1:length(list_Z)) {
    out = pi_fn(list_Z[[i]], list_w[[i]])
    out_all[[i]] = out
  }

  return(out_all)
}


mu_fn = function(data_input, Z, w) {
  d = dim(data_input)[2]
  G = dim(Z)[2]
  Z_new = Z*w # n by g
  out_all_mu = matrix(0, d, G) # d by g matrix where d is the dimension and g is number of groups

  for (g in 1:G) {
    val = data_input * Z_new[ , g]
    out_all_mu[ , g] = apply(val, 2, sum)
  }

  return(out_all_mu)
}

mu_fn_all = function(list_data, list_Z, list_w) {
    out_all = list()

  for (i in 1:length(list_Z)) {
    out = mu_fn(list_data[[i]], list_Z[[i]], list_w[[i]])
    out_all[[i]] = out
  }

  return(out_all)
}


sigma_fn = function(data_input, Z, w, mu) {
    d = ncol(data_input)
  G = ncol(Z)
  n = nrow(data_input)
  Z_new = Z*w # n by g

  out1_test1 = array(NA, dim=c(d, d, n))
  out1_test2 = array(NA, dim=c(d, d, n))
  sigma_star = array(NA, dim = c(d, d, G)) #p by p by g array for sigma

  for (g in 1:G) {
      for (i in 1:n){
      out1_test1[ , , i] = (data_input[i, ] - mu[ , g]) %*% t(data_input[i, ] - mu[ , g])
      out1_test2[ , , i] = out1_test1[ , , i]  * Z_new[i , g]
    }

    sigma_star[ , , g] = rowSums(out1_test2, dims = 2)
  }
  return(sigma_star)
}

sigma_fn_all = function(list_data, list_Z, list_w, mu){

  out_all = list()

  for (i in 1:length(list_Z)) {
    out = sigma_fn(list_data[[i]], list_Z[[i]], list_w[[i]], mu)
    out_all[[i]] = out
  }
  return(out_all)
}


M_fn = function(data_list, Z_hat, w_list) {
  pi_all = pi_fn_all(Z_hat, w_list)
  combined_pi = Reduce("+", pi_all)
  pi_update = combined_pi/sum(combined_pi)

  mu_all = mu_fn_all(data_list, Z_hat, w_list)
  combined_mu = Reduce("+", mu_all)
  mu_update = t(t(combined_mu)/combined_pi)

  # combining across all sets of data for sigma
  sigma_all = sigma_fn_all(data_list, Z_hat, w_list, mu_update) # sigma is updated last because it requires the update of mu
  set_no = length(data_list)
  d = ncol(data_list[[1]])
  G = dim(sigma_all[[1]])[3]
  final_sigma = array(NA, c(d, d, length(pi_update)))
  for (g in 1:G) {
    all_set_sigma = array(NA, c(d, d, set_no))
    for(set in 1:set_no) {
      all_set_sigma[ , , set] = sigma_all[[set]][ , , g]
    }
    final_sigma[ , , g] = rowSums(all_set_sigma, dims = 2)/combined_pi[g]
  }
  out = list("mu" = mu_update, "sigma" = final_sigma, "pi" = pi_update)
  return(out)
}









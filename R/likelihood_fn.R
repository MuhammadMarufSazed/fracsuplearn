liklihood_oneset = function(data, mu, sigma) {
  density_val = matrix(NA, nrow = nrow(data), ncol = ncol(mu))
  G = ncol(mu)
  for (g in 1:G) {
    density_val[ , g] = dmvnorm(data, mu[ , g], sigma[ , , g])
  }
  return(density_val)
}

liklihood_allset = function(data_list, mu, sigma) {
  set_no = length(data_list)
  all_density = list()
  for (s in 1:set_no) {
    all_density[[s]] = liklihood_oneset(data_list[[s]], mu, sigma)
  }
  return(all_density)
}


weighted_lik = function(data_list, mu, sigma, Z_list, w_list) {
  all_w = unlist(w_list)
  set_no = length(data_list)
  all_density = liklihood_allset(data_list, mu, sigma)
  val = numeric(set_no)
  for (s in 1:set_no) {
    out = all_densit[[s]] * Z_list[[s]] 
    sum_out = apply(out, 1, sum)
    val[s] = sum_out
  }
  all_L = val ^ all_w
  final_val = prod(all_L)
  return(final_val)
}

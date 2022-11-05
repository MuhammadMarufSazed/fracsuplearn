#' Creates Z matrix from given class labels
#'
#' @param input_vec the vector that contains the class labels
#' @param totaL_groups total number of groups; an integer
#' @export
z_from_label = function(input_vec, totaL_groups) {
  n = length(input_vec)
  g = (totaL_groups)
  Z_labelled = matrix(0, nrow = n, ncol = g) # n by g matrix
  for (i in 1:n) {
    Z_labelled[i, input_vec[i]] = 1
  }
  return(Z_labelled)
}

#' Creates Create a hard Z matrix from a soft Z matrix
#'
#' @param Zhat the Z matrix that is extimated and values between 0 and 1
#' @export
map_fn<-function(Zhat) {temp = apply(Zhat, 1, function(x) {
  map_est = numeric(length(x))
  ind = which.max(x)
  map_est[ind] = 1
  map_est[-ind] = 0
  return(map_est)
})
return(t(temp))
}

#' Creates Create labels from a hard Z matrix
#'
#' @param Z the hard Z matrix with elements only 0 or 1
#' @export
label_from_Z = function(Z) {

  class_labels = numeric(nrow(Z))
  for (i in 1:nrow(Z)) {
    class_labels[i] = which(Z[i, ] == 1)
  }

  return(class_labels)
}

within_grp_W = function(input_data_class_c) {
  d = ncol(input_data_class_c)
  n = nrow(input_data_class_c)
  dataX = input_data_class_c[ , d-1]
  xbar = apply(dataX, 2, mean)
  val = dataX - dataX
  all_mat = array(NA, c(d, d, n))
  for(i in 1:n) {
    all_mat[ , , i] = val[i, ] %*% t(val[i, ])
  }
  temp = rowSums(all_mat, dims = 2)
  return(temp)
}

#' Calculates the trace and the determinant of a within group matrix
#'
#' @param input_data_all_class the dataframe with the class labels
#' @export
criteria_W = function(input_data_all_class) {
  g = length(unique(input_data_all_class$class))
  class_list = unique(input_data_all_class$class)
  all_grp = array(NA, c(d, d, g))
  for (c in class_list) {
    data_class_c = input_data_all_class[input_data_all_class$class == c, ]
    all_grp[ , , c] = within_grp_W(data_class_c)
  }
  temp = rowSums(all_grp, dims = 2)
  det_val = det(temp)
  trace_val = sum(diag(temp))
  out = list("determinant" = det_val, "trace" = trace_val)
  return(out)
}

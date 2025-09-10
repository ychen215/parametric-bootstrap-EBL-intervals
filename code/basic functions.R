#--------- Estimation of A ---------------

#---------Prasad-Rao Method -----------
PR_fun = function(x, y, d, m, p){
  beta_ols = solve(t(x) %*% x) %*% t(x) %*% y
  h = diag(x %*% solve(t(x) %*% x) %*% t(x))
  A_pr = (t(y - x %*% beta_ols) %*% (y - x %*% beta_ols) - sum(d * (1 - h))) / (m - p)
  return(A_pr)
}

#---------Fay-Herriot Method -----------
F_fh = function(Ahat, y, d, x){
  Vi = diag(1/(c(Ahat) + d))
  yhat = x %*% solve(t(x) %*% Vi %*% x) %*% t(x) %*% Vi %*% y
  f_fh = sum((y - yhat)^2 / (Ahat + d))
  return(f_fh)
}
# Expectation of the first derivative of f
G_fh = function(Ahat, d, x){
  Vi = diag(1/(c(Ahat) + d))
  P = Vi - Vi %*% x %*% solve(t(x) %*% Vi %*% x) %*% t(x) %*% Vi
  return(-sum(diag(P)))
}

# A estimator with PR or FH method
estimate_A = function(method, x, y, d, m, p, tol = 1e-4, maxiter = 100, k0 = 0) {
  if (method == "PR") {
    A <- PR_fun(x, y, d, m, p)
    if (A < 0) {
      A = 0.01
      k0 = 1
    }
    return(list(A = A, k0 = k0))
  }
  
  if (method == "FH") {
    if ((F_fh(0, y, d, x) - (m - p)) < 0) return(list(A = 0.01, k0 = 1))
    A = 0
    k = 0
    diff = 10
    while (diff > tol && k < maxiter) {
      Anew = A + (m - p - F_fh(A, y, d, x)) / G_fh(A, d, x)
      Anew = ifelse(Anew < 0, -Anew / 10, Anew)
      diff = abs(Anew - A)
      A = Anew
      k = k + 1
    }
    if (A < 0 || k == maxiter) {
      A = 0.01
      k0 = 1
      }
    return(list(A = A, k0 = k0))
  }
}

# EBLUP and mse with PR or FH variance estimator
compute_eblup = function(A, x, y, d){
  Vi <- diag(1 / (c(A) + d))
  beta <- solve(t(x) %*% Vi %*% x) %*% t(x) %*% Vi %*% y
  B <- d / (d + c(A))
  eblup <- (1 - B) * y + B * (x %*% beta)
  list(A = A, beta = beta, B = B, eblup = eblup, Vi = Vi)
}

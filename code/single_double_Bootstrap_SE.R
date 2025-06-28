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

# library(flexsurv)
library(mvtnorm)
#-----Generate Data Parameters
m = 15  # number of small areas
# m = 50
Nosim = 1000
B = 400
DB = 100
mm = m/5
d = rep(c(4.0, 0.6, 0.5, 0.4, 0.2), rep(mm, 5))
x = rep(1, m)
beta = 0
A = 1.0
Alpha = c(0.2, 0.1, 0.05)
p = 1
df = 9
tol = 0.0001
maxiter = 100

## number of negative estimates
cc_fh = 0
cc_pr = 0
## number of negative estimates in single bootstrap
s_fh = 0
s_pr = 0
## number of negative estimates in double bootstrap
ss_fh = 0
ss_pr = 0

Int1 <- matrix(0, m, 9)
Int2 <- matrix(0, m, 9)
Int3 <- matrix(0, m, 9)
qdif1 <- matrix(0, m, 9)
qdif2 <- matrix(0, m, 9)
qdif3 <- matrix(0, m, 9)

system.time(
  for (s in 1:Nosim) {
    v = (rexp(m)-1)*sqrt(A)
    e = rmvnorm(1, sigma = diag(d))
    theta = x * beta + v 
    y = theta + matrix(e)
    
    Direct = y
    
    #------ theta_PR-----
    A_PR = PR_fun(x, y, d, m, p)
    if(A_PR < 0){
      A_PR = 0.01
      cc_pr = cc_pr + 1
    }
    
    Vi_PR = diag(1/(c(A_PR) + d))
    beta_PR = solve(t(x) %*% Vi_PR %*% x) %*% t(x) %*% Vi_PR %*% y
    B_PR = d / (d + c(A_PR))
    eblup.PR = (1 - B_PR) * y + B_PR * (x %*% beta_PR)
    
    g1_PR = d * (1 - B_PR)
    g2_PR = B_PR^2 * diag(x %*% solve(t(x) %*% Vi_PR %*% x) %*% t(x))
    g3_PR = d^2 * (c(A_PR) + d)^(-3) * 2 * sum((c(A_PR) + d)^2) / m^2
    
    mse_PR = g1_PR + g2_PR + 2*g3_PR

    #------ theta_FH-----
    if((F_fh(0, y, d, x)-(m-p))<0){
      A_FH = 0.01
      cc_fh = cc_fh +1
    }
    else{
      A_FH = 0
      k=0
      diff = 10
      while((diff > tol)&(k < maxiter)) {
        Anew_fh = A_FH + (m - p - F_fh(A_FH, y, d, x)) / G_fh(A_FH, d, x)
        Anew_fh = ifelse(Anew_fh < 0, -Anew_fh/10, Anew_fh)
        diff = abs(Anew_fh - A_FH)
        A_FH = Anew_fh
        k = k +1
      }
      if((A_FH < 0)|(k==maxiter)) {
        A_FH = 0.01
        cc_fh = cc_fh + 1
      }
    }
    Vi_FH = diag(1/(c(A_FH) + d))
    beta_FH = solve(t(x) %*% Vi_FH %*% x) %*% t(x) %*% Vi_FH %*% y
    B_FH = d / (d + c(A_FH))
    eblup.FH = (1 - B_FH) * y + B_FH * (x %*% beta_FH)
    
    g1_FH = d * (1 - B_FH)
    g2_FH = B_FH^2 * diag(x %*% solve(t(x) %*% Vi_FH %*% x) %*% t(x))
    g3_FH = 2 * m * d^2 * (c(A_FH) + d)^(-3) * (sum(Vi_FH))^(-2)
    g4_FH = 2 * B_FH^2 * (m * sum(diag(Vi_FH^2)) - (sum(diag(Vi_FH)))^2) / (sum(diag(Vi_FH)))^3
    
    mse_FH = g1_FH + g2_FH + 2*g3_FH - g4_FH
    mse_FH2 = g1_FH + g2_FH + 2*g3_FH
    if(sum(mse_FH<0)>0){mse_FH[mse_FH<0] = mse_FH2[mse_FH<0]}
    
    Z.j.PR = matrix(0, B, m)
    Z.j.FH = matrix(0, B, m)
    ####------Parametric Bootstrap------
    pivot_PR = matrix(0, B, m)
    HM_PR = matrix(0, B, m)
    pivot_FH = matrix(0, B, m)
    HM_FH = matrix(0, B, m)
    for (bb in 1:B) {
      bv_PR = (rexp(m)-1)*sqrt(A_PR)
      btheta_PR = c(beta_PR) + bv_PR
      
      bv_FH = (rexp(m)-1)*sqrt(A_FH)
      btheta_FH = c(beta_FH) + bv_FH
    
      be = rmvnorm(1, sigma = diag(d))
      
      by_FH = btheta_FH + matrix(be)
      by_PR = btheta_PR + matrix(be)
      
      #Single bootstrap for PR
      A.star.PR = PR_fun(x, by_PR, d, m, p)
      if(A.star.PR < 0){
        A.star.PR = 0.01
        s_pr = s_pr + 1
      }
      
      Vi.star.PR = diag(1/(c(A.star.PR) + d))
      beta.s.PR = solve(t(x) %*% Vi.star.PR %*% x) %*% t(x) %*% Vi.star.PR %*% by_PR
      B.star.PR = d / (d + c(A.star.PR))
      eblup.s.PR = (1 - B.star.PR) * by_PR + B.star.PR * (x %*% beta.s.PR)
      g1.s.PR = d * (1 - B.star.PR)
      
      pivot_PR[bb, ] = t((btheta_PR - eblup.s.PR) / sqrt(g1.s.PR))
      HM_PR[bb, ] = t((btheta_PR - x %*% beta.s.PR) / c(sqrt(A.star.PR)))
      
      #Single bootstrap for FH
      if((F_fh(0, by_FH, d, x)-(m-p))<0){
        A.star.FH = 0.01
        s_fh = s_fh +1
      }
      else{
        A.star.FH = 0
        diff = 10
        k = 0
        while((diff > tol) & (k < maxiter)) {
          Anew_fh = A.star.FH + (m - p - F_fh(A.star.FH, by_FH, d, x)) / G_fh(A.star.FH, d, x)
          Anew_fh = ifelse(Anew_fh < 0, -Anew_fh/10, Anew_fh)
          diff = abs(Anew_fh - A.star.FH)
          A.star.FH = Anew_fh
          k = k + 1
        }
        if((A.star.FH < 0)|(k==maxiter)) {
          A.star.FH = 0.01
          s_fh = s_fh +1
        }
      }
      Vi.star.FH = diag(1/(c(A.star.FH) + d))
      beta.s.FH = solve(t(x) %*% Vi.star.FH %*% x) %*% t(x) %*% Vi.star.FH %*% by_FH
      B.star.FH = d / (d + c(A.star.FH))
      eblup.s.FH = (1 - B.star.FH) * by_FH + B.star.FH * (x %*% beta.s.FH)
      g1.s.FH = d * (1 - B.star.FH)
      
      pivot_FH[bb, ] = t((btheta_FH - eblup.s.FH) / sqrt(g1.s.FH))
      
      HM_FH[bb, ] = t((btheta_FH - x %*% beta.s.FH) / c(sqrt(A.star.FH)))
      
      d_pivot_PR = matrix(0, DB, m)
      d_HM_PR = matrix(0, DB, m)
      d_pivot_FH = matrix(0, DB, m)
      d_HM_FH = matrix(0, DB, m)
      for (dd in 1:DB) {

        dv_FH = (rexp(m)-1)*sqrt(A.star.FH)
        dtheta_FH = c(beta.s.FH) + dv_FH

        dv_PR = (rexp(m)-1)*sqrt(A.star.PR)
        dtheta_PR = c(beta.s.PR) + dv_PR
        
        de = rmvnorm(1, sigma = diag(d))
        
        dy_FH = dtheta_FH + matrix(de)
        dy_PR = dtheta_PR + matrix(de)
        
        # Double bootstrap for PR
        A.ss.PR = PR_fun(x, dy_PR, d, m, p)
        if(A.ss.PR < 0){
          A.ss.PR = 0.01
          ss_pr = ss_pr + 1
        }
        
        Vi.ss.PR = diag(1/(c(A.ss.PR) + d))
        beta.ss.PR = solve(t(x) %*% Vi.ss.PR %*% x) %*% t(x) %*% Vi.ss.PR %*% dy_PR
        B.ss.PR = d / (d + c(A.ss.PR))
        eblup.ss.PR = (1 - B.ss.PR) * dy_PR + B.ss.PR * (x %*% beta.ss.PR)
        g1.ss.PR = d * (1 - B.ss.PR)
        
        d_pivot_PR[dd, ] = t((dtheta_PR - eblup.ss.PR) / sqrt(g1.ss.PR))
        d_HM_PR[dd, ] = t((dtheta_PR - x %*% beta.ss.PR) / c(sqrt(A.ss.PR)))
        
        # Double bootstrap for FH
        if((F_fh(0, dy_FH, d, x)-(m-p))<0){
          A.ss.FH = 0.01
          ss_fh = ss_fh +1
        }
        else{
          A.ss.FH = 0
          diff = 10
          k = 0
          while((diff > tol) & (k < maxiter)) {
            Anew_fh = A.ss.FH + (m - p - F_fh(A.ss.FH, dy_FH, d, x)) / G_fh(A.ss.FH, d, x)
            Anew_fh = ifelse(Anew_fh < 0, -Anew_fh/10, Anew_fh)
            diff = abs(Anew_fh - A.ss.FH)
            A.ss.FH = Anew_fh
            k = k+1
          }
          if((A.ss.FH < 0)|(k==maxiter)) {
            A.ss.FH = 0.01
            ss_fh = ss_fh +1
          }
        }
        
        Vi.ss.FH = diag(1/(c(A.ss.FH) + d))
        beta.ss.FH = solve(t(x) %*% Vi.ss.FH %*% x) %*% t(x) %*% Vi.ss.FH %*% dy_FH
        B.ss.FH = d / (d + c(A.ss.FH))
        eblup.ss.FH = (1 - B.ss.FH) * dy_FH + B.ss.FH * (x %*% beta.ss.FH)
        g1.ss.FH = d * (1 - B.ss.FH)
        
        d_pivot_FH[dd, ] = t((dtheta_FH - eblup.ss.FH) / sqrt(g1.ss.FH))
        d_HM_FH[dd, ] = t((dtheta_FH - x %*% beta.ss.FH) / c(sqrt(A.ss.FH)))
      }
      
      ###---Double bootstrap -------
      eblup.PR.mat = matrix(rep(pivot_PR[bb, ], rep(DB, m)), DB, m, byrow = F)
      eblup.FH.mat = matrix(rep(pivot_FH[bb, ], rep(DB, m)), DB, m, byrow = F)
      Z.j.PR[bb, ] = colSums(d_pivot_PR<=eblup.PR.mat)/DB
      Z.j.FH[bb, ] = colSums(d_pivot_FH<=eblup.FH.mat)/DB
      
      plot(s,ylim=c(0,Nosim),main=paste(s,"/",bb))
    
    }
    ################################################################################
    #####----Single Bootstrap------
    q.s.upper1 <- matrix(0, m, 4)
    q.s.lower1 <- matrix(0, m, 4)
    q.d.upper1 <- matrix(0, m, 2)
    q.d.lower1 <- matrix(0, m, 2)
    
    q.s.upper2 <- matrix(0, m, 4)
    q.s.lower2 <- matrix(0, m, 4)
    q.d.upper2 <- matrix(0, m, 2)
    q.d.lower2 <- matrix(0, m, 2)
    
    q.s.upper3 <- matrix(0, m, 4)
    q.s.lower3 <- matrix(0, m, 4)
    q.d.upper3 <- matrix(0, m, 2)
    q.d.lower3 <- matrix(0, m, 2)
    for(i in 1:m){
      q.s.lower1[i,]<-c(quantile(pivot_FH[, i], prob=Alpha[1]/2), quantile(HM_FH[, i], prob=Alpha[1]/2),
                       quantile(pivot_PR[, i], prob=Alpha[1]/2), quantile(HM_PR[, i], prob=Alpha[1]/2))
      q.s.upper1[i,]<-c(quantile(pivot_FH[, i], prob=1-Alpha[1]/2), quantile(HM_FH[, i], prob=1-Alpha[1]/2),
                       quantile(pivot_PR[, i], prob=1-Alpha[1]/2), quantile(HM_PR[, i], prob=1-Alpha[1]/2))
      q.s.lower2[i,]<-c(quantile(pivot_FH[, i], prob=Alpha[2]/2), quantile(HM_FH[, i], prob=Alpha[2]/2),
                        quantile(pivot_PR[, i], prob=Alpha[2]/2), quantile(HM_PR[, i], prob=Alpha[2]/2))
      q.s.upper2[i,]<-c(quantile(pivot_FH[, i], prob=1-Alpha[2]/2), quantile(HM_FH[, i], prob=1-Alpha[2]/2),
                        quantile(pivot_PR[, i], prob=1-Alpha[2]/2), quantile(HM_PR[, i], prob=1-Alpha[2]/2))
      q.s.lower3[i,]<-c(quantile(pivot_FH[, i], prob=Alpha[3]/2), quantile(HM_FH[, i], prob=Alpha[3]/2),
                        quantile(pivot_PR[, i], prob=Alpha[3]/2), quantile(HM_PR[, i], prob=Alpha[3]/2))
      q.s.upper3[i,]<-c(quantile(pivot_FH[, i], prob=1-Alpha[3]/2), quantile(HM_FH[, i], prob=1-Alpha[3]/2),
                        quantile(pivot_PR[, i], prob=1-Alpha[3]/2), quantile(HM_PR[, i], prob=1-Alpha[3]/2))
    }
    ###---Double bootstrap -------
    #####----CLL-----
    au.CLL.FH1 = apply(Z.j.FH, 2, quantile, prob=1-Alpha[1]/2)
    al.CLL.FH1 = apply(Z.j.FH, 2, quantile, prob=Alpha[1]/2)
    au.CLL.PR1 = apply(Z.j.PR, 2, quantile, prob=1-Alpha[1]/2)
    al.CLL.PR1 = apply(Z.j.PR, 2, quantile, prob=Alpha[1]/2)
    
    au.CLL.FH2 = apply(Z.j.FH, 2, quantile, prob=1-Alpha[2]/2)
    al.CLL.FH2 = apply(Z.j.FH, 2, quantile, prob=Alpha[2]/2)
    au.CLL.PR2 = apply(Z.j.PR, 2, quantile, prob=1-Alpha[2]/2)
    al.CLL.PR2 = apply(Z.j.PR, 2, quantile, prob=Alpha[2]/2)
    
    au.CLL.FH3 = apply(Z.j.FH, 2, quantile, prob=1-Alpha[3]/2)
    al.CLL.FH3 = apply(Z.j.FH, 2, quantile, prob=Alpha[3]/2)
    au.CLL.PR3 = apply(Z.j.PR, 2, quantile, prob=1-Alpha[3]/2)
    al.CLL.PR3 = apply(Z.j.PR, 2, quantile, prob=Alpha[3]/2)
    
    for(i in 1:m){
      q.d.lower1[i,]<-c(quantile(pivot_FH[, i], prob=al.CLL.FH1[i]), quantile(pivot_PR[, i], prob=al.CLL.PR1[i]))
      q.d.upper1[i,]<-c(quantile(pivot_FH[, i], prob=au.CLL.FH1[i]), quantile(pivot_PR[, i], prob=au.CLL.PR1[i]))
      
      q.d.lower2[i,]<-c(quantile(pivot_FH[, i], prob=al.CLL.FH2[i]), quantile(pivot_PR[, i], prob=al.CLL.PR2[i]))
      q.d.upper2[i,]<-c(quantile(pivot_FH[, i], prob=au.CLL.FH2[i]), quantile(pivot_PR[, i], prob=au.CLL.PR2[i]))
      
      q.d.lower3[i,]<-c(quantile(pivot_FH[, i], prob=al.CLL.FH3[i]), quantile(pivot_PR[, i], prob=al.CLL.PR3[i]))
      q.d.upper3[i,]<-c(quantile(pivot_FH[, i], prob=au.CLL.FH3[i]), quantile(pivot_PR[, i], prob=au.CLL.PR3[i]))
    }
    
    qdif1[, 1] <- qdif1[, 1] + (q.s.upper1[, 1] - q.s.lower1[, 1]) * sqrt(g1_FH)
    qdif1[, 2] <- qdif1[, 2] + (q.s.upper1[, 2] - q.s.lower1[, 2]) * sqrt(A_FH)
    qdif1[, 3] <- qdif1[, 3] + (q.d.upper1[, 1] - q.d.lower1[, 1]) * sqrt(g1_FH)
    
    qdif1[, 4] <- qdif1[, 4] + (q.s.upper1[, 3] - q.s.lower1[, 3]) * sqrt(g1_PR)
    qdif1[, 5] <- qdif1[, 5] + (q.s.upper1[, 4] - q.s.lower1[, 4]) * sqrt(A_PR)
    qdif1[, 6] <- qdif1[, 6] + (q.d.upper1[, 2] - q.d.lower1[, 2]) * sqrt(g1_PR)
    
    qdif1[, 7] <- qdif1[, 7] + 2 * qnorm(1-Alpha[1]/2) * sqrt(mse_FH)
    qdif1[, 8] <- qdif1[, 8] + 2 * qnorm(1-Alpha[1]/2) * sqrt(mse_PR)
    qdif1[, 9] <- qdif1[, 9] + 2 * qnorm(1-Alpha[1]/2) * sqrt(d)
    
    
    
    qdif2[, 1] <- qdif2[, 1] + (q.s.upper2[, 1] - q.s.lower2[, 1]) * sqrt(g1_FH)
    qdif2[, 2] <- qdif2[, 2] + (q.s.upper2[, 2] - q.s.lower2[, 2]) * sqrt(A_FH)
    qdif2[, 3] <- qdif2[, 3] + (q.d.upper2[, 1] - q.d.lower2[, 1]) * sqrt(g1_FH)
    
    qdif2[, 4] <- qdif2[, 4] + (q.s.upper2[, 3] - q.s.lower2[, 3]) * sqrt(g1_PR)
    qdif2[, 5] <- qdif2[, 5] + (q.s.upper2[, 4] - q.s.lower2[, 4]) * sqrt(A_PR)
    qdif2[, 6] <- qdif2[, 6] + (q.d.upper2[, 2] - q.d.lower2[, 2]) * sqrt(g1_PR)
    
    qdif2[, 7] <- qdif2[, 7] + 2 * qnorm(1-Alpha[2]/2) * sqrt(mse_FH)
    qdif2[, 8] <- qdif2[, 8] + 2 * qnorm(1-Alpha[2]/2) * sqrt(mse_PR)
    qdif2[, 9] <- qdif2[, 9] + 2 * qnorm(1-Alpha[2]/2) * sqrt(d)
    
    
    
    
    qdif3[, 1] <- qdif3[, 1] + (q.s.upper3[, 1] - q.s.lower3[, 1]) * sqrt(g1_FH)
    qdif3[, 2] <- qdif3[, 2] + (q.s.upper3[, 2] - q.s.lower3[, 2]) * sqrt(A_FH)
    qdif3[, 3] <- qdif3[, 3] + (q.d.upper3[, 1] - q.d.lower3[, 1]) * sqrt(g1_FH)
    
    qdif3[, 4] <- qdif3[, 4] + (q.s.upper3[, 3] - q.s.lower3[, 3]) * sqrt(g1_PR)
    qdif3[, 5] <- qdif3[, 5] + (q.s.upper3[, 4] - q.s.lower3[, 4]) * sqrt(A_PR)
    qdif3[, 6] <- qdif3[, 6] + (q.d.upper3[, 2] - q.d.lower3[, 2]) * sqrt(g1_PR)
    
    qdif3[, 7] <- qdif3[, 7] + 2 * qnorm(1-Alpha[3]/2) * sqrt(mse_FH)
    qdif3[, 8] <- qdif3[, 8] + 2 * qnorm(1-Alpha[3]/2) * sqrt(mse_PR)
    qdif3[, 9] <- qdif3[, 9] + 2 * qnorm(1-Alpha[3]/2) * sqrt(d)
    
    cov.CLL.FH1 <- ((theta <= (eblup.FH + q.s.upper1[,1] * sqrt(g1_FH))) & (theta >= (eblup.FH + q.s.lower1[,1] * sqrt(g1_FH))))
    cov.HM.FH1 <- ((theta <= (c(beta_FH) + q.s.upper1[,2] * sqrt(A_FH))) & (theta >= (c(beta_FH) + q.s.lower1[,2] * sqrt(A_FH))))
    cov.CLL.PR1 <- ((theta <= (eblup.PR + q.s.upper1[,3] * sqrt(g1_PR))) & (theta >= (eblup.PR + q.s.lower1[,3] * sqrt(g1_PR))))
    cov.HM.PR1 <- ((theta <= (c(beta_PR) + q.s.upper1[,4] * sqrt(A_PR))) & (theta >= (c(beta_PR) + q.s.lower1[,4] * sqrt(A_PR))))
    cov.FH1 <- ((theta <= (c(eblup.FH) + qnorm(1-Alpha[1]/2) * sqrt(mse_FH))) & (theta >= (c(eblup.FH) + qnorm(Alpha[1]/2) * sqrt(mse_FH))))
    cov.PR1 <- ((theta <= (c(eblup.PR) + qnorm(1-Alpha[1]/2) * sqrt(mse_PR))) & (theta >= (c(eblup.PR) + qnorm(Alpha[1]/2) * sqrt(mse_PR))))
    cov.DIR1 <- ((theta <= (c(Direct) + qnorm(1-Alpha[1]/2) * sqrt(d))) & (theta >= (c(Direct) + qnorm(Alpha[1]/2) * sqrt(d))))
    cov.DB.FH1 <- ((theta <= (eblup.FH + q.d.upper1[,1] * sqrt(g1_FH))) & (theta >= (eblup.FH + q.d.lower1[,1] * sqrt(g1_FH))))
    cov.DB.PR1 <- ((theta <= (eblup.PR + q.d.upper1[,2] * sqrt(g1_PR))) & (theta >= (eblup.PR + q.d.lower1[,2] * sqrt(g1_PR))))
    
    coverge1 <- matrix(0, m, 9)
    
    coverge1[, 1] <- ifelse(cov.CLL.FH1, 1, 0)
    coverge1[, 2] <- ifelse(cov.HM.FH1, 1, 0)
    coverge1[, 3] <- ifelse(cov.DB.FH1, 1, 0)
    
    coverge1[, 4] <- ifelse(cov.CLL.PR1, 1, 0)
    coverge1[, 5] <- ifelse(cov.HM.PR1, 1, 0)
    coverge1[, 6] <- ifelse(cov.DB.PR1, 1, 0)
    
    coverge1[, 7] <- ifelse(cov.FH1, 1, 0)
    coverge1[, 8] <- ifelse(cov.PR1, 1, 0)
    coverge1[, 9] <- ifelse(cov.DIR1, 1, 0)

    Int1 <- coverge1 + Int1
    
    
    
    cov.CLL.FH2 <- ((theta <= (eblup.FH + q.s.upper2[,1] * sqrt(g1_FH))) & (theta >= (eblup.FH + q.s.lower2[,1] * sqrt(g1_FH))))
    cov.HM.FH2 <- ((theta <= (c(beta_FH) + q.s.upper2[,2] * sqrt(A_FH))) & (theta >= (c(beta_FH) + q.s.lower2[,2] * sqrt(A_FH))))
    cov.CLL.PR2 <- ((theta <= (eblup.PR + q.s.upper2[,3] * sqrt(g1_PR))) & (theta >= (eblup.PR + q.s.lower2[,3] * sqrt(g1_PR))))
    cov.HM.PR2 <- ((theta <= (c(beta_PR) + q.s.upper2[,4] * sqrt(A_PR))) & (theta >= (c(beta_PR) + q.s.lower2[,4] * sqrt(A_PR))))
    cov.FH2 <- ((theta <= (c(eblup.FH) + qnorm(1-Alpha[2]/2) * sqrt(mse_FH))) & (theta >= (c(eblup.FH) + qnorm(Alpha[2]/2) * sqrt(mse_FH))))
    cov.PR2 <- ((theta <= (c(eblup.PR) + qnorm(1-Alpha[2]/2) * sqrt(mse_PR))) & (theta >= (c(eblup.PR) + qnorm(Alpha[2]/2) * sqrt(mse_PR))))
    cov.DIR2 <- ((theta <= (c(Direct) + qnorm(1-Alpha[2]/2) * sqrt(d))) & (theta >= (c(Direct) + qnorm(Alpha[2]/2) * sqrt(d))))
    
    cov.DB.FH2 <- ((theta <= (eblup.FH + q.d.upper2[,1] * sqrt(g1_FH))) & (theta >= (eblup.FH + q.d.lower2[,1] * sqrt(g1_FH))))
    cov.DB.PR2 <- ((theta <= (eblup.PR + q.d.upper2[,2] * sqrt(g1_PR))) & (theta >= (eblup.PR + q.d.lower2[,2] * sqrt(g1_PR))))
    
    coverge2 <- matrix(0, m, 9)
    
    coverge2[, 1] <- ifelse(cov.CLL.FH2, 1, 0)
    coverge2[, 2] <- ifelse(cov.HM.FH2, 1, 0)
    coverge2[, 3] <- ifelse(cov.DB.FH2, 1, 0)
    
    coverge2[, 4] <- ifelse(cov.CLL.PR2, 1, 0)
    coverge2[, 5] <- ifelse(cov.HM.PR2, 1, 0)
    coverge2[, 6] <- ifelse(cov.DB.PR2, 1, 0)
    
    coverge2[, 7] <- ifelse(cov.FH2, 1, 0)
    coverge2[, 8] <- ifelse(cov.PR2, 1, 0)
    coverge2[, 9] <- ifelse(cov.DIR2, 1, 0)
    
    Int2 <- coverge2 + Int2
  
    
    cov.CLL.FH3 <- ((theta <= (eblup.FH + q.s.upper3[,1] * sqrt(g1_FH))) & (theta >= (eblup.FH + q.s.lower3[,1] * sqrt(g1_FH))))
    cov.HM.FH3 <- ((theta <= (c(beta_FH) + q.s.upper3[,2] * sqrt(A_FH))) & (theta >= (c(beta_FH) + q.s.lower3[,2] * sqrt(A_FH))))
    cov.CLL.PR3 <- ((theta <= (eblup.PR + q.s.upper3[,3] * sqrt(g1_PR))) & (theta >= (eblup.PR + q.s.lower3[,3] * sqrt(g1_PR))))
    cov.HM.PR3 <- ((theta <= (c(beta_PR) + q.s.upper3[,4] * sqrt(A_PR))) & (theta >= (c(beta_PR) + q.s.lower3[,4] * sqrt(A_PR))))
    cov.FH3 <- ((theta <= (c(eblup.FH) + qnorm(1-Alpha[3]/2) * sqrt(mse_FH))) & (theta >= (c(eblup.FH) + qnorm(Alpha[3]/2) * sqrt(mse_FH))))
    cov.PR3 <- ((theta <= (c(eblup.PR) + qnorm(1-Alpha[3]/2) * sqrt(mse_PR))) & (theta >= (c(eblup.PR) + qnorm(Alpha[3]/2) * sqrt(mse_PR))))
    cov.DIR3 <- ((theta <= (c(Direct) + qnorm(1-Alpha[3]/2) * sqrt(d))) & (theta >= (c(Direct) + qnorm(Alpha[3]/2) * sqrt(d))))
    
    cov.DB.FH3 <- ((theta <= (eblup.FH + q.d.upper3[,1] * sqrt(g1_FH))) & (theta >= (eblup.FH + q.d.lower3[,1] * sqrt(g1_FH))))
    cov.DB.PR3 <- ((theta <= (eblup.PR + q.d.upper3[,2] * sqrt(g1_PR))) & (theta >= (eblup.PR + q.d.lower3[,2] * sqrt(g1_PR))))
    
    coverge3 <- matrix(0, m, 9) 
    
    coverge3[, 1] <- ifelse(cov.CLL.FH3, 1, 0)
    coverge3[, 2] <- ifelse(cov.HM.FH3, 1, 0)
    coverge3[, 3] <- ifelse(cov.DB.FH3, 1, 0)
    
    coverge3[, 4] <- ifelse(cov.CLL.PR3, 1, 0)
    coverge3[, 5] <- ifelse(cov.HM.PR3, 1, 0)
    coverge3[, 6] <- ifelse(cov.DB.PR3, 1, 0)
    
    coverge3[, 7] <- ifelse(cov.FH3, 1, 0)
    coverge3[, 8] <- ifelse(cov.PR3, 1, 0)
    coverge3[, 9] <- ifelse(cov.DIR3, 1, 0)
    
    Int3 <- coverge3 + Int3
    
    if(s%%100==0){print(s)}
  }
)



G1 <- colMeans(Int1[1:3,])
G2 <- colMeans(Int1[4:6,])
G3 <- colMeans(Int1[7:9,])
G4 <- colMeans(Int1[10:12,])
G5 <- colMeans(Int1[13:15,])
rbind(G1, G2, G3, G4, G5)

a1 = rbind(rbind(G1, G2, G3, G4, G5))/Nosim*100
colnames(a1) = c("CLL.FH","HM.FH","DB.FH","CLL.PR","HM.PR", "DB.PR", "FH","PR","DIRECT")

G1 <- colMeans(qdif1[1:3,])
G2 <- colMeans(qdif1[4:6,])
G3 <- colMeans(qdif1[7:9,])
G4 <- colMeans(qdif1[10:12,])
G5 <- colMeans(qdif1[13:15,])
rbind(G1, G2, G3, G4, G5)

b1 = rbind(rbind(G1, G2, G3, G4, G5))/Nosim
colnames(b1) = c("CLL.FH","HM.FH","DB.FH","CLL.PR","HM.PR", "DB.PR", "FH","PR","DIRECT")

G1 <- colMeans(Int2[1:3,])
G2 <- colMeans(Int2[4:6,])
G3 <- colMeans(Int2[7:9,])
G4 <- colMeans(Int2[10:12,])
G5 <- colMeans(Int2[13:15,])
rbind(G1, G2, G3, G4, G5)

a2 = rbind(rbind(G1, G2, G3, G4, G5))/Nosim*100
colnames(a2) = c("CLL.FH","HM.FH","DB.FH","CLL.PR","HM.PR", "DB.PR", "FH","PR","DIRECT")

G1 <- colMeans(qdif2[1:3,])
G2 <- colMeans(qdif2[4:6,])
G3 <- colMeans(qdif2[7:9,])
G4 <- colMeans(qdif2[10:12,])
G5 <- colMeans(qdif2[13:15,])
rbind(G1, G2, G3, G4, G5)

b2 = rbind(rbind(G1, G2, G3, G4, G5))/Nosim
colnames(b2) = c("CLL.FH","HM.FH","DB.FH","CLL.PR","HM.PR", "DB.PR", "FH","PR","DIRECT")

G1 <- colMeans(Int3[1:3,])
G2 <- colMeans(Int3[4:6,])
G3 <- colMeans(Int3[7:9,])
G4 <- colMeans(Int3[10:12,])
G5 <- colMeans(Int3[13:15,])
rbind(G1, G2, G3, G4, G5)

a3 = rbind(rbind(G1, G2, G3, G4, G5))/Nosim*100
colnames(a3) = c("CLL.FH","HM.FH","DB.FH","CLL.PR","HM.PR", "DB.PR", "FH","PR","DIRECT")

G1 <- colMeans(qdif3[1:3,])
G2 <- colMeans(qdif3[4:6,])
G3 <- colMeans(qdif3[7:9,])
G4 <- colMeans(qdif3[10:12,])
G5 <- colMeans(qdif3[13:15,])
rbind(G1, G2, G3, G4, G5)

b3 = rbind(rbind(G1, G2, G3, G4, G5))/Nosim
colnames(b3) = c("CLL.FH","HM.FH","DB.FH","CLL.PR","HM.PR", "DB.PR", "FH","PR","DIRECT")

tab2 = rbind(cbind(a1, b1),cbind(a2, b2),cbind(a3, b3))
xtable::xtable(tab2)


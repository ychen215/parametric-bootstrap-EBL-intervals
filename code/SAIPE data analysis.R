saipet = read.table("CPS89T.TXT", header = T)
data5_17 = saipet[, c(1, 9:15, 30)]
library(mvtnorm)
library(ggplot2)
source("basic functions.R")
B =  400
DB = 100
## number of negative estimates in single bootstrap
s_fh = 0
s_pr = 0
## number of negative estimates in double bootstrap
ss_fh = 0
ss_pr = 0
# Data in modeling
y = data5_17$cps89.1
d = (data5_17$fnlse.1)^2
m = nrow(data5_17)
Alpha = 0.1
df = 3
x = as.matrix(data5_17[, 3:6])
p = dim(x)[2]

#---------Fay-Herriot Method -----------
est_A = estimate_A(method = "FH", x, y, d, m, p)
A_FH = est_A$A
all_FH = compute_eblup(A_FH, x, y, d)
Vi_FH = all_FH$Vi
beta_FH = all_FH$beta
B_FH = all_FH$B
eblup.FH = all_FH$eblup
g1_FH = d * (1 - B_FH)

Z.j.FH = matrix(0, B, m)
pivot_FH = matrix(0, B, m)
####------Parametric Bootstrap under t3 distribution------
for (bb in 1:B){
  bv_FH = rt(m, df)*sqrt((df-2)*c(A_FH)/df)
  btheta_FH = x%*%beta_FH + bv_FH
  
  be = rmvnorm(1, sigma = diag(d))
  
  by_FH = btheta_FH + matrix(be)
  
  #Single bootstrap for FH
  est_A_star = estimate_A(method = "FH", x, by_FH, d, m, p)
  A.star.FH = est_A_star$A
  s_fh = s_fh + est_A_star$k0
  
  all_FH_star = compute_eblup(A.star.FH, x, by_FH, d)
  
  Vi.star.FH = all_FH_star$Vi
  beta.s.FH = all_FH_star$beta
  B.star.FH = all_FH_star$B
  eblup.s.FH = all_FH_star$eblup
  g1.s.FH = d * (1 - B.star.FH)
  
  pivot_FH[bb, ] = t((btheta_FH - eblup.s.FH) / sqrt(g1.s.FH))
  
  d_pivot_FH = matrix(0, DB, m)
  d_HM_FH = matrix(0, DB, m)
  
  for (dd in 1:DB) {
    dv_FH = rt(m, df)*sqrt((df-2)*c(A.star.FH)/df)
    dtheta_FH = x%*%beta.s.FH + dv_FH
    de = rmvnorm(1, sigma = diag(d))
    dy_FH = dtheta_FH + matrix(de)
    
    # Double bootstrap for FH
    est_A_ss = estimate_A(method = "FH", x, dy_FH, d, m, p)
    A.ss.FH = est_A_ss$A
    ss_fh = ss_fh + est_A_ss$k0
    
    all_FH_ss = compute_eblup(A.ss.FH, x, dy_FH, d)
    Vi.ss.FH = all_FH_ss$Vi
    beta.ss.FH = all_FH_ss$beta
    B.ss.FH = all_FH_ss$B
    eblup.ss.FH = all_FH_ss$eblup
    g1.ss.FH = d * (1 - B.ss.FH)
    
    d_pivot_FH[dd, ] = t((dtheta_FH - eblup.ss.FH) / sqrt(g1.ss.FH))
  }
  eblup.FH.mat = matrix(rep(pivot_FH[bb, ], rep(DB, m)), DB, m, byrow = F)
  Z.j.FH[bb, ] = colSums(d_pivot_FH<=eblup.FH.mat)/DB
  
  plot(1,ylim=c(0,1),main=paste(1,"/",bb))
}

q.s.lower1 = as.numeric(sapply(as.data.frame(pivot_FH), quantile, prob=Alpha/2))
q.s.upper1 = as.numeric(sapply(as.data.frame(pivot_FH), quantile, prob=(1-Alpha/2)))

au.SB.FH1 = apply(Z.j.FH, 2, quantile, prob=1-Alpha/2)
al.SB.FH1 = apply(Z.j.FH, 2, quantile, prob=Alpha/2)
q.d.upper1 <- matrix(0, m, 1)
q.d.lower1 <- matrix(0, m, 1)
for(i in 1:m){
  q.d.lower1[i,]<-c(quantile(pivot_FH[, i], prob=al.SB.FH1[i]))
  q.d.upper1[i,]<-c(quantile(pivot_FH[, i], prob=au.SB.FH1[i]))
}
SB.fh1 = cbind(eblup.FH + q.s.lower1 * sqrt(g1_FH), eblup.FH + q.s.upper1 * sqrt(g1_FH))
DB.FH1 = cbind(eblup.FH + q.d.lower1 * sqrt(g1_FH), eblup.FH + q.d.upper1 * sqrt(g1_FH))

intervals_t3 = data.frame(state = data5_17$fips,
                          EBLUP = eblup.FH,
                          SB.t3.l = eblup.FH + q.s.lower1 * sqrt(g1_FH),
                          SB.t3.u = eblup.FH + q.s.upper1 * sqrt(g1_FH),
                          DB.t3.l = eblup.FH + q.d.lower1 * sqrt(g1_FH),
                          DB.t3.u = eblup.FH + q.d.upper1 * sqrt(g1_FH))
intervals_t3$model = rep("t", 51)

# Single and double bootstrap intervals under t3--------------------------------
intervals_t3_SB = intervals_t3[, 1:4]
colnames(intervals_t3_SB) = c("state", "estimates", "lower", "upper")
intervals_t3_SB$model = rep("t,SB", 51)

intervals_t3_DB = intervals_t3[, c(1:2, 5:6)]
colnames(intervals_t3_DB) = c("state", "estimates", "lower", "upper")
intervals_t3_DB$model = rep("t,DB", 51)

# Direct intervals--------------------------------------------------------------
intervals_dir = data.frame(state = data5_17$fips,
                           estimates = data5_17$cps89.1,
                           lower = data5_17$cps89.1 - qnorm(0.95) * data5_17$fnlse.1,
                           upper = data5_17$cps89.1 + qnorm(0.95) * data5_17$fnlse.1)
intervals_dir$model = rep("DIRECT", 51)

# Parametric bootstrap under normal distribution--------------------------------
pivot_FH = matrix(0, B, m)
for (bb in 1:B){
  # bv_FH = rt(m, df)*sqrt((df-2)*c(A_FH)/df)
  bv_FH = rnorm(m, 0, sqrt(A_FH))
  btheta_FH = x%*%beta_FH + bv_FH
  
  be = rmvnorm(1, sigma = diag(d))
  
  by_FH = btheta_FH + matrix(be)
  
  #Single bootstrap with FH estimate
  est_A_star = estimate_A(method = "FH", x, by_FH, d, m, p)
  A.star.FH = est_A_star$A
  s_fh = s_fh + est_A_star$k0
  
  all_FH_star = compute_eblup(A.star.FH, x, by_FH, d)
  
  Vi.star.FH = all_FH_star$Vi
  beta.s.FH = all_FH_star$beta
  B.star.FH =all_FH_star$B
  eblup.s.FH = all_FH_star$eblup
  g1.s.FH = d * (1 - B.star.FH)
  
  pivot_FH[bb, ] = t((btheta_FH - eblup.s.FH) / sqrt(g1.s.FH))
  
  plot(1,ylim=c(0,1),main=paste(1,"/",bb))
}

q.s.lower2 = as.numeric(sapply(as.data.frame(pivot_FH), quantile, prob=Alpha/2))
q.s.upper2 = as.numeric(sapply(as.data.frame(pivot_FH), quantile, prob=(1-Alpha/2)))

intervals_norm_SB = data.frame(state = data5_17$fips,
                               estimates = eblup.FH,
                               lower = eblup.FH + q.s.lower2 * sqrt(g1_FH),
                               upper = eblup.FH + q.s.upper2 * sqrt(g1_FH))
intervals_norm_SB$model = rep("normal", 51)

fips = read.csv("us_state_fips_codes.csv")
final_ints = rbind(intervals_t3_SB, 
                   intervals_t3_DB, 
                   intervals_norm_SB,
                   intervals_dir)
final_ints$fips = final_ints$state

merge_df = merge(final_ints, fips, by.x = "fips")
merge_df$model = factor(merge_df$model, c("DIRECT", "normal", "t,SB", "t,DB"))


ggplot(merge_df, aes(x = factor(State), y = estimates, shape = model, colour = model)) +
  geom_point(position = position_dodge(width = 1), size = 2.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, linetype = model),
                position = position_dodge(width = 0.8),
                width = 1) +
  scale_shape_manual(values = c(15, 17, 19, 21)) +  # different shapes for methods
  theme_minimal() +
  xlab("State") +
  ylab("Estimate with 90% Prediction Interval")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())
save.image("SAIPE data analysis.RData")

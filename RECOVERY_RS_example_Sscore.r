
library(mvtnorm)
library(parallel)

#rm(list=ls())

#-------------------------------------------------------------------------------
conditional_error <- function (look_size, upper, observed_statistic){
  n.looks = length(look_size)- length(observed_statistic) + 2
  s_var = matrix(rep(NA, n.looks^2), ncol=n.looks)
  for (i in seq(1, n.looks)){
    for (j in seq(1,i)){
      s_var[i,j] = s_var[j,i] = look_size[j+length(observed_statistic)-2] 
    }
  }
  #print(s_var)
  counter = length(observed_statistic)-1
  CTI_err = c(rep(NA, 2))
  J.sigma = matrix(rep(NA, (n.looks-1)^2), ncol=(n.looks-1))
  J.mean = c(rep(0, (n.looks-1)))
  sigma_11 = matrix(rep(NA, 1), ncol=1)
  sigma_12 = matrix(NA, nrow=1, ncol=(n.looks-1))
  sigma_21 = matrix(NA, nrow=(n.looks-1), ncol=1)
  sigma_22 = matrix(rep(NA, (n.looks-1)^2), ncol=(n.looks-1))
  sigma_11 = s_var[1,1]
  sigma_12 = s_var[1,2:n.looks]
  sigma_21 = s_var[2:n.looks,1]
  sigma_22 = s_var[2:n.looks,2:n.looks]
  J.sigma = sigma_22 - sigma_21%*%solve(sigma_11)%*%sigma_12
  J.mean = sigma_21%*%solve(sigma_11)%*%observed_statistic[counter]
  
  if (length(observed_statistic) == length(look_size)) {
    CTI_err[2] = 1 - pnorm(upper[length(observed_statistic)], mean = J.mean, sd = sqrt(J.sigma))
  }
  else {
    u = c(upper[length(observed_statistic):length(look_size)])
    #l = c(rep(-Inf, (n.looks-counter)))
    CTI_err[2] = 1-pmvnorm(upper = u, mean = J.mean[1:(n.looks-1)], sigma = J.sigma[1:(n.looks-1),1:(n.looks-1)], abseps=0.00001,maxpts=10000)[1]
  }
  CTI_err[1] = observed_statistic[counter]
  #print(CTI_err_new)
  return(CTI_err)
  
}

#-------------------------------------------------------------------------------
search.function <- function(u, target, mean, sigma){
  ( 1 - pnorm(u, mean = mean, sd = sqrt(sigma)) - target)
}

#-------------------------------------------------------------------------------
CE_boundary <- function(look_size, CTI_err, observed_statistic){
  n.looks = 2#(length(observed_statistic)-1)+1 #same as CTI_err size
  s_var = matrix(rep(NA, n.looks^2), ncol=n.looks)
  for (i in seq(1, n.looks)){
    for (j in seq(1,i)){
      s_var[i,j] = s_var[j,i] = look_size[j+length(observed_statistic)-2] 
    }
  }
  #print(s_var)
  counter = length(observed_statistic)-1
  #CTI_err = c(rep(NA, counter+1))
  J.sigma = c(rep(NA, 1))
  J.mean = c(rep(NA, 1))
  sigma_11 = c(rep(NA, 1))#matrix(rep(NA, counter^2), ncol=counter)
  sigma_12 = c(rep(NA, 1))#matrix(NA, nrow=counter , ncol=(n.looks-counter))
  sigma_21 = c(rep(NA, 1))#matrix(NA, nrow=(n.looks-counter), ncol=counter)
  sigma_22 = c(rep(NA, 1))#matrix(rep(NA, (n.looks-counter)^2), ncol=(n.looks-counter))
  sigma_11 = s_var[1,1]
  sigma_12 = s_var[1,2]
  sigma_21 = s_var[2,1]
  sigma_22 = s_var[2,2]
  J.sigma = sigma_22 - sigma_21%*%solve(sigma_11)%*%sigma_12
  J.mean = sigma_21%*%solve(sigma_11)%*%observed_statistic[counter]
  #print(J.sigma)
  #print(J.mean)
  new_upper = uniroot(search.function, target = CTI_err[2], mean = J.mean, sigma = J.sigma, interval=c(-50,50)*sqrt(look_size[counter+1]))$root
  #print(new_upper_new)
  return(new_upper) 
}
#-------------------------------------------------------------------------------
MSN.search <- function(Max_size, target, Max_looks, org_z.boundaries){
  # calculating power under the alternative using z scores. see your notes in the booklet. 
  tail(Power(seq(1,Max_looks)/Max_looks * Max_size, org_z.boundaries, lower = rep(-Inf, Max_looks), theta = delta), n=1 ) - target
}

#-------------------------------------------------------------------------------

design <- function(look_size, alpha.star.u, alpha.star.l)
{
  n.looks = length(look_size)
  alpha.star.u.inc = increments(alpha.star.u)
  alpha.star.l.inc = increments(alpha.star.l)
  
  upper = rep(Inf, n.looks)
  lower = rep(-Inf, n.looks)
  for (look in seq(1, n.looks))
  {
    if (alpha.star.u.inc[look] > 0) {
      upper[look] = uniroot(upper_boundaries, c(-50,50)*sqrt(look_size[look]), look, look_size, upper, lower, target=alpha.star.u[look])$root
    }
    if (alpha.star.l.inc[look] > 0) {
      lower[look] = uniroot(lower_boundaries, c(-50,50)*sqrt(look_size[look]), look, look_size, upper, lower, target=alpha.star.l[look])$root
    }
  }
  data.frame(upper, lower)
}

#-------------------------------------------------------------------------------

increments <- function(v){
  c(v[1], diff(v))
}

#-------------------------------------------------------------------------------

upper_boundaries <- function(upper.value, look, look_size, upper, lower, target)
{
  look_size = look_size[1:look]
  lower = lower[1:look]
  upper = upper[1:look]
  upper[look] = upper.value
  Power(look_size, upper, lower, theta=0)[look] - target
}


lower_boundaries <- function(lower.value, look, look_size, upper, lower, target)
{
  look_size = look_size[1:look]
  upper = upper[1:look]
  lower = lower[1:look]
  lower[look] = lower.value
  Power(look_size, upper=-lower, lower=-upper, theta=0)[look] - target
}

#-------------------------------------------------------------------------------

Power <- function(look_size, upper, lower, theta){
  # calculating the z scores when interim looks are equally spaced does not depend on look sizes as they cancel out one another.
  n.looks = length(look_size)
  s_mean = theta*sqrt(look_size)
  s_var = matrix(rep(NA, n.looks^2), ncol=n.looks)
  
  for (i in seq(1, n.looks)){
    for (j in seq(1,i)){
      s_var[i,j] = s_var[j,i] = sqrt(j)/sqrt(i)
    }
  }
  look_power = 0
  pr_look = rep(0,n.looks)
  for (look in seq(1,n.looks)){
    pr_look[look] = Pr_stopping(look, upper, lower, s_mean, s_var)
  } 
  
  cumsum(pr_look)
}


#-------------------------------------------------------------------------------

Pr_stopping <- function(look, upper, lower, mean, Sigma){
  if (look == 1) {
    pr = 1 - pnorm(upper[1], mean = mean[1], sd = sqrt(Sigma[1,1]))
  }
  else {
    u = c(upper[1:(look-1)],Inf)
    l = c(lower[1:(look-1)],upper[look])
    pr = pmvnorm(lower = l, upper = u, mean = mean[1:look], sigma = Sigma[1:look,1:look], abseps=0.00001,maxpts=10000)[1]
  }
  pr
}
#-------------------------------------------------------------------------------

p_C = 0.15 # standard therapy
p_E = 0.10 # SPAP
#we are interested in reduction in intubation when it comes to the experimental arm.
p_mean = (p_E + p_C)/ 2
#print(p_mean)
Max_looks = 12
Max_power = 0.9
alpha = 0.05

delta = -log((p_E*(1-p_C))/(p_C*(1-p_E)))# due to pc and pe values
#print(delta)
look_space = seq(1,Max_looks)/Max_looks 
alpha.star.u = (alpha/2)*look_space
alpha.star.l = rep(0, Max_looks)#alpha.star.u
print(alpha.star.u)

get_design = design(look_space, alpha.star.u, alpha.star.l)
org_z.boundaries = get_design$upper
print(org_z.boundaries)


Max_info = uniroot(MSN.search, Max_power, Max_looks, org_z.boundaries, interval = c(0,10000))$root
print(Max_info)
Max_size = 4 * Max_info / (p_mean * (1 - p_mean))# this is less reliable as it depends on P_mean
#print(Max_size)
look_size = look_space*Max_info
print(look_size)
#-------------------------------------------------------------------------------

#so we prefer to calculate the info using observed data
#(s_E, s_C, f_E, f_C)-----> s_E|f_E
#                           -------
#                           s_C|f_C
data_1 = matrix(c(12, 11, 10, 12), ncol=2)
v_1 = (data_1[1,1] + data_1[1,2])*(data_1[2,1] + data_1[2,2])*(data_1[1,1] + data_1[2,1])*(data_1[1,2] + data_1[2,2])/(sum(data_1)^3)
data_2 = matrix(c(19, 15, 19, 13), ncol=2)
v_2 = (data_2[1,1] + data_2[1,2])*(data_2[2,1] + data_2[2,2])*(data_2[1,1] + data_2[2,1])*(data_2[1,2] + data_2[2,2])/(sum(data_2)^3)
data_3 = matrix(c(83, 81, 43, 41), ncol=2)
v_3 = (data_3[1,1] + data_3[1,2])*(data_3[2,1] + data_3[2,2])*(data_3[1,1] + data_3[2,1])*(data_3[1,2] + data_3[2,2])/(sum(data_3)^3)
data_4 = matrix(c(240, 198, 137, 158), ncol=2)
v_4 = (data_4[1,1] + data_4[1,2])*(data_4[2,1] + data_4[2,2])*(data_4[1,1] + data_4[2,1])*(data_4[1,2] + data_4[2,2])/(sum(data_4)^3)
#print(data_1)
#print(data_2)
#print(data_3)
#print(data_4)
#print(v_1)
#print(v_2)
#print(v_3)
#print(v_4)

#look_size = c(sum(data_1), sum(data_2), sum(data_3), sum(data_4))* Max_info /  Max_size
#new_alpha.star.u_1 = (alpha/2)*c(sum(data_1), sum(data_2), sum(data_3), sum(data_4)) / Max_size
#print(look_size_1)


#-------------------------------------------------------------------------------

new_look_size = c(v_1, v_2, v_3, v_4) 
print(new_look_size)
new_alpha.star.u = (alpha/2)*(c(v_1, v_2, v_3, v_4) / Max_info)
print(new_alpha.star.u)
get_new_design = design(new_look_size, new_alpha.star.u, alpha.star.l)
new_z.boundaries = get_new_design$upper
print(new_z.boundaries*sqrt(new_look_size))
#-------------------------------------------------------------------------------
observed_statistic = c(((data_1[2,1] + data_1[2,2])*data_1[1,1] - (data_1[1,1] + data_1[1,2])*data_1[2,1])/sum(data_1),
                       ((data_2[2,1] + data_2[2,2])*data_2[1,1] - (data_2[1,1] + data_2[1,2])*data_2[2,1])/sum(data_2), 
                       ((data_3[2,1] + data_3[2,2])*data_3[1,1] - (data_3[1,1] + data_3[1,2])*data_3[2,1])/sum(data_3), 
                       ((data_4[2,1] + data_4[2,2])*data_4[1,1] - (data_4[1,1] + data_4[1,2])*data_4[2,1])/sum(data_4)) #/ sqrt(new_look_size)
print(observed_statistic)
#-------------------------------------------------------------------------------
rem_alpha.star.u = c((alpha/2)*v_1/ Max_info, (alpha/2)*v_2/ Max_info, (alpha/2)*v_3/ Max_info, (alpha/2))
#print(rem_alpha.star.u)
get_rem_design = design(new_look_size, rem_alpha.star.u, alpha.star.l)
rem_z.boundaries = get_rem_design$upper
print(rem_z.boundaries*sqrt(new_look_size))

#-------------------------------------------------------------------------------

CE_look_size = c(rep(0, Max_looks))
CE_alpha.star.u = c(rep(0, Max_looks))
for (i in seq(1,Max_looks)){
  if (i < 4){
    CE_look_size[i] = new_look_size[i]
    CE_alpha.star.u[i] = new_alpha.star.u[i]
  }
  else{
    CE_look_size[i] = look_space[i]*Max_info
    CE_alpha.star.u[i] = alpha.star.u[i]
  }
}
print(CE_look_size)
#print(CE_alpha.star.u)
get_CE_design = design(CE_look_size, CE_alpha.star.u, alpha.star.l)
CE_z.boundaries = get_CE_design$upper
print(CE_z.boundaries)
CTIe = conditional_error(CE_look_size, CE_z.boundaries*sqrt(CE_look_size), observed_statistic)#*sqrt(new_look_size))#$CTI_err
print(CTIe)
CE_boundary = CE_boundary(new_look_size, CTIe, observed_statistic)
print(CE_boundary)#/sqrt(new_look_size[4]))

print(observed_statistic)#/sqrt(new_look_size))


#---------------------------------------------------------------------------------

#xdata <- c(look_size)
xdata <- c(2.809723, 4.026602,13.883484,44.032724)
#y1 <- c(org_z.boundaries)
y1 <- c(new_z.boundaries*sqrt(new_look_size))
y2 <- c(observed_statistic)#/sqrt(new_look_size))
y3 <- c(rem_z.boundaries[4]*sqrt(new_look_size[4]))
y4 <- c(CE_boundary)#/sqrt(new_look_size[4]))
par(mar=c(4.1, 4.1, 4.1, 4.1), xpd=TRUE)
plot(xdata[order(xdata)], y1[order(xdata)], type = "o", pch = 3, col = "darkgreen", lty=4, xlab = "V", ylab = "S", xlim=c(0,50), ylim=c(-5,15),  cex.lab=0.75, cex.axis=0.75)
lines(xdata[order(xdata)], y1[order(xdata)], type = "o", pch = 3, col = "darkgreen", lty=4)
#points(c(2.809723, 4.026602,13.883484,44.032724), y2, type = "o", pch = 20, col = "darkgreen",lty=2)
points(c(44.032724), y4, type = "o", pch = 2, lty=3, col = "orchid3")
points(c(44.032724), y3, type = "o", pch = 1, lty=6, col = "darkred")
points(c(2.809723, 4.026602,13.883484,44.032724), y2, pch = "*", col = "black", lty=1)
#lines(c(2.809723, 4.026602,13.883484,44.032724), y2, pch = "*", col = "black", lty=1)
#legend(x="topleft", legend=c("Using original spending function", "Using conditional error", "Spending all remaining alpha", "observed statistic"), col=c("darkgreen", "orchid3", "darkred", "black"),  pch=c(19,20,18,8), cex =0.4, ncol=1)




# Software (c) 2021, The University of Warwick (the "Software"), comprising of code written in R. 
#
# The Software remains the property of the University of Warwick ("the University"). 
# 
# The Software is distributed "AS IS" under this Licence solely for non-commercial use, 
# such as academic research, clinical research, and clinical trials, in the hope that it 
# will be useful, but in order that the University as a charitable foundation protects 
# its assets for the benefit of its educational and research purposes, the University  
# makes clear that no condition is made or to be implied, nor is any warranty given or to  
# be implied, as to the accuracy of the Software, or that it will be suitable for any  
# particular purpose or for use under any specific conditions. Furthermore, the University  
# disclaims all responsibility for the use which is made of the Software. It further  
# disclaims any liability for the outcomes arising from using the Software. 
# 
# The Licensee agrees to indemnify the University and hold the University harmless from  
# and against any and all claims, damages and liabilities asserted by third parties  
# (including claims for negligence) which arise directly or indirectly from the use of  
# the Software. No part of the Software may be reproduced, modified, transmitted or  
# transferred in any form or by any means, electronic or mechanical, without the  
# express permission of the University. The permission of the University is not required  
# if the said reproduction, modification, transmission or transference is done without  
# financial return, the conditions of this Licence are imposed upon the receiver of the  
# product, and all original and amended source code is included in any transmitted product.  
# You may be held legally responsible for any copyright infringement that is caused or  
# encouraged by your failure to abide by these terms and conditions.
# 
# You are not permitted under this Licence to use this Software commercially. Use for which  
# any financial return is received shall be defined as commercial use, and includes  
# (1) integration of all or part of the source code or the Software into a product for sale  
# or license by or on behalf of Licensee to third parties or (2) use of the Software or any  
# derivative of it for research with the final aim of developing software products for  
# sale or license to a third party or (3) use of the Software or any derivative of it for  
# research with the final aim of developing non-software products for sale or license to  
# a third party, or (4) use of the Software to provide any service to an external organisation  
# for which payment is received.   
# 
# If you are interested in using the Software commercially, please contact Warwick Ventures  
# Limited ("WVL"), the technology transfer company of the University, to negotiate a licence.  
# Contact details are: ventures@warwick.ac.uk.  
#
# If you are in any doubt if your use constitutes commercial use, contact WVL.

library(mvtnorm)
library(parallel)

#rm(list=ls())
#-------------------------------------------------------------------------------
do_simulation <- function(sim=1, look_size, org_uppers, sprt_org_uppers, remaining_uppers, sprt_remaining_uppers, seed){
  #if (!missing("seed")) set.seed(seed+sim)
  
  org_sim_result = c(rep(0, Max_looks))
  rem_sim_result = c(rep(0, Max_looks))
  CE_sim_result = c(rep(0, Max_looks))
  
  sprt_org_sim_result = c(rep(0, Max_looks))
  sprt_rem_sim_result = c(rep(0, Max_looks))
  sprt_CE_sim_result = c(rep(0, Max_looks))
  
  observed_statistic = statistic_calculator(look_size, length(look_size)-1, 0)
  #print(observed_statistic)#$observed_statistic
  
  org_boundary = org_uppers
  sprt_org_boundary = sprt_org_uppers
  
  rem_boundary = remaining_uppers
  sprt_rem_boundary = sprt_remaining_uppers
  
  
  stop_marker = 0
  sprt_stop_marker = 0
  
  for (j in seq(1, length(look_size))){
    if (stop_marker == 0){
      if(j == 1){
        if (((org_boundary[j] - observed_statistic[j]) > 0) == TRUE){
          org_sim_result[j] = 1
          rem_sim_result[j] = 1
          CE_sim_result[j] = 1
        }
        else {
          org_sim_result[j] = -1
          rem_sim_result[j] = -1
          CE_sim_result[j] = -1
          stop_marker = 1
        }
      }
      else{
        #print("hi")
        ter_CTIe = termination_conditional_error(look_size, remaining_uppers, observed_statistic[1:j])
        org_CTIe = conditional_error(look_size, org_boundary, observed_statistic[1:j])
        CE_boundary = CE_boundary(look_size, org_CTIe, observed_statistic[1:j])
        #print(ter_CTIe[2])
        #print(org_CTIe[2])
        #print(CE_boundary)
        
        #if(abs(ter_CTIe[2] - org_CTIe[2]) > 0){
        if(ter_CTIe[2] > org_CTIe[2]){
          for (i in 1:3){
            if (i == 1){
              if (((org_boundary[j] - observed_statistic[j]) > 0) == TRUE){
                org_sim_result[j] = 1
              }
              else {
                org_sim_result[j] = -1
              }
            }
            else if (i == 2){
              if (((rem_boundary[j] - observed_statistic[j]) > 0) == TRUE){
                rem_sim_result[j] = 1
              }
              else {
                rem_sim_result[j] = -1
              }
            }
            else{
              if (((CE_boundary - observed_statistic[j]) > 0) == TRUE){
                CE_sim_result[j] = 1
              }
              else {
                CE_sim_result[j] = -1
              }
            }
          }
          stop_marker = 1
        }
        else{
          if (((org_boundary[j] - observed_statistic[j]) > 0) == TRUE){
            org_sim_result[j] = 1
            rem_sim_result[j] = 1
            CE_sim_result[j] = 1
          }
          else {
            org_sim_result[j] = -1
            rem_sim_result[j] = -1
            CE_sim_result[j] = -1
            stop_marker = 1
          }
        }
      }
    }
  }
  for (j in seq(1, length(look_size))){
    if (sprt_stop_marker == 0){
      if(j == 1){
        if (((sprt_org_boundary[j] - observed_statistic[j]) > 0) == TRUE){
          sprt_org_sim_result[j] = 1
          sprt_rem_sim_result[j] = 1
          sprt_CE_sim_result[j] = 1
        }
        else {
          sprt_org_sim_result[j] = -1
          sprt_rem_sim_result[j] = -1
          sprt_CE_sim_result[j] = -1
          sprt_stop_marker = 1
        }
      }
      else{
        #print("hi")
        sprt_ter_CTIe = termination_conditional_error(look_size, sprt_remaining_uppers, observed_statistic[1:j])
        sprt_org_CTIe = sprt_CE(look_size, observed_statistic[1:j], j)#$CTI_err
        sprt_CE_boundary = CE_boundary(look_size, sprt_org_CTIe, observed_statistic[1:j])
        #print(ter_CTIe[2])
        #print(org_CTIe[2])
        #print(CE_boundary)
        
        #if(abs(ter_CTIe[2] - org_CTIe[2]) > 0){
        if(sprt_ter_CTIe[2] > sprt_org_CTIe[2]){
          for (i in 1:3){
            if (i == 1){
              if (((sprt_org_boundary[j] - observed_statistic[j]) > 0) == TRUE){
                sprt_org_sim_result[j] = 1
              }
              else {
                sprt_org_sim_result[j] = -1
              }
            }
            else if (i == 2){
              if (((sprt_rem_boundary[j] - observed_statistic[j]) > 0) == TRUE){
                sprt_rem_sim_result[j] = 1
              }
              else {
                sprt_rem_sim_result[j] = -1
              }
            }
            else{
              if (((sprt_CE_boundary - observed_statistic[j]) > 0) == TRUE){
                sprt_CE_sim_result[j] = 1
              }
              else {
                sprt_CE_sim_result[j] = -1
              }
            }
          }
          sprt_stop_marker = 1
        }
        else{
          if (((sprt_org_boundary[j] - observed_statistic[j]) > 0) == TRUE){
            sprt_org_sim_result[j] = 1
            sprt_rem_sim_result[j] = 1
            sprt_CE_sim_result[j] = 1
          }
          else {
            sprt_org_sim_result[j] = -1
            sprt_rem_sim_result[j] = -1
            sprt_CE_sim_result[j] = -1
            sprt_stop_marker = 1
          }
        }
      }
    }
  }
  sim_result = data.frame(org_sim_result, rem_sim_result, CE_sim_result, sprt_org_sim_result, sprt_rem_sim_result, sprt_CE_sim_result)
  #print(sim_result)
}

#-------------------------------------------------------------------------------

run_simulation <- function(simulation_fn, no_simulation, look_size, org_uppers, sprt_org_uppers, remaining_uppers, sprt_remaining_uppers, seed){# runs simulation.fn nsim times and gives output in matrix with one row per simulation
  # seed is optional for reproducible simulations
  sim_result_total = mclapply(seq(1, no_simulation), simulation_fn, look_size=look_size, org_uppers=org_uppers, sprt_org_uppers=sprt_org_uppers, remaining_uppers=remaining_uppers, sprt_remaining_uppers=sprt_remaining_uppers, seed=seed)
  sim_result_mat = matrix(unlist(sim_result_total), ncol=no_simulation)
  t(sim_result_mat)
  #print(t(sim_result_mat))
}
#--
#-------------------------------------------------------------------------------
statistic_calculator <- function(look_size, n.stage, theta){
  n.look = length(look_size)
  if (n.stage < n.look){
    observed_statistic = c(rep(NA, n.stage+1))
    for (look in seq(1, n.stage+1)){
      if (look == 1){
        s_mean = theta*look_size[look]
        s_var = look_size[look]
        observed_statistic[look] = rnorm(1, mean = s_mean, sd = sqrt(s_var))
      }
      else{
        s_mean = theta*(look_size[look] - look_size[look-1])
        s_var = look_size[look] - look_size[look-1]
        observed_statistic[look] = observed_statistic[look-1] + rnorm(1, mean = s_mean, sd = sqrt(s_var))#, mean = 0, sd = cumsum(look_size[1]:look_size[j]))        }
      }
    }
  }
  else{
    print("Oopss!")
  }
  #print(observed_statistic)
  return(observed_statistic)
}
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
termination_conditional_error <- function (look_size, upper, observed_statistic){
  n.looks = 2#length(look_size)- length(observed_statistic) + 2
  s_var = matrix(rep(NA, n.looks^2), ncol=n.looks)
  for (i in seq(1, n.looks)){
    for (j in seq(1,i)){
      s_var[i,j] = s_var[j,i] = look_size[j+length(observed_statistic)-2] 
    }
  }
  counter = length(observed_statistic)-1
  #print(s_var)
  #print(upper[length(observed_statistic),length(observed_statistic)])
  CTI_err = c(rep(NA, 2))
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
  #print(sigma_11)
  #print(sigma_12)
  #print(sigma_21)
  #print(sigma_22)
  
  CTI_err[2] = 1 - pnorm(upper[length(observed_statistic)], mean = J.mean, sd = sqrt(J.sigma))
  
  CTI_err[1] = observed_statistic[counter]
  #print(CTI_err)
  return(CTI_err)
}

#-------------------------------------------------------------------------------
search.function <- function(u, target, mean, sigma){
  ( 1 - pnorm(u, mean = mean, sd = sqrt(sigma))) - target
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
  tail(Power(seq(1,Max_looks)/Max_looks * Max_size, org_z.boundaries, lower = rep(-Inf, Max_looks), theta=delta), n=1 ) - target
}

#-------------------------------------------------------------------------------
compy_design <- function(look_size, alpha.star.u, alpha.star.l)
{
  Uppers =  matrix(rep(NA, length(look_size)^2), ncol=length(look_size))
  Lowers =  matrix(rep(NA, length(look_size)^2), ncol=length(look_size))
  for (i in seq(1, Max_looks)){
    Look = c(rep(NA, i))
    Alpha_u = c(rep(NA, i))
    Alpha_l = c(rep(NA, i))
    for (j in seq(1,i)){
      Look[j] = look_size[j]
    }
    if (i == 1){
      Alpha_u[1] = alpha/2
      Alpha_l[1] = alpha/2
    }
    else{
      Alpha_u[i] = alpha/2 #cumulative one as the increments is calculated in the design.
      Alpha_l[i] = alpha/2 #cumulative one
      Alpha_u[1:(i-1)] = alpha.star.u[1:(i-1)]
      Alpha_l[1:(i-1)] = alpha.star.l[1:(i-1)]
    }
    Uppers[i,1:i]= design(Look, Alpha_u, Alpha_l)$upper
    Lowers[i,1:i]= design(Look, Alpha_u, Alpha_l)$lower
  }
  data.frame(Uppers)
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

sprt_alpha_star <-function(look_size, intersect, slope){
  alpha.star = rep(NA, length(look_size))
  for(i in seq(1, length(look_size))){
    alpha.star[i] = 1 - pnorm(intersect*(1/sqrt(look_size[i])) + slope*sqrt(look_size[i])) + exp(-2*intersect*slope)*pnorm(-(intersect*(1/sqrt(look_size[i]))) + slope*sqrt(look_size[i]))
  }
  alpha.star
}
#-------------------------------------------------------------------------------
apx_sprt_uppers <-function(look_size, intersect, slope){
  upper_apx = rep(NA, length(look_size))
  for(i in seq(1, length(look_size))){
    upper_apx[i] = (intersect - 0.583*sqrt(look_size[1])) + look_size[i]*slope
  }
  upper_apx
}
#-------------------------------------------------------------------------------
sprt_CE <-function(look_size, observed_statistic, stopping_look){
  sprt_CE_val = (alpha/2) * exp(-((delta/2)*look_size[stopping_look-1] - observed_statistic[stopping_look-1])*delta)
  return(c(sprt_CE_val, sprt_CE_val))
}
#-------------------------------------------------------------------------------

delta = 0.5
Max_looks = 10
Max_power = 0.9
alpha = 0.05
B = 1 /(alpha/2)
intersect = log(B) / delta
slope = delta / 2

no_simulation = 100000
#-------------------------------------------------------------------------------

look_size = seq(1,Max_looks)/ Max_looks 
t = look_size
alpha.star.u = (alpha/2)*t
#print(alpha.star.u)
alpha.star.l = alpha.star.u

get_design = design(look_size, alpha.star.u, alpha.star.l)
org_z.boundaries = get_design$upper

Max_size = uniroot(MSN.search, Max_power, Max_looks, org_z.boundaries, interval = c(0,100))$root
#print(Max_size)
new_look_size = seq(1, Max_looks)/ Max_looks * Max_size

org_uppers = org_z.boundaries*sqrt(new_look_size)
#print(org_uppers)


sprt_alpha.star.u = sprt_alpha_star(new_look_size, intersect, slope)
sprt_alpha.star.l = 0.000001*seq(1, Max_looks)
sprt_get_design = design(look_size, sprt_alpha.star.u, sprt_alpha.star.l)
sprt_org_z.boundaries = sprt_get_design$upper
sprt_org_uppers = sprt_org_z.boundaries*sqrt(new_look_size)
#print(sprt_org_uppers)

design_output = compy_design(look_size, alpha.star.u, alpha.star.l)	
rem_z.boundaries = data.frame(design_output)

sprt_design_output = compy_design(look_size, sprt_alpha.star.u, sprt_alpha.star.l)	
sprt_rem_z.boundaries = data.frame(sprt_design_output)


remaining_uppers = rep(NA, length(look_size))
sprt_remaining_uppers = rep(NA, length(look_size))
for(i in seq(1,length(new_look_size))){
  remaining_uppers[i] = rem_z.boundaries[i,i]*sqrt(new_look_size[i])
  sprt_remaining_uppers[i] = sprt_rem_z.boundaries[i,i]*sqrt(new_look_size[i])
} 
#print(remaining_uppers)
#print(sprt_remaining_uppers)
#upper_apx = apx_sprt_uppers(new_look_size, intersect, slope)

start.time = Sys.time()
sim_output = run_simulation(do_simulation, no_simulation=no_simulation, look_size=new_look_size, org_uppers=org_uppers, sprt_org_uppers=sprt_org_uppers, remaining_uppers=remaining_uppers, sprt_remaining_uppers=sprt_remaining_uppers, seed=1)
simulation = data.frame(sim_output)
#print(simulation)
end.time = Sys.time()
difftime(end.time,start.time)



#------------------------------------------------------------------------------

org_counter = rep(0, Max_looks)
rem_counter = rep(0, Max_looks)
CE_counter = rep(0, Max_looks)
sprt_org_counter = rep(0, Max_looks)
sprt_rem_counter = rep(0, Max_looks)
sprt_CE_counter = rep(0, Max_looks)
for (i in seq(1, no_simulation)){
  for (j in 1:6){
    if (j == 1){
      for(k in seq(1,Max_looks)){
        if (simulation[i,k] == -1){
          org_counter[k] = org_counter[k] + 1
        } 
      }
    }
    else if (j == 2){
      for(k in seq(Max_looks+1,2*Max_looks)){
        if (simulation[i,k] == -1){
          rem_counter[k-Max_looks] = rem_counter[k-Max_looks] + 1
        } 
      }
    }
    else if(j == 3){
      for(k in seq(2*Max_looks+1,3*Max_looks)){
        if (simulation[i,k] == -1){
          CE_counter[k-(2*Max_looks)] = CE_counter[k-(2*Max_looks)] + 1
        }
      }
    }
    if (j == 4){
      for(k in seq(3*Max_looks+1,4*Max_looks)){
        if (simulation[i,k] == -1){
          sprt_org_counter[k-(3*Max_looks)] = sprt_org_counter[k-(3*Max_looks)] + 1
        } 
      }
    }
    else if (j == 5){
      for(k in seq(4*Max_looks+1,5*Max_looks)){
        if (simulation[i,k] == -1){
          sprt_rem_counter[k-(4*Max_looks)] = sprt_rem_counter[k-(4*Max_looks)] + 1
        } 
      }
    }
    else if (j == 6){
      for(k in seq(5*Max_looks+1,6*Max_looks)){
        if (simulation[i,k] == -1){
          sprt_CE_counter[k-(5*Max_looks)] = sprt_CE_counter[k-(5*Max_looks)] + 1
        }
      }
    }
  }
}
#print(org_counter / no_simulation)
#print(rem_counter / no_simulation)
#print(CE_counter / no_simulation)
org_type.I = rep(0, Max_looks)#org_counter / no_simulation
rem_type.I = rep(0, Max_looks)#rem_counter / no_simulation
CE_type.I = rep(0, Max_looks)#CE_counter / no_simulation
sprt_org_type.I = rep(0, Max_looks)#org_counter / no_simulation
sprt_rem_type.I = rep(0, Max_looks)#rem_counter / no_simulation
sprt_CE_type.I = rep(0, Max_looks)
for (i in 1:Max_looks){
  org_type.I[i] = sum(org_counter[1:i])/ no_simulation
  rem_type.I[i] = sum(rem_counter[1:i])/ no_simulation
  CE_type.I[i] = sum(CE_counter[1:i])/ no_simulation
  sprt_org_type.I[i] = sum(sprt_org_counter[1:i])/ no_simulation
  sprt_rem_type.I[i] = sum(sprt_rem_counter[1:i])/ no_simulation
  sprt_CE_type.I[i] = sum(sprt_CE_counter[1:i])/ no_simulation
}
print(org_type.I)
print(rem_type.I)
print(CE_type.I)
print(sprt_org_type.I)
print(sprt_rem_type.I)
print(sprt_CE_type.I)

#-------------------------------------------------------------------------------

xdata <- c(1:Max_looks)
y1 <- c(org_type.I)
y2 <- c(rem_type.I)
y3 <- c(CE_type.I)
y4 <- c(alpha.star.u)

y5 <- c(sprt_org_type.I)
y6 <- c(sprt_rem_type.I)
y7 <- c(sprt_CE_type.I)
y8 <- c(sprt_alpha.star.u)

plot(xdata[order(xdata)], y1[order(xdata)], type = "o", pch = 1, col = "blue", lty=1, xlab = "Look", ylab = "Type I error spent", ylim=c(0,0.04))
points(xdata[order(xdata)], y2[order(xdata)], pch = 3, col = "red", lty=1)
lines(xdata[order(xdata)], y2[order(xdata)], pch = 3, col = "red", lty=1)
points(xdata[order(xdata)], y3[order(xdata)], pch = 6, col = "black", lty=1)
lines(xdata[order(xdata)], y3[order(xdata)], pch = 6, col = "black", lty=1)
points(xdata[order(xdata)], y4[order(xdata)], pch = 2, col = "darkgreen", lty=1)
lines(xdata[order(xdata)], y4[order(xdata)], pch = 2, col = "darkgreen", lty=1)
points(xdata[order(xdata)], y5[order(xdata)], pch = 1, col = "blue", lty=2)
lines(xdata[order(xdata)], y5[order(xdata)], pch = 1, col = "blue", lty=2)
points(xdata[order(xdata)], y6[order(xdata)], pch = 3, col = "red", lty=2)
lines(xdata[order(xdata)], y6[order(xdata)], pch = 3, col = "red", lty=2)
points(xdata[order(xdata)], y7[order(xdata)], pch = 6, col = "black", lty=2)
lines(xdata[order(xdata)], y7[order(xdata)], pch = 6, col = "black", lty=2)
points(xdata[order(xdata)], y8[order(xdata)], pch = 2, col = "darkgreen", lty=2)
lines(xdata[order(xdata)], y8[order(xdata)], pch = 2, col = "darkgreen", lty=2)






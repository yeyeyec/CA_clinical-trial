# n = 400 (or 473, see below)
# DELTA = 0, 0.05, 0.1, 0.15
# K = 4 looks
# N = 5000 replications

##########################################################
#######                     base               ###########
##########################################################
library(statip)
library(gsDesign)
install.packages('essentials')
library(essentials)

# Read Data

base = read.csv("CSP1019_ASICSP1019_ASI.csv")

# Processing data

## Generate Dummies

### TREAT: 1- TREAT ; 0 - PLACEBO
base$TREAT_1 = rep(2,300)
for (i in 1:300){
  if (base[i,'TREAT']=='SELEGILINE'){
    base[i,'TREAT_1']=1
  }else{
    base[i,'TREAT_1']=0
  }
}

### GENDER: 1- MALE ; 0 - FEMALE
for (i in 1:300){
  if (base[i,'Sex']==2){
    base[i,'Sex']=0
  }
}

## Log-transformed ASI
base$LogASI = log(base$ASI)

## Generate a subdataset to conduct our project
base = base[c('LogASI','TREAT_1','Sex','ADHD','Pre.HAMD')]
colnames(base) = c('LogASI','Treat','Gender','ADHD','Pre.HAMD')

##########################################################

## Regression
lin.fit = lm('LogASI~Treat + Gender + ADHD + Pre.HAMD',data = base)

beta_hat = lin.fit$coefficients
sigma2_hat = sum(lin.fit$residuals^2)/lin.fit$df.residual
sigma_hat = sqrt(sigma2_hat)


##########################################################
#######               Simulate Data            ###########
##########################################################
n = 400
# Adjust sample size to make Pocock achieves power of 0.9 under alternative
# n = 473

## Generate a single dataset
GENERATE <- function(mu0, delta){
  indice = sample(1:300,n, replace = T)
  df = data.frame(matrix(NA, ncol=3, nrow = n))
  for (j in 1:n){
    df[j,] = base[indice[j],c('Gender','ADHD','Pre.HAMD')]
  }
  ti = rbern(n,1/2)
  df = cbind(df,ti)
  colnames(df) = c('Gender','ADHD','Pre.HAMD','Treat')
  epsilon = rnorm(n,0,sigma_hat)
  df$Y = df$Treat*(mu0+delta)+(1-df$Treat)*mu0+beta_hat[3]*df$Gender+beta_hat[4]*df$ADHD+beta_hat[5]*df$Pre.HAMD+epsilon
  df$Y1 = mu0+delta+beta_hat[3]*df$Gender+beta_hat[4]*df$ADHD+beta_hat[5]*df$Pre.HAMD+epsilon
  df$Y0 = mu0+beta_hat[3]*df$Gender+beta_hat[4]*df$ADHD+beta_hat[5]*df$Pre.HAMD+epsilon
  return(df)
}

##########################################################
#######              Numerical Study           ###########
##########################################################
K = 4

# Calculate T_WO
cal_T_W0 <- function(data){
  n1 = sum(data$Treat)
  nt = length(data$Treat)
  n0 = nt - n1
  y1bar = (1/n1)*(sum(data[,'Treat']*data[,'Y1']))
  y0bar = (1/n0)*(sum((1-data[,'Treat'])*data[,'Y0']))
  s12 = 1/(n1-1) * sum(data[,'Treat']*(data[,'Y'] - y1bar)^2)
  s02 = 1/(n0-1) * sum((1-data[,'Treat'])*(data[,'Y'] - y0bar)^2)
  WaldT = (y1bar-y0bar)/sqrt(s12/n1 + s02/n0)
  return(WaldT)
}

# Calculate T_W1
cal_T_W1 <- function(data){
  n1 = sum(data$Treat)
  nt = length(data$Treat)
  n0 = nt - n1
  ## Obtain the covariance matrix for x1,x2,x3,y1 and y0
  Sigma = cov(data[,c('Gender','ADHD','Pre.HAMD','Y1','Y0')])
  
  ## Obtain Sigma_xx, Sigma_xy1 and Sigma_xy0
  Sigma_xx = Sigma[c(1,2,3),c(1,2,3)]
  Sigma_xy1 = Sigma[c(1,2,3),4]
  Sigma_xy0 = Sigma[c(1,2,3),5]
  ## Test whether Sigma_xx is invertible or not.
  if (det(Sigma_xx)==0){
    return('Colinear.')
  }else{
    ## Obtain the estimated value of beta
    betahat = solve(Sigma_xx) %*% (n1/nt * Sigma_xy1 + n0/nt * Sigma_xy0)
    
    ## Obtain the estimated value of mu1 and mu0
    xbar1 = mean(data$Gender)
    xbar2 = mean(data$ADHD)
    xbar3 = mean(data$Pre.HAMD)
    
    mu1hat = 1/n1 * sum(data$Treat*(data$Y - (data$Gender - xbar1)*betahat[1,1]-(data$ADHD - xbar2)*betahat[2,1]-(data$Pre.HAMD - xbar3)*betahat[3,1]))
    mu0hat = 1/n0 * sum((1-data$Treat)*(data$Y - (data$Gender - xbar1)*betahat[1,1]-(data$ADHD - xbar2)*betahat[2,1]-(data$Pre.HAMD - xbar3)*betahat[3,1]))
    
    ## Obtain the estimated variance
    sigma12hat = 1/n1 * sum(data$Treat*(data$Y - mu1hat -(data$Gender - xbar1)*betahat[1,1]-(data$ADHD - xbar2)*betahat[2,1]-(data$Pre.HAMD - xbar3)*betahat[3,1] )^2)
    sigma02hat = 1/n0 * sum((1-data$Treat)*(data$Y - mu0hat -(data$Gender - xbar1)*betahat[1,1]-(data$ADHD - xbar2)*betahat[2,1]-(data$Pre.HAMD - xbar3)*betahat[3,1])^2)
    
    ## Obtain test statistics
    WaldT = (mu1hat - mu0hat)/sqrt(sigma12hat/n1 + sigma02hat/n0)
    return(WaldT)
  }
}

# Calculate T_W2
cal_T_W2 <- function(data){
  n1 = sum(data$Treat)
  nt = length(data$Treat)
  n0 = nt - n1
  ## Obtain the covariance matrix for x1,x2,x3,y1 and y0
  Sigma = cov(data[,c('Gender','ADHD','Pre.HAMD','Y1','Y0')])
  
  ## Obtain Sigma_xx, Sigma_xy1 and Sigma_xy0
  Sigma_xx = Sigma[c(1,2,3),c(1,2,3)]
  Sigma_xy1 = Sigma[c(1,2,3),4]
  Sigma_xy0 = Sigma[c(1,2,3),5]
  
  ## Test whether Sigma_xx is invertible or not.
  if (det(Sigma_xx)==0){
    return('Colinear.')
  }else{
    ## Obtain the estimated value of beta
    inv = solve(Sigma_xx)
    betahat1 = inv %*% Sigma_xy1
    betahat0 = inv %*% Sigma_xy0
    
    ## Obtain the estimated value of mu1 and mu0
    xbar1 = mean(data$Gender)
    xbar2 = mean(data$ADHD)
    xbar3 = mean(data$Pre.HAMD)
    mu1hat = 1/n1 * sum(data$Treat*(data$Y - (data$Gender - xbar1)*betahat1[1,1] - (data$ADHD - xbar2)*betahat1[2,1] - (data$Pre.HAMD - xbar3)*betahat1[3,1] ))
    mu0hat = 1/n0 * sum((1-data$Treat)*(data$Y - (data$Gender - xbar1)*betahat0[1,1] - (data$ADHD - xbar2)*betahat0[2,1] - (data$Pre.HAMD - xbar3)*betahat0[3,1] ))
    
    ## Obtain the estimated variance of tauhat
    sigma12hat = 1/n1 * sum(data$Treat * 
                              (data$Y - 
                                 data$Treat*(mu1hat+(data$Gender - xbar1)*betahat1[1,1]+(data$ADHD - xbar2)*betahat1[2,1]+(data$Pre.HAMD-xbar3)*betahat1[3,1])-
                                 (1-data$Treat)*(mu0hat+(data$Gender - xbar1)*betahat0[1,1]+(data$ADHD - xbar2)*betahat0[2,1]+(data$Pre.HAMD-xbar3)*betahat0[3,1]))^2)
    sigma02hat = 1/n0 * sum((1-data$Treat) * 
                              (data$Y - 
                                 data$Treat*(mu1hat+(data$Gender - xbar1)*betahat1[1,1]+(data$ADHD - xbar2)*betahat1[2,1]+(data$Pre.HAMD-xbar3)*betahat1[3,1])-
                                 (1-data$Treat)*(mu0hat+(data$Gender - xbar1)*betahat0[1,1]+(data$ADHD - xbar2)*betahat0[2,1]+(data$Pre.HAMD-xbar3)*betahat0[3,1]))^2)
    varhat = nt/n1 * sigma12hat + nt/n0 * sigma02hat +as.scalar(t(betahat1-betahat0) %*% inv %*% (betahat1 - betahat0))
    
    WaldT = sqrt(nt)*(mu1hat - mu0hat)/sqrt(varhat)
    return(WaldT)
  }
}

# A function to conduct a single simulation
# If type == 'W0' : Y ~ Ti *mu1 + (1-Ti)* mu0
# If type == 'W1' : Y ~ Ti *mu1 + (1-Ti)* mu0 + (xi-xbar)^T*beta
# If type == 'W2' : Y ~ Ti *mu1 + (1-Ti)* mu0 + Ti*(xi-xbar)^T*beta1+ (1-Ti)(xi-xbar)^T*beta0

SIMULATE <- function(mu0,delta, type, alpha = 0.05, beta = 0.1){
  x = GENERATE(mu0 = mu0,delta = delta)
  ## Calculate Wald T
  if (type == 'W0'){
    WaldT = cal_T_W0(x)
  }else if (type == 'W1'){
    WaldT = cal_T_W1(x)
  }else{
    WaldT = cal_T_W2(x)
  }
  if (WaldT == 'Colinear.'){
    return('Fail.')
  }else{
    ## Make Decision
    ### Fixed-size
    cv = abs(qnorm(1-alpha/2))
    if (abs(WaldT)>cv){
      result3 = 1
    }else{
      result3 = 0
    }
    
    ### GST
    #### Spliting the data
    X1 <- x[c(1:ceiling(n/4)),] 
    X2 <- x[c(1:ceiling(n/2)),] 
    X3 <- x[c(1:ceiling(3*n/4)),] 
    X4 <- x[c(1:n),]
    
    #### Calculate Z(T)
    Z_l = rep(9999,K)
    
    if (type == 'W0'){
      Z_l[1] = cal_T_W0(X1) 
      Z_l[2] = cal_T_W0(X2) 
      Z_l[3] = cal_T_W0(X3) 
      Z_l[4] = cal_T_W0(X4)
    }else if (type == 'W1'){
      Z_l[1] = cal_T_W1(X1) 
      Z_l[2] = cal_T_W1(X2) 
      Z_l[3] = cal_T_W1(X3) 
      Z_l[4] = cal_T_W1(X4)
    }else{
      Z_l[1] = cal_T_W2(X1) 
      Z_l[2] = cal_T_W2(X2) 
      Z_l[3] = cal_T_W2(X3) 
      Z_l[4] = cal_T_W2(X4)
    }
    if (Z_l[1] == 'Colinear.' ||Z_l[2] == 'Colinear.'||Z_l[3] == 'Colinear.'||Z_l[4] == 'Colinear.' ){
      return('Fail.')
    }else{
      #### Pocock
      X.PO <- gsDesign(k = K, test.type = 2, alpha = alpha/2, beta = beta, sfu = "Pocock") 
      PO.u <- X.PO$upper$bound 
      n_l1 <- c(ceiling(n/4),ceiling(n/2),ceiling(3*n/4),n) 
      flag1 = 0 
      for (i in c(1:K)){ 
        if (abs(Z_l[i])>PO.u[i]){
          n_final1 = n_l1[i]
          reject1 = 1
          break 
        }else{
          flag1 = flag1+1
        } 
        if (flag1==4){
          n_final1 = n_l1[4]
          reject1 = 0
        }
      } 
      decision1 = reject1 
      ss1 = n_final1
      
      #### OF
      X.OF <- gsDesign(k = K, test.type = 2, alpha = alpha/2, beta = beta, sfu = "OF") 
      OF.u <- X.OF$upper$bound 
      n_l2 <- c(ceiling(n/4),ceiling(n/2),ceiling(3*n/4),n) 
      flag2 = 0 
      for (i in c(1:K)){ 
        if (abs(Z_l[i])>OF.u[i]){
          n_final2 = n_l2[i]
          reject2 = 1
          break 
        }else{
          flag2 = flag2+1
        } 
        if (flag2==4){
          n_final2 = n_l2[4]
          reject2 = 0
        }
      }
      decision2 = reject2
      ss2 = n_final2
      return(list(fixd = result3,
                  fixn = n,
                  pod = decision1,
                  pon = ss1,
                  ofd = decision2,
                  ofn = ss2))
    }
  }
}
  
  


# Now conduct N=5000 times Simulation for each delta = 0, 0.05, 0.1 and 0.15
N = 5000
RESULT <- function(delta_v,type){
  decision_list.PO = rep(NA,5000)
  decision_list.OF = rep(NA,5000)
  decision_list.FX = rep(NA,5000)
  sample_size_list.PO = rep(0,5000)
  sample_size_list.OF = rep(0,5000)
  sample_size_list.FX = rep(0,5000)
  i = 1
  while (i <5001){ 
    
    y = SIMULATE(mu0 = 0,delta = delta_v,type = type)
    if (length(y) == 1){
      i = i - 1
      next
    }else{
      decision_list.PO[i] = y$pod
      decision_list.OF[i] = y$ofd
      decision_list.FX[i] = y$fixd
      
      sample_size_list.PO[i] = y$pon
      sample_size_list.OF[i] = y$ofn
      sample_size_list.FX[i] = y$fixn
      i = i + 1
    }
    
  
  }
  #Pr(Reject H0) 
  Pr.PO = mean(decision_list.PO) 
  Pr.OF = mean(decision_list.OF) 
  Pr.FX = mean(decision_list.FX)
  
  #EN 
  Esamplesize.PO = mean(sample_size_list.PO) 
  Esamplesize.OF = mean(sample_size_list.OF) 
  Esamplesize.FX = mean(sample_size_list.FX)
  
  #Sd of N 
  Std.dev.PO = sd(sample_size_list.PO) 
  Std.dev.OF = sd(sample_size_list.OF) 
  Std.dev.FX = sd(sample_size_list.FX) 
  return(list(Pr = c(Pr.PO,Pr.OF,Pr.FX),
              EN = c(Esamplesize.PO,Esamplesize.OF,Esamplesize.FX), 
              SD.N = c(Std.dev.PO,Std.dev.OF,Std.dev.FX)))

}

##########################################################
#######                 W0 AND GST             ###########
##########################################################

W0.0 = RESULT(0,'W0')
W0.05 = RESULT(0.05,'W0')
W0.1 = RESULT(0.1,'W0')
W0.15 = RESULT(0.15,'W0')

##########################################################
#######                 W1 AND GST             ###########
##########################################################

W1.0 = RESULT(0,'W1')
W1.05 = RESULT(0.05,'W1')
W1.1 = RESULT(0.1,'W1')
W1.15 = RESULT(0.15,'W1')

##########################################################
#######                 W2 AND GST             ###########
##########################################################

W2.0 = RESULT(0,'W2')
W2.05 = RESULT(0.05,'W2')
W2.1 = RESULT(0.1,'W2')
W2.15 = RESULT(0.15,'W2')
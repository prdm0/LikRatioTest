##########################################################################################
# Monte Carlo simulations for obtaining the null rejection rates of the LR and ILR tests #
# using the LikRatioTest package                                                         #
##########################################################################################

#Loading LikRatioTest and AdequacyModel packages
library(LikRatioTest)
library(AdequacyModel)

#Testing the shape parameter of the gamma distribution (Section 4, Table 1)

#Probability density function of the gamma distribution
pdf_gamma <- function(par, x, var = NULL){ 
  shape <- par[1]
  scale <- par[2]
  
  if (is.list(var)) eval(parse(text = paste(var[[1]], " <- ", unlist(var[[2]]), sep = "")))
  
  dgamma(x, shape , scale)
}

#Use of the mc() function from the LikRatioTest package
n <- 20 #(n = 20, 30, 40, 50, 100)
r <- 2
k <- 1
RNGkind("L'Ecuyer-CMRG")
mc(N = 15000L,
   n = 20,
   sig = 0.1, #(sig = 0.1, 0.05, 0.01)
   f = pdf_gamma,
   q = rgamma,
   kicks = c(1.01,0.04),
   par0 = list(c("shape"), c(1.0)),
   ncores = 4L,
   shape = 1.0,
   rate = 20,
   #Value of the parameter c (Subsection 3.1, eq.(6))
   c = log(1 - (0.9/(0.9 + 0.1*exp(-0.1*n/(r-k+(k/r))))))/log(0.1 - 0.1*exp(-0.1*n/(r-k+(k/r)))),
   bilateral = FALSE)

#Testing the parameter lambda of the exponentiated Weibull distribution (Section 4, Table 2)

#Probability density function of the exponentiated Weibull distribution
pdf_ew <- function(par, x, var = NULL){
  alpha <- par[1]
  beta <- par[2]
  lambda <- par[3]
  
  if (is.list(var)) eval(parse(text = paste(var[[1]], " <- ", unlist(var[[2]]), sep = "")))
  
  lambda*dweibull(x, shape = beta, scale = alpha)*(1 - exp(-(x/alpha)^beta))^(lambda-1)
}

#Random number generator for the exponentiated Weibull distribution
r_ew <- function(n,alpha,beta,lambda){
  u <- runif(n,0,1)
  alpha*(-log(1-u^(1/lambda)))^(1/beta)
}

#Use of the mc() function from the LikRatioTest package
n <- 30 #(n = 30, 40, 50, 60, 120)
r <- 3
k <- 1
RNGkind("L'Ecuyer-CMRG")
mc(N = 15000L,
   n = 30,
   sig = 0.1, #(sig = 0.1, 0.05, 0.01)
   f = pdf_ew,
   q = r_ew,
   kicks = c(89.9,7.9,1.01),
   par0 = list(c("lambda"),c(1)),
   ncores = 4L,
   alpha = 90,
   beta = 8,
   lambda = 1,
   c = log(1 - (0.9/(0.9 + 0.1*exp(-0.1*n/(r-k+(k/r))))))/log(0.1 - 0.1*exp(-0.1*n/(r-k+(k/r)))),
   bilateral = FALSE)

#Testing the parameters beta and lambda of the exponentiated Weibull distribution (Section 4, Table 3)

n <- 30 #(n = 30, 40, 50, 60, 120)
r <- 3
k <- 2
#data <- r_ew(n = 120L, alpha = 90, beta = 9, lambda = 1.3)
RNGkind("L'Ecuyer-CMRG")
mc(N = 15000L,
   n = 30,
   sig = 0.1, #(sig = 0.1, 0.05, 0.01)
   f = pdf_ew,
   q = r_ew,
   kicks = c(89.5,8.9,1.31),
   par0 = list(c("beta", "lambda"),c(9.0, 1.3)),
   ncores = 4L,
   alpha = 90,
   beta = 9,
   lambda = 1.3,
   c = log(1 - (0.9/(0.9 + 0.1*exp(-0.1*n/(r-k+(k/r))))))/log(0.1 - 0.1*exp(-0.1*n/(r-k+(k/r)))),
   bilateral = FALSE)

##########################################################################################
##########################################################################################

#APPLICATION

#Dataset 1 (Bhaumik et al. (2009)) 
#Successive failures of air-conditioning equipment in Boeing 720 aircraft (n = 29), see Table 5
data01 <- c(90, 14, 44, 310, 130, 10, 24, 59, 76, 208, 60, 56, 29, 26, 70, 186,
            20, 118, 44, 101, 61, 79, 25, 23, 208, 49, 84, 156, 62)

#value of the Wn statistic (MLEs: 1.67 and 0.02, see Table 7)
w <- lrt(f = pdf_gamma, data = data01, kicks = c(1.66, 0.021), par0 = list("shape", 1))

#LR p-values
pchisq(w, df = 1, lower.tail = FALSE)

#Chi-squared inf distribution
fdp_chisq_inf <- function(par, x) {
  k <- par[1]
  c <- par[2]
  
  indet <- 1 - pchisq(q = x, df = k)
  
  indet[indet == 0] <- .Machine$double.eps
  
  dchisq(x = x, df = k) * (1 - (1 - pchisq(q = x, df = k)) ^ c + c * pchisq(q = x, df = k) *
                             (indet) ^ (c - 1))
}

#ILR p-values
p_inf <- function(par) {
  k <- par[1]
  c <- par[2]
  integrate(
    f = fdp_chisq_inf,
    lower = w,
    upper = Inf,
    par = c(k, c)
  )$value
}
n <- 29
r <- 2
k <- 1
p_inf(c(k,log(1 - (0.9/(0.9 + 0.1*exp(-0.1*n/(r-k+(k/r))))))/log(0.1 - 0.1*exp(-0.1*n/(r-k+(k/r))))))

##########################################################################################
##########################################################################################
pdf_ew <- function(par, x, var = NULL) {
  alpha <- par[1]
  sigma <- par[2]
  theta <- par[3]

  if (is.list(var))
    eval(parse(text = paste(var[[1]], " <- ", unlist(var[[2]]), sep = "")))

  alpha * theta / sigma * (1 - exp(-(x / sigma) ^ alpha)) ^ (theta - 1) *
    exp(-(x / sigma) ^ alpha) * (x / sigma) ^ (alpha - 1)
}

rew <- function(n, alpha, sigma, theta){
  u <- runif(n, 0, 1)
  sigma * (-log(1 - u ^ (1 / theta))) ^ (1 / alpha)
}

set.seed(1L, kind = "L'Ecuyer-CMRG")

tictoc::tic()
power_test(N = 1000L,
           B = 250L,
           n = 150L,
           f = pdf_ew,
           sig = 0.05,
           q = rew,
           kicks = c(1, 1, 1),
           par0 = list(c("alpha", "theta"), c(1.5, 1.7)),
           ncores = 4L,
           alpha = 1.5,
           sigma = 1.5,
           theta = 1.7
)
tictoc::toc()

#bootstraping(B = 10, f = pdf_ew, sample_true = rew(100, 1, 1, 1))

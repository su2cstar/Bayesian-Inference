#library(splines)

n <- 300
sigma <- 0.5

## data setting
x <- rnorm(n) # a predictor
t <- runif(n) # time stamp
beta <- function(t) 3*exp(-200*(t-0.2)^2)+exp(-50*(t-0.6)^2) # varying coefficient function
set.seed(123)
y <- sapply(1:n, function(time) beta(t[time])*x[time]+rnorm(1,sd = sigma))
K <- 10 # make a basis spline
knots <- seq(0,1,length.out = K+2)[2:K+1] # Specify knots!

## Make an expanded design matrix!
b_mat <- bs(t, degree = 3, knots = knots,intercept = FALSE)
design_matrix <- x * b_mat

# Obtaining coefficients
fit <- lm(y~design_matrix-1)
phi <- fit$coefficients

# Drawing estimated varying plot!
var_coeff <- b_mat%*%phi
plot(t, var_coeff)

# Okay! Simulation circumstance is made!

# Let's do variational inference!


VI_regression <- function(y,X,max_iters){
  # # of coefficients
  N <- dim(X)[1] ; p <- dim(X)[2]
  # For weight precision
  a <- 1e-7 ; b <- 1e-7 
  # For error precision
  c <- 1e-7 ; d <- 1e-7
  # Fixed params
  a_tilde <- rep(a+0.5, p) ; b_tilde <- rep(b,p)
  c_tilde <- c+(N+1)/2  ; d_tilde <- d
  
  # Default values of params
  mu_coeffs <- rep(0,p)
  sigma_coeffs <- diag(p)
  
  for(i in 1:max_iters){
    # required contributions
    expected_coeffs <- mu_coeffs
    double_expected_coeffs <- sigma_coeffs + mu_coeffs%*%t(mu_coeffs)
    diagonal_sigma <- diag(sigma_coeffs)
    expected_alpha <- sapply(1:p, function(i) a_tilde[i]/b_tilde[i])
    log_expected_alpha <- sapply(1:p, function(i) digamma(a_tilde[i])-log(b_tilde[i]))
    expected_tau <- as.numeric(c_tilde / d_tilde)
    log_expected_tau <- digamma(c_tilde)-log(d_tilde)
    
    sigma_coeffs <- solve(diag(expected_alpha)+expected_tau*(t(X)%*%X))
    mu_coeffs <- expected_tau*sigma_coeffs%*%(t(X)%*%y)
    b_tilde <- sapply(1:p, function(i) (diagonal_sigma[i]+mu_coeffs[i]^2)/2 + b)
    d_tilde <- d+0.5*(t(y)%*%y) - t(expected_coeffs)%*%(t(X)%*%y) + 0.5*sum(diag(X%*%double_expected_coeffs%*%t(X)))
  }
  return(list(mu=mu_coeffs, sigma=sigma_coeffs))
}

## Version of Gibbs_sampler
library(rlist)
Gibbs_VI_regression <- function(y,X,max_iter=1000){
  
  N <- dim(X)[1] ; p <- dim(X)[2]
  
  a <- 1e-7 ; b <- 1e-7 ; c <- 1e-7 ; d <-1e-7
  
  alpha <- rep(0,p)
  tau <- 1
  update_phi <- list()
  update_tau <- c()
  for(iter in 1:max_iter){
    # for phi
    sigma_phi <- solve(tau*t(X)%*%X+diag(p)*alpha)
    mu_phi <- tau*sigma_phi%*%t(X)%*%y
    
    phi <- mvtnorm::rmvnorm(1,mean=mu_phi, sigma=sigma_phi)
    
    # for tau
    shape <- c+N/2 ; rate <- d+0.5*t(y-X%*%t(phi))%*%(y-X%*%t(phi))
    tau <- rgamma(1, shape=shape,rate=rate)
    
    # for alpha_k 
    alpha <- sapply(1:p, function(i) rgamma(1, a+0.5, 0.5*phi[i]^2+b))
    
    if(iter %% 5 == 0){
    update_phi <- list.append(update_phi,phi)
    update_tau <- append(update_tau, tau)
    }
  }
  return(list(phi=update_phi,tau=update_tau))
}

# 95% CI
show_VI_CI_plots <- function(timestamp, sim_data=200, mu, sigma){
  gen_data <- mvtnorm::rmvnorm(sim_data, mean=mu, sigma=sigma)
  basis_matrix <- bs(timestamp, degree = 3, knots = knots, intercept = FALSE)
  coeff_sets <- basis_matrix %*% t(gen_data)
  poly_range <- c(timestamp, rev(timestamp))
  var_matrix <- apply(coeff_sets, 1, function(vec) quantile(vec, probs = c(0.025, 0.5, 0.975)))
  poly_coef_UL <- c(var_matrix[1,], rev(var_matrix[3,]))
  plot(NULL,type="l",ylim=c(min(var_matrix),max(var_matrix)),xlim=c(min(timestamp),max(timestamp)), xlab="",ylab="")
  polygon(poly_range, poly_coef_UL,col=gray(0:9/9)[8],border=F)
  lines(timestamp, var_matrix[2,],lty=1)
  #mtext(paste0('cluster',which_clusters,';varying',j),side=2,line=2.3,cex=0.9)
  mtext("t",side=1,line=2.3,cex=0.9)
  
}

show_Gibbs_CI_plots <- function(timestamp, sim_data=200, coeffs_list){
  basis_matrix <- bs(timestamp, degree = 3, knots = knots, intercept = FALSE)
  coeff_sets <- lapply(coeffs_list, function(lst) basis_matrix %*% as.vector(lst))
  poly_range <- c(timestamp, rev(timestamp))
  coeff_sets <- do.call(cbind, coeff_sets)
  var_matrix <- apply(coeff_sets, 1, function(vec) quantile(vec, probs = c(0.025, 0.5, 0.975)))
  poly_coef_UL <- c(var_matrix[1,], rev(var_matrix[3,]))
  plot(NULL,type="l",ylim=c(min(var_matrix),max(var_matrix)),xlim=c(min(timestamp),max(timestamp)), xlab="",ylab="")
  polygon(poly_range, poly_coef_UL,col=gray(0:9/9)[8],border=F)
  lines(timestamp, var_matrix[2,],lty=1)
  #mtext(paste0('cluster',which_clusters,';varying',j),side=2,line=2.3,cex=0.9)
  mtext("t",side=1,line=2.3,cex=0.9)
}


# Drawing plots 
timestamp <- seq(min(t),max(t),length.out = 50)
out_VI <- VI_regression(y,design_matrix, max_iters = 100)

# Verification!
out_gibbs <- Gibbs_VI_regression(y,design_matrix,max_iter = 5000)
hist(sqrt(1/out_gibbs$tau))

graphics.off()
show_VI_CI_plots(timestamp,mu = out$mu, sigma = out$sigma)
lines(timestamp, beta(timestamp), lty=2,col='blue')

show_Gibbs_CI_plots(timestamp, coeffs_list = out_gibbs$phi)
lines(timestamp, beta(timestamp), lty=2, col='blue')

# 



# Draw the signal plot 

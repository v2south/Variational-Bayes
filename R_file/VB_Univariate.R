library(MCMCpack)
library(MASS)
# True Value: mu = 4, sigma^2 = 1
# For Hyper-parameters
# \mu_{\mu} = 0, \sigma_{\mu}^{2} = 10; A = 0.1, B=0.1
mu_mu = 0; sigma2_mu=10
A = 0.1; B = 0.1
n= 200
# Simulate the data 
x <- rnorm(mean= 4, sd=1, n)

######################Gibbs Sampling###################
M = 10000
samples <- matrix(0, nrow = M, ncol = 2)
colnames(samples) <- c("mu","sigma2")
# Initial values 
mu <-0.5
sigma2 <- 0.005

for( i in 1:M)
{
		mu  <- rnorm( mean = (n*mean (x)/sigma2+ mu_mu/sigma2_mu) / (n/sigma2 + 1/sigma2_mu), sd  = sqrt(1 / ( n/sigma2 + 1/ sigma2_mu)), 1)
		
		sigma2 <- rinvgamma(1, A + n/2, B + 0.5* sum( (x - mu)^2))
		
		samples[i, 1] <-  mu
		
		samples[i, 2] <- sigma2
}

print(apply(samples,2,mean))

print(apply(samples,2,sd))


############### Variaiotnal Bayesian Approximation###########
# ORMEROD  2.2.2 normal Random sample example
# Full derivation is written in seperate sheets
# Tolerence is e-10
x <- rnorm(mean= 4, sd=2, n)
tol <- 10^(-10)

log_delta <- 0.5
# Intial Value for B_{q(sigma^2)}
B_qsigma2 <- 0.000001
iter <- 0 
log_prev <-0
post <- matrix(0,12, 5)
colnames(post) <- c("mu","sigma2","B_qsigma2","Post_Sigma2","log(p(x,q))")

A_qsigma2 <- A+ n/2

while(log_delta > tol )
{
		iter <- iter +1  
		
		sigma2_qu <- 1/ ( n*(A+ n/2)/B_qsigma2 + 1 / sigma2_mu )

	    	u_qu <- ( n* mean(x)*(A+ n /2)/ B_qsigma2 + 0) * sigma2_qu

	    	B_qsigma2 <- B + (sum ((x - u_qu)^2) + n* sigma2_qu) / 2 
		
		log_p <- 0.5 -  n/2  * log( pi * 2) + 0.5* log( sigma2_qu /  sigma2_mu) - ((u_qu - mu_mu)^2 + sigma2_qu) / (2* sigma2_mu) + A*log(B) - (A + n /2)*log(B_qsigma2) + log(gamma(A + n/2)) - log(gamma(A))
		
		log_delta <- abs(log_p - log_prev)
		
		log_prev <- log_p
		
		post_sigma2 <- B_qsigma2 / (A_qsigma2 - 1)
		
		post[iter,] <- cbind(u_qu, sigma2_qu, B_qsigma2 ,post_sigma2,log_p)
}

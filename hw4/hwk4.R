main <- function(datafile,NumberOfIterations,clusterSize)
{
  library(snow)
  library(MASS)
  
  #read the data
  data = read.table(datafile, header = FALSE, sep = "", na.strings = "NA", stringsAsFactors = FALSE)

  newtonRaphsonUpdate <- function(data, response, explanatory, betaInit, maxIter = 100, tol = 1e-4) {
    beta <- betaInit  # Initial beta values (beta0, beta1)
    y <- data[, response]
    x <- data[, explanatory]
    
    # Logistic probability function
    logisticFunction <- function(beta, x) {
      p <- 1 / (1 + exp(-(beta[1] + beta[2] * x)))
      return(p)
    }
    
    # Gradient of the log-likelihood
    computeGradient <- function(beta, y, x) {
      p <- logisticFunction(beta, x)
      gradient_beta0 <- sum(y - p)
      gradient_beta1 <- sum((y - p) * x)
      return(c(gradient_beta0, gradient_beta1))
    }
    
    # Hessian matrix of the log-likelihood
    computeHessian <- function(beta, x) {
      p <- logisticFunction(beta, x)
      W <- p * (1 - p)  # Weight for the second derivative
      H11 <- -sum(W)
      H22 <- -sum(W * x^2)
      H12 <- -sum(W * x)
      Hessian <- matrix(c(H11, H12, H12, H22), nrow = 2)
      return(Hessian)
    }
    
    # Newton-Raphson iteration
    for (k in 1:maxIter) {
      grad <- computeGradient(beta, y, x)
      hessian <- computeHessian(beta, x)
      beta_new <- beta - solve(hessian, grad)  # Update beta using the Newton-Raphson formula
      # Check convergence
      if (max(abs(beta_new - beta)) < tol) {
        break
      }
      beta <- beta_new
    }
    return(beta)
  }
  
  ### PROBLEM 1
  getLaplaceApprox <- function(response, explanatory, data, betaMode) {
    # Extract the response and explanatory variables from data
    y <- data[, response]
    x <- data[, explanatory]
    print("x")
    print(x)
    
    # Calculate probabilities using logistic function
    p <- 1 / (1 + exp(- (betaMode[1] + betaMode[2] * x)))
    
    
    # Compute the log-likelihood
    log_likelihood <- sum(y * log(p) + (1 - y) * log(1 - p))
    print("log_likelihood")
    cat(log_likelihood)
    
    # Calculate the Hessian matrix
    pi_derivative <- p * (1 - p)
    H11 <- -sum(pi_derivative)  # Second derivative w.r.t beta0
    H22 <- -sum(x^2 * pi_derivative)  # Second derivative w.r.t beta1
    H12 <- -sum(x * pi_derivative)  # Mixed derivative
    Hessian <- matrix(c(H11, H12, H12, H22), nrow = 2)
    
    # Calculate the determinant of the Hessian matrix
    det_Hessian <- det(Hessian)
    print("det_Hessian")
    cat(det_Hessian)
    
    # Compute the log of the Laplace approximation of the marginal likelihood
    log_laplace_approx <- log(2 * pi) + log_likelihood - 0.5 * log(abs(det_Hessian))
    print("log_laplace_approx")
    cat(log_laplace_approx)
    
    return(log_laplace_approx)
  }
  
  ### PROBLEM 2
  getPosteriorMeans <- function(response, explanatory, data, betaMode, niter) {
    y <- data[, response]
    x <- data[, explanatory]
    
    # Initialize storage for beta samples
    beta_samples <- matrix(NA, nrow = niter, ncol = 2)
    colnames(beta_samples) <- c("beta0", "beta1")
    beta_samples[1, ] <- betaMode  # Start at the mode
    
    # Define logistic probability function
    logisticFunction <- function(beta, x) {
      p <- 1 / (1 + exp(-(beta[1] + beta[2] * x)))
      return(p)
    }
    # Define Hessian matrix
    computeHessian <- function(beta, x) {
      p <- logisticFunction(beta, x)
      W <- p * (1 - p)  # Weight for the second derivative
      H11 <- -sum(W)
      H22 <- -sum(W * x^2)
      H12 <- -sum(W * x)
      Hessian <- matrix(c(H11, H12, H12, H22), nrow = 2)
      return(Hessian)
    }
    
    # Define function to calculate the log-likelihood
    logLikelihood <- function(beta, y, x) {
      p <- logisticFunction(beta, x)
      return(sum(y * log(p) + (1 - y) * log(1 - p)))
    }
    
    # Metropolis-Hastings Algorithm
    for (i in 2:niter) {
      current_beta <- beta_samples[i - 1, ]
      
      # Propose new beta by drawing from the normal distribution around the current beta
      proposed_beta <- mvrnorm(1, mu = current_beta, Sigma = solve(-computeHessian(current_beta, x)))
      
      # Calculate log-likelihoods
      current_log_likelihood <- logLikelihood(current_beta, y, x)
      proposed_log_likelihood <- logLikelihood(proposed_beta, y, x)
      
      # Calculate acceptance probability
      accept_ratio <- exp(proposed_log_likelihood - current_log_likelihood)
      
      # Accept or reject the proposed move
      if (runif(1) < accept_ratio) {
        beta_samples[i, ] <- proposed_beta
      } else {
        beta_samples[i, ] <- current_beta
      }
    }
    
    # Calculate the sample means of beta0 and beta1
    beta_mean <- colMeans(beta_samples)
    return(beta_mean)
  }
  
  bayesLogistic <- function(apredictor, response, data, NumberOfIterations) {
    # Define initial values for the beta parameters (common to start at zero)
    initialBeta <- c(0, 0)
    
    # Find the mode of the posterior distribution using Newton-Raphson method
    betaMode <- newtonRaphsonUpdate(data, response, apredictor, initialBeta, maxIter = 100, tol = 1e-4)
    
    # Calculate the Laplace approximation of the marginal likelihood
    logmarglik <- getLaplaceApprox(response, apredictor, data, betaMode)
    
    # Generate posterior samples and calculate means
    betaPosteriorMeans <- getPosteriorMeans(response, apredictor, data, betaMode, NumberOfIterations)
    
    # Optionally calculate MLEs for comparison (using the mode as an approximation here)
    beta0mle <- betaMode[1]
    beta1mle <- betaMode[2]
    
    # Return a list of results
    return(list(
      apredictor = apredictor,
      logmarglik = logmarglik,
      beta0bayes = betaPosteriorMeans[1],
      beta1bayes = betaPosteriorMeans[2],
      beta0mle = beta0mle,
      beta1mle = beta1mle
    ))
  }
  
  
  #the sample size is 148 (number of rows)
  #the explanatory variables are the first 60 columns for '534binarydata.txt'
  #the last column is the binary response
  response = ncol(data);
  lastPredictor = ncol(data)-1;
  
  #initialize a cluster for parallel computing
  cluster <- makeCluster(clusterSize, type = "SOCK")

  # Load required libraries on each worker
  clusterEvalQ(cluster, library(MASS))
  
  #run the MC3 algorithm from several times
  results = clusterApply(cluster, 1:lastPredictor, bayesLogistic,
                         response,data,NumberOfIterations);
  
  #print out the results
  for(i in 1:lastPredictor)
  {
    cat('Regression of Y on explanatory variable ',results[[i]]$apredictor,
        ' has log marginal likelihood ',results[[i]]$logmarglik,
        ' with beta0 = ',results[[i]]$beta0bayes,' (',results[[i]]$beta0mle,')',
        ' and beta1 = ',results[[i]]$beta1bayes,' (',results[[i]]$beta1mle,')',
        '\n');    
  }
  
  #destroy the cluster
  stopCluster(cluster);  
}

#this is where the program starts
main("534binarydata.txt", 10000, 10)

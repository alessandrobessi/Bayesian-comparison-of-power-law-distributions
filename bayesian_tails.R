bayesian.tails <- function(x, y, q = 0.75, nsim = 50000, burnin = 5000) {
  require(poweRlaw)
  options(digits = 3)
  x <- x[x > 0]; y <- y[y > 0]
  ox <- displ(x); oy <- displ(y)
  xmin <- quantile(x, q); ymin <- quantile(y, q)
  message("xmin = ", xmin, ">>", q*100  ,"%\nymin = ",
      ymin, ">>", q*100, "%")
  ox$setXmin(xmin); oy$setXmin(ymin)
  alpha.mle.x <- estimate_pars(ox)$pars
  alpha.mle.y <- estimate_pars(oy)$pars

  likelihood <- function(distr, alpha, minim) {
    ll <- dpldis(distr, minim, alpha, log = T)
    return(sum(ll))
  }

  prior <- function(alpha, alpha.mle){
    alpha.prior <- dnorm(alpha, mean = alpha.mle, sd = 0.25, log = T)
    return(alpha.prior)
  }

  posterior <- function(distr, alpha, minim, alpha.mle){
    return(likelihood(distr, alpha, minim) + prior(alpha, alpha.mle))
  }

  proposal.function <- function(param) {
    return(rnorm(1, mean = param, sd = 0.1))
  }

  metropolis <- function(distr, minim, alpha.mle, startvalue, iterations) {
    chain <- array(dim = c(iterations + 1 , 1))
    chain[1,] <- startvalue
    for (i in 1:iterations) {
      proposal <- proposal.function(chain[i,])

      probab <- exp(posterior(distr, proposal, minim, alpha.mle) -
                      posterior(distr, chain[i,], minim, alpha.mle))
      if (runif(1) < probab) {
        chain[i + 1,] <- proposal
      } else {
        chain[i + 1,] <- chain[i,]
      }
    }
    return(chain)
  }

  # ALPHA X
  startvalue <- c(alpha.mle.x)
  chain <- metropolis(x, xmin, alpha.mle.x, startvalue, nsim)

  burn.in <- burnin
  acceptance.rate <- 1 - mean(duplicated(chain[-(1:burn.in),]))
  acceptance.rate

  par(mfrow = c(1,2))
  hist(chain[-(1:burn.in),1], nclass = 100, border = 0, col = rgb(0,0,0,.15),
       main = expression(paste("Posterior of ", alpha[x])), xlab = "MLE estimate = red line")
  abline(v = median(chain[-(1:burn.in),1]), col = "skyblue", lwd = 2)
  abline(v = alpha.mle.x, col = "red", lwd = 2)
  abline(v = quantile(chain[-(1:burn.in),1], 0.05), col = "skyblue", lwd = 2, lty = "dotted")
  abline(v = quantile(chain[-(1:burn.in),1], 0.95), col = "skyblue", lwd = 2, lty = "dotted")

  plot(chain[-(1:burn.in),1], type = "l",
       xlab = "MLE estimate = red line", ylab = expression(alpha[x]),
       main = "MCMC Trace")
  abline(h = alpha.mle.x, col = "red", lwd = 2)

  message("\nmedian alpha x = ", median(chain[-(1:burn.in),1]),
  "\nupper alpha x = ", quantile(chain[-(1:burn.in),1], 0.95),
  "\nlower alpha x = ", quantile(chain[-(1:burn.in),1], 0.05),
  "\nHDI 90% alpha x = ", quantile(chain[-(1:burn.in),1], 0.95) - quantile(chain[-(1:burn.in),1], 0.05))

  posterior.x <- chain[-(1:burn.in),1]

  # ALPHA Y
  startvalue <- c(alpha.mle.y)
  chain <- metropolis(y, ymin, alpha.mle.y, startvalue, nsim)

  burn.in <- burnin
  acceptance.rate <- 1 - mean(duplicated(chain[-(1:burn.in),]))
  acceptance.rate

  par(mfrow = c(1,2))
  hist(chain[-(1:burn.in),1], nclass = 100, border = 0, col = rgb(0,0,0,.15),
       main = expression(paste("Posterior of ", alpha[y])), xlab = "MLE estimate = red line")
  abline(v = median(chain[-(1:burn.in),1]), col = "skyblue", lwd = 2)
  abline(v = alpha.mle.y, col = "red", lwd = 2)
  abline(v = quantile(chain[-(1:burn.in),1], 0.05), col = "skyblue", lwd = 2, lty = "dotted")
  abline(v = quantile(chain[-(1:burn.in),1], 0.95), col = "skyblue", lwd = 2, lty = "dotted")

  plot(chain[-(1:burn.in),1], type = "l",
       xlab = "MLE estimate = red line", ylab = expression(alpha[y]),
       main = "MCMC Trace")
  abline(h = alpha.mle.y, col = "red", lwd = 2)

  message("\nmedian alpha y = ", median(chain[-(1:burn.in),1]),
      "\nupper alpha y = ", quantile(chain[-(1:burn.in),1], 0.95),
      "\nlower alpha y = ", quantile(chain[-(1:burn.in),1], 0.05),
      "\nHDI 90% alpha y = ", quantile(chain[-(1:burn.in),1], 0.95) - quantile(chain[-(1:burn.in),1], 0.05))

  posterior.y <- chain[-(1:burn.in),1]

  # POSTERIOR X - Y
  par(mfrow = c(1,1), mgp = c(3,1,0))
  posterior.diff <- posterior.x - posterior.y
  hist(posterior.diff, nclass = 100, border = 0, col = rgb(0,0,0,.15),
       xlab = expression(alpha[x] - alpha[y]),
       main = expression(paste("Posterior distribution of ", alpha[x] - alpha[y])))
  abline(v = median(posterior.diff), col = "skyblue", lwd = 2)
  abline(v = 0, col = "red", lwd = 2)
  abline(v = quantile(posterior.diff, 0.05), col = "skyblue", lwd = 2, lty = "dotted")
  abline(v = quantile(posterior.diff, 0.95), col = "skyblue", lwd = 2, lty = "dotted")

  message("\nmedian diff = ", median(posterior.diff),
      "\nupper diff = ", quantile(posterior.diff, 0.95),
      "\nlower diff = ", quantile(posterior.diff, 0.05),
      "\nHDI 90% diff = ", quantile(posterior.diff, 0.95) - quantile(posterior.diff, 0.05))
}

library(poweRlaw)
x <- rpldis(1e4, 10, 2)
y <- rpldis(1e4, 10, 2.1)

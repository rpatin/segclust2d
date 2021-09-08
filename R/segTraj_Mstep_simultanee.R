  # Mstep_simultanee
  #' Mstep_simultanee computes the MLE within the EM framework

  #' @param x the bivariate signal
  #' @param rupt the rupture dataframe
  #' @param phi the parameters of the mixture
  #' @param tau the K*P matrix containing posterior probabilities of membership
  #'   to clusters
  #' @param sameSigma TRUE if all segment have the same variance
  #' @return phi the updated value of the parameters

  Mstep_simultanee <- function(x, rupt, tau, phi, sameSigma = TRUE) {
    K <- nrow(tau)
    P <- ncol(tau)
    m <- matrix(nrow = 2, ncol = P)
    s <- matrix(nrow = 2, ncol = P)
    prop <- matrix(nrow = 1, ncol = P)
    Yk <- apply(rupt, 1, FUN = function(y) rowSums(x[, y[1]:y[2]]))
    nk <- rupt[, 2] - rupt[, 1] + 1
    n <- sum(nk)

    #
    np <- nk %*% tau
    m <- Yk %*% tau / rep(np, each = 2)
    if (!sameSigma) {
      for (i in 1:2) {
        s[i, ] <- colSums(
          tau * (vapply(seq_len(P), function(p) {
            apply(rupt, 1, FUN = function(y) sum((x[i, y[1]:y[2]] - m[i, p])^2))
          }))
        )
      }
      s <- sqrt(s / rep(np, each = 2))
    } else {
      for (i in 1:2) {
        s[i, ] <- rep(
          sum(
            tau * (vapply(1:P, function(p) {
              apply(rupt, 1, FUN = function(y) {
                sum((x[i, y[1]:y[2]] - m[i, p])^2)
              })
            }))
          ),
          P
        )
      }
      s <- sqrt(s / n)
    }

    # prop = apply(tau,2,sum)/K
    # emptyCluster = which(prop==0)
    # if(length(emptyCluster)>0){
    #   prop = pmax(prop, eps)
    #   prop = prop /sum(prop)
    #   for (d in emptyCluster){
    #     m[,d]=rep(0,2)
    #     s[,d]=rep(1e9,2)
    #   }
    # }

    prop <- apply(tau, 2, sum) / K
    b <- order(m[1, ])
    m <- m[, b]
    s <- s[, b]
    prop <- prop[b]
    phi <- list(mu = m, sigma = s, prop = prop)

    invisible(phi)
  }

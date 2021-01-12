#' Clustering with EM for Mixed Data
#'
#' @param X A dataframe | matrix
#' @param K A number
#' @param nIter A number
#' @param init A string
#' @return list(Log-likelihood value over all iterations, \code{z} predicted classes)
em_clust <- function(X, K, nIter, init = 'kmeans') {
  N <- nrow(X); J <- ncol(X); theta <- list() # paramètres du modèle

  #Variables quantitatives
  Xqt <- X %>% dplyr::select(where(is.numeric))
  Jqt <- ncol(Xqt)

  #Variables catégorielles
  Xcat <-
    X %>% dplyr::select(where(is.factor) | where(is.character))
  if (!is.character(Xcat))
    Xcat <-
    data.frame(lapply(Xcat, as.factor), stringsAsFactors = FALSE)
  Jcat <- ncol(Xcat)
  ## Encore
  if(Jcat == 0) print('Pas de variables catégorielles')

  # Le nombre de modalités de chaque variable j
  hj <- numeric(Jcat)
  for (j in 1:Jcat)
    hj[j] <- length(as.matrix(unique(Xcat[, j])))

  #Initialisation aléatoire de theta
  if (init == 'random') {
    theta <- RandomInit(Xqt, Xcat, K, N, hj)
  }

  #Initialisation kmeans de theta
  if (init == 'kmeans') {
    z <- kmeans(Xqt, K)$cluster
    theta <- KmeansInit(z, Xqt, Xcat, K, N, hj)
  }

  # Initialisation du loglik de la première itération
  loglik <- rep(0, nIter)
  # tmp <- 0
  # for(i in 1:N)
  #   tmp <- tmp + loglik_fun(Xcat[i, ], Xqt[i, ], K, theta)
  loglik[1] <- Reduce('+', lapply(1:N, function(i) loglik_fun(Xcat[i, ], Xqt[i, ], K, theta)))

  Delta <- 1; iter <- 2; tk <- matrix(0, N, K)

  while (abs(Delta) > 1e-6 && iter <= nIter) {

    # Expectation
    for(i in 1:N)
      tk[i,] <- pkfk_fun(Xcat[i, ], Xqt[i, ], K, theta)
    tk <- tk / rowSums(tk)

    # Maximisation
    for (k in 1:K) {
      nk <- sum(tk[, k])
      theta$mu[k, ] <- colSums(tk[, k] * Xqt) / nk
      theta$pk[k] <- nk / N
      for (j in 1:Jcat) {
        jbin <- vector2binary(Xcat[, j])
        for (h in  1:ncol(jbin))
          theta$alpha[h, j, k] <- colSums(tk[, k] * jbin)[h] / nk
      }
    }

    for (i in 1:Jqt) {
      for (j in 1:Jqt) {
        for (k in 1:K) {
          theta$mat_cov[i, j, k] <- (1 / sum(tk[, k])) * sum(tk[, k] * (Xqt[, i] - theta$mu[k, i]) * (Xqt[, j] - theta$mu[k, j]))
        }
      }
    }


    # Calcul de loglik
    # tmp <- 0
    # for(i in 1:N)
    #   tmp <- tmp + loglik_fun(Xcat[i, ], Xqt[i, ], K, theta)
    loglik[iter] <- Reduce('+', lapply(1:N, function(i) loglik_fun(Xcat[i, ], Xqt[i, ], K, theta)))

    # Mise à jour pour la vérification de la condition
    Delta <- loglik[iter] - loglik[iter - 1]
    iter <- iter + 1
  }

  # Résultats du clustering
  z <- max.col(tk)

  # Nombre de paramètres du modèle
  nu <-  K * (sum(hj - 1) + ((Jqt + 1)/2 + Jqt + 1) - 1) 

  # Critère BIC
  BICrit <- loglik[iter-1] - nu * log(N) /2

  # Critère ICL
  ICL <- BICrit + 2 * Reduce('+', lapply(1:N, function(i) log(tk[i, z[i]])))


  return(list(loglik = loglik,
              theta = theta,
              z = z,
              iter = iter,
              BIC = BICrit,
              ICL = ICL))
}

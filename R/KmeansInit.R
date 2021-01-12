#' Kmeans initialization for clustering with EM
#'
#' @param z A list
#' @param Xqt A matrix
#' @param Xcat A matrix
#' @param K A list
#' @param N A number
#' @param hj A list
#' @return Model parameters : theta = list(mu, mat_cov, pk, alpha)
KmeansInit <- function(z, Xqt, Xcat, K, N, hj) {
  ## Initialisation des paramètres de modèle par kmeans
  #  Xqt : Données quantitaives
  #  Xcat : Données catégorielles
  #  K : Nombre de classes
  #  N : Nombre d'observations
  #  Jcat : Nombre de variables catégorielles
  #  hj : Nombre de modalités pour chaque variable j

  # Déclaration
  Jcat <- ncol(Xcat)
  Jqt <- ncol(Xqt)
  mu <- matrix(0, K, Jqt)
  mat_cov <- array(0, dim = c(Jqt, Jqt, K))
  pk <- numeric(K)
  alpha <- array(0, dim = c(max(hj), Jcat, K))

  # Initialisation Kmeans
  for (iterK in 1:K) {
    mu[iterK, ] <- colMeans(Xqt[which(z == iterK), ])
    mat_cov[, , iterK] <- cov(Xqt[which(z == iterK), ])
    nk <- nrow(Xqt[which(z == iterK), ])
    pk[iterK] <- mean(z == iterK)

    for (j in 1:Jcat) {
      jbin <- vector2binary(Xcat[which(z == iterK), j])
      for (h in  1:ncol(jbin))
        alpha[h, j, iterK] <- colSums(jbin)[h] / nk
    }
  }

  return(theta = list(mu = mu, mat_cov = mat_cov, pk = pk, alpha = alpha))
}

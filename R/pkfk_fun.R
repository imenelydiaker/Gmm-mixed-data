#' Compute proportions for each cluster
#'
#' @param Vcat An array
#' @param Vqt An array
#' @param K An array
#' @param theta A list
#' @return loglik value
pkfk_fun <- function(Vcat, Vqt, K, theta) {
  ## retourne un vecteur de taille K qui calcule pk * fk pour chaque cluster
  # Xqt : vecteur colonne
  # Xcat : vecteur colonne
  # K : nombre de Clusters
  # theta : paramètres du modèle

  Jcat <- ncol(Vcat)
  Jqt <- ncol(Vqt)
  pkfk <- matrix(0, 1, K)
  loglik <- 0
  for (iterK in 1:K) {
    multinom <- 1
    for (j in 1:Jcat) {
      bin_data <- vector2binary(Vcat[,j])
      multinom <- multinom * theta$alpha[which(bin_data == 1), j, iterK]
    }
    pkfk[iterK] <-
      theta$pk[iterK] * multinom * dmvnorm(Vqt, mean = theta$mu[iterK,], sigma = theta$mat_cov[, , iterK])
  }
  return(pkfk)
}

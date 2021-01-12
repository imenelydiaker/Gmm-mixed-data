#' Compute loglikelihood value
#'
#' @param Vcat An array
#' @param Vqt An array
#' @param K An array
#' @param theta A list
#' @return loglik value
loglik_fun <- function(Vcat, Vqt, K, theta) {
  ## retourne la valeur de la log-vraisemblance pour une itération
  # Vqt : vecteur colonne
  # Vcat : vecteur colonne
  # K : nombre de Clusters
  # theta : paramètres du modèle

  Jcat <- ncol(Vcat)
  Jqt <- ncol(Vqt)
  tmp <- 1e-9
  loglik <- 0
  for (iterK in 1:K) {
    multinom <- 1
    for(j in 1:Jcat){
      bin_data <- vector2binary(Vcat[,j])
      multinom <- multinom * theta$alpha[which(bin_data == 1), j, iterK]
    }
    tmp <-
      tmp + theta$pk[iterK] * multinom * dmvnorm(Vqt, mean = theta$mu[iterK, ], sigma = theta$mat_cov[, , iterK])
  }
  loglik <- loglik + log(tmp)

  return(loglik)
}

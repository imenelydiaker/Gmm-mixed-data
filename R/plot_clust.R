#' Plot Clustering Results for Mixed Data
#'
#' @param X A dataframe | matrix
#' @param labels A list
#' @return plot of individuals colored in predicted classes by EM algorithm
plot_clust <- function(X, labels){
  res.famd <- FAMD(X, graph = F)
  plot(res.famd$ind$coord[,1],res.famd$ind$coord[,2],
       main = paste("Clustering K = ", k),
       col = labels, pch = c(17, 15, 18, 19)[labels])
}

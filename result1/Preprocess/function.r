calculate_centroid <- function(data, k) {
  kmeans_model <- kmeans(t(data), centers = k)
  centroids <- kmeans_model$centers
  overall_centroid <- colMeans(centroids)
  return(list(centroids = centroids, overall_centroid = overall_centroid))
}


euclidean_distance <- function(a, b) {
  sqrt(sum((a - b)^2))
}


cosine_similarity <- function(x, y) {
  as.numeric(crossprod(x, y) / sqrt(crossprod(x) * crossprod(y)))
}


calculate_centroid <- function(data, k) {

  kmeans_model <- kmeans(data, centers = k)

  centroids <- kmeans_model$centers
  

  overall_centroid <- colMeans(centroids)
  
  return(list(centroids = centroids, overall_centroid = overall_centroid))
}

mahalanobis_distance <- function(x, y, cov_matrix,...) {
  diff <- x - y
  sqrt(t(diff) %*% solve(cov_matrix,...) %*% diff)
}


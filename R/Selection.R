#' selection()
#' @param genes n by P matrix; P columns of genes with length n (composed of 1s and 0s)
#' @param fit_values list of doubles; list of fitness values, with the order of @param genes
#' @param type string; type could be "tournament" or "rank". tournament is selection by tuornament method, and rank is selection by rank
#' @description This function receives input from the output
#'  of evaluation function, and select the best P pairs of genes as parents,
#'  stored as column indexes for the matrix. The output will be a 2 by P matrix,
#'  with each column as the parents' index in the original matrix
#' @name selection
#' @author Jiayang Nie
#' @export
#'
selection <- function(genes, fit_values, type="tournament") {
  P <- dim(genes)[2]
  if (type == "rank") {
    fit_rank <- (P+1) - rank(fit_values)
    select_prob <- 2 * fit_rank / (P * (P + 1))
    parents <- replicate(P, sample(1:P, 2, prob = select_prob))
    return(parents)
  }
  else if (type == "tournament") {
    ind <- sample(1:P, 6*P, replace = TRUE) # Shuffle
    parent_pool <- matrix(ind, 6, P) # Partition into groups of 3
    parents <- apply(parent_pool, 2, tournament_vectorization, fit_values = fit_values)
    return(parents) # Output parents is 2 by P matrix, with indexes of the parents' position in the gene matrix
  }
}

#' tournament_vectorization()
#' @param x list of 6 indexes; each index can be mapped to a fitness value, and a gene
#' @param fit_values list of doubles; list of fitness values passed in, with the order of @param genes
#' @description This function is a helper function for the apply vectorization in selection function
#' @author Jiayang Nie
#'
tournament_vectorization <- function(x, fit_values) {
  group_a <- x[1:3]
  group_b <- x[4:6]
  parent_a <- group_a[which.min(fit_values[group_a])]
  parent_b <- group_a[which.min(fit_values[group_b])]
  return(c(parent_a, parent_b))
}

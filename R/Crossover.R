#' crossover()
#' @param genes n by P matrix; P columns of genes with length n (composed of 1s and 0s)
#' @param parents 2 by P marix; P columns of indexes
#' @description This function receives input from the output
#'  of selection function, and do crossovers on the pairs of parents
#' @name crossover
#' @author Fernando Machado, Jiayang Nie
#' @export

crossover <- function(genes, parents) {
  offsprings <- apply(parents, 2, crossover_helper, genes = genes)
  return(offsprings)
}

#' crossover_helper()
#' @param genes n by P matrix; P columns of genes with length n (composed of 1s and 0s)
#' @param parent.ind list of 2 integers; indexes of the parents
#' @description A helper function for crossover
#' @author Fernando Machado, Jiayang Nie
#'
crossover_helper <- function(parent.ind, genes) {
  parent_a <- genes[,parent.ind[1]]
  parent_b <- genes[,parent.ind[2]]
  n = length(parent_a)
  #Randomly select the cut-off point.
  cutoff <- sample(2:(length(parent_a)-1), 1)
  offspring <- c(parent_a[1:cutoff], parent_b[(cutoff+1):n])
  return(offspring)
}


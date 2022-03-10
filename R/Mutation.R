#' mutation()
#' @param genes n by P matrix; P columns of genes with length n (composed of 1s and 0s)
#' @param mutation_rate double; indicate the mutation rate
#' @description This function receives input from the output
#'  of crossover function, and do mutation on offsprings. By default,
#'  Mutation_rate is set to 1/C where C is the gene length.
#' @name mutation
#' @author Fernando Machado, Jiayang Nie
#' @export
#'
mutation <- function(genes, mutation_rate=NULL) {
  if (is.null(mutation_rate)) {
    mutation_rate = 1/dim(genes)[1]
  }
  P <- dim(genes)[2]
  mutants <- rbinom(P, 1, mutation_rate)
  mutant.inds <- which(mutants == 1)
  if (length(mutant.inds) > 0) {
    for (ind in mutant.inds) {
      genes[,ind] = random_mutate(genes[,ind])
    }
  }
  return(genes)
}

#' random_mutate()
#' @param gene list of integers; a single gene of genes
#' @description Randomly mutate a gene
#' @author Fernando Machado, Jiayang Nie
#'
#'
random_mutate <- function(gene) {
  allel <- sample(1:length(gene), 1)
  gene[allel] <- abs(1-gene[allel])
  return(gene)
}

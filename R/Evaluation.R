
#' evaluation()
#' @param X m by n matrix of X as covariate
#' @param y list of numbers; the variable to be regressed on
#' @param genes n by P matrix; P columns of genes with length n (composed of 1s and 0s)
#' @param fit.function a user defined function; this function will give the
#' fitness value of any given gene and X, and y variable. By defaut, it is
#' the AIC.
#' @name evaluation
#' @author Jiayang Nie
#' @export
evaluation <- function(X, y, genes, fit.function = AIC, glm.family,
                       HallOfFame.genes = NULL, HallOfFame.fit = NULL) {
  n <- dim(genes)[1]
  P <- dim(genes)[2]
  fit.values <- rep(0, P)

  if (is.null(HallOfFame.genes)){
    HallOfFame.genes <- matrix(0, n, P)
    HallOfFame.fit <- rep(Inf, P)
  }
  for (i in 1:P) {
    gene = genes[,i]
    # If all allels in the gene is zero, fit value should be Inf
    if (sum(gene) == 0) {
      fit.value = Inf
    }
    else {
      # Select the columns needed to the regression based on genes
      X.i = X[,which(gene==1)]
      model = glm(y~., family = glm.family, data = as.data.frame(X.i))
      fit.value = fit.function(model)
      fit.values[i] = fit.value
    }
    # If the fit value is better than the worst gene in HoF, then remove the
    # worst gene in HoF and add the new gene into it.
    if (fit.value < max(HallOfFame.fit)) {
      max.ind = which.max(HallOfFame.fit)
      HallOfFame.genes[,max.ind] = gene
      HallOfFame.fit[max.ind] = fit.value
    }
  }
  HallOfFame = list(genes = HallOfFame.genes, fit = HallOfFame.fit)
  toReturn <- list(genes = genes, fit.value = fit.values, HallOfFame = HallOfFame)
  return(toReturn)
}

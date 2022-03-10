# load the source code of the functions to be tested
#source("Selection.R")
#source("Crossover.R")
#source("Mutation.R")
#source("Evaluation.R")


#' select()
#' @param X m by n matrix of X as covariate
#' @param y list of numbers; the variable to be regressed on
#' @param fit.function a user-defined function; this function will give the
#' fitness value of any given gene and X, and y variable. By defaut, it is
#' the AIC.
#' @param glm.family String; user-defined regression type, e.g. Gaussian for lm,
#' binomial for logistic
#' @param operator a user-defined function; this function determines how to
#' mutate the genes. By default it is NULL.
#' @param pop.size integer; population size
#' @param select.type String; type for selection, rank or tournament
#' @param mutate.rate double; rate for mutation, default to be 0.01
#' @param min_iter integer; minimum number of iterations
#' @param max_iter integer; maximum number of iterations
#' @param converge.thres double; threshold to determine whether to stop iterating
#' early. If it is 0.01, that means if the ratio of the change of mean of fitness
#' value from last generation and the mean of fitness value from last generation
#' < 0.01, then we finish our iterations early, ignoring the max_iter.
#' @param verbose integer; determine the type of log to be printed.
#' @name select
#' @description Find the best covariates in X such that the regression on y by
#' a specific type has a low fitness value.
#' @author Jiayang Nie
#' @examples
#' library(GA)
#' X = matrix(rnorm(120*100, 10, 5), 120, 100)
#' y = rowSums(X)
#' X_noise = matrix(rnorm(120*100, -10, 1), 120, 100)
#' X = cbind(X, X_noise)
#' test1 <- select(X,y, max_iter = 100, verbose=2)
#' test2 <- select(X,y, max_iter = 100, select.type="rank", verbose=2)
#' test1$HallOfFame$genes[,1]
#' @export
select <- function(X, y, fit.function = AIC, glm.family = "gaussian",operator = NULL,
                   pop.size = 100,select.type = "tournament", mutate.rate = NULL,
                   min_iter = 5, max_iter = 100, converge.thres = 0.01, verbose=0) {
  library(testit)
  # Checking input parameters' validity
  assert("X is a matrix",
         is.matrix(X))
  assert("X has positive dimensions",
         dim(X)[1] > 0, dim(X)[2] > 0)
  assert("y is a vector",
         is.vector(y))
  assert("y has the same amount of rows that X has",
         length(y) == dim(X)[1])
  assert("The fitness function is a function",
         is.function(fit.function))
  glms<- c('binomial', 'gaussian','Gamma', 'inverse.gaussian', 'poisson',
           'quasibinomial', 'quasipoisson')
  assert("The glm family entry is a valid member of the generalize linear models family",
         glm.family %in% glms)
  assert("The operator is a function or null",
         is.function(operator)|is.null(operator))
  assert("The population size is an integer",
         pop.size %% 1 == 0)
  assert("The selection type is a character string saying either 'rank' or 'tournament'",
         is.character(select.type),
         identical(select.type,'rank')|identical(select.type,'tournament'))
  if (!is.null(mutate.rate)) {
    assert("The mutation rate is a double between 0.0 and 1.0",
         is.double(mutate.rate),
         mutate.rate > 0.0,
         mutate.rate < 1.0)
  }
  assert("The minimum and maximum iterations are positive integers",
         min_iter %% 1 == 0,
         max_iter %% 1 == 0,
         min_iter > 0,
         max_iter > 0 )
  assert("The convergence threshold is of type double",
         is.double(converge.thres))
  assert("The verbose parameter is set to either 1 or 2 or has the default of 0.",
         verbose %in% c(0,1,2))

  n = dim(X)[2]
  P = pop.size
  genes = matrix(rbinom(P*n, 1, 0.5), ncol = P, nrow = n)
  # First Generation generated outside of the loop
  eval <- evaluation(X, y, genes, glm.family = glm.family)
  fit.values <- eval$fit.value
  parents <- selection(genes, fit.values, type=select.type)
  offsprings <- crossover(genes, parents)
  offsprings <- mutation(offsprings, mutate.rate)
  if (!is.null(operator)) {
    offsprings <- operator(offsprings)
  }
  iter = 0
  is_converge = FALSE
  prev.fit.mean = Inf
  # Start iterating for more generations
  while (iter < max_iter && !is_converge) {
    eval <- evaluation(X, y, offsprings, glm.family = glm.family, HallOfFame.fit = eval$HallOfFame$fit,
                       HallOfFame.genes = eval$HallOfFame$genes)
    fit.values <- eval$fit.value
    parents <- selection(offsprings, fit.values, type=select.type)
    offsprings <- crossover(offsprings, parents)
    offsprings <- mutation(offsprings, mutate.rate)
    if (!is.null(operator)) {
      offsprings <- operator(offsprings)
    }
    if (verbose==1) {
      cat("current iteration: ", iter, "\n")
      cat("min of current population fit: ", min(fit.values), "\n")
      cat("25 percentile of current population fit: ", quantile(fit.values, 0.25), "\n")
      cat("median of current population fit: ", median(fit.values), "\n")
      cat("75 percentile of current population fit: ", quantile(fit.values, 0.75), "\n")
      cat("max of current population fit: ", max(fit.values), "\n")
      cat("mean of current population fit: ", mean(fit.values, na.rm = TRUE), "\n")
      cat("sd of current population fit: ", sd(fit.values, na.rm = TRUE), "\n")
      cat("---------------------------------------------------------", "\n")
    }
    else if (verbose==2) {
      cat("current iteration: ", iter, "\n")
      cat("mean of current population fit: ", mean(fit.values, na.rm = TRUE), "\n")
    }
    fit.mean <- mean(fit.values)
    if ((iter > min_iter) && (abs((fit.mean - prev.fit.mean)/fit.mean) < converge.thres)) {
      is_converge = TRUE
    }
    prev.fit.mean = fit.mean
    iter = iter + 1
  }
  HallOfFame = list(genes = eval$HallOfFame$genes, fit = eval$HallOfFame$fit)
  lastPopulation = list(genes = offsprings, fit = fit.values)
  return(list(bestModel = HallOfFame$genes[,which.min(HallOfFame$fit)],HallOfFame=HallOfFame, lastPopulation=lastPopulation))
}

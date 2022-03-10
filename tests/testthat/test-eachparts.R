


# Initializing the parameters for testing
P <- sample(5:100, 1)
n <- sample(5:20, 1)
X <- matrix(rnorm(2*n*n), 2*n, n)
y <- rnorm(2*n)
genes <- matrix(rbinom(P*n, 1, 0.5), ncol = P, nrow = n)
eval <- evaluation(X, y, genes, glm.family="gaussian")

# Unit tests

# Generation 0
test_that("Gene dimensions are valid", {
  expect_equal(dim(eval$genes)[1], n)
  expect_equal(dim(eval$genes)[2], P)
  expect_equal(dim(eval$HallOfFame$genes)[1], n)
  expect_equal(dim(eval$HallOfFame$genes)[2], P)
})

test_that("There are the same amount of fitness values as population values", {
  expect_equal(length(eval$fit.value), P)
  expect_equal(length(eval$HallOfFame$fit), P)
})

test_that("There are no NA values in the HallOfFame's genes and fitness values", {
  expect_false(any(is.na(eval$HallOfFame$genes)))
  expect_false(any(is.na(eval$HallOfFame$fit)))
})

test_that("There are no NAN values in the HallOfFame's genes and fitness values", {
  expect_false(any(is.nan(eval$HallOfFame$genes)))
  expect_false(any(is.nan(eval$HallOfFame$fit)))
})

test_that("There are no NULL values in the HallOfFame's genes and fitness values", {
  expect_false(any(is.null(eval$HallOfFame$genes)))
  expect_false(any(is.null(eval$HallOfFame$fit)))
})

parents <- selection(genes, eval$fit.value, type="rank")

test_that("Check that the 'rank' selection dimensions are valid", {
  expect_equal(dim(parents)[1], 2)
  expect_equal(dim(parents)[2], P)
})

parents <- selection(genes, eval$fit.value, type="tournament")

test_that("Check that the 'tournament' selection dimensions are valid", {
  parents <- selection(genes, eval$fit.value, type="tournament")
  expect_equal(dim(parents)[1], 2)
  expect_equal(dim(parents)[2], P)
})


# Generation 1
offsprings <- crossover(genes, parents)

test_that("Check dimensions of first generation", {
  expect_equal(dim(offsprings)[1], n)
  expect_equal(dim(offsprings)[2], P)
})

mutate1 <- mutation(offsprings)

test_that("Check that dimensions stay the same", {
  expect_equal(dim(mutate1)[1], n)
  expect_equal(dim(mutate1)[2], P)
})

test_that("Check that no mutation happens if we set the mutation rate to 0", {
  mutate1 <- mutation(offsprings, 0)
  expect_equal(mutate1, offsprings)
})

eval <- evaluation(X, y, mutate1, HallOfFame.fit = eval$HallOfFame$fit,
                   HallOfFame.genes = eval$HallOfFame$genes, glm.family="gaussian")

test_that("Check dimensions of output after evaluating first generation", {
  expect_equal(dim(eval$genes)[1], n)
  expect_equal(dim(eval$genes)[2], P)
  expect_equal(dim(eval$HallOfFame$genes)[1], n)
  expect_equal(dim(eval$HallOfFame$genes)[2], P)
})


test_that("There are the same amount of fitness values as population values in the first generation", {
  expect_equal(length(eval$fit.value), P)
  expect_equal(length(eval$HallOfFame$fit), P)
})


test_that("There are no NA values in the HallOfFame's genes and fitness values in the first generation", {
  expect_false(any(is.na(eval$HallOfFame$genes)))
  expect_false(any(is.na(eval$HallOfFame$fit)))
})

test_that("There are no NAN values in the HallOfFame's genes and fitness values in the first generation", {
  expect_false(any(is.nan(eval$HallOfFame$genes)))
  expect_false(any(is.nan(eval$HallOfFame$fit)))
})

test_that("There are no NULL values in the HallOfFame's genes and fitness values in the first generation", {
  expect_false(any(is.null(eval$HallOfFame$genes)))
  expect_false(any(is.null(eval$HallOfFame$fit)))
})

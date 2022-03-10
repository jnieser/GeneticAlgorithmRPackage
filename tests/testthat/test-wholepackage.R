# Testing the package in its entirety

X = matrix(rnorm(120*100, 10, 5), 120, 100)
y = rowSums(X)

# Overview of default values
fit.function <- AIC
glm.family = "gaussian"
operator = NULL
pop.size = 100
select.type = "tournament"
mutate.rate = NULL
min_iter = 5
max_iter = 100
converge.thres = 0.01
verbose=0

#Selection with 1 max iteration has a larger minimum fitness value in it's hall
# of fame than with 5 max iterations
test_that("Selection: 1 vs. 5 max iterations, tournament", {
  test1 <- select(X,y, max_iter = 1)
  test2 <- select(X,y, max_iter = 5)
  boo <- min(test1$HallOfFame$fit) > min(test2$HallOfFame$fit)
  expect_true(boo)
  })

test_that("Selection: 1 vs. 5 max iterations, rank", {
  test1 <- select(X, y, max_iter = 1, select.type = 'rank')
  test2 <- select(X, y, max_iter = 5, select.type = 'rank')
  boo <- min(test1$HallOfFame$fit) > min(test2$HallOfFame$fit)
  expect_true(boo)
})

test_that("Selection: 1 vs. 5 max iterations, mutation rate 0.2", {
  test1 <- select(X, y, max_iter = 1, mutate.rate = 0.2)
  test2 <- select(X, y, max_iter = 5, mutate.rate = 0.2)
  boo <- min(test1$HallOfFame$fit) > min(test2$HallOfFame$fit)
  expect_true(boo)
})

test_that("Selection: 1 vs. 5 max iterations, mutation rate 0.7", {
  test1 <- select(X, y, max_iter = 1, mutate.rate = 0.7)
  test2 <- select(X, y, max_iter = 5, mutate.rate = 0.7)
  boo <- min(test1$HallOfFame$fit) > min(test2$HallOfFame$fit)
  expect_true(boo)
})

test_that("Selection: 1 vs. 5 max iterations, convergence threshold 0.2", {
  test1 <- select(X, y, max_iter = 1, converge.thres = 0.2)
  test2 <- select(X, y, max_iter = 5, converge.thres = 0.2)
  boo <- min(test1$HallOfFame$fit) > min(test2$HallOfFame$fit)
  expect_true(boo)
})


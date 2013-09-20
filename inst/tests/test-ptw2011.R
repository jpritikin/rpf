library(testthat)
library(rpf)

context("goodness of fit")

# traditional Pearson X^2 goodness of fit test
pearson.gof <- function(observed, expected) {
  x2 <- sum((observed - expected)^2/expected)
  df <- (dim(observed)[1]-1) * (dim(observed)[2]-1)
  pchisq(x2, df, lower.tail=FALSE)
}

ms <- function(observed, expected, draws) {
  draws * sum((observed - expected)^2)
}

mc.gof.test <- function(observed, expected) {
  reps <- 10000
  cells <- prod(dim(observed))
  n <- sum(observed)
  sim.size <- min(n, 185)  # Perkins, Tygert, Ward (2011, p. 12)
  observed <- observed / n
  expected <- expected / n
  rms.ref <- ms(observed, expected, n)
  
  rms.mc <- rep(0,reps)

  for (h in 1:reps) {
    sim <- sample.int(cells, size=sim.size,
                      prob=c(expected), replace=TRUE)
    got <- rep(0L, cells)
    for (cell in 1:cells) got[cell] <- sum(sim==cell)
    got <- got / sim.size
    got <- matrix(got, nrow=dim(observed)[1], ncol=dim(observed)[2])
    rms.mc[h] <- ms(got, expected, sim.size)
  }
  
  sum(rms.mc >= rms.ref)/reps
}

geolen <- function(v) sqrt(sum(v * v))
angle <- function(v1, v2) acos(sum(v1 * v2) / (geolen(v1) * geolen(v2)))

test_that("mc", {
  bins <- 8L
  v1 <- matrix(runif(bins), nrow=2)
  v2 <- matrix(runif(bins), nrow=2)
  v1 <- v1 / sum(v1)
  v2 <- v2 / sum(v2)
  draws <- 300  # ptw2011 accuracy is proportional to draws
  
  got <- data.frame(gap=seq(.5, .05, -.025))
  for (xx in 1:dim(got)[1]) {
    gap <- got[xx,'gap']
    v3 <- v1 * (1-gap) + v2 * gap
    got$angle[xx] <- angle(v1, v3)
    observed <- v1 * draws
    expected <- v3 * draws
    got$mc[xx] <- mc.gof.test(observed, expected)
    got$ptw[xx] <- ptw2011.gof.test(observed, expected)
    got$x2[xx] <- pearson.gof(observed, expected)
  }

  expect_equal(got$mc, got$ptw, .01)
  expect_true(max(abs(got$x2 - got$mc)) > .25)   # not fair, but illustrative
})

test_that("crazy1", {
  # can obtain less than 0 without check
  observed <- structure(c(0L, 0L, 0L, 0L, 0L, 158L, 35L, 51L, 22L, 40L), .Dim = c(5L,  2L))
  expected <- structure(c(78.8628, 5.0351, 2.835, 0.5783, 0.4605, 108.4319,  23.5408, 35.5775, 17.2235, 33.4546), .Dim = c(5L, 2L))
  expect_true(ptw2011.gof.test(observed, expected) >= 0)
})

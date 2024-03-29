context("objective bayes")
library(breakaway)
set.seed(1)

test_that("objective bayes", {
  z <- rnbinomtable(20, 5, 0.5)
  while (z$Frequency[1] == 20) {
    z <- rnbinomtable(20, 0.5, 0.5)
  }
  expect_is(z, "data.frame")
  
  objective_bayes_geometric(z, 
                            tau = 3,
                            iterations = 20, 
                            burn.in = 2,
                            plot = F,
                            output = F,
                            answers = T)
  
  expect_is(objective_bayes_geometric(z, 
                                      iterations = 10, 
                                      burn.in = 2,
                                      plot = F,
                                      output = F,
                                      answers = T),
            "list")
  expect_is(objective_bayes_negbin(z, 
                                   iterations = 10, 
                                   burn.in = 2,
                                   plot = F,
                                   output = F,
                                   answers = T),
            "list")
  expect_is(objective_bayes_poisson(z, 
                                    iterations = 10, 
                                    burn.in = 2,
                                    plot = F,
                                    output = F,
                                    answers = T),
            "list")
  expect_is(objective_bayes_mixedgeo(z,
                                     iterations = 10,
                                     burn.in = 2,
                                     plot = F,
                                     output = F,
                                     answers = T),
            "list")
  
  z <- rnbinomtable(20, 0.5, 0.5)
  while (z$Frequency[1] == 20) {
    z <- rnbinomtable(20, 0.5, 0.5)
  }
  expect_is(z, "data.frame")
  
  objective_bayes_geometric(z, 
                            tau = 3,
                            iterations = 20, 
                            burn.in = 2,
                            plot = F,
                            output = F,
                            answers = T)
  
  expect_is(objective_bayes_geometric(z, 
                                      iterations = 10, 
                                      burn.in = 2,
                                      plot = F,
                                      output = F,
                                      answers = T),
            "list")
  expect_is(objective_bayes_negbin(z, 
                                   iterations = 10, 
                                   burn.in = 2,
                                   plot = F,
                                   output = F,
                                   answers = T),
            "list")
  expect_is(objective_bayes_poisson(z, 
                                    iterations = 10, 
                                    burn.in = 2,
                                    plot = F,
                                    output = F,
                                    answers = T),
            "list")
  expect_is(objective_bayes_mixedgeo(z,
                                     iterations = 10,
                                     burn.in = 2,
                                     plot = F,
                                     output = F,
                                     answers = T),
            "list")
})

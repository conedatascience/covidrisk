#' Converts Infected persons, Number of Contact and Population
#' to the associated probability of a contact
#'
#' Uses a Taylor approximation for repeated Binomial trials
#' Classically associated with the "Birthday Problem"
#'
#' P(A) = 1- P(A)'
#'
#'@param I integer, the number of persons infected
#'@param n integer, the number of contacts per day
#'@param pop integer, the total population
#'@export


calc_risk <- function(I, n, pop) {
  p_I <- I / pop
  r <- 1 - exp(-p_I*n)
  round(100 * r, 1)
}

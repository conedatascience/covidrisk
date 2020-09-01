#' Create the Number of Likely Infectious Persons
#' Given a contact rate and quarantine period
#'
#'This function fits a stochastic three compartment model for
#'S -> Q (quarantined) -> Recovered with transitions
#'from S->Q governed by a poisson distribution
#'
#' @param n_people integer, the number of people to consider in the Susceptible population
#' @param contact_rate integer, the average number of new cases per day
#' @param infectious_period integer, the average number of days someone is infectious
#' @param n_held_out integer, starting number of people infected
#' @export

run_infections_simulation <- function(n_people,
                                    contact_rate_per_day,
                                    infectious_period,
                                    n_held_out = 0,
                                    positive_test_rate = NULL){

  if(!is.null(positive_test_rate)){
    positivity_correction <- round(positive_test_rate*100)
  } else{
    positivity_correction <- 1
  }

  fit <- odin::odin({
    ## Core Transition Equations

    update(S) <- S - n_SQ
    update(Q) <- Q + n_SQ - n_QR
    update(R) <- R + n_QR

    ## Individual Transitions
    ## Probabilities
    p_QR <- 1- exp(-beta)

    # Poisson
    n_SQ <- rpois(epsilon)

    ## Recovery
    n_QR <- rbinom(Q, p_QR)

    ## Initial states
    initial(S) <- S_ini
    initial(Q) <- Q_ini
    initial(R) <- R_ini

    ## Parameters
    S_ini <- user(100)
    Q_ini <- user(10)
    R_ini <- user(0)

    epsilon <- user(6)
    beta <- user(.07)
  }, verbose = FALSE)

  beta_1 <- 1/infectious_period
  x <- fit(S_ini = n_people,
           beta =beta_1,
           epsilon = contact_rate_per_day,
           Q_ini = n_held_out)


  out <- as.data.frame(replicate(n = 100,
                                 x$run(step = 365),
                                 simplify = F))

  x_res <- data.table::as.data.table(out)

  x_comp <- x_res[, .SD, .SDcols = patterns("Q")]

  val_mu <- rowMeans(x_comp)

  val_range <- matrixStats::rowQuantiles(x = as.matrix(x_comp), probs = c(.05,.95))

  out <- data.frame(average_q=val_mu* positivity_correction,
                    lower_q =val_range[,1]* positivity_correction,
                    upper_q = val_range[,2]* positivity_correction)

  out$day <- 1:nrow(out)

  out

}

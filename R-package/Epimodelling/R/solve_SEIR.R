#' @export

solve_SEIR <- function(N, S0, E0, I0,  mu, beta, a, gamma, delta_t = 0.01, times = 1000){

  #setting up a container
  res <- matrix(0, nrow = times, ncol = 6,
                dimnames = list(NULL, c("N", "S", "E", "I", "R", "time")))

  #setting inital values
  res[1, ] <- c(N, S0, E0, I0, N-S0-E0-I0, delta_t)

  #differential equation for the change in S with respect to delta_t
  dS <- function(N, S, I){
    mu * N - mu * S - (beta * S * I)/N
  }

  #differential equation for the change in E with respect to delta_t
  dE <- function(N, S, I, E){
    (beta * S * I)/N - (mu + a) * E
  }

  #differential equation for the change in S with respect to delta_t
  dI <- function(N, E, I){
    a * E - (gamma + mu) * I
  }

  #solving the system using Euler's method
  for(i in 2:times){
    #saving in the previous state of the system
    S_prev <- res[i-1, 2]
    E_prev <- res[i-1, 3]
    I_prev <- res[i-1, 4]

    #calculating the next step in a small time interval
    S_next <- S_prev + delta_t * dS(N, S_prev, I_prev)
    E_next <- E_prev + delta_t * dE(N, S_prev, I_prev, E_prev)
    I_next <- I_prev + delta_t * dI(N, E_prev, I_prev)
    R_next <- N - S_next - E_next - I_next
    timer <- i * delta_t

    #saving the results
    res[i, ] <- c(N, S_next, E_next, I_next, R_next, timer)
  }

  return(res)

}

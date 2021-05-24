#define a generic method

solve <- function(obj, ...){
  UseMethod("solve", obj)
}

solve.Epimodel <- function(obj, delta_t = 0.01, times = 1000){
  
  #setting up a container
  res <- matrix(0, nrow = times, ncol = length(obj$variables) + 1, 
                dimnames = list(NULL, c(obj$variables , "time")))
  
  #setting inital values
  res[1, ] <- c(obj$var_values, delta_t)
  
  #solving the system using Euler's method
  for(i in 2:times){
    
    #calculating the next step in a small time interval based on the generalised equations 
    for(vars in 1:(length(obj$variables)-2)){
      res[i, vars] <- res[i-1, vars] + delta_t * obj$equations[[vars]](obj$variables, res[i-1, -(length(obj$variables)+1)], obj$parameters, obj$par_values)
    }
    
    #calculating the last variable in the system using the assumption of equal population
    sum_variables <- sum(res[i, 1:(length(obj$variables)-2)])
    res[i, (length(obj$variables)-1)] <- obj$var_values[obj$variables == "N"] - sum_variables
    res[i, length(obj$variables)] <- obj$var_values[obj$variables == "N"]
    res[i, (length(obj$variables)+1)] <- i * delta_t #calculating the time-step
  
  }
  
  class(res) <- c("Epimodelled", "matrix") 
  return(res)
  
}

  proba_res <- solve(proba)
  
  
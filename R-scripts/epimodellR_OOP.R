#OOP version of the 

epimodellR <- function(type = c("SIR", "SIS", "SIRD", "SEIR", "SEIS")){
  
  type <- match.arg(type)
  variables <- strsplit(type, "")
  variables <- unlist(variables)

  #SIR model
  if(type == "SIR"){
    
    #parameters of the model
    parameters <- c("beta", "gamma")
    
    #system of differential equations
    dS <- function(variables, var_values, parameters, par_values){
      
      S <- var_values[variables == "S"]
      I <- var_values[variables == "I"]
      N <- var_values[variables == "N"]
      beta <- par_values[parameters == "beta"]
      
      res <- -1 * (beta * I * S) / N
      return(res)
    }
    
    dI <- function(variables, var_values, parameters, par_values){
      
      S <- var_values[variables == "S"]
      I <- var_values[variables == "I"]
      N <- var_values[variables == "N"]
      beta <- par_values[parameters == "beta"]
      gamma <- par_values[parameters == "gamma"]
      
      res <- ((beta * I * S) / N) - (gamma * I)
      return(res)
    }
    
    dR <- function(variables, var_values, parameters, par_values){
      
      I <- var_values[variables == "I"]
      gamma <- par_values[parameters == "gamma"]
      
      res <- gamma * I
      return(res)
    }
    
    l <- list(type = type,
              variables = c(variables, "N"),
              var_values = c(1, 0, 0, 1),
              parameters = parameters,
              par_values = c(0,0),
              equations = list(dS = dS, dI = dI, dR = dR))
    class(l) <- "Epimodel"
  
  } else if(type == "SIS"){ #SIS model
    
    #parameters of the model
    parameters <- c("beta", "gamma")
    
    #system of differential equations
    dS <- function(variables, var_values, parameters, par_values){
      
      S <- var_values[variables == "S"]
      I <- var_values[variables == "I"]
      N <- var_values[variables == "N"]
      beta <- par_values[parameters == "beta"]
      gamma <- par_values[parameters == "gamma"]
      
      res <- (-1 * (beta * I * S) / N) + (gamma * I)
      return(res)
    }
    
    dI <- function(variables, var_values, parameters, par_values){
      
      S <- var_values[variables == "S"]
      I <- var_values[variables == "I"]
      N <- var_values[variables == "N"]
      beta <- par_values[parameters == "beta"]
      gamma <- par_values[parameters == "gamma"]
      
      res <- ((beta * I * S) / N) + (gamma * I)
      return(res)
    }
    
    
    l <- list(type = type,
              variables = c(unique(variables), "N"),
              var_values = c(1, 0, 1),
              parameters = parameters,
              par_values = c(0,0),
              equations = list(dS = dS, dI = dI))
    class(l) <- "Epimodel"
  
  } else if(type == "SIRD"){ #SIRD model
    
    #parameters of the model
    parameters <- c("beta", "gamma", "mu")
    
    #system of differential equations
    dS <- function(variables, var_values, parameters, par_values){
      
      S <- var_values[variables == "S"]
      I <- var_values[variables == "I"]
      N <- var_values[variables == "N"]
      beta <- par_values[parameters == "beta"]
      
      res <- (-1 * (beta * I * S) / N)
      return(res)
    }
    
    dI <- function(variables, var_values, parameters, par_values){
      
      S <- var_values[variables == "S"]
      I <- var_values[variables == "I"]
      N <- var_values[variables == "N"]
      beta <- par_values[parameters == "beta"]
      gamma <- par_values[parameters == "gamma"]
      mu <- par_values[parameters == "mu"]
      
      res <- ((beta * I * S) / N) - (gamma * I) - (mu *  I)
      return(res)
    }
    
    dR <- function(variables, var_values, parameters, par_values){
      
      
      I <- var_values[variables == "I"]
      gamma <- par_values[parameters == "gamma"]
      
      res <- gamma * I
      return(res)
    }
    
    dD <- function(variables, var_values, parameters, par_values){
      
      I <- var_values[variables == "I"]
      mu <- par_values[parameters == "mu"]
      
      res <- mu * I
      return(res)
    }
    
    
    l <- list(type = type,
              variables = c(variables, "N"),
              var_values = c(1, 0, 0, 0, 1),
              parameters = parameters,
              par_values = c(0,0,0),
              equations = list(dS = dS, dI = dI, dR = dR, dD = dD))
    class(l) <- "Epimodel"
  
  } else if(type == "SEIR"){
   
    #parameters of the model
    parameters <- c("beta", "gamma", "mu", "a")
    
    #system of differential equations
    dS <- function(variables, var_values, parameters, par_values){
      
      S <- var_values[variables == "S"]
      I <- var_values[variables == "I"]
      N <- var_values[variables == "N"]
      beta <- par_values[parameters == "beta"]
      gamma <- par_values[parameters == "gamma"]
      mu <- par_values[parameters == "mu"]  
      
     res <- mu * N - mu * S - (beta * S * I)/N
     return(res)
    }
    
    dE <- function(variables, var_values, parameters, par_values){
      
      S <- var_values[variables == "S"]
      I <- var_values[variables == "I"]
      N <- var_values[variables == "N"]
      E <- var_values[variables == "E"]
      beta <- par_values[parameters == "beta"]
      mu <- par_values[parameters == "mu"]
      a <- par_values[parameters == "a"]
      
     res <- (beta * S * I)/N - (mu + a) * E
     return(res)
    }
    
    dI <- function(variables, var_values, parameters, par_values){
      
      
      I <- var_values[variables == "I"]
      E <- var_values[variables == "E"]
      gamma <- par_values[parameters == "gamma"]
      mu <- par_values[parameters == "mu"]
      a <- par_values[parameters == "a"]
      
     res <- a * E - (gamma + mu) * I 
     return(res)
    }
    
    dR <- function(variables, var_values, parameters, par_values){
      
      I <- var_values[variables == "I"]
      R <- var_values[variables == "R"]
      gamma <- par_values[parameters == "gamma"]
      mu <- par_values[parameters == "mu"]
      
      res <- gamma * I - mu * R
      return(res)
    }
    
    
    l <- list(type = type,
              variables = c(variables, "N"),
              var_values = c(1, 0, 0, 0, 1),
              parameters = parameters,
              par_values = c(0,0,0,0),
              equations = list(dS = dS, dE = dE, dI = dI, dR = dR))
    class(l) <- "Epimodel"
  
  } else if(type == "SEIS"){
    
    
    #parameters of the model
    parameters <- c("beta", "gamma", "mu", "a", "lambda")
    
    #system of differential equations
    dS <- function(variables, var_values, parameters, par_values){
      
      S <- var_values[variables == "S"]
      I <- var_values[variables == "I"]
      N <- var_values[variables == "N"]
      beta <- par_values[parameters == "beta"]
      mu <- par_values[parameters == "mu"]
      gamma <- par_values[parameters == "gamma"]
      lambda <- par_values[parameters == "lambda"]
      
      res <- lambda - (beta * S * I)/N - mu * S + gamma * I
      return(res)
    }
    
    dE <- function(variables, var_values, parameters, par_values){
      
      S <- var_values[variables == "S"]
      I <- var_values[variables == "I"]
      N <- var_values[variables == "N"]
      E <- var_values[variables == "E"]
      beta <- par_values[parameters == "beta"]
      mu <- par_values[parameters == "mu"]
      a <- par_values[parameters == "a"]
      
      res <- (beta * S * I)/N - (mu + a) * E
      return(res)
    }
    
    dI <- function(variables, var_values, parameters, par_values){
      
      I <- var_values[variables == "I"]
      E <- var_values[variables == "E"]
      gamma <- par_values[parameters == "gamma"]
      mu <- par_values[parameters == "mu"]
      a <- par_values[parameters == "a"]
      
      res <- a * E - (gamma + mu) * I 
      return(res)
    }

    
    l <- list(type = type,
              variables = c(unique(variables), "N"),
              var_values = c(1, 0, 0, 1),
              parameters = parameters,
              par_values = c(0,0,0,0, 0),
              equations = list(dS = dS,dE = dE, dI = dI))
    class(l) <- "Epimodel"
  }
  
  return(l)
  
}

proba <- epimodellR("SIR")
proba$var_values <- c(0.8, 0.1, 0, 1)
proba$par_values <- c(0.8, 0.2)






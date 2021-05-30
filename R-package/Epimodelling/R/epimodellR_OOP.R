#'Set up a Compartmental Epidemiological Model
#'
#'
#'The function \emph{epimodellR()} helps you set up your compartmental epidemiological model
#'
#'@param type Type of the model. You can choose from c("SIR", "SIS", "SIRD", "SEIR", "SEIS")
#'
#'@return The function \emph{epimodellR()} return an \emph{Epimodel} object. First, you need to specify the variable and parameter values in the console!
#'The list contains the \emph{type}, the  \emph{variables}, the \emph{var_values}, the \emph{parameters}, the \emph{par_values} and the system of differential \emph{equations} of the given model.
#'
#'@examples
#'
#'mod <- epimodellR("SEIR")
#'  #parameter values
#'mod$par_values
#'
#'@export

epimodellR <- function(type = c("SIR", "SIS", "SIRD", "SEIR", "SEIS")){

  type <- match.arg(type)
  variables <- strsplit(type, "")
  variables <- unique(unlist(variables))

  #function for userinput
  userinputR <- function(variable, parameters){

    #container
    var_values <- numeric()
    par_values <- numeric()
    #var values
    for(vars in 1:length(variable)){
      var_values[vars] <- readline(prompt = paste("Please insert the", variable[vars], "value!"))

      testR <- function(var_values){
        if(var_values < 0){
          warning("A population variable can't be negative")
        }
      }

      tt_neg <- tryCatch(testR(var_values[vars]),error=function(e) e, warning=function(w) w)

      while(is(tt_neg, "warning")){
        cat("The parameters need to be non-negative!")
        var_values[vars] <- readline(prompt = paste("Please insert the", variable[vars], "value!"))

        tt_neg <- tryCatch(testR(var_values[vars]),error=function(e) e, warning=function(w) w)

      }

      tt <- tryCatch(as.integer(var_values[vars]),error=function(e) e, warning=function(w) w)

      while(is(tt, "warning")){
        cat("The parameters need to be in numerical format!")
        var_values[vars] <- readline(prompt = paste("Please insert the", variable[vars], "value!"))

        tt <- tryCatch(as.integer(var_values[vars]),error=function(e) e, warning=function(w) w)

      }

    }

    var_values <- as.numeric(var_values)

    #normalisation of the variables to sum to 1 (as 1 is the population size)
    if(sum(var_values) != 1){
      var_values <- var_values / sum(var_values)
    }

    #adding the constant population size
    var_values[length(var_values) + 1] <- 1


    for(par in 1:length(parameters)){
      par_values[par] <- readline(prompt = paste("Please insert the", parameters[par], "value!"))

      tt <- tryCatch(as.integer(par_values[par]),error=function(e) e, warning=function(w) w)

      while(is(tt, "warning")){
        cat("The parameters need to be in numerical format!")
        par_values[par] <- readline(prompt = paste("Please insert the", parameters[par], "value!"))

        tt <- tryCatch(as.integer(par_values[par]),error=function(e) e, warning=function(w) w)
      }

      par_values <- as.numeric(par_values)

    }


    #ADD A REASONABLE METRIC HERE FOR THE PARAMETERS
    return(list(var_values, par_values))

  }


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

    ##########################################################################
    #Adding user input functionality
    values <- userinputR(variables, parameters)
    var_values <- values[[1]]
    par_values <- values[[2]]


    l <- list(type = type,
              variables = c(variables, "N"),
              var_values = var_values,
              parameters = parameters,
              par_values = par_values,
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

    #userinput part
    values <- userinputR(variables, parameters)
    var_values <- values[[1]]
    par_values <- values[[2]]

    l <- list(type = type,
              variables = c(unique(variables), "N"),
              var_values = var_values,
              parameters = parameters,
              par_values = par_values,
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

    #userinput part
    values <- userinputR(variables, parameters)
    var_values <- values[[1]]
    par_values <- values[[2]]

    l <- list(type = type,
              variables = c(variables, "N"),
              var_values = var_values,
              parameters = parameters,
              par_values = par_values,
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

    #userinput part
    values <- userinputR(variables, parameters)
    var_values <- values[[1]]
    par_values <- values[[2]]


    l <- list(type = type,
              variables = c(variables, "N"),
              var_values = var_values,
              parameters = parameters,
              par_values = par_values,
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

    #userinput part
    values <- userinputR(variables, parameters)
    var_values <- values[[1]]
    par_values <- values[[2]]


    l <- list(type = type,
              variables = c(unique(variables), "N"),
              var_values = var_values,
              parameters = parameters,
              par_values = par_values,
              equations = list(dS = dS,dE = dE, dI = dI))
    class(l) <- "Epimodel"
  }

  return(l)

}

#The same function as before, but with manual input for var_values and par_values
#' @export
epimodellR_manual <- function(type = c("SIR", "SIS", "SIRD", "SEIR", "SEIS")){

  type <- match.arg(type)
  variables <- strsplit(type, "")
  variables <- unique(unlist(variables))


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
              var_values = numeric(),
              parameters = parameters,
              par_values = numeric(),
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
              var_values = numeric(),
              parameters = parameters,
              par_values = numeric(),
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
              var_values = numeric(),
              parameters = parameters,
              par_values = numeric(),
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
              var_values = numeric(),
              parameters = parameters,
              par_values = numeric(),
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
              var_values = numeric(),
              parameters = parameters,
              par_values = numeric(),
              equations = list(dS = dS,dE = dE, dI = dI))
    class(l) <- "Epimodel"
  }

  return(l)

}

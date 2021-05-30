#'Shiny application for simulating Compartmental Epidemiological Models
#'
#' The function \emph{EpiSimulator()} will start the Shiny application
#'
#' @return First, choose the model of your choice and the number of days you would like to simulate the model (also dt).
#' After choosing the model the slider outputs will update based on the model variables and the numeric input will update based on the model parameters.
#' Last, if you have set the starting values and parameters click the \emph{Solve and Plot!} button. The application will print the resulting graph.
#'
#'
#'
#' @export

EpiSimulator <- function(){
  require(shiny)
  require(tidyverse)

  epimodellR_shiny <- function(type = c("SIR", "SIS", "SIRD", "SEIR", "SEIS")){

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
                equations = list(dS = dS, dI = dI, dR = dR)
      )
      class(l) <- "Epimodel_shiny"

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
      class(l) <- "Epimodel_shiny"

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
      class(l) <- "Epimodel_shiny"

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
      class(l) <- "Epimodel_shiny"

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
      class(l) <- "Epimodel_shiny"
    }

    return(l)

  }

  solve <- function(obj, ...){
    UseMethod("solve", obj)
  }

  solve.Epimodel_shiny <- function(obj, variable_values, parameter_values, delta_t = 0.01, days = 10){


    times <- days / delta_t

    #setting up a container
    res <- matrix(0, nrow = times, ncol = length(obj$variables) + 1,
                  dimnames = list(NULL, c(obj$variables , "time")))

    #setting inital values
    res[1, ] <- c(variable_values, delta_t)

    #solving the system using Euler's method
    for(i in 2:times){

      #calculating the next step in a small time interval based on the generalised equations
      for(vars in 1:(length(obj$variables)-2)){
        res[i, vars] <- res[i-1, vars] + delta_t * obj$equations[[vars]](obj$variables, res[i-1, -(length(obj$variables)+1)], obj$parameters, parameter_values)
      }

      #calculating the last variable in the system using the assumption of equal population
      sum_variables <- sum(res[i, 1:(length(obj$variables)-2)])
      res[i, (length(obj$variables)-1)] <- variable_values[obj$variables == "N"] - sum_variables
      res[i, length(obj$variables)] <- variable_values[obj$variables == "N"]
      res[i, (length(obj$variables)+1)] <- i * delta_t #calculating the time-step

    }

    class(res) <- c("Epimodelled", "matrix")
    return(res)

  }

  plot.Epimodelled <- function(res, main = ""){

    vars <- colnames(res)
    vars <- vars[-(length(vars) : (length(vars) - 1)) ]

    time <- res[, ncol(res)]

    matplot(
      time,
      res[, 1:(ncol(res) - 2)],
      type = 'l',
      axes = FALSE,
      lty = 1, lwd = 2,
      ylab = 'Population %', xlab = 'Days', ylim = c(0, 1), xlim = c(0, tail(time, 1)),
      main = main, cex.main = 1.75, cex.lab = 1.25, font.main = 1,
      xaxs = 'i', yaxs = 'i'
    )

    axis(1, cex.axis = 1.5)
    axis(2, las = 2, cex.axis = 1.5)

    legend("topright",
           legend = vars,
           lty = 1, lwd = 2,
           bty = 'n', cex = 1.5,
           col = 1:length(vars)
    )


  }

  description <- list(beta = "The average number of connections per person per time interval",
                      gamma = "Recovery rate (1 / Days to recover)",
                      mu = "Fatality rate (1 / Number of fatal cases per day)",
                      lambda = "Birth rate (1 / Number of new members in the population per day)",
                      a = "Incubation period (1 / Average incubation in days)",
                      S = "Susceptible",
                      E = "Exposed",
                      I = "Infected",
                      R = "Recovered/Removed",
                      D = "Deceased")

  ui <- fluidPage(
    fluidRow(
      tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 0px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #CCFF66;
               z-index: 105;
             }
          "))),
    conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                     tags$div("Loading...",id="loadmessage")),


    fluidRow( style = "background-color: #F2F2F2;",
              column(3, div(style = "margin: 1.5em;
                             background-color: #F2F2F2;"),
                     p("EpiSimulator",
                       style = "font-weight: bold;
                           font-size: 2.5em;
                           "),
                     selectInput("modelselectR", label = "Please choose the model!", choices = c("SIR", "SIS", "SIRD", "SEIR", "SEIS")),
                     numericInput("Days", label = "How long do you want to simulate? (days)", value = 10),
                     numericInput("dt", label = "dt", value = 0.01),
                     actionButton("click", "Solve and Plot!",
                                  style = "background-color: #4f94d8;
                                   color: white;
                                   font-weight: bold;
                                   border: none;
                                   border-radius: 5em;
                                   margin-bottom: 1.5em;
                                   ")
              ),
              column(4, offset = 1, div(style = "margin: 1.5em;
                                        background-color: #F2F2F2;"),
                     uiOutput("vars")
              ),
              column(4, div(style = "margin: 1.5em;
                             background-color: #F2F2F2;"),
                     uiOutput("params")

              )
    ),
    fluidRow(
      plotOutput("plot")
    )
  )


  server <- function(input, output, session) {

    model <- reactive(input$modelselectR)
    model_object <- reactive(epimodellR_shiny(model()))
    var_names <- reactive(model_object()$variables)
    par_names <- reactive(model_object()$parameters)

    output$vars <- renderUI({
      map(var_names()[-length(var_names())], ~ sliderInput(.x,
                                                           label = paste("The proportion of ", description[[.x]], "in the population"),
                                                           min =  0, max = 1, value = 0, step =  0.1))
    })

    output$params <- renderUI({
      map(par_names(), ~ numericInput(.x, label = description[[.x]],
                                      min =  0, value = 0, step = 0.05))
    })

    var_values <- eventReactive(input$click, {
      c(unlist(map(var_names(), ~ input[[.x]])), 1)
    })


    par_values <- eventReactive(input$click, {
      unlist(map(par_names(), ~ input[[.x]]))
    })


    solved <- eventReactive(input$click, {
      solve(obj = model_object(),
            variable_values = var_values(),
            parameter_values = par_values(), delta_t = as.numeric(input$dt), days = as.numeric(input$Days))
    })

    #output$proba <- renderTable(var_values())


    output$plot <- renderPlot(plot(solved()))



  }



  # Run the application
  shinyApp(ui = ui, server = server)

}



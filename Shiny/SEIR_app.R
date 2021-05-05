################################################################################################
#This is the primitive version of the Shiny app based on the SEIR model
#Author: Bence Gergely
#Version: shiny_SEIR_05_05
#Language: R
#Language version: 4.0.5
#Reference: 
### https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology
### https://fabiandablander.com/r/Nonlinear-Infection.html
###############################################################################################

library(shiny)

ui <- fluidPage(
    sidebarPanel(
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
    
    #input buttons    ,
    numericInput("N", "Population size", value = 1,  min = 0),
    numericInput("S0", "Initial number of Susceptible", value = 1,  min = 0),
    numericInput("E0", "Initial number of Exposed", value = 1,  min = 0),
    numericInput("I0", "Initial number of Infectous", value = 1,  min = 0),
    numericInput("mu", "Vitality rate", value = 1,  min = 0),
    numericInput("beta", "Contact rate", value = 1,  min = 0),
    numericInput("a", "Incubation period", value = 1,  min = 0),
    numericInput("gamma", "Recovery rate", value = 1,  min = 0),
    numericInput("time", "Days", value = 1,  min = 0),
    actionButton("click", "Run!",
                 style = "background-color:#47CC83;
                      color:#000000;
                      border-color:#16400C;
                      border-style:solid;
                      border-width:2px;
                      font-size:18px;"),
    #output placeholders
    plotOutput("plot")
    
    
)


server <- function(input, output, session) {
    
    #setting the values - they are reactive to the click this can be changed
    N <- eventReactive(input$click, {input$N})
    S0 <- eventReactive(input$click, {input$S0})
    E0 <- eventReactive(input$click, {input$E0})
    I0 <- eventReactive(input$click, {input$I0})
    mu <- eventReactive(input$click, {input$mu})
    beta <- eventReactive(input$click, {input$beta})
    a <- eventReactive(input$click, {input$a})
    gamma <- eventReactive(input$click, {input$gamma})
    time <- eventReactive(input$click, {input$time})
    
    
    #change it to source, but otherwise the solver for the system
    solve_SEIR <- function(N, S0, E0, I0,  mu, beta, a, gamma, delta_t = 0.01, times = 100000){
        
        ####### setting up a container ##########
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
            
            res[i, ] <- c(N, S_next, E_next, I_next, R_next, timer) 
        }
        
        return(res)
        
    }
    
    
    #save the results of the solver
    
    res <- eventReactive(input$click, {solve_SEIR(N(), 
                                                  S0 = S0(), E0 = E0(), I0 = I0(),
                                                  mu = mu(), beta = beta(), a = a(), gamma = gamma(), 
                                                  times = time())})
    
    #change it to source something, but otherwise the code of the plot
    plot_SIRS <- function(res, main = '') {
        time <- res[, 6]
        
        matplot(
            time, res[, 2:5]/res[1,1], type = 'l', axes = FALSE, lty = 1, lwd = 2,
            ylab = 'Subpopulations', xlab = 'Days', ylim = c(0, 1),
            main = main, cex.main = 1.75, cex.lab = 1.25,
            font.main = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, tail(time, 1))
        )
        
        axis(1, cex.axis = 1.5)
        axis(2, las = 2, cex.axis = 1.5)
        legend(
            30, 0.95, legend = c('S', 'I', 'R'),
            lty = 1, lwd = 2, bty = 'n', cex = 1.5
        )
    }
    
  
    
    #render the plot based on the results
    output$plot <- renderPlot(plot_SIRS(res()))
    
}



# Run the application 
shinyApp(ui = ui, server = server)

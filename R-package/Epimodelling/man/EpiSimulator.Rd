% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interact_SEIR.R
\name{EpiSimulator}
\alias{EpiSimulator}
\title{Shiny application for simulating Compartmental Epidemiological Models}
\usage{
EpiSimulator()
}
\value{
First, choose the model of your choice and the number of days you would like to simulate the model (also dt).
After choosing the model the slider outputs will update based on the model variables and the numeric input will update based on the model parameters.
Last, if you have set the starting values and parameters click the \emph{Solve and Plot!} button. The application will print the resulting graph.
}
\description{
The function \emph{EpiSimulator()} will start the Shiny application
}

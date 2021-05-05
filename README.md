# Epimodelling-SEIR
R package and Shiny app for SEIR models

## The main aim of this project is twofold. 
First, I am creating a set of tools with which users can better understand the dynamics of an epidemic. These tools include R functions where the user can change the parameter values and inspect the results and a resulting simple Shiny application where these parameter changes can be made dynamically.

Second, the R package will be capable of estimating the parameters of the system when observed data is present.

## SEIR model

![Wikipedia_image](https://upload.wikimedia.org/wikipedia/commons/3/3d/SEIR.PNG)

![Model-Assumption](https://latex.codecogs.com/png.latex?\inline&space;S&space;&plus;&space;E&space;&plus;&space;I&space;&plus;&space;R&space;=&space;N)

![Model-S](https://latex.codecogs.com/png.latex?\inline&space;\frac{dS}{dt}&space;=&space;\mu&space;N&space;-&space;\mu&space;S&space;-&space;\frac{\beta&space;S&space;I}{N})

![Model-E](https://latex.codecogs.com/png.latex?\inline&space;\frac{dE}{dt}&space;=&space;\frac{\beta&space;S&space;I}{N}&space;-&space;(\mu&space;&plus;&space;a)E)

![Model-I](https://latex.codecogs.com/png.latex?\inline&space;\frac{dI}{dt}&space;=&space;a&space;E&space;-&space;(\gamma&space;&plus;&space;\mu)I)

![Model-R](https://latex.codecogs.com/png.latex?\inline&space;\frac{dR}{dt}&space;=&space;\gamma&space;I&space;-&space;\mu&space;R)

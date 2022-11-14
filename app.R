#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
#library(tidyverse)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library (deSolve)
library(lubridate)
library(plotly)

# https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf

# Source for the SEIQR model:
# https://www.medrxiv.org/content/10.1101/2020.04.17.20070086v1.full.pdf



# Define UI for application that draws a histogram
ui <- fluidPage(
  shiny::withMathJax(),
  tags$head(tags$style(type="text/css",
                       HTML('.container-fluid {  max-width: 1000px; padding:25px; margin-left: auto; margin-right: auto; }'))),
  
  # Application title
  titlePanel("COVID SEIR Simulation"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      sliderInput("f",
                  "% of normal interactions:",
                  min = 0,
                  max = 100,
                  value = 100),
      actionButton("go", "Run 30 Days"),
      actionButton("reset", "Reset"),
      h4("Simulation Credits:"),
      h5(tags$a("Christopher Belanger, PhD", href = "https://www.linkedin.com/in/christopherabelanger/")),
      h5(tags$a("Read the blog post.", href="https://cbelanger.netlify.app/post/interactive-covid-19-simulation-the-effects-of-physical-distancing-on-transmission/")),
      h5(tags$a("See the code on GitHub.", href="https://github.com/chris31415926535/covid_seir_shiny")),
      h4("Model Citation:"),
      h6(tags$a("Anderson, S. C., et al (2020). Estimating the impact of COVID-19 control measures using a Bayesian model of physical distancing [Preprint]. Epidemiology.", href="https://doi.org/10.1101/2020.04.17.20070086"))
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      #fluidRow(plotOutput("distPlot"))
      fluidRow(plotly::plotlyOutput("distPlot")),
      fluidRow(
        h2("Instructions"),
        p("This simulation is very simple: you choose the degree of physical distancing as a percentage of interactions (100% is business as usual, 0% is no interactions whatsoever) and click \"Run 30 Days,\" and the simulation shows you how many active COVID-19 infectiosn the model predicts over the next 30 days."),
        h2("About This App"),
        p("This simple Shiny app is an interactive simulation of how physical distancing affects the spread of COVID-19 at the population level. It puts you in the role of choosing how much physical distancing the whole society will do: choose the level of interaction compared to \"normal\" (100 is business as usual, 0 is no one interacts with anyone), then press the button to see what happens over the next 30 days."),
        p("It's based on the epidemiological model developed by Anderson SC et al. (2020) and presented in", tags$a("this paper", href = "https://www.medrxiv.org/content/10.1101/2020.04.17.20070086v1.full.pdf"), "which is used by the Public Health Agency of Canada in its", tags$a("official reports", href = "https://www.canada.ca/content/dam/phac-aspc/documents/services/diseases-maladies/coronavirus-disease-covid-19/epidemiological-economic-research-data/update-covid-19-canada-epidemiology-modelling-20210115-en.pdf"), "(p. 9). The authors have created", tags$a("covidseir, a wonderful R package that you can download and install from GitHub", href = "https://github.com/seananderson/covidseir"), "!"),
        
        h2("How does it work?"),
        p("The model is a variant of", tags$a("the epidemiological SEIR model", href = "https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SEIR_model"), "of how diseases spread. Please see", tags$a("the paper", href="https://www.medrxiv.org/content/10.1101/2020.04.17.20070086v1.full.pdf"), "for details, but at a high level it tracks how people move through six states over time:"),
        
        tags$ul(tags$li("$S$ : Susceptible to the virus;"),
                tags$li("$E_1$ : Exposed to the virus;"),
                tags$li("$E_2$ : Exposed, pre-symptomatic, and able to infect others;"),
                tags$li("$I$ : Symptomatic and able to infect others;"),
                tags$li("$Q$ : quarantined; and"),
                tags$li("$R$ Recovered or deceased.")
        ),
        p("There are also six variables for people practicing physical distancing: $S_d$, $E_1_d$, $E_2_d$, $I_d$, $Q_d$, and $R_d$. People are assumed to be immune once they enter state $R$."),
        p("The model describes how these groups change over time with differential equations. I won't repeat them here (you're welcome), but I *will*", tags$a("link to the paper one more time so you can see them yourselves on pages 6 and 7!", href = "https://www.medrxiv.org/content/10.1101/2020.04.17.20070086v1.full.pdf")),
        p("In brief, people start out $S$usceptible, then get $E_1$xposed, then after a while they are both $E_2$xposed and able to infect others, then they become visibly $I$nfected and symptomatic, then either go straight to $R$ecovery, or into $Q$uarantine and then into $R$ecovery."),
        p("At the same time, people can either start or stop physical distancing when they're in any of these states; the model assumes that there is some flow back and forth from distancing and not-distancing, but that in general more people start distancing than stop."),
        
        h2("Please note!"),
        tags$ul(
          tags$li("Policy changes in physical distancing in this simulation happen right away: we assume society as a whole turns on a dime, responsive to your whim."),
          tags$li("The simulation uses Anderson et al.'s (2020) initial conditions, which represent British Columbia."),
          tags$li("The simulation solves the differential equations  numerically (using the **deSolve** package) using a resolution of 0.1 days, and then keeps results for whole-numbered time-steps for plotting."),
          tags$li("**I'm not an epidemiologist** and am not affiliated with Anderson et al. in any way. I just like differential equations."),
          tags$li("This is for edutainment purposes only** to learn about how physical distancing can affect COVID-19 prevalence on the population level."),
          tags$li("It comes with no warranties, guarantees, or manatees!")
        )    
        
      )
    )
  )
)

# Define server logic for the simulation
server <- function(input, output) {
  
  # function for deSolve::ode(): defines the system of 12 differential equations that make up the model
  covid_seir <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      #non-physical distancing differential equations:
      dS <- -B * (I + E2 + f * (Id + E2d)) * (S/N) - ud * S + ur * Sd
      dE1 <- B * (I + E2 + f * (Id + E2d)) * (S/N) - k1 * E1 - ud * E1 + ur * E1d
      dE2 <- k1 * E1 - k2 * E2 - ud * E2 + ur * E2d
      dI <- k2 * E2 - q * I - I/D - ud * I + ur * Id
      dQ <- q * I - Q/D - ud * Q + ur * Qd
      dR <- I/D + Q/D - ud * R + ur * Rd
      
      # physical distancing differential equations:
      dSd <- -f * B * (I+E2 + f * (Id + E2d)) * Sd/N + ud * S - ur * Sd
      dE1d <- f * B * (I + E2 + f*(Id + E2d)) * Sd/N - k1 * E1d + ud * E1 - ur * E1d
      dE2d <- k1 * E1d - k2 * E2d + ud * E2 - ur * E2d
      dId <- k2 * E2d - q * Id - Id/D + ud * I - ur * Id
      dQd <- q * Id - Qd / D + ud * Q - ur * Qd
      dRd <- Id/D + Qd/D + ud * R - ur * Rd
      
      return (list(c(dS, dE1, dE2, dI, dQ, dR,
                     dSd, dE1d, dE2d, dId, dQd, dRd)))
    })
  }
  
  # parameters for the model: these are taken straight from Anderson et al. (2020) and represent British Columbia, Canada
  parameters <- c(
    B = 0.433,
    N = 5100000,    # population of British Columbia
    D = 5,          # Mean duration of infectious period in days
    k1 = 0.2,       # (time to infectiousness)-1 (E1 to E2)
    k2 = 1,         # (time period of pre-symptomatic transmissibility)-1 (E2 to I)
    q = 0.05,       # quarantine rate
    ud = 0.1,       # Rate of people moving to physical distancing
    ur = 0.02,      # Rate of people returning from physical distancing
    psir = 0.3,     # proportion of anticipated cases on day r that are tested and reported
    f=1            # proportion of regular contacts that remain (1 = no distancing, 0 = all hermits)
  )
  
  # Initial state for the model: these are taken straight from Anderson et al. (2020) and represent British Columbia, Canada
  state = c(
    S = 849999,
    E1 = 0.53,
    E2 = 0.13,
    I = 0.67,
    Q = 0,
    R = 0,
    Sd = 4249993,
    E1d = 2.67,
    E2d = 0.67,
    Id = 3.33,
    Qd = 0,
    Rd = 0
  )
  
  # convert the vector into a tibble that we can work with more easily
  state_init <- state %>%
    enframe() %>%
    pivot_wider(names_from = "name", values_from = "value") %>%
    mutate(time = 0, deaths = 0)
  
  # create a reactive dataframe 
  # hat tip to the RStudio community: https://community.rstudio.com/t/how-to-pass-reactive-data-frame-into-reactivevalues/24388/3
  vals <- reactiveValues(
    state = state_init
  )
  
  # function to use deSolve to evolve the model forward 30 days
  # it's called below when the button is clicked
  evolve_model <- function(state_now, parameters, duration, f){
    # we use all the same parameters as in the source paper, except we let the user choose "f",
    # the proportion of baseline interactions that continue
    parameters["f"] <- f
    
    # get the last day's state: these will be our inputs to the model
    last_day <- slice_tail(state_now, n=1)
    
    # extract the date/time, since we're going to trim it before we put the state into deSolve::ode()
    last_time <- last_day$time 
    
    # get a vector of times for the new values: start at the last time, and go ahead 30 days in 0.1-day increments
    times <- seq(last_time, last_time + duration, by = 0.1)
    
    # convert the last-day tibble to a named vector
    # note that we remove "time" and "deaths"--deSolve::ode() will crash if we include state variables that aren't defined
    # in the equations. I've decided not to include a deaths calculation because I haven't done the research to do it properly
    last_day <- last_day %>%
      select(-time, -deaths) %>%
      pivot_longer(cols = everything(), names_to = "name", values_to = "value") %>%
      deframe()
    
    # get next month's values
    out <- ode(y = last_day, times = times, func = covid_seir, parms = parameters) %>%
      as_tibble() %>%
      mutate(across(everything(), as.numeric))
    
    # add to our running total of values and calculate deaths
    # note! I'm not including the deaths because I haven't done the research to know what's appropriate.
    # it also seems like it would scale out of control if active cases passed a certain threshold, since
    # hospitals would be full and it would be out of control. This needs more than zero minutes of thought
    # to implement appropriately.
    state_now <- bind_rows(state_now, out) %>%
      mutate(deaths = (I + Id + E1 + E1d + E2 + E2d + Q + Qd + R + Rd) * 0.03)
    
    return (state_now)
  }
  
  # function to take the current state of the model and make a nice ggplot output
  plot_model <- function(state_now){
    # plot it
    # note that we take the generic variable "time" and convert it to a date starting in Feb 2020
    # note also that I've plotted everyone with an active case, excluding quarantined and recovered
    state_now %>%
      filter(time == floor(time)) %>%
      mutate(Date = lubridate::date("2020-02-01") + lubridate::days(time),
             `Active Cases` = round(I + Id + E1 + E1d + E2 + E2d, digits = 0)) %>%
      ggplot() + 
      geom_line(aes(x = Date, y = `Active Cases`)) +
      scale_y_continuous(labels = scales::comma) +
      labs(x = "Date",
           y = "Human Beings with Active COVID-19 Cases") +
      theme_minimal()
  }
  
  # when the button is clicked, evolve the model forward 30 days
  observeEvent(input$go, {
    vals$state <- evolve_model(vals$state, parameters, duration = 30, (input$f/100))
  })
  
  # if reset button is clicked, go back to original state
  observeEvent(input$reset, {
    vals$state <- state_init
  })  
  
  output$distPlot <- plotly::renderPlotly({#renderPlot({
    # update the plot whenever vals$state changes
    plot_model(vals$state) %>%
      plotly::ggplotly()
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

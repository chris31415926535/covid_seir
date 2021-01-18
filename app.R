#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library (deSolve)
library(lubridate)

# https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf
# https://www.medrxiv.org/content/10.1101/2020.04.17.20070086v1.full.pdf



# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("COVID SEIR Simulation"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("f",
                        "Percent of normal interactions:",
                        min = 0,
                        max = 100,
                        value = 100),
            actionButton("go", "Run 30 Days"),
            h3("Citation for Model:"),
            h5("Anderson, S. C., Edwards, A. M., Yerlanov, M., Mulberry, N., Stockdale, J. E., Iyaniwura, S. A., Falcao, R. C., Otterstatter, M. C., Irvine, M. A., Janjua, N. Z., Coombs, D., & Colijn, C. (2020). Estimating the impact of COVID-19 control measures using a Bayesian model of physical distancing [Preprint]. Epidemiology. https://doi.org/10.1101/2020.04.17.20070086")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot"),
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
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
    

    state_init <- state %>%
        enframe() %>%
        pivot_wider(names_from = "name", values_from = "value") %>%
        mutate(time = 0, deaths = 0)
    
    # reactive dataframe https://community.rstudio.com/t/how-to-pass-reactive-data-frame-into-reactivevalues/24388/3
    
    vals <- reactiveValues(
        state = state_init
    )
    

    evolve_model <- function(state_now, parameters, duration, f){
        
        parameters["f"] <- f
        
        last_day <- slice_tail(state_now, n=1)
        
        last_time <- last_day$time 
        
        times <- seq(last_time, last_time + duration, by = 0.1)
        
        last_day <- last_day %>%
            select(-time, -deaths) %>%
            pivot_longer(cols = everything(), names_to = "name", values_to = "value") %>%
            deframe()
        
        state <- last_day 
        
        # get next month's values
        out <- ode(y = state, times = times, func = covid_seir, parms = parameters) %>%
            as_tibble() %>%
            mutate(across(everything(), as.numeric))
        
        # add to our running total of values and calculate deaths
        state_now <- bind_rows(state_now, out) %>%
            mutate(deaths = (I + Id + E1 + E1d + E2 + E2d + Q + Qd + R + Rd) * 0.03)
        
        return (state_now)
    }
    
    plot_model <- function(state_now){
        # plot it
        state_now %>%
            filter(time == floor(time)) %>%
            mutate(days = lubridate::date("2020-02-01") + lubridate::days(time)) %>%
            #  select(days, I, Id, death)
            #ggplot(aes(x = time, y = I + Id + E1 + E1d + E2 + E2d + Q + Qd + R + Rd)) + 
            ggplot(aes(x = days, y = I + Id )) + #+ E1 + E1d + E2 + E2d)) + 
            geom_line() +
            scale_y_continuous(labels = scales::comma) +
            labs(x = "Date",
                 y = "Human Beings with Active COVID-19 Cases") +
            theme_minimal()
    }
    
    # when the button is clicked, evolve the model
    observeEvent(input$go, {
        vals$state <- evolve_model(vals$state, parameters, duration = 30, (input$f/100))
        
        
    })
    
    output$distPlot <- renderPlot({
        # # generate bins based on input$bins from ui.R
        # x    <- faithful[, 2]
        # bins <- seq(min(x), max(x), length.out = input$bins + 1)
        # 
        # # draw the histogram with the specified number of bins
        # hist(x, breaks = bins, col = 'darkgray', border = 'white')
        
        plot_model(vals$state)
        
        #state_now %>% select(time, I) %>% plot()
        #vals$state %>% select(time, I) %>% plot()
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

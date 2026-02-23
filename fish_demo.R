# app.R
# Shiny app: Logistic growth with constant harvest + slope field background
# dP/dt = r P (1 - P/K) - H

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(deSolve)

ui <- fluidPage(
  titlePanel("Logistic Growth with Harvest: dP/dt = r P (1 - P/K) - H"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("r", "r (1/time)", min = 0.01, max = 2.0, value = 0.6, step = 0.01),
      sliderInput("K", "K (carrying capacity)", min = 100, max = 20000, value = 10000, step = 100),
      sliderInput("H", "H (harvest per time)", min = 0, max = 4000, value = 1000, step = 10),
      sliderInput("P0", "Starting population P(0)", min = 0, max = 20000, value = 2000, step = 50),
      
      tags$hr(),
      
      sliderInput("tmax", "Time window (0 to tmax)", min = 5, max = 100, value = 30, step = 1),
      sliderInput("dt", "Time step (smaller = smoother)", min = 0.01, max = 0.5, value = 0.05, step = 0.01),
      
      tags$hr(),
      
      checkboxInput("show_eq", "Show equilibria (if they exist)", value = TRUE)
    ),
    
    mainPanel(
      plotOutput("plot", height = "650px")
    )
  )
)

server <- function(input, output, session) {
  
  # Right-hand side (slope function)
  f <- reactive({
    function(t, P) {
      input$r * P * (1 - P / input$K) - input$H
    }
  })
  
  # Trajectory using ODE solver
  traj_df <- reactive({
    times <- seq(0, input$tmax, by = input$dt)
    
    ode_fun <- function(t, state, parms) {
      P <- state[["P"]]
      dP <- input$r * P * (1 - P / input$K) - input$H
      list(c(dP))
    }
    
    out <- deSolve::ode(
      y = c(P = input$P0),
      times = times,
      func = ode_fun,
      parms = NULL,
      method = "rk4"
    )
    
    as.data.frame(out) |>
      mutate(P = pmax(P, 0))  # keep nonnegative for plotting
  })
  
  # Slope field grid (light grey segments)
  field_df <- reactive({
    # Plot window for P
    P_max <- max(input$K * 1.2, input$P0 * 1.2, 50)
    P_max <- min(P_max, 25000)
    
    # Denser grid helps visibility
    t_grid <- seq(0, input$tmax, length.out = 35)
    P_grid <- seq(0, P_max, length.out = 35)
    
    grid <- tidyr::expand_grid(t = t_grid, P = P_grid) |>
      mutate(
        slope = input$r * P * (1 - P / input$K) - input$H,
        
        # Fixed horizontal half-length (in time units)
        dt_seg = input$tmax / 70,
        
        # Raw vertical change from slope
        dP_raw = slope * dt_seg,
        
        # Cap vertical change so segments stay visible and on-plot
        P_cap = 0.06 * P_max,
        dP_seg = pmax(pmin(dP_raw, P_cap), -P_cap),
        
        t0 = t - dt_seg,
        t1 = t + dt_seg,
        P0 = P - dP_seg,
        P1 = P + dP_seg
      ) |>
      # keep segments inside plot window
      filter(t0 >= 0, t1 <= input$tmax, P0 >= 0, P1 >= 0, P0 <= P_max, P1 <= P_max)
    
    grid
  })
  
  # Equilibria for constant harvest: solve rP(1 - P/K) - H = 0
  eq_vals <- reactive({
    if (!isTRUE(input$show_eq)) return(numeric(0))
    
    r <- input$r
    K <- input$K
    H <- input$H
    
    # Handle trivial cases
    if (r <= 0 || K <= 0) return(numeric(0))
    
    # Solve: rP - (r/K)P^2 - H = 0
    # => (r/K)P^2 - rP + H = 0
    a <- r / K
    b <- -r
    c <- H
    
    disc <- b^2 - 4 * a * c
    if (disc < 0) return(numeric(0))
    
    P1 <- (-b - sqrt(disc)) / (2 * a)
    P2 <- (-b + sqrt(disc)) / (2 * a)
    
    # Keep real, nonnegative
    sort(unique(pmax(c(P1, P2), 0)))
  })
  
  output$plot <- renderPlot({
    tr <- traj_df()
    fd <- field_df()
    
    # plot limits
    P_max <- max(fd$P, tr$P, input$K * 1.1, input$P0 * 1.1, 50)
    P_max <- min(P_max, 25000)
    
    p <- ggplot() +
      # slope field in background
      geom_segment(
        data = fd,
        aes(x = t0, y = P0, xend = t1, yend = P1),
        color = "grey85",
        linewidth = 0.7
      ) +
      # trajectory
      geom_line(
        data = tr,
        aes(x = time, y = P),
        linewidth = 1.3
      ) +
      labs(
        x = "t",
        y = "P(t)",
        title = "Trajectory with Slope Field Background",
        subtitle = bquote(frac(dP, dt) == .(input$r) * P * (1 - P / .(input$K)) - .(input$H))
      ) +
      coord_cartesian(xlim = c(0, input$tmax), ylim = c(0, P_max)) +
      theme_minimal(base_size = 14)
    
    # Optional equilibrium lines
    ev <- eq_vals()
    if (length(ev) > 0) {
      p <- p +
        geom_hline(
          data = data.frame(Peq = ev),
          aes(yintercept = Peq),
          linetype = "dashed",
          linewidth = 0.9
        )
    }
    
    p
  })
}

shinyApp(ui, server)
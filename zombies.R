# app.R
# Zombie outbreak (logistic interaction) + optional zombie death
# Run: shiny::runApp()

library(shiny)
library(deSolve)
library(ggplot2)

# ---- ODE solvers ----
simulate_basic <- function(beta, N, Z0, t_end, dt) {
  # dZ/dt = beta * Z * (N - Z)
  model <- function(t, state, parms) {
    Z <- state[["Z"]]
    with(as.list(parms), {
      dZ <- beta * Z * (N - Z)
      list(c(dZ))
    })
  }

  times <- seq(0, t_end, by = dt)
  out <- ode(y = c(Z = Z0), times = times, func = model,
             parms = list(beta = beta, N = N), method = "lsoda")
  as.data.frame(out)
}

simulate_with_death <- function(beta, mu, N, Z0, t_end, dt) {
  # dZ/dt = beta * Z * (N - Z) - mu * Z
  model <- function(t, state, parms) {
    Z <- state[["Z"]]
    with(as.list(parms), {
      dZ <- beta * Z * (N - Z) - mu * Z
      list(c(dZ))
    })
  }

  times <- seq(0, t_end, by = dt)
  out <- ode(y = c(Z = Z0), times = times, func = model,
             parms = list(beta = beta, mu = mu, N = N), method = "lsoda")
  as.data.frame(out)
}

# ---- UI ----
ui <- fluidPage(
  titlePanel("Zombie Differential Equation Simulator"),
  tabsetPanel(
    tabPanel(
      "Basic interaction",
      sidebarLayout(
        sidebarPanel(
          sliderInput("beta1", "Beta (transmission/interaction rate)", min = 0, max = 0.01,
                      value = 0.002, step = 0.000001),
          numericInput("N1", "Starting population (total people, N)", value = 149, min = 2, step = 1),
          numericInput("Z01", "Initial zombies (Z(0))", value = 1, min = 1, step = 1),
          sliderInput("tend1", "Time horizon", min = 10, max = 500, value = 200, step = 10),
          sliderInput("dt1", "Time step (smaller = smoother)", min = 0.001, max = 5, value = 1, step = 0.1),
          helpText("Model: dZ/dt = beta * Z * (N - Z)")
        ),
        mainPanel(
          plotOutput("plot1", height = 450),
          verbatimTextOutput("summary1")
        )
      )
    ),

    tabPanel(
      "Interaction + zombie death",
      sidebarLayout(
        sidebarPanel(
          sliderInput("beta2", "Beta (transmission/interaction rate)", min = 0, max = 0.01,
                      value = 0.002, step = 0.000001),
          sliderInput("mu2", "Mu (zombie death/decay rate)", min = 0, max = 100,
                      value = 0.05, step = 0.01),
          numericInput("N2", "Starting population (total people, N)", value = 149, min = 2, step = 1),
          numericInput("Z02", "Initial zombies (Z(0))", value = 1, min = 1, step = 1),
          sliderInput("tend2", "Time horizon", min = 10, max = 500, value = 200, step = 10),
          sliderInput("dt2", "Time step (smaller = smoother)", min = 0.001, max = 5, value = 1, step = 0.1),
          helpText("Model: dZ/dt = beta * Z * (N - Z) - mu * Z")
        ),
        mainPanel(
          plotOutput("plot2", height = 450),
          verbatimTextOutput("summary2")
        )
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {

  # Keep Z0 <= N
  observeEvent(input$N1, {
    if (input$Z01 > input$N1) updateNumericInput(session, "Z01", value = input$N1)
  })
  observeEvent(input$Z01, {
    if (input$Z01 > input$N1) updateNumericInput(session, "Z01", value = input$N1)
  })
  observeEvent(input$N2, {
    if (input$Z02 > input$N2) updateNumericInput(session, "Z02", value = input$N2)
  })
  observeEvent(input$Z02, {
    if (input$Z02 > input$N2) updateNumericInput(session, "Z02", value = input$N2)
  })

  sim1 <- reactive({
    req(input$beta1, input$N1, input$Z01, input$tend1, input$dt1)
    simulate_basic(beta = input$beta1, N = input$N1, Z0 = input$Z01,
                   t_end = input$tend1, dt = input$dt1)
  })

  sim2 <- reactive({
    req(input$beta2, input$mu2, input$N2, input$Z02, input$tend2, input$dt2)
    simulate_with_death(beta = input$beta2, mu = input$mu2, N = input$N2, Z0 = input$Z02,
                        t_end = input$tend2, dt = input$dt2)
  })

  output$plot1 <- renderPlot({
    df <- sim1()
    ggplot(df, aes(x = time, y = Z)) +
      geom_line(linewidth = 1) +
      labs(x = "Time", y = "Zombies, Z(t)",
           title = "Zombie population over time (basic interaction)") +
      coord_cartesian(ylim = c(0, input$N1))
  })

  output$plot2 <- renderPlot({
    df <- sim2()
    ggplot(df, aes(x = time, y = Z)) +
      geom_line(linewidth = 1) +
      labs(x = "Time", y = "Zombies, Z(t)",
           title = "Zombie population over time (interaction + zombie death)") +
      coord_cartesian(ylim = c(0, input$N2))
  })

  output$summary1 <- renderText({
    df <- sim1()
    Z_end <- df$Z[nrow(df)]
    paste0("Final Z(t_end) = ", round(Z_end, 3), " out of N = ", input$N1)
  })

  output$summary2 <- renderText({
    df <- sim2()
    Z_end <- df$Z[nrow(df)]
    paste0("Final Z(t_end) = ", round(Z_end, 3), " out of N = ", input$N2)
  })
}

shinyApp(ui, server)


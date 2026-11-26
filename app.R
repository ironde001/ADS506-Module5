library(shiny)
library(fpp3)
library(tidyverse)
library(gt)

### Data prep
wine_raw <- read_csv("australianwines.csv", na = c("", "NA", "*"),
                     show_col_types = FALSE)

wine_ts <- wine_raw |>
  mutate(
    Month = yearmonth(my(Month)),
    across(-Month, as.numeric)) |>
  pivot_longer(-Month, names_to = "Varietal", values_to = "Sales") |>
  mutate(Varietal = str_to_title(Varietal)) |>
  as_tsibble(index = Month, key = Varietal) |>
  arrange(Varietal, Month) |>
  group_by_key() |>
  mutate(Sales = approx(seq_along(Sales), Sales, 
                        xout = seq_along(Sales), rule = 2)$y) |>
  ungroup()

varietals <- sort(unique(wine_ts$Varietal))
min_date  <- as.Date(min(wine_ts$Month))
max_date  <- as.Date(max(wine_ts$Month))
model_levels <- c("TSLM", "ETS", "ARIMA", "Ensemble")

### UI
ui <- fluidPage(
  titlePanel("Australian Wine Sales Forecaster"),
  sidebarLayout(
    sidebarPanel(
      selectizeInput("varietal", "Varietal(s)",
                     choices = varietals,
                     selected = varietals[1],
                     multiple = TRUE,
                     options = list(plugins = list("remove_button"))),
      dateRangeInput("dates", "Date range",
                     start = min_date, end = max_date,
                     min = min_date,   max = max_date,
                     format = "M yyyy",
                     startview = "year",
                     separator = " to ")),                                              
    mainPanel(
      tabsetPanel(
        tabPanel("Data", p("Monthly tsibble"), verbatimTextOutput("ts_print")),
        tabPanel("Overview", plotOutput("overview", height = "450px")),
        tabPanel("Seasonality", plotOutput("season_plot", height = "450px"),
                 plotOutput("stl_plot", height = "550px")),
        tabPanel("Models", br(), splitLayout(cellWidths = c("50%", "50%"),
                   numericInput("valid_n", "Validation window (months)",
                                value = 12, min = 6, max = 36, step = 6),
                   numericInput("h", "Forecast horizon (months)",
                                value = 24, min = 1, max = 60)),
                 textOutput("split_text"), br(),
                 tabsetPanel(
                   tabPanel("TSLM", tableOutput("tslm_table")),
                   tabPanel("ETS", tableOutput("ets_table")),
                   tabPanel("ARIMA", tableOutput("arima_table"),
                            helpText(HTML("Differencing is fixed at d = 1, D = 1;
                               p, q, P, and Q are chosen automatically."))),
                   tabPanel("Model Comparison",
                            h4("Accuracy: training vs. validation"),
                            gt_output("acc_table")),
                   tabPanel("Forecasts",
                            br(),
                            checkboxGroupInput("fc_model", "Models to plot",
                                               choices = model_levels,
                                               selected = model_levels, inline = TRUE),
                            plotOutput("fc_plot"))
                 )                                  
        )                                           
      )                                             
    )                                               
  )                                                 
)                                                   


### Server
server <- function(input, output, session) {
  
  # data
  start_ym <- reactive({ req(input$dates); yearmonth(input$dates[1]) })
  end_ym <- reactive({ req(input$dates); yearmonth(input$dates[2]) })
  
  filtered <- reactive({
    req(input$varietal)
    validate(need(start_ym() < end_ym(),
                  "'From' must be earlier than 'To'."))
    wine_ts |>
      filter(Varietal %in% input$varietal,
             Month >= start_ym(), Month <= end_ym())
  })
  
  output$ts_print <- renderPrint({
    print(filtered(), n = 24)
  })
  
  # overviw
  output$overview <- renderPlot({
    filtered() |>
      autoplot(Sales) +
      labs(title = "Wine Sales by Month", y = "Sales (in thousands of litres)", 
           x = NULL, colour = "Varietal")
  })
  
  # seasonality
  season_data <- reactive({
    validate(need(n_distinct(filtered()$Month) >= 24,
                  "Seasonality needs at least 2 years of data"))
    filtered()
  })
  
  output$season_plot <- renderPlot({
    season_data() |>
      gg_season(Sales, labels = "none") +
      facet_wrap(~ Varietal, scales = "free_y") +
      labs(title = "Seasonal Plot", y = "Sales (in thousands of litres)", x = NULL)
  })
  
  output$stl_plot <- renderPlot({
    season_data() |>
      model(STL(Sales ~ season(window = "periodic"))) |>
      components() |>
      as_tibble() |>
      select(Varietal, Month, Sales, trend, season_year, remainder) |>
      pivot_longer(c(Sales, trend, season_year, remainder),
                   names_to = "Component", values_to = "Value") |>
      mutate(Component = factor(Component,
                                levels = c("Sales", "trend", "season_year", "remainder"),
                                labels = c("Data", "Trend", "Seasonal", "Remainder"))) |>
      ggplot(aes(x = Month, y = Value)) +
      geom_line() +
      facet_grid(Component ~ Varietal, scales = "free_y") +
      labs(title = "STL Decomposition", x = NULL, y = NULL)
  })
  
  # comparison models
  valid_start <- reactive({
    req(input$valid_n)
    max(filtered()$Month) - input$valid_n + 1
  })
  
  train_data <- reactive({
    d <- filtered() |> filter(Month < valid_start())
    validate(need(n_distinct(d$Month) >= 24,
                  "Need at least 2 years of training data. Widen the date range or shorten the validation window."))
    d
  })
  
  # model fit
  fit_tslm <- reactive({
    train_data() |> model(TSLM = TSLM(Sales ~ trend() + season()))
  }) |>
    bindCache("TSLM", input$varietal, input$dates, input$valid_n)
  
  fit_ets <- reactive({
    train_data() |> model(ETS = ETS(Sales))
  }) |>
    bindCache("ETS", input$varietal, input$dates, input$valid_n)
  
  fit_arima <- reactive({
    train_data() |> model(ARIMA = ARIMA(Sales ~ pdq(d = 1) + PDQ(D = 1)))
  }) |>
    bindCache("ARIMA", input$varietal, input$dates, input$valid_n)
  
  fit_ens <- reactive({
    train_data() |>
      model(Ensemble = combination_model(
        ETS(Sales),
        ARIMA(Sales ~ pdq(d = 1) + PDQ(D = 1))
      ))
  }) |>
    bindCache("Ensemble", input$varietal, input$dates, input$valid_n)
  
  get_spec <- function(fit, name) {
    tibble(Varietal = fit$Varietal,
           Specification = str_remove_all(format(fit[[name]]), "[<>]"))
  }
  
  # tslm
  output$tslm_table <- renderTable({
    fit <- fit_tslm()
    coefs <- tidy(fit) |>
      filter(term %in% c("(Intercept)", "trend()")) |>
      select(Varietal, term, estimate) |>
      pivot_wider(names_from = term, values_from = estimate) |>
      rename(Intercept = `(Intercept)`, `Trend (per month)` = `trend()`)
    fitq <- glance(fit) |> select(Varietal, `Adj. R-squared` = adj_r_squared)
    left_join(coefs, fitq, by = "Varietal")
  }, digits = 3)
  
  # ets
  output$ets_table <- renderTable({
    fit <- fit_ets()
    params <- tidy(fit) |>
      filter(term %in% c("alpha", "beta", "gamma", "phi")) |>
      select(Varietal, term, estimate) |>
      pivot_wider(names_from = term, values_from = estimate)
    left_join(get_spec(fit, "ETS"), params, by = "Varietal")
  }, digits = 3)
  
  # arima
  output$arima_table <- renderTable({
    fit <- fit_arima()
    coefs <- tidy(fit) |>
      select(Varietal, term, estimate) |>
      pivot_wider(names_from = term, values_from = estimate)
    left_join(get_spec(fit, "ARIMA"), coefs, by = "Varietal")
  }, digits = 3, na = "")
  
  # comparison
  acc_data <- reactive({
    n <- input$valid_n
    fits <- list(fit_tslm(), fit_ets(), fit_arima(), fit_ens())
    
    train_acc <- map_dfr(fits, accuracy)
    valid_acc <- map_dfr(fits, \(f) f |> forecast(h = n) |> accuracy(filtered()))
    
    bind_rows(train_acc, valid_acc) |>
      select(Varietal, Model = .model, Window = .type, RMSE, MAE, MAPE) |>
      mutate(Window = if_else(Window == "Test", "Validation", Window)) |>
      pivot_wider(names_from = Window,
                  values_from = c(RMSE, MAE, MAPE),
                  names_glue = "{Window}_{.value}")
  })
  
  output$acc_table <- render_gt({
    acc_data() |>
      select(Varietal, Model, Training_RMSE, Training_MAE, Training_MAPE,
             Validation_RMSE, Validation_MAE, Validation_MAPE) |>
      mutate(Model = factor(Model, levels = model_levels)) |>
      arrange(Varietal, Model) |>
      group_by(Varietal) |>
      gt() |>
      tab_spanner(label = "Model", columns = Model, id = "model_spanner") |>
      tab_spanner(label = "Training", columns = starts_with("Training_")) |>
      tab_spanner(label = "Validation", columns = starts_with("Validation_")) |>
      cols_label(Model = "", Training_RMSE = "RMSE", Training_MAE = "MAE", 
                 Training_MAPE = "MAPE", Validation_RMSE = "RMSE", 
                 Validation_MAE = "MAE", Validation_MAPE = "MAPE") |>
      fmt_number(columns = where(is.numeric), decimals = 2) |>
      tab_style(style = cell_text(weight = "bold"), 
                locations = list(cells_column_spanners(), cells_column_labels(),
                                 cells_row_groups()))
  })
  
  # forecasts
  make_fc <- function(fit, h) {
    fit |>
      forecast(h = h) |>
      hilo(level = c(80, 95)) |>
      unpack_hilo(c(`80%`, `95%`)) |>
      as_tibble() |>
      select(Varietal, Model = .model, Month, Forecast = .mean,
             lo80 = `80%_lower`, hi80 = `80%_upper`,
             lo95 = `95%_lower`, hi95 = `95%_upper`)
  }
  
  fc_data <- reactive({
    req(input$h)
    validate(need(length(input$fc_model) > 0, "Check at least one model"))
    all_fits <- list(TSLM     = fit_tslm,
                     ETS      = fit_ets,
                     ARIMA    = fit_arima,
                     Ensemble = fit_ens)
    map_dfr(input$fc_model, \(m) make_fc(all_fits[[m]](), input$h)) |>
      mutate(Model = factor(Model, levels = model_levels))
  })
  
  output$fc_plot <- renderPlot({
    fc <- fc_data()
    actuals <- filtered() |>
      filter(Month >= valid_start() - 36) |>
      as_tibble()
    
    ggplot() +
      geom_ribbon(data = fc, aes(x = Month, ymin = lo95, ymax = hi95, fill = Model),
                  alpha = 0.15) +
      geom_ribbon(data = fc, aes(x = Month, ymin = lo80, ymax = hi80, fill = Model),
                  alpha = 0.30) +
      geom_line(data = actuals, aes(x = Month, y = Sales), colour = "black") +
      geom_line(data = fc, aes(x = Month, y = Forecast, colour = Model),
                linewidth = 0.8) +
      geom_vline(xintercept = as.Date(valid_start()), linetype = "dashed") +
      scale_colour_discrete(drop = FALSE, limits = model_levels) +
      scale_fill_discrete(drop = FALSE, limits = model_levels) +
      facet_grid(Varietal ~ Model, scales = "free_y", switch = "y") +
      labs(title = paste0(input$h, "-month forecasts"),
           y = "Sales (in thousands of litres)", x = NULL) +
      theme(strip.placement = "outside",
            legend.position = "none")
  }, height = function() 120 + 260 * length(input$varietal))
}

shinyApp(ui, server)
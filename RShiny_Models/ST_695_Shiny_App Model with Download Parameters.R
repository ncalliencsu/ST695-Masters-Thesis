#setwd("C:/Users/norma/OneDrive/Documents/NC State/ST 695/R Models")
#source

# Load required libraries
library(shiny)        # For building interactive web apps
library(lmerTest)     # For linear mixed-effects models
library(ggplot2)      # For plotting
library(dplyr)        # For data manipulation
library(grid)         # For grid graphics
library(DT)           # For interactive data tables



# Define the UI (User Interface) section
ui <- fluidPage(
  titlePanel("Simulation Parameter Entry (Dynamic)"),  # App title
  sidebarLayout(
    sidebarPanel(
      # Numeric inputs for simulation parameters
      numericInput("HTC", "HTC (positive integer):", value = 2, min = 1, step = 1),
      numericInput("ETC", "ETC (positive integer):", value = 2, min = 1, step = 1),
      numericInput("r", "r (positive integer):", value = 2, min = 1, step = 1),
      numericInput("nsim", "nsim (positive integer):", value = 100, min = 1, step = 1),
      numericInput("Intercept", "Intercept (real number):", value = 10),
      # Dynamic UI for w, s, and interaction inputs
      uiOutput("w_inputs"),
      uiOutput("s_inputs"),
      uiOutput("interaction_inputs"),
      # Numeric inputs for variance parameters
      numericInput("sigma2_W", "sigma2_W (positive real):", value = 1, min = 0, step = 0.01),
      numericInput("sigma2_S", "sigma2_S (positive real):", value = 5, min = 0, step = 0.01),
      # Button to run the simulation
      actionButton("submit", "Run Simulation"),
      # Model Inadequate UI (moved to bottom)
      checkboxInput("model_inadequate", "Model Inadequate", value = FALSE),
      conditionalPanel(
        condition = "input.model_inadequate",
        sliderInput("inadequate_scale", "Model Inadequate Scale", min = 0.5, max = 1.5, value = 1, step = 0.01)
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Summary & Plots",
          h4("Summary Table:"),
          h4("Overlaid Histograms:"),
          # Output area for histograms
          uiOutput("all_histograms"),
          # Download buttons for results
          downloadButton("download_summary", "Download Summary"),
          downloadButton("download_histograms", "Download Histograms (PDF)"),
          downloadButton("download_params", "Download Parameters")
          # Lenth plot output and download (REMOVED)
        ),
        tabPanel("True Positive Rate (TPR)",
          h4("True Positive Rate (TPR) Table for All Models"),
          dataTableOutput("tpr_table"),
          downloadButton("download_tpr", "Download TPR Table")
        ),
        tabPanel("F-Test Results",
          h4("F-Test: SW_Model vs S_Model"),
          p("Testing H0: SW_Model"),
          dataTableOutput("Ftest_table"),
          downloadButton("download_F_Test", "Download F-Test Results")
        )
      )
    )
  )
)

# Define the Server Logic section
server <- function(input, output, session) {
  # Dynamically generate numeric inputs for w parameters
  output$w_inputs <- renderUI({
    lapply(1:input$HTC, function(i) {
      numericInput(paste0("w", i), paste0("w", i, " (real number):"), value = 0)
    })
  })
  # Dynamically generate numeric inputs for s parameters
  output$s_inputs <- renderUI({
    lapply(1:input$ETC, function(i) {
      numericInput(paste0("s", i), paste0("s", i, " (real number):"), value = 0)
    })
  })
  # Dynamically generate numeric inputs for interaction terms
  output$interaction_inputs <- renderUI({
    w_names <- paste0("w", 1:input$HTC)
    s_names <- paste0("s", 1:input$ETC)
    interaction_names <- c(
      if (length(w_names) > 1) combn(w_names, 2, FUN = function(x) paste(x, collapse = ":")),
      if (length(s_names) > 1) combn(s_names, 2, FUN = function(x) paste(x, collapse = ":")),
      as.vector(outer(w_names, s_names, paste, sep = ":"))
    )
    # Replace ':' with '_' for input IDs
    interaction_ids <- gsub(":", "_", interaction_names)
    lapply(seq_along(interaction_names), function(i) {
      numericInput(interaction_ids[i], paste0(interaction_names[i], " (real number):"), value = 0)
    })
  })

  # Main simulation and output logic, triggered when 'Run Simulation' is clicked
  observeEvent(input$submit, {
    # Extract and isolate user inputs
    HTC <- isolate(input$HTC)
    ETC <- isolate(input$ETC)
    r <- isolate(input$r)
    nsim <- isolate(input$nsim)
    Sigma2_W <- isolate(input$sigma2_W)
    Sigma2_S <- isolate(input$sigma2_S)
    w_names <- paste0("w", 1:HTC)
    s_names <- paste0("s", 1:ETC)
    # Apply Model Inadequate scaling if enabled
    if (!is.null(input$model_inadequate) && input$model_inadequate) {
      for (n in w_names) {
        assign(n, isolate(input[[n]]) * input$inadequate_scale)
      }
      for (n in s_names) {
        assign(n, isolate(input[[n]]) * input$inadequate_scale)
      }
    }
    interaction_names <- c(
      if (length(w_names) > 1) combn(w_names, 2, FUN = function(x) paste(x, collapse = ":")),
      if (length(s_names) > 1) combn(s_names, 2, FUN = function(x) paste(x, collapse = ":")),
      as.vector(outer(w_names, s_names, paste, sep = ":"))
    )
    interaction_ids <- gsub(":", "_", interaction_names)
    # Collect all effect parameters from user input
    effects <- c(
      isolate(input$Intercept),
      sapply(w_names, function(n) isolate(input[[n]])),
      sapply(s_names, function(n) isolate(input[[n]])),
      sapply(interaction_ids, function(n) isolate(input[[n]]))
    )
    # Create all combinations of w and s factors
    combined_list <- c(
      setNames(replicate(ETC, c(-1, 1), simplify = FALSE), s_names),
      setNames(replicate(HTC, c(-1, 1), simplify = FALSE), w_names)
    )
    df_S <- expand.grid(combined_list)
    df_rep_S <- df_S[rep(seq_len(nrow(df_S)), r), ]
    row.names(df_rep_S) <- NULL
    linear_terms_ws <- c(w_names, s_names)
    # Build model matrix for S model
    formula_S <- as.formula(
      paste("~ (", paste(linear_terms_ws, collapse = " + "), ")^2")
    )
    X_S <- model.matrix(formula_S, data = df_rep_S)
    colnames(X_S)[1] <- '0'
    # Build model matrix for W model
    wpf_terms <- setNames(replicate(HTC, c(-1,1), simplify = FALSE), w_names)
    df_W <- expand.grid(wpf_terms)
    df_rep_W <- df_W[rep(seq_len(nrow(df_W)), r), ]
    row.names(df_rep_W) <- NULL
    formula_W <- as.formula(
      paste("~ (", paste(w_names, collapse = " + "), ")^2")
    )
    # Extract relevant columns for W and SW models
    X_W <- X_S[, !grepl("s", colnames(X_S), ignore.case = TRUE)]
    X_SW <- X_S[, grep('s|0', colnames(X_S), value = TRUE)]
    colnames(X_SW) <- gsub(":", "_", colnames(X_SW))
    set.seed(1)  # For reproducibility

    #creates a factor variable WP for the whole-plot structure in the split-plot design.
    #repeats each whole-plot ID (from 1 to r*2^HTC) for each subplot (2^ETC times)
    WP <- as.factor(rep(seq(1, r*2^HTC), each = 2^ETC))

    # Prepare names for fixed effects
    fixed_effect_names_S <- paste0("beta", colnames(X_S))
    fixed_effect_names_W <- paste0("beta", colnames(X_W))

    # After building X_S, ensure effect assignment matches column order
    effect_names <- colnames(X_S)

    # Collect all effect parameters from user input in the order of effect_names
    # Build a named list of all possible effect values from user input
    effect_input_map <- list(Intercept = isolate(input$Intercept))
    for (n in w_names) effect_input_map[[n]] <- isolate(input[[n]])
    for (m in s_names) effect_input_map[[m]] <- isolate(input[[m]])
    for (i in seq_along(interaction_names))
      effect_input_map[[interaction_names[i]]] <- isolate(input[[interaction_ids[i]]])

    # For each effect_name, get the corresponding value from effect_input_map
    effects <- sapply(effect_names, function(nm) {
      # Remove spaces and replace ':' with '_' for matching input IDs
      nm_clean <- gsub(":", "_", nm)
      if (!is.null(effect_input_map[[nm]])) {
        effect_input_map[[nm]]
      } else if (!is.null(effect_input_map[[nm_clean]])) {
        effect_input_map[[nm_clean]]
      } else if (nm == "0") {
        effect_input_map[["Intercept"]]
      } else {
        0  # default to 0 if not found
      }
    })
    beta_S <- setNames(effects, effect_names)
    names(beta_S)[1] <- "beta*0"
    beta_W <- beta_S[!grepl("s", colnames(X_S), ignore.case = TRUE)]
    # Calculate means for S and W models
    mu_S <- X_S %*% beta_S
    mu_W <- X_W %*% beta_W

    # Simulation function for S model (mixed effects)
    sim_model_S <- function(X_S, mu_S, mu_W, Sigma2_S, Sigma2_W, WP) {
      n_groups <- length(unique(WP))
      epsilon_W <- rnorm(n_groups, mean = 0, sd = sqrt(Sigma2_W))
      epsilon_W <- rep(epsilon_W, times = table(WP))
      epsilon_S <- rnorm(nrow(X_S), mean = 0, sd = sqrt(Sigma2_S))
      # Model Inadequate logic
      if (!is.null(input$model_inadequate) && input$model_inadequate) {
        scalefactor <- input$inadequate_scale
        Y_W_star <- scalefactor * mu_W + epsilon_W
        Y_S <- mu_S + epsilon_W + epsilon_S + Y_W_star
        Y_W <- Y_W_star
      } else {
        Y_S <- mu_S + epsilon_W + epsilon_S
        Y_W <- mu_W + epsilon_W
      }

      df_S <- data.frame(epsilon_W, df_rep_S, Y_S)

      formula_S <- as.formula(
        paste("Y_S ~ (", paste(linear_terms_ws, collapse = " + "), ")^2 + (1|WP)"))
      rmodel_S <- lmer(formula_S, data = df_S)
      rmodel_S_summary <- summary(rmodel_S)
      rownames(rmodel_S_summary$coefficients) <-
        sub("\\(Intercept\\)", "Intercept", rownames(rmodel_S_summary$coefficients))


      fixed_effects_S <- rmodel_S_summary$coefficients[, "Estimate"]
      pvalues_S <- rmodel_S_summary$coefficients[, "Pr(>|t|)"]
      Sigma2_W <- (as.numeric(rmodel_S_summary$varcor$WP))
      Sigma2_S <- (rmodel_S_summary$sigma)^2
      dffe_S <- as.data.frame(t(fixed_effects_S))
      colnames(dffe_S) <- paste0("Fixed_Eff_Model_S (", names(fixed_effects_S), ")")
      dffp_S <- as.data.frame(t(pvalues_S))
      colnames(dffp_S) <- paste0("pvalue_Model_S (", names(pvalues_S), ")")
      df_out_S <- cbind(dffe_S, dffp_S)
      df_out_S$Sigma2_W <- Sigma2_W
      df_out_S$Sigma2_S <- Sigma2_S
      colnames(df_out_S)[colnames(df_out_S) == "Sigma2_W"] <- "Sigma2_W (Model S)"
      colnames(df_out_S)[colnames(df_out_S) == "Sigma2_S"] <- "Sigma2_S (Model S)"
      return(list(df_out_S = df_out_S, Y_W = Y_W, Y_S = Y_S, model_S = rmodel_S))
    }
    # Simulation function for W model (fixed effects)
    sim_model_W <- function(Y_W) {
      Y_W_unique <- unique(Y_W)
      df_W <- data.frame(df_rep_W, Y_W = Y_W_unique)
      formula_W <- as.formula(
        paste("Y_W ~ (", paste(w_names, collapse = " + "), ")^2"))
      rmodel_W <- summary(lm(formula_W, data = df_W))
      rownames(rmodel_W$coefficients) <- sub("\\(Intercept\\)", "Intercept", rownames(rmodel_W$coefficients))
      fixed_effects_W <- rmodel_W$coefficients[, "Estimate"]
      pvalues_W <- rmodel_W$coefficients[, "Pr(>|t|)"]
      Sigma2_W <- (rmodel_W$sigma)^2
      dffe_W <- as.data.frame(t(fixed_effects_W))
      colnames(dffe_W) <- paste0("Fixed_Eff_Model_W (", names(fixed_effects_W),")")
      dffp_W <- as.data.frame(t(pvalues_W))
      colnames(dffp_W) <- paste0("pvalue_Model_W (", names(fixed_effects_W), ")")
      df_out_W <- cbind(dffe_W, dffp_W)
      df_out_W$Sigma2_W <- Sigma2_W
      colnames(df_out_W)[colnames(df_out_W) == "Sigma2_W"] <- "Sigma2_W (Model W)"
      return(df_out_W)
    }

    # Simulation function for SW Models
    # SW1 is the null model built by subtracting the W Model from the S Model.  SW2 is the alternative model generated by fitting the S Model to SW responses.
    sim_model_SW <- function(X_SW, Y_S, Y_W, X_S) {
      X_SW <- X_SW[,-1]
      df_SW1 <- data.frame(Y_SW1 = Y_S - Y_W, X_SW)

      X_S <- X_S[,-1]
      df_SW2 <- data.frame(Y_SW2 = Y_S - Y_W, X_S)

      #Model SW1 - (S Model - W Model)
      formula_SW1 <- as.formula(
        paste("Y_SW1 ~", paste(colnames(X_SW), collapse = " + "))
      )

      rmodel_SW1 <- lm(formula_SW1, data = df_SW1)
      rmodel_SW1_summary <- summary(rmodel_SW1)
      rownames(rmodel_SW1_summary$coefficients) <-
        sub("\\(Intercept\\)", "Intercept", rownames(rmodel_SW1_summary$coefficients))
      fixed_effects_SW1 <- rmodel_SW1_summary$coefficients[, "Estimate"]
      pvalues_SW1 <- rmodel_SW1_summary$coefficients[, "Pr(>|t|)"]
      Sigma2_S <- (rmodel_SW1_summary$sigma)^2
      dffe_SW1 <- as.data.frame(t(fixed_effects_SW1))
      colnames(dffe_SW1) <- paste0("Fixed_Eff_Model_SW1 (", names(fixed_effects_SW1),")")
      dffp_SW1 <- as.data.frame(t(pvalues_SW1))
      colnames(dffp_SW1) <- paste0("pvalue_Model_SW1 (", names(pvalues_SW1), ")")
      df_out_SW1 <- cbind(dffe_SW1, dffp_SW1)
      df_out_SW1$Sigma2_S <- Sigma2_S
      colnames(df_out_SW1)[colnames(df_out_SW1) == "Sigma2_S"] <- "Sigma2_S (Model SW1)"

      #Model SW2: - Fitting S Model to SW Responses
      # formula_SW2 <- as.formula(
      #   paste("Y_SW2 ~ (", paste(colnames(X_S), collapse = " + "), ") + (1|WP)"))
      # rmodel_SW2 <- lmer(formula_SW2, data = df_SW2)
      formula_SW2 <- as.formula(paste("Y_SW1 ~ (", paste(colnames(X_SW), collapse = " + "), ") + WP"))
      rmodel_SW2 <- lm(formula_SW2, data = df_SW1)
      rmodel_SW2_summary <- summary(rmodel_SW2)
      rownames(rmodel_SW2_summary$coefficients) <-
        sub("\\(Intercept\\)", "Intercept", rownames(rmodel_SW2_summary$coefficients))

      return(list(df_out_SW1 = df_out_SW1, model_SW1 = rmodel_SW1, model_SW2 = rmodel_SW2))
    }

    # Function to perform F-test comparing S_Model vs SW_Model
    perform_F_test <- function(model_SW1, model_SW2) {
      tryCatch({

        # Extract relevant statistics
        SSE_SW1 <- sum(residuals(model_SW1)^2)
        SSE_SW2 <- sum(residuals(model_SW2)^2)
        dfE_SW1 <- model_SW1$df.residual
        dfE_SW2 <- model_SW2$df.residual

        f_statistic <- ((SSE_SW1 - SSE_SW2) / (dfE_SW1 - dfE_SW2)) / (SSE_SW2 / dfE_SW2)
        Fcrit <- qf(0.95, dfE_SW1 - dfE_SW2, dfE_SW2) # upper 5% critical value

        Significant <- FALSE
        if (f_statistic > Fcrit) Significant <- TRUE

        data.frame(
          F_statistic = f_statistic,
          Fcrit = Fcrit,
          Significant = Significant
        )
      }, error = function(e) {
        # Return NA values if test fails
        data.frame(
          F_statistic = NA,
          Fcrit = NA,
          Significant = NA
        )
      })
    }

    # Helper to organize model outputs for summary and plotting
    model_gen <- function(sim_df_S, sim_df_W, sim_df_SW1) {
      S_Model <- sim_df_S[, !grepl("pvalue", colnames(sim_df_S), ignore.case = TRUE)]
      W_Model <- sim_df_W[, !grepl("pvalue", colnames(sim_df_W), ignore.case = TRUE)]
      SW_Model <- sim_df_SW1[, !grepl("pvalue", colnames(sim_df_SW1), ignore.case = TRUE)]
      colnames(S_Model) <- c(colnames(X_S), rep(NA, ncol(S_Model) - length(colnames(X_S))))
      colnames(S_Model)[1] <- "Intercept"
      colnames(S_Model)[(ncol(S_Model)-1):ncol(S_Model)] <- c("Sigma2_W", "Sigma2_S")
      colnames(W_Model) <- c(colnames(X_W), rep(NA, ncol(W_Model) - length(colnames(X_W))))
      colnames(W_Model)[1] <- "Intercept"
      colnames(W_Model)[ncol(W_Model)] <- "Sigma2_W"
      colnames(X_SW) <- gsub("_", ":", colnames(X_SW))
      colnames(SW_Model) <- c(colnames(X_SW), rep(NA, ncol(SW_Model) - length(colnames(X_SW))))
      colnames(SW_Model)[1] <- "Intercept"
      colnames(SW_Model)[ncol(SW_Model)] <- "Sigma2_S"
      return(list(S_Model = S_Model, W_Model = W_Model, SW_Model = SW_Model))
    }
    # Function to plot overlaid histograms for each parameter
    plot_overlaid_histograms <- function(S_Model, W_Model, SW_Model, bins = 20) {
      plots <- list()
      # Build a map of all effect values, including interactions, using effect_names and effect_input_map
      param_map <- list()
      for (nm in colnames(S_Model)) {
        # Clean up for matching input IDs
        nm_clean <- gsub(":", "_", nm)
        if (!is.null(effect_input_map[[nm]])) {
          param_map[[nm]] <- effect_input_map[[nm]]
        } else if (!is.null(effect_input_map[[nm_clean]])) {
          param_map[[nm]] <- effect_input_map[[nm_clean]]
        } else if (nm == "0" || nm == "Intercept") {
          param_map[[nm]] <- effect_input_map[["Intercept"]]
        } else if (nm == "Sigma2_W") {
          param_map[[nm]] <- input$sigma2_W
        } else if (nm == "Sigma2_S") {
          param_map[[nm]] <- input$sigma2_S
        } else {
          param_map[[nm]] <- NA
        }
      }
      for (col in colnames(S_Model)) {
        df_long <- bind_rows(
          data.frame(value = S_Model[[col]], model = "S_Model"),
          if (col %in% colnames(W_Model)) data.frame(value = W_Model[[col]], model = "W_Model"),
          if (col %in% colnames(SW_Model)) data.frame(value = SW_Model[[col]], model = "SW_Model")
        )
        param_value <- param_map[[col]]
        label_text <- if (!is.null(param_value) && !is.na(param_value)) paste0(col, " = ", param_value) else col
        p <- ggplot(df_long, aes(x = value, fill = model, color = model)) +
          geom_histogram(position = "identity", alpha = 0.4, bins = bins) +
          labs(
            title = paste("Comparison:", col),
            x = col,
            subtitle = label_text
          ) +
          theme_minimal() +
          theme(plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5))
        plots[[col]] <- p
      }
      return(plots)
    }
    # Function to summarize statistics for each parameter
    summary_table <- function(S_Model, W_Model, SW_Model) {
      all_stats <- list()
      for (col in colnames(S_Model) ) {
        stats_list <- list()
        stats_list[["S_Model"]] <- summary(S_Model[[col]])
        if (col %in% colnames(W_Model)) {
          stats_list[["W_Model"]] <- summary(W_Model[[col]])
        }
        if (col %in% colnames(SW_Model)) {
          stats_list[["SW_Model"]] <- summary(SW_Model[[col]])
        }
        stats_df <- tryCatch({
          as.data.frame(stats_list)
        }, error = function(e) {
          stats_list
        })
        all_stats[[col]] <- stats_df
      }
      return(all_stats)
    }
    # Run the simulation nsim times and collect results
    sim_df_S <- data.frame()
    sim_df_W <- data.frame()
    sim_df_SW1 <- data.frame()
    Ftest_result <- data.frame()



    for(i in 1:nsim){
      sim_model_S_out <- sim_model_S(X_S, mu_S, mu_W, Sigma2_S = Sigma2_S, Sigma2_W = Sigma2_W, WP)
      Y_W <- sim_model_S_out$Y_W
      Y_S <- sim_model_S_out$Y_S

      model_S <- sim_model_S_out$model_S
      rmodel_S_df <- as.data.frame(sim_model_S_out$df_out_S)
      sim_df_S <- rbind(sim_df_S, rmodel_S_df)

      rmodel_W_df <- as.data.frame(sim_model_W(Y_W))
      sim_df_W <- rbind(sim_df_W, rmodel_W_df)

      sim_model_SW_out <- sim_model_SW(X_SW, Y_S, Y_W, X_S)
      model_SW1 <- sim_model_SW_out$model_SW1
      rmodel_SW1_df <- as.data.frame(sim_model_SW_out$df_out_SW1)
      sim_df_SW1 <- rbind(sim_df_SW1, rmodel_SW1_df)

      model_SW2 <- sim_model_SW_out$model_SW2

      F_Test <- perform_F_test(model_SW1, model_SW2)
      Ftest_result <- rbind(Ftest_result, F_Test)


  }

    # Organize model outputs
    models <- model_gen(sim_df_S, sim_df_W, sim_df_SW1)
    S_Model <- models$S_Model
    W_Model <- models$W_Model
    SW_Model <- models$SW_Model

    # Always plot all histograms (Lenth logic removed)
    plots <- plot_overlaid_histograms(S_Model, W_Model, SW_Model, bins = 20)

    #Download Handler for Summary(CSV)
    output$download_summary <- downloadHandler(
      filename = function() {
        paste0("summary_table_", Sys.Date(), ".csv")
      },
      content = function(file) {
        summary_stats <- summary_table(S_Model, W_Model, SW_Model)

        # Simple approach - just take the first element and write it
        if (length(summary_stats) > 0) {
          write.csv(summary_stats[[1]], file, row.names = TRUE)
        } else {
          # Fallback - create empty CSV
          write.csv(data.frame(Note = "No data available"), file, row.names = FALSE)
        }
      }
    )

    # Download handler for histograms (PDF)
    output$download_histograms <- downloadHandler(
      filename = function() {
        paste0("histograms_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        pdf(file)
        for (p in plots) {
          print(p)
        }
        dev.off()
      }
    )
    # Download handler for simulation parameters
    output$download_summary <- downloadHandler(
      filename = function() {
        paste0("summary_table_", Sys.Date(), ".csv")
      },
      content = function(file) {
        summary_stats <- summary_table(S_Model, W_Model, SW_Model)

        # Create a more robust data frame conversion
        combined_list <- list()

        for (param_name in names(summary_stats)) {
          param_data <- summary_stats[[param_name]]

          # Handle different data structures that might be returned
          if (is.data.frame(param_data)) {
            # If it's already a data frame, add parameter column
            param_data$Parameter <- param_name
            combined_list[[param_name]] <- param_data
          } else if (is.list(param_data)) {
            # If it's a list of summaries, convert to data frame
            df_rows <- list()
            for (model_name in names(param_data)) {
              if (is.numeric(param_data[[model_name]])) {
                # Convert summary statistics to a single row
                summary_row <- as.data.frame(t(param_data[[model_name]]))
                summary_row$Model <- model_name
                summary_row$Parameter <- param_name
                df_rows[[model_name]] <- summary_row
              }
            }
            if (length(df_rows) > 0) {
              combined_list[[param_name]] <- do.call(rbind, df_rows)
            }
          }
        }

        # Combine all parameter data
        if (length(combined_list) > 0) {
          final_df <- do.call(rbind, combined_list)
          rownames(final_df) <- NULL
        } else {
          # Fallback: create a simple summary
          final_df <- data.frame(
            Parameter = names(summary_stats),
            Note = "Summary data structure not compatible for CSV export"
          )
        }

        write.csv(final_df, file, row.names = FALSE)
      }
    )


    # Render all histogram plots in a grid layout
    output$all_histograms <- renderUI({
      n_cols <- 3
      n_plots <- length(plots)
      rows <- ceiling(n_plots / n_cols)
      plot_output_list <- list()
      for (row in 1:rows) {
        cols <- lapply(1:n_cols, function(col) {
          idx <- (row - 1) * n_cols + col
          if (idx <= n_plots) {
            plotOutput(paste0("plot", idx), height = "250px")
          } else {
            NULL
          }
        })
        plot_output_list[[row]] <- fluidRow(cols)
      }
      do.call(tagList, plot_output_list)
    })

    # Render each individual plot
    for (i in seq_along(plots)) {
      local({
        my_i <- i
        plotname <- paste0("plot", my_i)
        output[[plotname]] <- renderPlot({
          plots[[my_i]]
        })
      })
    }

    # --- True Positive Rate (TPR) Table and Download ---
    output$tpr_table <- renderDataTable({
      req(nsim)
      # Extract all p-value columns for each model
      pval_cols_S <- grep("^pvalue_Model_S", colnames(sim_df_S), value = TRUE)
      pval_cols_W <- grep("^pvalue_Model_W", colnames(sim_df_W), value = TRUE)
      pval_cols_SW <- grep("^pvalue_Model_SW", colnames(sim_df_SW1), value = TRUE)
      # Calculate TPR for each effect (proportion of p-values <= 0.05)
      tpr_S <- colMeans(sim_df_S[, pval_cols_S, drop = FALSE] <= 0.05, na.rm = TRUE)
      tpr_W <- colMeans(sim_df_W[, pval_cols_W, drop = FALSE] <= 0.05, na.rm = TRUE)
      tpr_SW <- colMeans(sim_df_SW1[, pval_cols_SW, drop = FALSE] <= 0.05, na.rm = TRUE)
      # Combine into a single table
      tpr_table <- data.frame(
        Parameter = c(names(tpr_S), names(tpr_W), names(tpr_SW)),
        Model = c(rep("S", length(tpr_S)), rep("W", length(tpr_W)), rep("SW", length(tpr_SW))),
        TPR = c(tpr_S, tpr_W, tpr_SW)
      )
      tpr_table <- tpr_table[order(tpr_table$Model, tpr_table$Parameter), ]
      DT::datatable(tpr_table, options = list(pageLength = 15, autoWidth = TRUE))
    })
    # Download handler for TPR table
    output$download_tpr <- downloadHandler(
      filename = function() {
        paste0("tpr_table_", Sys.Date(), ".csv")
      },
      content = function(file) {
        pval_cols_S <- grep("^pvalue_Model_S", colnames(sim_df_S), value = TRUE)
        pval_cols_W <- grep("^pvalue_Model_W", colnames(sim_df_W), value = TRUE)
        pval_cols_SW <- grep("^pvalue_Model_SW", colnames(sim_df_SW1), value = TRUE)
        tpr_S <- colMeans(sim_df_S[, pval_cols_S, drop = FALSE] <= 0.05, na.rm = TRUE)
        tpr_W <- colMeans(sim_df_W[, pval_cols_W, drop = FALSE] <= 0.05, na.rm = TRUE)
        tpr_SW <- colMeans(sim_df_SW1[, pval_cols_SW, drop = FALSE] <= 0.05, na.rm = TRUE)
        tpr_table <- data.frame(
          Parameter = c(names(tpr_S), names(tpr_W), names(tpr_SW)),
          Model = c(rep("S", length(tpr_S)), rep("W", length(tpr_W)), rep("SW", length(tpr_SW))),
          TPR = c(tpr_S, tpr_W, tpr_SW)
        )
        tpr_table <- tpr_table[order(tpr_table$Model, tpr_table$Parameter), ]
        write.csv(tpr_table, file, row.names = FALSE)
      }
    )

    # --- F-Test Results Table and Download ---

    output$Ftest_table <- renderDataTable({
       req(nrow(Ftest_result) > 0)

      # Create summary statistics for F Test results
      Ftest_summary <- data.frame(
        Statistic = c(
          "Mean F-Statistic Value",
          "Mean F-critical Value",
          "Proportion Significant",
          "Number of Successful Tests",
          "Number of Failed Tests"
        ),
        Value = c(
          round(mean(Ftest_result$F_statistic, na.rm = TRUE), 4),
          round(mean(Ftest_result$Fcrit, na.rm = TRUE), 4),
          round(mean(Ftest_result$Significant, na.rm = TRUE), 4),
          sum(!is.na(Ftest_result$F_statistic)),
          sum(is.na(Ftest_result$F_statistic))
        )
      )

      DT::datatable(Ftest_summary,
                    options = list(pageLength = 15, autoWidth = TRUE, searching = FALSE),
                    caption = "Summary of F-Test Results: S_Model vs SW_Model")
    })

    # Download handler for F Test results
    output$download_F_Test <- downloadHandler(
      filename = function() {
        paste0("Ftest_results_", Sys.Date(), ".csv")
      },
      content = function(file) {
        # Create detailed results with summary
        Ftest_summary <- data.frame(
          Statistic = c(
            "Mean F-Statistic Value",
            "Mean F-Critical Value",
            "Proportion Significant",
            "Number of Successful Tests",
            "Number of Failed Tests"
          ),
          Value = c(
            mean(Ftest_result$F_statistic, na.rm = TRUE),
            mean(Ftest_result$Fcrit, na.rm = TRUE),
            mean(Ftest_result$Significant, na.rm = TRUE),
            sum(!is.na(Ftest_result$F_statistic)),
            sum(is.na(Ftest_result$F_statistic))
          )
        )

        # Write summary and detailed results to file
        writeLines("F-Test Summary: S_Model vs SW_Model", file)
        writeLines("", file)
        write.table(Ftest_summary, file, append = TRUE, sep = ",", row.names = FALSE)
        writeLines("", file)
        writeLines("Detailed Results by Simulation Run:", file)
        write.table(Ftest_result, file, append = TRUE, sep = ",", row.names = FALSE)
      }
    )
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)

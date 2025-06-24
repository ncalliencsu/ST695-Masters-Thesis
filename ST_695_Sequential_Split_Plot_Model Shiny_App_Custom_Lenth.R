#setwd("C:/Users/norma/OneDrive/Documents/NC State/ST 695/R Models")
#source

# Load required libraries
library(shiny)        # For building interactive web apps
library(lmerTest)     # For linear mixed-effects models
library(ggplot2)      # For plotting
library(dplyr)        # For data manipulation
library(grid)         # For grid graphics
library(BsMD)         # For Lenth Plot



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
      actionButton("submit", "Run Simulation")
    ),
    mainPanel(
      h4("Summary Table:"),
      h4("Overlaid Histograms:"),
      # Output area for histograms
      uiOutput("all_histograms"),
      # Download buttons for results
      downloadButton("download_summary", "Download Summary"),
      downloadButton("download_histograms", "Download Histograms (PDF)"),
      downloadButton("download_params", "Download Parameters"),
      # Lenth plot output and download
      uiOutput("lenth_plot_ui")
    )
  )
)

# Custom Lenth half-normal plot function for Shiny app
custom_lenth_plot <- function(mean_effects, ME, SME) {
  require(ggplot2)
  # Remove NA and zero effects
  effects <- mean_effects[!is.na(mean_effects) & mean_effects != 0]
  abs_effects <- abs(effects)
  # Half-normal quantiles
  n <- length(abs_effects)
  q <- qnorm((1:n - 0.5) / n)
  df <- data.frame(
    abs_effect = sort(abs_effects),
    quantile = sort(q),
    label = names(sort(abs_effects))
  )
  p <- ggplot(df, aes(x = quantile, y = abs_effect, label = label)) +
    geom_point(size = 2) +
    geom_text(nudge_y = 0.05 * max(df$abs_effect), size = 3, hjust = 0, check_overlap = TRUE) +
    geom_hline(yintercept = ME, color = 'blue', linetype = 'dashed', size = 1) +
    geom_hline(yintercept = SME, color = 'red', linetype = 'dotted', size = 1) +
    annotate('text', x = min(df$quantile), y = ME, label = paste0('ME = ', signif(ME, 3)), vjust = -1, hjust = 0, color = 'blue', size = 4) +
    annotate('text', x = min(df$quantile), y = SME, label = paste0('SME = ', signif(SME, 3)), vjust = -1, hjust = 0, color = 'red', size = 4) +
    labs(
      title = 'Lenth Half-Normal Plot',
      x = 'Half-Normal Quantile',
      y = 'Absolute Effect'
    ) +
    theme_minimal()
  return(p)
}

# Define the server logic section
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
    for (n in s_names) effect_input_map[[n]] <- isolate(input[[n]])
    for (i in seq_along(interaction_names)) effect_input_map[[interaction_names[i]]] <- isolate(input[[interaction_ids[i]]])
    # Now, for each effect_name, get the corresponding value from effect_input_map
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
      Y_S <- mu_S + epsilon_W + epsilon_S
      Y_W <- mu_W + epsilon_W
      df_S <- data.frame(epsilon_W, df_rep_S, Y_S)
      formula_S <- as.formula(
        paste("Y_S ~ (", paste(linear_terms_ws, collapse = " + "), ")^2 + (1|WP)"))
      rmodel_S <- summary(lmer(formula_S, data = df_S))
      rownames(rmodel_S$coefficients) <- sub("\\(Intercept\\)", "Intercept", rownames(rmodel_S$coefficients))
      
    
      fixed_effects_S <- rmodel_S$coefficients[, "Estimate"]
      pvalues_S <- rmodel_S$coefficients[, "Pr(>|t|)"]
      Sigma2_W <- (as.numeric(rmodel_S$varcor$WP))
      Sigma2_S <- (rmodel_S$sigma)^2
      dffe_S <- as.data.frame(t(fixed_effects_S))
      colnames(dffe_S) <- paste0("Fixed_Eff_Model_S (", names(fixed_effects_S), ")")
      dffp_S <- as.data.frame(t(pvalues_S))
      colnames(dffp_S) <- paste0("pvalue_Model_S (", names(pvalues_S), ")")
      df_out_S <- cbind(dffe_S, dffp_S)
      df_out_S$Sigma2_W <- Sigma2_W
      df_out_S$Sigma2_S <- Sigma2_S
      colnames(df_out_S)[colnames(df_out_S) == "Sigma2_W"] <- "Sigma2_W (Model S)"
      colnames(df_out_S)[colnames(df_out_S) == "Sigma2_S"] <- "Sigma2_S (Model S)"
      return(list(df_out_S = df_out_S, Y_W = Y_W, Y_S = Y_S))
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
    # Simulation function for SW model (fixed effects on S terms)
    sim_model_SW <- function(X_SW, Y_S, Y_W) {
      X_SW <- X_SW[,-1]
      df_SW <- data.frame(Y_SW = Y_S - Y_W, X_SW)
      formula_SW <- as.formula(
        paste("Y_SW ~", paste(colnames(X_SW), collapse = " + "))
      )
      rmodel_SW <- summary(lm(formula_SW, data = df_SW))
      rownames(rmodel_SW$coefficients) <- sub("\\(Intercept\\)", "Intercept", rownames(rmodel_SW$coefficients))
      fixed_effects_SW <- rmodel_SW$coefficients[, "Estimate"]
      pvalues_SW <- rmodel_SW$coefficients[, "Pr(>|t|)"]
      Sigma2_S <- (rmodel_SW$sigma)^2
      dffe_SW <- as.data.frame(t(fixed_effects_SW))
      colnames(dffe_SW) <- paste0("Fixed_Eff_Model_SW (", names(fixed_effects_SW),")")
      dffp_SW <- as.data.frame(t(pvalues_SW))
      colnames(dffp_SW) <- paste0("pvalue_Model_SW (", names(pvalues_SW), ")")
      df_out_SW <- cbind(dffe_SW, dffp_SW)
      df_out_SW$Sigma2_S <- Sigma2_S
      colnames(df_out_SW)[colnames(df_out_SW) == "Sigma2_S"] <- "Sigma2_S (Model SW)"
      return(df_out_SW)
    }
    # Helper to organize model outputs for summary and plotting
    model_gen <- function(sim_df_S, sim_df_W, sim_df_SW) {
      S_Model <- sim_df_S[, !grepl("pvalue", colnames(sim_df_S), ignore.case = TRUE)]
      W_Model <- sim_df_W[, !grepl("pvalue", colnames(sim_df_W), ignore.case = TRUE)]
      SW_Model <- sim_df_SW[, !grepl("pvalue", colnames(sim_df_SW), ignore.case = TRUE)]
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
    sim_df_SW <- data.frame()
    for(i in 1:nsim){
      sim_model_S_out <- sim_model_S(X_S, mu_S, mu_W, Sigma2_S = Sigma2_S, Sigma2_W = Sigma2_W, WP)
      Y_W <- sim_model_S_out$Y_W
      Y_S <- sim_model_S_out$Y_S
      rmodel_S <- as.data.frame(sim_model_S_out$df_out_S)
      sim_df_S <- rbind(sim_df_S, rmodel_S)
      rmodel_W <- as.data.frame(sim_model_W(Y_W))
      sim_df_W <- rbind(sim_df_W, rmodel_W)
      rmodel_SW <- as.data.frame(sim_model_SW(X_SW, Y_S, Y_W))
      sim_df_SW <- rbind(sim_df_SW, rmodel_SW)
    }
    # Organize model outputs
    models <- model_gen(sim_df_S, sim_df_W, sim_df_SW)
    S_Model <- models$S_Model
    W_Model <- models$W_Model
    SW_Model <- models$SW_Model

    # Function to run Lenth method, plot, and return significant effect names and thresholds
    run_lenth <- function(S_Model, r) {
      if (r != 1) return(list(sig_effect_names = NULL, ME = NA, SME = NA, lenth_result = NULL))
      effect_cols <- !(colnames(S_Model) %in% c("Intercept", "Sigma2_W", "Sigma2_S"))
      mean_effects <- colMeans(S_Model[, effect_cols], na.rm = TRUE)
      lenth_result <- BsMD::LenthPlot(mean_effects)  # will plot
      ME <- lenth_result["ME"]
      SME <- lenth_result["SME"]
      significant_ME <- abs(mean_effects) > ME
      sig_effect_names <- names(mean_effects)[significant_ME]
      cat("Lenth Method Thresholds:\n")
      cat("ME  =", ME, "\n")
      cat("SME =", SME, "\n\n")
      cat("Significant effects by Lenth method (|effect| > ME):\n")
      print(sig_effect_names)
      return(list(sig_effect_names = sig_effect_names, ME = ME, SME = SME, lenth_result = lenth_result, mean_effects = mean_effects))
    }

    # Run Lenth method and get significant effect names and thresholds (only if r == 1)
    lenth_out <- run_lenth(S_Model, r)
    sig_effect_names <- lenth_out$sig_effect_names
    ME <- lenth_out$ME
    SME <- lenth_out$SME
    lenth_result <- lenth_out$lenth_result
    mean_effects <- lenth_out$mean_effects
    # Optionally, you can display ME, SME, and coefficients in the UI or console

    # If significant effects found, always include Sigma2_W and Sigma2_S for both models
    if (!is.null(sig_effect_names) && length(sig_effect_names) > 0) {
      cols_to_plot <- c("Intercept", sig_effect_names, "Sigma2_W", "Sigma2_S")
      cols_to_plot <- intersect(cols_to_plot, colnames(S_Model))
      S_Model_sig <- S_Model[, cols_to_plot, drop = FALSE]
      w_cols <- unique(c(intersect(cols_to_plot, colnames(W_Model)), "Sigma2_W"))
      w_cols <- w_cols[w_cols %in% colnames(W_Model)]
      W_Model_sig <- W_Model[, w_cols, drop = FALSE]
      SW_Model_sig <- SW_Model[, intersect(cols_to_plot, colnames(SW_Model)), drop = FALSE]
      plots <- plot_overlaid_histograms(S_Model_sig, W_Model_sig, SW_Model_sig, bins = 20)
    } else {
      plots <- plot_overlaid_histograms(S_Model, W_Model, SW_Model, bins = 20)
    }
    # Download handler for summary statistics
    output$download_summary <- downloadHandler(
      filename = function() {
        paste0("summary_table_", Sys.Date(), ".csv")
      },
      content = function(file) {
        summary_stats <- summary_table(S_Model, W_Model, SW_Model)
        write.csv(summary_stats[[1]], file)
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
    output$download_params <- downloadHandler(
      filename = function() {
        paste0("simulation_parameters_", Sys.Date(), ".txt")
      },
      content = function(file) {
        param_lines <- c(
          "Simulation Parameters:",
          paste("HTC =", input$HTC),
          paste("ETC =", input$ETC),
          paste("r =", input$r),
          paste("nsim =", input$nsim),
          paste("Intercept =", input$Intercept),
          sapply(w_names, function(n) paste(n, "=", input[[n]])),
          sapply(s_names, function(n) paste(n, "=", input[[n]])),
          sapply(seq_along(interaction_names), function(i) paste(interaction_names[i], "=", input[[interaction_ids[i]]])),
          paste("sigma2_W =", input$sigma2_W),
          paste("sigma2_S =", input$sigma2_S)
        )
        writeLines(param_lines, file)
      }
    )
    # --- Lenth Plot Output (UI and Server) ---
      output$lenth_plot <- renderPlot({
        if (r == 1 && !is.null(mean_effects) && !is.na(ME) && !is.na(SME)) {
          custom_lenth_plot(mean_effects, ME, SME)
        }
      })
      output$download_lenth_plot <- downloadHandler(
        filename = function() {
          paste0("lenth_plot_", Sys.Date(), ".pdf")
        },
        content = function(file) {
          pdf(file)
          if (r == 1 && !is.null(mean_effects) && !is.na(ME) && !is.na(SME)) {
            print(custom_lenth_plot(mean_effects, ME, SME))
          }
          dev.off()
        }
      )
      # --- End Lenth Plot Output ---

      # Add Lenth plot and download button to UI
      output$lenth_plot_ui <- renderUI({
        if (r == 1) {
          tagList(
            h4("Lenth Half-Normal Plot (r = 1 only):"),
            plotOutput("lenth_plot", height = "350px"),
            downloadButton("download_lenth_plot", "Download Lenth Plot (PDF)")
          )
        }
      })
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
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)

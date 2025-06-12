## app.R: Shiny App for Mixture Fitting of Expert‑Elicited Summaries

library(shiny)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(openxlsx)
library(ggplot2)
library(plotly)
library(DT)
library(tibble)

#––– 0) Settings –––

default_seed   <- 1234
default_n_boot <- 25
default_maxit <- 50000
probs <- c(
  p1  = 0.01,
  p25 = 0.25,
  p75 = 0.75,
  p99 = 0.99
)

#––– 1) Objective functions –––

obj_beta <- function(par, q_obs, mean_obs) {
  a <- par[1]
  b <- par[2]
  qf <- suppressWarnings(qbeta(probs, a, b))
  sum((qf - q_obs)^2) +
    (a/(a + b) - mean_obs)^2
}

obj_gamma <- function(par, q_obs, mean_obs) {
  shape <- par[1]
  rate  <- par[2]
  qf <- suppressWarnings(
    qgamma(probs, shape = shape, rate = rate)
  )
  sum((qf - q_obs)^2) +
    ((shape/rate) - mean_obs)^2
}

#––– 2) Fit one expert block (with convergence checks) –––
fit_one <- function(stats, distn, maxit = default_maxit) {
  if (length(unique(stats)) == 1) {
    return(list(dist = "point", value = stats[1],
                conv1 = NA_integer_, conv2 = NA_integer_))
  }
  
  q_obs <- stats[c("p1","p25","p75","p99")]
  m_obs <- stats["mean"]
  distn <- tolower(distn)
  
  # choose objective
  obj       <- if (distn=="gamma") obj_gamma else obj_beta
  par_names <- if (distn=="gamma") c("shape","rate") else c("alpha","beta")
  
  # 2.1) rough search
  par0 <- runif(2, 0.1, 5)
  s0   <- optim(par0, obj, q_obs = q_obs, mean_obs = m_obs,
                method = "SANN",
                control = list(maxit = maxit, reltol = 1e-10))
  conv1 <- s0$convergence
  if (conv1 != 0) {
    warning(sprintf(
      "fit_one[%s]: SANN failed to converge (code=%d)",
      distn, conv1
    ))
  }
  
  # 2.2) refine via L-BFGS-B, but never abort on error
  f0 <- tryCatch({
    optim(s0$par, obj, q_obs = q_obs, mean_obs = m_obs,
          method = "L-BFGS-B", lower = c(1e-8,1e-8), control = list(
            maxit = maxit,   # allow up to fifty thousand iterations
            factr = 1e6,     # tighter tolerance
            pgtol = 1e-8     # stop when ∥proj grad∥ < 1e-8
          ))
  }, error = function(e) {
    warning(sprintf(
      "fit_one[%s]: L-BFGS-B error: %s
  falling back to SANN params",
      distn, e$message
    ))
    # fabricate a "converged" result
    list(par = s0$par, value = NA, convergence = 1L)
  })
  conv2 <- f0$convergence
  if (conv2 != 0) {
    warning(sprintf(
      "fit_one[%s]: L-BFGS-B did not converge (code=%d); using SANN result",
      distn, conv2
    ))
  }
  
  # choose final parameters
  par_used <- if (conv2 == 0) f0$par else s0$par
  names(par_used) <- par_names
  
  # assemble output
  if (distn == "gamma") {
    list(
      dist  = "gamma",
      shape = par_used["shape"],
      rate  = par_used["rate"],
      conv1 = conv1,
      conv2 = conv2
    )
  } else {
    list(
      dist  = "beta",
      alpha = par_used["alpha"],
      beta  = par_used["beta"],
      conv1 = conv1,
      conv2 = conv2
    )
  }
}

#––– 3) Build mixture functions –––

make_mixture <- function(comp_list) {
  k <- length(comp_list)
  
  comp_means <- map_dbl(comp_list, function(cmp) {
    switch(
      cmp$dist,
      gamma = cmp$shape / cmp$rate,
      beta  = cmp$alpha / (cmp$alpha + cmp$beta),
      point = cmp$value
    )
  })
  
  is_pt <- map_lgl(comp_list, ~ .x$dist == "point")
  pts   <- comp_list[is_pt]
  vals  <- sort(unique(map_dbl(pts, "value")))
  
  if (length(vals)) {
    weights <- sapply(vals, function(v) {
      sum(map_dbl(pts, "value") == v) / k
    })
    cum_w   <- cumsum(weights)
  } else {
    weights <- cum_w <- numeric(0)
  }
  
  p_mix <- function(x) {
    sapply(x, function(xi) {
      cont <- sum(map_dbl(comp_list, function(cmp) {
        switch(
          cmp$dist,
          gamma = pgamma(xi, cmp$shape, cmp$rate),
          beta  = pbeta(xi, cmp$alpha, cmp$beta),
          point = 0
        )
      })) / k
      
      disc <- if (length(vals)) {
        sum(weights[vals <= xi])
      } else {
        0
      }
      cont + disc
    })
  }
  
  q_mix <- function(ps) {
    sapply(ps, function(pi) {
      if (length(vals)) {
        hit <- which(pi <= cum_w)[1]
        if (!is.na(hit)) return(vals[hit])
      }
      
      lo <- if (length(vals)) max(vals) else 0
      hi <- max(c(
        map_dbl(comp_list, function(cmp) {
          switch(
            cmp$dist,
            gamma = cmp$shape / cmp$rate +
              3 * sqrt(cmp$shape/(cmp$rate^2)),
            beta  = 1,
            point = cmp$value
          )
        }), lo + 1
      ))
      
      fl <- p_mix(lo) - pi
      fu <- p_mix(hi) - pi
      
      if (fl >= 0) return(lo)
      if (fu <= 0) return(hi)
      
      uniroot(
        function(x) p_mix(x) - pi,
        lower = lo, upper = hi
      )$root
    })
  }
  
  d_mix <- function(x) {
    sapply(x, function(xi) {
      sum(map_dbl(comp_list, function(cmp) {
        switch(
          cmp$dist,
          gamma = dgamma(xi, cmp$shape, cmp$rate),
          beta  = dbeta(xi, cmp$alpha, cmp$beta),
          point = 0
        )
      })) / k
    })
  }
  
  r_mix <- function(n) {
    idx <- sample.int(k, n, TRUE)
    map_dbl(idx, function(i) {
      cmp <- comp_list[[i]]
      switch(
        cmp$dist,
        gamma = rgamma(1, cmp$shape, cmp$rate),
        beta  = rbeta(1, cmp$alpha, cmp$beta),
        point = cmp$value
      )
    })
  }
  
  list(
    density   = d_mix,
    cdf       = p_mix,
    quantile  = q_mix,
    random    = r_mix,
    comp_means= comp_means
  )
}

#––– 4) Fit single distribution to mixture (with convergence checks) –––
fit_to_mixture <- function(mix, distn, maxit = default_maxit) {
  mix_q <- mix$quantile(probs)
  mix_m <- mean(mix$comp_means)
  distn <- tolower(distn)
  
  # pick objective + parameter names
  obj       <- if (distn == "gamma") obj_gamma else obj_beta
  par_names <- if (distn == "gamma") c("shape","rate") else c("alpha","beta")
  
  # 1) Simulated annealing rough fit
  par0  <- runif(2, 0.1, 5)
  s0    <- optim(par0, obj, q_obs = mix_q, mean_obs = mix_m,
                 method = "SANN",
                 control = list(maxit = maxit, reltol = 1e-10))
  conv1 <- s0$convergence
  if (conv1 != 0) {
    warning(sprintf("fit_to_mixture[%s]: SANN failed (code=%d)",
                    distn, conv1))
  }
  
  # 2) Refine with L-BFGS-B
  f0 <- tryCatch({
    optim(s0$par, obj, q_obs = mix_q, mean_obs = mix_m,
          method = "L-BFGS-B",
          lower = c(1e-8, 1e-8),
          control = list(maxit = maxit, factr = 1e6, pgtol = 1e-8))
  }, error = function(e) {
    warning(sprintf("fit_to_mixture[%s]: L-BFGS-B error: %s; using SANN params",
                    distn, e$message))
    list(par = s0$par, convergence = 1L)
  })
  conv2 <- f0$convergence
  if (conv2 != 0) {
    warning(sprintf("fit_to_mixture[%s]: L-BFGS-B did not converge (code=%d)",
                    distn, conv2))
  }
  
  # choose parameters
  par_used <- if (conv2 == 0) f0$par else s0$par
  names(par_used) <- par_names
  
  # assemble return, including convergence codes
  if (distn == "gamma") {
    shape <- par_used["shape"]
    rate  <- par_used["rate"]
    m     <- shape / rate
    se    <- sqrt(shape / (rate^2))
    list(dist  = "gamma",
         alpha = shape,
         beta  = rate,
         mean  = m,
         SE    = se,
         conv1 = conv1,
         conv2 = conv2)
  } else {
    a  <- par_used["alpha"]
    b  <- par_used["beta"]
    m  <- a / (a + b)
    se <- sqrt((a * b) / (((a + b)^2) * (a + b + 1)))
    list(dist  = "beta",
         alpha = a,
         beta  = b,
         mean  = m,
         SE    = se,
         conv1 = conv1,
         conv2 = conv2)
  }
}

#––– 5) UI Definition –––

ui <- fluidPage(
  titlePanel("Mixture Fitting App"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload .xlsx", accept = ".xlsx"),
      numericInput("seed", "Seed", value = default_seed),
      numericInput("n_boot", "Boot reps", value = default_n_boot, min = 1),
      actionButton("run", "Run"),
      br(),
      actionButton("prev_param", "Prev"),
      actionButton("next_param", "Next"),
      selectInput("param", "Param", choices = NULL, width = "100%"),
      numericInput("maxit", "Max iterations", value = default_maxit, min = 1),
      downloadButton("download","Download .xlsx"),
      conditionalPanel(
        condition = "input.main_tabs == 'Validation'",
        uiOutput("val_series_ui")
      ),
      conditionalPanel(
        condition = "input.main_tabs == 'Quantiles'",
        uiOutput("quant_series_ui")
      )
    ),
    mainPanel(
      tabsetPanel(
        id = "main_tabs",
        tabPanel("Median",    tableOutput("medTable")),
        tabPanel("Bootstrap", DTOutput("bootTable")),
        tabPanel("Histogram", plotOutput("histPlots")),
        tabPanel(
          "Validation",
          plotlyOutput("valPlots"),
          br(),
          DTOutput("valTable")
        ),
        tabPanel(
          "Quantiles",
          plotlyOutput("quantPlots")
        )
      )
    )
  )
)

#––– 6) Server Logic –––
generate_palette <- function(experts) {
  expert_colors <- c(
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
  )
  pal_ex <- setNames(
    rep(expert_colors, length.out = length(experts)),
    experts
  )
  c(pal_ex,
    Mixture    = "grey70",
    `Final Fit` = "black"
  )
}

server <- function(input, output, session) {
  vals <- reactiveValues(
    med = NULL,
    boot = NULL,
    expert_input = NULL,
    param_fits = NULL,
    mixtures = NULL
  )
  
  # Run estimation
  observeEvent(input$run, {
    req(input$file)
    set.seed(input$seed)
    
    withProgress(message = "Running estimation...", value = 0, {
      # 10%: Load & preprocess data
      incProgress(0.1, detail = "Loading data...")
      df <- read_excel(input$file$datapath) %>%
        mutate(
          idx = rep(1:(n()/5), each = 5),
          across(SR:BO, as.numeric)
        ) %>%
        group_by(idx) %>%
        mutate(
          Parameter   = first(Parameter),
          Parameter2  = first(`...2`),
          Distributie = first(Distributie)
        ) %>%
        ungroup() %>%
        mutate(
          Parameter = paste(Parameter, Parameter2),
          across(SR:BO, ~ if_else(
            tolower(Distributie) %in% c("beta","beta?"), .x/100, .x
          ))
        )
      
      # 20%: Fit individual experts and store stats
      incProgress(0.2, detail = "Fitting expert blocks...")
      long <- df %>%
        mutate(rid = rep(1:(n()/5), each = 5)) %>%
        pivot_longer(SR:BO, names_to = "comp", values_to = "val") %>%
        group_by(Parameter, Distributie, comp, rid) %>%
        summarise(vals = list(val), .groups = "drop_last") %>%
        mutate(stats = map(vals, ~ set_names(.x, c("p1","p99","mean","p25","p75")))) %>%
        ungroup()
      
      # store expert inputs
      vals$expert_input <- long %>%
        select(Parameter, comp, stats) %>%
        unnest_wider(stats) %>%
        select(Parameter, Expert = comp, p1, p25, p75, p99, mean) %>%
        split(.$Parameter)
      
      # 30%: Build and store fits
      fits_raw <- long %>%
        rowwise() %>%
        mutate(fit = list(fit_one(stats, Distributie, maxit = input$maxit))) %>%
        ungroup()
      
      # 1) Expert‐fits table (for Excel)
      #    this unnests fit and captures alpha/beta/conv1/conv2 per expert block
      vals$expert_fits <- fits_raw %>%
        unnest_wider(fit)
      
      # 2) Mixture component list (pf) for later plotting
      pf <- fits_raw %>%
        group_by(Parameter) %>%
        summarise(
          components = list(set_names(fit, comp)),
          .groups    = "drop"
        ) %>%
        deframe()
      
      vals$param_fits <- pf
      
      # 40%: Build mixtures
      incProgress(0.4, detail = "Building mixtures...")
      mixtures <- map(pf, make_mixture)
      vals$mixtures <- mixtures
      
      # 50-90%: Bootstrap
      incProgress(0.5, detail = "Bootstrapping...")
      boot_list <- vector("list", input$n_boot)
      for (b in seq_len(input$n_boot)) {
        incProgress(1/input$n_boot, detail = paste0("Bootstrap ", b, " of ", input$n_boot))
        fb <- map2(mixtures, distinct(long, Parameter, Distributie)$Distributie, ~ fit_to_mixture(.x, .y, maxit = input$maxit))
        boot_list[[b]] <- map2_dfr(
          names(fb), fb,
          ~ tibble(
            Parameter = .x,
            dist      = .y$dist,
            alpha     = .y$alpha,
            beta      = .y$beta,
            mean      = .y$mean,
            SE        = .y$SE,
            conv1     = .y$conv1,
            conv2     = .y$conv2
          )
        ) %>% mutate(rep = b)
        
      }
      boot <- bind_rows(boot_list)
      vals$boot <- boot
      
      # 95%: Summarise medians
      incProgress(0.95, detail = "Summarising...")
      med <- boot %>%
        group_by(Parameter) %>%
        summarise(
          dist  = first(dist),
          alpha = median(alpha, na.rm = TRUE),
          beta  = median(beta,  na.rm = TRUE),
          mean  = median(mean,  na.rm = TRUE),
          .groups = "drop"
        )
      vals$med <- med
      
      # Update UI
      updateSelectInput(
        session, "param",
        choices  = vals$med$Parameter,
        selected = vals$med$Parameter[1]
      )
    })
  })

# Navigation
observeEvent(input$prev_param, {
  req(vals$med)
  sel <- vals$med$Parameter
  i   <- match(input$param, sel)
  updateSelectInput(session, "param",
                    selected = sel[max(1, i-1)]
  )
})

observeEvent(input$next_param, {
  req(vals$med)
  sel <- vals$med$Parameter
  i   <- match(input$param, sel)
  updateSelectInput(session, "param",
                    selected = sel[min(length(sel), i+1)]
  )
})

# Toggle UIs
output$val_series_ui <- renderUI({
  req(vals$param_fits, input$param)
  ex <- names(vals$param_fits[[input$param]])
  checkboxGroupInput(
    "val_items", "Choose series:",
    choices = c(ex, "Mixture", "Final Fit"),
    selected = c(ex, "Mixture", "Final Fit")
  )
})

output$quant_series_ui <- renderUI({
  req(vals$param_fits, input$param)
  ex <- names(vals$param_fits[[input$param]])
  checkboxGroupInput(
    "quant_items", "Choose series:",
    choices = c(ex, "Mixture", "Final Fit"),
    selected = c(ex, "Mixture", "Final Fit")
  )
})

# Output: Median
output$medTable <- renderTable({
  req(vals$med)
  vals$med
})

# Output: Bootstrap
output$bootTable <- renderDT({
  req(vals$boot, input$param)
  datatable(
    vals$boot %>% filter(Parameter == input$param),
    rownames = FALSE
  )
})

# Output: Histogram
output$histPlots <- renderPlot({
  req(vals$boot, input$param)
  hist_df <- vals$boot %>%
    filter(Parameter == input$param) %>%
    pivot_longer(
      c(alpha, beta),
      names_to  = "param",
      values_to = "val"
    )
  ggplot(hist_df, aes(x = val)) +
    geom_histogram(bins = 30) +
    facet_wrap(~ param, scales = "free") +
    labs(
      x     = "Value",
      y     = "Count",
      title = paste("Bootstrap –", input$param)
    ) +
    theme_minimal()
})

# Output: Validation plot
output$valPlots <- renderPlotly({
  req(
    vals$expert_input,
    vals$param_fits,
    vals$mixtures,
    input$param,
    input$val_items
  )
  
  expert_vals <- unlist(
    vals$expert_input[[input$param]][c(
      'p1','p25','p75','p99','mean'
    )]
  )
  max_val <- max(expert_vals, na.rm = TRUE) * 1.1
  
  ex_df <- map2_dfr(
    vals$param_fits[[input$param]],
    names(vals$param_fits[[input$param]]),
    function(cmp, nm) {
      xs <- seq(
        0,
        if (cmp$dist == 'beta') 1 else max_val,
        length.out = 200
      )
      ys <- switch(
        cmp$dist,
        beta  = dbeta(xs, cmp$alpha, cmp$beta),
        gamma = dgamma(xs, cmp$shape, cmp$rate)
      )
      tibble(
        x      = xs,
        y      = ys,
        series = nm
      )
    }
  )
  
  xs_mix <- seq(0, max_val, length.out = 300)
  mix_df <- tibble(
    x      = xs_mix,
    y      = vals$mixtures[[input$param]]$density(xs_mix),
    series = 'Mixture'
  )
  
  fr <- vals$med %>% filter(Parameter == input$param)
  fin_df <- tibble(
    x      = xs_mix,
    y      = switch(
      fr$dist,
      beta  = dbeta(xs_mix, fr$alpha, fr$beta),
      gamma = dgamma(xs_mix, fr$alpha, fr$beta)
    ),
    series = 'Final Fit'
  )
  
  all_df <- bind_rows(ex_df, mix_df, fin_df) %>%
    filter(series %in% input$val_items)
  
  # Manual palette
  ex <- names(vals$param_fits[[input$param]])
  pal <- generate_palette(ex)

  #pal <- c(
  #  setNames(
  #    rep('steelblue', length(unique(exp_df$series))),
  #    unique(exp_df$series)
  #  ),
  #  Mixture   = 'grey70',
  #  Final Fit = 'black'
  #)
  
  p <- ggplot(
    all_df,
    aes(x = x, y = y, group = series, color = series)
  ) +
    geom_line(linewidth = 1, alpha = 0.8) +
    scale_color_manual(values = pal) +
    labs(
      x     = 'Value',
      y     = 'Density',
      title = paste('Validation –', input$param)
    ) +
    theme_minimal()
  
  ggplotly(p) %>%
    layout(legend = list(itemclick = 'toggle'))
})

# Output: Expert input table
output$valTable <- renderDT({
  req(vals$expert_input, input$param)
  datatable(
    vals$expert_input[[input$param]] %>% select(-Parameter),
    rownames = FALSE,
    options  = list(pageLength = 5)
  )
})

# Output: Quantiles plot
output$quantPlots <- renderPlotly({
  req(
    vals$expert_input,
    vals$param_fits,
    input$param,
    input$quant_items
  )
  
  raw <- vals$expert_input[[input$param]]
  exp_df <- raw %>%
    pivot_longer(
      c(p1, p25, mean, p75, p99),
      names_to  = 'prob',
      values_to = 'value'
    ) %>%
    rename(series = Expert)
  
  mix_q <- vals$mixtures[[input$param]]$quantile(probs)
  mix_m <- median(raw$mean, na.rm = TRUE)
  mix_df <- tibble(
    series = 'Mixture',
    prob   = names(mix_q),
    value  = mix_q
  ) %>%
    bind_rows(
      tibble(series='Mixture', prob='mean', value=mix_m)
    )
  
  fr      <- vals$med %>% filter(Parameter == input$param)
  final_q <- switch(
    fr$dist,
    beta  = qbeta(probs, fr$alpha, fr$beta),
    gamma = qgamma(probs, fr$alpha, fr$beta)
  )
  fin_df  <- tibble(
    series = 'Final Fit',
    prob   = names(final_q),
    value  = final_q
  ) %>%
    bind_rows(
      tibble(series='Final Fit', prob='mean', value=fr$mean)
    )
  
  all_q <- bind_rows(exp_df, mix_df, fin_df) %>%
    mutate(
      prob = factor(prob, levels = c('p1','p25','mean','p75','p99'))
    ) %>%
    filter(series %in% input$quant_items)
  
  # Manual palette
  ex <- names(vals$param_fits[[input$param]])
  pal <- generate_palette(ex)
  
  p2 <- ggplot(
    all_q,
    aes(x = prob, y = value, group = series, color = series)
  ) +
    geom_point(size=2, alpha = 0.8) +
    geom_line(linewidth=1, alpha = 0.8) +
    scale_color_manual(values=pal) +
    labs(
      x     = 'Quantile/Mean',
      y     = 'Value',
      title = paste('Quantiles –', input$param)
    ) +
    theme_minimal()
  
  ggplotly(p2)
})

# Download handler
output$download <- downloadHandler(
  filename = function() 'results.xlsx',
  content = function(file) {
    wb <- createWorkbook()
    
    # Prepare Median sheet: set alpha=shape, beta=scale for gamma (1/rate), keep beta for beta dist, then mean and sd
    med_df <- vals$med %>%
      rename(alpha_orig = alpha, beta_orig = beta) %>%
      mutate(
        shape = alpha_orig,
        beta  = case_when(
          dist == 'gamma' ~ 1/beta_orig,
          dist == 'beta'  ~ beta_orig,
          TRUE            ~ NA_real_
        ),
        mean = mean,
        sd   = case_when(
          dist == 'gamma' ~ beta * sqrt(shape), # sd = scale * sqrt(shape)
          dist == 'beta'  ~ sqrt((alpha_orig * beta_orig) / (((alpha_orig + beta_orig)^2) * (alpha_orig + beta_orig + 1))),
          TRUE            ~ NA_real_
        )
      ) %>%
      select(Parameter, dist, shape, beta, mean, sd)
    
    addWorksheet(wb, 'Median')
    writeData(wb, 'Median', med_df)
    
    # Add Bootstrap sheet
    addWorksheet(wb, 'Bootstrap')
    writeData(wb, 'Bootstrap', vals$boot)
    
    # Add Expert‐fits sheet
    addWorksheet(wb, "ExpertFits")
    
    # Build an export table from vals$expert_fits
    fits_export <- vals$expert_fits %>%
      # vals$expert_fits came from unnest_wider(fit), so it has:
      #  • dist (“gamma” / “beta” / “point”)
      #  • shape + rate for gamma
      #  • alpha + beta for beta
      #  • value        for point
      #  • conv1, conv2
      mutate(
        alpha = case_when(
          dist == "beta"  ~ alpha,
          dist == "gamma" ~ shape,
          dist == "point" ~ value
        ),
        beta  = case_when(
          dist == "beta"  ~ beta,
          dist == "gamma" ~ 1/rate,    # convert rate→scale
          dist == "point" ~ NA_real_
        ),
        # zero‐fill convergence codes for point‐masses
        conv1 = if_else(dist == "point", 0L, conv1),
        conv2 = if_else(dist == "point", 0L, conv2)
      ) %>%
      select(
        Parameter,
        Expert = comp,
        dist,
        alpha,
        beta,
        conv1,
        conv2
      )
    
    writeData(wb, "ExpertFits", fits_export)
    
    saveWorkbook(wb, file, overwrite = TRUE)
  }
)
}

#––– Run the app –––

shinyApp(ui, server)

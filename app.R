library(shiny)

# =====================================================================
# The Death Correction: protecting relational duration in an open
# population. An interactive companion to the SISMID "Network Modeling
# for Epidemics" course (Module 7).
#
# In a closed population an edge ends only one way: the two partners
# break up, at hazard 1/duration per step. Open the population and an
# edge gains a second way to end: a partner departs (dies), adding a
# hazard of about 2 * d.rate. These are competing risks. The death
# correction lowers the breakup hazard just enough that the two together
# sum back to 1/duration, restoring the realized relational duration.
#
# The analytic core is exact and instant; the simulation panel shows the
# same result with run-to-run variability. Plot and text builders are
# pure functions so they can be tested outside Shiny.
# =====================================================================


# ---- Analytic core --------------------------------------------------

# Per-step probability that at least one of the two partners departs.
hazard_dep <- function(d_rate) 1 - (1 - d_rate)^2

# Per-step breakup hazard. Uncorrected it is simply 1/duration. Corrected
# it is lowered so breakup and departure together return the realized
# ending hazard to 1/duration. Clamped at 0 when departure alone already
# ends edges faster than the target (the correction ceiling).
breakup_hazard <- function(duration, h_d, corrected) {
  if (!corrected) return(1 / duration)
  hb <- 1 - (1 - 1 / duration) / (1 - h_d)
  max(hb, 0)
}

# Realized mean relational duration under the two competing hazards.
realized_duration <- function(duration, h_d, corrected) {
  hb <- breakup_hazard(duration, h_d, corrected)
  1 / (1 - (1 - hb) * (1 - h_d))
}

# Dissolution (persistence) coefficient for offset(edges): logit(1 - hb).
diss_coef <- function(hb) log((1 - hb) / hb)

# Share of edge endings attributable to departure vs breakup.
departure_share <- function(hb, h_d) {
  lam_b <- -log(1 - hb)
  lam_d <- -log(1 - h_d)
  if (lam_b + lam_d == 0) return(0)
  lam_d / (lam_b + lam_d)
}

# Departure rate at which the correction can no longer hold the target.
ceiling_drate <- function(duration) 1 - sqrt(1 - 1 / duration)

# Stationary simulation of an ongoing edge population. Every ended edge is
# replaced by a fresh one, and we track the mean age of active edges each
# step. For geometric durations the mean age of extant ties equals the
# mean total duration, so this converges to realized_duration.
sim_mean_age <- function(duration, h_d, corrected,
                         nsteps = 250, n_edges = 400, reps = 5) {
  hb <- breakup_hazard(duration, h_d, corrected)
  p_end <- 1 - (1 - hb) * (1 - h_d)
  m <- matrix(NA_real_, nsteps, reps)
  for (r in seq_len(reps)) {
    ages <- rep(0, n_edges)
    for (t in seq_len(nsteps)) {
      ended <- runif(n_edges) < p_end
      ages[ended] <- 0
      ages <- ages + 1
      m[t, r] <- mean(ages)
    }
  }
  m
}


# ---- Palette --------------------------------------------------------
col_corr   <- "#1a6b54"   # teal: with correction
col_none   <- "#c0392b"   # red: without correction
col_steel  <- "#4682B4"   # steel blue: departure component
col_target <- "#333333"


# ---- Plot and text builders (pure; tested outside Shiny) ------------

plot_duration_curve <- function(duration, d_rate, drate_max = 0.01) {
  D <- duration
  dseq <- seq(0, drate_max, length.out = 240)
  hd <- hazard_dep(dseq)
  dur_c <- vapply(hd, function(x) realized_duration(D, x, TRUE),  numeric(1))
  dur_n <- vapply(hd, function(x) realized_duration(D, x, FALSE), numeric(1))
  hnow <- hazard_dep(d_rate)

  par(family = "sans", mar = c(4.5, 4.5, 2.5, 1), bg = "#ffffff")
  plot(dseq, dur_c, type = "n", ylim = c(0, D * 1.12),
       xlab = "Departure rate (d.rate)", ylab = "Realized mean duration",
       main = "Realized duration as departures rise",
       col.main = col_corr, font.main = 2, cex.main = 1.25,
       col.lab = "#555", cex.lab = 1.05, col.axis = "#777")
  abline(h = D, lty = 2, lwd = 2, col = col_target)
  lines(dseq, dur_n, lwd = 3, col = col_none)
  lines(dseq, dur_c, lwd = 3, col = col_corr)
  abline(v = d_rate, lty = 3, lwd = 2, col = "#999999")
  points(d_rate, realized_duration(D, hnow, FALSE), pch = 19, cex = 1.4, col = col_none)
  points(d_rate, realized_duration(D, hnow, TRUE),  pch = 19, cex = 1.4, col = col_corr)
  legend("bottomleft", bty = "n", lwd = c(3, 3, 2), lty = c(1, 1, 2),
         col = c(col_corr, col_none, col_target),
         legend = c("With correction", "Without correction", paste0("Target = ", D)),
         text.col = "#555")
}

plot_hazard_decomp <- function(duration, d_rate) {
  D <- duration
  hd <- hazard_dep(d_rate)
  hb_n <- breakup_hazard(D, hd, FALSE)
  hb_c <- breakup_hazard(D, hd, TRUE)
  mat <- rbind(breakup = c(hb_n, hb_c), departure = c(hd, hd))

  par(family = "sans", mar = c(4, 4.5, 2.5, 1), bg = "#ffffff")
  bp <- barplot(mat, names.arg = c("Without correction", "With correction"),
                col = c(col_none, col_steel), border = "#ffffff",
                ylim = c(0, max(colSums(mat)) * 1.25),
                ylab = "Per-step ending hazard",
                main = "Where edge endings come from",
                col.main = col_corr, font.main = 2, cex.main = 1.25,
                col.lab = "#555", cex.lab = 1.05, col.axis = "#777",
                cex.names = 1.05)
  abline(h = 1 / D, lty = 2, lwd = 2, col = col_target)
  text(bp[2], 1 / D, labels = "target hazard = 1 / duration",
       pos = 3, offset = 0.4, col = col_target, cex = 0.95)
  legend("topright", bty = "n", fill = c(col_none, col_steel), border = "#ffffff",
         legend = c("breakup", "departure"), text.col = "#555")
}

plot_sim <- function(sim_corr, sim_none, duration) {
  D <- duration
  par(family = "sans", mar = c(4.5, 4.5, 2.5, 1), bg = "#ffffff")
  matplot(sim_corr, type = "l", lty = 1, lwd = 1.5, col = adjustcolor(col_corr, 0.6),
          ylim = c(0, D * 1.2), xlab = "Time step",
          ylab = "Mean age of active partnerships",
          main = "Simulated realized duration, with and without correction",
          col.main = col_corr, font.main = 2, cex.main = 1.2,
          col.lab = "#555", cex.lab = 1.05, col.axis = "#777")
  matlines(sim_none, lty = 1, lwd = 1.5, col = adjustcolor(col_none, 0.6))
  lines(rowMeans(sim_corr), lwd = 3, col = col_corr)
  lines(rowMeans(sim_none), lwd = 3, col = col_none)
  abline(h = D, lty = 2, lwd = 2, col = col_target)
  legend("bottomright", bty = "n", lwd = c(3, 3, 2), lty = c(1, 1, 2),
         col = c(col_corr, col_none, col_target),
         legend = c("With correction", "Without correction", paste0("Target = ", D)),
         text.col = "#555")
}

# Returns a list of recap items (character or HTML) for the results card.
recap_items <- function(duration, d_rate) {
  D <- duration
  hd <- hazard_dep(d_rate)
  hb_n <- breakup_hazard(D, hd, FALSE)
  hb_c <- breakup_hazard(D, hd, TRUE)
  dur_n <- realized_duration(D, hd, FALSE)
  dur_c <- realized_duration(D, hd, TRUE)
  items <- list(
    sprintf("Departure adds a per-step hazard of %.4f (about 2 x %.4f).", hd, d_rate),
    HTML(sprintf("Realized duration <strong>without</strong> the correction: <strong>%.1f</strong> steps (%.0f%% of the %d target).",
                 dur_n, 100 * dur_n / D, D)),
    HTML(sprintf("Realized duration <strong>with</strong> the correction: <strong>%.1f</strong> steps (on target).",
                 dur_c)),
    sprintf("Dissolution coefficient rises from %.3f (crude) to %.3f (corrected).",
            diss_coef(hb_n), diss_coef(hb_c)),
    sprintf("With the correction, about %.0f%% of endings are departures and %.0f%% are breakups.",
            100 * departure_share(hb_c, hd), 100 * (1 - departure_share(hb_c, hd))),
    HTML(sprintf("Mean degree comes along for the ride: without the correction it sags to about %.0f%% of target; with it, it holds.",
                 100 * dur_n / D))
  )
  if (d_rate > ceiling_drate(D)) {
    items <- c(items, list(HTML(sprintf(
      "<span class='ceiling-note'>Above the ceiling: at this duration the correction can hold departures only up to d.rate = %.4f. Beyond that, departures alone end edges faster than the target, so even zero breakup cannot reach %d.</span>",
      ceiling_drate(D), D))))
  }
  items
}


# ---- UI -------------------------------------------------------------
ui <- fluidPage(

  tags$head(
    tags$link(
      href = "https://fonts.googleapis.com/css2?family=Open+Sans:wght@400;600;700&display=swap",
      rel = "stylesheet"
    ),
    tags$style(HTML("
      html { font-size: 16px; }
      body {
        font-family: 'Open Sans', -apple-system, BlinkMacSystemFont, sans-serif;
        font-size: 1rem; color: #444; background-color: #f0f5f9;
        -webkit-font-smoothing: antialiased; -moz-osx-font-smoothing: grayscale;
      }
      .container-fluid { max-width: 1200px; padding: 0 2rem; }
      .app-header {
        background: linear-gradient(135deg, #2f5d47, #1a6b54);
        color: #ffffff; padding: 2rem 2.5rem; margin: -15px -15px 1.5rem -15px;
        border-radius: 0 0 8px 8px;
      }
      .app-header h2 { font-weight: 700; font-size: 2rem; margin: 0 0 0.4rem 0; }
      .app-header p { font-size: 1.15rem; margin: 0; opacity: 0.9; }
      .card {
        background: #ffffff; border: 1px solid #dee2e6; border-radius: 8px;
        padding: 1.5rem; box-shadow: 0 2px 6px rgba(0,0,0,0.05); margin-bottom: 1.25rem;
      }
      .card h4 {
        font-weight: 600; font-size: 1.25rem; color: #1a6b54; margin-top: 0;
        margin-bottom: 0.75rem; padding-bottom: 0.5rem; border-bottom: 2px solid #e8eef3;
      }
      .card p, .card li { font-size: 1.05rem; line-height: 1.65; color: #555; }
      .card strong { color: #1a6b54; }
      .well {
        background: #ffffff !important; border: 1px solid #dee2e6 !important;
        border-radius: 8px !important; box-shadow: 0 2px 6px rgba(0,0,0,0.05) !important;
        padding: 1.5rem !important;
      }
      .control-label { font-weight: 600; font-size: 1.0rem; color: #2f5d47; margin-bottom: 0.25rem; }
      .irs--shiny .irs-bar { background: #1a6b54; border-top: 1px solid #2f5d47; border-bottom: 1px solid #2f5d47; }
      .irs--shiny .irs-single, .irs--shiny .irs-from, .irs--shiny .irs-to { background-color: #1a6b54; }
      .irs--shiny .irs-handle { border: 1px solid #1a6b54; }
      .btn-primary { background-color: #1a6b54; border-color: #1a6b54; font-weight: 600; }
      .btn-primary:hover { background-color: #14543f; border-color: #14543f; }
      .plot-card, .results-card {
        background: #ffffff; border: 1px solid #dee2e6; border-radius: 8px;
        padding: 1.25rem; box-shadow: 0 2px 6px rgba(0,0,0,0.05); margin-bottom: 1.25rem;
      }
      .results-card .results-title {
        font-weight: 600; font-size: 1.25rem; color: #1a6b54; margin-top: 0;
        margin-bottom: 0.75rem; padding-bottom: 0.5rem; border-bottom: 2px solid #e8eef3;
      }
      .results-card li {
        font-size: 1.05rem; line-height: 1.6; color: #555; padding: 0.4rem 0;
        border-bottom: 1px solid #f0f0f0; list-style: none; padding-left: 0.2rem;
      }
      .results-card li:last-child { border-bottom: none; }
      .results-card ul { padding-left: 0; margin: 0; }
      .ceiling-note { color: #b34700; font-weight: 600; }
      .app-footer { text-align: center; color: #777; font-size: 0.95rem; padding: 1rem 0 1.5rem 0; }
      .app-footer a { color: #1a6b54; text-decoration: none; }
      .app-footer a:hover { text-decoration: underline; }
    "))
  ),

  div(class = "app-header",
    h2("The Death Correction"),
    p("Why an open population needs an adjustment to hold relational duration on target")
  ),

  fluidRow(
    column(6,
      div(class = "card",
        h4("Two ways a partnership ends"),
        p("In a", strong("closed"), "population an edge ends only one way: the two
           partners break up, at a hazard of", strong("1 / duration"), "per step.
           Open the population and an edge gains a", strong("second"), "way to end:
           a partner leaves the population, most often by dying. Either one is
           enough, so departures add a hazard of about", strong("2 x d.rate."),
          "These are", strong("competing risks."))
      )
    ),
    column(6,
      div(class = "card",
        h4("What the correction does"),
        p("The dissolution model was fit knowing nothing about departure. Layer
           departures on top and the two hazards stack, so realized durations come
           out short. The", strong("death correction"), "lowers the breakup hazard
           just enough that the two together sum back to 1 / duration, restoring the
           duration. That is why it makes edges", strong("more"), "persistent."),
        p("Drag the", strong("departure rate"), "and watch the uncorrected duration
           fall away from target while the corrected one holds.")
      )
    )
  ),

  sidebarLayout(
    sidebarPanel(width = 3,
      tags$h4(style = "font-weight:600; font-size:1.25rem; color:#1a6b54; margin-top:0;
                       margin-bottom:0.75rem; padding-bottom:0.5rem; border-bottom:2px solid #e8eef3;",
              "Parameters"),
      sliderInput("duration", "Target mean duration (steps)",
                  min = 10, max = 100, value = 40, step = 1),
      sliderInput("d_rate", "Departure rate (per partner, per step)",
                  min = 0, max = 0.01, value = 0.005, step = 0.0005),
      helpText("The course uses 0.001 (realistic) and 0.005 (amplified for a visible effect)."),
      tags$hr(),
      actionButton("run_sim", "Run simulation", class = "btn-primary"),
      helpText("Simulates an ongoing edge population to show the same result with
                run-to-run variability.")
    ),
    mainPanel(width = 9,
      div(class = "plot-card", plotOutput("dur_curve", height = "360px")),
      div(class = "results-card",
        div(class = "results-title", "At the current departure rate"),
        htmlOutput("recap")
      ),
      div(class = "plot-card", plotOutput("hazard_bar", height = "340px")),
      div(class = "plot-card",
        tags$p(style = "color:#777; font-size:0.95rem; margin:0 0 0.5rem 0.25rem;",
               "Simulation: mean age of active partnerships over time (press Run)."),
        plotOutput("sim_plot", height = "380px")
      )
    )
  ),

  div(class = "app-footer",
    p("Part of the", tags$a(href = "https://epimodel.org", "EpiModel"),
      "ecosystem • Powered by",
      tags$a(href = "https://posit-dev.github.io/r-shinylive/", "shinylive"))
  )
)


# ---- Server ---------------------------------------------------------
server <- function(input, output, session) {

  output$dur_curve <- renderPlot(plot_duration_curve(input$duration, input$d_rate))
  output$hazard_bar <- renderPlot(plot_hazard_decomp(input$duration, input$d_rate))
  output$recap <- renderUI(
    tags$ul(lapply(recap_items(input$duration, input$d_rate), function(x) tags$li(x)))
  )

  sim_data <- eventReactive(input$run_sim, {
    hd <- hazard_dep(input$d_rate)
    list(D = input$duration,
         corr = sim_mean_age(input$duration, hd, TRUE),
         none = sim_mean_age(input$duration, hd, FALSE))
  })

  output$sim_plot <- renderPlot({
    d <- sim_data()
    plot_sim(d$corr, d$none, d$D)
  })
}

shinyApp(ui = ui, server = server)

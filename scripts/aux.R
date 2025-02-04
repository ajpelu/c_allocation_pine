# Auxiliary scripts with colours,etc 
colours_elev <- c("low" = "#ef7159",
                  "low-Dec" ="#FFCC3F",
                  "medium" = "#1B1F47",
                  "high" = "#53A049")

colours_elev2 <- c("low-Dec" = "#8c510a",
                  "low" ="#fc8d59",
                  "medium" = "#2166ac",
                  "high" = "#72A55A")


colours_temp <- c("tmin" = "#a6bddb", "tmed" = "#67a9cf", "tmax" = "#016c59")


lines_lm <- c("sig" = "solid", "no sig" = "11") # Used this pattern see https://ggplot2.tidyverse.org/reference/aes_linetype_size_shape.html
# lines_lm <- c("sig" = "solid", "no sig" = "dashed")

lines_lm_size <- c("sig" = .8, "no sig" = 1.3)

# colours_sp <- c("pinaster" = "#f9e900", "nigra" = "darkgreen", "halepensis" = "orange", "sylvestris" = "#67a9cf")
colours_sp <- c("halepensis" = "#c1666b",
                "pinaster" = "#ecb42e",
                "nigra" = "#43ba85",
                "sylvestris" = "#006494")


colours_Specie <- c("P. halepensis" = "#c1666b",
                    "P. pinaster" = "#ecb42e",
                    "P. nigra" = "#43ba85",
                    "P. sylvestris" = "#006494")


colours_Specie_dark <- c(
  "P. halepensis" = "#8e4a4e", # Más oscuro que "#c1666b"
  "P. pinaster"   = "#b48d21", # Más oscuro que "#ecb42e"
  "P. nigra"      = "#2e8c63", # Más oscuro que "#43ba85"
  "P. sylvestris" = "#004266"  # Más oscuro que "#006494"
)

# pinaster, "#ffffbf", "beige" or "#f9e900"
# nigra, "#abdda4", "lightgreen"
# halepensis, "orange", "orange"
# sylvestris, "lightblue", "lightblue"

# shape
# shape_elev <- c("low" = 16, "low2" = 18, "medium" = 15, "high" = 17)
shape_elev <- c("low" = 21, "low-Dec" = 25, "medium" = 22, "high" = 24)
shape_elev2 <- c("low" = 19, "low-Dec" = 18, "medium" = 15, "high" = 17)





## Base theme
base_theme <- list(
  theme_bw(),
  theme(text = element_text(family = "Helvetica"),
        axis.title.x = element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold", size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 8))
) 









### Custom functions 


## Custom function to compute drought events ----------------------
# Compare the severe drought events.
# A severe drought event starts when SPEI falls below the threshold of −1.28 (Páscoa et al 2017; Spinoni et al. 2018). A drought event is considered only when SPEI values fall below that threshold for at least two consecutive months. For each drought event, we computed:
#
#   - the duration as the number of consecutive months with the SPEI lower than a certain threshold
# - the severity as the sum of the absolute SPEI values during the drought event
# - the intensity and the Lowest SPEI refer to the mean and lowest value of SPEI, respectively, during the drought event.
#
# We computed the severe drought events (below -1.28) by site and per each index (i.e. SPEI-06 and SPEI-12) since 1960.

droughtIndicators <- function(df, vname, threshold) {
  require(data.table)
  
  # Add a new column indicating if the value is below the threshold and the following month is also below the threshold
  out <- df |>
    mutate(is_drought = ifelse(
      .data[[vname]] < threshold & lead(.data[[vname]], default = .data[[vname]][n()]) < threshold, 1, 0
    )) |>
    mutate(date = lubridate::make_date(year, month))
  
  # Compute the drought duration of the events
  out2 <- out |>
    group_by(index_events = data.table::rleid(is_drought)) |>
    mutate(drought_duration = sum(is_drought)) |>
    as.data.frame()
  
  # Filter events with drought duration > 1
  out3 <- out2 |>
    filter(drought_duration > 1) |>
    as.data.frame()
  
  # Compute several indicators (drought assessments)
  da <- out3 |>
    group_by(index_events) |>
    summarise(
      d_duration = unique(drought_duration),
      d_intensity = mean(.data[[vname]], na.rm = TRUE),
      d_severity = sum(abs(.data[[vname]]), na.rm = TRUE),
      lowest_spei = min(.data[[vname]]),
      month_peak = month[which.min(.data[[vname]])],
      minyear = min(year),
      maxyear = max(year),
      rangeDate = paste(
        lubridate::month(min(date, na.rm = TRUE), label = TRUE), "-",
        (lubridate::month(max(date, na.rm = TRUE), label = TRUE))
      )
    ) |>
    as.data.frame()
  
  return(list(data = out2, drought_events = out3, drought_assessment = da))
}


#### Functions to extract results from MLE
# Define custom functions to extract data from models
MLE_results_format <- function(x, yvar) {
  
  # compute RMSE
  x$source_data$residual <- x$source_data[[yvar]] - x$source_data[["predicted"]]
  rmse <- sqrt(sum(x$source_data$residual^2)/(length(x$source_data$residual)-1))
  
  out <- data.frame(
    max_likeli = x$max_likeli,
    n_params = length(x$best_pars),
    aic_cor = x$aic_corr,
    aic = x$aic,
    R2 = x$R2,
    slope = x$slope,
    RMSE = rmse
  )
  return(out)
}

MLE_results_params <- function(x) {
  
  params <- list(
    value = as.data.frame(x$best_pars),
    support_interval_lower = as.data.frame(x$lower_limits),
    support_interval_upper = as.data.frame(x$upper_limits),
    std_error = as.data.frame(x$std_errs),
    bound_lower = as.data.frame(x$par_lo),
    bound_upper = as.data.frame(x$par_hi)
  )
  
  out <- bind_rows(params, .id = "var") |>
    pivot_longer(cols = -var, names_to = "parameter", values_to = "value") |>
    pivot_wider(names_from = var, values_from = value)
  
  return(out)
}



# Function to compute deltaAIC
dAIC <- function(x) {
  x <- x[!is.na(x)]
  delta.aic <- x - min(x, na.rm = TRUE)
  return(delta.aic)
}

wAIC <- function(delta) {
  delta <- delta[!is.na(delta)]
  aux_w <- exp(-0.5*delta)
  wi <- aux_w / sum(aux_w, na.rm = TRUE)
  return(wi)
}


# Combined function to plot observed vs. predicted or residuals vs. predicted

plot_mle <- function(x, yvar, title = NULL, plot_type = c("observed", "residuals")) {
  plot_type <- match.arg(plot_type)
  
  data <- x$source_data |>
    mutate(residuals = !!sym(yvar)  - predicted) |>
    rename(observed = !!sym(yvar))
  
  p <- ggplot(data, aes(x = predicted)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    ggtitle(title) +
    labs(x = "Predicted")
  
  if (plot_type == "observed") {
    p <- p + geom_point(aes(y = observed)) +
      geom_abline(slope = 1, intercept = 0, color = "red") +
      labs(y = "Observed")
  } else if (plot_type == "residuals") {
    p <- p + geom_point(aes(y = residuals)) +
      geom_hline(yintercept = 0, color = "red") +
      labs(y = "Residuals")
  }
  
  return(p)
}


# plot_mle(x = result_mtabi1, yvar = "abi",  plot_type = "observed")
# plot_mle(x = result_mtabi1, yvar = "abi",  plot_type = "residuals")




######## Geom text 
# Default geom text is 3.88 
# 1 pt = 0.35 mm 
# size = 12 / .pt 





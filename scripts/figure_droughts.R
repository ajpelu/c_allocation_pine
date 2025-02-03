# title: "Drought events. Compute events and plot "

library(tidyverse)
library(patchwork)
library(scico)


source("scripts/aux.R")


# Read spei data
## annual pet
spei_data <- read_csv("data/spei_climate.csv")


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
# -----------------

# Select spei06 y spei12 data
s <- spei_data |>
  dplyr::select(sp_elev, year, month, spei06, spei12) |>
  pivot_longer(c(spei06, spei12), names_to = "spei") |>
  filter(!is.na(value))

spei_events <- s |>
  group_by(spei, sp_elev) |>
  nest() |>
  mutate(drought_results = purrr::map(data, ~ droughtIndicators(., "value", -1.28)$drought_assessment)) |>
  unnest(drought_results) |>
  select(-data)

d <- spei_events |>
  separate(rangeDate, into = c("minmonth", "maxmonth"), remove = FALSE) |>
  unite("start_date", minyear, minmonth, sep = "-", remove = FALSE) |>
  unite("end_date", maxyear, maxmonth, sep = "-", remove = FALSE) |>
  mutate(start_date = as.Date(paste0(start_date, "-01"), format = "%Y-%b-%d")) |>
  mutate(end_date = as.Date(paste0(end_date, "-01"), format = "%Y-%b-%d")) |>
  rowwise() |>
  mutate(middate = (start_date + ((end_date - start_date) / 2)))


#### Drought events SPEI-12
d12 <- d |>
  filter(spei == "spei12") |>
  ungroup() |>
  dplyr::select(sp_elev, d_duration, d_intensity, d_severity, lowest_spei, month_peak, start_date, end_date, middate) |>
  separate(sp_elev, c("sp_code", "elev_code"), sep = "_", remove = FALSE) |>
  mutate(species = paste0("P. ", sp_code))

image_data <- data.frame(
  species = c("P. sylvestris", "P. nigra", "P. pinaster", "P. halepensis"),
  elev_code = c("medium", "medium", "low", "medium"),
  middate = ymd("1970-01-01"),
  img = c("data/images/sylvestris.png", "data/images/nigra.png", 
          "data/images/pinaster.png", "data/images/halepensis.png") 
)

plot_droughts <- d12 |>
  mutate(elev_code = fct_relevel(elev_code, "low-Dec", "low", "medium", "high")) |>
  ggplot(aes(y = elev_code)) +
  geom_point(aes(x = middate, y = elev_code, size = d_severity, color = d_duration)) +
  labs(size = "Drought severity", x = "Year", y = "Elevation") +
  facet_wrap(~ factor(species, levels = c("P. sylvestris", "P. nigra", "P. pinaster", "P. halepensis")), 
             ncol = 1, scales = "free_y") +
  theme_bw() +
  ggimage::geom_image(
    data = image_data,
    aes(x = middate, y = elev_code, image = img), size = 0.8, asp = 0.5
  ) +
  scale_color_scico(
    palette = "vik", direction = 1,
    name = "Drought duration (months)",
    limits = c(2, 20)
  ) +
  scale_x_date(
    limits = c(ymd("1970-01-01"), ymd("2023-01-01")),
    breaks = seq.Date(from = ymd("1970-01-01"), to = ymd("2023-01-01"), by = "5 years"),
    date_labels = "%Y",
    date_minor_breaks = "1 year"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "italic", hjust = 0),
    legend.position = "bottom"
  ) +
  guides(
    color = guide_colorbar(title.hjust = 0.5), #
    size = guide_legend(title.hjust = 0.5)
  )

plot_droughts

ggsave(
  plot_droughts,
  file = "output/plot_drouhgts.jpg",
  dpi = 400,
  width = 7.09 * 1, height = 7.09 * 1
)

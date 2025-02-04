# title: "Figure NPP EVI ABI vs climate"


library(tidyverse)
library(ggh4x)
library(patchwork)
source("scripts/aux.R")

# Read data
## annual pet
annual_pet <- read_csv("data/spei_climate.csv") |> 
  dplyr::select(sp_elev, year, monthly_pet, monthly_tmed, monthly_prec) |> 
  group_by(sp_elev, year) |>
  summarise(pet = sum(monthly_pet), 
            prec = sum(monthly_prec),
            tmed = mean(monthly_tmed, na.rm = TRUE)) |> 
  rowwise() |> 
  mutate(water_balance = prec - pet)

## abi 
abi <- read_csv("data/abi.csv") |> 
  rename(mean = IBT_ag_m2) |> 
  mutate(se = NA, sd = NA, variable = "abi") 

## evi Landsat
evi_landsat <- read_csv("data/iv_landsat.csv") |> 
  filter(iv == "evi") |> 
  dplyr::select(year, sp_code, elev_code, sp_elev, mean, sd, se) |> 
  mutate(variable = "evi_landsat")

npp <- read_csv("data/npp_modis.csv") |> 
  rename(mean = npp) |> 
  mutate(se = NA, sd = NA, variable = "npp") |> 
  dplyr::select(year, sp_code, elev_code, sp_elev, mean, sd, se, variable) 


df_index <- bind_rows(
  abi, evi_landsat, npp) |> 
  mutate(elev_code = fct_relevel(elev_code, "low-Dec","low", "medium", "high")) |> 
  mutate(Specie = paste0("P. ", sp_code)) |> 
  rename(mean_y = mean, y_variable = variable)

df <- df_index |> 
  inner_join(annual_pet) |> 
  pivot_longer(pet:water_balance, values_to = "mean_climate", 
               names_to = "climate_variable")

df_plot <- df |> 
  # filter(year > 1986) |> 
  mutate(y_variable2 = case_when(
    y_variable == "abi" ~ "ABI~(g~C~m^2~year^{-1})",
    y_variable == "evi_landsat" ~ "EVI[Landsat]",
    y_variable == "npp" ~"NPP[MODIS]~(g~C~m^2~year^{-1})")) |> 
  mutate(y_variable2 = fct_relevel(y_variable2, 
                                   "EVI[Landsat]", 
                                   "NPP[MODIS]~(g~C~m^2~year^{-1})",
                                   "ABI~(g~C~m^2~year^{-1})"))

## Linear regressions for plot
l <- df_plot |> 
  filter(climate_variable != "pet") |> 
  group_by(y_variable,climate_variable, sp_code) |> 
  group_modify(
    ~ broom::tidy(lm(mean_y ~ mean_climate, data = .x))
  ) |>
  filter(term != "(Intercept)") |> 
  mutate(p = case_when(
    p.value < 0.001 ~ "<0.001",
    p.value >= 0.001 ~ as.character(round(p.value, 3))
  )) |> 
  mutate(sig = case_when(
    p.value < 0.001 ~ "sig",
    p.value >= 0.001 ~ "no sig"  
  )) 


df_plot <- 
  df_plot |> inner_join(
    l |> dplyr::select(y_variable, climate_variable, sp_code, p, sig)) |> 
  mutate(sp_code = fct_relevel(sp_code, "halepensis","pinaster", "nigra", "sylvestris")) |> 
  mutate(Specie = fct_relevel(Specie, "P. halepensis","P. pinaster", "P. nigra", "P. sylvestris"))


# Get time range periods 
df_plot |> group_by(y_variable) |> 
  summarise(start = min(year),
            end = max(year))


# Plot 
## general_parameters

alpha_points <- 0.3
main_line_width <- 1
main_line_color <- "black"
partial_lines_width <- 1.1
size_points <- 1.7
stroke_points <- 0.1 

label_landsat <- "EVI[Landsat]"
label_npp <- "NPP[MODIS]~(g~C~m^2~year^{-1})"
label_abi <- "ABI~(g~C~m^2~year^{-1})"
label_wb <- "P-PET (mm)"
label_prec <- "Precipitation (mm)"

y_scales <- list(
  scale_y_continuous(limits = c(0, 0.5)),
  scale_y_continuous(limits = c(0,1250)),
  scale_y_continuous(limits = c(0, 700))
)

## Custom themes 
custom_scales <- list(
  scale_colour_manual(
    values = colours_Specie, 
    labels = expression(italic("P. halepensis"), italic("P. pinaster"), italic("P. nigra"), italic("P. sylvestris")),
    name = "Species"),
  scale_fill_manual(
    values = colours_Specie, 
    labels = expression(italic("P. halepensis"), italic("P. pinaster"), italic("P. nigra"), italic("P. sylvestris")),
    name = "Species"),
  scale_shape_manual(values = shape_elev, name = "Elevation")
)

base_theme <- list(
  theme_bw(),
  theme(text = element_text(family = "Helvetica"),
        axis.title.x = element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold", size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 8))
) 

custom_theme_figure <- list(
  theme(
    panel.grid = element_blank(), 
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y = element_text(vjust = 0)
  )
) 


## Letters for label 
letras <- data.frame(
  tipo = rep(c("prec", "temp", "wb"), each = 3),  
  y_variable2 = rep(
    c("EVI[Landsat]", "NPP[MODIS]~(g~C~m^2~year^{-1})", "ABI~(g~C~m^2~year^{-1})"),
    times = 3
  ),
  label = c("a", "b", "c", "d", "e", "f", "f", "g", "h"),  
  x = -Inf,  
  y = Inf    
)


# ----------

plot_prec <- 
  df_plot|> 
  filter(climate_variable == "prec") |> 
  ggplot(aes(x=mean_climate, y = mean_y, colour = Specie)) + 
  geom_point(aes(shape = elev_code, fill = Specie), 
             alpha = alpha_points, size = size_points, stroke = stroke_points) + 
  geom_smooth(aes(linetype = sig), linewidth = partial_lines_width, method = "lm", se = FALSE) +
  scale_linetype_manual(values = lines_lm, guide = "none") + 
  geom_smooth(aes(), 
              method = "nls", 
              formula=y~SSlogis(x, Asym, xmid, scal),
              se =  FALSE, # this is important 
              linewidth = main_line_width, colour = main_line_color) +
  ylab("") + 
  xlab("Precipitation (mm)") +
  custom_scales + 
  base_theme + 
  custom_theme_figure + 
  theme(
    legend.position = "top"
    # legend.box = "vertical"
  ) + 
  scale_x_continuous(limits = c(0, 1750)) +
  facet_wrap(~factor(y_variable2, c("EVI[Landsat]", 
                                    "NPP[MODIS]~(g~C~m^2~year^{-1})",
                                    "ABI~(g~C~m^2~year^{-1})")),
             scales = "free_y", 
             labeller = label_parsed, 
             strip.position = "left") +
  ggh4x::facetted_pos_scales(y = y_scales) +
  geom_text(data = (letras |> filter(tipo == "prec")), 
            aes(x = x, y = y, label = label),
            vjust = 1.2, hjust = -1.2, size = 5, 
            fontface = "bold", inherit.aes = FALSE) +
  guides(
    colour = guide_legend(override.aes = list(size = 4)), 
    shape = guide_legend(override.aes = list(size = 4)), 
  )



plot_tmed <- 
  df_plot |> 
  # filter(year > 1988) |> 
  filter(climate_variable == "tmed") |> 
  ggplot(aes(x=mean_climate, y = mean_y, colour = Specie)) + 
  geom_point(aes(shape = elev_code, fill = Specie), 
             alpha = alpha_points, size = size_points, stroke = stroke_points) +
  geom_smooth(aes(linetype = sig), linewidth = partial_lines_width, method = "lm", se = FALSE) +
  scale_linetype_manual(values = lines_lm, guide = "none") +
  geom_smooth(aes(),
              method = "nls", 
              formula = y ~ SSlogis(x, Asym, xmid, scal),
              se =  FALSE, # this is important 
              linewidth = main_line_width, colour = main_line_color, 
              control = nls.control(maxiter = 10000)) +
  ylab("") + 
  xlab("Annual Mean Temperature (ºC)") +
  custom_scales + 
  base_theme + 
  custom_theme_figure + 
  theme(
    legend.position = "none", 
  ) + 
  scale_x_continuous(limits = c(8, 16)) +
  facet_wrap(~factor(y_variable2, c("EVI[Landsat]", 
                                    "NPP[MODIS]~(g~C~m^2~year^{-1})",
                                    "ABI~(g~C~m^2~year^{-1})")),
             scales = "free_y", 
             labeller = label_parsed, 
             strip.position = "left") +
  facetted_pos_scales(y = y_scales) +
  geom_text(data = (letras |> filter(tipo == "temp")), 
            aes(x = x, y = y, label = label),
            vjust = 1.2, hjust = -1.2, size = 5, 
            fontface = "bold", inherit.aes = FALSE)


plot_tmed 

plot_wb <- 
  df_plot|> 
  filter(climate_variable == "water_balance") |> 
  ggplot(aes(x=mean_climate, y = mean_y, colour = Specie)) + 
  geom_point(aes(shape = elev_code, fill = Specie), alpha = alpha_points, size = size_points, stroke = stroke_points) +
  geom_smooth(aes(linetype = sig), linewidth = partial_lines_width, method = "lm", se = FALSE) +
  scale_linetype_manual(values = lines_lm, guide = "none") +
  geom_smooth(aes(), 
              method = "nls", 
              formula=y~SSlogis(x, Asym, xmid, scal),
              se =  FALSE, # this is important 
              linewidth = main_line_width, colour = main_line_color) + 
  ylab("") + 
  xlab("Water balance (mm)") +
  custom_scales + 
  base_theme + 
  custom_theme_figure + 
  theme(
    legend.position = "none", 
  ) + 
  scale_x_continuous(limits = c(-900, 1200)) +
  facet_wrap(~factor(y_variable2, c("EVI[Landsat]", 
                                    "NPP[MODIS]~(g~C~m^2~year^{-1})",
                                    "ABI~(g~C~m^2~year^{-1})")),
             scales = "free_y", 
             labeller = label_parsed, 
             strip.position = "left") +
  facetted_pos_scales(y = y_scales) +
  geom_text(data = (letras |> filter(tipo == "wb")), 
            aes(x = x, y = y, label = label),
            vjust = 1.2, hjust = -1.2, size = 5, 
            fontface = "bold", inherit.aes = FALSE)


### Export 

wb_prec_tmed <- plot_prec / plot_tmed /plot_wb  &
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 3)),
    colour = guide_legend(order = 1),
    shape = guide_legend(override.aes = list(stroke = .5, size = 3, colour = "black")), 
  )



# 180 mm figure size (width)
# 1 inch = 25.4 mm
# 180 mm = 7.09 inches

ggsave(
  wb_prec_tmed, 
  file = "output/plot_combined_wb_prec_tmed.png",
  dpi = 400,
  width = 7.09*1.3, height = 7.09*1.3
)

# ggsave(
#   wb_prec_tmed, 
#   file = "output/plot_combined_wb_prec_tmed.pdf",
#   dpi = 400,
#   width = 7.09*1.2, height = 7.09*1.3
# )

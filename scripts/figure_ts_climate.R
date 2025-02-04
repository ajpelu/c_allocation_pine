# title: Figure climate ts 

library(tidyverse)
source("scripts/aux.R")
library(Kendall)
library(trend)

avg_yearly <- read_csv("data/climate_avg_yearly.csv") |> 
  mutate(elev_code = fct_recode(elev_code, `low-Dec` = "low2")) |> 
  mutate(elev_code = fct_relevel(elev_code, "low-Dec","low", "medium", "high"))

aux <- avg_yearly |> 
  dplyr::select(-annual_sd, -annual_se) |> 
  pivot_wider(names_from = var, values_from = annual_value) 


# Tmed 1971, Prec 1970 
mk_tmed <- aux |> 
  dplyr::select(-prec) |> 
  na.omit() |> 
  ungroup() |> 
  group_by(sp_code, elev_code, sp_elev)|> 
  summarise(across(c(tmed), ~MannKendall(.)$tau, .names ="tau"),
            across(c(tmed), ~MannKendall(.)$sl, .names ="pvalue_mk"),
            across(c(tmed), ~trend::sens.slope(.)$estimate, .names ="senslope"),
            across(c(tmed), ~trend::sens.slope(.)$p.value, .names ="pvalue_sen")) |> 
  mutate(ypos = 
           case_when(
             elev_code == "low" ~ 14,
             elev_code == "low-Dec" ~ 8,
             elev_code == "medium" ~ 15,
             elev_code == "high" ~ 16), 
         p_value_string = symnum(pvalue_mk, corr = FALSE,
                                 cutpoints = c(0,  .001,.01,.05, 1),
                                 symbols = c("***","**","*","")), 
         variable = "tmed") |> 
  mutate(taulabel = paste(expression(tau), "==", paste0('"', round(tau, 3), p_value_string, '"'))) |> 
  mutate(
    Specie = paste0("P. ", sp_code)
  ) |> 
  mutate(taulabel = gsub("-", "\U2212", taulabel))


mk_tmed <- mk_tmed |> 
  mutate(ypos = 
           case_when(
             sp_code == "halepensis" & elev_code == "low" ~ 8,
             sp_code == "halepensis" & elev_code == "medium" ~ 9,
             sp_code == "halepensis" & elev_code == "high" ~ 10,
             TRUE ~ ypos)) 

mk_prec <- aux |> 
  ungroup() |> 
  filter(year > 1969) |> 
  group_by(sp_code, elev_code, sp_elev)|> 
  summarise(across(c(prec), ~MannKendall(.)$tau, .names ="tau"),
            across(c(prec), ~MannKendall(.)$sl, .names ="pvalue_mk"),
            across(c(prec), ~trend::sens.slope(.)$estimate, .names ="senslope"),
            across(c(prec), ~trend::sens.slope(.)$p.value, .names ="pvalue_sen")) |> 
  mutate(ypos = 
           case_when(
             elev_code == "low" ~ 1200,
             elev_code == "low-Dec" ~ 1050,
             elev_code == "medium" ~ 1350,
             elev_code == "high" ~ 1500), 
         p_value_string = symnum(pvalue_mk, corr = FALSE,
                                 cutpoints = c(0,  .001,.01,.05, 1),
                                 symbols = c("***","**","*","")),
         variable = "prec") |> 
  mutate(taulabel = paste(expression(tau), "==", paste0('"', round(tau, 3), p_value_string, '"'))) |> 
  mutate(
    Specie = paste0("P. ", sp_code)
  ) |> 
  mutate(taulabel = gsub("-", "\U2212", taulabel))



mk <- bind_rows(mk_prec, mk_tmed)

clima <- aux |> 
  filter(year > 1969) |>
  dplyr::select(-tmax, -tmin) |> 
  pivot_longer(c(prec, tmed), names_to = "variable", values_to = "mean") |> 
  mutate(Specie = paste0("P. ", sp_code)) |> 
  mutate(elev_code = fct_relevel(elev_code, "low-Dec","low", "medium", "high")) 


to_label <- as_labeller(c(
  "tmed" = "Annual Mean Temperature (ºC)", 
  "prec" = "Precipitation (mm)"))

colours_elev <- c("low-Dec" = "#8c510a",
                  "low" ="#fc8d59",
                  "medium" = "#2166ac",
                  "high" = "#72A55A")

figura_ts_clima <- clima |> 
  ggplot(aes(x = year, y = mean, group = elev_code, colour = elev_code)) +
  geom_line() +
  facet_grid(variable~factor(Specie, levels = c("P. halepensis", "P. pinaster", "P. nigra", "P. sylvestris")), 
             scales = "free_y", 
             switch = "y",
             labeller = labeller(variable = to_label)) +
  scale_colour_manual(values = colours_elev, name = "",
                      guide = guide_legend(override.aes = list(linewidth = 2))) + 
  scale_fill_manual(values = colours_elev, name="") + 
  theme_bw() +
  ylab("") +  xlab("") + 
  theme(
    text = element_text(family = "Helvetica"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.text.x = element_text(face = "italic", size = 13),
    strip.text.y = element_text(face = "bold", size = 12),
    strip.background = element_blank(),
    strip.placement = "outside", 
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 13)
  ) +
  scale_x_continuous(limits = c(1970, 2022), breaks = seq(1971, 2021, by = 10)) +
  geom_text(x = 1971, aes(y = ypos, label = taulabel), data = mk, parse = TRUE, show.legend = FALSE, hjust = "left", size = 4.5)

figura_ts_clima 


ggsave(
  figura_ts_clima, 
  file = "output/plot_ts_clima.jpg",
  dpi = 400,
  width = 7.09*1.15*1.5, height = 7.09*1.15*0.8
)

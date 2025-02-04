# title: Figure EVI vs LAI 

library(tidyverse)
source("scripts/aux.R")


stand <- read_csv("data/stand_summary.csv") |>
  dplyr::select(-locality, -site_code, -elev_category) |> 
  separate(sp_elev, into = c("sp_code", "elev_code"), remove = FALSE) |> 
  mutate(elev_code = fct_recode(elev_code, 
                                `low-Dec` = "low2",
                                "medium" = "med")) |> 
  # mutate(elev_code = fct_recode(elev_code, "medium" = "med")) |>
  unite("sp_elev", c(sp_code, elev_code), remove = FALSE) 

evi_landsat <- read_csv("data/iv_landsat.csv") |>
  filter(iv == "evi") |>
  dplyr::select(year, sp_code, elev_code, sp_elev, evi = mean) |>
  mutate(variable = "evi_landsat")  |> 
  mutate(species = paste0("P. ", sp_code))

d <- evi_landsat |> 
  left_join(
    (stand |> dplyr::select(sp_elev, lai, dentreec, bac)), by = "sp_elev") |> 
  mutate(elev_code = fct_relevel(elev_code, "low-Dec", "low", "medium", "high")) |>
  mutate(sp_code = fct_relevel(sp_code, "halepensis","pinaster", "nigra", "sylvestris")) |> 
  mutate(species = fct_relevel(species, "P. halepensis","P. pinaster", "P. nigra", "P. sylvestris")) 

summary_evi <- evi_landsat |> 
  group_by(sp_elev, elev_code, species, sp_code) |>
  summarise(
    evi_avg = mean(evi), 
    evi_avg_sd = sd(evi), 
    evi_avg_se = sqrt(evi_avg_sd/length(evi))) |> 
  left_join(
    (stand |> dplyr::select(sp_elev, lai, dentreec, bac)), by = "sp_elev")


p <- ggplot(summary_evi, 
            aes(x = lai, y = evi_avg)) +
  geom_point(data = d, 
             aes(y = evi, x = lai, fill = species, color = species), 
             size = .5, 
             alpha = 0.7) + 
  geom_point(aes(shape = elev_code, 
                 fill = species, color = species, size = bac), 
             alpha = 0.7, 
             stroke = 1) + 
  scale_shape_manual(
    values = shape_elev, 
    name = "Elevation") + 
  scale_colour_manual(
    values = colours_Specie, 
    labels = expression(italic("P. halepensis"), italic("P. pinaster"), italic("P. nigra"), italic("P. sylvestris")),
    name = "Species") +
  scale_fill_manual(
    values = colours_Specie, 
    labels = expression(italic("P. halepensis"), italic("P. pinaster"), italic("P. nigra"), italic("P. sylvestris")),
    name = "Species") +
  scale_size_continuous(name = expression(BA~(m^2~ha^{-1}))) + 
  #scale_size_continuous(name = "Tree Density\n(ind/ha)") +
  theme_bw() +
  ylab(expression(EVI[Landsat])) +
  xlab(expression(LAI~(m^2~m^{-2}))) +
  theme(
    panel.grid = element_blank()
  ) +
  scale_y_continuous(limits = c(0,0.5)) +
  scale_x_continuous(limits = c(1, 4)) 


p_evi_lai <- p + theme(
  panel.grid = element_blank(),
  legend.text = element_text(size =  7),
  legend.title = element_text(size =  8),
  legend.box = "horizontal",        
  legend.position = c(0.98, 0.02),
  legend.justification = c(1, 0),   
  legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
  axis.title = element_text(size = 14, face = "bold")
) +
  guides(
    colour = guide_legend(order = 1),
    fill = guide_legend(order = 1),
    shape = guide_legend(order = 2, 
                         override.aes = list(stroke = .5, size = 3, colour = "black")), 
    size = guide_legend(order = 3))  

nls_model_evi_lai <- nls(evi_avg ~ a + b * log(lai), data = summary_evi, 
                         start = list(a = 0.2, b = 0.1))

observed <- summary_evi$evi_avg
# Sum of Squared Residuals (SSR)
SSR <- sum((observed - predict(nls_model_evi_lai))^2)
# Total Sum of Squares (TSS)
TSS <- sum((observed - mean(observed))^2)
# R-squared
r2_evi_lai <- 1 - (SSR / TSS)
r2_evi_lai

p_evi_lai <- p_evi_lai + geom_smooth(data = summary_evi,
                                     method = "nls", 
                                     formula = y ~ a + b * log(x), 
                                     method.args = list(start = list(a = 0.2, b = 0.1)), 
                                     se = FALSE, colour = "black", linetype = "dashed") +
  annotate("text", x = 1.4, y = 0.5, label= sprintf("R^2 == %.3g", round(r2_evi_lai, 2)),
           color="black",
           parse = TRUE, 
           size = 6)

ggsave(
  p_evi_lai, 
  file = "output/plot_evi_lai.jpg",
  dpi = 400,
  width = 7.09*.8, height = 7.09*.53
)   






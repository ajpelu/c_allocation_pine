###################################### 
library(tidyverse)
source("scripts/aux.R")
library(ggnewscale)
library(likelihood)
library(nlstools)

annual_pet <- read_csv("data/spei_climate.csv")  |> 
  dplyr::select(sp_elev, year, monthly_pet, monthly_tmed, monthly_prec) |> 
  group_by(sp_elev, year) |>
  summarise(pet = sum(monthly_pet), 
            prec = sum(monthly_prec),
            tmed = mean(monthly_tmed, na.rm = TRUE)) |> 
  rowwise() |> 
  mutate(water_balance = prec - pet)

ratio <- read_csv("data/ratio_abinpp.csv") |> 
  dplyr::select(year, sp_code, elev_code, sp_elev, ratio_abi_npp = ratio) 

df <- ratio |> 
  inner_join(annual_pet) |> 
  pivot_longer(pet:water_balance, values_to = "mean_climate", 
               names_to = "climate_variable")

d <- df |> filter(climate_variable %in% c("tmed", "prec")) |> 
  pivot_wider(values_from = mean_climate, names_from = climate_variable) |> 
  dplyr::rename(rat = ratio_abi_npp) |> 
  mutate(species = paste0("P. ", sp_code)) |> 
  mutate(elev_code = fct_relevel(elev_code, "low-Dec", "low", "medium", "high")) |>
  mutate(sp_code = fct_relevel(sp_code, "halepensis","pinaster", "nigra", "sylvestris")) |> 
  mutate(species = fct_relevel(species, "P. halepensis","P. pinaster", "P. nigra", "P. sylvestris")) 


# 
# Read data from modelling # see modelling ratio selected.rmd
teffect <- read_csv("data/final_model_teffect.csv") 
peffect <- read_csv("data/final_model_peffect.csv")
summary_models <- read_csv("data/final_model_results.csv")

# Plots parameters 
line_colour_model <- "black"
alpha_points <- 0.35 
size_points <- 2.5

label_npp <- "Annual~NPP[MODIS]~(g~C~m^2~year^{-1})"
label_ratio <- "ABI:NPP"
label_prec <- "Precipitation (mm)"
label_tmed <- "Annual Mean Temperature (ºC)"
label_wb <- "Water Balance (mm)"


custom_options <- list(
  scale_shape_manual(
    values = shape_elev, 
    name = "Elevation"), 
  scale_colour_manual(
    values = colours_Specie, 
    labels = expression(italic("P. halepensis"), italic("P. pinaster"), italic("P. nigra"), italic("P. sylvestris")),
    name = "Species"),
  scale_fill_manual(
    values = colours_Specie, 
    labels = expression(italic("P. halepensis"), italic("P. pinaster"), italic("P. nigra"), italic("P. sylvestris")),
    name = "Species"),
  theme_bw(),
  theme(
    panel.grid = element_blank(), 
    axis.title = element_text(size = 15, face = "bold"), 
    axis.text = element_text(size = 12),
    axis.ticks.length = unit(.2, "cm"),
    legend.title=element_text(size=15, face = "bold"), 
    legend.text=element_text(size=15)
  )
)

# - model 
var=list(mean = "predicted", x = "rat", tmed = "tmed", prec = "prec", log = TRUE)

# temp Logistic_mod
set.seed(1234)

mtrat6 <- function(MR,k,t0, b){
  tmed_effect <- (1-b)*(1 / (1 + exp(- k *( d$tmed - t0)))) + b 
  MR*tmed_effect
}


result_mtrat6 <- anneal(
  model = mtrat6, 
  var = var, 
  source_data = d,
  par = list(MR = 0.7, b = 0.09, k = -1, t0 = 10),
  par_lo = list(MR = 0, b = 0, k = -2.5, t0 = 0),
  par_hi = list(MR = 1.2, b = 0.2, k = 2, t0 = 20),
  pdf = dnorm,
  dep_var = "rat",
  initial_temp = 5, 
  temp_red = 0.975, 
  max_iter = 200000, 
  show_display = FALSE
)

best_par <- unlist(result_mtrat6$best_par, use.names = TRUE)

plot_t <- teffect |>
  filter(type_value == "predicted") |> 
  filter(modelo == "mtrat6") |> 
  ggplot(
    aes(x = tmed, y = value))  + 
  geom_point(data = d, 
             aes(x = tmed, y = rat, shape = elev_code, fill = species, colour = species), 
             alpha = alpha_points, size = size_points) +
  geom_line(
    data = (teffect |>
              filter(type_value == "predicted") |> 
              filter(modelo == "mtrat6") |> 
              filter(tmed > min(d$tmed)) |> 
              filter(tmed < max(d$tmed))),
    colour = line_colour_model, 
    linewidth = 2.1) +
  scale_y_continuous(limits = c(0,0.63), expand = expansion(add = c(0.005, 0.005))) + 
  scale_x_continuous(limits = c(7,18), 
                     expand = expansion(add = c(0.1, 0.1)),
                     breaks = seq(7, 18, 1)) + 
  xlab(label_tmed) +
  ylab(label_ratio) +
  custom_options + 
  annotate(geom="text", 
           x=15, y=.6, 
           label= sprintf("R^2 == %.3g",round(summary_models |> filter(modelo == "mtrat6") |> dplyr::select(R2), 3)), 
           color="black",
           parse = TRUE, 
           size = 6)


# Asíntotas: Inferior es MR * b; superior MR 

asintota_inferior <- best_par["MR"] * best_par["b"]
asintota_superior <- best_par["MR"]


# Para calcular el punto de inflexión, asintota inferior + 5 % del rango 
pi <- asintota_inferior + 0.025 * (asintota_superior - asintota_inferior)

plot_t + 
  geom_hline(yintercept = asintota_inferior, linetype = "dashed", color = "red") +
  geom_hline(yintercept = asintota_superior, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = pi, linetype = "dashed", color = "black") 

### Ahora computo el valor de t* tal que f_mtrat6(t*) = pi 


# esta funcion es la inversa de f_mtrat6; la he calculado manualmente 
t_target <- function(MR, k, t0, b, y) {
  MR <- as.numeric(MR)
  k <- as.numeric(k)
  t0 <- as.numeric(t0)
  b <- as.numeric(b)
  y <- as.numeric(y)
  
  (log(1 /((y / MR - b) / (1 - b)) - 1))*(-1 / k) + t0
  
}

t_hat <- t_target(MR = best_par["MR"], 
                  k = best_par["k"],
                  b = best_par["b"], 
                  t0 = best_par["t0"],
                  y = pi)

# IC95% para ese valor de t_hat. Uso nlsBoot 
# formula 
form <- rat ~ MR * ( b + (1 - b) / (1 + exp(-k * (tmed - t0))) )

# ajusto
fit <- nls(form, data = d, 
           start = list(MR = best_par["MR"], k= best_par["k"], 
                        t0 = best_par["t0"], b  = best_par["b"]),
           control = nls.control(maxiter = 5000, minFactor = 1e-8, warnOnly = TRUE))

### bootstrap con nlsBoot (IC95% de t_hat)
set.seed(1)
boot <- nlstools::nlsBoot(fit, niter = 1000)   # re-muestreo 
bpars <- as.data.frame(boot$coefboot) 
t_boot <- mapply(
  function(MR, k, t0, b) t_target(MR, k, t0, b, y = pi),
  MR = bpars[,"MR"], k = bpars[,"k"], t0 = bpars[,"t0"], b = bpars[,"b"]
)

temp_CI <- quantile(t_boot, c(.025, .975), na.rm = TRUE)
t_hat 
temp_CI


## Final plot 
plot_t + 
  geom_hline(yintercept = asintota_inferior, linetype = "dashed", color = "red") +
  geom_hline(yintercept = asintota_superior, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = pi, linetype = "dashed", color = "black") +
  geom_vline(xintercept = t_hat, linetype = "solid", color = "black", size=1, alpha = .4) +
  annotate("rect",
           xmin = temp_CI[1], xmax = temp_CI[2],
           ymin = -Inf, ymax = Inf,
           fill = "black", alpha = 0.2) 



##############################

mprecrat3 <- function(MR, Popt, Pb)
{
  prec_effect <- exp(-0.5*(log(d$prec/Popt)/Pb)^2)
  MR * prec_effect 
}


result_mprecrat3 <- anneal(
  model = mprecrat3, 
  var = var, 
  source_data = d,
  par = list(MR = 0.7, Popt = 500, Pb = 2),
  par_lo = list(MR = 0, Popt = 300, Pb = 0),
  par_hi = list(MR = 1, Popt = 1000, Pb = 5),
  pdf = dnorm,
  dep_var = "rat",
  initial_temp = 5, 
  temp_red = 0.975, 
  max_iter = 200000, 
  show_display = FALSE
)

best_par <- unlist(result_mprecrat3$best_par, use.names = TRUE)

# mprecrat3 
plot_p <- peffect |> 
  filter(type_value == "predicted") |> 
  ggplot(
    aes(x = prec, y = value)) + 
  geom_point(data = d, 
             aes(x = prec, 
                 y = rat, 
                 shape = elev_code, fill = species, colour = species), 
             alpha = alpha_points, size = size_points) +
  geom_line(
    data = (peffect |>
              filter(type_value == "predicted") |> 
              filter(prec > min(d$prec)) |> 
              filter(prec < max(d$prec))), 
    colour = line_colour_model, 
    linewidth = 2.1) + 
  # scale_y_continuous(limits = c(0,0.63), expand = expansion(add = c(0.005, 0.005)), position = "right") + 
  scale_x_continuous(limits = c(50, 1550), 
                     breaks = c(100, 300, 500, 700, 900, 1100, 1300, 1500)) + 
  xlab(label_prec) + 
  ylab("") + 
  custom_options +
  annotate(geom="text", 
           x=1050, y=.6, 
           label = sprintf("R^2 == %.3g",round(summary_models |> filter(modelo == "mprecrat3") |> dplyr::select(R2), 3)), 
           color="black",
           parse = TRUE, 
           size = 6) 


# Asintotas horizontales, cuando prec -> +infinto, rat -> 0; y cuando prec -> 0, rat ->0 
# Máximo. El máximo de la función es en Popt 

# Para calcular el IC95% de Popt, uso nlsBoot
# formula
form <- rat ~ MR * exp(-0.5 * (log(prec / Popt) / Pb)^2)
# ajusto
fit <- nls(form, data = d, 
           start = list(MR = best_par["MR"], Popt= best_par["Popt"], Pb = best_par["Pb"]),
           control = nls.control(maxiter = 5000, minFactor = 1e-8, warnOnly = TRUE))
# bootstrap con nlsBoot (IC95% de Popt)


# maximo de ratio al 97.5 % 
ymax <- best_par["MR"] - best_par["MR"]*0.05

prec_target_left <- function(MR, Popt, Pb, y) {
  MR <- as.numeric(MR)
  Popt <- as.numeric(Popt)
  Pb <- as.numeric(Pb)
  y <- as.numeric(y)
  
  Popt * exp(-Pb * sqrt(-2 * log(y / MR))) 
  
}

prec_hat <- prec_target_left(MR = best_par["MR"], 
                             Popt = best_par["Popt"],
                             Pb = best_par["Pb"],
                             y = best_par["MR"]*.95)


form <- rat ~ MR * exp(-0.5 * (log(prec / Popt) / Pb)^2)

# Ajuste con nls; partimos de best_par
fit <- nls(
  formula = form,
  data    = d,
  start   = list(MR = as.numeric(best_par["MR"]),
                 Popt = as.numeric(best_par["Popt"]),
                 Pb   = as.numeric(best_par["Pb"])),
  control = nls.control(maxiter = 5000, minFactor = 1e-8, warnOnly = TRUE)
)

# bootstrap con nlsBoot (IC95% de Popt)
set.seed(1)
boot <- nlstools::nlsBoot(fit, niter = 1000)  # remuestreo 
bpars <- as.data.frame(boot$coefboot) # extrae matrix de parámetros


p_left_boot <- mapply(
  function(MR, Popt, Pb) prec_target_left(MR, Popt, Pb, y = 0.95 * MR),
  MR = bpars[,"MR"], Popt = bpars[,"Popt"], Pb = bpars[,"Pb"]
)

CI_left <- quantile(p_left_boot, c(0.025, 0.975), na.rm = TRUE)
prec_hat
CI_left



plot_p +
  geom_hline(yintercept = 0.95 * best_par["MR"], linetype = "dashed", color = "black") +
  geom_vline(xintercept = prec_hat,  linetype = "solid", color = "red", size = 1, alpha = .6, na.rm = TRUE) +
  annotate("rect",
           xmin = CI_left[1], xmax = 1500, # El IC superior sale enorme 
           ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.2) 
 









# prec al % del máximo para el modelo mprecrat3
prec_at_pct_of_max <- function(Popt, Pb, pct = 0.95, side = c("both","left","right")) {
  Popt <- as.numeric(Popt); Pb <- as.numeric(Pb); pct <- as.numeric(pct)
  stopifnot(Popt > 0, Pb > 0, pct > 0, pct < 1)
  s <- sqrt(-2 * log(pct))
  left  <- Popt * exp(-Pb * s)  # rama "seca" (prec < Popt)
  right <- Popt * exp( Pb * s)  # rama "húmeda" (prec > Popt)
  side <- match.arg(side)
  switch(side,
         both  = c(left = left, right = right),
         left  = left,
         right = right
  )
}


# como calculo los valors IC del prec_at_pct_of_max rama seca
prec_at_pct_of_max(Popt = best_par["Popt"], Pb = best_par["Pb"], pct = 0.95, side = "left")













































library(tidyverse)
source("scripts/aux.R")
library(ggnewscale)

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


# La asintota es: 
plot_t + 
  geom_hline(yintercept = 0.11, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 11.05441, linetype = "dashed", color = "blue") 


######## Detectar valores umbral y IC 

### Option 1. Usar segemented regression 
# modelo lineal base 
m0 <- lm(rat ~ tmed, data = d)

# ¿existe cambio de pendiente? 
# test Davies H0 -> no hay cambio de pendiente; H1 -> hay cambio de pendiente 
dt <- davies.test(m0, ~ tmed, k = 100)
# p < 0.00001 Hay cambio de pendiente 

# Regression segmentada 
# calcula usa punto inicial (puedes meter el "best" del Davies si lo tienes)
m.seg <- segmented(m0, seg.Z = ~ tmed, psi = list(tmed = as.numeric(dt$statistic)))

# calculo los IC de dos formas: 
### sencillo 
summary(m.seg)
confint(m.seg)

# Est.                 CI(95%).low  CI(95%).up
# psi1.tmed 11.0544     10.7723     11.3365

### bootstrap no paramétrico del breakpoint
set.seed(123)
nrep <- 1000
bp <- numeric(nrep) # object para guardar los breakpoints 

for (i in 1:nrep) {
  # remuestreo 
  idx <- sample(nrow(d), replace = TRUE)
  d_b <- d[idx, ]
  m0_b <- lm(rat ~ tmed, data = d_b)
  m.seg_b <- segmented(m0_b, seg.Z = ~tmed, psi = as.numeric(dt$statistic))
  bp[i] <- m.seg_b$psi[2]  # guardo el valor estimado de breakpoint  
}

# IC
quantile(bp, c(0.025, 0.975))

# 2.5%    97.5% 
# 10.77581 12.49154 


### Option 2. Usar el modelo NLS, estimar a que valor ocurre 
# la estabilización del ratio (invirtiendo la función) y hacer bootstrap de ese valor
## mas complejo 

y_star <- 0.13   # <-- nivel objetivo de ratio (asintota inferior del modelo)

# modelo 
f_mtrat6 <- function(t, MR, k, t0, b) {
  MR * ( b + (1 - b) / (1 + exp(-k * (t - t0))) )
}

# Inversa: temp* tal que f_mtrat6(temp*) = y
t_from_y_mtrat6 <- function(MR, k, t0, b, y) {
  z <- y / MR
  r <- (z - b) / (1 - b)
  t0 - (1 / k) * log(1 / r - 1)
}


## valores iniciales (vienen del ajuste )


best_par <- unlist(result_mtrat6$best_par, use.names = TRUE)
storage.mode(best_par) <- "numeric"
start_list <- list(MR = best_par["MR"],
                   k  = best_par["k"],
                   t0 = best_par["t0"],
                   b  = best_par["b"])

## ajusta una nls, si falla usa nlsLM 
form <- rat ~ MR * ( b + (1 - b) / (1 + exp(-k * (tmed - t0))) )

fit <- nls(form, data = d, start = start_list,
           control = nls.control(maxiter = 500, warnOnly = TRUE))

# library(minpack.lm)
# fit <- nlsLM(form, data = d, start = start_list,
#                    control = nls.lm.control(maxiter = 1000))

# chequear rango
y_min <- coef(fit)["MR"] * coef(fit)["b"]
y_max <- coef(fit)["MR"]


# th <- coef(fit)  # c(MR, k, t0, b)
# MRhat <- th["MR"]; khat <- th["k"]; t0hat <- th["t0"]; bhat <- th["b"]
# y_min <- MRhat * bhat
# y_max <- MRhat

# if (!(y_star > min(y_min, y_max) && y_star < max(y_min, y_max))) {
#   warning(sprintf("y* = %.3f está fuera del rango [%.3f, %.3f]; t* no existe.", 
#                   y_star, min(y_min, y_max), max(y_min, y_max)))
# }

target_temp <- t_from_y_mtrat6(MR = coef(fit)["MR"], 
                               k = coef(fit)["k"],
                               b = coef(fit)["b"], 
                               t0 = coef(fit)["t0"],
                               y = y_star)

target_temp # t* estimada 

### bootstrap con nlsBoot (IC95% de t*)
set.seed(1)
boot <- nlstools::nlsBoot(fit, niter = 1000)   # re-muestreo de residuales y re-estimación

# Asegura nombres correctos de columnas (MR,k,t0,b):
bpars <- as.data.frame(boot$coefboot)


t_boot <- mapply(
  function(MR, k, t0, b) t_from_y_mtrat6(MR, k, t0, b, y = y_star),
  MR = bpars[,"MR"], k = bpars[,"k"], t0 = bpars[,"t0"], b = bpars[,"b"]
)

temp_CI <- quantile(t_boot, c(.025, .5, .975), na.rm = TRUE)


umbral_temp <- plot_t + 
  geom_vline(xintercept = target_temp, linetype = "dotted", color = "darkgreen", size=1) +
  annotate("rect",
           xmin = temp_CI[1], xmax = temp_CI[3],
           ymin = -Inf, ymax = Inf,
           fill = "darkgreen", alpha = 0.2) 


### Precipita 

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



mprecrat3 <- function(MR, Popt, Pb)
{
  prec_effect <- exp(-0.5*(log(d$prec/Popt)/Pb)^2)
  MR * prec_effect 
}


### Option 1. Usar segemented regression 
# modelo lineal base 
m0 <- lm(rat ~ prec, data = d)

# ¿existe cambio de pendiente? 
# test Davies H0 -> no hay cambio de pendiente; H1 -> hay cambio de pendiente 
dt <- davies.test(m0, ~ prec, k = 100)
# p < 0.00001 Hay cambio de pendiente 

# Regression segmentada 
# calcula usa punto inicial (puedes meter el "best" del Davies si lo tienes)
m.seg <- segmented(m0, seg.Z = ~ prec, psi = list(prec = as.numeric(dt$statistic)))

# calculo los IC de dos formas: 
### sencillo 
summary(m.seg)
confint(m.seg)

# > confint(m.seg)
# Est. CI(95%).low CI(95%).up
# psi1.prec 622.5     462.491    782.509

### bootstrap no paramétrico del breakpoint
set.seed(123)
nrep <- 1000
bp <- numeric(nrep) # object para guardar los breakpoints 

for (i in 1:nrep) {
  # remuestreo 
  idx <- sample(nrow(d), replace = TRUE)
  d_b <- d[idx, ]
  m0_b <- lm(rat ~ prec, data = d_b)
  m.seg_b <- segmented(m0_b, seg.Z = ~prec, psi = as.numeric(dt$statistic))
  bp[i] <- m.seg_b$psi[2]  # guardo el valor estimado de breakpoint  
}

# IC
quantile(bp, c(0.025, 0.975))
# 2.5%       97.5% 
# 373.0976   758.4000 


### Option 2.

# buscar el pto maximo
Pq_humedo <- function(MR, Popt, Pb, q = 0.95) {
  Popt * exp(Pb * sqrt(-2 * log(q)))
}

bp <- unlist(result_mprecrat3$best_par, use.names = TRUE)
storage.mode(bp) <- "numeric"
start_list <- list(MR = bp["MR"], Popt = bp["Popt"], Pb = bp["Pb"])

form_p <- rat ~ MR * exp(-0.5 * (log(prec / Popt) / Pb)^2)

fit <- nls(form_p, data = d, start = start_list, control = nls.control(maxiter = 500, warnOnly = TRUE))

coef_fit <- coef(fit)
MR_hat    <- unname(coef_fit["MR"])
Popt_hat  <- unname(coef_fit["Popt"])
Pb_hat    <- unname(coef_fit["Pb"])


set.seed(1)
boot <- nlsBoot(fit, niter = 1000)                # bootstrap de residuales
bpars <- as.data.frame(boot$coefboot)
# asegura nombres
if (is.null(colnames(bpars)) || any(colnames(bpars) == "")) {
  colnames(bpars) <- names(coef_fit)
}

CI_Popt <- quantile(bpars$Popt, c(.025,.5,.975), na.rm = TRUE)
CI_MR   <- quantile(bpars$MR,   c(.025,.5,.975), na.rm = TRUE)

## 6) Punto de “meseta” al 95% del máximo (flanco húmedo) -----
q_star <- 0.95
P95_hat <- Pq_humedo(MR_hat, Popt_hat, Pb_hat, q = q_star)

P95_boot <- Pq_humedo(bpars$MR, bpars$Popt, bpars$Pb, q = q_star)
CI_P95   <- quantile(P95_boot, c(.025,.5,.975), na.rm = TRUE)

resultado <- list(
  params_hat          = c(MR = MR_hat, Popt = Popt_hat, Pb = Pb_hat),
  CI95_MR             = setNames(unname(CI_MR),   c("lower","median","upper")),
  CI95_Popt           = setNames(unname(CI_Popt), c("lower","median","upper")),
  q_star              = q_star,
  P_q_humedo_hat_mm   = P95_hat,
  CI95_P_q_humedo_mm  = setNames(unname(CI_P95),  c("lower","median","upper")),
  prop_NA_boot_Pq     = mean(is.na(P95_boot))
)
print(resultado)

























# Asegura nombres correctos de columnas (MR,k,t0,b):
bpars <- as.data.frame(boot$coefboot)


t_boot <- mapply(
  function(MR, k, t0, b) t_from_y_mtrat6(MR, k, t0, b, y = y_star),
  MR = bpars[,"MR"], k = bpars[,"k"], t0 = bpars[,"t0"], b = bpars[,"b"]
)


umbral_temp <- plot_t + 
  geom_vline(xintercept = target_temp, linetype = "dotted", color = "darkgreen", size=1) +
  annotate("rect",
           xmin = temp_CI[1], xmax = temp_CI[3],
           ymin = -Inf, ymax = Inf,
           fill = "darkgreen", alpha = 0.2) 








######## Detectar valores umbral y IC 

### Option 1. Usar segemented regression 
# modelo lineal base 
m0 <- lm(rat ~ tmed, data = d)

# ¿existe cambio de pendiente? 
# test Davies H0 -> no hay cambio de pendiente; H1 -> hay cambio de pendiente 
dt <- davies.test(m0, ~ tmed, k = 100)
# p < 0.00001 Hay cambio de pendiente 

# Regression segmentada 
# calcula usa punto inicial (puedes meter el "best" del Davies si lo tienes)
m.seg <- segmented(m0, seg.Z = ~ tmed, psi = list(tmed = as.numeric(dt$statistic)))

# calculo los IC de dos formas: 
### sencillo 
summary(m.seg)
confint(m.seg)

# Est.                 CI(95%).low  CI(95%).up
# psi1.tmed 11.0544     10.7723     11.3365

### bootstrap no paramétrico del breakpoint
set.seed(123)
nrep <- 1000
bp <- numeric(nrep) # object para guardar los breakpoints 

for (i in 1:nrep) {
  # remuestreo 
  idx <- sample(nrow(d), replace = TRUE)
  d_b <- d[idx, ]
  m0_b <- lm(rat ~ tmed, data = d_b)
  m.seg_b <- segmented(m0_b, seg.Z = ~tmed, psi = as.numeric(dt$statistic))
  bp[i] <- m.seg_b$psi[2]  # guardo el valor estimado de breakpoint  
}

# IC
quantile(bp, c(0.025, 0.975))

# 2.5%    97.5% 
# 10.77581 12.49154 


### Option 2. Usar el modelo NLS, estimar a que valor ocurre 
# la estabilización del ratio (invirtiendo la función) y hacer bootstrap de ese valor
## mas complejo 

y_star <- 0.13   # <-- nivel objetivo de ratio (asintota inferior del modelo)

# modelo 
f_mtrat6 <- function(t, MR, k, t0, b) {
  MR * ( b + (1 - b) / (1 + exp(-k * (t - t0))) )
}

# Inversa: temp* tal que f_mtrat6(temp*) = y
t_from_y_mtrat6 <- function(MR, k, t0, b, y) {
  z <- y / MR
  r <- (z - b) / (1 - b)
  t0 - (1 / k) * log(1 / r - 1)
}


## valores iniciales (vienen del ajuste )


best_par <- unlist(result_mtrat6$best_par, use.names = TRUE)
storage.mode(best_par) <- "numeric"
start_list <- list(MR = best_par["MR"],
                   k  = best_par["k"],
                   t0 = best_par["t0"],
                   b  = best_par["b"])

## ajusta una nls, si falla usa nlsLM 
form <- rat ~ MR * ( b + (1 - b) / (1 + exp(-k * (tmed - t0))) )

fit <- nls(form, data = d, start = start_list,
           control = nls.control(maxiter = 500, warnOnly = TRUE))

# library(minpack.lm)
# fit <- nlsLM(form, data = d, start = start_list,
#                    control = nls.lm.control(maxiter = 1000))

# chequear rango
y_min <- coef(fit)["MR"] * coef(fit)["b"]
y_max <- coef(fit)["MR"]


# th <- coef(fit)  # c(MR, k, t0, b)
# MRhat <- th["MR"]; khat <- th["k"]; t0hat <- th["t0"]; bhat <- th["b"]
# y_min <- MRhat * bhat
# y_max <- MRhat

# if (!(y_star > min(y_min, y_max) && y_star < max(y_min, y_max))) {
#   warning(sprintf("y* = %.3f está fuera del rango [%.3f, %.3f]; t* no existe.", 
#                   y_star, min(y_min, y_max), max(y_min, y_max)))
# }

target_temp <- t_from_y_mtrat6(MR = coef(fit)["MR"], 
                               k = coef(fit)["k"],
                               b = coef(fit)["b"], 
                               t0 = coef(fit)["t0"],
                               y = y_star)

target_temp # t* estimada 

### bootstrap con nlsBoot (IC95% de t*)
set.seed(1)
boot <- nlstools::nlsBoot(fit, niter = 1000)   # re-muestreo de residuales y re-estimación

# Asegura nombres correctos de columnas (MR,k,t0,b):
bpars <- as.data.frame(boot$coefboot)


t_boot <- mapply(
  function(MR, k, t0, b) t_from_y_mtrat6(MR, k, t0, b, y = y_star),
  MR = bpars[,"MR"], k = bpars[,"k"], t0 = bpars[,"t0"], b = bpars[,"b"]
)

temp_CI <- quantile(t_boot, c(.025, .5, .975), na.rm = TRUE)


umbral_temp <- plot_t + 
  geom_vline(xintercept = target_temp, linetype = "dotted", color = "darkgreen", size=1) +
  annotate("rect",
           xmin = temp_CI[1], xmax = temp_CI[3],
           ymin = -Inf, ymax = Inf,
           fill = "darkgreen", alpha = 0.2) 



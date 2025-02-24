title: "Modelling Ratio Selected ABI:NPP"


library(tidyverse)
source("scripts/aux.R")
library(likelihood)
library(kableExtra)

library(patchwork)


library(patchwork)
library(ggnewscale)
source("scripts/aux.R")
library(sf)
library(ggh4x)
library(geomtextpath)
library(metR)
library(ggsci)
library(fields)
library(splines)


annual_pet <- read_csv("data/spei_climate.csv") |> 
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
  mutate(species = paste0("P. ", sp_code))


# Model temperature 

# # Temperature 
# - Tmed: Log-normal and Log-modified
# 
# $$(ABI:NPP)_i=MR \times g(T_{med}) + \epsilon_i$$
#   
#   - Log-normal (**mtrat3**)
# $$g(T_{med})=\exp \left[ -\frac{1}{2} \left(\frac{\log{(\frac{T_{med}}{OT_{med}}}) }{b}\right)^2 \right]$$ 
#   siendo $OT_{med}$ el valor óptimo de $T_{med}$ al cual ocurre el máximo Ratio $(ABI:NPP)_i$; y $b$ es la desviación estándar o amplitud de la función. 
# 
# - Logistic modified (**mtrat6*)
# 
# $$g(T_{med})=\left(1 - b \right) \times \left[ \frac{1}{1+\exp \left(-k \times (T_{med} - t0) \right) } \right] + b$$ 
#   
#   
#   
  
 
var=list(mean = "predicted", x = "rat", tmed = "tmed", prec = "prec", log = TRUE)

# temp Log-normal
set.seed(123456)

mtrat3 <- function(MR, OTmed, Tb)
{
  tmed_effect <- exp(-0.5*(log(d$tmed/OTmed)/Tb)^2)
  MR * tmed_effect
}

# change model 
result_mtrat3 <- anneal(
  model = mtrat3, 
  var = var, 
  source_data = d,
  par = list(MR = 0.7, OTmed = 10.5, Tb = 1.5),
  par_lo = list(MR = 0, OTmed = 0, Tb = 0.1),
  par_hi = list(MR = 1.2, OTmed = 20, Tb = 5),
  pdf = dnorm,
  dep_var = "rat",
  initial_temp = 5, 
  temp_red = 0.975, 
  max_iter = 300000, 
  show_display = FALSE
)

# result_mtrat3 |> likelihood::write_results("output/models/univariate_temp_lognormal.txt")

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

# result_mtrat6 |> likelihood::write_results("output/models/univariate_temp_logisticmod.txt")



# Prec 
# - Precipitation: Log-normal 
# 
# $$(ABI:NPP)_i=MR \times g(P) + \epsilon_i$$
#   - Log-normal (**mprecrat3**)
# $$g(P)=\exp \left[ -\frac{1}{2} \left(\frac{\log{(\frac{P}{P{opt}}}) }{Pb}\right)^2 \right]$$ 
#   

 
## Log-normal 
set.seed(1234)

mprecrat3 <- function(MR, Popt, Pb)
{
  prec_effect <- exp(-0.5*(log(d$prec/Popt)/Pb)^2)
  MR * prec_effect 
}


# change model 
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

# result_mprecrat3 |> likelihood::write_results("output/models/univariate_prec_lognormal.txt")



## Multiplicative 

# - ***Model combrat1***
#   - log-normal temperature
# - log-normal prec


mcombrat1 <- function(MR, OTmed, Tb, Popt, Pb){
  # log-normal temperature
  # log-normal prec 
  tmed_effect <- exp(-0.5*(log(d$tmed/OTmed)/Tb)^2) 
  prec_effect <- exp(-0.5*(log(d$prec/Popt)/Pb)^2)
  MR * tmed_effect * prec_effect
}

set.seed(2023)

# change model 
result_mcombrat1 <- anneal(
  model = mcombrat1, 
  var = var, 
  source_data = d,
  par = list(MR = 0.7, OTmed = 10.5, Tb = 1.5, Popt = 500, Pb = 2), 
  par_lo = list(MR = 0, OTmed = 5, Tb = 0.01, Popt = 300, Pb = 0),
  par_hi = list(MR = 1, OTmed = 20, Tb = 5, Popt = 1000, Pb = 5),
  pdf = dnorm,
  dep_var = "rat",
  initial_temp = 5, 
  temp_red = 0.975, 
  max_iter = 200000, 
  show_display = FALSE
)

# result_mcombrat1 |> likelihood::write_results("output/models/multiplicative_tmed_lognormal_prec_lognormal.txt")



# 
# - ***Model combrat7***
#   - log-mod  temperature
# - log-normal prec


mcombrat7 <- function(MR, Tb, k, Topt, Popt, Pb){
  # log-mod temperature
  # log-normal prec 
  tmed_effect <- (1-Tb)*(1 / (1 + exp(- k *( d$tmed - Topt)))) + Tb
  prec_effect <- exp(-0.5*(log(d$prec/Popt)/Pb)^2)
  MR * tmed_effect * prec_effect
}

set.seed(123456)

result_mcombrat7 <- anneal(
  model = mcombrat7, 
  var = var, 
  source_data = d,
  par = list(MR = 0.7, Tb = 0.5, k = -1, Topt = 10, Popt = 500, Pb = 2), 
  par_lo = list(MR = 0, Tb = 0, k = -2, Topt = 0, Popt = 300, Pb = 0),
  par_hi = list(MR = 1, Tb = 1, k = 2, Topt = 20, Popt = 1000, Pb = 5),
  pdf = dnorm,
  dep_var = "rat",
  initial_temp = 5, 
  temp_red = 0.975, 
  max_iter = 200000, 
  show_display = FALSE
)


# result_mcombrat7 |> likelihood::write_results("output/models/multiplicative_tmed_logisticmod_prec_lognormal.txt")



# Compare models 
name_models <- data.frame(
  modelo = c(paste0("mtrat",c(3,6)), "mprecrat3", "mcombrat1", "mcombrat7"), 
  name_modelo = c("Temp (Log-normal)", "Temp (Logistic_Mod)", "Prec (Log-normal)",
                  "Temp (Log-normal) | Prec (Log-normal)", 
                  "Temp (Logistic_Mod) | Prec (Log-normal)"), 
  type = c(rep("univariate",3), rep("multiplicative",2))
)

all_results <-
  mget(paste0("result_", name_models$modelo)) |>
  purrr::imap_dfr(~ MLE_results_format(.x, yvar = "rat"), .id = "modelo") |>
  arrange(aic_cor) |>
  mutate(modelo = str_remove(modelo, "result_")) |> 
  inner_join(name_models) 

all_results |> 
  mutate(deltaAIC = dAIC(aic_cor)) |> 
  mutate(w = wAIC(deltaAIC)) |> 
  kbl(digits = c(0,0,0,2,0,2,2,3,3,2,3,3)) |> 
  kable_styling()



univariate <- all_results |> 
  filter(type == "univariate") |> 
  mutate(deltaAIC = dAIC(aic_cor)) |> 
  mutate(w = wAIC(deltaAIC)) 

univariate |> 
  kbl(digits = c(0,0,0,2,0,2,2,3,3,2,3,3)) |> 
  kable_styling()



multiplicative  <- all_results |> 
  filter(type == "multiplicative") |> 
  mutate(deltaAIC = dAIC(aic_cor)) |> 
  mutate(w = wAIC(deltaAIC))

multiplicative |> 
  kbl(digits = c(0,0,0,2,0,2,2,3,3,2,3,3)) |> 
  kable_styling()


#   Cuando el likelihood es similar, aplicar el Ockham’s razor (mas parsimonioso)
# 
# - Usar el modelo con menor valor de AIC (mas cercano a la "verdad")
# - Ojo el AIC, identifica el mejor modelo de un conjunto, aún cuando los modelos sean pobres 
# - Los modelos incluidos en el conjunto de modelos, tienen que ser modelos realistas, con sentido ecológico (no cualquier modelo) 
# - Usar AICc si n/k (size / nparams) es < 40 (ojo usar siempre o AIC o AICc, pero no mezclar)
# - BIC 
# 
# - Differences in AIC (Δi’s) can be used to interpret strength of evidence for one model vs. another.
# - A model with a Δ value within 1-2 of the best model has substantial support in the data, and should be considered along with the best model.
# - A Δ value within only 4-7 units of the best model has considerably less support.
# - A Δ value > 10 indicates that the worse model has virtually no support and can be omitted from further consideration.
# 
# - Akaike weights 

  

# Predicted 

pcomb1 <- result_mcombrat1$best_pars |> as.data.frame()
pcomb7 <- result_mcombrat7$best_pars |> as.data.frame()



pred_mcombrat7<- expand.grid(tmed = seq(5,20, by =.5),
                             prec = seq(100, 1500, by = 5)) |> 
  mutate(
    ratio_pred =   pcomb7$MR * ((1 - pcomb7$Tb)*(1 / (1 + exp(- pcomb7$k *(tmed - pcomb7$Topt)))) + pcomb7$Tb) * (exp(-0.5*(log(prec/pcomb7$Popt)/pcomb7$Pb)^2) )
  ) 

write_csv(pred_mcombrat7, "data/dendroapadtamed_final_model_pred_mcombrat7.csv")

breaks <- seq(0, .75, by=0.05)





### Model validation
plot_observed <- function(x, yvar, annotate = TRUE, 
                          lab_x = "Observed", 
                          lab_y = "Predicted", 
                          size_annotate = 5, ...) { 
  
  d <- x$source_data |> 
    mutate(residuals = !!sym(yvar)  - predicted) |> 
    rename(observed = !!sym(yvar))
  
  model_results <- MLE_results_format(x, yvar = "rat")
  
  max_value <- max(max(d$predicted, na.rm = TRUE), max(d$observed, na.rm = TRUE))
  max_range <- max(0, max_value)
  
  p <- ggplot(d, aes(x = predicted, y = observed)) +
    labs(x = lab_x, y = lab_y) + 
    geom_point() +
    geom_abline() +
    theme_bw() + 
    theme(panel.grid = element_blank()) +
    xlim(0, max_range) + 
    ylim(0, max_range)
  
  if (annotate) { 
    out_p <- p + 
      annotate("text", 
               x = 0.01 * max_range, y = 0.95 * max_range, 
               label = paste("R^2~'='", sprintf("%.3f", model_results$R2), sep="~"), 
               hjust = 0, vjust = 0, parse = TRUE, size = size_annotate) +
      annotate("text", 
               x = 0.01 * max_range, y = 0.85 * max_range, 
               label = paste("Slope~'='", sprintf("%.4f", model_results$slope), sep="~"), 
               hjust = 0, vjust = 0, parse = TRUE, size = size_annotate)
  } else {
    out_p <- p 
  }
  
  print(out_p)
}


p_observed <- plot_observed(result_mcombrat7, yvar = "rat", 
                            lab_y = "Observed ABI:NPP", 
                            lab_x= "Predicted ABI:NPP") 



plot_residuals <- function(x, yvar, lab_residuals = "Residuals", 
                           lab_predicted = "Predicted", ...) { 
  
  d <- x$source_data |> 
    mutate(residuals = !!sym(yvar)  - predicted) |> 
    rename(observed = !!sym(yvar))
  
  
  out <- ggplot(d, aes(x = predicted, y = residuals)) +
    labs(x = lab_predicted, y = lab_residuals) + 
    geom_point() +
    geom_hline(yintercept = 0) +  
    xlim(0, NA) +
    theme_bw() + 
    theme(panel.grid = element_blank()) 
  
  print(out)
}


p_residuos <- plot_residuals(result_mcombrat7, yvar = "rat", lab_residuals = "Residuals ABI:NPP", 
                             lab_predicted = "Predicted ABI:NPP") 



prepare_data <- function(x, yvar) {
  x$source_data |> 
    mutate(residuals = !!sym(yvar)  - predicted) |> 
    rename(observed = !!sym(yvar)) 
}

h <- prepare_data(result_mcombrat7, yvar = "rat")

dist_residuos <- ggplot(h, aes(x = residuals)) +
  geom_histogram(aes(y = after_stat(density), x = residuals), fill = "gray", col = "white") +
  stat_function(fun = dnorm, 
                args = list(mean = mean(h$residuals, na.rm = TRUE), 
                            sd = sd(h$residuals, na.rm = TRUE)), 
                color = "black", size = .85) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("Residuals") +
  ylab("Density")

plot_validation <- (p_observed + p_residuos + dist_residuos +  plot_annotation(tag_levels = "a")) & 
  theme(axis.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 14),
        plot.tag = element_text(size = 18, face = "bold"))

ggsave(
  plot_validation,
  file = "output/plot_model_validation.jpg",
  dpi = 500,
  width = 13, height = 5.5
)

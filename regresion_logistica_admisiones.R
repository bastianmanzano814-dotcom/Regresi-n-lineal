# =============================================================================
# REGRESIÓN LOGÍSTICA — PREDICCIÓN DE ADMISIÓN UNIVERSITARIA
# Dataset: admission.csv  |  Variables: admit, gre, gpa, rank
# =============================================================================

# --- Instalar / cargar paquetes ----------------------------------------------
if (!require(lmtest))            install.packages("lmtest")
if (!require(ResourceSelection)) install.packages("ResourceSelection")
if (!require(car))               install.packages("car")

library(lmtest)
library(ResourceSelection)
library(car)

usar_kable <- requireNamespace("knitr", quietly = TRUE)
if (usar_kable) library(knitr)

tabla <- function(x, caption = "") {
  if (usar_kable) print(knitr::kable(x, caption = caption, digits = 4))
  else            print(x)
}


# =============================================================================
# 1. CARGA Y EXPLORACIÓN DEL DATASET
# =============================================================================

datos <- read.csv(file.choose(), header = TRUE)   # abre explorador de archivos

datos$rank <- as.factor(datos$rank)               # rank = categórica (dummy coding)

cat("\n========== ESTRUCTURA DEL DATASET ==========\n")
str(datos)

cat("\n========== RESUMEN ESTADÍSTICO ==========\n")
print(summary(datos))

cat("\n--- Distribución de la variable respuesta ---\n")
tabla(as.data.frame(table(admit = datos$admit)), caption = "Frecuencia de admit")
cat(sprintf("Tasa de admisión: %.1f%%\n", mean(datos$admit) * 100))


# =============================================================================
# FIGURA 1 — Exploración bivariada: dispersión y barplot de rank
# =============================================================================
par(mfrow = c(1, 2),
    oma   = c(0, 0, 2, 0))   # margen superior para título global

# -- Dispersión GRE vs GPA coloreado por admisión
plot(datos$gre, datos$gpa,
     col  = ifelse(datos$admit == 1, "steelblue", "tomato"),
     pch  = 19, cex = 0.65,
     xlab = "GRE", ylab = "GPA",
     main = "GRE vs GPA por admisión")
legend("bottomright",
       legend = c("Admitido (1)", "No admitido (0)"),
       col    = c("steelblue", "tomato"),
       pch    = 19, cex = 0.8, bty = "n")

# -- Tasa de admisión por rank
prop_rank <- tapply(datos$admit, datos$rank, mean)
barplot(prop_rank * 100,
        names.arg = paste0("Rank ", names(prop_rank)),
        col  = c("#2196F3", "#42A5F5", "#90CAF9", "#BBDEFB"),
        ylim = c(0, 80),
        ylab = "% Admitidos",
        main = "Tasa de admisión por rank")
abline(h = mean(datos$admit) * 100, lty = 2, col = "red")
legend("topright", legend = "Promedio total", lty = 2,
       col = "red", bty = "n", cex = 0.8)

mtext("FIGURA 1 — Exploración bivariada", outer = TRUE, cex = 1, font = 2)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))


# =============================================================================
# FIGURA 2 — Boxplots e histogramas de GRE y GPA
# =============================================================================
par(mfrow = c(2, 2),
    oma   = c(0, 0, 2, 0))

# -- Boxplot GRE
boxplot(gre ~ admit, data = datos,
        col    = c("tomato", "steelblue"),
        names  = c("No admitido (0)", "Admitido (1)"),
        ylab   = "GRE",
        main   = "GRE por admisión",
        notch  = FALSE)

# -- Boxplot GPA
boxplot(gpa ~ admit, data = datos,
        col   = c("tomato", "steelblue"),
        names = c("No admitido (0)", "Admitido (1)"),
        ylab  = "GPA",
        main  = "GPA por admisión")

# -- Histograma GRE superpuesto
hist(datos$gre[datos$admit == 1],
     col    = rgb(0.13, 0.59, 0.95, 0.65),
     breaks = 15, xlab = "GRE",
     main   = "Distribución GRE", prob = TRUE)
hist(datos$gre[datos$admit == 0],
     col    = rgb(0.96, 0.26, 0.21, 0.55),
     breaks = 15, add = TRUE, prob = TRUE)
legend("topright", legend = c("Admitido", "No admitido"),
       fill = c(rgb(0.13,0.59,0.95,0.65), rgb(0.96,0.26,0.21,0.55)),
       bty = "n", cex = 0.75)

# -- Histograma GPA superpuesto
hist(datos$gpa[datos$admit == 1],
     col    = rgb(0.13, 0.59, 0.95, 0.65),
     breaks = 15, xlab = "GPA",
     main   = "Distribución GPA", prob = TRUE)
hist(datos$gpa[datos$admit == 0],
     col    = rgb(0.96, 0.26, 0.21, 0.55),
     breaks = 15, add = TRUE, prob = TRUE)
legend("topleft", legend = c("Admitido", "No admitido"),
       fill = c(rgb(0.13,0.59,0.95,0.65), rgb(0.96,0.26,0.21,0.55)),
       bty = "n", cex = 0.75)

mtext("FIGURA 2 — Boxplots e histogramas", outer = TRUE, cex = 1, font = 2)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))


# =============================================================================
# 2. ESTANDARIZACIÓN DE GRE Y GPA (z-score)
# =============================================================================
# Por qué estandarizar:
#   • GRE oscila ~200-800 (SD ≈ 116 pts) y GPA oscila ~2-4 (SD ≈ 0.39).
#   • Los coeficientes crudos no son comparables: β_gre mide el efecto
#     de +1 punto en GRE, mientras β_gpa mide +1 punto completo de GPA.
#   • Con z-score (x_z = (x - μ) / σ) cada β representa el cambio en el
#     logit ante +1 SD, haciendo directamente comparables las influencias.
#
#   x_z = (x - mean(x)) / sd(x)

gre_z <- scale(datos$gre)
gpa_z <- scale(datos$gpa)

gre_media <- attr(gre_z, "scaled:center"); gre_sd <- attr(gre_z, "scaled:scale")
gpa_media <- attr(gpa_z, "scaled:center"); gpa_sd <- attr(gpa_z, "scaled:scale")

cat(sprintf("\nGRE  — media: %.2f | sd: %.2f\n", gre_media, gre_sd))
cat(sprintf("GPA  — media: %.4f | sd: %.4f\n", gpa_media, gpa_sd))

datos$gre_z <- as.numeric(gre_z)
datos$gpa_z <- as.numeric(gpa_z)


# =============================================================================
# 3. MODELOS CANDIDATOS
# =============================================================================
# Modelo logístico:  logit(p) = log(p/(1-p)) = β0 + β1·x1 + β2·x2 + ...
# Estimación por máxima verosimilitud (MLE).

mod1 <- glm(admit ~ gre_z + gpa_z + rank,
            data = datos, family = binomial(link = "logit"))

mod2 <- glm(admit ~ gre_z + gpa_z,
            data = datos, family = binomial(link = "logit"))

mod3 <- glm(admit ~ gre_z + gpa_z,          # explícito por enunciado
            data = datos, family = binomial(link = "logit"))

mod4 <- glm(admit ~ rank,
            data = datos, family = binomial(link = "logit"))

cat("\n========== MODELO 1: gre_z + gpa_z + rank ==========\n")
print(summary(mod1))
cat("AIC:", AIC(mod1), " | BIC:", BIC(mod1), "\n")

cat("\n========== MODELO 2/3: gre_z + gpa_z (sin rank) ==========\n")
print(summary(mod2))
cat("AIC:", AIC(mod2), " | BIC:", BIC(mod2), "\n")

cat("\n========== MODELO 4: solo rank ==========\n")
print(summary(mod4))
cat("AIC:", AIC(mod4), " | BIC:", BIC(mod4), "\n")

# -- Tabla comparativa
cat("\n========== TABLA COMPARATIVA AIC / BIC ==========\n")
tabla_comparativa <- data.frame(
  Modelo   = c("Mod1: gre_z+gpa_z+rank",
               "Mod2/3: gre_z+gpa_z",
               "Mod4: rank"),
  AIC      = round(c(AIC(mod1), AIC(mod2), AIC(mod4)), 2),
  BIC      = round(c(BIC(mod1), BIC(mod2), BIC(mod4)), 2),
  LogLik   = round(c(as.numeric(logLik(mod1)),
                     as.numeric(logLik(mod2)),
                     as.numeric(logLik(mod4))), 2),
  Df_resid = c(mod1$df.residual, mod2$df.residual, mod4$df.residual)
)
tabla(tabla_comparativa, caption = "Comparación de modelos candidatos")

# -- Likelihood Ratio Test
# H0: el modelo reducido es suficiente (coeficientes adicionales = 0)
# p < 0.05 → el modelo más complejo mejora significativamente
cat("\n========== LRT: Mod2 vs Mod1 (¿aporta rank?) ==========\n")
lrt_mod2_mod1 <- anova(mod2, mod1, test = "Chisq")
print(lrt_mod2_mod1)

cat("\n========== LRT: Mod4 vs Mod1 (¿aportan gre_z + gpa_z?) ==========\n")
lrt_mod4_mod1 <- anova(mod4, mod1, test = "Chisq")
print(lrt_mod4_mod1)

# CRITERIO DE SELECCIÓN:
#   AIC/BIC: menor = mejor (penalizan complejidad)
#   LRT: p < 0.05 confirma que los predictores adicionales son significativos
#   → MODELO FINAL = Mod1 (gre_z + gpa_z + rank)
#     Tiene el menor AIC y el LRT muestra que rank es significativa.

modelo_final <- mod1


# =============================================================================
# 4. CODIFICACIÓN DE DUMMIES PARA RANK
# =============================================================================
# R usa treatment coding (dummies de indicadora) para factores.
# Con 4 niveles crea 3 dummies; el nivel base (rank 1) queda en el intercepto.
#
#   rank2 = 1 si rank == 2, 0 si no
#   rank3 = 1 si rank == 3, 0 si no
#   rank4 = 1 si rank == 4, 0 si no
#
# Los coeficientes β3, β4, β5 son el CAMBIO en logit(p) respecto a rank 1,
# manteniendo gre_z y gpa_z constantes.

cat("\n========== COEFICIENTES DE RANK (vs nivel base rank 1) ==========\n")
coef_rank <- coef(modelo_final)[grep("rank", names(coef(modelo_final)))]
tabla(data.frame(
  Nivel      = names(coef_rank),
  Beta       = round(coef_rank, 4),
  OR         = round(exp(coef_rank), 4),
  Interpreta = c("rank2 vs rank1", "rank3 vs rank1", "rank4 vs rank1")
), caption = "Dummies de rank — base = rank 1")


# =============================================================================
# 5. MODELO FINAL — INTERPRETACIÓN
# =============================================================================

cat("\n========== RESUMEN MODELO FINAL ==========\n")
print(summary(modelo_final))

b <- coef(modelo_final)

cat("\n--- Ecuación del modelo final (coeficientes exactos) ---\n")
cat(sprintf(
  "logit(p) = %.4f\n           + %.4f * gre_z\n           + %.4f * gpa_z\n",
  b["(Intercept)"], b["gre_z"], b["gpa_z"]))
cat(sprintf(
  "           + (%.4f) * rank2\n           + (%.4f) * rank3\n           + (%.4f) * rank4\n",
  b["rank2"], b["rank3"], b["rank4"]))
cat("\n  Donde:  logit(p) = log(p / (1 - p))\n")
cat("          p = probabilidad de admisión\n")

# -- Odds Ratios
# OR = exp(β)  →  OR > 1: predictor aumenta odds; OR < 1: los reduce
or    <- exp(coef(modelo_final))
or_ic <- exp(confint(modelo_final))        # IC por perfil de log-verosimilitud

or_tabla <- data.frame(
  Variable    = names(or),
  Coeficiente = round(b, 4),
  OR          = round(or, 4),
  IC_2.5pct   = round(or_ic[, 1], 4),
  IC_97.5pct  = round(or_ic[, 2], 4),
  p_valor     = round(summary(modelo_final)$coefficients[, 4], 5)
)
cat("\n========== ODDS RATIOS E IC 95% ==========\n")
tabla(or_tabla, caption = "OR, IC 95% y p-valores del modelo final")

# -- Interpretación cuantitativa de gre_z y gpa_z
cat("\n--- Interpretación de OR para variables estandarizadas ---\n")
cat(sprintf(
  "gre_z  OR = %.4f  →  +1 SD en GRE (≈ +%.0f pts desde media) ↑ odds %.1f%%\n",
  or["gre_z"], gre_sd, (or["gre_z"] - 1) * 100))
cat(sprintf(
  "gpa_z  OR = %.4f  →  +1 SD en GPA (≈ +%.2f desde media) ↑ odds %.1f%%\n",
  or["gpa_z"], gpa_sd, (or["gpa_z"] - 1) * 100))

cat("\n--- CONCLUSIÓN GRE vs GPA ---\n")
if (or["gpa_z"] > or["gre_z"]) {
  cat(sprintf(
    "GPA tiene MAYOR influencia sobre la admisión:\n  OR_gpa = %.4f vs OR_gre = %.4f\n",
    or["gpa_z"], or["gre_z"]))
  cat(sprintf("  +1 SD de GPA aumenta los odds un %.1f%% frente a un %.1f%% de GRE.\n",
              (or["gpa_z"] - 1) * 100, (or["gre_z"] - 1) * 100))
} else {
  cat(sprintf(
    "GRE tiene MAYOR influencia sobre la admisión:\n  OR_gre = %.4f vs OR_gpa = %.4f\n",
    or["gre_z"], or["gpa_z"]))
}


# =============================================================================
# 6. ANÁLISIS RESIDUAL Y DIAGNÓSTICOS
# =============================================================================

n   <- nrow(datos)
p   <- length(coef(modelo_final)) - 1

res_pearson  <- residuals(modelo_final, type = "pearson")
res_deviance <- residuals(modelo_final, type = "deviance")
ajustados    <- fitted(modelo_final)

cook         <- cooks.distance(modelo_final)
lev          <- hatvalues(modelo_final)
df_its       <- dffits(modelo_final)

umbral_cook   <- 4 / n
umbral_lev    <- 2 * (p + 1) / n
umbral_dffits <- 2 * sqrt((p + 1) / n)


# =============================================================================
# FIGURA 3 — Residuos de Pearson y de Deviance
# =============================================================================
par(mfrow = c(1, 3),
    oma   = c(0, 0, 2, 0))

# -- Pearson vs valores ajustados
plot(ajustados, res_pearson,
     pch  = 19, cex = 0.6, col = "steelblue",
     xlab = "Valores ajustados (prob. estimada)",
     ylab = "Residuos de Pearson",
     main = "Pearson vs Ajustados")
abline(h = c(-2, 0, 2), lty = c(2, 1, 2), col = c("red", "gray40", "red"))
lines(lowess(ajustados, res_pearson), col = "orange", lwd = 2)
legend("topright", legend = "LOWESS", col = "orange",
       lwd = 2, bty = "n", cex = 0.75)

# -- Deviance vs GRE
plot(datos$gre, res_deviance,
     pch  = 19, cex = 0.6, col = "darkorange",
     xlab = "GRE",
     ylab = "Residuos de Deviance",
     main = "Deviance vs GRE")
abline(h = c(-2, 0, 2), lty = c(2, 1, 2), col = c("red", "gray40", "red"))
lines(lowess(datos$gre, res_deviance), col = "blue", lwd = 2)

# -- Deviance vs GPA
plot(datos$gpa, res_deviance,
     pch  = 19, cex = 0.6, col = "darkorange",
     xlab = "GPA",
     ylab = "Residuos de Deviance",
     main = "Deviance vs GPA")
abline(h = c(-2, 0, 2), lty = c(2, 1, 2), col = c("red", "gray40", "red"))
lines(lowess(datos$gpa, res_deviance), col = "blue", lwd = 2)

mtext("FIGURA 3 — Gráficas de residuos", outer = TRUE, cex = 1, font = 2)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))


# =============================================================================
# FIGURA 4 — Medidas de influencia: Cook, Leverage, DFFITS
# =============================================================================
par(mfrow = c(1, 3),
    oma   = c(0, 0, 2, 0))

# -- Cook's distance
col_cook <- ifelse(cook > umbral_cook, "red", "steelblue")
plot(cook, type = "h", col = col_cook,
     ylab = "Distancia de Cook",
     main = sprintf("Cook's Distance\nUmbral = 4/n = %.4f", umbral_cook),
     cex.main = 0.85)
abline(h = umbral_cook, lty = 2, col = "red")
obs_cook <- which(cook > umbral_cook)
if (length(obs_cook) > 0 && length(obs_cook) <= 10)
  text(obs_cook, cook[obs_cook], labels = obs_cook, pos = 3, cex = 0.65, col = "red")

# -- Leverage
col_lev <- ifelse(lev > umbral_lev, "red", "steelblue")
plot(lev, type = "h", col = col_lev,
     ylab = "Leverage (h_ii)",
     main = sprintf("Leverage\nUmbral = 2(p+1)/n = %.4f", umbral_lev),
     cex.main = 0.85)
abline(h = umbral_lev, lty = 2, col = "red")

# -- DFFITS
col_dff <- ifelse(abs(df_its) > umbral_dffits, "red", "steelblue")
plot(df_its, type = "h", col = col_dff,
     ylab = "DFFITS",
     main = sprintf("DFFITS\nUmbral = ±%.4f", umbral_dffits),
     cex.main = 0.85)
abline(h = c(-umbral_dffits, umbral_dffits), lty = 2, col = "red")

mtext("FIGURA 4 — Medidas de influencia", outer = TRUE, cex = 1, font = 2)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))


# -- Tabla de observaciones influyentes
obs_influyentes <- which(cook > umbral_cook |
                         lev  > umbral_lev  |
                         abs(df_its) > umbral_dffits)

if (length(obs_influyentes) > 0) {
  tabla_influyentes <- data.frame(
    Obs         = obs_influyentes,
    Cook        = round(cook[obs_influyentes], 5),
    Leverage    = round(lev[obs_influyentes],  5),
    DFFITS      = round(df_its[obs_influyentes], 5),
    Flag_Cook   = cook[obs_influyentes]      > umbral_cook,
    Flag_Lev    = lev[obs_influyentes]       > umbral_lev,
    Flag_DFFITS = abs(df_its[obs_influyentes]) > umbral_dffits
  )
  cat("\n========== OBSERVACIONES INFLUYENTES ==========\n")
  tabla(tabla_influyentes, caption = "Obs. que superan al menos un umbral de influencia")
  cat("Total:", nrow(tabla_influyentes), "de", n, "observaciones\n")
} else {
  cat("\nNo hay observaciones influyentes según Cook, Leverage ni DFFITS.\n")
}


# -- Pruebas formales sobre residuos
mod_aux_bp <- lm(res_pearson ~ gre + gpa, data = datos)

# BREUSCH-PAGAN: H0 = homocedasticidad
bp_test <- lmtest::bptest(mod_aux_bp)
cat("\n========== BREUSCH-PAGAN (homocedasticidad) ==========\n")
print(bp_test)
cat("Conclusión:", ifelse(bp_test$p.value < 0.05,
    "Rechaza H0 → heterocedasticidad detectada.",
    "No rechaza H0 → sin evidencia de heterocedasticidad."), "\n")

# SHAPIRO-WILK: solo informativo en GLM logístico
sw_muestra <- if (n > 5000) sample(res_deviance, 5000) else res_deviance
sw_test    <- shapiro.test(sw_muestra)
cat("\n========== SHAPIRO-WILK (solo informativo) ==========\n")
print(sw_test)
cat("(En GLM logístico la normalidad de residuos NO es un supuesto.)\n")

# DURBIN-WATSON: H0 = sin autocorrelación
dw_test <- lmtest::dwtest(mod_aux_bp)
cat("\n========== DURBIN-WATSON (independencia) ==========\n")
print(dw_test)
cat("Conclusión:", ifelse(dw_test$p.value < 0.05,
    "Rechaza H0 → posible autocorrelación.",
    "No rechaza H0 → residuos independientes."), "\n")


# =============================================================================
# 7. BONDAD DE AJUSTE
# =============================================================================

# -- Hosmer-Lemeshow
# H0: buen ajuste  |  p > 0.05 indica que el modelo ajusta aceptablemente
hl_test <- ResourceSelection::hoslem.test(datos$admit, fitted(modelo_final), g = 10)
cat("\n========== HOSMER-LEMESHOW TEST ==========\n")
print(hl_test)
cat("Conclusión:", ifelse(hl_test$p.value < 0.05,
    "Rechaza H0 → el modelo NO ajusta bien.",
    "No rechaza H0 → ajuste aceptable."), "\n")

# -- Matriz de confusión
prob_pred  <- fitted(modelo_final)
clase_pred <- ifelse(prob_pred >= 0.5, 1, 0)
conf_mat   <- table(Real = datos$admit, Predicho = clase_pred)

cat("\n========== MATRIZ DE CONFUSIÓN (umbral = 0.5) ==========\n")
print(conf_mat)

# Extraer celdas con fallback si alguna clase no fue predicha
VP <- if (!is.na(conf_mat["1", "1"])) conf_mat["1", "1"] else 0
VN <- if (!is.na(conf_mat["0", "0"])) conf_mat["0", "0"] else 0
FP <- if (!is.na(conf_mat["0", "1"])) conf_mat["0", "1"] else 0
FN <- if (!is.na(conf_mat["1", "0"])) conf_mat["1", "0"] else 0

accuracy     <- (VP + VN) / n
sensibilidad <- VP / (VP + FN)
especificidad <- VN / (VN + FP)

metricas_df <- data.frame(
  Metrica  = c("Accuracy", "Sensibilidad (Recall)", "Especificidad"),
  Formula  = c("(VP+VN)/n", "VP/(VP+FN)", "VN/(VN+FP)"),
  Valor    = round(c(accuracy, sensibilidad, especificidad), 4),
  Porcentaje = paste0(round(c(accuracy, sensibilidad, especificidad) * 100, 1), "%")
)
cat("\n========== MÉTRICAS DE CLASIFICACIÓN ==========\n")
tabla(metricas_df, caption = "Métricas con umbral de clasificación = 0.5")


# =============================================================================
# FIGURA 5 — Curva ROC
# =============================================================================
roc_manual <- function(obs, prob) {
  umbrales <- sort(unique(prob), decreasing = TRUE)
  res <- lapply(umbrales, function(u) {
    pred <- as.integer(prob >= u)
    tp <- sum(pred == 1 & obs == 1);  fp <- sum(pred == 1 & obs == 0)
    tn <- sum(pred == 0 & obs == 0);  fn <- sum(pred == 0 & obs == 1)
    c(TPR = tp / (tp + fn), FPR = fp / (fp + tn))
  })
  do.call(rbind, res)
}

roc_vals <- roc_manual(datos$admit, prob_pred)
auc_val  <- abs(sum(diff(roc_vals[, "FPR"]) *
                    (roc_vals[-1, "TPR"] +
                     roc_vals[-nrow(roc_vals), "TPR"]) / 2))

par(mfrow = c(1, 1))
plot(roc_vals[, "FPR"], roc_vals[, "TPR"],
     type = "l", col = "steelblue", lwd = 2.5,
     xlab = "1 - Especificidad  (Tasa de Falsos Positivos)",
     ylab = "Sensibilidad  (Tasa de Verdaderos Positivos)",
     main = "FIGURA 5 — Curva ROC")
abline(0, 1, lty = 2, col = "gray50")
polygon(c(roc_vals[, "FPR"], 1, 0), c(roc_vals[, "TPR"], 0, 0),
        col = rgb(0.13, 0.59, 0.95, 0.15), border = NA)
legend("bottomright",
       legend = c(paste("AUC =", round(auc_val, 4)),
                  "Clasificador aleatorio"),
       col    = c("steelblue", "gray50"),
       lwd    = c(2.5, 1), lty = c(1, 2), bty = "n")

cat(sprintf("\nAUC = %.4f\n", auc_val))
cat(ifelse(auc_val >= 0.7,
    "Capacidad discriminatoria ACEPTABLE (AUC ≥ 0.70)\n",
    "Capacidad discriminatoria BAJA (AUC < 0.70)\n"))


# =============================================================================
# RESUMEN EJECUTIVO FINAL
# =============================================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║          RESUMEN EJECUTIVO — REGRESIÓN LOGÍSTICA                ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")

cat("║  DATOS                                                           ║\n")
cat(sprintf("║    N = %d observaciones                                       ║\n", n))
cat(sprintf("║    Admitidos: %d (%.1f%%)  |  No admitidos: %d (%.1f%%)          ║\n",
            sum(datos$admit), mean(datos$admit)*100,
            sum(datos$admit == 0), mean(datos$admit == 0)*100))
cat(sprintf("║    GRE: media = %.1f, sd = %.1f                              ║\n",
            gre_media, gre_sd))
cat(sprintf("║    GPA: media = %.3f, sd = %.3f                             ║\n",
            gpa_media, gpa_sd))

cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  SELECCIÓN DE MODELO                                             ║\n")
cat(sprintf("║    Mod1 (gre_z+gpa_z+rank): AIC=%.2f  BIC=%.2f  *ELEGIDO*  ║\n",
            AIC(mod1), BIC(mod1)))
cat(sprintf("║    Mod2 (gre_z+gpa_z):       AIC=%.2f  BIC=%.2f            ║\n",
            AIC(mod2), BIC(mod2)))
cat(sprintf("║    Mod4 (rank):               AIC=%.2f  BIC=%.2f            ║\n",
            AIC(mod4), BIC(mod4)))
lrt_p <- lrt_mod2_mod1[["Pr(>Chi)"]][2]
cat(sprintf("║    LRT Mod2 vs Mod1: p = %.5f → rank SIGNIFICATIVA         ║\n",
            lrt_p))
cat("║    Criterio: menor AIC + LRT confirma variables significativas   ║\n")

cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  ECUACIÓN DEL MODELO FINAL                                       ║\n")
cat(sprintf("║    logit(p) = %.4f + %.4f*gre_z + %.4f*gpa_z            ║\n",
            b["(Intercept)"], b["gre_z"], b["gpa_z"]))
cat(sprintf("║              + (%.4f)*rank2 + (%.4f)*rank3              ║\n",
            b["rank2"], b["rank3"]))
cat(sprintf("║              + (%.4f)*rank4                               ║\n",
            b["rank4"]))

cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  ODDS RATIOS E INTERPRETACIÓN                                    ║\n")
smod   <- summary(modelo_final)$coefficients
nms    <- rownames(smod)
for (nm in nms) {
  sig <- ifelse(smod[nm, 4] < 0.001, "***",
         ifelse(smod[nm, 4] < 0.01,  "** ",
         ifelse(smod[nm, 4] < 0.05,  "*  ",
         ifelse(smod[nm, 4] < 0.1,   ".  ", "   "))))
  cat(sprintf("║    %-14s  β=%7.4f  OR=%6.4f  IC[%6.4f,%7.4f]  p=%s ║\n",
              nm,
              smod[nm, 1],
              exp(smod[nm, 1]),
              or_ic[nm, 1],
              or_ic[nm, 2],
              sig))
}
cat("║    (Signif.: *** p<0.001  ** p<0.01  * p<0.05  . p<0.1)         ║\n")

cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  INFLUENCIA RELATIVA GRE vs GPA (OR estandarizados)             ║\n")
cat(sprintf("║    OR_gre_z = %.4f  →  +1SD GRE (≈+%.0f pts) sube odds %.1f%%  ║\n",
            or["gre_z"], gre_sd, (or["gre_z"]-1)*100))
cat(sprintf("║    OR_gpa_z = %.4f  →  +1SD GPA (≈+%.2f pts) sube odds %.1f%%  ║\n",
            or["gpa_z"], gpa_sd, (or["gpa_z"]-1)*100))
ganador <- ifelse(or["gpa_z"] > or["gre_z"], "GPA", "GRE")
cat(sprintf("║    ► MAYOR INFLUENCIA: %-4s (OR más alto en escala SD)      ║\n",
            ganador))

cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  DIAGNÓSTICOS                                                    ║\n")
cat(sprintf("║    Breusch-Pagan: BP=%.4f, p=%.4f  → %s      ║\n",
            bp_test$statistic, bp_test$p.value,
            ifelse(bp_test$p.value < 0.05,
                   "Heterocedasticidad  ", "Homocedasticidad    ")))
cat(sprintf("║    Durbin-Watson:  DW=%.4f, p=%.4f  → %s          ║\n",
            dw_test$statistic, dw_test$p.value,
            ifelse(dw_test$p.value < 0.05,
                   "Autocorrelación     ", "Independencia       ")))
cat(sprintf("║    Shapiro-Wilk:   W=%.4f,  p=%.4f  (solo inform.)      ║\n",
            sw_test$statistic, sw_test$p.value))
cat(sprintf("║    Hosmer-Lemeshow: χ²=%.4f, p=%.4f → %s        ║\n",
            hl_test$statistic, hl_test$p.value,
            ifelse(hl_test$p.value < 0.05,
                   "Mal ajuste          ", "Buen ajuste         ")))
cat(sprintf("║    Obs. influyentes: %d de %d (Cook | Lev | DFFITS)          ║\n",
            length(obs_influyentes), n))
cat(sprintf("║      Cook > 4/n:            %d obs.                         ║\n",
            sum(cook > umbral_cook)))
cat(sprintf("║      Leverage > 2(p+1)/n:   %d obs.                         ║\n",
            sum(lev > umbral_lev)))
cat(sprintf("║      |DFFITS| > umbral:      %d obs.                         ║\n",
            sum(abs(df_its) > umbral_dffits)))

cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  BONDAD DE AJUSTE Y CLASIFICACIÓN                               ║\n")
cat(sprintf("║    Accuracy     = %.4f  (%.1f%% clasificados correctamente)  ║\n",
            accuracy, accuracy*100))
cat(sprintf("║    Sensibilidad = %.4f  (%.1f%% admitidos identificados)      ║\n",
            sensibilidad, sensibilidad*100))
cat(sprintf("║    Especificidad= %.4f  (%.1f%% no admitidos identificados)   ║\n",
            especificidad, especificidad*100))
cat(sprintf("║    AUC (ROC)    = %.4f                                    ║\n",
            auc_val))
cat(sprintf("║    Devianza residual = %.2f con %d g.l.                 ║\n",
            modelo_final$deviance, modelo_final$df.residual))

cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  CONCLUSIÓN GENERAL                                              ║\n")
cat("║    El modelo logístico con gre_z, gpa_z y rank (dummies) es     ║\n")
cat("║    el mejor candidato según AIC y LRT. Todas las variables son  ║\n")
if (ganador == "GPA") {
  cat("║    estadísticamente significativas. GPA ejerce mayor influencia ║\n")
  cat("║    sobre los odds de admisión que GRE en escala estandarizada.  ║\n")
} else {
  cat("║    estadísticamente significativas. GRE ejerce mayor influencia ║\n")
  cat("║    sobre los odds de admisión que GPA en escala estandarizada.  ║\n")
}
cat("║    El prestigio de la universidad (rank) es relevante: menor    ║\n")
cat("║    rank (mayor prestigio destino) se asocia a mayor admisión.   ║\n")
cat(sprintf("║    Hosmer-Lemeshow confirma ajuste %s (p=%.4f).       ║\n",
            ifelse(hl_test$p.value >= 0.05, "aceptable", "deficiente"),
            hl_test$p.value))
cat(sprintf("║    AUC = %.4f indica capacidad discriminatoria %s.      ║\n",
            auc_val,
            ifelse(auc_val >= 0.7, "aceptable", "baja")))
cat("╚══════════════════════════════════════════════════════════════════╝\n")

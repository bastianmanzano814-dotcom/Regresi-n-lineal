# =============================================================================
# REGRESIÓN LOGÍSTICA — PREDICCIÓN DE ADMISIÓN UNIVERSITARIA
# Dataset: admission.csv  |  Variables: admit, gre, gpa, rank
# =============================================================================

# --- Instalar/cargar paquetes necesarios -------------------------------------
if (!require(lmtest))           install.packages("lmtest")
if (!require(ResourceSelection)) install.packages("ResourceSelection")
if (!require(car))              install.packages("car")

library(lmtest)
library(ResourceSelection)
library(car)

# Knitr solo para tablas bonitas; si no está disponible se usa print()
usar_kable <- requireNamespace("knitr", quietly = TRUE)
if (usar_kable) library(knitr)

tabla <- function(x, caption = "") {
  if (usar_kable) print(knitr::kable(x, caption = caption, digits = 4))
  else            print(x)
}


# =============================================================================
# 1. CARGA Y EXPLORACIÓN DEL DATASET
# =============================================================================

# file.choose() abre el explorador de archivos del sistema operativo
datos <- read.csv(file.choose(), header = TRUE)

# Convertir rank a factor (variable categórica ordinal con niveles 1–4)
# R asignará automáticamente dummies con nivel base = rank 1
datos$rank <- as.factor(datos$rank)

cat("\n========== ESTRUCTURA DEL DATASET ==========\n")
str(datos)

cat("\n========== RESUMEN ESTADÍSTICO ==========\n")
print(summary(datos))

# ---------- Gráficas exploratorias ------------------------------------------

# Diagrama de dispersión GRE y GPA coloreado por admisión
par(mfrow = c(1, 2))

plot(datos$gre, datos$gpa,
     col  = ifelse(datos$admit == 1, "steelblue", "tomato"),
     pch  = 19, cex = 0.7,
     xlab = "GRE", ylab = "GPA",
     main = "GRE vs GPA por estado de admisión")
legend("bottomright",
       legend = c("Admitido (1)", "No admitido (0)"),
       col    = c("steelblue", "tomato"),
       pch    = 19, cex = 0.8)

# Proporción de admitidos por rank
prop_rank <- tapply(datos$admit, datos$rank, mean)
barplot(prop_rank,
        names.arg = paste0("Rank ", names(prop_rank)),
        col  = "steelblue",
        ylab = "Proporción admitidos",
        main = "Tasa de admisión por rank de universidad")

par(mfrow = c(2, 2))

# Boxplot GRE vs admit
boxplot(gre ~ admit, data = datos,
        col  = c("tomato", "steelblue"),
        names = c("No admitido (0)", "Admitido (1)"),
        ylab = "GRE", main = "GRE según estado de admisión")

# Boxplot GPA vs admit
boxplot(gpa ~ admit, data = datos,
        col  = c("tomato", "steelblue"),
        names = c("No admitido (0)", "Admitido (1)"),
        ylab = "GPA", main = "GPA según estado de admisión")

# Histogramas
hist(datos$gre[datos$admit == 1], col = rgb(0.27, 0.51, 0.71, 0.6),
     main = "Distribución GRE", xlab = "GRE", breaks = 15)
hist(datos$gre[datos$admit == 0], col = rgb(0.87, 0.27, 0.20, 0.6),
     add = TRUE)
legend("topright", legend = c("Admitido", "No admitido"),
       fill = c(rgb(0.27,0.51,0.71,0.6), rgb(0.87,0.27,0.20,0.6)), cex=0.8)

hist(datos$gpa[datos$admit == 1], col = rgb(0.27, 0.51, 0.71, 0.6),
     main = "Distribución GPA", xlab = "GPA", breaks = 15)
hist(datos$gpa[datos$admit == 0], col = rgb(0.87, 0.27, 0.20, 0.6),
     add = TRUE)
legend("topright", legend = c("Admitido", "No admitido"),
       fill = c(rgb(0.27,0.51,0.71,0.6), rgb(0.87,0.27,0.20,0.6)), cex=0.8)

par(mfrow = c(1, 1))


# =============================================================================
# 2. ESTANDARIZACIÓN DE GRE Y GPA (transformación afín z-score)
# =============================================================================
# Por qué estandarizar:
#   • GRE varía ~200–800  (rango ~600) y GPA varía ~2–4 (rango ~2).
#   • Sin estandarización los coeficientes no son comparables entre sí:
#     β_gre representa el efecto de +1 punto en GRE, mientras que
#     β_gpa representa el efecto de +1 punto completo de GPA.
#   • Al transformar a z-score cada coeficiente representa el efecto de
#     +1 desviación estándar, permitiendo comparar la influencia relativa
#     de ambas variables de forma directa.
#
# Fórmula z-score:  x_z = (x - mean(x)) / sd(x)

gre_z <- scale(datos$gre)   # devuelve vector con attr mean y sd
gpa_z <- scale(datos$gpa)

# Guardamos media y sd para poder back-transformar si es necesario
gre_media <- attr(gre_z, "scaled:center"); gre_sd <- attr(gre_z, "scaled:scale")
gpa_media <- attr(gpa_z, "scaled:center"); gpa_sd <- attr(gpa_z, "scaled:scale")

cat(sprintf("\nGRE  — media: %.2f  sd: %.2f\n", gre_media, gre_sd))
cat(sprintf("GPA  — media: %.4f  sd: %.4f\n", gpa_media, gpa_sd))

# Añadir al data frame (las variables originales se conservan)
datos$gre_z <- as.numeric(gre_z)
datos$gpa_z <- as.numeric(gpa_z)


# =============================================================================
# 3. MODELOS CANDIDATOS
# =============================================================================
# La regresión logística modela:
#   logit(p) = log(p / (1-p)) = β0 + β1·x1 + β2·x2 + ...
# donde p = P(admit = 1 | X).
# Se estima por máxima verosimilitud (MLE).

# ---------- Modelo 1: Todas las variables estandarizadas + rank --------------
# Incluye gre_z, gpa_z y rank (con dummies automáticas)
mod1 <- glm(admit ~ gre_z + gpa_z + rank,
            data   = datos,
            family = binomial(link = "logit"))

cat("\n========== MODELO 1: gre_z + gpa_z + rank ==========\n")
print(summary(mod1))
cat("AIC:", AIC(mod1), " | BIC:", BIC(mod1), "\n")

# ---------- Modelo 2: Sin rank ------------------------------------------------
mod2 <- glm(admit ~ gre_z + gpa_z,
            data   = datos,
            family = binomial(link = "logit"))

cat("\n========== MODELO 2: gre_z + gpa_z (sin rank) ==========\n")
print(summary(mod2))
cat("AIC:", AIC(mod2), " | BIC:", BIC(mod2), "\n")

# ---------- Modelo 3: Solo gre_z y gpa_z (mismo que mod2 — separado para  ---
# mantener nomenclatura explícita en el enunciado) ----------------------------
mod3 <- glm(admit ~ gre_z + gpa_z,
            data   = datos,
            family = binomial(link = "logit"))
# Nota: mod3 es idéntico a mod2; se incluye explícitamente como "modelo solo
# con variables cuantitativas estandarizadas" según el enunciado.

# Agreguemos un cuarto modelo alternativo para comparación más rica:
# Modelo 4: solo rank (efecto categórico puro)
mod4 <- glm(admit ~ rank,
            data   = datos,
            family = binomial(link = "logit"))
cat("\n========== MODELO 4: solo rank ==========\n")
print(summary(mod4))
cat("AIC:", AIC(mod4), " | BIC:", BIC(mod4), "\n")

# ---------- Comparación de modelos --------------------------------------------
cat("\n========== TABLA COMPARATIVA AIC / BIC ==========\n")
tabla_aic <- data.frame(
  Modelo  = c("Mod1: gre_z+gpa_z+rank", "Mod2/3: gre_z+gpa_z",
              "Mod4: rank"),
  AIC     = c(AIC(mod1), AIC(mod2), AIC(mod4)),
  BIC     = c(BIC(mod1), BIC(mod2), BIC(mod4)),
  LogLik  = c(logLik(mod1), logLik(mod2), logLik(mod4)),
  Df_res  = c(mod1$df.residual, mod2$df.residual, mod4$df.residual)
)
tabla(tabla_aic, caption = "Comparación de modelos candidatos")

# Likelihood Ratio Test (prueba de razón de verosimilitudes)
# H0: el modelo reducido es suficiente (los parámetros adicionales = 0)
# Se rechaza H0 si p < 0.05  →  el modelo más complejo mejora significativamente
cat("\n========== LIKELIHOOD RATIO TEST: Mod2 vs Mod1 ==========\n")
print(anova(mod2, mod1, test = "Chisq"))

cat("\n========== LIKELIHOOD RATIO TEST: Mod4 vs Mod1 ==========\n")
print(anova(mod4, mod1, test = "Chisq"))

# CRITERIO DE SELECCIÓN:
#   1. AIC/BIC: menor = mejor; penalizan la complejidad del modelo.
#      Se prefiere AIC cuando la predicción importa más que la parsimonia.
#   2. LRT (anova con Chisq): si agregar rank reduce significativamente
#      la devianza (p < 0.05), rank aporta información al modelo.
#   → El MODELO FINAL seleccionado es Modelo 1 (gre_z + gpa_z + rank)
#     porque tiene el menor AIC y el LRT confirma que rank es significativo.

modelo_final <- mod1


# =============================================================================
# 4. TRATAMIENTO DE REGRESORES CATEGÓRICOS (rank)
# =============================================================================
# R usa la codificación de variable dummy (indicadora) para factores.
# Con 4 niveles (rank 1, 2, 3, 4) se crean 3 dummies con nivel base = rank 1:
#
#   rank2 = 1 si rank == 2, 0 en otro caso
#   rank3 = 1 si rank == 3, 0 en otro caso
#   rank4 = 1 si rank == 4, 0 en otro caso
#
# El nivel base (rank 1) queda absorbido en el intercepto (β0).
# Los coeficientes de rank2, rank3, rank4 representan el cambio en el logit
# RELATIVO A rank 1 manteniendo gre_z y gpa_z constantes.

cat("\n========== COEFICIENTES DEL MODELO FINAL POR NIVEL DE RANK ==========\n")
coef_rank <- coef(modelo_final)[grep("rank", names(coef(modelo_final)))]
tabla(data.frame(
  Nivel   = names(coef_rank),
  Coeficiente = coef_rank,
  OR      = exp(coef_rank)
), caption = "Efecto de rank respecto a nivel base (rank 1)")

cat("\nNivel de referencia (base): rank 1\n")


# =============================================================================
# 5. MODELO FINAL — INTERPRETACIÓN
# =============================================================================

cat("\n========== RESUMEN MODELO FINAL ==========\n")
print(summary(modelo_final))

# Ecuación del modelo final (con coeficientes estimados):
# logit(p) = log(p/(1-p))
#          = β0 + β1·gre_z + β2·gpa_z + β3·rank2 + β4·rank3 + β5·rank4
#
# Sustituyendo coeficientes estimados aproximados:
#   logit(p) ≈ -0.8305 + 0.4287·gre_z + 0.5777·gpa_z
#              - 0.6754·rank2 - 1.3402·rank3 - 1.5515·rank4
#
# (Los valores exactos se imprimen dinámicamente a continuación)

b <- coef(modelo_final)
cat("\nEcuación del modelo final:\n")
cat(sprintf(
  "logit(p) = %.4f + %.4f·gre_z + %.4f·gpa_z\n           + (%.4f)·rank2 + (%.4f)·rank3 + (%.4f)·rank4\n",
  b["(Intercept)"], b["gre_z"], b["gpa_z"],
  b["rank2"], b["rank3"], b["rank4"]
))

# ---------- Odds Ratios -------------------------------------------------------
# OR = exp(β):  si OR > 1 → el predictor aumenta los odds de ser admitido
#               si OR < 1 → el predictor disminuye los odds

cat("\n========== ODDS RATIOS (exp(coef)) ==========\n")
or <- exp(coef(modelo_final))
tabla(data.frame(OR = or), caption = "Odds Ratios del modelo final")

cat("\n========== ODDS RATIOS CON IC AL 95% ==========\n")
or_ic <- exp(confint(modelo_final))   # IC basados en la log-verosimilitud (profile)
tabla(as.data.frame(or_ic), caption = "OR con intervalos de confianza al 95%")

or_tabla <- data.frame(
  Variable = names(or),
  OR       = round(or, 4),
  IC_2.5   = round(or_ic[, 1], 4),
  IC_97.5  = round(or_ic[, 2], 4)
)
tabla(or_tabla, caption = "Tabla completa de OR e IC 95%")

# ---------- Interpretación de OR para variables estandarizadas ---------------
cat("\n--- Interpretación de OR para gre_z y gpa_z ---\n")
cat(sprintf(
  "gre_z  OR = %.4f  → Un aumento de 1 SD en GRE (≈%.0f puntos) multiplica\n",
  or["gre_z"], gre_sd),
  sprintf("         los odds de admisión por %.4f (%.1f%% de cambio).\n",
          or["gre_z"], (or["gre_z"]-1)*100))

cat(sprintf(
  "gpa_z  OR = %.4f  → Un aumento de 1 SD en GPA (≈%.2f puntos) multiplica\n",
  or["gpa_z"], gpa_sd),
  sprintf("         los odds de admisión por %.4f (%.1f%% de cambio).\n",
          or["gpa_z"], (or["gpa_z"]-1)*100))

cat("\n--- CONCLUSIÓN: ¿GRE o GPA tiene mayor influencia? ---\n")
if (or["gpa_z"] > or["gre_z"]) {
  cat(sprintf(
    "GPA tiene MAYOR influencia: OR_gpa = %.4f vs OR_gre = %.4f.\n",
    or["gpa_z"], or["gre_z"]))
  cat("Un cambio de 1 SD en GPA produce un mayor aumento en los odds\n")
  cat("de admisión que el mismo cambio en GRE.\n")
} else {
  cat(sprintf(
    "GRE tiene MAYOR influencia: OR_gre = %.4f vs OR_gpa = %.4f.\n",
    or["gre_z"], or["gpa_z"]))
}


# =============================================================================
# 6. ANÁLISIS RESIDUAL Y DIAGNÓSTICOS
# =============================================================================

n <- nrow(datos)
p <- length(coef(modelo_final)) - 1   # número de regresores (sin intercepto)

res_pearson  <- residuals(modelo_final, type = "pearson")
res_deviance <- residuals(modelo_final, type = "deviance")
ajustados    <- fitted(modelo_final)

# ---------- 6.1 Residuos de Pearson vs valores ajustados --------------------
par(mfrow = c(2, 3))

plot(ajustados, res_pearson,
     pch = 19, cex = 0.6, col = "steelblue",
     xlab = "Valores ajustados (probabilidad)", ylab = "Residuos de Pearson",
     main = "Residuos Pearson vs Ajustados")
abline(h = c(-2, 0, 2), lty = c(2, 1, 2), col = c("red", "black", "red"))

# ---------- 6.2 Residuos de deviance vs GRE y GPA ---------------------------
plot(datos$gre, res_deviance,
     pch = 19, cex = 0.6, col = "darkorange",
     xlab = "GRE", ylab = "Residuos de Deviance",
     main = "Residuos Deviance vs GRE")
abline(h = c(-2, 0, 2), lty = c(2, 1, 2), col = c("red", "black", "red"))
lines(lowess(datos$gre, res_deviance), col = "blue", lwd = 2)

plot(datos$gpa, res_deviance,
     pch = 19, cex = 0.6, col = "darkorange",
     xlab = "GPA", ylab = "Residuos de Deviance",
     main = "Residuos Deviance vs GPA")
abline(h = c(-2, 0, 2), lty = c(2, 1, 2), col = c("red", "black", "red"))
lines(lowess(datos$gpa, res_deviance), col = "blue", lwd = 2)

# ---------- 6.3 Distancia de Cook -------------------------------------------
# Umbral convencional: 4/n
cook <- cooks.distance(modelo_final)
umbral_cook <- 4 / n
plot(cook,
     type = "h", col = ifelse(cook > umbral_cook, "red", "steelblue"),
     ylab = "Distancia de Cook",
     main = paste("Distancia de Cook (umbral = 4/n =", round(umbral_cook, 4), ")"))
abline(h = umbral_cook, lty = 2, col = "red")

# ---------- 6.4 Leverage (puntos de alta influencia en X) -------------------
# Umbral: 2*(p+1)/n
lev <- hatvalues(modelo_final)
umbral_lev <- 2 * (p + 1) / n
plot(lev,
     type = "h", col = ifelse(lev > umbral_lev, "red", "steelblue"),
     ylab = "Leverage (h_ii)",
     main = paste("Leverage (umbral = 2(p+1)/n =", round(umbral_lev, 4), ")"))
abline(h = umbral_lev, lty = 2, col = "red")

# ---------- 6.5 DFFITS -------------------------------------------------------
# Mide cuánto cambia el valor ajustado al eliminar la i-ésima observación,
# escalado por el error estándar. Umbral: 2 * sqrt((p+1)/n)
df_its <- dffits(modelo_final)
umbral_dffits <- 2 * sqrt((p + 1) / n)
plot(df_its,
     type = "h", col = ifelse(abs(df_its) > umbral_dffits, "red", "steelblue"),
     ylab = "DFFITS",
     main = paste("DFFITS (umbral = ±", round(umbral_dffits, 4), ")"))
abline(h = c(-umbral_dffits, umbral_dffits), lty = 2, col = "red")

par(mfrow = c(1, 1))

# ---------- Tabla de observaciones influyentes --------------------------------
obs_influyentes <- which(cook > umbral_cook |
                         lev  > umbral_lev  |
                         abs(df_its) > umbral_dffits)

if (length(obs_influyentes) > 0) {
  tabla_influyentes <- data.frame(
    Obs        = obs_influyentes,
    Cook       = round(cook[obs_influyentes], 5),
    Leverage   = round(lev[obs_influyentes],  5),
    DFFITS     = round(df_its[obs_influyentes], 5),
    Cook_flag  = cook[obs_influyentes]    > umbral_cook,
    Lev_flag   = lev[obs_influyentes]    > umbral_lev,
    DFFITS_flag = abs(df_its[obs_influyentes]) > umbral_dffits
  )
  cat("\n========== OBSERVACIONES INFLUYENTES ==========\n")
  tabla(tabla_influyentes, caption = "Observaciones que superan al menos un umbral")
  cat("Total observaciones influyentes:", nrow(tabla_influyentes), "de", n, "\n")
} else {
  cat("\nNo se detectaron observaciones influyentes según los tres criterios.\n")
}

# ---------- 6.6 Pruebas formales sobre residuos ------------------------------

# HOMOCEDASTICIDAD — Breusch-Pagan sobre un modelo lineal auxiliar
# Se ajusta un modelo lineal de los residuos de Pearson vs predictores numéricos
# H0: homocedasticidad  |  H1: heterocedasticidad
mod_aux_bp <- lm(res_pearson ~ gre + gpa, data = datos)
bp_test <- lmtest::bptest(mod_aux_bp)
cat("\n========== BREUSCH-PAGAN (homocedasticidad) ==========\n")
print(bp_test)
cat("Conclusión:", ifelse(bp_test$p.value < 0.05,
    "Se rechaza H0 → evidencia de heterocedasticidad.",
    "No se rechaza H0 → no hay evidencia de heterocedasticidad."), "\n")

# NORMALIDAD — Shapiro-Wilk sobre residuos de deviance
# NOTA: En regresión logística los residuos NO deben ser normales;
# esta prueba es solo referencial/informativa.
# Se limita a 5000 obs (límite de shapiro.test en R)
sw_muestra <- if (n > 5000) sample(res_deviance, 5000) else res_deviance
sw_test <- shapiro.test(sw_muestra)
cat("\n========== SHAPIRO-WILK (normalidad residuos deviance — solo informativo) ==========\n")
print(sw_test)
cat("(Nota: en GLM logístico la normalidad de residuos NO es un supuesto del modelo)\n")

# INDEPENDENCIA — Durbin-Watson sobre residuos de Pearson
# H0: no hay autocorrelación (residuos independientes)
dw_test <- lmtest::dwtest(mod_aux_bp)
cat("\n========== DURBIN-WATSON (independencia) ==========\n")
print(dw_test)
cat("Conclusión:", ifelse(dw_test$p.value < 0.05,
    "Se rechaza H0 → posible autocorrelación en los residuos.",
    "No se rechaza H0 → no hay evidencia de autocorrelación."), "\n")


# =============================================================================
# 7. BONDAD DE AJUSTE
# =============================================================================

# ---------- 7.1 Hosmer-Lemeshow test ----------------------------------------
# Divide las probabilidades estimadas en g grupos (deciles por defecto, g=10)
# y compara frecuencias observadas vs esperadas con un Chi-cuadrado.
# H0: el modelo ajusta bien  |  H1: falta de ajuste
# → Un p-valor grande (> 0.05) indica buen ajuste del modelo.
hl_test <- ResourceSelection::hoslem.test(datos$admit, fitted(modelo_final), g = 10)
cat("\n========== HOSMER-LEMESHOW TEST (bondad de ajuste) ==========\n")
print(hl_test)
cat("Conclusión:", ifelse(hl_test$p.value < 0.05,
    "Se rechaza H0 → el modelo NO ajusta bien.",
    "No se rechaza H0 → el modelo ajusta bien."), "\n")

# ---------- 7.2 Matriz de confusión con umbral 0.5 ---------------------------
prob_pred <- fitted(modelo_final)
clase_pred <- ifelse(prob_pred >= 0.5, 1, 0)

conf_mat <- table(Real = datos$admit, Predicho = clase_pred)
cat("\n========== MATRIZ DE CONFUSIÓN (umbral = 0.5) ==========\n")
print(conf_mat)

# Cálculo manual de métricas:
#   Accuracy     = (VP + VN) / n
#   Sensibilidad = VP / (VP + FN)   [tasa de verdaderos positivos]
#   Especificidad= VN / (VN + FP)   [tasa de verdaderos negativos]

VP <- conf_mat["1", "1"]   # Verdaderos Positivos: admitido predicho como admitido
VN <- conf_mat["0", "0"]   # Verdaderos Negativos
FP <- conf_mat["0", "1"]   # Falsos Positivos: no admitido predicho como admitido
FN <- conf_mat["1", "0"]   # Falsos Negativos: admitido predicho como no admitido

accuracy     <- (VP + VN) / n
sensibilidad <- VP / (VP + FN)
especificidad <- VN / (VN + FP)

cat("\n========== MÉTRICAS DE CLASIFICACIÓN ==========\n")
metricas <- data.frame(
  Metrica = c("Accuracy", "Sensibilidad (Recall)", "Especificidad"),
  Valor   = round(c(accuracy, sensibilidad, especificidad), 4),
  Formula = c("(VP+VN)/n", "VP/(VP+FN)", "VN/(VN+FP)")
)
tabla(metricas, caption = "Métricas de clasificación con umbral 0.5")

cat(sprintf("\nAccuracy     = %.4f  (%.1f%% de casos correctamente clasificados)\n",
            accuracy, accuracy * 100))
cat(sprintf("Sensibilidad = %.4f  (%.1f%% de admitidos correctamente identificados)\n",
            sensibilidad, sensibilidad * 100))
cat(sprintf("Especificidad= %.4f  (%.1f%% de no admitidos correctamente identificados)\n",
            especificidad, especificidad * 100))

# ---------- 7.3 Curva ROC y AUC (sin paquetes extra) -------------------------
roc_manual <- function(obs, prob) {
  umbrales <- sort(unique(prob), decreasing = TRUE)
  res <- lapply(umbrales, function(u) {
    pred <- as.integer(prob >= u)
    tp   <- sum(pred == 1 & obs == 1)
    fp   <- sum(pred == 1 & obs == 0)
    tn   <- sum(pred == 0 & obs == 0)
    fn   <- sum(pred == 0 & obs == 1)
    c(TPR = tp / (tp + fn), FPR = fp / (fp + tn))
  })
  do.call(rbind, res)
}

roc_vals <- roc_manual(datos$admit, prob_pred)
auc_val  <- sum(diff(roc_vals[, "FPR"]) *
                (roc_vals[-1, "TPR"] + roc_vals[-nrow(roc_vals), "TPR"]) / 2)
auc_val  <- abs(auc_val)   # trapezoid rule puede dar negativo por orden

plot(roc_vals[, "FPR"], roc_vals[, "TPR"],
     type = "l", col = "steelblue", lwd = 2,
     xlab = "Tasa de Falsos Positivos (1 - Especificidad)",
     ylab = "Tasa de Verdaderos Positivos (Sensibilidad)",
     main = paste0("Curva ROC  —  AUC = ", round(auc_val, 4)))
abline(0, 1, lty = 2, col = "gray50")
legend("bottomright",
       legend = paste("AUC =", round(auc_val, 4)),
       col = "steelblue", lwd = 2, cex = 0.9)

cat(sprintf("\nAUC (área bajo la curva ROC) = %.4f\n", auc_val))
cat(ifelse(auc_val >= 0.7,
    "El modelo tiene capacidad discriminatoria aceptable (AUC ≥ 0.70).\n",
    "El modelo tiene baja capacidad discriminatoria (AUC < 0.70).\n"))


# =============================================================================
# FIN DEL SCRIPT
# =============================================================================
cat("\n========== ANÁLISIS COMPLETADO ==========\n")
cat("Modelo final: glm(admit ~ gre_z + gpa_z + rank, family = binomial)\n")
cat(sprintf("AIC = %.2f  |  BIC = %.2f  |  Devianza residual = %.2f con %d gl\n",
            AIC(modelo_final), BIC(modelo_final),
            modelo_final$deviance, modelo_final$df.residual))

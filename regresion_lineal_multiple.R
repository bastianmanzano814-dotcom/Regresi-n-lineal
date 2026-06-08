# =============================================================================
# REGRESIÓN LINEAL MÚLTIPLE — SCRIPT GENÉRICO
# El usuario solo especifica Y y X; el script construye y analiza todo solo.
# =============================================================================

# --- Paquetes ----------------------------------------------------------------
if (!require(lmtest))   install.packages("lmtest")
if (!require(car))      install.packages("car")
if (!require(nortest))  install.packages("nortest")   # Anderson-Darling

library(lmtest)
library(car)
library(nortest)

usar_kable <- requireNamespace("knitr", quietly = TRUE)
if (usar_kable) library(knitr)
tabla <- function(x, caption = "") {
  if (usar_kable) print(knitr::kable(x, caption = caption, digits = 4))
  else { cat(if (nchar(caption)) paste0("\n-- ", caption, " --\n")); print(x) }
}


# =============================================================================
#  ► CONFIGURACIÓN — ÚNICA SECCIÓN QUE EL USUARIO DEBE EDITAR ◄
# =============================================================================

NOMBRE_Y  <- "admit"          # nombre exacto de la columna respuesta (Y)

NOMBRES_X <- NULL             # vector de nombres de columnas predictoras (X)
                              # Ej: c("gre", "gpa", "rank")
                              # NULL = usar TODAS las columnas excepto Y


# =============================================================================
# 1. CARGA DE DATOS
# =============================================================================

datos <- read.csv(file.choose(), header = TRUE)

cat("\n========== DATASET CARGADO ==========\n")
cat("Dimensiones:", nrow(datos), "filas x", ncol(datos), "columnas\n")
cat("Columnas disponibles:", paste(names(datos), collapse = ", "), "\n")

# Validar que Y existe
if (!NOMBRE_Y %in% names(datos))
  stop(sprintf("La columna '%s' no existe en el dataset.", NOMBRE_Y))

# Si X es NULL, usar todas las columnas menos Y
if (is.null(NOMBRES_X)) {
  NOMBRES_X <- setdiff(names(datos), NOMBRE_Y)
  cat("X inferidas automáticamente:", paste(NOMBRES_X, collapse = ", "), "\n")
} else {
  faltantes <- setdiff(NOMBRES_X, names(datos))
  if (length(faltantes) > 0)
    stop(sprintf("Columnas no encontradas: %s", paste(faltantes, collapse = ", ")))
}

# Subconjunto de trabajo (solo columnas usadas, sin NAs)
df <- na.omit(datos[, c(NOMBRE_Y, NOMBRES_X)])
cat(sprintf("Obs. disponibles (sin NAs): %d\n", nrow(df)))

# Variables de trabajo
y_vec <- df[[NOMBRE_Y]]
n     <- nrow(df)
k     <- length(NOMBRES_X)    # número de predictores

cat(sprintf("\nY  = '%s'\n",   NOMBRE_Y))
cat(sprintf("X  = [%s]\n",    paste(NOMBRES_X, collapse = ", ")))
cat(sprintf("n  = %d  |  k = %d predictores\n", n, k))


# =============================================================================
# 2. EXPLORACIÓN INICIAL
# =============================================================================

cat("\n========== RESUMEN ESTADÍSTICO ==========\n")
print(summary(df))

cat("\n========== CORRELACIONES CON Y ==========\n")
cor_y <- sapply(NOMBRES_X, function(x) cor(df[[x]], y_vec, use = "complete.obs"))
tabla(data.frame(
  Predictor   = names(cor_y),
  Correlacion = round(cor_y, 4),
  Abs_cor     = round(abs(cor_y), 4)
)[order(abs(cor_y), decreasing = TRUE), ],
caption = paste("Correlación de Pearson con", NOMBRE_Y))


# =============================================================================
# FIGURA 1 — Matriz de dispersión
# =============================================================================
pares_cols <- c(NOMBRE_Y, NOMBRES_X)
pares_max  <- min(length(pares_cols), 6)   # limitar a 6 para legibilidad

pairs(df[, pares_cols[seq_len(pares_max)]],
      pch  = 19, cex = 0.4, col = rgb(0.2, 0.4, 0.8, 0.5),
      main = "FIGURA 1 — Matriz de dispersión (primeras 6 variables)")


# =============================================================================
# FIGURA 2 — Boxplots e histograma de Y
# =============================================================================
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))

boxplot(y_vec,
        col  = "steelblue", horizontal = TRUE,
        xlab = NOMBRE_Y,
        main = paste("Boxplot de", NOMBRE_Y))

hist(y_vec,
     col  = rgb(0.13, 0.59, 0.95, 0.7),
     breaks = 20, prob = TRUE,
     xlab = NOMBRE_Y, main = paste("Distribución de", NOMBRE_Y))
lines(density(y_vec), col = "navy", lwd = 2)

mtext("FIGURA 2 — Variable respuesta", outer = TRUE, cex = 1, font = 2)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))


# =============================================================================
# 3. AJUSTE DEL MODELO DE REGRESIÓN LINEAL MÚLTIPLE
# =============================================================================
# Modelo:  Y = β0 + β1·X1 + β2·X2 + ... + βk·Xk + ε
#          ε ~ N(0, σ²)  (supuesto de normalidad e independencia)

formula_str <- paste(NOMBRE_Y, "~", paste(NOMBRES_X, collapse = " + "))
formula_mod <- as.formula(formula_str)

cat(sprintf("\n========== MODELO: %s ==========\n", formula_str))

modelo <- lm(formula_mod, data = df)

cat("\n--- Resumen del modelo ---\n")
smod <- summary(modelo)
print(smod)

# Guardar valores derivados
coefs       <- coef(modelo)
fitted_vals <- fitted(modelo)
res_std     <- rstandard(modelo)      # residuos estandarizados
res_stud    <- rstudent(modelo)       # residuos studentizados (jackknife)
res_raw     <- residuals(modelo)
cook        <- cooks.distance(modelo)
lev         <- hatvalues(modelo)
df_its      <- dffits(modelo)

p_total <- k + 1                      # parámetros (β0 + k βs)
umbral_cook   <- 4 / n
umbral_lev    <- 2 * p_total / n
umbral_dffits <- 2 * sqrt(p_total / n)

cat(sprintf("\nR²          = %.4f\n",    smod$r.squared))
cat(sprintf("R² ajustado = %.4f\n",    smod$adj.r.squared))
cat(sprintf("Error std.  = %.4f\n",    smod$sigma))
cat(sprintf("F-estadístico = %.4f  (p = %.6f)\n",
            smod$fstatistic[1],
            pf(smod$fstatistic[1], smod$fstatistic[2],
               smod$fstatistic[3], lower.tail = FALSE)))

# -- Tabla de coeficientes enriquecida
coef_tab <- as.data.frame(smod$coefficients)
names(coef_tab) <- c("Estimado", "Std_Error", "t_valor", "p_valor")
coef_tab$IC_2.5  <- confint(modelo)[, 1]
coef_tab$IC_97.5 <- confint(modelo)[, 2]
coef_tab$Signif  <- symnum(coef_tab$p_valor,
                            cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                            symbols   = c("***", "**", "*", ".", " "))
cat("\n========== COEFICIENTES CON IC 95% ==========\n")
tabla(coef_tab, caption = "Coeficientes, IC 95% y significancia")


# =============================================================================
# 4. MULTICOLINEALIDAD — Factor de Inflación de Varianza (VIF)
# =============================================================================
# VIF > 5 (o 10): posible multicolinealidad problemática
# VIF = 1 / (1 - R²_j) donde R²_j proviene de regresar Xj vs resto de X

if (k >= 2) {
  vif_vals <- car::vif(modelo)
  cat("\n========== VIF — MULTICOLINEALIDAD ==========\n")
  vif_df <- if (is.matrix(vif_vals)) {
    data.frame(Predictor = rownames(vif_vals),
               VIF       = round(vif_vals[, 1], 4),
               Alerta    = ifelse(vif_vals[, 1] > 10, "ALTA",
                           ifelse(vif_vals[, 1] > 5, "MEDIA", "OK")))
  } else {
    data.frame(Predictor = names(vif_vals),
               VIF       = round(vif_vals, 4),
               Alerta    = ifelse(vif_vals > 10, "ALTA",
                           ifelse(vif_vals > 5, "MEDIA", "OK")))
  }
  tabla(vif_df, caption = "VIF por predictor (umbral problemático > 5 o > 10)")
} else {
  cat("\n(VIF requiere al menos 2 predictores)\n")
}


# =============================================================================
# FIGURA 3 — Diagnóstico global del modelo (4 gráficas estándar de lm)
# =============================================================================
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
plot(modelo, which = 1:4, col = "steelblue", pch = 19, cex = 0.6)
mtext("FIGURA 3 — Diagnósticos estándar del modelo", outer = TRUE, cex = 1, font = 2)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))


# =============================================================================
# 5. ANÁLISIS RESIDUAL POR REGRESOR
# =============================================================================
# Para cada predictor Xi se generan:
#   a) Residuos estandarizados vs Xi  (detecta no linealidad y outliers por Xi)
#   b) Gráfica de regresión parcial (added-variable plot)

# Calcular layout dinámico según k
cols_fig  <- min(k, 3)
filas_fig <- ceiling(k / cols_fig)

# ---- 5a. Residuos estandarizados vs cada Xi --------------------------------
cat("\n--- Residuos estandarizados vs cada predictor ---\n")

par(mfrow = c(filas_fig, cols_fig),
    oma   = c(0, 0, 2, 0),
    mar   = c(4, 4, 3, 1))

for (xi in NOMBRES_X) {
  x_vals <- df[[xi]]

  plot(x_vals, res_std,
       pch  = 19, cex  = 0.6,
       col  = rgb(0.13, 0.59, 0.95, 0.6),
       xlab = xi,
       ylab = "Residuos estand.",
       main = paste("Res. estand. vs", xi))

  abline(h = c(-2, 0, 2), lty = c(2, 1, 2),
         col = c("red", "gray30", "red"), lwd = c(1, 1, 1))

  # Curva LOWESS para detectar tendencias no lineales
  lw <- lowess(x_vals, res_std)
  lines(lw, col = "orange", lwd = 2)

  # Marcar residuos fuera de ±2
  out_idx <- which(abs(res_std) > 2)
  if (length(out_idx) > 0 && length(out_idx) <= 15)
    text(x_vals[out_idx], res_std[out_idx],
         labels = out_idx, pos = 3, cex = 0.55, col = "red")
}

# Rellenar celdas vacías si el layout no es perfecto
n_vacias <- filas_fig * cols_fig - k
if (n_vacias > 0) for (i in seq_len(n_vacias)) plot.new()

mtext(paste("FIGURA 4 — Residuos estand. vs cada predictor"),
      outer = TRUE, cex = 1, font = 2)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(5, 4, 4, 2) + 0.1)


# ---- 5b. Added-Variable Plots (regresión parcial) --------------------------
# Muestra la relación lineal NETA de cada Xi con Y después de controlar
# el efecto de todos los demás predictores.
if (k >= 2) {
  par(oma = c(0, 0, 2, 0))
  car::avPlots(modelo,
               pch     = 19,
               cex     = 0.6,
               col     = rgb(0.2, 0.6, 0.2, 0.5),
               col.lines = "red",
               main    = "")
  mtext("FIGURA 5 — Added-Variable Plots (efecto parcial de cada Xi)",
        outer = TRUE, cex = 1, font = 2)
  par(oma = c(0, 0, 0, 0))
}


# ---- 5c. Residuos estandarizados vs valores ajustados ----------------------
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))

plot(fitted_vals, res_std,
     pch  = 19, cex = 0.6,
     col  = rgb(0.13, 0.59, 0.95, 0.6),
     xlab = "Valores ajustados (Ŷ)",
     ylab = "Residuos estand.",
     main = "Res. estand. vs Ŷ")
abline(h = c(-2, 0, 2), lty = c(2, 1, 2),
       col = c("red", "gray30", "red"))
lines(lowess(fitted_vals, res_std), col = "orange", lwd = 2)
legend("topright", legend = "LOWESS", col = "orange",
       lwd = 2, bty = "n", cex = 0.8)

# Residuos vs orden de observación (independencia)
plot(seq_along(res_std), res_std,
     type = "b", pch = 19, cex = 0.5,
     col  = rgb(0.13, 0.59, 0.95, 0.7),
     xlab = "Índice de observación",
     ylab = "Residuos estand.",
     main = "Res. estand. vs Orden")
abline(h = c(-2, 0, 2), lty = c(2, 1, 2),
       col = c("red", "gray30", "red"))

mtext("FIGURA 6 — Residuos vs ajustados y vs orden",
      outer = TRUE, cex = 1, font = 2)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))


# ---- 5d. Q-Q plot y densidad de residuos (normalidad) ----------------------
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))

qqnorm(res_std,
       pch  = 19, cex = 0.6,
       col  = "steelblue",
       main = "Q-Q Normal — Residuos estand.")
qqline(res_std, col = "red", lwd = 2)

plot(density(res_std),
     col  = "steelblue", lwd = 2,
     main = "Densidad — Residuos estand.",
     xlab = "Residuos estand.")
x_seq <- seq(min(res_std) - 1, max(res_std) + 1, length = 200)
lines(x_seq, dnorm(x_seq), col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("Empírica", "N(0,1) teórica"),
       col = c("steelblue", "red"), lwd = 2, lty = c(1, 2),
       bty = "n", cex = 0.8)

mtext("FIGURA 7 — Normalidad de residuos",
      outer = TRUE, cex = 1, font = 2)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))


# =============================================================================
# FIGURA 8 — Medidas de influencia: Cook, Leverage, DFFITS
# =============================================================================
par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))

# -- Cook's distance
col_cook <- ifelse(cook > umbral_cook, "red", "steelblue")
plot(cook, type = "h", col = col_cook,
     ylab = "Distancia de Cook",
     main = sprintf("Cook's D\nUmbral = 4/n = %.4f", umbral_cook),
     cex.main = 0.85)
abline(h = umbral_cook, lty = 2, col = "red")
obs_c <- which(cook > umbral_cook)
if (length(obs_c) > 0 && length(obs_c) <= 10)
  text(obs_c, cook[obs_c], labels = obs_c, pos = 3, cex = 0.6, col = "red")

# -- Leverage
col_lev <- ifelse(lev > umbral_lev, "red", "steelblue")
plot(lev, type = "h", col = col_lev,
     ylab = "Leverage (h_ii)",
     main = sprintf("Leverage\nUmbral = 2p/n = %.4f", umbral_lev),
     cex.main = 0.85)
abline(h = umbral_lev, lty = 2, col = "red")
obs_l <- which(lev > umbral_lev)
if (length(obs_l) > 0 && length(obs_l) <= 10)
  text(obs_l, lev[obs_l], labels = obs_l, pos = 3, cex = 0.6, col = "red")

# -- DFFITS
col_dff <- ifelse(abs(df_its) > umbral_dffits, "red", "steelblue")
plot(df_its, type = "h", col = col_dff,
     ylab = "DFFITS",
     main = sprintf("DFFITS\nUmbral = ±%.4f", umbral_dffits),
     cex.main = 0.85)
abline(h = c(-umbral_dffits, umbral_dffits), lty = 2, col = "red")
obs_d <- which(abs(df_its) > umbral_dffits)
if (length(obs_d) > 0 && length(obs_d) <= 10)
  text(obs_d, df_its[obs_d], labels = obs_d, pos = 3, cex = 0.6, col = "red")

mtext("FIGURA 8 — Medidas de influencia", outer = TRUE, cex = 1, font = 2)
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))

# -- Tabla de observaciones influyentes
obs_infl <- sort(unique(c(obs_c, obs_l, obs_d)))
if (length(obs_infl) > 0) {
  inf_df <- data.frame(
    Obs         = obs_infl,
    Cook        = round(cook[obs_infl], 5),
    Leverage    = round(lev[obs_infl],  5),
    DFFITS      = round(df_its[obs_infl], 5),
    Flag_Cook   = cook[obs_infl]        > umbral_cook,
    Flag_Lev    = lev[obs_infl]         > umbral_lev,
    Flag_DFFITS = abs(df_its[obs_infl]) > umbral_dffits
  )
  cat("\n========== OBSERVACIONES INFLUYENTES ==========\n")
  tabla(inf_df, caption = "Obs. que superan al menos un umbral de influencia")
  cat("Total:", nrow(inf_df), "de", n, "\n")
} else {
  cat("\nNo hay observaciones influyentes según Cook, Leverage ni DFFITS.\n")
}


# =============================================================================
# 6. PRUEBAS FORMALES SOBRE RESIDUOS
# =============================================================================

# -- BREUSCH-PAGAN: H0 = homocedasticidad (varianza constante de ε)
bp_test <- lmtest::bptest(modelo)
cat("\n========== BREUSCH-PAGAN (homocedasticidad) ==========\n")
print(bp_test)
cat("Conclusión:", ifelse(bp_test$p.value < 0.05,
    "Rechaza H0 → heterocedasticidad detectada.",
    "No rechaza H0 → sin evidencia de heterocedasticidad."), "\n")

# -- SHAPIRO-WILK: H0 = normalidad de residuos
# (límite de 5000 obs de la función shapiro.test)
sw_muestra <- if (n > 5000) sample(res_raw, 5000) else res_raw
sw_test <- shapiro.test(sw_muestra)
cat("\n========== SHAPIRO-WILK (normalidad de residuos) ==========\n")
print(sw_test)
cat("Conclusión:", ifelse(sw_test$p.value < 0.05,
    "Rechaza H0 → residuos NO normales.",
    "No rechaza H0 → residuos compatibles con normalidad."), "\n")

# -- ANDERSON-DARLING (más potente que SW para n grande)
ad_test <- nortest::ad.test(res_raw)
cat("\n========== ANDERSON-DARLING (normalidad) ==========\n")
print(ad_test)
cat("Conclusión:", ifelse(ad_test$p.value < 0.05,
    "Rechaza H0 → residuos NO normales.",
    "No rechaza H0 → normalidad aceptable."), "\n")

# -- DURBIN-WATSON: H0 = sin autocorrelación de residuos
dw_test <- lmtest::dwtest(modelo)
cat("\n========== DURBIN-WATSON (independencia / autocorrelación) ==========\n")
print(dw_test)
cat("Conclusión:", ifelse(dw_test$p.value < 0.05,
    "Rechaza H0 → posible autocorrelación en residuos.",
    "No rechaza H0 → residuos independientes."), "\n")

# -- RESET TEST de Ramsey: H0 = forma funcional correcta (linealidad)
reset_test <- lmtest::resettest(modelo, power = 2:3)
cat("\n========== RAMSEY RESET (linealidad) ==========\n")
print(reset_test)
cat("Conclusión:", ifelse(reset_test$p.value < 0.05,
    "Rechaza H0 → posible mala especificación (no linealidad).",
    "No rechaza H0 → forma funcional lineal aceptable."), "\n")


# =============================================================================
# 7. BONDAD DE AJUSTE — ANOVA DE REGRESIÓN
# =============================================================================

cat("\n========== TABLA ANOVA DE REGRESIÓN ==========\n")
print(anova(modelo))

# SCR, SCE, SCT manual para verificación
SCT <- sum((y_vec - mean(y_vec))^2)
SCR <- sum((fitted_vals - mean(y_vec))^2)
SCE <- sum(res_raw^2)

cat(sprintf("\nSCT = %.4f  (variación total de Y)\n", SCT))
cat(sprintf("SCR = %.4f  (variación explicada por el modelo = %.1f%%)\n",
            SCR, SCR/SCT*100))
cat(sprintf("SCE = %.4f  (variación no explicada / residual)\n", SCE))
cat(sprintf("R²  = SCR/SCT = %.4f\n", SCR/SCT))


# =============================================================================
# FIGURA 9 — Y observado vs Y predicho
# =============================================================================
rango <- range(c(y_vec, fitted_vals))

plot(fitted_vals, y_vec,
     pch  = 19, cex = 0.7,
     col  = rgb(0.13, 0.59, 0.95, 0.6),
     xlim = rango, ylim = rango,
     xlab = paste("Ŷ — Valores ajustados"),
     ylab = paste(NOMBRE_Y, "— Valores observados"),
     main = "FIGURA 9 — Observado vs Ajustado")
abline(0, 1, col = "red", lwd = 2, lty = 2)
abline(lm(y_vec ~ fitted_vals), col = "orange", lwd = 1.5)
legend("topleft",
       legend = c("Línea perfecta (y=ŷ)", "Regresión obs~pred"),
       col    = c("red", "orange"),
       lwd    = 2, lty = c(2, 1), bty = "n", cex = 0.8)


# =============================================================================
# RESUMEN EJECUTIVO FINAL
# =============================================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║        RESUMEN EJECUTIVO — REGRESIÓN LINEAL MÚLTIPLE            ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")

cat(sprintf("║  DATOS                                                           ║\n"))
cat(sprintf("║    N = %d obs.  |  k = %d predictor(es)  |  p = %d parámetros    ║\n",
            n, k, p_total))
cat(sprintf("║    Y = %-15s  rango: [%.3f, %.3f]               ║\n",
            NOMBRE_Y, min(y_vec), max(y_vec)))
for (xi in NOMBRES_X)
  cat(sprintf("║    X: %-15s  rango: [%.3f, %.3f]               ║\n",
              xi, min(df[[xi]]), max(df[[xi]])))

cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  AJUSTE DEL MODELO                                               ║\n")
cat(sprintf("║    R²          = %.4f                                         ║\n",
            smod$r.squared))
cat(sprintf("║    R² ajustado = %.4f                                         ║\n",
            smod$adj.r.squared))
cat(sprintf("║    Error std.  = %.4f                                         ║\n",
            smod$sigma))
f_stat  <- smod$fstatistic
f_pval  <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
cat(sprintf("║    F(%d,%d)     = %.4f  p = %.6f  → modelo %s       ║\n",
            f_stat[2], f_stat[3], f_stat[1], f_pval,
            ifelse(f_pval < 0.05, "SIGNIFICATIVO  ", "NO significativo")))

cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  ECUACIÓN AJUSTADA                                               ║\n")
cat(sprintf("║    Ŷ = %.4f\n", coefs["(Intercept)"]))
for (xi in NOMBRES_X)
  cat(sprintf("║        %s %.4f * %s\n",
              ifelse(coefs[xi] >= 0, "+", "-"), abs(coefs[xi]), xi))
cat("║                                                                  ║\n")

cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  COEFICIENTES                                                    ║\n")
cat(sprintf("║    %-16s  %8s  %8s  %8s  %-6s ║\n",
            "Variable", "Beta", "t", "p-val", "Signif"))
for (nm in rownames(smod$coefficients)) {
  cc  <- smod$coefficients[nm, ]
  sig <- symnum(cc[4], cutpoints = c(0,.001,.01,.05,.1,1),
                symbols = c("***","** ","*  ",".  ","   "))
  cat(sprintf("║    %-16s  %8.4f  %8.4f  %8.5f  %-6s ║\n",
              nm, cc[1], cc[3], cc[4], as.character(sig)))
}
cat("║    (*** p<.001  ** p<.01  * p<.05  . p<.1)                      ║\n")

if (k >= 2) {
  cat("╠══════════════════════════════════════════════════════════════════╣\n")
  cat("║  MULTICOLINEALIDAD (VIF)                                         ║\n")
  for (xi in names(vif_vals)) {
    vv <- if (is.matrix(vif_vals)) vif_vals[xi, 1] else vif_vals[xi]
    cat(sprintf("║    %-16s  VIF = %6.3f  → %s                     ║\n",
                xi, vv,
                ifelse(vv > 10, "ALTA COLINEALIDAD ",
                ifelse(vv > 5,  "COLINEALIDAD MEDIA", "OK               "))))
  }
}

cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  SUPUESTOS DEL MODELO                                            ║\n")
cat(sprintf("║    Breusch-Pagan : BP=%.4f  p=%.5f  → %s      ║\n",
            bp_test$statistic, bp_test$p.value,
            ifelse(bp_test$p.value < 0.05,
                   "Heterocedasticidad  ", "Homocedasticidad    ")))
cat(sprintf("║    Shapiro-Wilk  :  W=%.4f  p=%.5f  → %s      ║\n",
            sw_test$statistic, sw_test$p.value,
            ifelse(sw_test$p.value < 0.05,
                   "No normalidad       ", "Normalidad ok       ")))
cat(sprintf("║    Anderson-Darling: A=%.4f p=%.5f  → %s      ║\n",
            ad_test$statistic, ad_test$p.value,
            ifelse(ad_test$p.value < 0.05,
                   "No normalidad       ", "Normalidad ok       ")))
cat(sprintf("║    Durbin-Watson :  DW=%.4f  p=%.5f  → %s      ║\n",
            dw_test$statistic, dw_test$p.value,
            ifelse(dw_test$p.value < 0.05,
                   "Autocorrelación     ", "Independencia ok    ")))
cat(sprintf("║    Ramsey RESET  :   F=%.4f  p=%.5f  → %s      ║\n",
            reset_test$statistic, reset_test$p.value,
            ifelse(reset_test$p.value < 0.05,
                   "Posible no lineal   ", "Linealidad ok       ")))

cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  OBSERVACIONES INFLUYENTES                                       ║\n")
cat(sprintf("║    Cook > 4/n:            %3d obs.  (umbral = %.4f)          ║\n",
            sum(cook > umbral_cook), umbral_cook))
cat(sprintf("║    Leverage > 2p/n:       %3d obs.  (umbral = %.4f)          ║\n",
            sum(lev > umbral_lev), umbral_lev))
cat(sprintf("║    |DFFITS| > 2√(p/n):   %3d obs.  (umbral = %.4f)          ║\n",
            sum(abs(df_its) > umbral_dffits), umbral_dffits))
cat(sprintf("║    Total influyentes:     %3d obs. de %d                      ║\n",
            length(obs_infl), n))

cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  VARIACIÓN EXPLICADA                                             ║\n")
cat(sprintf("║    SCT = %.2f  SCR = %.2f  SCE = %.2f             ║\n",
            SCT, SCR, SCE))
cat(sprintf("║    El modelo explica el %.1f%% de la variabilidad total de Y.   ║\n",
            smod$r.squared * 100))

cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║  CONCLUSIÓN GENERAL                                              ║\n")

if (f_pval < 0.05) {
  cat("║    El modelo es estadísticamente significativo en su conjunto.   ║\n")
} else {
  cat("║    El modelo NO es significativo globalmente (p ≥ 0.05).         ║\n")
}

# Predictores significativos
sig_preds <- rownames(smod$coefficients)[
  smod$coefficients[, 4] < 0.05 &
  rownames(smod$coefficients) != "(Intercept)"]
if (length(sig_preds) > 0) {
  cat(sprintf("║    Predictores significativos (p<.05): %s\n",
              paste(sig_preds, collapse = ", ")))
} else {
  cat("║    Ningún predictor es individualmente significativo (p<.05).   ║\n")
}

cat(sprintf("║    R² ajustado = %.4f: el modelo explica el %.1f%% de Y.      ║\n",
            smod$adj.r.squared, smod$adj.r.squared * 100))

cat(sprintf("║    Supuestos: %s | %s | %s | %s\n",
            ifelse(bp_test$p.value    >= 0.05, "Homo ✓", "Hetero ✗"),
            ifelse(sw_test$p.value    >= 0.05, "Normal ✓", "No-normal ✗"),
            ifelse(dw_test$p.value    >= 0.05, "Indep ✓", "Autocorr ✗"),
            ifelse(reset_test$p.value >= 0.05, "Lineal ✓", "No-lineal ✗")))

cat("╚══════════════════════════════════════════════════════════════════╝\n")

# =============================================================================
# REGRESIÓN LINEAL MÚLTIPLE SIN INTERCEPTO — ELIMINACIÓN HACIA ATRÁS
# Modelo: Y = β₁X₁ + β₂X₂ + ... + βₖXₖ + ε
# =============================================================================

# --- Paquetes ---
instalar_y_cargar <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
  library(p, character.only = TRUE)
}
instalar_y_cargar("lmtest")
instalar_y_cargar("car")

# =============================================================================
# SECCIÓN 1: LECTURA DE DATOS
# =============================================================================
datos <- read.csv(file.choose(), header = TRUE, sep = ",", dec = ".")

if (!"Y" %in% names(datos)) stop("El CSV debe tener una columna llamada 'Y'.")

regresores_iniciales <- setdiff(names(datos)[sapply(datos, is.numeric)], "Y")
if (length(regresores_iniciales) == 0) stop("No se encontraron regresores numéricos.")

n <- nrow(datos)

cat("Observaciones:", n, "| Regresores detectados:", paste(regresores_iniciales, collapse = ", "), "\n\n")

# =============================================================================
# SECCIÓN 2: MODELO COMPLETO SIN INTERCEPTO
# =============================================================================
# "Y ~ 0 +" suprime el intercepto; todos los regresores entran sin β₀.
formula_completa <- as.formula(paste("Y ~ 0 +", paste(regresores_iniciales, collapse = " + ")))

modelo_completo  <- lm(formula_completa, data = datos)
resumen_completo <- summary(modelo_completo)

cat("====================================================\n")
cat("  MODELO COMPLETO (sin intercepto)\n")
cat("====================================================\n")
print(resumen_completo)
cat("NOTA: R² no centrado (SST = ΣYᵢ²). Usar solo para comparar modelos sin intercepto.\n\n")

# =============================================================================
# SECCIÓN 3: ELIMINACIÓN HACIA ATRÁS (criterio: p-value > 0.05)
# =============================================================================
cat("====================================================\n")
cat("  ELIMINACIÓN HACIA ATRÁS — α = 0.05\n")
cat("====================================================\n\n")

alfa               <- 0.05
paso               <- 0
modelo_actual      <- modelo_completo
regresores         <- regresores_iniciales
variables_eliminadas <- c()

repeat {
  resumen_actual <- summary(modelo_actual)
  coef_tabla     <- resumen_actual$coefficients

  # Sin intercepto: todas las filas de coef_tabla son regresores.
  p_vals <- coef_tabla[, "Pr(>|t|)"]
  if (is.null(names(p_vals))) names(p_vals) <- rownames(coef_tabla)

  max_p    <- max(p_vals)
  var_elim <- names(which.max(p_vals))

  paso <- paso + 1
  if (paso > length(regresores_iniciales))
    stop("ERROR: Se superó el máximo de eliminaciones. Revise el ciclo.")
  if (!(var_elim %in% regresores))
    stop(paste("ERROR: Variable a eliminar no encontrada:", var_elim))

  cat(sprintf("--- PASO %d ---\n", paso))
  cat("Modelo actual: Y ~ 0 +", paste(regresores, collapse = " + "), "\n")
  print(coef_tabla)
  cat(sprintf("R² ajustado: %.6f  (no centrado)\n", resumen_actual$adj.r.squared))
  cat(sprintf("Mayor p-value: %s = %.6f\n", var_elim, max_p))

  if (max_p > alfa) {
    cat(sprintf("→ ELIMINAR '%s'  (p = %.6f > 0.05)\n\n", var_elim, max_p))
    variables_eliminadas <- c(variables_eliminadas, var_elim)
    regresores           <- setdiff(regresores, var_elim)
    modelo_actual        <- lm(
      as.formula(paste("Y ~ 0 +", paste(regresores, collapse = " + "))),
      data = datos
    )
  } else {
    cat(sprintf("→ TODOS los p-values ≤ 0.05. Proceso concluido.\n\n"))
    paso <- paso - 1
    break
  }

  if (length(regresores) == 0) { cat("Sin regresores restantes.\n\n"); break }
}

cat(sprintf("Pasos de eliminación realizados: %d\n\n", paso))

# =============================================================================
# SECCIÓN 4: MODELO FINAL
# =============================================================================
cat("====================================================\n")
cat("  MODELO FINAL (sin intercepto)\n")
cat("====================================================\n\n")

modelo_final  <- modelo_actual
resumen_final <- summary(modelo_final)
coef_finales  <- coef(modelo_final)

# Ecuación: Ŷ = β̂₁X₁ + β̂₂X₂ + ...  (sin término independiente)
terminos_eq <- sprintf("%.4f*%s", coef_finales[1], names(coef_finales)[1])
if (length(coef_finales) > 1) {
  for (i in 2:length(coef_finales)) {
    v <- coef_finales[i]
    terminos_eq <- paste0(terminos_eq,
      ifelse(v >= 0,
        sprintf(" + %.4f*%s", v, names(coef_finales)[i]),
        sprintf(" - %.4f*%s", abs(v), names(coef_finales)[i])))
  }
}
cat("Ecuación ajustada:\n  Ŷ =", terminos_eq, "\n\n")

cat("Tabla de coeficientes:\n")
print(resumen_final$coefficients)

# Medidas de ajuste
# df.residual() da gl = n - k (sin intercepto) automáticamente.
gl_error <- df.residual(modelo_final)
SSE      <- sum(residuals(modelo_final)^2)
MSE      <- SSE / gl_error
RMSE     <- sqrt(MSE)

cat(sprintf("\nR² (no centrado)  : %.6f\n", resumen_final$r.squared))
cat(sprintf("R² ajustado       : %.6f\n",   resumen_final$adj.r.squared))
cat(sprintf("SSE               : %.6f\n",   SSE))
cat(sprintf("MSE               : %.6f\n",   MSE))
cat(sprintf("RMSE              : %.6f\n",   RMSE))
cat(sprintf("gl error          : %d  (n - k = %d - %d)\n", gl_error, n, length(regresores)))

f_stat    <- resumen_final$fstatistic
p_valor_F <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
cat(sprintf("\nPrueba F: F(%d,%d) = %.4f,  p-value = %.6f\n", f_stat[2], f_stat[3], f_stat[1], p_valor_F))
cat(ifelse(p_valor_F <= alfa,
  "→ Modelo GLOBALMENTE SIGNIFICATIVO (p ≤ 0.05)\n",
  "→ Modelo NO globalmente significativo (p > 0.05)\n"))

cat("\nVariables en el modelo final :", paste(regresores, collapse = ", "), "\n")
cat("Variables eliminadas          :",
  ifelse(length(variables_eliminadas) > 0, paste(variables_eliminadas, collapse = ", "), "ninguna"), "\n\n")

# =============================================================================
# SECCIÓN 5: ANÁLISIS RESIDUAL — GRÁFICAS
# =============================================================================
residuos  <- residuals(modelo_final)
ajustados <- fitted(modelo_final)

total_graficas <- 4 + length(regresores)
par(mfrow = c(ceiling(total_graficas / 2), 2), mar = c(4, 4, 3, 1))

plot(ajustados, datos$Y, main = "Observados vs Ajustados",
     xlab = "Ŷᵢ", ylab = "Yᵢ", pch = 19, col = "steelblue")
abline(0, 1, col = "red", lwd = 2, lty = 2)

plot(ajustados, residuos, main = "Residuos vs Ajustados",
     xlab = "Ŷᵢ", ylab = "eᵢ", pch = 19, col = "darkorange")
abline(h = 0, col = "red", lwd = 2, lty = 2)

qqnorm(residuos, main = "QQ-Plot de Residuos", pch = 19, col = "purple")
qqline(residuos, col = "red", lwd = 2)

hist(residuos, main = "Histograma de Residuos", xlab = "eᵢ",
     col = "lightblue", border = "white", breaks = "Sturges")
xseq <- seq(min(residuos), max(residuos), length.out = 100)
lines(xseq, dnorm(xseq, mean(residuos), sd(residuos)) *
        length(residuos) * diff(hist(residuos, plot = FALSE)$breaks)[1],
      col = "red", lwd = 2)

for (var in regresores) {
  plot(datos[[var]], residuos, main = paste("Residuos vs", var),
       xlab = var, ylab = "eᵢ", pch = 19, col = "darkgreen")
  abline(h = 0, col = "red", lwd = 2, lty = 2)
}
par(mfrow = c(1, 1))

# =============================================================================
# SECCIÓN 6: NORMALIDAD — SHAPIRO-WILK
# =============================================================================
cat("====================================================\n")
cat("  NORMALIDAD (Shapiro-Wilk)\n")
cat("  H₀: residuos ~ Normal  |  H₁: residuos ≁ Normal\n")
cat("====================================================\n")
prueba_shapiro <- shapiro.test(residuos)
print(prueba_shapiro)
cat(ifelse(prueba_shapiro$p.value > alfa,
  sprintf("→ No se rechaza H₀ (p = %.6f > 0.05). Se puede ASUMIR NORMALIDAD.\n\n",
          prueba_shapiro$p.value),
  sprintf("→ Se rechaza H₀ (p = %.6f ≤ 0.05). NO se puede asumir normalidad.\n\n",
          prueba_shapiro$p.value)))

# =============================================================================
# SECCIÓN 7: HOMOCEDASTICIDAD — BREUSCH-PAGAN
# =============================================================================
cat("====================================================\n")
cat("  HOMOCEDASTICIDAD (Breusch-Pagan)\n")
cat("  H₀: varianza constante  |  H₁: heterocedasticidad\n")
cat("====================================================\n")
prueba_bp <- bptest(modelo_final)
print(prueba_bp)
cat(ifelse(prueba_bp$p.value > alfa,
  sprintf("→ No se rechaza H₀ (p = %.6f > 0.05). Se puede ASUMIR HOMOCEDASTICIDAD.\n\n",
          prueba_bp$p.value),
  sprintf("→ Se rechaza H₀ (p = %.6f ≤ 0.05). NO se puede asumir homocedasticidad.\n\n",
          prueba_bp$p.value)))

# =============================================================================
# SECCIÓN 8: INDEPENDENCIA — DURBIN-WATSON
# =============================================================================
cat("====================================================\n")
cat("  INDEPENDENCIA (Durbin-Watson)\n")
cat("  H₀: sin autocorrelación  |  H₁: autocorrelación\n")
cat("====================================================\n")
prueba_dw <- dwtest(modelo_final)
print(prueba_dw)
cat(ifelse(prueba_dw$p.value > alfa,
  sprintf("→ No se rechaza H₀ (p = %.6f > 0.05). Se puede ASUMIR INDEPENDENCIA.\n\n",
          prueba_dw$p.value),
  sprintf("→ Se rechaza H₀ (p = %.6f ≤ 0.05). NO se puede asumir independencia.\n\n",
          prueba_dw$p.value)))

# =============================================================================
# SECCIÓN 9: MULTICOLINEALIDAD — VIF
# =============================================================================
cat("====================================================\n")
cat("  MULTICOLINEALIDAD (VIF)\n")
cat("  VIF ≤ 5: sin problema  |  5-10: moderada  |  >10: fuerte\n")
cat("====================================================\n")
if (length(regresores) >= 2) {
  vif_valores <- vif(modelo_final)
  print(vif_valores)
  if (any(vif_valores > 10)) {
    cat("→ Multicolinealidad FUERTE (VIF > 10) en:",
        paste(names(vif_valores)[vif_valores > 10], collapse = ", "), "\n\n")
  } else if (any(vif_valores > 5)) {
    cat("→ Multicolinealidad MODERADA (VIF > 5) en:",
        paste(names(vif_valores)[vif_valores > 5], collapse = ", "), "\n\n")
  } else {
    cat("→ Sin multicolinealidad significativa (todos los VIF ≤ 5).\n\n")
  }
} else {
  vif_valores <- NULL
  cat("→ Un solo regresor: VIF no aplica.\n\n")
}

# =============================================================================
# SECCIÓN 10: OBSERVACIONES INFLUYENTES — DISTANCIA DE COOK
# =============================================================================
cat("====================================================\n")
cat("  OBSERVACIONES INFLUYENTES (Distancia de Cook)\n")
cat(sprintf("  Umbral: 4/n = 4/%d = %.6f\n", n, 4/n))
cat("====================================================\n")

cook_dist   <- cooks.distance(modelo_final)
umbral_cook <- 4 / n

plot(cook_dist, type = "h",
     main = "Distancia de Cook (sin intercepto)",
     xlab = "Observación", ylab = "Dᵢ",
     col = ifelse(cook_dist > umbral_cook, "red", "steelblue"), lwd = 1.5)
abline(h = umbral_cook, col = "red", lwd = 2, lty = 2)

obs_influyentes <- which(cook_dist > umbral_cook)
if (length(obs_influyentes) > 0) {
  cat(sprintf("→ Observaciones influyentes (Dᵢ > %.6f): %s\n\n",
              umbral_cook, paste(obs_influyentes, collapse = ", ")))
  print(cook_dist[obs_influyentes])
} else {
  cat(sprintf("→ Ninguna observación supera el umbral (%.6f).\n\n", umbral_cook))
}

# =============================================================================
# SECCIÓN 11: CONCLUSIÓN GENERAL
# =============================================================================
multi_str <- if (!is.null(vif_valores)) {
  if (any(vif_valores > 10)) "multicolinealidad FUERTE detectada"
  else if (any(vif_valores > 5)) "multicolinealidad MODERADA detectada"
  else "sin multicolinealidad significativa"
} else "VIF no aplica (un regresor)"

inf_str <- ifelse(length(obs_influyentes) > 0,
  paste("observaciones influyentes:", paste(obs_influyentes, collapse = ", ")),
  "sin observaciones influyentes")

cat("====================================================\n")
cat("  CONCLUSIÓN GENERAL\n")
cat("  MODELO SIN INTERCEPTO — Eliminación Hacia Atrás\n")
cat("====================================================\n\n")
cat(sprintf(paste0(
  "Ecuación final  : Ŷ = %s\n",
  "Variables final : %s\n",
  "Eliminadas      : %s\n",
  "R² ajustado*    : %.4f  (* no centrado — interpretar con cautela)\n",
  "Prueba F global : F(%d,%d) = %.4f, p = %.6f  →  %s\n",
  "Normalidad      : Shapiro-Wilk p = %.6f  →  %s\n",
  "Homocedasticidad: Breusch-Pagan p = %.6f  →  %s\n",
  "Independencia   : Durbin-Watson p = %.6f  →  %s\n",
  "Multicolinealidad: %s\n",
  "Cook            : %s\n\n",
  "INTERPRETACIÓN:\n",
  "El modelo SIN INTERCEPTO asume E(Y) = 0 cuando todos los regresores\n",
  "son cero (plano de regresión obligado a pasar por el origen).\n",
  "Use este modelo solo si esa restricción tiene sentido teórico.\n",
  "Cada β̂ representa el cambio promedio en Y por unidad del regresor,\n",
  "manteniendo constantes los demás regresores del modelo.\n"
),
terminos_eq,
paste(regresores, collapse = ", "),
ifelse(length(variables_eliminadas) > 0, paste(variables_eliminadas, collapse = ", "), "ninguna"),
round(resumen_final$adj.r.squared, 4),
f_stat[2], f_stat[3], f_stat[1], p_valor_F,
ifelse(p_valor_F <= alfa, "SIGNIFICATIVO", "no significativo"),
prueba_shapiro$p.value,
ifelse(prueba_shapiro$p.value > alfa, "normalidad asumida", "normalidad NO asumida"),
prueba_bp$p.value,
ifelse(prueba_bp$p.value > alfa, "homocedasticidad asumida", "homocedasticidad NO asumida"),
prueba_dw$p.value,
ifelse(prueba_dw$p.value > alfa, "independencia asumida", "independencia NO asumida"),
multi_str,
inf_str
))

cat("====================================================\n")
cat("  FIN DEL ANÁLISIS\n")
cat("====================================================\n")


# =============================================================================
# =============================================================================
# GUÍA DE INTERPRETACIÓN DE RESULTADOS PARA EL REPORTE
# (MODELO SIN INTERCEPTO)
# =============================================================================
# =============================================================================
#
# Adapte los valores numéricos entre corchetes [] según sus resultados.
#
# ----------------------------------------------------------------------------
# 1. ELIMINACIÓN HACIA ATRÁS (SIN INTERCEPTO)
# ----------------------------------------------------------------------------
# "Se aplicó el método de eliminación hacia atrás partiendo del modelo
#  completo sin intercepto:  Y = β₁X₁ + β₂X₂ + ... + βₖXₖ + ε.
#  En el Paso 1, el regresor [X_j] presentó el mayor p-value
#  (p = [valor] > 0.05), por lo que fue eliminado. Este proceso se repitió
#  hasta que todos los regresores restantes resultaron estadísticamente
#  significativos al nivel α = 0.05. En total se realizaron [m] pasos de
#  eliminación, obteniéndose un modelo final con [q] regresores."
#
# ----------------------------------------------------------------------------
# 2. MODELO FINAL (SIN INTERCEPTO)
# ----------------------------------------------------------------------------
# "El modelo final obtenido es:
#      Ŷ = β̂₁X₁ + β̂₂X₂
#  donde β̂₁ = [valor] representa el cambio promedio en Y por cada unidad
#  de aumento en X₁, manteniendo constante X₂, y β̂₂ = [valor] representa
#  el cambio promedio en Y por cada unidad de aumento en X₂, manteniendo
#  constante X₁. El modelo no incluye intercepto, por lo que se asume que
#  cuando todos los regresores son iguales a cero, Y = 0."
#
# ----------------------------------------------------------------------------
# 3. R² AJUSTADO EN MODELO SIN INTERCEPTO
# ----------------------------------------------------------------------------
# "El coeficiente de determinación ajustado del modelo final fue
#  R²_adj = [valor]. Dado que el modelo fue ajustado sin intercepto, el R²
#  reportado por R corresponde a un R² no centrado, calculado con
#  SST = Σ(Yᵢ²) en lugar de Σ(Yᵢ - Ȳ)². Por ello, este valor puede ser
#  artificialmente elevado y debe interpretarse con cautela, utilizándose
#  principalmente para comparar modelos sin intercepto entre sí."
#
# ----------------------------------------------------------------------------
# 4. SIGNIFICANCIA INDIVIDUAL DE LOS REGRESORES
# ----------------------------------------------------------------------------
# "De acuerdo con la tabla de coeficientes del modelo final, el regresor X₁
#  resultó estadísticamente significativo (t = [valor], p-value = [valor] <
#  0.05), lo que indica que contribuye significativamente a explicar la
#  variabilidad de Y, dado el resto de regresores en el modelo."
#
# ----------------------------------------------------------------------------
# 5. PRUEBA F GLOBAL
# ----------------------------------------------------------------------------
# "La prueba F global del modelo final arrojó F([gl1], [gl2]) = [valor],
#  con p-value = [valor] < 0.05. Se rechaza la hipótesis nula de que todos
#  los coeficientes son simultáneamente iguales a cero, por lo que el modelo
#  en su conjunto es estadísticamente significativo."
#
# ----------------------------------------------------------------------------
# 6. NORMALIDAD (Shapiro-Wilk)
# ----------------------------------------------------------------------------
# "El QQ-plot de los residuos muestra que los puntos se distribuyen
#  aproximadamente sobre la línea diagonal de referencia. La prueba de
#  Shapiro-Wilk arrojó W = [valor], p-value = [valor]. Dado que el p-value
#  [es mayor / es menor] que α = 0.05, [no se rechaza / se rechaza] H₀.
#  Por lo tanto, [se puede / no se puede] asumir normalidad de los residuos."
#
# ----------------------------------------------------------------------------
# 7. HOMOCEDASTICIDAD (Breusch-Pagan)
# ----------------------------------------------------------------------------
# "La gráfica de residuos vs valores ajustados no muestra un patrón
#  sistemático en la dispersión. La prueba de Breusch-Pagan arrojó
#  BP = [valor], p-value = [valor]. Dado que el p-value [es mayor / es menor]
#  que α = 0.05, [no se rechaza / se rechaza] H₀. En consecuencia,
#  [se puede / no se puede] asumir homocedasticidad."
#
# ----------------------------------------------------------------------------
# 8. INDEPENDENCIA (Durbin-Watson)
# ----------------------------------------------------------------------------
# "La prueba de Durbin-Watson arrojó DW = [valor], p-value = [valor].
#  Dado que el p-value [es mayor / es menor] que α = 0.05,
#  [no se rechaza / se rechaza] H₀. Por lo tanto, [se puede / no se puede]
#  asumir que los errores del modelo son independientes entre sí."
#
# ----------------------------------------------------------------------------
# 9. MULTICOLINEALIDAD (VIF)
# ----------------------------------------------------------------------------
# "Se calculó el VIF para cada regresor del modelo final, obteniendo
#  [X₁: VIF = valor, X₂: VIF = valor]. Dado que [todos los valores son ≤ 5 /
#  algunos valores superan 5 o 10], [no se detecta multicolinealidad
#  significativa / se detecta multicolinealidad moderada o fuerte] entre los
#  regresores del modelo final."
#
# ----------------------------------------------------------------------------
# 10. OBSERVACIONES INFLUYENTES (Distancia de Cook)
# ----------------------------------------------------------------------------
# "Se calculó la distancia de Cook para cada observación, con umbral
#  Dᵢ > 4/n = [umbral]. [No se identificó ninguna observación influyente.]
#  [Alternativamente:] Las observaciones [i₁, i₂, ...] presentaron
#  Dᵢ > [umbral], por lo que podrían ejercer una influencia considerable
#  sobre los coeficientes estimados. Se recomienda investigar si son errores
#  de medición, valores atípicos legítimos o puntos extremos en el espacio
#  de los regresores."
#
# =============================================================================
# FIN DE LA GUÍA DE INTERPRETACIÓN
# =============================================================================


# =============================================================================
# =============================================================================
# RESUMEN DE CAMBIOS RESPECTO AL CÓDIGO CON INTERCEPTO
# =============================================================================
# =============================================================================
#
# 1. FÓRMULA "Y ~ 0 + ..." EN LUGAR DE "Y ~ ..."
#    Al escribir "Y ~ 0 +" se le indica explícitamente a R que elimine el
#    término (Intercept) del modelo. La notación alternativa "Y ~ X - 1"
#    produce el mismo resultado matemático, pero "0 +" es más clara e intuitiva
#    para el lector. Este cambio afecta la fórmula del modelo completo (Sección 2)
#    y todas las fórmulas dentro del ciclo de eliminación (Sección 3).
#
# 2. EXTRACCIÓN DE P-VALUES SIN EXCLUIR "(Intercept)"
#    En el modelo con intercepto, la tabla de coeficientes incluye una fila
#    "(Intercept)" que no es un regresor y no debe participar en la selección
#    de variables a eliminar, de ahí el filtro:
#      coef_tabla[rownames(coef_tabla) != "(Intercept)", 4]
#    En el modelo sin intercepto, R simplemente no genera esa fila. Filtrar
#    por "(Intercept)" no causaría error, pero sería código innecesario y
#    confuso. Se usa directamente:
#      coef_tabla[, "Pr(>|t|)"]
#    que toma todos los p-values de la tabla, que en este caso son solo los
#    de los regresores.
#
# 3. ECUACIÓN FINAL SIN TÉRMINO INDEPENDIENTE
#    El código original construía la ecuación iniciando con el intercepto:
#      terminos_eq <- sprintf("%.4f", coef_finales[1])   # β₀
#      for (i in 2:length(coef_finales)) { ... }          # β₁X₁, β₂X₂, ...
#    En el modelo sin intercepto, coef_finales[1] es el coeficiente del
#    primer regresor (β̂₁), no un intercepto. Por ello la ecuación comienza
#    directamente con el primer término:
#      terminos_eq <- sprintf("%.4f*%s", coef_finales[1], names(coef_finales)[1])
#      for (i in 2:length(coef_finales)) { ... }
#    Esto garantiza que la ecuación impresa tenga la forma correcta:
#      Ŷ = β̂₁X₁ + β̂₂X₂ + ...   (sin β̂₀)
#
# 4. GRADOS DE LIBERTAD CON df.residual()
#    En un modelo CON intercepto se estiman k+1 parámetros (k regresores
#    + 1 intercepto), por lo que gl_error = n - k - 1.
#    En un modelo SIN intercepto solo se estiman k parámetros, por lo que
#    gl_error = n - k.
#    Usar df.residual(modelo_final) hace que R calcule automáticamente los
#    grados de libertad correctos para cualquier tipo de modelo, evitando
#    el error de restar 1 de más cuando no hay intercepto.
#
# 5. R² NO CENTRADO EN MODELOS SIN INTERCEPTO
#    Cuando se ajusta un modelo sin intercepto, R calcula el R² como:
#      R² = 1 - SSE / SST,  donde  SST = Σ(Yᵢ²)
#    en lugar del R² centrado habitual:
#      R² = 1 - SSE / SST,  donde  SST = Σ(Yᵢ - Ȳ)²
#    El SST no centrado (Σ(Yᵢ²)) es siempre ≥ Σ(Yᵢ - Ȳ)², lo que hace que
#    el cociente SSE/SST sea más pequeño y el R² resultante sea artificialmente
#    mayor. Esto puede llevar a pensar que el modelo ajusta muy bien cuando en
#    realidad solo "se beneficia" de no centrar la variabilidad total. Por ello
#    el R² de un modelo sin intercepto NO debe compararse con el de un modelo
#    con intercepto: son métricas con denominadores distintos.
#
# =============================================================================
# FIN DEL SCRIPT
# =============================================================================

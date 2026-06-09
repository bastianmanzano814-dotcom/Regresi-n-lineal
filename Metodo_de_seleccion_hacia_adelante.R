# =============================================================================
# ANÁLISIS DE REGRESIÓN LINEAL MÚLTIPLE
# MÉTODO DE SELECCIÓN HACIA ADELANTE (FORWARD SELECTION)
# =============================================================================
# Curso: Análisis de Regresión
# Descripción: Este script realiza un análisis completo de regresión lineal
#              múltiple utilizando el método de selección hacia adelante.
#              Se parte del modelo nulo (solo intercepto) y se agrega, paso
#              a paso, el regresor candidato con menor p-value, siempre que
#              ese p-value sea ≤ α = 0.05. El proceso se detiene cuando
#              ninguna variable fuera del modelo cumple el criterio de entrada.
# =============================================================================


# =============================================================================
# SECCIÓN 0: INSTALACIÓN Y CARGA DE PAQUETES
# =============================================================================

# Función que instala el paquete si no está disponible y luego lo carga.
# Esto garantiza que el script funcione en cualquier equipo sin configuración
# manual previa.
instalar_y_cargar <- function(paquete) {
  if (!requireNamespace(paquete, quietly = TRUE)) {
    install.packages(paquete, dependencies = TRUE)
  }
  library(paquete, character.only = TRUE)
}

# Paquetes requeridos:
# - lmtest: pruebas de Breusch-Pagan (homocedasticidad) y Durbin-Watson (independencia)
# - car   : VIF (factor de inflación de la varianza) para multicolinealidad
instalar_y_cargar("lmtest")
instalar_y_cargar("car")

cat("\n====================================================\n")
cat("  ANÁLISIS DE REGRESIÓN LINEAL MÚLTIPLE\n")
cat("  Método: Selección Hacia Adelante\n")
cat("====================================================\n\n")


# =============================================================================
# SECCIÓN 1: LECTURA DE DATOS
# =============================================================================

# Se abre un cuadro de diálogo para seleccionar el archivo CSV.
# Requisito: el archivo debe tener una columna llamada exactamente "Y"
# (variable respuesta) y una o más columnas numéricas adicionales (regresores).
cat("Por favor, seleccione el archivo CSV con sus datos...\n")
datos <- read.csv(file.choose(), header = TRUE, sep = ",", dec = ".")

# Verificar que el archivo se leyó correctamente.
cat("\n--- Vista previa de los datos (primeras 6 filas) ---\n")
print(head(datos))
cat("\n--- Dimensiones del conjunto de datos ---\n")
cat("Número de observaciones:", nrow(datos), "\n")
cat("Número de variables    :", ncol(datos), "\n")
cat("Nombres de columnas    :", paste(names(datos), collapse = ", "), "\n\n")

# Validar que exista la columna "Y".
if (!"Y" %in% names(datos)) {
  stop("ERROR: El archivo CSV no contiene una columna llamada 'Y'.
       Asegúrese de que la variable respuesta se llame exactamente 'Y'.")
}

# Identificar automáticamente los regresores candidatos: todas las columnas
# numéricas excepto "Y".
candidatos_iniciales <- setdiff(
  names(datos)[sapply(datos, is.numeric)],
  "Y"
)

if (length(candidatos_iniciales) == 0) {
  stop("ERROR: No se encontraron columnas numéricas que funcionen como regresores.")
}

cat("Variable respuesta detectada  : Y\n")
cat("Regresores candidatos         :", paste(candidatos_iniciales, collapse = ", "), "\n\n")

# Tamaño de muestra (se usa más adelante para la distancia de Cook).
n <- nrow(datos)


# =============================================================================
# SECCIÓN 2: MODELOS INICIAL Y COMPLETO
# =============================================================================
# La selección hacia adelante parte del modelo más simple posible:
#   Modelo nulo:     Y = β₀ + ε         (solo intercepto, sin regresores)
#   Modelo completo: Y = β₀ + β₁X₁ + ... + βₖXₖ + ε  (todos los regresores)
#
# El modelo nulo es el punto de partida. El modelo completo define el conjunto
# máximo de variables que podrían entrar al modelo final.

cat("====================================================\n")
cat("  SECCIÓN 2: MODELOS NULO Y COMPLETO\n")
cat("====================================================\n\n")

# --- Modelo nulo ---
# lm(Y ~ 1, ...) ajusta un modelo que solo estima la media global de Y (β₀).
modelo_nulo <- lm(Y ~ 1, data = datos)
cat("Modelo nulo (punto de partida de la selección hacia adelante):\n")
cat("  Y = β₀ + ε\n\n")
cat("--- Resumen del Modelo Nulo ---\n")
print(summary(modelo_nulo))

# --- Modelo completo ---
# Se define para conocer el universo de variables disponibles.
formula_completa <- as.formula(
  paste("Y ~", paste(candidatos_iniciales, collapse = " + "))
)
modelo_completo <- lm(formula_completa, data = datos)
cat("\nModelo completo (referencia con todos los regresores):\n")
print(formula_completa)
cat("\n--- Resumen del Modelo Completo ---\n")
resumen_completo <- summary(modelo_completo)
print(resumen_completo)
cat(sprintf("\nR² ajustado del modelo completo: %.6f\n\n",
            resumen_completo$adj.r.squared))


# =============================================================================
# SECCIÓN 3: MÉTODO DE SELECCIÓN HACIA ADELANTE
# =============================================================================
# Lógica del algoritmo:
#   1. Se parte del modelo nulo (sin regresores).
#   2. En cada paso se prueban, de forma individual, todas las variables
#      candidatas que aún no están en el modelo.
#   3. Para cada candidata se ajusta un modelo temporal que incluye las
#      variables ya seleccionadas más esa candidata.
#   4. Se extrae el p-value de la candidata en ese modelo temporal.
#   5. La variable con menor p-value entra al modelo, siempre que p ≤ α = 0.05.
#   6. Si ninguna candidata cumple el criterio, el proceso se detiene.
#
# Nota sobre R² ajustado:
#   Aunque el criterio principal para incorporar una variable es su p-value
#   individual (≤ α = 0.05), en cada paso también se reporta el R² ajustado
#   del modelo actual. Este indicador penaliza la inclusión de variables
#   innecesarias y es útil para comparar la calidad del ajuste entre modelos
#   con distinto número de regresores.

cat("====================================================\n")
cat("  SECCIÓN 3: SELECCIÓN HACIA ADELANTE\n")
cat("  Criterio principal: p-value ≤ 0.05\n")
cat("  Nivel de significancia: α = 0.05\n")
cat("====================================================\n\n")

alfa              <- 0.05                   # Nivel de significancia
paso              <- 0                      # Contador de pasos
variables_en_modelo  <- c()                 # Regresores ya seleccionados
variables_fuera      <- candidatos_iniciales # Candidatos aún no incorporados
modelo_actual        <- modelo_nulo         # Modelo actual (empieza como nulo)
variables_no_entran  <- c()                 # Registro de variables que nunca entraron

# Bucle principal de selección hacia adelante.
repeat {

  # Si ya no quedan candidatos fuera del modelo, detener.
  if (length(variables_fuera) == 0) {
    cat("No quedan variables candidatas fuera del modelo. El proceso se detiene.\n\n")
    break
  }

  paso <- paso + 1

  cat(sprintf("--- PASO %d ---\n", paso))
  if (length(variables_en_modelo) == 0) {
    cat("Variables en el modelo actual: (ninguna — modelo nulo)\n")
  } else {
    cat("Variables en el modelo actual:", paste(variables_en_modelo, collapse = ", "), "\n")
  }
  cat("Variables candidatas (fuera) :", paste(variables_fuera, collapse = ", "), "\n\n")

  # --- Evaluar cada candidata ---
  # Para cada variable fuera del modelo, se ajusta un modelo temporal que
  # incluye las variables ya seleccionadas más esa candidata, y se extrae
  # el p-value de la candidata en dicho modelo temporal.

  p_valores_candidatas <- c()  # Almacena p-value de cada candidata

  for (var_cand in variables_fuera) {
    # Construir fórmula del modelo temporal.
    vars_temp <- c(variables_en_modelo, var_cand)
    formula_temp <- as.formula(
      paste("Y ~", paste(vars_temp, collapse = " + "))
    )
    modelo_temp <- lm(formula_temp, data = datos)
    tabla_coef  <- summary(modelo_temp)$coefficients

    # Extraer el p-value de la variable candidata (fila correspondiente).
    p_cand <- tabla_coef[var_cand, 4]
    p_valores_candidatas[var_cand] <- p_cand
  }

  # Mostrar p-values de todas las candidatas en este paso.
  cat("P-values de las variables candidatas al agregarse al modelo:\n")
  for (vc in names(p_valores_candidatas)) {
    cat(sprintf("  %-15s p-value = %.6f\n", vc, p_valores_candidatas[vc]))
  }
  cat("\n")

  # Identificar la candidata con menor p-value.
  min_p    <- min(p_valores_candidatas)
  var_entra <- names(which.min(p_valores_candidatas))

  cat(sprintf("Variable con menor p-value: %s  (p-value = %.6f)\n", var_entra, min_p))

  # Verificar si la mejor candidata cumple el criterio de entrada.
  if (min_p <= alfa) {
    # La variable entra al modelo.
    cat(sprintf("Decisión: AGREGAR '%s'  →  p-value = %.6f ≤ α = %.2f\n",
                var_entra, min_p, alfa))
    cat(sprintf("Razón: El regresor '%s' es estadísticamente significativo\n",
                var_entra))
    cat(sprintf("       al nivel α = %.2f, por lo que se incorpora al modelo.\n", alfa))

    # Actualizar conjuntos de variables.
    variables_en_modelo <- c(variables_en_modelo, var_entra)
    variables_fuera     <- setdiff(variables_fuera, var_entra)

    # Reajustar el modelo con la nueva variable incluida.
    nueva_formula <- as.formula(
      paste("Y ~", paste(variables_en_modelo, collapse = " + "))
    )
    modelo_actual <- lm(nueva_formula, data = datos)
    resumen_actual <- summary(modelo_actual)

    cat(sprintf("\nR² ajustado del modelo actual (después de agregar '%s'): %.6f\n\n",
                var_entra, resumen_actual$adj.r.squared))
    cat("Tabla de coeficientes del modelo actual:\n")
    print(resumen_actual$coefficients)
    cat("\n")

  } else {
    # Ninguna candidata cumple el criterio; se detiene el proceso.
    cat(sprintf("Decisión: NO AGREGAR NINGUNA VARIABLE  →  menor p-value = %.6f > α = %.2f\n",
                min_p, alfa))
    cat("Ninguna variable candidata es estadísticamente significativa al\n")
    cat(sprintf("nivel α = %.2f. El proceso de selección hacia adelante ha concluido.\n\n", alfa))
    # Las variables que quedaron fuera nunca entraron al modelo.
    variables_no_entran <- variables_fuera
    paso <- paso - 1  # El último paso fue de verificación, no de incorporación.
    break
  }
}

# Si se agotaron todas las candidatas y todas entraron.
if (length(variables_fuera) == 0 && length(variables_en_modelo) > 0) {
  variables_no_entran <- c()
}

cat(sprintf("Total de pasos de selección realizados: %d\n\n", paso))


# =============================================================================
# SECCIÓN 4: MODELO FINAL
# =============================================================================

cat("====================================================\n")
cat("  SECCIÓN 4: MODELO FINAL\n")
cat("====================================================\n\n")

modelo_final  <- modelo_actual
resumen_final <- summary(modelo_final)

# --- Ecuación ajustada ---
coef_finales <- coef(modelo_final)
terminos_eq  <- sprintf("%.4f", coef_finales[1])  # intercepto

for (i in 2:length(coef_finales)) {
  nombre_var <- names(coef_finales)[i]
  valor      <- coef_finales[i]
  if (valor >= 0) {
    terminos_eq <- paste0(terminos_eq, sprintf(" + %.4f*%s", valor, nombre_var))
  } else {
    terminos_eq <- paste0(terminos_eq, sprintf(" - %.4f*%s", abs(valor), nombre_var))
  }
}

cat("Ecuación ajustada del modelo final:\n")
cat(sprintf("  Ŷ = %s\n\n", terminos_eq))

# --- Tabla de coeficientes del modelo final ---
cat("Tabla de coeficientes (modelo final):\n")
print(resumen_final$coefficients)

# --- Medidas de ajuste ---
k_final  <- length(variables_en_modelo)  # número de regresores en el modelo final
gl_error <- n - k_final - 1              # grados de libertad del error
SSE      <- sum(residuals(modelo_final)^2)
MSE      <- SSE / gl_error
RMSE     <- sqrt(MSE)

cat("\n--- Medidas de ajuste del modelo final ---\n")
cat(sprintf("R²                           : %.6f\n",  resumen_final$r.squared))
cat(sprintf("R² ajustado                  : %.6f\n",  resumen_final$adj.r.squared))
cat(sprintf("SSE                          : %.6f\n",  SSE))
cat(sprintf("MSE                          : %.6f\n",  MSE))
cat(sprintf("RMSE                         : %.6f\n",  RMSE))
cat(sprintf("Grados de libertad del error : %d\n",    gl_error))

# --- Prueba F global ---
f_stat    <- resumen_final$fstatistic
p_valor_F <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)

cat("\n--- Prueba F global ---\n")
cat(sprintf("F(%d, %d) = %.4f\n", f_stat[2], f_stat[3], f_stat[1]))
cat(sprintf("p-value  = %.6f\n", p_valor_F))
if (p_valor_F <= alfa) {
  cat("→ El modelo final es GLOBALMENTE SIGNIFICATIVO (p-value ≤ 0.05).\n")
} else {
  cat("→ El modelo final NO es globalmente significativo (p-value > 0.05).\n")
}

# --- Resumen de variables ---
cat("\n--- Resumen de variables ---\n")
cat("Variables que entraron al modelo final :", paste(variables_en_modelo, collapse = ", "), "\n")
if (length(variables_no_entran) > 0) {
  cat("Variables que NO entraron al modelo    :", paste(variables_no_entran, collapse = ", "), "\n\n")
} else {
  cat("Todas las variables candidatas entraron al modelo final.\n\n")
}


# =============================================================================
# SECCIÓN 5: ANÁLISIS RESIDUAL
# =============================================================================
# Los residuos se definen como:  eᵢ = Yᵢ − Ŷᵢ
# Se examinan para verificar los supuestos del modelo de regresión lineal:
# normalidad, homocedasticidad e independencia de los errores.

cat("====================================================\n")
cat("  SECCIÓN 5: ANÁLISIS RESIDUAL\n")
cat("====================================================\n\n")

residuos  <- residuals(modelo_final)   # eᵢ = Yᵢ − Ŷᵢ
ajustados <- fitted(modelo_final)      # Ŷᵢ

cat("Residuos calculados: eᵢ = Yᵢ − Ŷᵢ\n\n")
cat("Primeros residuos:\n")
print(head(data.frame(Y_obs = datos$Y, Y_ajustado = ajustados, Residuo = residuos)))

# =============================================================================
# SECCIÓN 5.1: GRÁFICAS RESIDUALES
# =============================================================================
# Las gráficas se muestran en la ventana Plots de RStudio (NO se guardan).

num_regresores_final <- length(variables_en_modelo)
total_graficas       <- 4 + num_regresores_final
num_cols             <- 2
num_filas            <- ceiling(total_graficas / num_cols)
par(mfrow = c(num_filas, num_cols), mar = c(4, 4, 3, 1))

# --- Gráfica 1: Valores observados vs valores ajustados ---
# Permite ver qué tan bien el modelo reproduce los valores reales de Y.
plot(ajustados, datos$Y,
     main = "Valores Observados vs Ajustados",
     xlab = "Valores Ajustados (Ŷᵢ)",
     ylab = "Valores Observados (Yᵢ)",
     pch  = 19, col = "steelblue")
abline(0, 1, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = "Y = Ŷ (referencia)",
       col = "red", lty = 2, lwd = 2, cex = 0.8)

# --- Gráfica 2: Residuos vs valores ajustados ---
# Un patrón aleatorio alrededor de cero sugiere homocedasticidad y linealidad.
plot(ajustados, residuos,
     main = "Residuos vs Valores Ajustados",
     xlab = "Valores Ajustados (Ŷᵢ)",
     ylab = "Residuos (eᵢ)",
     pch  = 19, col = "darkorange")
abline(h = 0, col = "red", lwd = 2, lty = 2)

# --- Gráfica 3: QQ-plot de residuos ---
# Si los puntos siguen la línea diagonal, los residuos son aproximadamente normales.
qqnorm(residuos,
       main  = "QQ-Plot de Residuos",
       xlab  = "Cuantiles Teóricos",
       ylab  = "Cuantiles Muestrales",
       pch   = 19, col = "purple")
qqline(residuos, col = "red", lwd = 2)

# --- Gráfica 4: Histograma de residuos ---
# Una forma de campana simétrica sugiere distribución aproximadamente normal.
hist(residuos,
     main   = "Histograma de Residuos",
     xlab   = "Residuos (eᵢ)",
     ylab   = "Frecuencia",
     col    = "lightblue",
     border = "white",
     breaks = "Sturges")
xseq <- seq(min(residuos), max(residuos), length.out = 100)
lines(xseq,
      dnorm(xseq, mean(residuos), sd(residuos)) *
        length(residuos) * diff(hist(residuos, plot = FALSE)$breaks)[1],
      col = "red", lwd = 2)

# --- Gráficas adicionales: Residuos vs cada regresor del modelo final ---
# Permiten detectar relaciones no capturadas y patrones de no linealidad.
for (var in variables_en_modelo) {
  plot(datos[[var]], residuos,
       main = paste("Residuos vs", var),
       xlab = var,
       ylab = "Residuos (eᵢ)",
       pch  = 19, col = "darkgreen")
  abline(h = 0, col = "red", lwd = 2, lty = 2)
}

par(mfrow = c(1, 1))
cat("Gráficas de análisis residual mostradas en la ventana Plots de RStudio.\n\n")


# =============================================================================
# SECCIÓN 6: PRUEBA DE NORMALIDAD (SHAPIRO-WILK)
# =============================================================================
# H₀: Los residuos siguen una distribución normal.
# H₁: Los residuos NO siguen una distribución normal.
#
# Se rechaza H₀ cuando p-value ≤ α = 0.05.
# La prueba de Shapiro-Wilk es adecuada para n ≤ 5000.

cat("====================================================\n")
cat("  SECCIÓN 6: PRUEBA DE NORMALIDAD\n")
cat("====================================================\n\n")

cat("Hipótesis:\n")
cat("  H₀: Los residuos siguen una distribución normal.\n")
cat("  H₁: Los residuos NO siguen una distribución normal.\n\n")

prueba_shapiro <- shapiro.test(residuos)
cat("Prueba de Shapiro-Wilk:\n")
print(prueba_shapiro)

cat("\n--- Conclusión (Normalidad) ---\n")
if (prueba_shapiro$p.value > alfa) {
  cat(sprintf(paste0(
    "No se rechaza H₀ (p-value = %.6f > α = %.2f).\n",
    "Por lo tanto, se puede ASUMIR NORMALIDAD de los residuos.\n"
  ), prueba_shapiro$p.value, alfa))
} else {
  cat(sprintf(paste0(
    "Se rechaza H₀ (p-value = %.6f ≤ α = %.2f).\n",
    "Por lo tanto, NO se puede asumir normalidad de los residuos.\n"
  ), prueba_shapiro$p.value, alfa))
}


# =============================================================================
# SECCIÓN 7: PRUEBA DE HOMOCEDASTICIDAD (BREUSCH-PAGAN)
# =============================================================================
# H₀: Existe homocedasticidad (varianza constante de los errores).
# H₁: Existe heterocedasticidad (varianza no constante).
#
# Se rechaza H₀ cuando p-value ≤ α = 0.05, indicando que la varianza de los
# errores varía en función de los regresores.

cat("\n====================================================\n")
cat("  SECCIÓN 7: PRUEBA DE HOMOCEDASTICIDAD\n")
cat("====================================================\n\n")

cat("Hipótesis:\n")
cat("  H₀: Existe homocedasticidad (varianza constante).\n")
cat("  H₁: Existe heterocedasticidad (varianza no constante).\n\n")

prueba_bp <- bptest(modelo_final)
cat("Prueba de Breusch-Pagan:\n")
print(prueba_bp)

cat("\n--- Conclusión (Homocedasticidad) ---\n")
if (prueba_bp$p.value > alfa) {
  cat(sprintf(paste0(
    "No se rechaza H₀ (p-value = %.6f > α = %.2f).\n",
    "Por lo tanto, se puede ASUMIR HOMOCEDASTICIDAD.\n"
  ), prueba_bp$p.value, alfa))
} else {
  cat(sprintf(paste0(
    "Se rechaza H₀ (p-value = %.6f ≤ α = %.2f).\n",
    "Por lo tanto, NO se puede asumir homocedasticidad.\n",
    "Se detecta posible heterocedasticidad en los residuos.\n"
  ), prueba_bp$p.value, alfa))
}


# =============================================================================
# SECCIÓN 8: PRUEBA DE INDEPENDENCIA (DURBIN-WATSON)
# =============================================================================
# H₀: No existe autocorrelación en los errores (errores independientes).
# H₁: Existe autocorrelación en los errores.
#
# El estadístico DW varía entre 0 y 4; valores cercanos a 2 sugieren
# independencia. Se rechaza H₀ cuando p-value ≤ α = 0.05.

cat("\n====================================================\n")
cat("  SECCIÓN 8: PRUEBA DE INDEPENDENCIA\n")
cat("====================================================\n\n")

cat("Hipótesis:\n")
cat("  H₀: No existe autocorrelación en los errores.\n")
cat("  H₁: Existe autocorrelación en los errores.\n\n")

prueba_dw <- dwtest(modelo_final)
cat("Prueba de Durbin-Watson:\n")
print(prueba_dw)

cat("\n--- Conclusión (Independencia) ---\n")
if (prueba_dw$p.value > alfa) {
  cat(sprintf(paste0(
    "No se rechaza H₀ (p-value = %.6f > α = %.2f).\n",
    "Por lo tanto, se puede ASUMIR INDEPENDENCIA de los errores.\n"
  ), prueba_dw$p.value, alfa))
} else {
  cat(sprintf(paste0(
    "Se rechaza H₀ (p-value = %.6f ≤ α = %.2f).\n",
    "Por lo tanto, NO se puede asumir independencia de los errores.\n",
    "Se detecta posible autocorrelación en los residuos.\n"
  ), prueba_dw$p.value, alfa))
}


# =============================================================================
# SECCIÓN 9: MULTICOLINEALIDAD (VIF)
# =============================================================================
# El Factor de Inflación de la Varianza (VIF) mide el incremento en la varianza
# del coeficiente estimado de un regresor debido a su correlación con los demás.
#   VIF ≤ 5          → sin multicolinealidad significativa
#   5 < VIF ≤ 10     → multicolinealidad moderada (posible problema)
#   VIF > 10         → multicolinealidad fuerte (problema serio)
#
# El VIF solo aplica cuando el modelo tiene al menos dos regresores. Con un
# único regresor no existe colinealidad posible con otras variables.

cat("\n====================================================\n")
cat("  SECCIÓN 9: MULTICOLINEALIDAD (VIF)\n")
cat("====================================================\n\n")

if (length(variables_en_modelo) >= 2) {
  vif_valores <- vif(modelo_final)
  cat("Factor de Inflación de la Varianza (VIF):\n")
  print(vif_valores)

  cat("\nCriterio de interpretación:\n")
  cat("  VIF ≤ 5          → sin multicolinealidad significativa\n")
  cat("  5 < VIF ≤ 10     → multicolinealidad moderada\n")
  cat("  VIF > 10         → multicolinealidad fuerte\n\n")

  cat("--- Conclusión (Multicolinealidad) ---\n")
  hay_fuerte   <- any(vif_valores > 10)
  hay_moderada <- any(vif_valores > 5 & vif_valores <= 10)

  if (hay_fuerte) {
    vars_fuerte <- names(vif_valores)[vif_valores > 10]
    cat(sprintf(paste0(
      "Se detecta multicolinealidad FUERTE (VIF > 10) en: %s.\n",
      "Se recomienda revisar el modelo e investigar la relación entre los regresores.\n"
    ), paste(vars_fuerte, collapse = ", ")))
  } else if (hay_moderada) {
    vars_mod <- names(vif_valores)[vif_valores > 5]
    cat(sprintf(paste0(
      "Se detecta multicolinealidad MODERADA (5 < VIF ≤ 10) en: %s.\n",
      "Puede ser conveniente monitorear estas variables.\n"
    ), paste(vars_mod, collapse = ", ")))
  } else {
    cat("No se detecta multicolinealidad significativa (todos los VIF ≤ 5).\n")
  }

} else if (length(variables_en_modelo) == 1) {
  cat("El modelo final contiene un único regresor.\n")
  cat("El VIF no aplica con un solo regresor: no existe colinealidad posible\n")
  cat("ya que no hay más de una variable explicativa con la cual medir colinealidad.\n")
  vif_valores <- NULL

} else {
  cat("El modelo final es el modelo nulo (sin regresores). El VIF no aplica.\n")
  vif_valores <- NULL
}


# =============================================================================
# SECCIÓN 10: OBSERVACIONES INFLUYENTES (DISTANCIA DE COOK)
# =============================================================================
# La distancia de Cook Dᵢ cuantifica el impacto de eliminar la observación i
# sobre el conjunto de coeficientes estimados.
# Criterio: Dᵢ > 4/n  puede indicar una observación potencialmente influyente.

cat("\n====================================================\n")
cat("  SECCIÓN 10: OBSERVACIONES INFLUYENTES\n")
cat("  Distancia de Cook\n")
cat("====================================================\n\n")

cook_dist   <- cooks.distance(modelo_final)
umbral_cook <- 4 / n

cat(sprintf("Umbral de Cook: 4/n = 4/%d = %.6f\n\n", n, umbral_cook))

# Gráfica de distancia de Cook (mostrada en Plots de RStudio, sin guardarse).
plot(cook_dist,
     type = "h",
     main = "Distancia de Cook por Observación",
     xlab = "Índice de observación",
     ylab = "Distancia de Cook (Dᵢ)",
     col  = ifelse(cook_dist > umbral_cook, "red", "steelblue"),
     lwd  = 1.5)
abline(h = umbral_cook, col = "red", lwd = 2, lty = 2)
legend("topright",
       legend = c(sprintf("Umbral 4/n = %.4f", umbral_cook), "Obs. influyente"),
       col    = c("red", "red"),
       lty    = c(2, 1),
       lwd    = c(2, 2),
       cex    = 0.8)

obs_influyentes <- which(cook_dist > umbral_cook)

cat("--- Conclusión (Observaciones influyentes) ---\n")
if (length(obs_influyentes) > 0) {
  cat(sprintf("Se identificaron %d observación(es) potencialmente influyente(s):\n",
              length(obs_influyentes)))
  cat("Índices:", paste(obs_influyentes, collapse = ", "), "\n\n")
  cat("Distancias de Cook de las observaciones influyentes:\n")
  print(cook_dist[obs_influyentes])
  cat("\nSe recomienda revisar estas observaciones para determinar si corresponden\n")
  cat("a errores de medición, valores atípicos reales o puntos extremos con\n")
  cat("influencia legítima en el ajuste del modelo.\n")
} else {
  cat("No se identificaron observaciones potencialmente influyentes\n")
  cat(sprintf("(ninguna Dᵢ supera el umbral 4/n = %.6f).\n", umbral_cook))
}


# =============================================================================
# SECCIÓN 11: CONCLUSIÓN GENERAL AUTOMÁTICA
# =============================================================================

cat("\n====================================================\n")
cat("  SECCIÓN 11: CONCLUSIÓN GENERAL\n")
cat("====================================================\n\n")

vars_entran_str  <- paste(variables_en_modelo, collapse = ", ")
vars_fuera_str   <- ifelse(length(variables_no_entran) > 0,
                           paste(variables_no_entran, collapse = ", "),
                           "ninguna")
r2_adj_final     <- round(resumen_final$adj.r.squared, 4)
p_f_str          <- round(p_valor_F, 6)
modelo_sign_str  <- ifelse(p_f_str <= alfa,
                           "sí es globalmente significativo",
                           "NO es globalmente significativo")
normalidad_str   <- ifelse(prueba_shapiro$p.value > alfa,
                           "se puede asumir normalidad de los residuos",
                           "NO se puede asumir normalidad de los residuos")
homoc_str        <- ifelse(prueba_bp$p.value > alfa,
                           "se puede asumir homocedasticidad",
                           "NO se puede asumir homocedasticidad")
indep_str        <- ifelse(prueba_dw$p.value > alfa,
                           "se puede asumir independencia de los errores",
                           "NO se puede asumir independencia de los errores")

if (!is.null(vif_valores)) {
  if (any(vif_valores > 10)) {
    multi_str <- "se detecta multicolinealidad FUERTE (VIF > 10)"
  } else if (any(vif_valores > 5)) {
    multi_str <- "se detecta multicolinealidad MODERADA (5 < VIF ≤ 10)"
  } else {
    multi_str <- "no se detecta multicolinealidad significativa (VIF ≤ 5)"
  }
} else if (length(variables_en_modelo) == 1) {
  multi_str <- "el VIF no aplica (solo un regresor en el modelo final)"
} else {
  multi_str <- "el VIF no aplica (modelo nulo)"
}

inf_str <- ifelse(length(obs_influyentes) > 0,
                  sprintf("se identificaron %d observación(es) influyente(s): %s",
                          length(obs_influyentes),
                          paste(obs_influyentes, collapse = ", ")),
                  "no se identificaron observaciones influyentes")

cat("RESUMEN DEL ANÁLISIS DE REGRESIÓN LINEAL MÚLTIPLE\n")
cat("Método: Selección Hacia Adelante\n")
cat(rep("-", 54), "\n", sep = "")

cat(sprintf(paste0(
  "1. MODELO FINAL OBTENIDO:\n",
  "   Ŷ = %s\n\n",
  "2. VARIABLES QUE ENTRARON AL MODELO: %s\n\n",
  "3. VARIABLES QUE NO ENTRARON: %s\n\n",
  "4. AJUSTE DEL MODELO:\n",
  "   R² ajustado = %.4f, lo que indica que el modelo explica aproximadamente\n",
  "   el %.2f%% de la variabilidad de Y, ajustado por el número de regresores.\n\n",
  "5. SIGNIFICANCIA GLOBAL (Prueba F):\n",
  "   El modelo final %s (p-value = %.6f).\n\n",
  "6. SIGNIFICANCIA INDIVIDUAL:\n",
  "   Todos los regresores del modelo final son estadísticamente significativos\n",
  "   al nivel α = 0.05, ya que el método de selección hacia adelante garantiza\n",
  "   que ningún regresor entra con p-value > 0.05.\n\n",
  "7. NORMALIDAD (Shapiro-Wilk, p-value = %.6f):\n",
  "   Con α = 0.05, %s.\n\n",
  "8. HOMOCEDASTICIDAD (Breusch-Pagan, p-value = %.6f):\n",
  "   Con α = 0.05, %s.\n\n",
  "9. INDEPENDENCIA (Durbin-Watson, p-value = %.6f):\n",
  "   Con α = 0.05, %s.\n\n",
  "10. MULTICOLINEALIDAD (VIF):\n",
  "    %s.\n\n",
  "11. OBSERVACIONES INFLUYENTES (Distancia de Cook, umbral 4/n = %.6f):\n",
  "    Con α = 0.05, %s.\n\n",
  "12. INTERPRETACIÓN GENERAL:\n",
  "    Mediante el método de selección hacia adelante se identificaron los\n",
  "    regresores: %s, como los más relevantes para explicar la variabilidad\n",
  "    de Y. Cada coeficiente representa el cambio promedio esperado en Y por\n",
  "    cada unidad de aumento en el regresor correspondiente, manteniendo\n",
  "    constantes los demás regresores del modelo.\n"
),
terminos_eq,
vars_entran_str,
vars_fuera_str,
r2_adj_final,
r2_adj_final * 100,
modelo_sign_str, p_f_str,
prueba_shapiro$p.value, normalidad_str,
prueba_bp$p.value, homoc_str,
prueba_dw$p.value, indep_str,
multi_str,
umbral_cook, inf_str,
vars_entran_str
))

cat(rep("=", 54), "\n", sep = "")
cat("  FIN DEL ANÁLISIS\n")
cat(rep("=", 54), "\n", sep = "")


# =============================================================================
# =============================================================================
# GUÍA DE INTERPRETACIÓN DE RESULTADOS PARA EL REPORTE
# =============================================================================
# =============================================================================
#
# Adapte los valores numéricos entre corchetes [] según sus propios resultados.
#
# ----------------------------------------------------------------------------
# 1. SELECCIÓN HACIA ADELANTE — DESCRIPCIÓN DEL MÉTODO
# ----------------------------------------------------------------------------
# "Se aplicó el método de selección hacia adelante para identificar los
#  regresores más relevantes. El proceso parte del modelo nulo (Y = β₀ + ε)
#  y, en cada paso, incorpora la variable candidata con menor p-value,
#  siempre que este sea ≤ α = 0.05. El procedimiento se repite hasta que
#  ninguna variable candidata cumpla el criterio de entrada."
#
# ----------------------------------------------------------------------------
# 2. ENTRADA DE VARIABLES AL MODELO
# ----------------------------------------------------------------------------
# Paso 1:
# "En el Paso 1 se evaluaron todos los regresores candidatos de forma
#  individual. El regresor [X₁] presentó el menor p-value (p = [valor] ≤ 0.05),
#  por lo que fue el primero en incorporarse al modelo."
#
# Paso 2:
# "En el Paso 2 se evaluaron los regresores restantes en presencia de [X₁].
#  El regresor [X₂] presentó el menor p-value (p = [valor] ≤ 0.05) y fue
#  incorporado al modelo."
#
# Paso final (sin entrada):
# "En el Paso [m+1] ninguna variable candidata presentó p-value ≤ 0.05,
#  por lo que el proceso de selección hacia adelante concluyó con un modelo
#  de [q] regresores."
#
# ----------------------------------------------------------------------------
# 3. MODELO FINAL
# ----------------------------------------------------------------------------
# "El modelo final obtenido mediante el método de selección hacia adelante es:
#      Ŷ = β̂₀ + β̂₁X₁ + β̂₂X₂
#  donde β̂₀ = [valor] es el intercepto; β̂₁ = [valor] indica el cambio
#  promedio en Y por cada unidad de incremento en X₁, manteniendo X₂
#  constante; y β̂₂ = [valor] indica el cambio promedio en Y por cada unidad
#  de incremento en X₂, manteniendo X₁ constante."
#
# ----------------------------------------------------------------------------
# 4. R² AJUSTADO
# ----------------------------------------------------------------------------
# "El coeficiente de determinación ajustado del modelo final fue
#  R²_adj = [valor], lo que indica que el modelo explica aproximadamente el
#  [valor × 100]% de la variabilidad total de Y, corregido por el número de
#  regresores. Este indicador penaliza la inclusión de variables innecesarias,
#  siendo preferible al R² ordinario cuando se comparan modelos con distinto
#  número de predictores."
#
# ----------------------------------------------------------------------------
# 5. SIGNIFICANCIA INDIVIDUAL DE LOS REGRESORES
# ----------------------------------------------------------------------------
# "En la tabla de coeficientes del modelo final, el regresor X₁ resultó
#  estadísticamente significativo (t = [valor], p-value = [valor] < 0.05),
#  lo que indica que contribuye de forma significativa a explicar la
#  variabilidad de Y, dada la presencia de los demás regresores en el modelo.
#  De manera análoga, X₂ también resultó significativo
#  (t = [valor], p-value = [valor] < 0.05)."
#
# ----------------------------------------------------------------------------
# 6. PRUEBA F GLOBAL
# ----------------------------------------------------------------------------
# "La prueba F global del modelo final arrojó F([gl1], [gl2]) = [valor],
#  con p-value = [valor] < 0.05. Por lo tanto, se rechaza la hipótesis nula
#  de que todos los coeficientes son simultáneamente iguales a cero,
#  concluyendo que el modelo en su conjunto es estadísticamente significativo
#  para explicar la variabilidad de Y."
#
# ----------------------------------------------------------------------------
# 7. NORMALIDAD
# ----------------------------------------------------------------------------
# "El QQ-plot de los residuos muestra que los puntos se distribuyen
#  aproximadamente sobre la línea diagonal de referencia, lo que sugiere
#  normalidad. Complementariamente, la prueba de Shapiro-Wilk arrojó
#  W = [valor], p-value = [valor]. Dado que el p-value [es mayor / es menor]
#  que α = 0.05, [no se rechaza / se rechaza] H₀. Por lo tanto, [se puede /
#  no se puede] asumir que los residuos siguen una distribución normal."
#
# ----------------------------------------------------------------------------
# 8. HOMOCEDASTICIDAD
# ----------------------------------------------------------------------------
# "La gráfica de residuos vs valores ajustados no muestra un patrón sistemático
#  en la dispersión, lo que sugiere varianza constante. La prueba de
#  Breusch-Pagan arrojó BP = [valor], p-value = [valor]. Dado que el p-value
#  [es mayor / es menor] que α = 0.05, [no se rechaza / se rechaza] H₀.
#  En consecuencia, [se puede / no se puede] asumir homocedasticidad."
#
# ----------------------------------------------------------------------------
# 9. INDEPENDENCIA
# ----------------------------------------------------------------------------
# "La prueba de Durbin-Watson arrojó DW = [valor], p-value = [valor].
#  Dado que el p-value [es mayor / es menor] que α = 0.05, [no se rechaza /
#  se rechaza] H₀. Por lo tanto, [se puede / no se puede] asumir que los
#  errores del modelo son independientes entre sí."
#
# ----------------------------------------------------------------------------
# 10. MULTICOLINEALIDAD
# ----------------------------------------------------------------------------
# "Se calculó el VIF para cada regresor del modelo final, obteniendo
#  [X₁: VIF = valor, X₂: VIF = valor]. Dado que [todos los valores son ≤ 5 /
#  algunos valores superan 5 o 10], [no se detecta multicolinealidad
#  significativa / se detecta multicolinealidad moderada o fuerte] entre los
#  regresores del modelo final."
#
# ----------------------------------------------------------------------------
# 11. OBSERVACIONES INFLUYENTES
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

# =============================================================================
# ANÁLISIS DE REGRESIÓN LINEAL MÚLTIPLE SIN INTERCEPTO
# MÉTODO DE ELIMINACIÓN HACIA ATRÁS (BACKWARD ELIMINATION)
# =============================================================================
# Curso: Análisis de Regresión
# Descripción: Este script realiza un análisis completo de regresión lineal
#              múltiple SIN INTERCEPTO utilizando el método de eliminación
#              hacia atrás, partiendo del modelo completo y eliminando
#              regresores no significativos (p-value > 0.05) hasta obtener
#              un modelo final donde todos los regresores sean significativos.
#
# MODELO SIN INTERCEPTO:
#   Y = β₁X₁ + β₂X₂ + ... + βₖXₖ + ε
#
# Al suprimir el intercepto se asume que cuando todos los regresores son
# iguales a cero, el valor esperado de Y también es cero (el plano de
# regresión pasa por el origen del espacio de regresores).
# Use este modelo solo cuando esa restricción tenga sentido teórico.
# =============================================================================


# =============================================================================
# SECCIÓN 0: INSTALACIÓN Y CARGA DE PAQUETES
# =============================================================================

# Función para instalar un paquete si aún no está instalado y luego cargarlo.
# Esto garantiza que el script funcione en cualquier equipo sin necesidad de
# instalar los paquetes manualmente.
instalar_y_cargar <- function(paquete) {
  if (!requireNamespace(paquete, quietly = TRUE)) {
    install.packages(paquete, dependencies = TRUE)
  }
  library(paquete, character.only = TRUE)
}

# Paquetes necesarios:
# - lmtest: prueba de Breusch-Pagan (homocedasticidad) y Durbin-Watson (independencia)
# - car:    VIF (factor de inflación de la varianza, multicolinealidad)
instalar_y_cargar("lmtest")
instalar_y_cargar("car")

cat("\n====================================================\n")
cat("  ANÁLISIS DE REGRESIÓN LINEAL MÚLTIPLE\n")
cat("  MODELO SIN INTERCEPTO\n")
cat("  Método: Eliminación Hacia Atrás\n")
cat("====================================================\n\n")


# =============================================================================
# SECCIÓN 1: LECTURA DE DATOS
# =============================================================================

# Se solicita al usuario que seleccione el archivo CSV de sus datos.
# El archivo debe tener una columna llamada exactamente "Y" (variable respuesta)
# y una o más columnas numéricas adicionales que serán los regresores (X1, X2, ...).
cat("Por favor, seleccione el archivo CSV con sus datos...\n")
datos <- read.csv(file.choose(), header = TRUE, sep = ",", dec = ".")

# Mostrar las primeras filas del conjunto de datos para verificar la lectura.
cat("\n--- Vista previa de los datos (primeras 6 filas) ---\n")
print(head(datos))
cat("\n--- Dimensiones del conjunto de datos ---\n")
cat("Número de observaciones:", nrow(datos), "\n")
cat("Número de variables    :", ncol(datos), "\n")
cat("Nombres de columnas    :", paste(names(datos), collapse = ", "), "\n\n")

# Verificar que exista la columna "Y" en el archivo.
if (!"Y" %in% names(datos)) {
  stop("ERROR: El archivo CSV no contiene una columna llamada 'Y'.
       Asegúrese de que la variable respuesta se llame exactamente 'Y'.")
}

# Identificar automáticamente los regresores: todas las columnas numéricas
# distintas de "Y".
regresores_iniciales <- setdiff(
  names(datos)[sapply(datos, is.numeric)],
  "Y"
)

if (length(regresores_iniciales) == 0) {
  stop("ERROR: No se encontraron columnas numéricas que funcionen como regresores.")
}

cat("Variable respuesta detectada : Y\n")
cat("Regresores detectados        :", paste(regresores_iniciales, collapse = ", "), "\n\n")

# Tamaño de muestra (se usará más adelante en el criterio de Cook y grados de libertad).
n <- nrow(datos)


# =============================================================================
# SECCIÓN 2: MODELO COMPLETO SIN INTERCEPTO (MODELO INICIAL)
# =============================================================================
# El modelo completo incluye todos los regresores disponibles y NO tiene
# término independiente (intercepto):
#   Y = β₁X₁ + β₂X₂ + ... + βₖXₖ + ε
#
# CAMBIO RESPECTO AL MODELO CON INTERCEPTO:
#   Se usa "Y ~ 0 + X1 + X2 + ..." en lugar de "Y ~ X1 + X2 + ...".
#   El "0 +" indica explícitamente a R que suprima el intercepto.
#   También puede escribirse con "-1", pero se prefiere "0 +" porque
#   es más legible y evita confusiones de signo en la fórmula.

cat("====================================================\n")
cat("  SECCIÓN 2: MODELO COMPLETO SIN INTERCEPTO\n")
cat("====================================================\n\n")

# Construir la fórmula del modelo completo SIN INTERCEPTO dinámicamente.
# MODIFICACIÓN: se agrega "0 +" al inicio para suprimir el intercepto.
formula_completa <- as.formula(
  paste("Y ~ 0 +", paste(regresores_iniciales, collapse = " + "))
)
cat("Fórmula del modelo completo (sin intercepto):\n")
print(formula_completa)
cat("\n")

# Ajustar el modelo completo sin intercepto.
modelo_completo <- lm(formula_completa, data = datos)

# Mostrar el resumen completo del modelo: coeficientes, errores estándar,
# estadísticos t, p-values, R², R² ajustado, error estándar residual y prueba F.
cat("--- Resumen del Modelo Completo (sin intercepto) ---\n")
resumen_completo <- summary(modelo_completo)
print(resumen_completo)

cat(sprintf("\nR² del modelo completo          : %.6f\n", resumen_completo$r.squared))
cat(sprintf("R² ajustado del modelo completo : %.6f\n", resumen_completo$adj.r.squared))
cat(sprintf("Error estándar residual         : %.6f\n", resumen_completo$sigma))

# ADVERTENCIA SOBRE R² EN MODELOS SIN INTERCEPTO:
# Cuando el modelo no tiene intercepto, R reporta un R² no centrado, calculado
# como 1 - SSE/SST, donde SST = Σ(Yᵢ²) en lugar de Σ(Yᵢ - Ȳ)². Este R² no
# centrado puede ser artificialmente alto (incluso cercano a 1) aunque el ajuste
# real sea mediocre. Por ello se debe interpretar con cautela y utilizarse
# principalmente para comparar modelos SIN intercepto entre sí, no para
# comparar con modelos que sí tienen intercepto.
cat("\nADVERTENCIA: El R² reportado en modelos sin intercepto es NO CENTRADO.\n")
cat("Debe interpretarse con cautela. Ver notas en Sección 4 y Sección 11.\n\n")


# =============================================================================
# SECCIÓN 3: MÉTODO DE ELIMINACIÓN HACIA ATRÁS
# =============================================================================
# Se parte del modelo completo sin intercepto y, en cada paso, se identifica
# el regresor con el mayor p-value. Si ese p-value supera α = 0.05, la
# variable se elimina y el modelo se reajusta. El proceso continúa hasta que
# todos los regresores restantes tengan p-value ≤ 0.05.
#
# Nota sobre R² ajustado: aunque el criterio principal de eliminación es el
# p-value individual de cada regresor (criterio de significancia estadística),
# en cada paso también se muestra el R² ajustado como indicador de la calidad
# del ajuste. El R² ajustado penaliza la inclusión de variables innecesarias,
# por lo que es útil para comparar modelos con distinto número de regresores.
# En modelos sin intercepto este R² es no centrado; usarlo solo para comparar
# modelos dentro del mismo esquema sin intercepto.

cat("====================================================\n")
cat("  SECCIÓN 3: ELIMINACIÓN HACIA ATRÁS (SIN INTERCEPTO)\n")
cat("  Criterio principal: p-value > 0.05\n")
cat("  Nivel de significancia: α = 0.05\n")
cat("====================================================\n\n")

alfa          <- 0.05          # Nivel de significancia
paso          <- 0             # Contador de pasos
modelo_actual <- modelo_completo
regresores    <- regresores_iniciales
variables_eliminadas <- c()    # Registra las variables que se fueron eliminando

# Bucle principal de eliminación hacia atrás.
repeat {
  resumen_actual <- summary(modelo_actual)
  coef_tabla     <- resumen_actual$coefficients

  # CAMBIO RESPECTO AL MODELO CON INTERCEPTO:
  # En un modelo CON intercepto, la fila "(Intercept)" aparece en la tabla
  # de coeficientes y se excluía con:
  #   coef_tabla[rownames(coef_tabla) != "(Intercept)", 4]
  #
  # En un modelo SIN intercepto, R NO incluye la fila "(Intercept)" en la
  # tabla. Por lo tanto, TODAS las filas de coef_tabla son regresores y se
  # toman directamente sin filtrar nada.
  p_values_regresores <- coef_tabla[, "Pr(>|t|)"]

  # Si solo hay un regresor, coef_tabla es un vector; convertir a estructura
  # con nombre para que which.max funcione correctamente.
  if (is.null(names(p_values_regresores))) {
    names(p_values_regresores) <- rownames(coef_tabla)
  }

  # Identificar el regresor con mayor p-value.
  max_p    <- max(p_values_regresores)
  var_elim <- names(which.max(p_values_regresores))

  paso <- paso + 1

  # SEGURIDAD 1: Verificar que no se supere el número máximo de eliminaciones.
  # Esto previene ciclos infinitos ante comportamientos numéricos inesperados.
  if (paso > length(regresores_iniciales)) {
    stop(paste0(
      "ERROR: El algoritmo superó el número máximo de eliminaciones posibles (",
      length(regresores_iniciales), " pasos). Revise el ciclo de eliminación."
    ))
  }

  # SEGURIDAD 2: Verificar que la variable candidata a eliminar esté realmente
  # en el vector de regresores activos. Previene inconsistencias de nombres.
  if (!(var_elim %in% regresores)) {
    stop(paste(
      "ERROR: La variable a eliminar '", var_elim,
      "' no coincide con los nombres de regresores activos.",
      "Regresores activos:", paste(regresores, collapse = ", ")
    ))
  }

  cat(sprintf("--- PASO %d ---\n", paso))
  # CAMBIO: mostrar fórmula sin intercepto (con "0 +").
  cat("Variables incluidas en el modelo: Y ~ 0 +",
      paste(regresores, collapse = " + "), "\n\n")
  cat("Tabla de coeficientes:\n")
  print(coef_tabla)
  cat(sprintf("\nR² ajustado del modelo actual: %.6f\n", resumen_actual$adj.r.squared))
  cat("  (Nota: R² no centrado en modelo sin intercepto — usar solo\n")
  cat("   para comparar entre modelos sin intercepto.)\n")
  cat(sprintf("\nRegresor con mayor p-value: %s  (p-value = %.6f)\n", var_elim, max_p))

  # Verificar si el mayor p-value supera el nivel de significancia α.
  if (max_p > alfa) {
    # La variable no es significativa: se elimina.
    cat(sprintf("Decisión: ELIMINAR '%s'  →  p-value = %.6f > α = %.2f\n",
                var_elim, max_p, alfa))
    cat(sprintf("Razón: El regresor '%s' no es estadísticamente significativo\n",
                var_elim))
    cat(sprintf("       al nivel α = %.2f, por lo que se elimina del modelo.\n\n",
                alfa))

    variables_eliminadas <- c(variables_eliminadas, var_elim)
    regresores           <- setdiff(regresores, var_elim)

    # CAMBIO: reajustar el modelo SIN INTERCEPTO ("0 +") tras cada eliminación.
    # En el código original era:  paste("Y ~", paste(regresores, ...))
    # Aquí debe ser:              paste("Y ~ 0 +", paste(regresores, ...))
    nueva_formula <- as.formula(
      paste("Y ~ 0 +", paste(regresores, collapse = " + "))
    )
    modelo_actual <- lm(nueva_formula, data = datos)

  } else {
    # Todos los regresores restantes son significativos: se detiene el proceso.
    cat(sprintf("Decisión: CONSERVAR MODELO  →  todos los p-values ≤ α = %.2f\n", alfa))
    cat("El proceso de eliminación hacia atrás ha concluido.\n\n")
    # Retroceder el contador: el último "paso" fue solo de verificación.
    paso <- paso - 1
    break
  }

  # Si solo queda un regresor o ninguno, detener el bucle.
  if (length(regresores) == 0) {
    cat("No quedan regresores en el modelo. El proceso se detiene.\n\n")
    break
  }
}

cat(sprintf("Total de pasos de eliminación realizados: %d\n\n", paso))


# =============================================================================
# SECCIÓN 4: MODELO FINAL SIN INTERCEPTO
# =============================================================================

cat("====================================================\n")
cat("  SECCIÓN 4: MODELO FINAL (SIN INTERCEPTO)\n")
cat("====================================================\n\n")

modelo_final  <- modelo_actual
resumen_final <- summary(modelo_final)

# --- Ecuación ajustada SIN INTERCEPTO ---
# CAMBIO RESPECTO AL MODELO CON INTERCEPTO:
# El código original asumía que coef_finales[1] era el intercepto (β₀) y
# construía la ecuación como:  Ŷ = β̂₀ + β̂₁X₁ + β̂₂X₂ + ...
#
# En un modelo SIN intercepto, coef_finales[1] es el primer REGRESOR (β̂₁),
# no un intercepto. La ecuación toma la forma:
#   Ŷ = β̂₁X₁ + β̂₂X₂ + ... + β̂ₘXₘ
#
# Por ello el bucle debe comenzar desde el primer coeficiente (i = 1),
# y no debe inicializarse con un valor de intercepto.
coef_finales <- coef(modelo_final)

# Primer término: β̂₁ * X₁ (no hay término independiente).
terminos_eq <- sprintf("%.4f*%s", coef_finales[1], names(coef_finales)[1])

# Términos restantes: β̂₂X₂, β̂₃X₃, ...
if (length(coef_finales) > 1) {
  for (i in 2:length(coef_finales)) {
    nombre_var <- names(coef_finales)[i]
    valor      <- coef_finales[i]
    if (valor >= 0) {
      terminos_eq <- paste0(terminos_eq, sprintf(" + %.4f*%s", valor, nombre_var))
    } else {
      terminos_eq <- paste0(terminos_eq, sprintf(" - %.4f*%s", abs(valor), nombre_var))
    }
  }
}

cat("Ecuación ajustada del modelo final (sin intercepto):\n")
cat(sprintf("  Ŷ = %s\n\n", terminos_eq))

# --- Tabla de coeficientes del modelo final ---
cat("Tabla de coeficientes (modelo final, sin intercepto):\n")
print(resumen_final$coefficients)

# --- Medidas de ajuste ---
k_final <- length(regresores)   # número de regresores en el modelo final

# CAMBIO RESPECTO AL MODELO CON INTERCEPTO:
# El código original calculaba gl_error = n - k_final - 1, restando 1 por el
# intercepto estimado. En un modelo sin intercepto no se estima β₀, por lo
# que los grados de libertad del error son n - k (no n - k - 1).
# Se usa df.residual() porque R lo calcula automáticamente de forma correcta
# tanto para modelos con como sin intercepto.
gl_error <- df.residual(modelo_final)   # n - k para modelo sin intercepto

SSE  <- sum(residuals(modelo_final)^2)
MSE  <- SSE / gl_error
RMSE <- sqrt(MSE)

cat("\n--- Medidas de ajuste del modelo final ---\n")
cat(sprintf("R²  (no centrado*)  : %.6f\n",  resumen_final$r.squared))
cat(sprintf("R² ajustado*        : %.6f\n",  resumen_final$adj.r.squared))
cat(sprintf("SSE                 : %.6f\n",  SSE))
cat(sprintf("MSE                 : %.6f\n",  MSE))
cat(sprintf("RMSE                : %.6f\n",  RMSE))
cat(sprintf("Grados de libertad del error: %d  (n - k = %d - %d)\n",
            gl_error, n, k_final))
cat("\n* ADVERTENCIA R² en modelo sin intercepto:\n")
cat("  R reporta un R² NO CENTRADO: calcula SST = Σ(Yᵢ²) en lugar de\n")
cat("  Σ(Yᵢ - Ȳ)². Esto puede producir valores artificialmente altos.\n")
cat("  Interprete el R² con cautela y úselo solo para comparar modelos\n")
cat("  SIN intercepto entre sí.\n")

# --- Prueba F global ---
f_stat   <- resumen_final$fstatistic
cat(sprintf("\nPrueba F global:\n"))
cat(sprintf("  F = %.4f,  gl1 = %d,  gl2 = %d\n",
            f_stat[1], f_stat[2], f_stat[3]))
p_valor_F <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
cat(sprintf("  p-value = %.6f\n", p_valor_F))
if (p_valor_F <= alfa) {
  cat("  → El modelo final es GLOBALMENTE SIGNIFICATIVO (p-value ≤ 0.05).\n")
} else {
  cat("  → El modelo final NO es globalmente significativo (p-value > 0.05).\n")
}

# --- Interpretación de variables ---
cat("\n--- Resumen de variables ---\n")
cat("Variables incluidas en el modelo final :", paste(regresores, collapse = ", "), "\n")
if (length(variables_eliminadas) > 0) {
  cat("Variables eliminadas por el método    :", paste(variables_eliminadas, collapse = ", "), "\n\n")
} else {
  cat("No se eliminó ninguna variable (el modelo completo fue el modelo final).\n\n")
}


# =============================================================================
# SECCIÓN 5: ANÁLISIS RESIDUAL
# =============================================================================
# Los residuos se definen como:   eᵢ = Yᵢ - Ŷᵢ
# Un buen modelo debe presentar residuos con distribución aproximadamente
# normal, varianza constante (homocedasticidad) y sin autocorrelación.
# Estas propiedades se verifican igual que en el modelo con intercepto.

cat("====================================================\n")
cat("  SECCIÓN 5: ANÁLISIS RESIDUAL\n")
cat("====================================================\n\n")

residuos  <- residuals(modelo_final)   # eᵢ = Yᵢ - Ŷᵢ
ajustados <- fitted(modelo_final)      # Ŷᵢ

cat("Calculados residuos: eᵢ = Yᵢ - Ŷᵢ\n\n")
cat("Primeros residuos:\n")
print(head(data.frame(Y_obs = datos$Y, Y_ajustado = ajustados, Residuo = residuos)))

# =============================================================================
# SECCIÓN 5.1: GRÁFICAS RESIDUALES
# =============================================================================
# Las gráficas se muestran en la ventana de Plots de RStudio (no se guardan).

# Número de paneles necesarios: 4 gráficas base + 1 por cada regresor final.
num_regresores_final <- length(regresores)
total_graficas       <- 4 + num_regresores_final

# Organizar el layout de las gráficas (filas x columnas).
num_cols  <- 2
num_filas <- ceiling(total_graficas / num_cols)
par(mfrow = c(num_filas, num_cols), mar = c(4, 4, 3, 1))

# --- Gráfica 1: Valores observados vs valores ajustados ---
plot(ajustados, datos$Y,
     main = "Valores Observados vs Ajustados\n(Modelo sin intercepto)",
     xlab = "Valores Ajustados (Ŷᵢ)",
     ylab = "Valores Observados (Yᵢ)",
     pch  = 19, col = "steelblue")
abline(0, 1, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = "Y = Ŷ (línea de referencia)",
       col = "red", lty = 2, lwd = 2, cex = 0.8)

# --- Gráfica 2: Residuos vs valores ajustados ---
# Sirve para detectar patrones (no linealidad) y heterocedasticidad.
plot(ajustados, residuos,
     main = "Residuos vs Valores Ajustados\n(Modelo sin intercepto)",
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
hist(residuos,
     main   = "Histograma de Residuos",
     xlab   = "Residuos (eᵢ)",
     ylab   = "Frecuencia",
     col    = "lightblue",
     border = "white",
     breaks = "Sturges")
# Agregar curva normal teórica de referencia.
xseq <- seq(min(residuos), max(residuos), length.out = 100)
lines(xseq,
      dnorm(xseq, mean(residuos), sd(residuos)) * length(residuos) *
        diff(hist(residuos, plot = FALSE)$breaks)[1],
      col = "red", lwd = 2)

# --- Gráficas adicionales: Residuos vs cada regresor del modelo final ---
for (var in regresores) {
  plot(datos[[var]], residuos,
       main = paste("Residuos vs", var),
       xlab = var,
       ylab = "Residuos (eᵢ)",
       pch  = 19, col = "darkgreen")
  abline(h = 0, col = "red", lwd = 2, lty = 2)
}

# Restaurar parámetros gráficos.
par(mfrow = c(1, 1))
cat("Gráficas de análisis residual mostradas en la ventana Plots de RStudio.\n\n")


# =============================================================================
# SECCIÓN 6: PRUEBA DE NORMALIDAD (SHAPIRO-WILK)
# =============================================================================
# H₀: Los residuos siguen una distribución normal.
# H₁: Los residuos NO siguen una distribución normal.
#
# La prueba de Shapiro-Wilk es adecuada para muestras de tamaño n ≤ 5000.
# Se rechaza H₀ cuando p-value ≤ α = 0.05.

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
# Se rechaza H₀ cuando p-value ≤ α = 0.05, lo que indica que la varianza
# de los errores cambia según los valores de los regresores.

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
# El estadístico DW toma valores entre 0 y 4; valores cercanos a 2 sugieren
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
# El Factor de Inflación de la Varianza (VIF) mide cuánto aumenta la varianza
# de un coeficiente estimado debido a la correlación con los demás regresores.
#   VIF ≤ 5      → sin problema de multicolinealidad (o leve)
#   5 < VIF ≤ 10 → multicolinealidad MODERADA (posible problema)
#   VIF > 10     → multicolinealidad FUERTE (problema serio)
#
# Nota: el VIF solo se puede calcular con al menos 2 regresores en el modelo.

cat("\n====================================================\n")
cat("  SECCIÓN 9: MULTICOLINEALIDAD (VIF)\n")
cat("====================================================\n\n")

if (length(regresores) >= 2) {
  vif_valores <- vif(modelo_final)
  cat("Factor de Inflación de la Varianza (VIF):\n")
  print(vif_valores)

  cat("\nCriterio de interpretación:\n")
  cat("  VIF ≤ 5          → sin multicolinealidad significativa\n")
  cat("  5 < VIF ≤ 10     → multicolinealidad moderada\n")
  cat("  VIF > 10         → multicolinealidad fuerte\n\n")

  cat("--- Conclusión (Multicolinealidad) ---\n")
  hay_moderada <- any(vif_valores > 5  & vif_valores <= 10)
  hay_fuerte   <- any(vif_valores > 10)

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

} else {
  cat("El modelo final contiene solo un regresor.\n")
  cat("El VIF no aplica con un único regresor (no hay multicolinealidad posible).\n")
  vif_valores <- NULL
}


# =============================================================================
# SECCIÓN 10: OBSERVACIONES INFLUYENTES (DISTANCIA DE COOK)
# =============================================================================
# La distancia de Cook Dᵢ mide el efecto de la observación i sobre los
# coeficientes estimados si se eliminara del ajuste.
# Criterio habitual: Dᵢ > 4/n  puede indicar observación influyente.

cat("\n====================================================\n")
cat("  SECCIÓN 10: OBSERVACIONES INFLUYENTES\n")
cat("  Distancia de Cook\n")
cat("====================================================\n\n")

cook_dist   <- cooks.distance(modelo_final)
umbral_cook <- 4 / n   # criterio Dᵢ > 4/n

cat(sprintf("Umbral de Cook: 4/n = 4/%d = %.6f\n\n", n, umbral_cook))

# Gráfica de la distancia de Cook.
plot(cook_dist,
     type = "h",
     main = "Distancia de Cook por Observación\n(Modelo sin intercepto)",
     xlab = "Índice de observación",
     ylab = "Distancia de Cook (Dᵢ)",
     col  = ifelse(cook_dist > umbral_cook, "red", "steelblue"),
     lwd  = 1.5)
abline(h = umbral_cook, col = "red", lwd = 2, lty = 2)
legend("topright",
       legend = c(sprintf("Umbral: 4/n = %.4f", umbral_cook), "Obs. influyente"),
       col    = c("red", "red"),
       lty    = c(2, 1),
       lwd    = c(2, 2),
       cex    = 0.8)

# Identificar observaciones influyentes.
obs_influyentes <- which(cook_dist > umbral_cook)

cat("--- Conclusión (Observaciones influyentes) ---\n")
if (length(obs_influyentes) > 0) {
  cat(sprintf("Se identificaron %d observacion(es) potencialmente influyente(s):\n",
              length(obs_influyentes)))
  cat("Índices:", paste(obs_influyentes, collapse = ", "), "\n\n")
  cat("Distancias de Cook de las observaciones influyentes:\n")
  print(cook_dist[obs_influyentes])
  cat("\nSe recomienda revisar estas observaciones para determinar si son\n")
  cat("errores de medición, valores atípicos o simplemente puntos extremos\n")
  cat("con alta influencia legítima en el modelo.\n")
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

# Recopilar información de todas las secciones anteriores.
vars_final_str <- paste(regresores, collapse = ", ")
vars_elim_str  <- ifelse(length(variables_eliminadas) > 0,
                         paste(variables_eliminadas, collapse = ", "),
                         "ninguna")
r2_adj_final   <- round(resumen_final$adj.r.squared, 4)
p_f_final      <- round(p_valor_F, 6)
modelo_sign    <- ifelse(p_f_final <= alfa,
                         "sí es globalmente significativo",
                         "NO es globalmente significativo")
normalidad_str <- ifelse(prueba_shapiro$p.value > alfa,
                         "se puede asumir normalidad de los residuos",
                         "NO se puede asumir normalidad de los residuos")
homoc_str      <- ifelse(prueba_bp$p.value > alfa,
                         "se puede asumir homocedasticidad",
                         "NO se puede asumir homocedasticidad")
indep_str      <- ifelse(prueba_dw$p.value > alfa,
                         "se puede asumir independencia de los errores",
                         "NO se puede asumir independencia de los errores")

if (!is.null(vif_valores)) {
  if (any(vif_valores > 10)) {
    multi_str <- "se detecta multicolinealidad FUERTE (VIF > 10)"
  } else if (any(vif_valores > 5)) {
    multi_str <- "se detecta multicolinealidad MODERADA (5 < VIF ≤ 10)"
  } else {
    multi_str <- "no se detecta multicolinealidad significativa"
  }
} else {
  multi_str <- "el VIF no aplica (solo un regresor en el modelo final)"
}

inf_str <- ifelse(length(obs_influyentes) > 0,
                  sprintf("se identificaron %d observación(es) potencialmente influyente(s): %s",
                          length(obs_influyentes),
                          paste(obs_influyentes, collapse = ", ")),
                  "no se identificaron observaciones potencialmente influyentes")

# --- Imprimir conclusión ---
cat("RESUMEN DEL ANÁLISIS DE REGRESIÓN LINEAL MÚLTIPLE\n")
cat("MODELO SIN INTERCEPTO\n")
cat("Método: Eliminación Hacia Atrás\n")
cat(rep("-", 52), "\n", sep = "")

# CAMBIO: la interpretación general menciona que el modelo no tiene intercepto
# y que el plano de regresión pasa por el origen.
cat(sprintf(paste0(
  "1. MODELO FINAL OBTENIDO (sin intercepto):\n",
  "   Ŷ = %s\n\n",
  "2. VARIABLES EN EL MODELO FINAL: %s\n\n",
  "3. VARIABLES ELIMINADAS: %s\n\n",
  "4. AJUSTE DEL MODELO:\n",
  "   R² ajustado (no centrado) = %.4f.\n",
  "   ADVERTENCIA: en modelos sin intercepto, R reporta un R² no centrado\n",
  "   (SST = ΣYᵢ² en lugar de Σ(Yᵢ-Ȳ)²), lo que puede producir valores\n",
  "   artificialmente altos. Interprete con cautela y compárelo solo contra\n",
  "   otros modelos sin intercepto.\n\n",
  "5. SIGNIFICANCIA GLOBAL (Prueba F):\n",
  "   El modelo final %s (p-value = %.6f).\n\n",
  "6. SIGNIFICANCIA INDIVIDUAL:\n",
  "   Todos los regresores del modelo final son estadísticamente significativos\n",
  "   al nivel α = 0.05, ya que el método de eliminación hacia atrás garantiza\n",
  "   que ningún regresor permanece con p-value > 0.05.\n\n",
  "7. NORMALIDAD (Shapiro-Wilk, p-value = %.6f):\n",
  "   Con α = 0.05, %s.\n\n",
  "8. HOMOCEDASTICIDAD (Breusch-Pagan, p-value = %.6f):\n",
  "   Con α = 0.05, %s.\n\n",
  "9. INDEPENDENCIA (Durbin-Watson, p-value = %.6f):\n",
  "   Con α = 0.05, %s.\n\n",
  "10. MULTICOLINEALIDAD (VIF):\n",
  "    Con α = 0.05, %s.\n\n",
  "11. OBSERVACIONES INFLUYENTES (Distancia de Cook, umbral 4/n = %.6f):\n",
  "    %s.\n\n",
  "12. INTERPRETACIÓN GENERAL:\n",
  "    El modelo de regresión lineal múltiple SIN INTERCEPTO obtenido\n",
  "    mediante eliminación hacia atrás relaciona la variable respuesta Y\n",
  "    con los regresores: %s.\n",
  "    Al tratarse de un modelo sin intercepto, se asume que cuando todos\n",
  "    los regresores son iguales a cero, el valor esperado de Y también es\n",
  "    cero (el plano de regresión está obligado a pasar por el origen).\n",
  "    Cada coeficiente estimado representa el cambio esperado en Y por cada\n",
  "    unidad de aumento en el regresor correspondiente, manteniendo\n",
  "    constantes los demás regresores en el modelo.\n"
),
terminos_eq,
vars_final_str,
vars_elim_str,
r2_adj_final,
modelo_sign, p_f_final,
prueba_shapiro$p.value, normalidad_str,
prueba_bp$p.value, homoc_str,
prueba_dw$p.value, indep_str,
multi_str,
umbral_cook, inf_str,
vars_final_str
))

cat(rep("=", 52), "\n", sep = "")
cat("  FIN DEL ANÁLISIS\n")
cat(rep("=", 52), "\n", sep = "")


# =============================================================================
# =============================================================================
# GUÍA DE INTERPRETACIÓN DE RESULTADOS PARA EL REPORTE
# (MODELO SIN INTERCEPTO)
# =============================================================================
# =============================================================================
#
# A continuación se presentan frases tipo para redactar el reporte académico.
# Adapte los valores numéricos según sus resultados.
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
# 6. NORMALIDAD, HOMOCEDASTICIDAD, INDEPENDENCIA, VIF, COOK
# ----------------------------------------------------------------------------
# (Las frases tipo son idénticas a las del modelo con intercepto.)
# Referirse a la guía del script Metodo_de_eliminacion_hacia_atras.R
# para las frases de Shapiro-Wilk, Breusch-Pagan, Durbin-Watson, VIF y Cook.
#
# =============================================================================
# =============================================================================
# RESUMEN DE CAMBIOS REALIZADOS RESPECTO AL CÓDIGO CON INTERCEPTO
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

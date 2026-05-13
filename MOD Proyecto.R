### PROYECTO MOD: Diagnóstico y monitoreo multivariante de la eficiencia #####
### operativa en proyectos ERP mediante análisis de tickets. ###
###### Integrantes: Maria Laura Montenegro, Saraí Torres. ###

### En este script abordamos los objetivos de modelar:
### PLS-DA para priority
#### Para predecir worked hours
### PLS para predecir worked hours
### PCR para predecir worked hpurs
### MLR para predecir worked hours
### RF para predecir worked hours
 

#------------------------- Cargamos librerias -------------------------
library("readxl")       # leer archivos Excel
library("dplyr")        # manipulación de datos
library("ggplot2")      # gráficos
library("pls")          # PCR y PLS 
library("randomForest") # Random Forest
library("caret")        # framework unificado de validación y partición
library("ggrepel")      # graficos
library("gridExtra")    # combinar múltiples gráficos en un panel
library("pROC")         # curvas ROC y AUROC por clase

set.seed(42)



# ------------------------- Cargamos los datos -------------------------

ruta    <- "C:/Users/52999/OneDrive - UPV/Escritorio/Master ADMPTD/Análisis, Monitorización y Diagnóstico de Procesos Multivariantes/Proyecto"
archivo <- file.path(ruta, "Muestra de tickets.xlsx")

df_raw <- read_excel(archivo, sheet = "Muestra dragonet")

cat("Dimensiones:", nrow(df_raw), "filas x", ncol(df_raw), "columnas\n")
print(names(df_raw))


# ------------------------- Seleccion de variables -------------------------
# Seleccion segun el modelo :
# Para modelos MLR y PLS-DA: solo 5 variables desagregadas ya que MLR no maneja colinealidad y PLS-DA podria
# replicar informacion
# Y para PCR / PLS / RF: 8 variables completas ya que manejan bien la colinealidad por diseño

vars_modelo <- c( "Priority level", "Worked hours", "Inward Tickets", "Outward Tickets",
  "Inward Project", "Outward Project", "Inward", "Outward", "Total Link" )

# Eliminar Worked hours = 0 
df <- df_raw %>%
  select(all_of(vars_modelo)) %>%
  filter(`Worked hours` > 0)

cat("\nFilas tras eliminar Worked hours = 0:", nrow(df), "\n")


# Limpieza de datos (replicando lo que se hzio en Dragonet)
# Buscamos replicar el proceso que se hizo en Dragonet, es decir:
#   - Centrado y escalado (media=0, sd=1) sobre las variables X
#   - PCA con A=4 componentes (seleccionado por 7-fold CV)
#   - Eliminación iterativa: SPE > 18 Y T² en percentil 99

# Se replica exactamente este proceso para partir del mismo dataset limpio.
# Todos los modelos posteriores (MLR, PCR, PLS, RF, PLS-DA) trabajan
# sobre df_limpio con lo que podemos tener entonces modelos comparables

vars_X_pca <- setdiff(vars_modelo, "Worked hours")
X_scaled   <- scale(as.matrix(df[, vars_X_pca]))
A_pca      <- 4  # número de componentes usados en Dragonet

#Limpieza iterativa de outliers

UMBRAL_SPE <- 18
idx        <- 1:nrow(df)
iteracion  <- 1

repeat {
  X_iter    <- X_scaled[idx, ]
  p_iter    <- prcomp(X_iter, center = FALSE, scale. = FALSE)
  T_iter    <- p_iter$x[, 1:A_pca]
  P_iter    <- p_iter$rotation[, 1:A_pca]
  E_iter    <- X_iter - T_iter %*% t(P_iter)
  SPE_iter  <- rowSums(E_iter^2)
  lam_iter  <- p_iter$sdev[1:A_pca]^2
  T2_iter   <- rowSums(sweep(T_iter^2, 2, lam_iter, "/"))
  umbral_T2 <- quantile(T2_iter, 0.99)
  idx_malos <- which(SPE_iter > UMBRAL_SPE & T2_iter > umbral_T2)
  
  cat("  Iteración", iteracion,
      "| obs. activas:", length(idx),
      "| outliers detectados:", length(idx_malos), "\n")
  
  if (length(idx_malos) == 0) break
  idx       <- idx[-idx_malos]
  iteracion <- iteracion + 1
}

df_limpio <- df[idx, ]
cat("\n Dataset limpio:", nrow(df_limpio), "observaciones\n")
cat(" Eliminadas:", nrow(df) - nrow(df_limpio), "observaciones\n")



# ------------------------- Creamos laa funciones auxiliares -------------------------

#Creamos funciones de métricas para modelos de regresion
funcRMSE <- function(Ypred, Yobs, ndec = 4) {
  round(sqrt(mean((Ypred - Yobs)^2)), ndec)
}

funcMAE <- function(Ypred, Yobs, ndec = 4) {
  round(mean(abs(Ypred - Yobs)), ndec)
}

funcNMSE <- function(Ypred, Yobs, ndec = 4) {
  round(sum((Ypred - Yobs)^2) / sum((mean(Yobs) - Yobs)^2), ndec)
}

funcNMAE <- function(Ypred, Yobs, ndec = 4) {
  round(sum(abs(Ypred - Yobs)) / sum(abs(mean(Yobs) - Yobs)), ndec)
}

# medidasER
medidasER <- function(Ypred, Yobs, ndec = 4) {
  c(RMSE = funcRMSE(Ypred, Yobs, ndec),
    MAE  = funcMAE(Ypred, Yobs, ndec),
    NMSE = funcNMSE(Ypred, Yobs, ndec),
    NMAE = funcNMAE(Ypred, Yobs, ndec))
}


# VIP (Variable Importance in Projection)
# El paquete 'pls' no lo calcula directamente por eso hacemos la implementación manual.

calcular_VIP <- function(modelo, ncomp) {
  W        <- loading.weights(modelo)[, 1:ncomp, drop = FALSE]
  Q        <- Yloadings(modelo)[, 1:ncomp, drop = FALSE]
  T_       <- scores(modelo)[, 1:ncomp, drop = FALSE]
  SS       <- colSums(T_^2) * colSums(Q^2)
  SS_total <- sum(SS)
  J        <- nrow(W)
  VIP      <- sqrt(J * rowSums(sweep(W^2, 2, SS / SS_total, "*")))
  return(VIP)
}


# MCC multiclase 
calcular_MCC_multiclase <- function(cm) {
  n   <- sum(cm)
  pk  <- colSums(cm)
  tk  <- rowSums(cm)
  c_  <- sum(diag(cm))
  num <- c_ * n - sum(pk * tk)
  den <- sqrt((n^2 - sum(pk^2)) * (n^2 - sum(tk^2)))
  if (den == 0) return(0)
  return(round(num / den, 4))
}


# Usamos el criterio de parsimonia para selección de A (componentes)
# which.min() selecciona el mínimo absoluto del RMSE-CV.
# Esto puede llevar a sobreajuste, es decir  componentes adicionales con mejora
# marginal mínima

seleccionar_A_parsimonia <- function(rmse_cv_vector, verbose = TRUE) {
  rmse_vec <- rmse_cv_vector[-1]  # eliminar A=0 (intercepto)
  A_vals   <- seq_along(rmse_vec)
  min_rmse <- min(rmse_vec)
  idx_min  <- which.min(rmse_vec)
  se_rmse  <- sd(rmse_vec)
  umbral   <- min_rmse + se_rmse
  A_optimo <- A_vals[which(rmse_vec <= umbral)[1]]
  
  if (verbose) {
    cat("  RMSE-CV mínimo:", round(min_rmse, 4), "en A =", idx_min, "\n")
    cat("  Umbral parsimonia:", round(umbral, 4), "\n")
    cat("  A óptimo:", A_optimo,
        "(modelo más parsimonioso dentro del umbral)\n")
  }
  return(A_optimo)
}

######################### MODELOS #########################

######################### Modelo de clasificación para Priority ##########################

# ------------------------- PLS-DA (para priority) -------------------------

# Variable respuesta
y_clase <- factor(df_limpio$`Priority level`,
                  levels = c(1, 2, 3, 4),
                  labels = c("P4", "P3", "P2", "P1"))

cat("Distribución de prioridades:\n")
print(table(y_clase))
prop_tabla        <- prop.table(table(y_clase))
prop_mayoritaria  <- max(prop_tabla)
clase_mayoritaria <- names(which.max(prop_tabla))

cat("\nProporciones:\n")
print(round(prop_tabla, 4))
cat("\nClase mayoritaria:", clase_mayoritaria,
    "(", round(prop_mayoritaria * 100, 1), "%)\n")
cat("Umbral clasificador trivial:", round(prop_mayoritaria, 4), "\n")
cat(" El PLS-DA debe superar este umbral para ser útil\n")

# Partición train/test 

idx_train <- createDataPartition(y_clase, p = 0.70, list = FALSE)
idx_test  <- setdiff(seq_along(y_clase), idx_train)

cat("\nPartición: Train =", length(idx_train),
    "obs | Test =", length(idx_test), "obs\n")

# Variables X para PLS-DA 
vars_X_da <- c("Worked hours", "Inward Tickets", "Outward Tickets",
               "Inward Project", "Outward Project")

X_da       <- df_limpio %>% select(all_of(vars_X_da))
X_da_train <- X_da[idx_train, ]
X_da_test  <- X_da[idx_test, ]
y_da_train <- y_clase[idx_train]
y_da_test  <- y_clase[idx_test]

# Escalamos sin data leakage o sea parámetros calculados solo sobre train
prep_da    <- preProcess(X_da_train, method = c("center", "scale"))
X_da_tr_sc <- predict(prep_da, X_da_train)
X_da_ts_sc <- predict(prep_da, X_da_test)

# Codificamos como dummy a Y (variable de respuesta)
# PLS-DA requiere Y como matriz (N × K) donde K = número de clases.
clases        <- levels(y_da_train)
K             <- length(clases)
Y_dummy_train <- model.matrix(~ y_da_train - 1)
colnames(Y_dummy_train) <- clases

# Ajuste PLS-DA completo para selección de A
pls_da_full <- plsr(
  Y_dummy_train ~ as.matrix(X_da_tr_sc),
  method     = "kernelpls",
  scale      = FALSE,
  validation = "CV",
  segments   = 7
)

# Selección de A 
rmsep_da   <- RMSEP(pls_da_full, estimate = "CV")
rmse_medio <- apply(rmsep_da$val[1, , ], 2, mean)

cat("\n--- SELECCIÓN DE COMPONENTES PLS-DA ---\n")
cat("RMSE-CV medio sobre las", K, "clases por componente:\n")
print(round(rmse_medio, 4))

A_da_whichmin <- which.min(rmse_medio[-1])
A_da          <- seleccionar_A_parsimonia(rmse_medio)

cat("\nComparación de criterios:\n")
cat("  which.min (mínimo absoluto):", A_da_whichmin,
    "comps | RMSE =", round(rmse_medio[A_da_whichmin + 1], 4), "\n")
cat("  Criterio parsimonia:        ", A_da,
    "comps | RMSE =", round(rmse_medio[A_da + 1], 4), "\n")
cat("  Diferencia de RMSE:",
    round(rmse_medio[A_da + 1] - rmse_medio[A_da_whichmin + 1], 4),
    "→ mejora marginal despreciable\n")
cat("  → Se selecciona A =", A_da, "para evitar sobreajuste\n")

# Reajustamos el modelo PLS-DA con A óptimo
# Se reajusta con A_da para que scores(), loadings() y Yloadings()
# devuelvan exactamente las A componentes seleccionadas.
pls_da <- plsr(
  Y_dummy_train ~ as.matrix(X_da_tr_sc),
  ncomp  = A_da,
  method = "kernelpls",
  scale  = FALSE
)

# Predicción y clasificación
Y_pred_da  <- predict(pls_da_full,
                      newdata = as.matrix(X_da_ts_sc),
                      ncomp   = A_da)[, , 1]
clase_pred <- clases[apply(Y_pred_da, 1, which.max)]
clase_real <- as.character(y_da_test)

# Hacemos el ajuste bayesiano por probabilidades a priori ya que tenemos clases desbalanceadas
# Dividir por la proporción de cada clase en train penaliza las clases
# frecuentes (como puede ser P3) y da oportunidad real a P1 y P2 — equivalente al umbral
# bayesiano óptimo: clasificar según P(clase|x) en lugar de P(x|clase).
prior            <- as.numeric(prop.table(table(y_da_train)))
names(prior)     <- clases
Y_pred_bayes     <- sweep(Y_pred_da, 2, prior, "/")
clase_pred_bayes <- clases[apply(Y_pred_bayes, 1, which.max)]

# Métricas de clasificación

acc_da       <- mean(clase_pred       == clase_real)
acc_da_bayes <- mean(clase_pred_bayes == clase_real)

cat("Accuracy sin ajuste:          ", round(acc_da,           4), "\n")
cat("Accuracy con ajuste bayesiano:", round(acc_da_bayes,     4), "\n")
cat("Umbral trivial (", clase_mayoritaria, "):   ",
    round(prop_mayoritaria, 4), "\n")
cat("Mejora sin ajuste:            ",
    round((acc_da       - prop_mayoritaria) * 100, 2), "pp\n")
cat("Mejora con ajuste bayesiano:  ",
    round((acc_da_bayes - prop_mayoritaria) * 100, 2), "pp\n")

# Matrices de confusión
cm_da       <- table(Real = clase_real, Predicho = clase_pred)
cm_da_bayes <- table(Real = clase_real, Predicho = clase_pred_bayes)

cat("\nMatriz de confusión — sin ajuste:\n")
print(cm_da)
cat("\nMatriz de confusión — con ajuste bayesiano:\n")
print(cm_da_bayes)

# Sensibilidad y especificidad por clase
cat("\nSin ajuste bayesiano:\n")
cat("Sensibilidad (TP / (TP + FN)):\n")
for (k in clases) {
  if (!(k %in% rownames(cm_da))) { cat("  ", k, ": NA\n"); next }
  tp  <- if (k %in% colnames(cm_da)) cm_da[k, k] else 0
  fn  <- sum(cm_da[k, ]) - tp
  sen <- ifelse((tp + fn) == 0, NA, round(tp / (tp + fn), 4))
  cat("  ", k, ":", sen, "\n")
}
cat("Especificidad (TN / (TN + FP)):\n")
for (k in clases) {
  if (!(k %in% rownames(cm_da))) { cat("  ", k, ": NA\n"); next }
  tp  <- if (k %in% colnames(cm_da)) cm_da[k, k] else 0
  fp  <- if (k %in% colnames(cm_da)) sum(cm_da[, k]) - tp else 0
  fn  <- sum(cm_da[k, ]) - tp
  tn  <- sum(cm_da) - tp - fp - fn
  esp <- ifelse((tn + fp) == 0, NA, round(tn / (tn + fp), 4))
  cat("  ", k, ":", esp, "\n")
}

cat("\nCon ajuste bayesiano:\n")
cat("Sensibilidad (TP / (TP + FN)):\n")
for (k in clases) {
  if (!(k %in% rownames(cm_da_bayes))) { cat("  ", k, ": NA\n"); next }
  tp  <- if (k %in% colnames(cm_da_bayes)) cm_da_bayes[k, k] else 0
  fn  <- sum(cm_da_bayes[k, ]) - tp
  sen <- ifelse((tp + fn) == 0, NA, round(tp / (tp + fn), 4))
  cat("  ", k, ":", sen, "\n")
}
cat("Especificidad (TN / (TN + FP)):\n")
for (k in clases) {
  if (!(k %in% rownames(cm_da_bayes))) { cat("  ", k, ": NA\n"); next }
  tp  <- if (k %in% colnames(cm_da_bayes)) cm_da_bayes[k, k] else 0
  fp  <- if (k %in% colnames(cm_da_bayes)) sum(cm_da_bayes[, k]) - tp else 0
  fn  <- sum(cm_da_bayes[k, ]) - tp
  tn  <- sum(cm_da_bayes) - tp - fp - fn
  esp <- ifelse((tn + fp) == 0, NA, round(tn / (tn + fp), 4))
  cat("  ", k, ":", esp, "\n")
}

# MCC
mcc_da       <- calcular_MCC_multiclase(cm_da)
mcc_da_bayes <- calcular_MCC_multiclase(cm_da_bayes)
cat("\nMCC sin ajuste:          ", mcc_da,       "\n")
cat("MCC con ajuste bayesiano:", mcc_da_bayes, "\n")

# AUROC por clase (one-vs-rest)
cat("\nAUROC por clase (one-vs-rest):\n")
auroc_por_clase <- numeric(K)
names(auroc_por_clase) <- clases

for (k in seq_len(K)) {
  clase_k            <- clases[k]
  y_bin              <- as.integer(clase_real == clase_k)
  score_k            <- if (is.matrix(Y_pred_da)) Y_pred_da[, k] else Y_pred_da
  roc_k              <- roc(y_bin, score_k, quiet = TRUE)
  auroc_por_clase[k] <- round(auc(roc_k), 4)
  cat("  ", clase_k, ":", auroc_por_clase[k], "\n")
}
cat("  AUROC medio:", round(mean(auroc_por_clase), 4), "\n")

# Permutation test
# Si permutar las etiquetas produce accuracy similar al modelo real,
# entonces el modelo NO captura estructura real ;(
# p-valor = proporción de permutaciones con accuracy >= accuracy real.
cat("\n--- PERMUTATION TEST (B=200) ---\n")
B        <- 200
acc_perm <- numeric(B)

for (b in seq_len(B)) {
  y_perm       <- sample(y_da_train)
  Y_perm_dummy <- model.matrix(~ y_perm - 1)
  colnames(Y_perm_dummy) <- clases
  mod_perm     <- plsr(Y_perm_dummy ~ as.matrix(X_da_tr_sc),
                       ncomp = A_da, method = "kernelpls", scale = FALSE)
  Y_pred_p     <- predict(mod_perm,
                          newdata = as.matrix(X_da_ts_sc),
                          ncomp   = A_da)[, , 1]
  acc_perm[b]  <- mean(clases[apply(Y_pred_p, 1, which.max)] == clase_real)
}

p_valor_perm <- mean(acc_perm >= acc_da)
cat("Accuracy real:     ", round(acc_da,                    4), "\n")
cat("Media permutada:   ", round(mean(acc_perm),            4), "\n")
cat("Percentil 95 nulo: ", round(quantile(acc_perm, 0.95),  4), "\n")
cat("p-valor:           ", round(p_valor_perm,              4), "\n")

if (p_valor_perm < 0.05) {
  cat("El modelo captura estructura real (p < 0.05)\n")
} else {
  cat(" No se puede descartar resultado aleatorio (p >= 0.05)\n")
}

# Interpretación
cat("\n--- INTERPRETACIÓN SISTEMA DE PRIORIZACIÓN ---\n")
if (acc_da > prop_mayoritaria && p_valor_perm < 0.05) {
  cat("PLS-DA SUPERA al clasificador trivial Y es estadísticamente\n")
  cat("significativo → existe coherencia PARCIAL entre comportamiento\n")
  cat("y priorización.\n")
} else if (acc_da > prop_mayoritaria && p_valor_perm >= 0.05) {
  cat("PLS-DA supera al clasificador trivial PERO el permutation test\n")
  cat("no es significativo → resultado posiblemente espurio.\n")
} else {
  cat("PLS-DA NO supera al clasificador trivial.\n")
  cat(" La prioridad NO es predecible desde el comportamiento real.\n")
}

# Gráficos PLS-DA
ruta_out <- ruta

# G8: RMSE-CV PLS-DA — justificación visual de A
rmse_vec_plot <- rmse_medio[-1]
A_max_plot    <- length(rmse_vec_plot)
umbral_plot   <- min(rmse_vec_plot) + sd(rmse_vec_plot)

df_rmse_da <- data.frame(A = seq_len(A_max_plot), RMSE = rmse_vec_plot)

g8 <- ggplot(df_rmse_da, aes(x = A, y = RMSE)) +
  geom_line(color = "#9C27B0", linewidth = 1.2) +
  geom_point(color = "#9C27B0", size = 2.5) +
  geom_hline(yintercept = umbral_plot, linetype = "dashed",
             color = "gray50", linewidth = 0.8) +
  geom_vline(xintercept = A_da, linetype = "dashed",
             color = "#F44336", linewidth = 0.8) +
  geom_vline(xintercept = A_da_whichmin, linetype = "dotted",
             color = "gray40", linewidth = 0.8) +
  annotate("text", x = A_da + 0.15, y = max(rmse_vec_plot) * 0.97,
           label = paste0("A = ", A_da, " (seleccionado)"),
           color = "#F44336", hjust = 0, size = 3.8) +
  annotate("text", x = A_da_whichmin + 0.15, y = max(rmse_vec_plot) * 0.91,
           label = paste0("A = ", A_da_whichmin, " (mín.)"),
           color = "gray40", hjust = 0, size = 3.5) +
  annotate("text", x = A_max_plot * 0.6, y = umbral_plot + 0.001,
           label = "Umbral parsimonia", color = "gray50",
           hjust = 0, size = 3.3) +
  scale_x_continuous(breaks = seq_len(A_max_plot)) +
  labs(title    = "Selección de componentes PLS-DA",
       subtitle = paste0("RMSE-CV medio sobre ", K,
                         " clases | Criterio parsimonia → A = ", A_da),
       x = "Número de componentes (A)", y = "RMSE-CV medio") +
  theme_minimal(base_size = 12)
print(g8)

# G9: Score plot PLS-DA
sc_da    <- scores(pls_da)[, 1:min(2, A_da), drop = FALSE]
df_sc_da <- data.frame(
  LV1       = sc_da[, 1],
  LV2       = if (ncol(sc_da) >= 2) sc_da[, 2] else 0,
  Prioridad = y_da_train
)

g9 <- ggplot(df_sc_da, aes(x = LV1, y = LV2, color = Prioridad)) +
  geom_point(alpha = 0.4, size = 1.5) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_manual(values = c("P1" = "#9C27B0", "P2" = "#F44336",
                                "P3" = "#2196F3",  "P4" = "#4CAF50")) +
  labs(title    = paste0("Score plot PLS-DA (A = ", A_da, ")"),
       subtitle = "Elipses 95% | Solapamiento → baja separabilidad entre clases",
       x = "Score LV1", y = "Score LV2") +
  theme_minimal(base_size = 12)
print(g9)

# G10: W*C plot PLS-DA
wstar_da <- as.data.frame(
  loading.weights(pls_da)[, 1:min(2, A_da), drop = FALSE])
colnames(wstar_da) <- c("LV1", if (ncol(wstar_da) >= 2) "LV2" else NULL)
if (ncol(wstar_da) < 2) wstar_da$LV2 <- 0
wstar_da$Variable <- gsub("as.matrix\\(X_da_tr_sc\\)", "",
                          rownames(wstar_da))
wstar_da$Tipo <- "X (W*)"

c_da <- as.data.frame(
  Yloadings(pls_da)[, 1:min(2, A_da), drop = FALSE])
colnames(c_da) <- c("LV1", if (ncol(c_da) >= 2) "LV2" else NULL)
if (ncol(c_da) < 2) c_da$LV2 <- 0
c_da$Variable <- rownames(c_da)
c_da$Tipo     <- "Y — Clase (C)"

df_wc <- rbind(wstar_da, c_da)

g10 <- ggplot(df_wc, aes(x = LV1, y = LV2, label = Variable, color = Tipo)) +
  geom_point(size = 3.5) +
  geom_text_repel(size = 3.5, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  scale_color_manual(values = c("X (W*)"        = "#2196F3",
                                "Y — Clase (C)" = "#F44336")) +
  labs(title    = paste0("W*C plot — PLS-DA (A = ", A_da, ")"),
       subtitle = "Azul: variables X (W*) | Rojo: clases Y (C) | Proximidad = asociación",
       x = "LV1", y = "LV2", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
print(g10)

# G11: Curvas ROC por clase
png(file.path(ruta_out, "g11_roc_plsda.png"),
    width = 1000, height = 900, res = 130)
par(mfrow = c(2, 2))
colores_roc <- c("#9C27B0", "#F44336", "#2196F3", "#4CAF50")
for (k in seq_len(K)) {
  clase_k <- clases[k]
  y_bin   <- as.integer(clase_real == clase_k)
  score_k <- if (is.matrix(Y_pred_da)) Y_pred_da[, k] else Y_pred_da
  roc_k   <- roc(y_bin, score_k, quiet = TRUE)
  plot(roc_k,
       main = paste0("ROC — ", clase_k,
                     " | AUROC = ", round(auc(roc_k), 3)),
       col  = colores_roc[k], lwd = 2)
  abline(a = 0, b = 1, lty = 2, col = "gray60")
}
par(mfrow = c(1, 1))
dev.off()

# G12: Permutation test
df_perm <- data.frame(acc_permutada = acc_perm)

g12 <- ggplot(df_perm, aes(x = acc_permutada)) +
  geom_histogram(bins = 30, fill = "#90A4AE", color = "white", alpha = 0.8) +
  geom_vline(xintercept = acc_da, color = "#F44336", linewidth = 1.2) +
  geom_vline(xintercept = quantile(acc_perm, 0.95),
             color = "gray30", linewidth = 0.8, linetype = "dashed") +
  annotate("text",
           x     = acc_da + diff(range(acc_perm)) * 0.02,
           y     = B * 0.07,
           label = paste0("Acc. real\n", round(acc_da, 3)),
           color = "#F44336", hjust = 0, size = 3.5) +
  annotate("text",
           x     = quantile(acc_perm, 0.95) + diff(range(acc_perm)) * 0.02,
           y     = B * 0.12,
           label = "p = 0.05",
           color = "gray30", hjust = 0, size = 3.3) +
  labs(title    = paste0("Permutation test PLS-DA (B = ", B, ")"),
       subtitle = paste0("p-valor = ", round(p_valor_perm, 4),
                         " | Rojo = accuracy modelo real"),
       x = "Accuracy bajo permutación aleatoria", y = "Frecuencia") +
  theme_minimal(base_size = 12)
print(g12)

# G13: Matrices de confusión comparativas
df_cm <- as.data.frame(cm_da)
colnames(df_cm) <- c("Real", "Predicho", "Freq")

g13a <- ggplot(df_cm, aes(x = Predicho, y = Real, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5, color = "white", fontface = "bold") +
  scale_fill_gradient(low = "#BBDEFB", high = "#1565C0") +
  labs(title    = "Matriz de confusión — sin ajuste",
       subtitle = paste0("Accuracy = ", round(acc_da, 3),
                         " | MCC = ", mcc_da),
       x = "Clase predicha", y = "Clase real", fill = "N obs.") +
  theme_minimal(base_size = 12)

df_cm_bayes <- as.data.frame(cm_da_bayes)
colnames(df_cm_bayes) <- c("Real", "Predicho", "Freq")

g13b <- ggplot(df_cm_bayes, aes(x = Predicho, y = Real, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5, color = "white", fontface = "bold") +
  scale_fill_gradient(low = "#BBDEFB", high = "#1565C0") +
  labs(title    = "Matriz de confusión — ajuste bayesiano",
       subtitle = paste0("Accuracy = ", round(acc_da_bayes, 3),
                         " | MCC = ", mcc_da_bayes),
       x = "Clase predicha", y = "Clase real", fill = "N obs.") +
  theme_minimal(base_size = 12)

print(g13a)
print(g13b)

# Guardar gráficos PLS-DA
ggsave(file.path(ruta_out, "g8_rmse_cv_plsda.png"),
       g8,  width = 7, height = 5, dpi = 150)
ggsave(file.path(ruta_out, "g9_score_plsda.png"),
       g9,  width = 7, height = 6, dpi = 150)
ggsave(file.path(ruta_out, "g10_wc_plsda.png"),
       g10, width = 7, height = 6, dpi = 150)
ggsave(file.path(ruta_out, "g12_permutation_test.png"),
       g12, width = 8, height = 5, dpi = 150)
ggsave(file.path(ruta_out, "g13_confusion_plsda.png"),
       grid.arrange(g13a, g13b, ncol = 2,
                    top = "Sin ajuste vs Ajuste bayesiano"),
       width = 12, height = 5, dpi = 150)



##################### Modelos de regresion para Worked Hours ###############################

# ------------------------- Preparacion de variables -------------------------

# Se reutilizan idx_train e idx_test generados anteriormente
# esto garantiza que todas las comparaciones son sobre el mismo test.


# Variables Y e X para regresión
y <- df_limpio$`Worked hours`

# X para MLR: sin variables derivadas
X_mlr <- df_limpio %>%
  select("Priority level", "Inward Tickets", "Outward Tickets",
         "Inward Project", "Outward Project")

# X completa: para PLS, PCR y RF ya que absorben colinealidad por diseño
X_full <- df_limpio %>% select(-`Worked hours`)

cat("Variables X_mlr (", ncol(X_mlr), "):", names(X_mlr), "\n")
cat("Variables X_full (", ncol(X_full), "):", names(X_full), "\n")

# Subconjuntos train/test
X_full_train <- X_full[idx_train, ]
X_full_test  <- X_full[idx_test, ]
X_mlr_train  <- X_mlr[idx_train, ]
X_mlr_test   <- X_mlr[idx_test, ]
y_train      <- y[idx_train]
y_test       <- y[idx_test]

cat("Train:", length(y_train), "obs | Test:", length(y_test), "obs\n")

# Escalado con parámetros SOLO del train para evitar data leakage
prep_full    <- preProcess(X_full_train, method = c("center", "scale"))
X_full_tr_sc <- predict(prep_full, X_full_train)
X_full_ts_sc <- predict(prep_full, X_full_test)

prep_mlr    <- preProcess(X_mlr_train, method = c("center", "scale"))
X_mlr_tr_sc <- predict(prep_mlr, X_mlr_train)
X_mlr_ts_sc <- predict(prep_mlr, X_mlr_test)


# ------------------------- MLR -------------------------

#Ajustamos el modelo
mlr_mod  <- lm(y_train ~ ., data = as.data.frame(X_mlr_tr_sc))
pred_mlr <- predict(mlr_mod, newdata = as.data.frame(X_mlr_ts_sc))

cat("\nResumen MLR:\n")
print(summary(mlr_mod))

cat("\nMétricas MLR en test:\n")
res_mlr <- medidasER(pred_mlr, y_test)
print(res_mlr)

cat("\nCoeficientes ordenados por |valor|:\n")
coefs_mlr <- sort(abs(coef(mlr_mod)[-1]), decreasing = TRUE)
print(round(coefs_mlr, 4))

# Verificación de supuestos
# 1. Residuals vs Fitted: verifica linealidad
# 2. Q-Q Residuals: verifica normalidad de residuos
# 3. Scale-Location: verifica homocedasticidad
# 4. Residuals vs Leverage: detecta observaciones influyentes
par(mfrow = c(2, 2))
plot(mlr_mod, main = "MLR — Verificación de supuestos")
par(mfrow = c(1, 1))
# ------------------------- PCR -------------------------
# Ajustamos el modelo
pcr_mod <- pcr(
  y_train ~ as.matrix(X_full_tr_sc),
  scale      = FALSE,
  validation = "CV",
  segments   = 7
)

rmsep_pcr <- RMSEP(pcr_mod, estimate = "CV")

cat("Selección de componentes PCR:\n")
A_pcr <- seleccionar_A_parsimonia(as.numeric(rmsep_pcr$val[1, 1, ]))

cat("\nRMSE-CV por componente:\n")
print(round(rmsep_pcr$val[1, 1, ], 3))

pred_pcr <- predict(pcr_mod,
                    newdata = as.matrix(X_full_ts_sc),
                    ncomp   = A_pcr)[, 1, 1]

cat("\nMétricas PCR en test:\n")
res_pcr <- medidasER(pred_pcr, y_test)
print(res_pcr)
# ------------------------- PLS -------------------------
# Ajustamos el modelo
pls_mod <- plsr(
  y_train ~ as.matrix(X_full_tr_sc),
  method     = "kernelpls",
  scale      = FALSE,
  validation = "CV",
  segments   = 7
)

rmsep_pls <- RMSEP(pls_mod, estimate = "CV")

cat("Selección de componentes PLS:\n")
A_pls <- seleccionar_A_parsimonia(as.numeric(rmsep_pls$val[1, 1, ]))

cat("\nRMSE-CV por componente:\n")
print(round(rmsep_pls$val[1, 1, ], 3))

pred_pls <- predict(pls_mod,
                    newdata = as.matrix(X_full_ts_sc),
                    ncomp   = A_pls)[, 1, 1]

cat("\nMétricas PLS en test:\n")
res_pls <- medidasER(pred_pls, y_test)
print(res_pls)

cat("\nPesos W* (componente 1):\n")
wstar <- loading.weights(pls_mod)[, 1:min(2, A_pls), drop = FALSE]
rownames(wstar) <- gsub("as.matrix\\(X_full_tr_sc\\)", "", rownames(wstar))
print(round(wstar, 4))

# VIP
vip_vals <- calcular_VIP(pls_mod, A_pls)
names(vip_vals) <- gsub("as.matrix\\(X_full_tr_sc\\)", "", names(vip_vals))
cat("\nVIP (A =", A_pls, "):\n")
print(round(sort(vip_vals, decreasing = TRUE), 4))
cat("Variables con VIP > 1:", names(vip_vals[vip_vals > 1]), "\n")


# Gráficos PLS
# G1: RMSE-CV PLS vs PCR, comparación de selección de A
df_rmse <- rbind(
  data.frame(A      = seq_along(rmsep_pls$val[1, 1, ]) - 1,
             RMSE   = as.numeric(rmsep_pls$val[1, 1, ]),
             Modelo = "PLS"),
  data.frame(A      = seq_along(rmsep_pcr$val[1, 1, ]) - 1,
             RMSE   = as.numeric(rmsep_pcr$val[1, 1, ]),
             Modelo = "PCR")
) %>% filter(A > 0)

g1 <- ggplot(df_rmse, aes(x = A, y = RMSE, color = Modelo)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = A_pls, linetype = "dashed",
             color = "#2196F3", linewidth = 0.8) +
  geom_vline(xintercept = A_pcr, linetype = "dashed",
             color = "#F44336", linewidth = 0.8) +
  annotate("text", x = A_pls + 0.15, y = max(df_rmse$RMSE) * 0.98,
           label = paste0("A_PLS=", A_pls), color = "#2196F3",
           hjust = 0, size = 3.5) +
  annotate("text", x = A_pcr + 0.15, y = max(df_rmse$RMSE) * 0.93,
           label = paste0("A_PCR=", A_pcr), color = "#F44336",
           hjust = 0, size = 3.5) +
  scale_color_manual(values = c("PLS" = "#2196F3", "PCR" = "#F44336")) +
  labs(title    = "Selección de componentes: PLS vs PCR",
       subtitle = "RMSE-CV 7-fold | Criterio parsimonia | Líneas = A óptimo",
       x = "Número de componentes (A)", y = "RMSE-CV (horas)") +
  theme_minimal(base_size = 12)
print(g1)

# G2: VIP plot
names(vip_vals) <- gsub("as.matrix\\(X_full_tr_sc\\)", "", names(vip_vals))
df_vip <- data.frame(
  Variable  = names(vip_vals),
  VIP       = as.numeric(vip_vals),
  Relevante = ifelse(vip_vals > 1, "VIP > 1", "VIP ≤ 1")
) %>% arrange(desc(VIP))

g2 <- ggplot(df_vip, aes(x = VIP, y = reorder(Variable, VIP),
                         fill = Relevante)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "black", linewidth = 0.8) +
  scale_fill_manual(values = c("VIP > 1" = "#F44336",
                               "VIP ≤ 1" = "#90A4AE")) +
  labs(title    = paste0("VIP — PLS (A = ", A_pls, ")"),
       subtitle = "VIP > 1 → variable relevante para predecir Worked hours",
       x = "VIP", y = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
print(g2)

# G3: W* plot
wstar_df <- as.data.frame(loading.weights(pls_mod)[, 1:min(2, A_pls)])
colnames(wstar_df) <- c("LV1", if (A_pls >= 2) "LV2" else NULL)
if (ncol(wstar_df) < 2) wstar_df$LV2 <- 0
wstar_df$Variable <- gsub("as.matrix\\(X_full_tr_sc\\)", "",
                          rownames(wstar_df))

g3 <- ggplot(wstar_df, aes(x = LV1, y = LV2, label = Variable)) +
  geom_point(color = "#F44336", size = 3.5) +
  geom_text_repel(size = 3.5, max.overlaps = 20) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  labs(title    = paste0("W* plot — PLS (A = ", A_pls, ")"),
       subtitle = "Pesos W* = W(P'W)⁻¹ | Relación neta X original → scores",
       x = "W* LV1", y = "W* LV2") +
  theme_minimal(base_size = 12)
print(g3)

# G4: t vs u
t1     <- scores(pls_mod)[, 1]
u1     <- Yscores(pls_mod)[, 1]
df_tu  <- data.frame(t = t1, u = u1)
cor_tu <- round(cor(t1, u1), 3)
cat("Correlación t1 vs u1:", cor_tu, "\n")

g4 <- ggplot(df_tu, aes(x = t, y = u)) +
  geom_point(alpha = 0.4, size = 1.5, color = "#2196F3") +
  geom_smooth(method = "lm", color = "#F44336", se = TRUE, linewidth = 1) +
  labs(title    = "Relación interna PLS — t₁ vs u₁",
       subtitle = paste0("Cor(t₁, u₁) = ", cor_tu,
                         " | Alta correlación → modelo X→Y robusto"),
       x = "Score t₁ (espacio X)", y = "Score u₁ (espacio Y)") +
  theme_minimal(base_size = 12)
print(g4)

# Guardamos los graficos
ggsave(file.path(ruta_out, "g1_rmse_cv_pls_pcr.png"),
       g1, width = 8, height = 5, dpi = 150)
ggsave(file.path(ruta_out, "g2_vip_pls.png"),
       g2, width = 7, height = 5, dpi = 150)
ggsave(file.path(ruta_out, "g3_wstar_pls.png"),
       g3, width = 7, height = 6, dpi = 150)
ggsave(file.path(ruta_out, "g4_t_vs_u_pls.png"),
       g4, width = 7, height = 6, dpi = 150)




# ------------------------- RF -------------------------

#Ajustamos el mdoelo
rf_mod <- randomForest(
  x          = as.matrix(X_full[idx_train, ]),
  y          = y_train,
  ntree      = 300,
  mtry       = floor(sqrt(ncol(X_full))),
  nodesize   = 5,
  importance = TRUE
)

cat("Error OOB (RMSE interno):", round(sqrt(rf_mod$mse[300]), 3), "h\n")

pred_rf <- predict(rf_mod, newdata = as.matrix(X_full[-idx_train, ]))

cat("\nMétricas RF en test:\n")
res_rf <- medidasER(pred_rf, y_test)
print(res_rf)

imp_rf <- importance(rf_mod, type = 1)
cat("\nImportancia variables RF (%IncMSE):\n")
print(round(sort(imp_rf[, 1], decreasing = TRUE), 3))

# varImpPlot
varImpPlot(rf_mod,
           main = "Random Forest — Importancia de variables (%IncMSE)",
           col  = "darkblue")

# G7: RF %IncMSE vs VIP PLS
imp_rf_norm <- imp_rf[, 1] / max(abs(imp_rf[, 1]))
vip_norm    <- vip_vals    / max(vip_vals)

df_imp_comp <- rbind(
  data.frame(Variable    = names(imp_rf_norm),
             Importancia = as.numeric(imp_rf_norm),
             Metodo      = "RF (%IncMSE)"),
  data.frame(Variable    = names(vip_norm),
             Importancia = as.numeric(vip_norm),
             Metodo      = "PLS (VIP)")
)

g7 <- ggplot(df_imp_comp,
             aes(x = Importancia, y = reorder(Variable, Importancia),
                 fill = Metodo)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("RF (%IncMSE)" = "#9C27B0",
                               "PLS (VIP)"    = "#F44336")) +
  labs(title    = "Importancia: RF vs PLS (normalizado)",
       subtitle = "Convergencia → importancia robusta | Divergencia → no linealidades",
       x = "Importancia normalizada", y = NULL, fill = "Método") +
  theme_minimal(base_size = 12)
print(g7)

# Guardamos el grafico
ggsave(file.path(ruta_out, "g7_rf_vs_vip.png"),
       g7, width = 8, height = 5, dpi = 150)


#------------------------- Comparacion de modelos de regresión -------------------------

tabla_reg <- rbind(
  MLR             = res_mlr,
  PCR             = res_pcr,
  PLS             = res_pls,
  `Random Forest` = res_rf
)
rownames(tabla_reg)[2] <- paste0("PCR (A=", A_pcr, ")")
rownames(tabla_reg)[3] <- paste0("PLS (A=", A_pls, ")")
tabla_reg  <- tabla_reg[order(tabla_reg[, "RMSE"]), ]
mejor_reg  <- rownames(tabla_reg)[1]

print(tabla_reg)
cat("\n Mejor modelo por RMSE:", mejor_reg, "\n")
cat("  NMSE < 1: supera predecir siempre la media\n")
cat("  Si PLS ≈ RF: relación esencialmente lineal pero PLS preferible\n")

# G5: Observed vs Predicted de todos los modelos
modelos_pred <- list(pred_mlr, pred_pcr, pred_pls, pred_rf)
names(modelos_pred) <- c(
  "MLR",
  paste0("PCR (A=", A_pcr, ")"),
  paste0("PLS (A=", A_pls, ")"),
  "Random Forest"
)

plots_ovp <- lapply(names(modelos_pred), function(nm) {
  yp   <- modelos_pred[[nm]]
  rmse <- round(funcRMSE(yp, y_test), 2)
  nmse <- round(funcNMSE(yp, y_test), 3)
  ggplot(data.frame(Pred = yp, Obs = y_test), aes(x = Pred, y = Obs)) +
    geom_point(alpha = 0.3, size = 0.9, color = "#2196F3") +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "black") +
    labs(title    = nm,
         subtitle = paste0("RMSE=", rmse, "h  NMSE=", nmse),
         x = "Predicho (h)", y = "Observado (h)") +
    theme_minimal(base_size = 9)
})

# G6: Comparación RMSE y NMSE
df_tabla <- data.frame(
  Modelo = rownames(tabla_reg),
  RMSE   = tabla_reg[, "RMSE"],
  NMSE   = tabla_reg[, "NMSE"]
)
df_tabla$Modelo <- factor(df_tabla$Modelo,
                          levels = df_tabla$Modelo[order(df_tabla$RMSE)])

g6a <- ggplot(df_tabla, aes(x = RMSE, y = Modelo, fill = Modelo)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f h", RMSE)), hjust = -0.1, size = 3.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "RMSE — Test", x = "RMSE (horas)", y = NULL) +
  xlim(0, max(df_tabla$RMSE) * 1.2) +
  theme_minimal(base_size = 11)

g6b <- ggplot(df_tabla, aes(x = NMSE, y = Modelo, fill = Modelo)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", NMSE)), hjust = -0.1, size = 3.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  scale_fill_brewer(palette = "Set2") +
  labs(title    = "NMSE — Test",
       subtitle = "NMSE < 1 = supera predecir la media",
       x = "NMSE", y = NULL) +
  xlim(0, max(df_tabla$NMSE) * 1.2) +
  theme_minimal(base_size = 11)

print(g6a)
print(g6b)

# G14: Residuos del mejor modelo
pred_mejor <- modelos_pred[[mejor_reg]]
df_resid   <- data.frame(
  Predicho = pred_mejor,
  Residuo  = y_test - pred_mejor
)

g14 <- ggplot(df_resid, aes(x = Predicho, y = Residuo)) +
  geom_point(alpha = 0.3, size = 0.9, color = "#2196F3") +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "black", linewidth = 1) +
  geom_smooth(method = "loess", color = "#F44336",
              se = FALSE, linewidth = 1) +
  labs(title    = paste("Análisis de residuos —", mejor_reg),
       subtitle = "Residuos sin estructura → modelo adecuado",
       x = "Predicho (h)", y = "Residuo (h)") +
  theme_minimal(base_size = 12)
print(g14)

# Guardar gráficos regresión
ggsave(file.path(ruta_out, "g6_comparacion_modelos.png"),
       grid.arrange(g6a, g6b, ncol = 2),
       width = 12, height = 5, dpi = 150)
ggsave(file.path(ruta_out, "g14_residuos.png"),
       g14, width = 8, height = 5, dpi = 150)

png(file.path(ruta_out, "g5_obs_vs_pred.png"),
    width = 1200, height = 900, res = 130)
grid.arrange(grobs = plots_ovp, ncol = 2,
             top   = "Observed vs Predicted — Test (30%)")
dev.off()



# ------------------------- Resumen final :) -------------------------

cat("\n", strrep("=", 60), "\n")
cat("RESUMEN FINAL\n")
cat(strrep("=", 60), "\n")

cat("\n--- PLS-DA (Priority level) ---\n")
cat("Componentes seleccionados (parsimonia): A =", A_da, "\n")
cat("  (which.min daría A =", A_da_whichmin,
    "— mejora de solo",
    round(rmse_medio[A_da_whichmin + 1] - rmse_medio[A_da + 1], 4),
    "en RMSE-CV)\n")
cat("Accuracy sin ajuste:          ", round(acc_da,              4), "\n")
cat("Accuracy con ajuste bayesiano:", round(acc_da_bayes,        4), "\n")
cat("Umbral trivial:               ", round(prop_mayoritaria,    4), "\n")
cat("MCC sin ajuste:               ", mcc_da,                        "\n")
cat("MCC con ajuste bayesiano:     ", mcc_da_bayes,                  "\n")
cat("AUROC medio:                  ", round(mean(auroc_por_clase), 4), "\n")
cat("p-valor permutación:          ", round(p_valor_perm,         4), "\n")

cat("\n--- REGRESIÓN (Worked hours) ---\n")
print(tabla_reg)
cat("\n▶ Mejor modelo por RMSE:", mejor_reg, "\n")

#' Media (promedio) con NA ignorados
#' @param x Vector numérico.
#' @return Valor medio de `x`.
#' @examples
#' media(1:5)
#' @export
media <- function(x){
  stopifnot(is.numeric(x))
  mean(x, na.rm = TRUE)
}

#' Varianza con NA ignorados
#' @param x Vector numérico.
#' @return Varianza muestral de `x`.
#' @examples
#' varianza(1:5)
#' @export
varianza <- function(x){
  stopifnot(is.numeric(x))
  stats::var(x, na.rm = TRUE)
}

#' Desviación típica con NA ignorados
#' @param x Vector numérico.
#' @return Desviación estándar muestral de `x`.
#' @examples
#' desviacion(1:5)
#' @export
desviacion <- function(x){
  stopifnot(is.numeric(x))
  stats::sd(x, na.rm = TRUE)
}

#' Cuantiles
#' @param x Vector numérico.
#' @param probs Probabilidades (0-1). Por defecto 0.25, 0.5 y 0.75.
#' @return Vector con cuantiles.
#' @examples
#' cuantiles(1:9)
#' @export
cuantiles <- function(x, probs = c(0.25, 0.5, 0.75)){
  stopifnot(is.numeric(x))
  stats::quantile(x, probs = probs, na.rm = TRUE, names = TRUE)
}

#' Normalización [0,1]
#' @param x Vector numérico.
#' @return Vector reescalado a [0,1]. Si min==max devuelve ceros.
#' @examples
#' normalizar(1:10)
#' @export
normalizar <- function(x){
  stopifnot(is.numeric(x))
  rng <- range(x, na.rm = TRUE)
  if (isTRUE(all.equal(rng[1], rng[2]))) return(ifelse(is.na(x), NA_real_, 0))
  (x - rng[1]) / (rng[2] - rng[1])
}

#' Estandarización (z-score)
#' @param x Vector numérico.
#' @return (x - media)/sd con NA ignorados.
#' @examples
#' estandarizar(1:10)
#' @export
estandarizar <- function(x){
  stopifnot(is.numeric(x))
  m <- mean(x, na.rm = TRUE)
  s <- stats::sd(x, na.rm = TRUE)
  if (isTRUE(all.equal(s, 0))) return(ifelse(is.na(x), NA_real_, 0))
  (x - m) / s
}

#' Correlación rápida
#' @param x,y Vectores numéricos.
#' @param metodo "pearson", "spearman" o "kendall".
#' @return Correlación entre `x` y `y`.
#' @examples
#' correlacion(1:10, 1:10)
#' @export
correlacion <- function(x, y, metodo = c("pearson","spearman","kendall")){
  metodo <- match.arg(metodo)
  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y))
  stats::cor(x, y, use = "complete.obs", method = metodo)
}

#' Regresión lineal simple
#' @param x Predictor numérico.
#' @param y Respuesta numérica.
#' @return Objeto `lm`.
#' @examples
#' set.seed(1); x <- 1:10; y <- x + rnorm(10); m <- regresion_lineal(x,y); summary(m)
#' @export
regresion_lineal <- function(x, y){
  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y))
  stats::lm(y ~ x)
}

#' Gráfica rápida (dispersión o línea / histograma)
#' @param x Vector numérico.
#' @param y Opcional: si se pasa, hace plot(x,y). Si no, hist(x).
#' @param tipo Tipo de gráfico base R para `plot` ("p" puntos, "l" líneas, etc.).
#' @examples
#' grafica(1:10, (1:10)^2)
#' grafica(rnorm(100))
#' @export
grafica <- function(x, y = NULL, tipo = "p"){
  if (is.null(y)){
    graphics::hist(x, main = "Histograma", xlab = deparse(substitute(x)))
  } else {
    graphics::plot(x, y, type = tipo, xlab = deparse(substitute(x)), ylab = deparse(substitute(y)))
  }
}

#' Diagnóstico sencillo de residuos de un modelo lineal
#' @param modelo Objeto `lm`.
#' @examples
#' m <- lm(mpg ~ wt, data = mtcars); grafica_residuos(m)
#' @export
grafica_residuos <- function(modelo){
  stopifnot(inherits(modelo, "lm"))
  op <- graphics::par(mfrow = c(1,2))
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(stats::fitted(modelo), stats::residuals(modelo),
                 xlab = "Ajustados", ylab = "Residuos", main = "Residuos vs Ajustados")
  graphics::abline(h = 0, lty = 2)
  stats::qqnorm(stats::rstandard(modelo), main = "QQ de residuos estandarizados")
  stats::qqline(stats::rstandard(modelo))
}

#' Generador de datos de ejemplo (lineales con ruido)
#' @param n Número de observaciones.
#' @param seed Semilla opcional para reproducibilidad.
#' @return data.frame con columnas x, y.
#' @examples
#' df <- datos_ejemplo(20, seed = 123)
#' @export
datos_ejemplo <- function(n = 50, seed = NULL){
  stopifnot(is.numeric(n), length(n) == 1, n > 1)
  if (!is.null(seed)) set.seed(seed)
  x <- seq_len(n)
  y <- 3 + 1.5*x + stats::rnorm(n, sd = sd(x)/2)
  data.frame(x = x, y = y)
}

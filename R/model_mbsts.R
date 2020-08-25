#' Multivariate structural time series model definition
#'
#' @param y t x d matrix or data.frame of observations.
#' @param components Character vector specifying the components of the
#'   multivariate structural time series model. Possible values in c("trend",
#'   "slope", "seasonal", "cycle").
#' @param seas.period Length of the seasonal pattern.
#' @param cycle.period Length of the cycle pattern.
#'
#' @return An object of class 'SSModel'.
#' @export
#' @keywords internal
#'
#' @examples
#' # Example 1 : local level + seasonal
#' y <- cbind(seq(0.5,200,by=0.5)*0.1 + rnorm(400),
#'            seq(100.25,200,by=0.25)*0.05 + rnorm(400),
#'            rnorm(400, 5,1))
#' model.1 <- model(y, components = c("trend","seasonal"), seas.period = 7)
#'
#' # Example 2: local level  + cycle
#' t <- seq(from = 0,to = 4*pi, length.out=300)
#' y <- cbind(3*sin(2*t)+rnorm(300), 2*cos(2*t) + rnorm(300))
#' model.2 <- model(y, components = c("trend", "cycle"), cycle.period = 75)

model <- function(y, components, seas.period = NULL, cycle.period = NULL){

  y <- as.matrix(y)

  if("trend" %in% components & "slope" %in% components){

    mt<-paste("SSM","trend","(", "degree = 2", " , ", "Q = list(matrix(NA), matrix(NA))",")", sep="")
  } else if ("trend" %in% components & !("slope" %in% components)) {
    mt<-paste("SSM","trend","(", "degree = 1", " , ", "Q = matrix(NA)",")", sep="")
  } else { mt <- NULL }

  if("seasonal" %in% components){
    if(is.null(seas.period)){ stop("A seasonal model needs a seas.period") }
    ms <- paste("SSM","seasonal", "(", "period = ", seas.period, " , ", "Q = matrix(NA)", ")", sep="")
  } else { ms <- NULL }

  if("cycle" %in% components){
    if(is.null(cycle.period)){ stop("A cyclical model needs a cycle.period") }
    mc <- paste("SSM","cycle", "(", "period = ", cycle.period, " , ", "Q = matrix(NA)", ")", sep="")
  } else { mc <- NULL }

  formula <- as.formula(paste("y ~ ", paste(unlist(list(mt,ms,mc)), collapse = " + ")))

  m <- SSModel(formula)

  return(m)

}

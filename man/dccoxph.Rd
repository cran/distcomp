% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survivalstuff.R
\name{dccoxph}
\alias{dccoxph}
\alias{dccoxph.fit}
\title{Functions copied and modified from survival package}
\usage{
dccoxph(
  formula,
  data,
  weights,
  subset,
  na.action,
  init,
  control,
  ties = c("efron", "breslow", "exact"),
  singular.ok = TRUE,
  robust = FALSE,
  model = FALSE,
  x = FALSE,
  y = TRUE,
  tt,
  method = ties,
  ...
)

dccoxph.fit(x, y, strata, offset, init, control, weights, method, rownames)
}
\description{
The distcomp package makes use of code from the survival package,
with the permission of the original authors. This includes R code as well as C code.
That is, the underlying Cox model code is derived from that in the R survival package.
The original copyrights are retained for these files and the notices preserved.
However, these are for internal use and future implementations may change how we use them.
In order to avoid confusion and any name collision, the names of these functions have
been modified to include a prefix "dc".
}
\keyword{internal}

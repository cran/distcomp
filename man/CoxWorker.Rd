% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coxstuff.R
\name{CoxWorker}
\alias{CoxWorker}
\title{R6 class for object to use as a worker with \code{CoxMaster} master objects}
\description{
\code{CoxWorker} objects are worker objects at each data site of
a distributed Cox model computation
}
\seealso{
\code{CoxMaster} which goes hand-in-hand with this object
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{CoxWorker$new()}}
\item \href{#method-getP}{\code{CoxWorker$getP()}}
\item \href{#method-getStateful}{\code{CoxWorker$getStateful()}}
\item \href{#method-logLik}{\code{CoxWorker$logLik()}}
\item \href{#method-var}{\code{CoxWorker$var()}}
\item \href{#method-kosher}{\code{CoxWorker$kosher()}}
\item \href{#method-clone}{\code{CoxWorker$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\subsection{Method \code{new()}}{
Create a new \code{CoxWorker} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CoxWorker$new(defn, data, stateful = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{defn}}{the computation definition}

\item{\code{data}}{the local data}

\item{\code{stateful}}{a boolean flag indicating if state needs to be preserved between REST calls}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a new \code{CoxWorker} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getP"></a>}}
\subsection{Method \code{getP()}}{
Return the dimension of the parameter vector.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CoxWorker$getP(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{other args ignored}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
the dimension of the parameter vector
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getStateful"></a>}}
\subsection{Method \code{getStateful()}}{
Return the stateful status of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CoxWorker$getStateful()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
the stateful flag, \code{TRUE} or \code{FALSE}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-logLik"></a>}}
\subsection{Method \code{logLik()}}{
Return the partial log likelihood on local data for given \code{beta} parameter.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CoxWorker$logLik(beta, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{beta}}{the parameter vector}

\item{\code{...}}{further arguments, currently unused}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a named list with three components: \code{value} contains the value of the
log likelihood, \code{gradient} contains the score vector, and \code{hessian} contains
the estimated hessian matrix
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-var"></a>}}
\subsection{Method \code{var()}}{
Return the variance of estimate for given \code{beta} parameter on local data.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CoxWorker$var(beta, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{beta}}{the parameter vector}

\item{\code{...}}{further arguments, currently unused}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
variance vector
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-kosher"></a>}}
\subsection{Method \code{kosher()}}{
Check if inputs and state of object are sane. For future use
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CoxWorker$kosher()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
\code{TRUE} or \code{FALSE}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{CoxWorker$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}

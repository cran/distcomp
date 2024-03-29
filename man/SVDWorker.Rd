% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svdstuff.R
\name{SVDWorker}
\alias{SVDWorker}
\title{R6 class for a SVD worker object to use with master objects generated by \code{\link[=SVDMaster]{SVDMaster()}}}
\description{
\code{SVDWorker} objects are worker objects at each site of a distributed SVD model computation
}
\seealso{
\code{\link[=SVDMaster]{SVDMaster()}} which goes hand-in-hand with this object
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-SVDWorker-new}{\code{SVDWorker$new()}}
\item \href{#method-SVDWorker-reset}{\code{SVDWorker$reset()}}
\item \href{#method-SVDWorker-dimX}{\code{SVDWorker$dimX()}}
\item \href{#method-SVDWorker-updateV}{\code{SVDWorker$updateV()}}
\item \href{#method-SVDWorker-updateU}{\code{SVDWorker$updateU()}}
\item \href{#method-SVDWorker-normU}{\code{SVDWorker$normU()}}
\item \href{#method-SVDWorker-fixU}{\code{SVDWorker$fixU()}}
\item \href{#method-SVDWorker-getN}{\code{SVDWorker$getN()}}
\item \href{#method-SVDWorker-getP}{\code{SVDWorker$getP()}}
\item \href{#method-SVDWorker-getStateful}{\code{SVDWorker$getStateful()}}
\item \href{#method-SVDWorker-kosher}{\code{SVDWorker$kosher()}}
\item \href{#method-SVDWorker-clone}{\code{SVDWorker$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SVDWorker-new"></a>}}
\if{latex}{\out{\hypertarget{method-SVDWorker-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{SVDWorker} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SVDWorker$new(defn, data, stateful = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{defn}}{the computation definition}

\item{\code{data}}{the local \code{x} matrix}

\item{\code{stateful}}{a boolean flag indicating if state needs to be preserved between REST calls, \code{TRUE} by default}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a new \code{SVDWorker} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SVDWorker-reset"></a>}}
\if{latex}{\out{\hypertarget{method-SVDWorker-reset}{}}}
\subsection{Method \code{reset()}}{
Reset the computation state by initializing work matrix and set up starting values for iterating
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SVDWorker$reset()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SVDWorker-dimX"></a>}}
\if{latex}{\out{\hypertarget{method-SVDWorker-dimX}{}}}
\subsection{Method \code{dimX()}}{
Return the dimensions of the matrix
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SVDWorker$dimX(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{other args ignored}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
the dimension of the matrix
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SVDWorker-updateV"></a>}}
\if{latex}{\out{\hypertarget{method-SVDWorker-updateV}{}}}
\subsection{Method \code{updateV()}}{
Return an updated value for the \code{V} vector, normalized by \code{arg}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SVDWorker$updateV(arg, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{arg}}{the normalizing value}

\item{\code{...}}{other args ignored}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
updated \code{V}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SVDWorker-updateU"></a>}}
\if{latex}{\out{\hypertarget{method-SVDWorker-updateU}{}}}
\subsection{Method \code{updateU()}}{
Update \code{U} and return the updated norm of \code{U}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SVDWorker$updateU(arg, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{arg}}{the initial value}

\item{\code{...}}{other args ignored}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
updated norm of \code{U}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SVDWorker-normU"></a>}}
\if{latex}{\out{\hypertarget{method-SVDWorker-normU}{}}}
\subsection{Method \code{normU()}}{
Normalize \code{U} vector
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SVDWorker$normU(arg, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{arg}}{the normalizing value}

\item{\code{...}}{other args ignored}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
\code{TRUE} invisibly
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SVDWorker-fixU"></a>}}
\if{latex}{\out{\hypertarget{method-SVDWorker-fixU}{}}}
\subsection{Method \code{fixU()}}{
Construct residual matrix using \code{arg}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SVDWorker$fixU(arg, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{arg}}{the value to use for residualizing}

\item{\code{...}}{other args ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SVDWorker-getN"></a>}}
\if{latex}{\out{\hypertarget{method-SVDWorker-getN}{}}}
\subsection{Method \code{getN()}}{
Getthe number of rows of \code{x} matrix
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SVDWorker$getN()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
the number of rows of \code{x} matrix
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SVDWorker-getP"></a>}}
\if{latex}{\out{\hypertarget{method-SVDWorker-getP}{}}}
\subsection{Method \code{getP()}}{
Getthe number of columnsof \code{x} matrix
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SVDWorker$getP()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
the number of columns of \code{x} matrix
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SVDWorker-getStateful"></a>}}
\if{latex}{\out{\hypertarget{method-SVDWorker-getStateful}{}}}
\subsection{Method \code{getStateful()}}{
Return the stateful status of the object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SVDWorker$getStateful()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
the stateful flag, \code{TRUE} or \code{FALSE}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SVDWorker-kosher"></a>}}
\if{latex}{\out{\hypertarget{method-SVDWorker-kosher}{}}}
\subsection{Method \code{kosher()}}{
Check if inputs and state of object are sane. For future use
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SVDWorker$kosher()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
\code{TRUE} or \code{FALSE}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-SVDWorker-clone"></a>}}
\if{latex}{\out{\hypertarget{method-SVDWorker-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SVDWorker$clone(deep = FALSE)}\if{html}{\out{</div>}}
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

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UndClassFunctions.R
\name{summary.UndPCDs}
\alias{summary.UndPCDs}
\title{Return a summary of an \code{UndPCDs} \code{object}}
\usage{
\method{summary}{UndPCDs}(object, ...)
}
\arguments{
\item{object}{An \code{object} of class \code{UndPCDs}.}

\item{\dots}{Additional parameters for \code{summary}.}
}
\value{
The \code{call} of the \code{object} of class \code{"UndPCDs"},
the \code{type} (i.e. the description) of the underlying or reflexivity graphs
of the proximity catch digraph (PCD),
some of the partition
(i.e. intervalization in the 1D case and triangulation
in the 2D case) points
(i.e., vertices of the intervals or the triangles),
parameters of the underlying or reflexivity graph of the PCD,
and various quantities
(number of vertices,
number of edges and edge density of the underlying
or reflexivity graphs of the PCDs,
number of vertices for the partition
and number of partition cells
(i.e., intervals or triangles)).
}
\description{
Returns the below information about the \code{object}:

\code{call} of the function defining the \code{object},
the \code{type} (i.e. the description) of the underlying or reflexivity graph
of the proximity catch digraph (PCD),
some of the partition
(i.e. intervalization in the 1D case and triangulation
in the 2D case) points
(i.e., vertices of the intervals or the triangles),
parameter of the underlying or reflexivity graphs of the PCD,
and various quantities
(number of vertices,
number of edges and edge density of the underlying
or reflexivity graphs of the PCDs,
number of vertices for the partition and number of partition cells
(i.e., intervals or triangles)).
}
\examples{
#\donttest{
nx<-20; ny<-5;
set.seed(1)
Xp<-cbind(runif(nx,0,1),runif(nx,0,1))
Yp<-cbind(runif(ny,0,.25),runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
M<-c(1,1,1)  #try also M<-c(1,2,3)
r<-1.5
Edges<-edgesPE(Xp,Yp,r,M)
Edges
summary(Edges)
#}

}
\seealso{
\code{\link{print.UndPCDs}}, \code{\link{print.summary.UndPCDs}},
and \code{\link{plot.UndPCDs}}
}

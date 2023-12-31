% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UndClassFunctions.R
\name{summary.NumEdges}
\alias{summary.NumEdges}
\title{Return a summary of a \code{NumEdges} \code{object}}
\usage{
\method{summary}{NumEdges}(object, ...)
}
\arguments{
\item{object}{An \code{object} of class \code{NumEdges}.}

\item{\dots}{Additional parameters for \code{summary}.}
}
\value{
The \code{call} of the \code{object} of class \code{"NumEdges"},
the \code{desc} of the output,
the type of the graph as "underlying" or "reflexivity",
total number of edges in the underlying or reflexivity graph.
Moreover, in the one Delaunay cell case, the function also provides
vertices of Delaunay cell, and
indices of target points in the Delaunay cell;
and in the multiple Delaunay cell case, it also provides
number of edges for the induced subgraphs for points in the Delaunay cells,
vertices of Delaunay cells or indices of points that form the the Delaunay cells,
indices of target points in the convex hull of nontarget points,
indices of Delaunay cells in which points reside,
and area or length of the the Delaunay cells.
}
\description{
Returns the below information about the \code{object}:

\code{call} of the function defining the \code{object},
the description of the output, \code{desc}, and
type of the graph as "underlying" or "reflexivity",
number of edges in the underlying or reflexivity graph of
the proximity catch digraph (PCD) and
related quantities in the induced subgraphs for points in the Delaunay cells.
In the one Delaunay cell case, the function provides
the total number of edges in the underlying or reflexivity graph,
vertices of Delaunay cell, and
indices of target points in the Delaunay cell.

In the multiple Delaunay cell case, the function provides
total number of edges in the underlying or reflexivity graph,
number of edges for the induced subgraphs for points in the Delaunay cells,
vertices of Delaunay cells or indices of points that form the the Delaunay cells,
indices of target points in the convex hull of nontarget points,
indices of Delaunay cells in which points reside,
and area or length of the the Delaunay cells.
}
\examples{
#\donttest{
nx<-15; ny<-5;
set.seed(1)
Xp<-cbind(runif(nx),runif(nx))
Yp<-cbind(runif(ny,0,.25),
runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))

M<-c(1,1,1)  #try also M<-c(1,2,3)

Nedges = num.edgesAS(Xp,Yp,M)
Nedges
summary(Nedges)
#}

}
\seealso{
\code{\link{print.NumEdges}}, \code{\link{print.summary.NumEdges}},
and \code{\link{plot.NumEdges}}
}

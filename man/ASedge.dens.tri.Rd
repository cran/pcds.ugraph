% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ArcSliceUGFunctions.R
\name{ASedge.dens.tri}
\alias{ASedge.dens.tri}
\title{Edge density of the underlying or reflexivity graph of
Arc Slice Proximity Catch Digraphs (AS-PCDs) -
one triangle case}
\usage{
ASedge.dens.tri(
  Xp,
  tri,
  M = "CC",
  ugraph = c("underlying", "reflexivity"),
  in.tri.only = FALSE
)
}
\arguments{
\item{Xp}{A set of 2D points
which constitute the vertices of the underlying or reflexivity graph of the AS-PCD.}

\item{tri}{A \eqn{3 \times 2} matrix with each row
representing a vertex of the triangle.}

\item{M}{The center of the triangle.
\code{"CC"} stands for circumcenter of the triangle \code{tri}
or a 2D point in Cartesian coordinates or
a 3D point in barycentric coordinates
which serves as a center in the interior of \code{tri};
default is \code{M="CC"}, i.e., the circumcenter of \code{tri}.}

\item{ugraph}{The type of the graph based on AS-PCDs,
\code{"underlying"} is for the underlying graph,
and \code{"reflexivity"} is for
the reflexivity graph (default is \code{"underlying"}).}

\item{in.tri.only}{A logical argument (default is \code{in.tri.only=FALSE})
for computing the edge density for only the points inside the triangle, \code{tri}.
That is,
if \code{in.tri.only=TRUE} edge density of the induced subgraph with the vertices
inside \code{tri} is computed, otherwise
otherwise edge density of the entire graph
(i.e., graph with all the vertices) is computed.}
}
\value{
Edge density of the underlying
or reflexivity graphs based on the AS-PCD
whose vertices are the 2D numerical data set, \code{Xp};
AS proximity regions are defined
with respect to the triangle \code{tri} and \code{M}-vertex regions.
}
\description{
Returns the edge density
of the underlying or reflexivity graph of
Arc Slice Proximity Catch Digraphs (AS-PCDs)
whose vertex set is the given 2D numerical data set, \code{Xp},
(some of its members are) in the triangle \code{tri}.

AS proximity regions are defined with respect to \code{tri}
and vertex regions are defined with the center \code{M="CC"}
for circumcenter of \code{tri};
or \eqn{M=(m_1,m_2)} in Cartesian coordinates
or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
interior of the triangle \code{tri};
default is \code{M="CC"}, i.e., circumcenter of \code{tri}.
For the number of edges,
loops are not allowed so edges are only possible for points inside \code{tri}
for this function.

\code{in.tri.only} is a logical argument (default is \code{FALSE})
for considering only the points
inside the triangle or all the points as the vertices of the digraph.
if \code{in.tri.only=TRUE}, edge density is computed only for
the points inside the triangle
(i.e., edge density of the subgraph of the underlying or reflexivity graph
induced by the vertices in the triangle is computed),
otherwise edge density of the entire graph
(i.e., graph with all the vertices) is computed.

See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph}).
}
\examples{
#\donttest{
A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
Tr<-rbind(A,B,C);
n<-10

set.seed(1)
Xp<-pcds::runif.tri(n,Tr)$g

M<-as.numeric(pcds::runif.tri(1,Tr)$g)

#For the underlying graph
(num.edgesAStri(Xp,Tr,M)$num.edges)/(n*(n-1)/2)
ASedge.dens.tri(Xp,Tr,M)
ASedge.dens.tri(Xp,Tr,M,in.tri.only = TRUE)

#For the reflexivity graph
(num.edgesAStri(Xp,Tr,M,ugraph="r")$num.edges)/(n*(n-1)/2)
ASedge.dens.tri(Xp,Tr,M,ugraph="r")
ASedge.dens.tri(Xp,Tr,M,in.tri.only = TRUE,ugraph="r")
#}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{PEedge.dens.tri}}, \code{\link{CSedge.dens.tri}},
and \code{\link[pcds]{ASarc.dens.tri}}
}
\author{
Elvan Ceyhan
}

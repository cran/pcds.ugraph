% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ArcSliceUGFunctions.R
\name{inci.mat.undAStri}
\alias{inci.mat.undAStri}
\title{Incidence matrix for the underlying or reflexivity graph of
Arc Slice Proximity Catch Digraphs (AS-PCDs) -
one triangle case}
\usage{
inci.mat.undAStri(Xp, tri, M = "CC", ugraph = c("underlying", "reflexivity"))
}
\arguments{
\item{Xp}{A set of 2D points
which constitute the vertices of the underlying
or reflexivity graph of the AS-PCD.}

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
}
\value{
Incidence matrix for the underlying or reflexivity graph
of the AS-PCD whose vertices are the 2D data set, \code{Xp}
in the triangle \code{tri} with vertex regions based on the center \code{M}
}
\description{
Returns the incidence matrix
for the underlying or reflexivity graph of the AS-PCD
whose vertices are the given 2D numerical data set, \code{Xp},
in the triangle \code{tri}\eqn{=T(v=1,v=2,v=3)}.

AS proximity regions are defined
with respect to the triangle \code{tri}\eqn{=T(v=1,v=2,v=3)} and
vertex regions are based on the center, \eqn{M=(m_1,m_2)} in Cartesian coordinates
or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
in the interior of the triangle \code{tri}
or based on circumcenter of  \code{tri};
default is \code{M="CC"}, i.e., circumcenter of \code{tri}.
Loops are allowed, so the diagonal entries are all equal to 1.

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
(IM<-inci.mat.undAStri(Xp,Tr,M))
pcds::dom.num.greedy(IM)
pcds::Idom.num.up.bnd(IM,3)

(IM<-inci.mat.undAStri(Xp,Tr,M,ugraph="r"))
pcds::dom.num.greedy(IM)
pcds::Idom.num.up.bnd(IM,3)
#}

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{inci.mat.undAS}}, \code{\link{inci.mat.undPEtri}},
\code{\link{inci.mat.undCStri}}, and \code{\link[pcds]{inci.matAStri}}
}
\author{
Elvan Ceyhan
}

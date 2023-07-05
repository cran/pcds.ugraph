#ArcSliceUGFunctions.R;
#Functions for Underlying or Reflexivity Graphs of AS-PCD in R^2
#################################################################

#' @import interp
#' @importFrom Rdpack reprompt

#' @title The indicator for the presence of an edge from a point to another
#' for the underlying or reflexivity graph of
#' Arc Slice Proximity Catch Digraphs (AS-PCDs) -
#' standard basic triangle case
#'
#' @description Returns \eqn{I(}\code{p1p2} is an edge
#' in the underlying or reflexivity graph of AS-PCDs \eqn{)}
#' for points \code{p1} and \code{p2} in the standard basic triangle.
#'
#' More specifically, when the argument \code{ugraph="underlying"},
#' it returns the edge indicator for the AS-PCD underlying graph,
#' that is, returns 1 if \code{p2} is
#' in \eqn{N_{AS}(p1)} **or** \code{p1} is in \eqn{N_{AS}(p2)},
#' returns 0 otherwise.
#' On the other hand,
#' when \code{ugraph="reflexivity"}, it returns
#' the edge indicator for the AS-PCD reflexivity graph,
#' that is, returns 1 if \code{p2} is
#' in \eqn{N_{AS}(p1)} **and** \code{p1} is in \eqn{N_{AS}(p2)},
#' returns 0 otherwise.
#'
#' AS proximity region is constructed in the standard basic triangle
#' \eqn{T_b=T((0,0),(1,0),(c_1,c_2))}
#' where \eqn{c_1} is in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \leq 1}.
#'
#' Vertex regions are based on the center, \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the standard basic triangle \eqn{T_b}
#' or based on circumcenter of \eqn{T_b};
#' default is \code{M="CC"}, i.e., circumcenter of \eqn{T_b}.
#'
#' If \code{p1} and \code{p2} are distinct
#' and either of them are outside \eqn{T_b}, it returns 0,
#' but if they are identical,
#' then it returns 1 regardless of their locations (i.e., it allows loops).
#'
#' Any given triangle can be mapped to the standard basic triangle
#' by a combination of rigid body motions
#' (i.e., translation, rotation and reflection) and scaling,
#' preserving uniformity of the points in the original triangle.
#' Hence, standard basic triangle is useful for simulation
#' studies under the uniformity hypothesis.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param p1 A 2D point whose AS proximity region is constructed.
#' @param p2 A 2D point. The function determines
#' whether there is an edge from \code{p1} to \code{p1} or not
#' in the underlying or reflexivity graph of AS-PCDs.
#' @param c1,c2 Positive real numbers
#' which constitute the vertex of the standard basic triangle
#' adjacent to the shorter edges;
#' \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M The center of the triangle. \code{"CC"} stands for
#' circumcenter of the triangle \eqn{T_b}
#' or a 2D point in Cartesian coordinates or
#' a 3D point in barycentric coordinates
#' which serves as a center in the interior of \eqn{T_b};
#' default is \code{M="CC"}, i.e., the circumcenter of \eqn{T_b}.
#' @param ugraph The type of the graph based on AS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return Returns 1 if there is an edge between points \code{p1} and \code{p2}
#' in the underlying or reflexivity graph of AS-PCDs
#' in the standard basic triangle, and 0 otherwise.
#'
#' @seealso \code{\link{IedgeAStri}}, \code{\link{IedgeCSbasic.tri}},
#' \code{\link{IedgePEbasic.tri}} and \code{\link[pcds]{IarcASbasic.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #\donttest{
#' c1<-.4; c2<-.6
#' A<-c(0,0); B<-c(1,0); C<-c(c1,c2);
#' Tb<-rbind(A,B,C);
#'
#' M<-as.numeric(pcds::runif.basic.tri(1,c1,c2)$g)
#' set.seed(4)
#' P1<-as.numeric(pcds::runif.basic.tri(1,c1,c2)$g)
#' P2<-as.numeric(pcds::runif.basic.tri(1,c1,c2)$g)
#' IedgeASbasic.tri(P1,P2,c1,c2,M)
#' IedgeASbasic.tri(P1,P2,c1,c2,M,ugraph = "reflexivity")
#'
#' P1<-c(.4,.2)
#' P2<-c(.5,.26)
#' IedgeASbasic.tri(P1,P2,c1,c2,M)
#' IedgeASbasic.tri(P1,P2,c1,c2,M,ugraph="r")
#' #}
#'
#' @export IedgeASbasic.tri
IedgeASbasic.tri <- function(p1,p2,c1,c2,M="CC",
                             ugraph=c("underlying", "reflexivity"))
{
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  arc12 = pcds::IarcASbasic.tri(p1,p2,c1,c2,M)
  arc21 = pcds::IarcASbasic.tri(p2,p1,c1,c2,M)

  edge = ifelse(ugraph == "underlying",
                max(arc12,arc21),
                arc12*arc21)
  edge
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an edge from a point to another
#' for the underlying or reflexivity graph of
#' Arc Slice Proximity Catch Digraphs (AS-PCDs) -
#' one triangle case
#'
#' @description Returns \eqn{I(}\code{p1p2} is an edge
#' in the underlying or reflexivity graph of AS-PCDs \eqn{)}
#' for points \code{p1} and \code{p2} in a given triangle.
#'
#' More specifically, when the argument \code{ugraph="underlying"},
#' it returns the edge indicator for the AS-PCD underlying graph,
#' that is, returns 1 if \code{p2} is
#' in \eqn{N_{AS}(p1)} **or** \code{p1} is in \eqn{N_{AS}(p2)},
#' returns 0 otherwise.
#' On the other hand,
#' when \code{ugraph="reflexivity"}, it returns
#' the edge indicator for the AS-PCD reflexivity graph,
#' that is, returns 1 if \code{p2} is
#' in \eqn{N_{AS}(p1)} **and** \code{p1} is in \eqn{N_{AS}(p2)},
#' returns 0 otherwise.
#'
#' In both cases AS proximity region is constructed
#' with respect to the triangle \code{tri} and
#' vertex regions are based on the center, \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri}
#' or based on circumcenter of  \code{tri};
#' default is \code{M="CC"}, i.e., circumcenter of \code{tri}.
#'
#' If \code{p1} and \code{p2} are distinct
#' and either of them are outside \code{tri}, it returns 0,
#' but if they are identical,
#' then it returns 1 regardless of their locations
#' (i.e., it allows loops).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param p1 A 2D point whose AS proximity region is constructed.
#' @param p2 A 2D point. The function determines
#' whether there is an edge from \code{p1} to \code{p1} or not
#' in the underlying or reflexivity graph of AS-PCDs.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param M The center of the triangle. \code{"CC"} stands for
#' circumcenter of the triangle \code{tri}
#' or a 2D point in Cartesian coordinates or
#' a 3D point in barycentric coordinates
#' which serves as a center in the interior of \code{tri};
#' default is \code{M="CC"}, i.e., the circumcenter of \code{tri}.
#' @param ugraph The type of the graph based on AS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return Returns 1 if there is an edge between points \code{p1} and \code{p2}
#' in the underlying or reflexivity graph of AS-PCDs
#' in a given triangle \code{tri}, and 0 otherwise.
#'
#' @seealso \code{\link{IedgeASbasic.tri}}, \code{\link{IedgePEtri}},
#' \code{\link{IedgeCStri}} and \code{\link[pcds]{IarcAStri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @importFrom graphics abline points polygon segments text
#' @importFrom stats pnorm qnorm runif
#'
#' @examples
#' #\donttest{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#' M<-as.numeric(pcds::runif.tri(1,Tr)$g)
#'
#' n<-3
#' set.seed(1)
#' Xp<-pcds::runif.tri(n,Tr)$g
#'
#' IedgeAStri(Xp[1,],Xp[3,],Tr,M)
#' IedgeAStri(Xp[1,],Xp[3,],Tr,M,ugraph = "reflexivity")
#'
#' set.seed(1)
#' P1<-as.numeric(pcds::runif.tri(1,Tr)$g)
#' P2<-as.numeric(pcds::runif.tri(1,Tr)$g)
#' IedgeAStri(P1,P2,Tr,M)
#' IedgeAStri(P1,P2,Tr,M,ugraph="r")
#' #}
#'
#' @export IedgeAStri
IedgeAStri <- function(p1,p2,tri,M="CC",ugraph=c("underlying", "reflexivity"))
{
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  arc12 = pcds::IarcAStri(p1,p2,tri,M)
  arc21 = pcds::IarcAStri(p2,p1,tri,M)

  edge = ifelse(ugraph == "underlying",
                max(arc12,arc21),
                arc12*arc21)
  edge
} #end of the function
#'

#################################################################

#' @title Number of edges of the underlying or reflexivity graph of
#' Arc Slice Proximity Catch Digraphs (AS-PCDs) -
#' one triangle case
#'
#' @description
#' An object of class \code{"NumEdges"}.
#' Returns the number of edges of
#' the underlying or reflexivity graph of
#' Arc Slice Proximity Catch Digraphs (AS-PCDs)
#' whose vertices are the
#' given 2D numerical data set, \code{Xp}
#' in a given triangle \code{tri}.
#' It also provides number of vertices
#' (i.e., number of data points inside the triangle)
#' and indices of the data points that reside in the triangle.
#'
#' AS proximity regions are defined
#' with respect to the triangle \code{tri} and vertex regions are
#' based on the center, \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri}
#' or based on circumcenter of  \code{tri};
#' default is \code{M="CC"}, i.e., circumcenter of \code{tri}.
#' For the number of edges, loops are not allowed,
#' so edges are only possible for points inside the triangle, \code{tri}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graph of the AS-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' @param ugraph The type of the graph based on AS-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of edges
#' and quantities related to the triangle}
#' \item{und.graph}{Type of the graph as "Underlying"
#' or "Reflexivity" for the AS-PCD}
#' \item{num.edges}{Number of edges of the underlying
#' or reflexivity graphs of the AS-PCD
#' with vertices in the given triangle \code{tri}}
#' \item{num.in.tri}{Number of \code{Xp} points in the triangle, \code{tri}}
#' \item{ind.in.tri}{The vector of indices of the \code{Xp} points
#' that reside in the triangle}
#' \item{tess.points}{Tessellation points,
#' i.e., points on which the tessellation of the study region is performed,
#' here, tessellation is the support triangle.}
#' \item{vertices}{Vertices of the underlying or reflexivity graph, \code{Xp}.}
#'
#' @seealso \code{\link{num.edgesAS}}, \code{\link{num.edgesPEtri}},
#' \code{\link{num.edgesCStri}}, and \code{\link[pcds]{num.arcsAStri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #\donttest{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#'
#' n<-10
#' set.seed(1)
#' Xp<-pcds::runif.tri(n,Tr)$g
#'
#' M<-as.numeric(pcds::runif.tri(1,Tr)$g)
#'
#' Nedges = num.edgesAStri(Xp,Tr,M)
#' Nedges
#' summary(Nedges)
#' plot(Nedges)
#' #}
#'
#' @export num.edgesAStri
num.edgesAStri <- function(Xp,tri,M="CC",
                           ugraph=c("underlying", "reflexivity"))
{
  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (pcds::is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('The triangle is degenerate!')}

  if (!(pcds::is.point(M) || pcds::is.point(M,3) || identical(M,"CC")))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates or the circumcenter "CC" ')}

  CC = pcds::circumcenter.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (pcds::dimension(M)==3)
  {M<-pcds::bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) ||
        pcds::in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('M is not the circumcenter or not a center in the interior of the triangle')}

  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  n<-nrow(Xp)
  tot.edges<-edges.in.tri<-0
  ind.in.tri = c()
  if (n<=0)
  {
    tot.edges<-edges.in.tri<-0
  } else
  {
    for (i in 1:n)
    {
      if (pcds::in.triangle(Xp[i,],tri,boundary=TRUE)$in.tri)
      { ind.in.tri = c(ind.in.tri,i)
      for (k in (i:n)[-1]) #to avoid loops
      {
        edges.in.tri<-edges.in.tri+IedgeAStri(Xp[i,],Xp[k,],tri,M,ugraph=ugraph)
      }
      }

      for (j in (i:n)[-1]) #to avoid loops
      {
        tot.edges<-tot.edges+IedgeAStri(Xp[i,],Xp[j,],tri,M,ugraph=ugraph)
      }
    }
  }

  NinTri = length(ind.in.tri)

  und.graph = ifelse(ugraph=="underlying",
                     "Underlying",
                     "Reflexivity")
  desc<-paste("Number of Edges of the ",und.graph,
              " Graphs of the AS-PCD and the Related Quantities with vertices Xp in One Triangle",sep="")

  res<-list(desc=desc, #description of the output
            und.graph = und.graph, #"Underlying" or "Reflexivity"
            num.edges=tot.edges,
            #total number of edges in the entire underlying or reflexivity graph
            tri.num.edges=edges.in.tri, #vector of number of edges for the triangle
            num.in.tri=NinTri, # number of Xp points in CH of Yp points
            ind.in.tri=ind.in.tri, #indices of data points inside the triangle
            tess.points=tri, #tessellation points
            vertices=Xp #vertices of the underlying or reflexivity graph
  )

  class(res) <- "NumEdges"
  res$call <-match.call()

  res
} #end of the function
#'

################################################################

#' @title Edge density of the underlying or reflexivity graph of
#' Arc Slice Proximity Catch Digraphs (AS-PCDs) -
#' one triangle case
#'
#' @description Returns the edge density
#' of the underlying or reflexivity graph of
#' Arc Slice Proximity Catch Digraphs (AS-PCDs)
#' whose vertex set is the given 2D numerical data set, \code{Xp},
#' (some of its members are) in the triangle \code{tri}.
#'
#' AS proximity regions are defined with respect to \code{tri}
#' and vertex regions are defined with the center \code{M="CC"}
#' for circumcenter of \code{tri};
#' or \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of the triangle \code{tri};
#' default is \code{M="CC"}, i.e., circumcenter of \code{tri}.
#' For the number of edges,
#' loops are not allowed so edges are only possible for points inside \code{tri}
#' for this function.
#'
#' \code{in.tri.only} is a logical argument (default is \code{FALSE})
#' for considering only the points
#' inside the triangle or all the points as the vertices of the digraph.
#' if \code{in.tri.only=TRUE}, edge density is computed only for
#' the points inside the triangle
#' (i.e., edge density of the subgraph of the underlying or reflexivity graph
#' induced by the vertices in the triangle is computed),
#' otherwise edge density of the entire graph
#' (i.e., graph with all the vertices) is computed.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying or reflexivity graph of the AS-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param M The center of the triangle.
#' \code{"CC"} stands for circumcenter of the triangle \code{tri}
#' or a 2D point in Cartesian coordinates or
#' a 3D point in barycentric coordinates
#' which serves as a center in the interior of \code{tri};
#' default is \code{M="CC"}, i.e., the circumcenter of \code{tri}.
#' @param ugraph The type of the graph based on AS-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#' @param in.tri.only A logical argument (default is \code{in.tri.only=FALSE})
#' for computing the edge density for only the points inside the triangle, \code{tri}.
#' That is,
#' if \code{in.tri.only=TRUE} edge density of the induced subgraph with the vertices
#' inside \code{tri} is computed, otherwise
#' otherwise edge density of the entire graph
#' (i.e., graph with all the vertices) is computed.
#'
#' @return Edge density of the underlying
#' or reflexivity graphs based on the AS-PCD
#' whose vertices are the 2D numerical data set, \code{Xp};
#' AS proximity regions are defined
#' with respect to the triangle \code{tri} and \code{M}-vertex regions.
#'
#' @seealso \code{\link{PEedge.dens.tri}}, \code{\link{CSedge.dens.tri}},
#' and \code{\link[pcds]{ASarc.dens.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #\donttest{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#' n<-10
#'
#' set.seed(1)
#' Xp<-pcds::runif.tri(n,Tr)$g
#'
#' M<-as.numeric(pcds::runif.tri(1,Tr)$g)
#'
#' #For the underlying graph
#' (num.edgesAStri(Xp,Tr,M)$num.edges)/(n*(n-1)/2)
#' ASedge.dens.tri(Xp,Tr,M)
#' ASedge.dens.tri(Xp,Tr,M,in.tri.only = TRUE)
#'
#' #For the reflexivity graph
#' (num.edgesAStri(Xp,Tr,M,ugraph="r")$num.edges)/(n*(n-1)/2)
#' ASedge.dens.tri(Xp,Tr,M,ugraph="r")
#' ASedge.dens.tri(Xp,Tr,M,in.tri.only = TRUE,ugraph="r")
#' #}
#'
#' @export ASedge.dens.tri
ASedge.dens.tri <- function(Xp,tri,M="CC",
                            ugraph=c("underlying", "reflexivity"),in.tri.only=FALSE)
{
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  if (pcds::is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  nx<-nrow(Xp)
  if (nx<=1)
  {stop('The graph is void or has only one vertex!
    So, there are not enough Xp points to compute the edge density!')}

  nedges<-num.edgesAStri(Xp,tri,M,ugraph)$num.edges

  if (in.tri.only==TRUE)
  {
    ind.it<-c() #index of in triangle points
    for (i in 1:nx)
    {
      ind.it<-c(ind.it,pcds::in.triangle(Xp[i,],tri)$in.tri)
    }
    dat.it<-Xp[ind.it,] #Xp points inside the triangle
    NinTri<-nrow(dat.it)
    if (NinTri<=1)
    {stop('Induced subgraph for points in the triangle is void or has only one vertex!
    So, there are not enough Xp points in the triangle, tri, to compute the (corrected) edge density!')}
    n<-NinTri
  } else
  {
    n<-nx
  }
  rho<-nedges/(n*(n-1)/2)
  rho
} #end of the function
#'

#################################################################

#' @title Number of edges of
#' the underlying or reflexivity graph of
#' Arc Slice Proximity Catch Digraphs (AS-PCDs) -
#' multiple triangle case
#'
#' @description
#' An object of class \code{"NumEdges"}.
#' Returns the number of edges of
#' the underlying or reflexivity graph of
#' Arc Slice Proximity Catch Digraph (AS-PCD)
#' and various other quantities and vectors such as
#' the vector of number of vertices (i.e., number of data points)
#' in the Delaunay triangles,
#' number of data points in the convex hull of \code{Yp} points,
#' indices of the Delaunay triangles for the data points, etc.
#'
#' AS proximity regions are defined with respect to the
#' Delaunay triangles based on \code{Yp} points
#' and vertex regions in each triangle are based on the center \code{M="CC"}
#' for circumcenter of each Delaunay triangle
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of each Delaunay triangle;
#' default is \code{M="CC"}, i.e., circumcenter of each triangle.
#' Each Delaunay triangle is first converted to
#' a (nonscaled) basic triangle so that \code{M} will be the same
#' type of center for each Delaunay triangle
#' (this conversion is not necessary when \code{M} is \eqn{CM}).
#'
#' Convex hull of \code{Yp} is partitioned
#' by the Delaunay triangles based on \code{Yp} points
#' (i.e., multiple triangles are the set of these Delaunay triangles
#' whose union constitutes the
#' convex hull of \code{Yp} points).
#' For the number of edges,
#' loops are not allowed so edges are only possible
#' for points inside the convex hull of \code{Yp} points.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph})
#' for more on AS-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds.ugraph})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graph of the AS-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param M The center of the triangle.
#' \code{"CC"} stands for circumcenter of each Delaunay triangle
#' or 3D point in barycentric
#' coordinates which serves as a center
#' in the interior of each Delaunay triangle;
#' default is \code{M="CC"}, i.e., the circumcenter of each triangle.
#' @param ugraph The type of the graph based on AS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of edges
#' and related quantities for the induced subgraphs of the underlying
#' or reflexivity graphs (of AS-PCD) in the Delaunay triangles}
#' \item{und.graph}{Type of the graph as "Underlying" or "Reflexivity" for the AS-PCD}
#' \item{num.edges}{Total number of edges in all triangles,
#' i.e., the number of edges for the entire underlying
#' or reflexivity graphs of the AS-PCD}
#' \item{num.in.conv.hull}{Number of \code{Xp} points
#' in the convex hull of \code{Yp} points}
#' \item{num.in.tris}{The vector of number of \code{Xp} points
#' in the Delaunay triangles based on \code{Yp} points}
#' \item{weight.vec}{The \code{vector} of the areas of
#' Delaunay triangles based on \code{Yp} points}
#' \item{tri.num.edges}{The \code{vector} of the number of edges
#' of the components of the AS-PCD in the
#' Delaunay triangles based on \code{Yp} points}
#' \item{del.tri.ind}{A matrix of indices of vertices of
#' the Delaunay triangles based on \code{Yp} points,
#' each column corresponds to the vector of
#' indices of the vertices of one triangle.}
#' \item{data.tri.ind}{A \code{vector} of indices of vertices of
#' the Delaunay triangles in which data points reside,
#' i.e., column number of \code{del.tri.ind} for each \code{Xp} point.}
#' \item{tess.points}{Points on which the tessellation of the study region is performed,
#' here, tessellation is the Delaunay triangulation based on \code{Yp} points.}
#' \item{vertices}{Vertices of the underlying or reflexivity graph, \code{Xp}.}
#'
#' @seealso \code{\link{num.edgesAStri}}, \code{\link{num.edgesPE}},
#' \code{\link{num.edgesCS}}, and \code{\link[pcds]{num.arcsAS}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #\donttest{
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-15; ny<-5;
#'
#' set.seed(1)
#' Xp<-cbind(runif(nx),runif(nx))
#' Yp<-cbind(runif(ny,0,.25),
#' runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#'
#' pcds::plotDelaunay.tri(Xp,Yp,xlab="",ylab="")
#'
#' M<-c(1,1,1)
#'
#' Nedges = num.edgesAS(Xp,Yp,M)
#' Nedges
#' summary(Nedges)
#' plot(Nedges)
#' #}
#'
#' @export num.edgesAS
num.edgesAS <- function(Xp,Yp,M="CC",ugraph=c("underlying", "reflexivity"))
{
  if (!is.numeric(as.matrix(Xp)) || !is.numeric(as.matrix(Yp)))
  {stop('Xp and Yp must be numeric')}

  if (pcds::is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  Yp<-as.matrix(Yp)
  if (ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('Yp must be of dimension kx2 with k>=3')}

  if ((!pcds::is.point(M,3) && M!="CC") || !all(M>0))
  {stop('M must be a numeric 3D point with positive barycentric coordinates or
        "CC" for circumcenter')}

  nx<-nrow(Xp); ny<-nrow(Yp)

  #Delaunay triangulation of Yp points
  Ytrimesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
  Ytri<-matrix(interp::triangles(Ytrimesh)[,1:3],ncol=3);
  #the indices of the vertices of the Delaunay triangles (row-wise)
  ndt<-nrow(Ytri)  #number of Delaunay triangles

  inCH<-interp::in.convex.hull(Ytrimesh,Xp[,1],Xp[,2],strict=FALSE)
  NinCH<-sum(inCH)  #number of points in the convex hull

  Wvec=vector()
  for (i in 1:ndt)
  {
    ifelse(ndt==1,
           Tri<-Yp[Ytri,],
           Tri<-Yp[Ytri[i,],]) #vertices of ith triangle
    Wvec<-c(Wvec,pcds::area.polygon(Tri))
  }

  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  und.graph = ifelse(ugraph=="underlying",
                     "Underlying",
                     "Reflexivity")

  if (ny==3)
  {# tri<-pcds::as.basic.tri(Yp)$tri
 # NumEdges = num.edgesAStri(Xp,tri,M,ugraph) #for the entire digraph
  NinTri<-NinCH #NumEdges$num.in.tri #number of points in the triangle

    if (NinTri==0)
    {Tot.Edges<-0;
    ni.vec<-edges<-rep(0,ndt)
    data.tri.ind = ind.in.CH =  NA
    } else
    {
      Xdt<-matrix(Xp[inCH,],ncol=2)
      tri<-pcds::as.basic.tri(Yp)$tri
      #convert the triangle Yp into an nonscaled basic triangle, see as.basic.tri help page
      NumEdges = num.edgesAStri(Xdt,tri,M,ugraph) #for the vertices inside the triangle
      # Wvec<-pcds::area.polygon(tri)
      Tot.Edges <- NumEdges$num.edges #Tot.Edges<-num.edgesAStri(Xdt,tri,M,ugraph)$num.edges
      #number of edges in the triangle Yp
      ni.vec = NumEdges$num.in.tri
      Tri.Ind = NumEdges$ind.in.tri
      #returns 1's if the points Xp[i,]'s are inside triangle based on Yp, NA otherwise
      data.tri.ind = rep(NA,nx)
      data.tri.ind[Tri.Ind] = 1
      edges = NumEdges$num.edges
      ind.in.CH = which(inCH) #which(!is.na(Tri.Ind))
    }

  Tot.Edges = Tot.Edges + sum(duplicated(Xp[!inCH,]))

  desc<-paste("Number of Edges of the ",und.graph,
              " Graphs of the AS-PCD with vertices Xp and the Related Quantities for the Induced Subdigraph for the Points in the Delaunay Triangle",sep="")

    res<-list(desc=desc, #description of the output
              und.graph = und.graph, #"Underlying" or "Reflexivity"
              num.edges=Tot.Edges,
              tri.num.edges=edges,
              num.in.conv.hull=NinTri,
              ind.in.conv.hull= ind.in.CH, #indices of Xp points in the triangle
              num.in.tris=ni.vec,
              weight.vec=Wvec,
              del.tri.ind=t(Ytri),
              data.tri.ind=data.tri.ind,
              tess.points=Yp, #tessellation points
              vertices=Xp #vertices of the digraph
    )

  } else
  {
    if (NinCH==0)
    { Tot.Edges<-0;
    ni.vec<-edges<-rep(0,ndt)
    data.tri.ind = ind.in.CH =  NA
    } else
    {
      Tri.Ind<-pcds::indices.delaunay.tri(Xp,Yp,Ytrimesh)
      #indices of triangles in which the points in the data fall
      ind.in.CH = which(!is.na(Tri.Ind))

      #calculation of the total number of edges
      ni.vec<-edges<-vector()
      data.tri.ind = rep(NA,nx)
      for (i in 1:ndt)
      {
        dt.ind = which(Tri.Ind==i)
        #which indices of data points residing in ith Delaunay triangle
        Xpi<-Xp[dt.ind,] #points in ith Delaunay triangle
        data.tri.ind[dt.ind] = i
        #assigning the index of the Delaunay triangle that contains the data point
        ifelse(ndt==1,
               Tri<-Yp[Ytri,],
               Tri<-Yp[Ytri[i,],])
        #vertices of ith triangle
        tri<-pcds::as.basic.tri(Tri)$tri
        #convert the triangle Tri into an nonscaled basic triangle,
        #see as.basic.tri help page
        ni.vec<-c(ni.vec,length(Xpi)/2)
        #number of points in ith Delaunay triangle

        ifelse(identical(M,"CC"),
               cent<-pcds::circumcenter.tri(tri),
               cent<-M)
        num.edges<-num.edgesAStri(Xpi,tri,cent,ugraph)$num.edges  #bura
        #number of edges in ith triangle
        edges<-c(edges,num.edges)
        #number of edges in all triangles as A \code{vector}

      }

      Tot.Edges<-sum(edges)  #the total number of edges in all triangles
    }

    Tot.Edges = Tot.Edges + sum(duplicated(Xp[!inCH,]))

     desc<-paste("Number of Edges of the ",und.graph,
                " Graphs of the AS-PCD with vertices Xp and the Related Quantities for the Induced Subdigraphs for the Points in the Delaunay Triangles",sep="")

    res<-list(desc=desc, #description of the output
              und.graph = und.graph, #"Underlying" or "Reflexivity"
              num.edges=Tot.Edges, #number of edges for the entire PCD
              tri.num.edges=edges,
              #vector of number of edges for the Delaunay triangles
              num.in.conv.hull=NinCH,
              # number of Xp points in CH of Yp points
              ind.in.conv.hull= ind.in.CH, #indices of Xp points in CH of Yp points
              num.in.tris=ni.vec,
              # vector of number of Xp points in the Delaunay triangles
              weight.vec=Wvec, #areas of Delaunay triangles
              del.tri.ind=t(Ytri),
              # indices of the Delaunay triangles, each column corresponds to the vector of
              #indices of the vertices of one triangle
              data.tri.ind=data.tri.ind, #indices of the Delaunay triangles in which data
              #points reside, i.e., column number of del.tri.ind for each Xp point
              tess.points=Yp, #tessellation points
              vertices=Xp #vertices of the digraph
    )
  }

  class(res) <- "NumEdges"
  res$call <-match.call()

  res
} #end of the function
#'

#################################################################

#' @title Incidence matrix for the underlying or reflexivity graph of
#' Arc Slice Proximity Catch Digraphs (AS-PCDs) -
#' one triangle case
#'
#' @description Returns the incidence matrix
#' for the underlying or reflexivity graph of the AS-PCD
#' whose vertices are the given 2D numerical data set, \code{Xp},
#' in the triangle \code{tri}\eqn{=T(v=1,v=2,v=3)}.
#'
#' AS proximity regions are defined
#' with respect to the triangle \code{tri}\eqn{=T(v=1,v=2,v=3)} and
#' vertex regions are based on the center, \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri}
#' or based on circumcenter of  \code{tri};
#' default is \code{M="CC"}, i.e., circumcenter of \code{tri}.
#' Loops are allowed, so the diagonal entries are all equal to 1.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graph of the AS-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param M The center of the triangle.
#' \code{"CC"} stands for circumcenter of the triangle \code{tri}
#' or a 2D point in Cartesian coordinates or
#' a 3D point in barycentric coordinates
#' which serves as a center in the interior of \code{tri};
#' default is \code{M="CC"}, i.e., the circumcenter of \code{tri}.
#' @param ugraph The type of the graph based on AS-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return Incidence matrix for the underlying or reflexivity graph
#' of the AS-PCD whose vertices are the 2D data set, \code{Xp}
#' in the triangle \code{tri} with vertex regions based on the center \code{M}
#'
#' @seealso \code{\link{inci.mat.undAS}}, \code{\link{inci.mat.undPEtri}},
#' \code{\link{inci.mat.undCStri}}, and \code{\link[pcds]{inci.matAStri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #\donttest{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#' n<-10
#'
#' set.seed(1)
#' Xp<-pcds::runif.tri(n,Tr)$g
#'
#' M<-as.numeric(pcds::runif.tri(1,Tr)$g)
#' (IM<-inci.mat.undAStri(Xp,Tr,M))
#' pcds::dom.num.greedy(IM)
#' pcds::Idom.num.up.bnd(IM,3)
#'
#' (IM<-inci.mat.undAStri(Xp,Tr,M,ugraph="r"))
#' pcds::dom.num.greedy(IM)
#' pcds::Idom.num.up.bnd(IM,3)
#' #}
#'
#' @export inci.mat.undAStri
inci.mat.undAStri <- function(Xp,tri,M="CC",
                              ugraph=c("underlying", "reflexivity"))
{
  if (!is.numeric(as.matrix(Xp)))
  {stop('Xp must be numeric')}

  if (pcds::is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('The triangle is degenerate')}

  if (!(pcds::is.point(M) || pcds::is.point(M,3) || identical(M,"CC")))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  CC = pcds::circumcenter.tri(tri)
  if (identical(M,"CC") )
  { M<-CC }

  if (pcds::dimension(M)==3)
  {M<-pcds::bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) ||
        pcds::in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('M is not the circumcenter or not a center in the interior of the triangle')}

  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  n<-nrow(Xp)
  inc.mat<-matrix(0, nrow=n, ncol=n)
 # if (n>=1)
 # {
    for (i in 1:n)
    {p1<-Xp[i,]

    for (j in i:n )
    {p2<-Xp[j,]
    inc.mat[i,j]<-inc.mat[j,i]<-IedgeAStri(p1,p2,tri,M,ugraph)
    }
    }
  #}
  inc.mat
} #end of the function
#'

#################################################################

#' @title The edges of the underlying or reflexivity graph of
#' the Arc Slice Proximity Catch Digraph (AS-PCD) for 2D data -
#' one triangle case
#'
#' @description
#' An object of class \code{"UndPCDs"}.
#' Returns edges of the underlying or reflexivity graph of AS-PCD
#' as left and right end points
#' and related parameters and the quantities of these graphs.
#' The vertices of these graphs are the data points in \code{Xp}
#' in the multiple triangle case.
#'
#' AS proximity regions are constructed
#' with respect to the triangle \code{tri}, i.e.,
#' edges may exist only for points inside \code{tri}.
#' It also provides various descriptions
#' and quantities about the edges of
#' the underlying or reflexivity graph of the AS-PCD
#' such as number of edges, edge density, etc.
#'
#' Vertex regions are based on the center, \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri}
#' or based on circumcenter of  \code{tri};
#' default is \code{M="CC"}, i.e., circumcenter of \code{tri}.
#' The different consideration of circumcenter vs
#' any other interior center of the triangle is because
#' the projections from circumcenter are orthogonal to the edges,
#' while projections of \code{M} on the edges are on the extensions
#' of the lines connecting \code{M} and the vertices.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying or reflexivity graph of the AS-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param M The center of the triangle.
#' \code{"CC"} stands for circumcenter of the triangle \code{tri}
#' or a 2D point in Cartesian coordinates or a 3D point in
#' barycentric coordinates
#' which serves as a center in the interior of the triangle \eqn{T_b};
#' default is \code{M="CC"}, i.e., the circumcenter of \code{tri}.
#' @param ugraph The type of the graph based on AS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the underlying
#' or reflexivity graph of the digraph}
#' \item{parameters}{Parameters of the underlying or reflexivity graph of the digraph,
#' here, it is only the center \code{M} used to construct the vertex regions.}
#' \item{tess.points}{Tessellation points, i.e., points on which the tessellation of
#' the study region is performed,
#' here, tessellation is the support triangle.}
#' \item{tess.name}{Name of the tessellation points \code{tess.points}}
#' \item{vertices}{Vertices of the underlying
#' or reflexivity graph of the digraph, \code{Xp} points}
#' \item{vert.name}{Name of the data set
#' which constitutes the vertices of
#' the underlying or reflexivity graph of the digraph}
#' \item{LE, RE}{Left and right end points of the edges of
#' the underlying or reflexivity graph of AS-PCD for 2D data set \code{Xp}
#' as vertices of the underlying or reflexivity graph of the digraph}
#' \item{mtitle}{Text for \code{"main"} title
#' in the plot of the underlying or reflexivity graph of the digraph}
#' \item{quant}{Various quantities for the underlying or reflexivity graph of the digraph:
#' number of vertices, number of partition points,
#' number of intervals, number of edges, and edge density.}
#'
#' @seealso \code{\link{edgesAS}}, \code{\link{edgesPEtri}},
#' \code{\link{edgesCStri}}, and \code{\link[pcds]{arcsAStri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #\donttest{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#' n<-10
#'
#' set.seed(1)
#' Xp<-pcds::runif.tri(n,Tr)$g
#'
#' M<-as.numeric(pcds::runif.tri(1,Tr)$g)
#'
#' #for underlying graph
#' Edges<-edgesAStri(Xp,Tr,M)
#' Edges
#' summary(Edges)
#' plot(Edges)
#'
#' #for reflexivity graph
#' Edges<-edgesAStri(Xp,Tr,M,ugraph="r")
#' Edges
#' summary(Edges)
#' plot(Edges)
#'
#' #can add vertex regions, but we first need to determine center is the circumcenter or not,
#' #see the description for more detail.
#' CC<-pcds::circumcenter.tri(Tr)
#' if (isTRUE(all.equal(M,CC)))
#' {cent<-CC
#' D1<-(B+C)/2; D2<-(A+C)/2; D3<-(A+B)/2;
#' Ds<-rbind(D1,D2,D3)
#' cent.name<-"CC"
#' } else
#' {cent<-M
#' cent.name<-"M"
#' Ds<-pcds::prj.cent2edges(Tr,M)
#' }
#' L<-rbind(cent,cent,cent); R<-Ds
#' segments(L[,1], L[,2], R[,1], R[,2], lty=2)
#'
#' #now we can add the vertex names and annotation
#' txt<-rbind(Tr,cent,Ds)
#' xc<-txt[,1]+c(-.02,.02,.02,.02,.03,-.03,-.01)
#' yc<-txt[,2]+c(.02,.02,.03,.06,.04,.05,-.07)
#' txt.str<-c("A","B","C","M","D1","D2","D3")
#' text(xc,yc,txt.str)
#' #}
#'
#' @export edgesAStri
edgesAStri <- function(Xp,tri,M="CC",
                       ugraph=c("underlying", "reflexivity"))
{
  xname <-deparse(substitute(Xp))
  yname <-deparse(substitute(tri))

  if (!is.numeric(as.matrix(Xp)) )
  {stop('Xp must be numeric')}

  if (pcds::is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  tri<-as.matrix(tri)
  if (!is.numeric(tri) || nrow(tri)!=3 || ncol(tri)!=2)
  {stop('tri must be numeric and of dimension 3x2')}

  vec1<-rep(1,3);
  D0<-det(matrix(cbind(tri,vec1),ncol=3))
  if (round(D0,14)==0)
  {stop('The triangle is degenerate')}

  CC <- pcds::circumcenter.tri(tri)
  if (identical(M,"CC"))
  {M<-CC
  } else
  { if (!pcds::is.point(M) && !pcds::is.point(M,3))
  {stop('M must be the circumcenter "CC" or a numeric 2D point for Cartesian coordinates or
          3D point for barycentric coordinates')}

  if (pcds::dimension(M)==3)
  {M<-pcds::bary2cart(M,tri)}

  if (!(isTRUE(all.equal(M,CC)) ||
        pcds::in.triangle(M,tri,boundary=FALSE)$in.tri))
  {stop('M is not the circumcenter or not a center in the interior of the triangle')}

}
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  n<-nrow(Xp)
  in.tri<-rep(0,n)
  for (i in 1:n)
    in.tri[i]<-pcds::in.triangle(Xp[i,],tri,boundary=TRUE)$in.tri
  #indices the Xp points inside the triangle

  Xtri<-Xp[in.tri==1,] #the Xp points inside the triangle
  n2<-length(Xtri)/2

  #the edges of the underlying or reflexivity graph of AS-PCDs
  lep<-rep<-NULL #left and right end points for the edges
  if (n2>1)
  {
    for (j in 1:(n2-1))
    {
      p1<-Xtri[j,];
      for (k in (j+1):n2)  #to avoid loops
      {
        p2<-Xtri[k,];
        if (IedgeAStri(p1,p2,tri,M,ugraph)==1)
        {
          lep <-rbind(lep,p1); rep <-rbind(rep,p2);
        }
      }
    }
  }

  param<-list(M)
  Mr<-round(M,2)

  if (identical(M,"CC") || isTRUE(all.equal(M,CC)))
  {
    cname <-"CC"
    names(param)<-c("circumcenter")
    main.txt<-paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"),
                    " Graph of AS-PCD with CC-Vertex Regions",sep="")
    typ<-paste(ifelse(ugraph=="underlying","Underlying", "Reflexivity"),
               " Graph of Arc Slice Proximity Catch Digraph (AS-PCD) for 2D Points in the Triangle with CC-Vertex Regions",sep="")
  } else
  {
    cname <-"M"
    names(param)<-c("center")
    main.txt<-paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"),
                    " Graph of AS-PCD with Center ", cname," = (",Mr[1],",",Mr[2],")",sep="")
    typ<-paste(ifelse(ugraph=="underlying","Underlying", "Reflexivity"),
               " Graph of Arc Slice Proximity Catch Digraph (AS-PCD) for 2D Points in the Triangle with Center ",
               cname," = (",Mr[1],",",Mr[2],")",sep="")
  }

  nvert<-n2; ny<-3; ntri<-1;
  nedges=ifelse(!is.null(lep),
                nrow(rep),
                0);
  edge.dens<-ifelse(nvert>1,
                    nedges/(nvert*(nvert-1)),
                    NA)

  quantities<-c(nvert,ny,ntri,nedges,edge.dens)
  names(quantities)<-c("number of vertices", "number of partition points",
                       "number of triangles","number of edges", "edge density")

  res<-list(
    type=typ,
    parameters=param,
    tess.points=tri, tess.name=yname, #tessellation points
    vertices=Xp, vert.name=xname, #vertices of the digraph
    LE=lep, RE=rep,
    mtitle=main.txt,
    quant=quantities,
    und.graph = ugraph
  )

  class(res) <- "UndPCDs"
  res$call <-match.call()
  res
} #end of the function
#'

#################################################################

#' @title The plot of the edges of the underlying or reflexivity graph of
#' the Arc Slice Proximity Catch Digraph
#' (AS-PCD) for 2D data - one triangle case
#'
#' @description Plots the edges of the underlying or reflexivity graph of
#' the Arc Slice Proximity Catch Digraph (AS-PCD)
#' whose vertices are the data points, \code{Xp}
#' and also the triangle \code{tri}.
#' AS proximity regions are constructed
#' with respect to the triangle \code{tri},
#' only for \code{Xp} points inside the triangle \code{tri}.
#' i.e., edges may exist only for \code{Xp} points inside the triangle \code{tri}.
#'
#' Vertex regions are based on the center, \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri}
#' or based on circumcenter of  \code{tri};
#' default is \code{M="CC"}, i.e., circumcenter of \code{tri}.
#' When the center is the circumcenter, \code{CC},
#' the vertex regions are constructed based on the
#' orthogonal projections to the edges,
#' while with any interior center \code{M},
#' the vertex regions are constructed using the extensions
#' of the lines combining vertices with \code{M}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graphs of the AS-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param M The center of the triangle.
#' \code{"CC"} stands for circumcenter of the triangle \code{tri}
#' or a 2D point in Cartesian coordinates or a 3D point in
#' barycentric coordinates
#' which serves as a center in the interior of the triangle \eqn{T_b};
#' default is \code{M="CC"}, i.e., the circumcenter of \code{tri}.
#' @param ugraph The type of the graph based on AS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#' @param asp A \code{numeric} value,
#' giving the aspect ratio \eqn{y/x} (default is \code{NA}),
#' see the official help page for \code{asp} by
#' typing "\code{? asp}".
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes,
#' respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2,
#' giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param vert.reg A logical argument to add vertex regions to the plot,
#' default is \code{vert.reg=FALSE}.
#' @param \dots	Additional \code{plot} parameters.
#'
#' @return A plot of the edges of the underlying
#' or reflexivity graphs of the AS-PCD
#' whose vertices are the points in data set \code{Xp}
#' and also the triangle \code{tri}
#'
#' @return A plot of the edges of the underlying
#' or reflexivity graphs of the AS-PCD
#' whose vertices are the points in data set \code{Xp}
#' where AS proximity regions
#' are defined with respect to the triangle \code{tri};
#' also plots the triangle \code{tri}
#'
#' @seealso \code{\link{plotASedges}}, \code{\link{plotPEedges.tri}},
#' \code{\link{plotCSedges.tri}}, and \code{\link[pcds]{plotASarcs.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #\donttest{
#' A<-c(1,1); B<-c(2,0); C<-c(1.5,2);
#' Tr<-rbind(A,B,C);
#' n<-10
#'
#' set.seed(1)
#' Xp<-pcds::runif.tri(n,Tr)$g
#'
#' M<-as.numeric(pcds::runif.tri(1,Tr)$g)
#' plotASedges.tri(Xp,Tr,M,vert.reg = TRUE,xlab="",ylab="")
#' plotASedges.tri(Xp,Tr,M,ugraph="r",vert.reg = TRUE,xlab="",ylab="")
#'
#' #can add vertex labels and text to the figure (with vertex regions)
#' ifelse(isTRUE(all.equal(M,pcds::circumcenter.tri(Tr))),
#' {Ds<-rbind((B+C)/2,(A+C)/2,(A+B)/2); cent.name="CC"},
#' {Ds<-pcds::prj.cent2edges(Tr,M); cent.name="M"})
#'
#' txt<-rbind(Tr,M,Ds)
#' xc<-txt[,1]+c(-.02,.02,.02,.02,.04,-0.03,-.01)
#' yc<-txt[,2]+c(.02,.02,.02,.07,.02,.04,-.06)
#' txt.str<-c("A","B","C",cent.name,"D1","D2","D3")
#' text(xc,yc,txt.str)
#' #}
#'
#' @export plotASedges.tri
plotASedges.tri <- function(Xp,tri,M="CC",
                   ugraph=c("underlying", "reflexivity"), asp=NA,
                   main=NULL,xlab=NULL,ylab=NULL, xlim=NULL,ylim=NULL,
                   vert.reg=FALSE,...)
{
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  if( any(duplicated(as.data.frame(Xp))) ) #if there are duplicates for Xp values, only one is taken for each
  {Xp = unique(as.data.frame(Xp))
  warning("There were duplicate Xp values;
          only one value is kept for each duplicate Xp value (to avoid edges of zero length)!")}

  EdgesAS<-edgesAStri(Xp,tri,M,ugraph)
  lep<-EdgesAS$LE; rep<-EdgesAS$RE
  #lep, rep are left and right end points of the edges of the graph
  cent = (EdgesAS$param)$c

  Xp<-matrix(Xp,ncol=2)
  if (is.null(xlim))
  {xlim<-range(tri[,1],Xp[,1],cent[1])}
  if (is.null(ylim))
  {ylim<-range(tri[,2],Xp[,2],cent[2])}

  if ( isTRUE(all.equal(cent,pcds::circumcenter.tri(tri))) )
  {M="CC"}

  if (is.null(main))
  {if (identical(M,"CC")){
  main=paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"),
             " Graph of AS-PCD\n with CC-Vertex Regions",sep="")
  } else {Mr=round(cent,2)
  Mvec= paste(Mr, collapse=",")
  main=paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"),
             " Graph of AS-PCD\n with M = (",Mvec,")",sep="")}
  }

  if (vert.reg)
  {main=c(main,"\n (vertex regions added)")}

  plot(Xp,main=main,asp=asp, xlab=xlab, ylab=ylab,
       xlim=xlim,ylim=ylim,pch=".",cex=3,...)
  polygon(tri,...)
  if (!is.null(lep)) {segments(lep[,1], lep[,2], rep[,1], rep[,2], col= 4)}

  if (vert.reg){
    ifelse(isTRUE(all.equal(cent,pcds::circumcenter.tri(tri))),
           {A=tri[1,];B=tri[2,];C=tri[3,];
           Ds<-rbind((B+C)/2,(A+C)/2,(A+B)/2)},
           Ds<-pcds::prj.cent2edges(tri,M))
    L<-rbind(cent,cent,cent); R<-Ds
    segments(L[,1], L[,2], R[,1], R[,2], lty=2)
  }
} #end of the function
#'

#################################################################

#' @title The edges of the underlying or reflexivity graph of
#' the Arc Slice Proximity Catch Digraph
#' (AS-PCD) for 2D data - multiple triangle case
#'
#' @description
#' An object of class \code{"UndPCDs"}.
#' Returns edges of the underlying or reflexivity graph of AS-PCD
#' as left and right end points
#' and related parameters and the quantities of these graphs.
#' The vertices of these graphs are the data points in \code{Xp}
#' in the multiple triangle case.
#'
#' AS proximity regions are defined
#' with respect to the Delaunay triangles based on
#' \code{Yp} points, i.e.,
#' AS proximity regions are defined only for \code{Xp} points
#' inside the convex hull of \code{Yp} points.
#' That is, edges may exist for points only
#' inside the convex hull of \code{Yp} points.
#' It also provides various descriptions
#' and quantities about the edges of the AS-PCD
#' such as number of edges, edge density, etc.
#'
#' Vertex regions are based on the center \code{M="CC"}
#' for circumcenter of each Delaunay triangle
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of each Delaunay triangle;
#' default is \code{M="CC"}, i.e., circumcenter of each triangle.
#' \code{M} must be entered in barycentric coordinates
#' unless it is the circumcenter.
#' The different consideration of circumcenter vs
#' any other interior center of the triangle is because
#' the projections from circumcenter are orthogonal to the edges,
#' while projections of \code{M} on the edges are on the extensions
#' of the lines connecting \code{M} and the vertices.
#' Each Delaunay triangle is first converted to
#' an (nonscaled) basic triangle so that \code{M} will be the same
#' type of center for each Delaunay triangle
#' (this conversion is not necessary when \code{M} is \eqn{CM}).
#'
#' Convex hull of \code{Yp} is partitioned
#' by the Delaunay triangles based on \code{Yp} points
#' (i.e., multiple triangles are the set of these Delaunay triangles
#' whose union constitutes the
#' convex hull of \code{Yp} points).
#' For the number of edges, loops are not allowed so edges are only possible
#' for points inside the convex hull of \code{Yp} points.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph})
#' for more on the AS-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds.ugraph})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graph of the AS-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param M The center of the triangle.
#' \code{"CC"} represents the circumcenter of each Delaunay triangle
#' or 3D point in barycentric coordinates
#' which serves as a center in the interior of each Delaunay triangle;
#' default is \code{M="CC"}, i.e., the circumcenter of each triangle.
#' \code{M} must be entered in barycentric coordinates
#' unless it is the circumcenter.
#' @param ugraph The type of the graph based on AS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the underlying
#' or reflexivity graph of the digraph}
#' \item{parameters}{Parameters of the underlying or reflexivity graph of the digraph,
#' here, it is only the center \code{M} used to
#' construct the vertex regions.}
#' \item{tess.points}{Tessellation points, i.e., points on which the tessellation
#' of the study region is performed,
#' here, tessellation is the Delaunay triangulation based on \code{Yp} points.}
#' \item{tess.name}{Name of the tessellation points \code{tess.points}}
#' \item{vertices}{Vertices of the digraph, \code{Xp} points}
#' \item{vert.name}{Name of the data set
#' which constitute the vertices of the digraph}
#' \item{LE, RE}{Left and right end points of the edges of
#' the underlying or reflexivity graph of AS-PCD for 2D data set \code{Xp}
#' as vertices of the underlying or reflexivity graph of the digraph}
#' \item{mtitle}{Text for \code{"main"} title
#' in the plot of the underlying or reflexivity graph of the digraph}
#' \item{quant}{Various quantities for the underlying
#' or reflexivity graph of the digraph:
#' number of vertices, number of partition points,
#' number of intervals, number of edges, and edge density.}
#'
#' @seealso \code{\link{edgesAStri}}, \code{\link{edgesPE}},
#' \code{\link{edgesCS}}, and \code{\link[pcds]{arcsAS}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #\donttest{
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-20; ny<-5;
#'
#' set.seed(1)
#' Xp<-cbind(runif(nx,0,1),runif(nx,0,1))
#' Yp<-cbind(runif(ny,0,.25),runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#'
#' M<-c(1,1,1)
#'
#' Edges<-edgesAS(Xp,Yp,M)
#' Edges
#' summary(Edges)
#' plot(Edges)
#' #}
#'
#' @export edgesAS
edgesAS <- function(Xp,Yp,M="CC",ugraph=c("underlying", "reflexivity"))
{
  xname <-deparse(substitute(Xp))
  yname <-deparse(substitute(Yp))

  if (!is.numeric(as.matrix(Xp)) || !is.numeric(as.matrix(Yp)))
  {stop('Xp and Yp must be numeric')}

  if (pcds::is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension n x 2')}
  }

  Yp<-as.matrix(Yp)
  if (ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('Yp must be of dimension kx2 with k>=3')}

  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  if (nrow(Yp)==3)
  {
    res<-edgesAStri(Xp,Yp,M,ugraph)
  } else
  {
    if ((!pcds::is.point(M,3) && M!="CC") || !all(M>0))
    {stop('M must be a numeric 3D point with positive barycentric coordinates or
          "CC" for circumcenter')}

    DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")

    nx<-nrow(Xp)
    ch<-rep(0,nx)
    for (i in 1:nx)
      ch[i]<-interp::in.convex.hull(DTmesh,Xp[i,1],Xp[i,2],strict=FALSE)

    Xch<-matrix(Xp[ch==1,],ncol=2)
    #the Xp points inside the convex hull of Yp

    DTr<-matrix(interp::triangles(DTmesh)[,1:3],ncol=3)  #Delaunay triangles
    nt<-nrow(DTr)
    nx2<-nrow(Xch) #number of Xp points inside the convex hull of Yp points

    lep<-rep<-NULL #left and right end points for the edges
    if (nx2>1)
    {
      i.tr<-rep(0,nx2)  #the vector of indices for the triangles that contain the Xp points
      for (i in 1:nx2)
        for (j in 1:nt)
        {
          tri<-Yp[DTr[j,],]
          if (pcds::in.triangle(Xch[i,],tri,boundary=TRUE)$in.tri )
            i.tr[i]<-j
        }

      for (i in 1:nt)
      {
        Xl<-matrix(Xch[i.tr==i,],ncol=2)
        if (nrow(Xl)>1)
        {
          Yi.Tri<-Yp[DTr[i,],] #vertices of the ith triangle
          Yi.tri<-pcds::as.basic.tri(Yi.Tri)$tri
          #convert the triangle Yi.Tri into an nonscaled basic triangle,
          #see as.basic.tri help page

          ifelse(identical(M,"CC"),
                cent<-pcds::circumcenter.tri(Yi.tri),
                cent<-M)
          nl<-nrow(Xl)
          for (j in 1:(nl-1))
          { for (k in (j+1):nl)  #to avoid loops
            if (IedgeAStri(Xl[j,],Xl[k,],Yi.tri,cent,ugraph)==1 )
            {
              lep <-rbind(lep,Xl[j,]); rep <-rbind(rep,Xl[k,]);
            }
          }
        }
      }
    }

    cname <-ifelse(identical(M,"CC"),"CC","M")
    param<-list(M)
    names(param)<-c("center")

    if (identical(M,"CC")){
      main.txt=paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"),
                     " Graph of AS-PCD\n with Circumcenter",sep="")
      typ<-paste(ifelse(ugraph=="underlying","Underlying", "Reflexivity"),
                 " Graph of Arc Slice Proximity Catch Digraph (AS-PCD) for 2D points in Multiple Triangles with Circumcenter",sep="")
    } else {Mvec= paste(M, collapse=",")
    main.txt=paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"),
                   " Graph of AS-PCD\n with Center ", cname," = (",Mvec,")",sep="")
    typ<-paste(ifelse(ugraph=="underlying","Underlying", "Reflexivity"),
               " Graph of Arc Slice Proximity Catch Digraph (AS-PCD) for 2D points in Multiple Triangles with Center ",
               cname," = (",Mvec,")",sep="")}

    nvert<-nx2; ny<-nrow(Yp); ntri<-nt;
    nedges=ifelse(!is.null(lep),
                  nrow(lep),
                  0);
    edge.dens<-ifelse(nvert>1,
                      nedges/(nvert*(nvert-1)),
                      NA)

    quantities<-c(nvert,ny,ntri,nedges,edge.dens)
    names(quantities)<-c("number of vertices", "number of partition points",
                         "number of triangles","number of edges", "edge density")

    res<-list(
      type=typ,
      parameters=param,
      tess.points=Yp, tess.name=yname, #tessellation points
      vertices=Xp, vert.name=xname, #vertices of the digraph
      LE=lep, RE=rep,
      mtitle=main.txt,
      quant=quantities,
      und.graph = ugraph
    )

    class(res) <- "UndPCDs"
    res$call <-match.call()
  }
res
} #end of the function
#'

#################################################################

#' @title Incidence matrix for the underlying or reflexivity graph of
#' Arc Slice Proximity Catch Digraphs (AS-PCDs) -
#' multiple triangle case
#'
#' @description Returns the incidence matrix
#' for the underlying or reflexivity graph of the AS-PCD
#' whose vertices are the data points in \code{Xp}
#' in the multiple triangle case.
#'
#' AS proximity regions are defined
#' with respect to the Delaunay triangles based on \code{Yp} points and
#' vertex regions are based on the center \code{M="CC"}
#' for circumcenter of each Delaunay triangle or \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the
#' interior of each Delaunay triangle;
#' default is \code{M="CC"}, i.e., circumcenter of each triangle.
#' Loops are allowed, so the diagonal entries are all equal to 1.
#'
#' Each Delaunay triangle is first converted to
#' an (nonscaled) basic triangle so that \code{M} will be the same
#' type of center for each Delaunay triangle
#' (this conversion is not necessary when \code{M} is \eqn{CM}).
#'
#' Convex hull of \code{Yp} is partitioned
#' by the Delaunay triangles based on \code{Yp} points
#' (i.e., multiple triangles are the set of these Delaunay triangles
#' whose union constitutes the
#' convex hull of \code{Yp} points).
#' For the incidence matrix loops are allowed,
#' so the diagonal entries are all equal to 1.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph})
#' for more on the AS-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds.ugraph})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graph of the AS-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param M The center of each triangle.
#' \code{"CC"} stands for circumcenter of each Delaunay triangle
#' or 3D point in barycentric
#' coordinates which serves as a center
#' in the interior of each Delaunay triangle;
#' default is \code{M="CC"}, i.e., the circumcenter of each triangle.
#' @param ugraph The type of the graph based on AS-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return Incidence matrix for the underlying or reflexivity graph
#' of the AS-PCD whose vertices are the 2D data set, \code{Xp}.
#' AS proximity regions are constructed
#' with respect to the Delaunay triangles and \code{M}-vertex regions.
#'
#' @seealso \code{\link{inci.mat.undAStri}}, \code{\link{inci.mat.undPE}},
#' \code{\link{inci.mat.undCS}}, and \code{\link[pcds]{inci.matAS}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #\donttest{
#' nx<-20; ny<-5;
#'
#' set.seed(1)
#' Xp<-cbind(runif(nx,0,1),runif(nx,0,1))
#' Yp<-cbind(runif(ny,0,.25),
#' runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#'
#' M<-c(1,1,1)
#'
#' IM<-inci.mat.undAS(Xp,Yp,M)
#' IM
#' pcds::dom.num.greedy(IM)
#' #}
#'
#' @export inci.mat.undAS
inci.mat.undAS <- function(Xp,Yp,M="CC",ugraph=c("underlying", "reflexivity"))
{
  if (!is.numeric(as.matrix(Xp)) || !is.numeric(as.matrix(Yp)))
  {stop('Xp and Yp must be numeric')}

  if (pcds::is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  Yp<-as.matrix(Yp)
  if (ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('Yp must be of dimension kx2 with k>=3')}

  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  if (nrow(Yp)==3)
  {
    inc.mat<-inci.mat.undAStri(Xp,Yp,M,ugraph)
  } else
  {
    if ((!pcds::is.point(M,3) && M!="CC") || !all(M>0))
    {stop('M must be a numeric 3D point with positive barycentric coordinates
          or "CC" for circumcenter')}

    DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")

    nx<-nrow(Xp)
    ch<-rep(0,nx)
    for (i in 1:nx)
      ch[i]<-interp::in.convex.hull(DTmesh,Xp[i,1],Xp[i,2],strict=FALSE)

    inc.mat<-matrix(0, nrow=nx, ncol=nx)

    DTr<-matrix(interp::triangles(DTmesh)[,1:3],ncol=3)
    nt<-nrow(DTr)  #number of Delaunay triangles

    #if (nx>1)
    #{
      i.tr<-rep(0,nx)
      #the vector of indices for the triangles that contain the Xp points
      for (i in 1:nx)
        for (j in 1:nt)
        {
          tri<-Yp[DTr[j,],]
          if (pcds::in.triangle(Xp[i,],tri,boundary=TRUE)$in.tri )
            i.tr[i]<-j
        }

      for (i in 1:nx)
      { p1<-Xp[i,]
      Yi.tri<-Yp[DTr[1,],]

      if (i.tr[i]!=0)
      {
        Yi.Tri<-Yp[DTr[i.tr[i],],] #vertices of the ith triangle
        Yi.tri<-pcds::as.basic.tri(Yi.Tri)$tri
        #convert the triangle Yi.Tri into an nonscaled basic triangle,
        #see as.basic.tri help page

      }
      ifelse(identical(M,"CC"),
             cent<-pcds::circumcenter.tri(Yi.tri),
             cent<-M)
        for (j in i:nx )
        {p2<-Xp[j,]
        inc.mat[i,j]<-inc.mat[j,i]<-IedgeAStri(p1,p2,Yi.tri,cent,ugraph)
        }
      #}
      }
    #}

   # diag(inc.mat)<-1
  }
  inc.mat
} #end of the function
#'

#################################################################

#' @title The plot of the edges of the underlying or reflexivity graph of
#' the Arc Slice Proximity Catch Digraph
#' (AS-PCD) for 2D data - multiple triangle case
#'
#' @description Plots the edges of the underlying or reflexivity graph of
#' the Arc Slice Proximity Catch Digraph
#' (AS-PCD) whose vertices are the data
#' points in \code{Xp} in the multiple triangle case
#' and the Delaunay triangles based on \code{Yp} points.
#'
#' AS proximity regions are constructed
#' with respect to the Delaunay triangles based on \code{Yp} points, i.e.,
#' AS proximity regions are defined only for \code{Xp} points
#' inside the convex hull of \code{Yp} points.
#' That is, edges may exist for \code{Xp} points
#' only inside the convex hull of \code{Yp} points.
#'
#' Vertex regions are based on the center \code{M="CC"}
#' for circumcenter of each Delaunay triangle
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates in the
#' interior of each Delaunay triangle;
#' default is \code{M="CC"}, i.e., circumcenter of each triangle.
#' When the center is the circumcenter, \code{CC},
#' the vertex regions are constructed based on the
#' orthogonal projections to the edges,
#' while with any interior center \code{M},
#' the vertex regions are constructed using the extensions
#' of the lines combining vertices with \code{M}.
#'
#' Convex hull of \code{Yp} is partitioned by
#' the Delaunay triangles based on \code{Yp} points
#' (i.e., multiple triangles are the set of these Delaunay triangles
#' whose union constitutes the
#' convex hull of \code{Yp} points).
#' Loops are not allowed so edges are only possible
#' for points inside the convex hull of \code{Yp} points.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph})
#' for more on the AS-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds.ugraph})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graphs of the AS-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param M The center of the triangle.
#' \code{"CC"} stands for circumcenter of each Delaunay triangle
#' or 3D point in barycentric coordinates
#' which serves as a center in the interior of each Delaunay triangle;
#' default is \code{M="CC"}, i.e., the circumcenter of each triangle.
#' @param ugraph The type of the graph based on AS-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#' @param asp A \code{numeric} value,
#' giving the aspect ratio \eqn{y/x} (default is \code{NA}),
#' see the official help page for \code{asp} by typing "\code{? asp}".
#' @param main An overall title for the plot (default=\code{NULL}).
#' @param xlab,ylab Titles for the \eqn{x} and \eqn{y} axes,
#' respectively (default=\code{NULL} for both).
#' @param xlim,ylim Two \code{numeric} vectors of length 2,
#' giving the \eqn{x}- and \eqn{y}-coordinate ranges
#' (default=\code{NULL} for both).
#' @param \dots Additional \code{plot} parameters.
#'
#' @return A plot of the edges of the underlying
#' or reflexivity graphs of the AS-PCD for a 2D data set \code{Xp}
#' where AS proximity regions
#' are defined with respect to the Delaunay triangles based on \code{Yp} points;
#' also plots the Delaunay triangles
#' based on \code{Yp} points.
#'
#' @seealso \code{\link{plotASedges.tri}}, \code{\link{plotPEedges}},
#' \code{\link{plotCSedges}}, and \code{\link[pcds]{plotASarcs}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #\donttest{
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-20; ny<-5;
#'
#' set.seed(1)
#' Xp<-cbind(runif(nx,0,1),runif(nx,0,1))
#' Yp<-cbind(runif(ny,0,.25),
#' runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#'
#' M<-c(1,1,1)
#'
#' plotASedges(Xp,Yp,M,xlab="",ylab="")
#' plotASedges(Xp,Yp,M,xlab="",ylab="",ugraph="r")
#' #}
#'
#' @export plotASedges
plotASedges <- function(Xp,Yp,M="CC",ugraph=c("underlying", "reflexivity"),
                           asp=NA,main=NULL,xlab=NULL,ylab=NULL,
                           xlim=NULL,ylim=NULL,...)
{
  Yp<-as.matrix(Yp)
  if (ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('Yp must be of dimension kx2 with k>=3')}

  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  if( any(duplicated(as.data.frame(Xp))) )
    #if there are duplicates for Xp values, only one is taken for each
  {Xp = unique(as.data.frame(Xp))
  warning("There were duplicate Xp values;
          only one value is kept for each duplicate Xp value (to avoid edges of zero length)!")}

  if (nrow(Yp)==3)
  {
    plotASedges.tri(Xp,Yp,M,ugraph,asp,main,xlab,ylab,xlim,ylim)
  } else
  {
    EdgesAS<-edgesAS(Xp,Yp,M,ugraph)
    lep<-EdgesAS$LE; rep<-EdgesAS$RE

    DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
    Xch<-pcds::Xin.convex.hullY(Xp,Yp)

    if (is.null(main))
    {if (identical(M,"CC")){
      main=paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"),
                 " Graph of AS-PCD\n with Circumcenter",sep="")
    } else {Mvec= paste(M, collapse=",")
    main=paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"),
               " Graph of AS-PCD\n with M = (",Mvec,")",sep="")}
    }

    Xlim<-xlim; Ylim<-ylim
    if (is.null(xlim))
    {xlim<-range(Yp[,1],Xp[,1])
    xr<-xlim[2]-xlim[1]
    xlim<-xlim+xr*c(-.05,.05)
    }

    if (is.null(ylim))
    {ylim<-range(Yp[,2],Xp[,2])
    yr<-ylim[2]-ylim[1]
    ylim<-ylim+yr*c(-.05,.05)
    }
    plot(rbind(Xp),asp=asp,main=main, xlab=xlab, ylab=ylab,xlim=xlim,
         ylim=ylim,pch=".",cex=3,...)
    interp::plot.triSht(DTmesh, add=TRUE, do.points = TRUE)
    if (!is.null(lep)) {segments(lep[,1], lep[,2], rep[,1], rep[,2], col= 4)}
  }
} #end of the function


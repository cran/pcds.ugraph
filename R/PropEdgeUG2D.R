#PropEdgeUG2D.R;
#Functions for Underlying or Reflexivity Graphs of PE-PCD in R^2
#################################################################

#' @title The indicator for the presence of an edge from a point to another
#' for the underlying or reflexivity graph of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' standard basic triangle case
#'
#' @description Returns \eqn{I(}\code{p1p2} is an edge
#' in the underlying or reflexivity graph of PE-PCDs \eqn{)}
#' for points \code{p1} and \code{p2} in the standard basic triangle.
#'
#' More specifically, when the argument \code{ugraph="underlying"}, it returns
#' the edge indicator for the PE-PCD underlying graph,
#' that is, returns 1 if \code{p2} is
#' in \eqn{N_{PE}(p1,r)} **or** \code{p1} is in \eqn{N_{PE}(p2,r)},
#' returns 0 otherwise.
#' On the other hand,
#' when \code{ugraph="reflexivity"}, it returns
#' the edge indicator for the PE-PCD reflexivity graph,
#' that is, returns 1 if \code{p2} is
#' in \eqn{N_{PE}(p1,r)} **and** \code{p1} is in \eqn{N_{PE}(p2,r)},
#' returns 0 otherwise.
#'
#' In both cases \eqn{N_{PE}(x,r)} is the PE proximity region for point \eqn{x}
#' with expansion parameter \eqn{r \ge 1}.
#' PE proximity region is defined
#' with respect to the standard basic triangle \eqn{T_b=T((0,0),(1,0),(c_1,c_2))}
#' where \eqn{c_1} is
#' in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' Vertex regions are based on the center, \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the standard basic triangle \eqn{T_b}
#' or based on circumcenter of \eqn{T_b};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_b}.
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
#' See also
#' (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param p1 A 2D point whose PE proximity region is constructed.
#' @param p2 A 2D point. The function determines
#' whether there is an edge from \code{p1} to \code{p2} or not
#' in the underlying or reflexivity graph of PE-PCDs.
#' @param r A positive real number
#' which serves as the expansion parameter
#' in PE proximity region; must be \eqn{\ge 1}
#' @param c1,c2 Positive real numbers
#' which constitute the vertex of the standard basic triangle
#' adjacent to the shorter edges;
#' \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the standard basic triangle
#' or circumcenter of \eqn{T_b}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_b}.
#' @param ugraph The type of the graph based on PE-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return Returns 1 if there is an edge between points \code{p1} and \code{p2}
#' in the underlying or reflexivity graph of PE-PCDs
#' in the standard basic triangle, and 0 otherwise.
#'
#' @seealso \code{\link{IedgePEtri}}, \code{\link{IedgeASbasic.tri}},
#' \code{\link{IedgeCSbasic.tri}} and \code{\link[pcds]{IarcPEbasic.tri}}
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
#'
#' r<-1.5
#' set.seed(4)
#'
#' P1<-as.numeric(pcds::runif.basic.tri(1,c1,c2)$g)
#' P2<-as.numeric(pcds::runif.basic.tri(1,c1,c2)$g)
#' IedgePEbasic.tri(P1,P2,r,c1,c2,M)
#' IedgePEbasic.tri(P1,P2,r,c1,c2,M,ugraph = "reflexivity")
#'
#' P1<-c(.4,.2)
#' P2<-c(.5,.26)
#' IedgePEbasic.tri(P1,P2,r,c1,c2,M)
#' IedgePEbasic.tri(P1,P2,r,c1,c2,M,ugraph = "reflexivity")
#' #}
#'
#' @export IedgePEbasic.tri
IedgePEbasic.tri <- function(p1,p2,r,c1,c2,M=c(1,1,1),
                             ugraph=c("underlying", "reflexivity"))
{
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  arc12 = pcds::IarcPEbasic.tri(p1,p2,r,c1,c2,M)
  arc21 = pcds::IarcPEbasic.tri(p2,p1,r,c1,c2,M)

  edge = ifelse(ugraph == "underlying",
                max(arc12,arc21),
                arc12*arc21)
  edge
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an edge from a point to another
#' for the underlying or reflexivity graph of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' standard equilateral triangle case
#'
#' @description Returns \eqn{I(}\code{p1p2} is an edge
#' in the underlying or reflexivity graph of PE-PCDs \eqn{)}
#' for points \code{p1} and \code{p2} in the standard equilateral triangle.
#'
#' More specifically,
#' when the argument \code{ugraph="underlying"}, it returns
#' the edge indicator for points \code{p1} and \code{p2}
#' in the standard equilateral triangle,
#' for the PE-PCD underlying graph,
#' that is, returns 1 if \code{p2} is
#' in \eqn{N_{PE}(p1,r)} **or** \code{p1} is in \eqn{N_{PE}(p2,r)},
#' returns 0 otherwise.
#' On the other hand,
#' when \code{ugraph="reflexivity"}, it returns
#' the edge indicator for points \code{p1} and \code{p2}
#' in the standard equilateral triangle,
#' for the PE-PCD reflexivity graph,
#' that is, returns 1 if \code{p2} is
#' in \eqn{N_{PE}(p1,r)} **and** \code{p1} is in \eqn{N_{PE}(p2,r)},
#' returns 0 otherwise.
#'
#' In both cases \eqn{N_{PE}(x,r)} is the PE proximity region
#' for point \eqn{x} with expansion parameter \eqn{r \ge 1}.
#' PE proximity region is defined
#' with respect to the standard equilateral triangle
#' \eqn{T_e=T(v=1,v=2,v=3)=T((0,0),(1,0),(1/2,\sqrt{3}/2))}
#' and vertex regions are based on the center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_e}.
#'
#' If \code{p1} and \code{p2} are distinct
#' and either of them are outside \eqn{T_e}, it returns 0,
#' but if they are identical,
#' then it returns 1 regardless of their locations (i.e., it allows loops).
#'
#' See also
#' (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param p1 A 2D point whose PE proximity region is constructed.
#' @param p2 A 2D point. The function determines
#' whether there is an edge from \code{p1} to \code{p2} or not
#' in the underlying or reflexivity graph of PE-PCDs.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center
#' in the interior of the standard equilateral triangle \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \eqn{T_e}.
#' @param ugraph The type of the graph based on PE-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return Returns 1 if there is an edge between points \code{p1} and \code{p2}
#' in the underlying or reflexivity graph of PE-PCDs
#' in the standard equilateral triangle, and 0 otherwise.
#'
#' @seealso \code{\link{IedgePEbasic.tri}}, \code{\link{IedgePEtri}},
#' and \code{\link[pcds]{IarcPEstd.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #\donttest{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C)
#' n<-3
#'
#' set.seed(1)
#' Xp<-pcds::runif.std.tri(n)$gen.points
#'
#' M<-as.numeric(pcds::runif.std.tri(1)$g)
#'
#' IedgePEstd.tri(Xp[1,],Xp[3,],r=1.5,M)
#' IedgePEstd.tri(Xp[1,],Xp[3,],r=1.5,M,ugraph="reflexivity")
#'
#' P1<-c(.4,.2)
#' P2<-c(.5,.26)
#' r<-2
#' IedgePEstd.tri(P1,P2,r,M)
#' IedgePEstd.tri(P1,P2,r,M,ugraph = "reflexivity")
#' #}
#'
#' @export IedgePEstd.tri
IedgePEstd.tri <- function(p1,p2,r,M=c(1,1,1),ugraph=c("underlying", "reflexivity"))
{
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  arc12 = pcds::IarcPEstd.tri(p1,p2,r,M)
  arc21 = pcds::IarcPEstd.tri(p2,p1,r,M)

  edge = ifelse(ugraph == "underlying",
                max(arc12,arc21),
                arc12*arc21)
  edge
} #end of the function
#'

#################################################################

#' @title Number of edges in the underlying or reflexivity graph of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' standard equilateral triangle case
#'
#' @description
#' An object of class \code{"NumEdges"}.
#' Returns the number of edges of
#' the underlying or reflexivity graph of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs)
#' whose vertices are the
#' given 2D numerical data set, \code{Xp}
#' in the standard equilateral triangle.
#' It also provides number of vertices
#' (i.e., number of data points inside the triangle)
#' and indices of the data points that reside in the triangle.
#'
#' PE proximity region \eqn{N_{PE}(x,r)} is defined
#' with respect to the standard equilateral triangle
#' \eqn{T_e=T(v=1,v=2,v=3)=T((0,0),(1,0),(1/2,\sqrt{3}/2))}
#' with expansion parameter \eqn{r \ge 1}
#' and vertex regions are based on the center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of \eqn{T_e};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \eqn{T_e}.
#' For the number of edges, loops are not allowed so
#' edges are only possible for points inside \eqn{T_e} for this function.
#'
#' See also (\insertCite{ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graphs based on the PE-PCD.
#' @param r A positive real number
#' which serves as the expansion parameter for PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center
#' in the interior of the standard equilateral triangle \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \eqn{T_e}.
#' @param ugraph The type of the graph based on PE-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of edges
#' and quantities related to the standard equilateral triangle}
#' \item{und.graph}{Type of the graph as "Underlying"
#' or "Reflexivity" for the PE-PCD}
#' \item{num.edges}{Number of edges of the underlying
#' or reflexivity graphs based on the PE-PCD
#' with vertices in the standard equilateral triangle \eqn{T_e}}
#' \item{num.in.tri}{Number of \code{Xp} points
#' in the standard equilateral triangle, \eqn{T_e}}
#' \item{ind.in.tri}{The vector of indices of the \code{Xp} points
#' that reside in \eqn{T_e}}
#' \item{tess.points}{Tessellation points, i.e., points on which
#' the tessellation of the study region is performed,
#' here, tessellation is the support triangle \eqn{T_e}.}
#' \item{vertices}{Vertices of the underlying or reflexivity graph, \code{Xp}.}
#'
#' @seealso \code{\link{num.edgesPEtri}}, \code{\link{num.edgesPE}},
#' and \code{\link[pcds]{num.arcsPEstd.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #\donttest{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' n<-10
#'
#' set.seed(1)
#' Xp<-pcds::runif.std.tri(n)$gen.points
#'
#' M<-c(.6,.2)
#'
#' Nedges = num.edgesPEstd.tri(Xp,r=1.25,M)
#' Nedges
#' summary(Nedges)
#' plot(Nedges)
#' #}
#'
#' @export num.edgesPEstd.tri
num.edgesPEstd.tri <- function(Xp,r,M=c(1,1,1),ugraph=c("underlying", "reflexivity"))
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

  if (!pcds::is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!pcds::is.point(M) && !pcds::is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates ')}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C);

  if (pcds::dimension(M)==3)
  {M<-pcds::bary2cart(M,Te)}

  if (pcds::in.triangle(M,Te,boundary=FALSE)$in.tri==F)
  {stop('M is not a center in the interior of the triangle')}

  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  n<-nrow(Xp)
  edges<-0
  ind.in.tri = NULL
  if (n<=0)
  {
    edges<-0
  } else
  {
    for (i in 1:n)
    {p1<-Xp[i,]
    if (!pcds::in.triangle(p1,Te,boundary = TRUE)$in.tri)
    {edges<-edges+0
    } else
    {
      ind.in.tri = c(ind.in.tri,i)
      for (j in i:n )
      {p2<-Xp[j,]
      if (!pcds::in.triangle(p2,Te,boundary = TRUE)$in.tri)
      {edges<-edges+0
      } else
      {
        edges<-edges+IedgePEstd.tri(p1,p2,r,M,ugraph)
      }
      }
    }
    }
  }

  NinTri = length(ind.in.tri)

  und.graph = ifelse(ugraph=="underlying","Underlying", "Reflexivity")
  desc<-paste("Number of Edges of the ",und.graph,
              " Graphs of the PE-PCD and the Related Quantities with vertices Xp in the Standard Equilateral Triangle",sep="")

  res<-list(desc=desc, #description of the output
            und.graph = und.graph, #"Underlying" or "Reflexivity"
            num.edges=edges-n, # -n is to avoid loops
            #number of edges for the underlying or reflexivity graph of the PE-PCD
            num.in.tri=NinTri, # number of Xp points in CH of Yp points
            ind.in.tri=ind.in.tri, #indices of data points inside the triangle
            tess.points=Te, #tessellation points
            vertices=Xp #vertices of the underlying or reflexivity graph
  )

  class(res) <- "NumEdges"
  res$call <-match.call()

  res
} #end of the function
#'

#################################################################

#' @title Incidence matrix for the underlying or reflexivity graph of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' standard equilateral triangle case
#'
#' @description Returns the incidence matrix
#' for the underlying or reflexivity graph of the PE-PCD
#' whose vertices are the given 2D numerical data set, \code{Xp},
#' in the standard equilateral triangle
#' \eqn{T_e=T(v=1,v=2,v=3)=T((0,0),(1,0),(1/2,\sqrt{3}/2))}.
#'
#' PE proximity region is constructed
#' with respect to the standard equilateral triangle \eqn{T_e} with
#' expansion parameter \eqn{r \ge 1} and vertex regions are based on
#' the center \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of \eqn{T_e}; default is \eqn{M=(1,1,1)},
#' i.e., the center of mass of \eqn{T_e}.
#' Loops are allowed,
#' so the diagonal entries are all equal to 1.
#'
#' See also
#' (\insertCite{ceyhan:Phd-thesis,ceyhan:comp-geo-2010,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying or
#' reflexivity graph of the PE-PCD.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center
#' in the interior of the standard equilateral triangle \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \eqn{T_e}.
#' @param ugraph The type of the graph based on PE-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return Incidence matrix for the underlying or reflexivity graph
#' of the PE-PCD with vertices
#' being 2D data set, \code{Xp}
#' in the standard equilateral triangle where PE proximity
#' regions are defined with \code{M}-vertex regions.
#'
#' @seealso \code{\link{inci.mat.undPEtri}}, \code{\link{inci.mat.undPE}},
#' and \code{\link[pcds]{inci.matPEstd.tri}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #\donttest{
#' A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
#' Te<-rbind(A,B,C)
#' n<-10
#'
#' set.seed(1)
#' Xp<-pcds::runif.std.tri(n)$gen.points
#'
#' M<-as.numeric(pcds::runif.std.tri(1)$g)
#'
#' inc.mat<-inci.mat.undPEstd.tri(Xp,r=1.25,M)
#' inc.mat
#' (sum(inc.mat)-n)/2
#' num.edgesPEstd.tri(Xp,r=1.25,M)$num.edges
#'
#' pcds::dom.num.greedy(inc.mat)
#' pcds::Idom.num.up.bnd(inc.mat,2)
#' #}
#'
#' @export inci.mat.undPEstd.tri
inci.mat.undPEstd.tri <- function(Xp,r,M=c(1,1,1),
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

  if (!pcds::is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!pcds::is.point(M) && !pcds::is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates ')}

  A<-c(0,0); B<-c(1,0); C<-c(1/2,sqrt(3)/2);
  Te<-rbind(A,B,C);

  if (pcds::dimension(M)==3)
  {M<-pcds::bary2cart(M,Te)}

  if (pcds::in.triangle(M,Te,boundary=FALSE)$in.tri==F)
  {stop('M is not a center in the interior of the triangle')}

  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  n<-nrow(Xp)
  inc.mat<-matrix(0, nrow=n, ncol=n)
  for (i in 1:n)
  {p1<-Xp[i,]
  for (j in i:n )
  {p2<-Xp[j,]
  inc.mat[i,j]<-inc.mat[j,i]<-IedgePEstd.tri(p1,p2,r,M,ugraph)
  }
  }
  inc.mat
} #end of the function
#'

#################################################################

# funsMuVarUndPE2D
#'
#' @title Returns the mean and (asymptotic) variance of edge density of
#' underlying or reflexivity graph of
#' Proportional Edge Proximity Catch Digraph (PE-PCD)
#' for 2D uniform data in one triangle
#'
#' @description
#' The mean and (asymptotic) variance functions
#' for the underlying or reflexivity graph of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs):
#' \code{muOrPE2D} and \code{asy.varOrPE2D} for the underlying graph
#' and
#' \code{muAndPE2D} and \code{asy.varAndPE2D} for the reflexivity graph.
#'
#' \code{muOrPE2D} and \code{muAndPE2D} return the mean of the (edge) density of
#' the underlying or reflexivity graph of PE-PCDs, respectively,
#' for 2D uniform data in a triangle.
#' Similarly,
#' \code{asy.varOrPE2D} and \code{asy.varAndPE2D} return the asymptotic variance
#' of the edge density of the underlying or reflexivity graph of PE-PCDs,
#' respectively, for 2D uniform data in a triangle.
#'
#' PE proximity regions are defined with expansion parameter \eqn{r \ge 1}
#' with respect to the triangle in which the points reside and
#' vertex regions are based on center of mass, \eqn{CM} of the triangle.
#'
#' See also (\insertCite{ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param r A positive real number which serves
#' as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param ugraph The type of the graph based on PE-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return \code{mu.undPE2D} returns the mean
#' and \code{asy.varUndOrPE2D} returns the (asymptotic) variance of the
#' edge density of the underlying graph of the PE-PCD for
#' uniform data in any triangle
#' if \code{ugraph="underlying"}, and those of the reflexivity graph
#' if \code{ugraph="reflexivity"}.
#' The functions \code{muOrPE2D}, \code{muAndPE2D}, \code{asy.varOrPE2D},
#' and \code{asy.varAndPE2D} are the corresponding mean
#' and asymptotic variance functions
#' for the edge density of the reflexivity graph of the PE-PCD,
#' respectively, for uniform data in any triangle.
#'
#' @name funsMuVarUndPE2D
NULL
#'
#' @seealso \code{\link{mu.undCS2D}}, \code{\link{asy.var.undCS2D}},
#' \code{\link[pcds]{muPE2D}}, \code{\link[pcds]{asy.varPE2D}},
#' \code{\link{muAndCS2D}}, and \code{\link{asy.varAndCS2D}}
#'
#' @rdname funsMuVarUndPE2D
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #\donttest{
#' mu.undPE2D(1.2)
#' mu.undPE2D(1.2,ugraph="r")
#'
#' rseq<-seq(1.01,5,by=.05)
#' lrseq<-length(rseq)
#'
#' muOR = muAND <- vector()
#' for (i in 1:lrseq)
#' {
#'   muOR<-c(muOR,mu.undPE2D(rseq[i]))
#'   muAND<-c(muAND,mu.undPE2D(rseq[i],ugraph="r"))
#' }
#'
#' plot(rseq, muOR,type="l",xlab="r",ylab=expression(mu(r)),lty=1,
#'      xlim=range(rseq),ylim=c(0,1))
#' lines(rseq,muAND,type="l",lty=2,col=2)
#' legend("bottomright", inset=.02,
#'        legend=c(expression(mu[or](r)),expression(mu[and](r))),
#'        lty=1:2,col=1:2)
#' #}
#'
#' @export muOrPE2D
muOrPE2D<-function(r)
{
  if (!pcds::is.point(r,1) || r<1)
  {stop('The argument must be a scalar greater than 1')}

  mn<-0;
  if (r < 4/3)
  {
    mn<- (-108*r^2+860*r^4-846*r^3+720*r-195*r^5-256+47*r^6)/(108*r^2*(r+2)*(r+1));
  } else {
    if (r < 3/2)
    {
      mn<- (175*r^5-579*r^4-536*r+1450*r^3-732*r^2+672)/(216*r*(r+2)*(r+1));
    } else {
      if (r<2)
      {
        mn<- -(3*r^8-7*r^7-30*r^6+84*r^5-264*r^4+304*r^3+144*r^2-368*r+96)/(8*r^4*(r+2)*(r+1));
      } else {
        mn<- (r^5+r^4-6*r+2)/(r^4*(r+1));
      }}}
  mn
} #end of the function
#'
#' @rdname funsMuVarUndPE2D
#'
#' @export muAndPE2D
muAndPE2D <- function(r)
{
  if (!pcds::is.point(r,1) || r<1)
  {stop('The argument must be a scalar greater than 1')}

  mn<-0;
  if (r < 4/3)
  {
    mn<- -(r-1)*(5*r^5-148*r^4+245*r^3-178*r^2-232*r+128)/(54*r^2*(r+2)*(r+1));
  } else {
    if (r < 3/2)
    {
      mn<- -(101*r^5-801*r^4+1302*r^3-732*r^2-536*r+672)/(216*r*(r+2)*(r+1));
    } else {
      if (r<2)
      {
        mn<- (r^8-13*r^7+30*r^6+148*r^5-448*r^4+264*r^3+288*r^2-368*r+96)/(8*r^4*(r+2)*(r+1));
      } else {
        mn<- (r^3+3*r^2-2+2*r)*(r-1)^2/(r^4*(r+1));
      }}}
  mn
} #end of the function
#'
#' @rdname funsMuVarUndPE2D
#'
#' @export mu.undPE2D
mu.undPE2D <-function(r,ugraph=c("underlying", "reflexivity"))
{
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  mn = ifelse(ugraph == "underlying",muOrPE2D(r),muAndPE2D(r))
  mn
}
#' @references
#' \insertAllCited{}
#'
#' @rdname funsMuVarUndPE2D
#'
#' @examples
#' #\donttest{
#' asy.var.undPE2D(1.2)
#' asy.var.undPE2D(1.2,ugraph="r")
#'
#' rseq<-seq(1.01,5,by=.05)
#' lrseq<-length(rseq)
#'
#' avarOR<-avarAND<-vector()
#' for (i in 1:lrseq)
#' {
#'   avarOR<-c(avarOR,asy.var.undPE2D(rseq[i]))
#'   avarAND<-c(avarAND,asy.var.undPE2D(rseq[i],ugraph="r"))
#' }
#'
#' oldpar <- par(mar=c(5,5,4,2))
#' plot(rseq, avarAND,type="l",lty=2,col=2,xlab="r",
#'      ylab=expression(paste(sigma^2,"(r)")),xlim=range(rseq))
#' lines(rseq,avarOR,type="l")
#' legend(3.75,.02,
#'        legend=c(expression(paste(sigma["underlying"]^"2","(r)")),
#'                  expression(paste(sigma["reflexivity"]^"2","(r)")) ),
#'        lty=1:2,col=1:2)
#'
#' par(oldpar)
#' #}
#'
#' @export asy.varOrPE2D
asy.varOrPE2D<-function(r)
{
  if (!pcds::is.point(r,1) || r<1)
  {stop('The argument must be a scalar greater than 1')}

  #internal OR underlying functions for covariance (i.e. asymptotic variance)
  PEcov.or0 <-function(r)
  {-(1458*r^22+13122*r^21+50731*r^20-84225*r^19-19193*r^18-1823223*r^17+5576151*r^16+2978697*r^15-
       33432692*r^14+37427862*r^13+15883834*r^12-60944766*r^11+49876417*r^10-1754523*r^9-36606859*r^8+
       32338215*r^7-10290256*r^6-2234754*r^5+7085471*r^4-5608569*r^3+1645826*r^2-132876*r+
       30824)/(58320*(2*r^2+1)*(r+2)^3*(r^2+1)*(r+1)^3*r^6)}

  PEcov.or1 <-function(r)
  {-(1458*r^22+13122*r^21+62825*r^20-175011*r^19+156014*r^18-3300900*r^17+11053023*r^16+5055135*r^15-
       67685050*r^14+75243552*r^13+33155180*r^12-120628524*r^11+99831906*r^10-4883958*r^9-74801558*r^8+
       64360782*r^7-19812000*r^6-3667716*r^5+14541630*r^4-11254002*r^3+3070468*r^2-413208*r+
       28880)/(r^2+1)/(116640*(2*r^2+1)*(r+2)^3*(r+1)^3*r^6)}

  PEcov.or2 <-function(r)
  {-(1458*r^22+13122*r^21+62825*r^20-175011*r^19+156014*r^18-3300900*r^17+11053023*r^16+5055135*r^15-
       67685050*r^14+75243552*r^13+33155180*r^12-120628524*r^11+99831906*r^10-4883958*r^9-74801558*r^8+
       64360782*r^7-19812000*r^6-3667716*r^5+14541630*r^4-11254002*r^3+3070468*r^2-413208*r+
       28880)/(116640*(r^2+1)*(2*r^2+1)*(r+2)^3*(r+1)^3*r^6)}

  PEcov.or3 <-function(r)
  {-(972*r^24+8748*r^23+29590*r^22-149106*r^21-36820*r^20-986280*r^19+5942884*r^18+2883672*r^17-
       47189711*r^16+43450125*r^15+85975304*r^14-156173934*r^13+27378901*r^12+123606417*r^11-152209261*r^10+
       64653597*r^9+56621894*r^8-88962768*r^7+43754559*r^6-5940597*r^5-13006396*r^4+17019366*r^3-7037340*r^2+
       413208*r-28880)/(58320*(r^2+1)*(2*r^2+1)*(r^2-2)*(r+2)^3*(r+1)^3*r^6)}

  PEcov.or4 <-function(r)
  {-(972*r^22+8748*r^21+31534*r^20-131610*r^19+261546*r^18-1552026*r^17+3745643*r^16+4573731*r^15-
       29416804*r^14+26163354*r^13+19600850*r^12-43126062*r^11+31497249*r^10-7381467*r^9-22237963*r^8+
       26778663*r^7-9107024*r^6-115074*r^5+3136927*r^4-5055609*r^3+2292994*r^2+14580*r-
       1944)/(58320*(r^2+1)*(2*r^2+1)*(r+1)^3*(r+2)^3*r^6)}

  PEcov.or5 <-function(r)
  {(486*r^22-7290*r^21-181459*r^20+1024401*r^19-2691213*r^18+3921057*r^17+1844321*r^16-33347697*r^15+
      80028903*r^14-29292735*r^13-98093906*r^12+125034492*r^11-46658244*r^10-57216612*r^9+88057996*r^8-
      26383068*r^7-12851392*r^6+14179848*r^5-8656508*r^4+1593828*r^3+134136*r^2-
      58320*r+7776)/(233280*(r^2+1)*(2*r^2+1)*(r+1)^3*(r+2)^3*r^6)}

  PEcov.or6 <-function(r)
  {(486*r^23-7776*r^22-174169*r^21+1205860*r^20-4656806*r^19+8763566*r^18+7460036*r^17-63559490*r^16+
      91134324*r^15+18516450*r^14-122708655*r^13+18577230*r^12+80410332*r^11-19357704*r^10-39129236*r^9+
      75311048*r^8-77449360*r^7+4053376*r^6+48283912*r^5-40690240*r^4+17736336*r^3-4315680*r^2+
      544320*r-31104)/(233280*(r+2)^3*(r^2+1)*(2*r^2+1)*(r+1)^3*(r-1)*r^6)}

  PEcov.or7 <-function(r)
  {(2*r^24-30*r^23-161*r^22+107*r^21+4137*r^20-10685*r^19+8367*r^18+78713*r^17-450859*r^16+697707*r^15+
      517846*r^14-3723120*r^13+6565124*r^12-1468692*r^11-8695792*r^10+9535720*r^9-6773160*r^8+526744*r^7+
      10691376*r^6-7797264*r^5+1137696*r^4+523712*r^3-2687872*r^2+1701888*r-
      245760)/(960*(2*r^2+1)*(r+2)^3*(r^2+1)*(r+1)^3*r^8)}

  PEcov.or8 <-function(r)
  {(2*r^25-32*r^24-129*r^23+236*r^22+4157*r^21-15610*r^20+21289*r^19+67536*r^18-511355*r^17+1161830*r^16-
      634128*r^15-3001568*r^14+9512164*r^13-11014136*r^12+2344968*r^11+7126240*r^10-13850504*r^9+
      14466592*r^8-3823216*r^7-4018976*r^6+5155776*r^5-4633984*r^4+1959808*r^3-244480*r^2-
      3584*r-1024)/(960*(2*r^2+1)*(r+1)^2*(r+2)^3*(r^2+1)*r^10)}

  PEcov.or9 <-function(r)
  {(2*r^24-34*r^23-101*r^22+433*r^21+5400*r^20-26982*r^19+23049*r^18+166787*r^17-717366*r^16+1196092*r^15+
      89468*r^14-5130844*r^13+12748688*r^12-11274744*r^11-12243496*r^10+33980568*r^9-14886656*r^8-19910592*r^7+
      20667776*r^6-1262208*r^5-5402752*r^4+2217088*r^3-235776*r^2-2560*r-
      1024)/(960*(2*r^2-1)*(r+2)^3*(r-1)*(r+1)^2*r^10)}

  PEcov.or10 <-function(r)
  {2*(180*r^8-48*r^7-648*r^6+396*r^5+214*r^4-190*r^3+39*r^2-4*r+1)/(15*(2*r^2-1)*(r+1)^2*r^10)}
  ###
  asy.var<-0;
  if (r < 2*sqrt(3)/3)
  {
    asy.var<-PEcov.or0(r);
  } else {
    if (r < 6/5)
    {
      asy.var<- PEcov.or1(r);
    } else {
      if (r < sqrt(5)-1)
      {
        asy.var<- PEcov.or2(r);
      } else {
        if (r < (6+2*sqrt(2))/7)
        {
          asy.var<- PEcov.or3(r);
        } else {
          if (r < 4/3)
          {
            asy.var<- PEcov.or4(r);
          } else {
            if (r < (6+sqrt(15))/7)
            {
              asy.var<- PEcov.or5(r);
            } else {
              if (r < 3/2)
              {
                asy.var<- PEcov.or6(r);
              } else {
                if (r < (1+sqrt(5))/2)
                {
                  asy.var<- PEcov.or7(r);
                } else {
                  if (r < 1+sqrt(2)/2)
                  {
                    asy.var<- PEcov.or8(r);
                  } else {
                    if (r < 2)
                    {
                      asy.var<- PEcov.or9(r);
                    } else {
                      asy.var<- PEcov.or10(r);
                    }}}}}}}}}}
  asy.var #need to multiply this by 4 in the asymptotic approximation
} #end of the function
#'
#' @rdname funsMuVarUndPE2D
#'
#' @export asy.varAndPE2D
asy.varAndPE2D<-function(r)
{
  if (!pcds::is.point(r,1) || r<1)
  {stop('The argument must be a scalar greater than 1')}

  #internal AND underlying functions for covariance (i.e. asymptotic variance)
  PEcov.and0<-function(r)
  {-(r-1)^2*(972*r^19+8748*r^18+44456*r^17+140328*r^16+121371*r^15-412117*r^14-27145*r^13-4503501*r^12+
               1336147*r^11+10640999*r^10-982009*r^9-6677105*r^8-2274458*r^7-1150162*r^6+249126*r^5+
               1232530*r^4+1234372*r^3+226776*r^2-184944*r-81920)/(58320*(2*r^2+1)*(r+2)^2*(r+1)^3*r^6)}

  PEcov.and1<-function(r)
  {-(486*r^21+3402*r^20-269*r^19-45155*r^18-118850*r^17+443518*r^16+3251855*r^15-13836295*r^14+
       13434672*r^13+11140788*r^12-27667544*r^11+13293088*r^10+7159710*r^9-13013598*r^8+4185440*r^7+
       3262952*r^6+586636*r^5-1616444*r^4-680120*r^3-55952*r^2+219936*r+
       49152)/(116640*(2*r^2+1)*(r+2)^2*(r+1)^3*r^6)}

  PEcov.and2<-function(r)
  {-(486*r^21+3402*r^20-269*r^19-45155*r^18-118850*r^17+443518*r^16+2751855*r^15-13736295*r^14+
       18084672*r^13+8770788*r^12-43009544*r^11+24604048*r^10+27137438*r^9-30889822*r^8-2832544*r^7+
       11101160*r^6-4168820*r^5+2364868*r^4+2305864*r^3-3041936*r^2+219936*r+
       49152)/(116640*(r+2)^2*(2*r^2+1)*(r+1)^3*r^6)}

  PEcov.and3<-function(r)
  {-(3632*r^22+25632*r^21-60328*r^20-441888*r^19+1353430*r^18-297666*r^17-4791125*r^16+12849927*r^15-
       10894618*r^14-26295324*r^13+62283823*r^12-2280753*r^11-81700012*r^10+32551926*r^9+39974410*r^8-
       11284026*r^7-5806580*r^6-9167580*r^5-2004944*r^4+4646688*r^3+1931776*r^2-489024*r-
       98304)/(58320*(r+2)^3*(r^2-2)*(2*r^2+1)*(r+1)^3*r^6)}

  PEcov.and4<-function(r)
  {-(3632*r^22+25632*r^21-49432*r^20-364992*r^19+958940*r^18-1167012*r^17+1200518*r^16+5424126*r^15-
       23566328*r^14+23837088*r^13+11797395*r^12-41623065*r^11+39261953*r^10-8239197*r^9-30178496*r^8+
       27901506*r^7-4936170*r^6+61038*r^5+4719720*r^4-5513952*r^3+340736*r^2+23328*r+
       65536)/(58320*(r+2)^3*(2*r^2+1)*(r^2+1)*(r+1)^3*r^6)}

  PEcov.and5<-function(r)
  {(1562*r^21-11142*r^20-103099*r^19+2105697*r^18-9774118*r^17+10220280*r^16+27825711*r^15-69243129*r^14+
      81624200*r^13-76052574*r^12-65530400*r^11+262451196*r^10-178092280*r^9-69106464*r^8+158439568*r^7-
      97568688*r^6+12246288*r^5+17591952*r^4-21111616*r^3+15628032*r^2-2545664*r+
      993024)/(466560*(r+2)^3*(2*r^2+1)*(r^2+1)*(r+1)^3*r^5)}

  PEcov.and6<-function(r)
  {(1562*r^21-11142*r^20-103099*r^19+2105697*r^18-9774118*r^17+10220280*r^16+27825711*r^15-69243129*r^14+
      81624200*r^13-76052574*r^12-65530400*r^11+262451196*r^10-178092280*r^9-69106464*r^8+158439568*r^7-
      97568688*r^6+12246288*r^5+17591952*r^4-21111616*r^3+15628032*r^2-2545664*r+
      993024)/(466560*(r+2)^3*(2*r^2+1)*(r^2+1)*(r+1)^3*r^5)}

  PEcov.and7<-function(r)
  {-(2*r^26-30*r^25+281*r^24-2395*r^23+8770*r^22+29528*r^21-268053*r^20+245667*r^19+2066216*r^18-
       5313494*r^17-1589216*r^16+18512684*r^15-18946136*r^14-2665248*r^13+22789584*r^12-32987760*r^11+
       20482512*r^10+13109584*r^9-28084416*r^8+17326976*r^7-3864576*r^6-4579328*r^5+6666240*r^4-
       3576320*r^3+635904*r^2-116736*r+61440)/(1920*(2*r^2+1)*(r+2)^3*(r^2+1)*(r+1)^3*r^10)}

  PEcov.and8<-function(r)
  {-(2*r^26-30*r^25+281*r^24-2395*r^23+8258*r^22+31064*r^21-262677*r^20+225443*r^19+2052136*r^18-
       5219030*r^17-1608928*r^16+18337836*r^15-18837080*r^14-2598688*r^13+22736336*r^12-32858736*r^11+
       20384720*r^10+12930896*r^9-27988416*r^8+17416832*r^7-3862784*r^6-4575488*r^5+6638848*r^4-
       3603200*r^3+640512*r^2-107520*r+63488)/(1920*(2*r^2+1)*(r+2)^3*(r+1)^3*(r^2+1)*r^10)}

  PEcov.and9<-function(r)
  {-(2*r^25-32*r^24+307*r^23-2612*r^22+11572*r^21+21934*r^20-328867*r^19+524994*r^18+2446870*r^17-
       8676180*r^16-437020*r^15+36944680*r^14-40677696*r^13-44860384*r^12+106256352*r^11-15515040*r^10-
       98636848*r^9+66358080*r^8+27142272*r^7-42614272*r^6+7781120*r^5+7327232*r^4-3388672*r^3+
       430592*r^2-171008*r+63488)/(1920*(r+2)^3*(r-1)*(r+1)^3*(2*r^2-1)*r^10)}

  PEcov.and10<-function(r)
  {(30*r^13+90*r^12-127*r^11-621*r^10+320*r^9+1568*r^8-858*r^7-1370*r^6+909*r^5+295*r^4-292*r^3+44*r^2+
      6*r-2)/(15*(2*r^2-1)*(r+1)^3*r^10)}
  ###
  asy.var<-0;
  if (r < 2*sqrt(3)/3)
  {
    asy.var<-PEcov.and0(r);
  } else {
    if (r < 6/5)
    {
      asy.var<- PEcov.and1(r);
    } else {
      if (r < sqrt(5)-1)
      {
        asy.var<- PEcov.and2(r);
      } else {
        if (r < (6+2*sqrt(2))/7)
        {
          asy.var<- PEcov.and3(r);
        } else {
          if (r < 4/3)
          {
            asy.var<- PEcov.and4(r);
          } else {
            if (r < (6+sqrt(15))/7)
            {
              asy.var<- PEcov.and5(r);
            } else {
              if (r < 3/2)
              {
                asy.var<- PEcov.and6(r);
              } else {
                if (r < (1+sqrt(5))/2)
                {
                  asy.var<- PEcov.and7(r);
                } else {
                  if (r < 1+sqrt(2)/2)
                  {
                    asy.var<- PEcov.and8(r);
                  } else {
                    if (r < 2)
                    {
                      asy.var<- PEcov.and9(r);
                    } else {
                      asy.var<- PEcov.and10(r);
                    }}}}}}}}}}
  asy.var #need to multiply this by 4 in the asymptotic approximation
} #end of the function
#'
#' @rdname funsMuVarUndPE2D
#'
#' @export asy.var.undPE2D
asy.var.undPE2D <-function(r,ugraph=c("underlying", "reflexivity"))
{
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  asy.var = ifelse(ugraph == "underlying",asy.varOrPE2D(r),asy.varAndPE2D(r))
  asy.var #need to multiply this by 4 in the asymptotic approximation
}

#################################################################

#' @title The indicator for the presence of an edge from a point to another
#' for the underlying or reflexivity graph of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' one triangle case
#'
#' @description Returns \eqn{I(}\code{p1p2} is an edge
#' in the underlying or reflexivity graph of PE-PCDs \eqn{)}
#' for points \code{p1} and \code{p2} in a given triangle.
#'
#' More specifically, when the argument \code{ugraph="underlying"}, it returns
#' the edge indicator for the PE-PCD underlying graph,
#' that is, returns 1 if \code{p2} is
#' in \eqn{N_{PE}(p1,r)} or \code{p1} is in \eqn{N_{PE}(p2,r)},
#' returns 0 otherwise.
#' On the other hand,
#' when \code{ugraph="reflexivity"}, it returns
#' the edge indicator for the PE-PCD reflexivity graph,
#' that is, returns 1 if \code{p2} is
#' in \eqn{N_{PE}(p1,r)} and \code{p1} is in \eqn{N_{PE}(p2,r)},
#' returns 0 otherwise.
#'
#' In both cases PE proximity region is constructed
#' with respect to the triangle \code{tri} and
#' vertex regions are based on the center, \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of \code{tri}
#' or based on the circumcenter of \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#'
#' If \code{p1} and \code{p2} are distinct
#' and either of them are outside \code{tri}, it returns 0,
#' but if they are identical,
#' then it returns 1 regardless of their locations
#' (i.e., it allows loops).
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param p1 A 2D point whose PE proximity region is constructed.
#' @param p2 A 2D point. The function determines
#' whether there is an edge from \code{p1} to \code{p2} or not
#' in the underlying or reflexivity graph of PE-PCDs.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' @param ugraph The type of the graph based on PE-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return Returns 1 if there is an edge between points \code{p1} and \code{p2}
#' in the underlying or reflexivity graph of PE-PCDs
#' in a given triangle \code{tri}, and 0 otherwise.
#'
#' @seealso \code{\link{IedgePEbasic.tri}}, \code{\link{IedgeAStri}},
#' \code{\link{IedgeCStri}} and \code{\link[pcds]{IarcPEtri}}
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
#' M<-as.numeric(pcds::runif.tri(1,Tr)$g)
#'
#' r<-1.5
#' n<-3
#' set.seed(1)
#' Xp<-pcds::runif.tri(n,Tr)$g
#'
#' IedgePEtri(Xp[1,],Xp[2,],Tr,r,M)
#' IedgePEtri(Xp[1,],Xp[2,],Tr,r,M,ugraph = "reflexivity")
#'
#' P1<-as.numeric(pcds::runif.tri(1,Tr)$g)
#' P2<-as.numeric(pcds::runif.tri(1,Tr)$g)
#' IedgePEtri(P1,P2,Tr,r,M)
#' IedgePEtri(P1,P2,Tr,r,M,ugraph="r")
#' #}
#'
#' @export IedgePEtri
IedgePEtri <- function(p1,p2,tri,r,M=c(1,1,1),
                       ugraph=c("underlying", "reflexivity"))
{
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  arc12 = pcds::IarcPEtri(p1,p2,tri,r,M)
  arc21 = pcds::IarcPEtri(p2,p1,tri,r,M)

  edge = ifelse(ugraph == "underlying",
                max(arc12,arc21),
                arc12*arc21)
  edge
} #end of the function
#'

#################################################################

#' @title Number of edges in the underlying or reflexivity graph of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' one triangle case
#'
#' @description
#' An object of class \code{"NumEdges"}.
#' Returns the number of edges of
#' the underlying or reflexivity graph of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs)
#' whose vertices are the
#' given 2D numerical data set, \code{Xp}
#' in a given triangle.
#' It also provides number of vertices
#' (i.e., number of data points inside the triangle)
#' and indices of the data points that reside in the triangle.
#'
#' PE proximity region \eqn{N_{PE}(x,r)} is defined
#' with respect to the triangle, \code{tri}
#' with expansion parameter \eqn{r \ge 1} and vertex regions are
#' based on the center \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri} or
#' based on circumcenter of \code{tri};
#' default is \eqn{M=(1,1,1)},
#' i.e., the center of mass of \code{tri}.
#' For the number of edges, loops are not allowed,
#' so edges are only possible for points
#' inside the triangle \code{tri} for this function.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of PE-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' @param ugraph The type of the graph based on PE-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of edges
#' and quantities related to the triangle}
#' \item{und.graph}{Type of the graph as "Underlying"
#' or "Reflexivity" for the PE-PCD}
#' \item{num.edges}{Number of edges of the underlying
#' or reflexivity graphs based on the PE-PCD
#' with vertices in the given triangle \code{tri}}
#' \item{num.in.tri}{Number of \code{Xp} points in the triangle, \code{tri}}
#' \item{ind.in.tri}{The vector of indices of the \code{Xp} points
#' that reside in the triangle}
#' \item{tess.points}{Tessellation points, i.e., points on which
#' the tessellation of the study region is performed,
#' here, tessellation is the support triangle.}
#' \item{vertices}{Vertices of the underlying or reflexivity graph, \code{Xp}.}
#'
#' @seealso \code{\link{num.edgesPE}}, \code{\link{num.edgesAStri}},
#' \code{\link{num.edgesCStri}}, and \code{\link[pcds]{num.arcsPEtri}}
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
#' Nedges = num.edgesPEtri(Xp,Tr,r=1.25,M)
#' Nedges
#' summary(Nedges)
#' plot(Nedges)
#' #}
#'
#' @export num.edgesPEtri
num.edgesPEtri <- function(Xp,tri,r,M=c(1,1,1),ugraph=c("underlying", "reflexivity"))
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

  if (!pcds::is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

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
        edges.in.tri<-edges.in.tri+IedgePEtri(Xp[i,],Xp[k,],tri,r,M,ugraph=ugraph)
      }
      }

      for (j in (i:n)[-1]) #to avoid loops
      {
        tot.edges<-tot.edges+IedgePEtri(Xp[i,],Xp[j,],tri,r,M,ugraph=ugraph)
      }
    }
  }

  NinTri = length(ind.in.tri)

  und.graph = ifelse(ugraph=="underlying",
                     "Underlying",
                     "Reflexivity")
  desc<-paste("Number of Edges of the ",und.graph,
              " Graphs of the PE-PCD and the Related Quantities with vertices Xp in One Triangle",sep="")

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

#################################################################

#' @title Edge density of the underlying or reflexivity graph of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' one triangle case
#'
#' @description Returns the edge density
#' of the underlying or reflexivity graph of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs)
#' whose vertex set is the given 2D numerical data set, \code{Xp},
#' (some of its members are) in the triangle \code{tri}.
#'
#' PE proximity regions is defined with respect to \code{tri} with
#' expansion parameter \eqn{r \ge 1} and vertex regions are
#' based on center \eqn{M=(m_1,m_2)} in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri} or based on
#' circumcenter of \code{tri}; default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' The function also provides edge density standardized
#' by the mean and asymptotic variance of the edge density
#' of the underlying or reflexivity graph of PE-PCD
#' for uniform data in the triangle \code{tri}
#' only when \code{M} is the center of mass.
#' For the number of edges, loops are not allowed.
#'
#' \code{in.tri.only} is a logical argument (default is \code{FALSE})
#' for considering only the points
#' inside the triangle or all the points as the vertices of the digraph.
#' if \code{in.tri.only=TRUE}, edge density is computed only for
#' the points inside the triangle (i.e., edge density of the subgraph of
#' the underlying or reflexivity graph
#' induced by the vertices in the triangle is computed),
#' otherwise edge density of the entire graph
#' (i.e., graph with all the vertices) is computed.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graph of the PE-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' @param ugraph The type of the graph based on PE-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#' @param in.tri.only A logical argument (default is \code{in.tri.only=FALSE})
#' for computing the edge density for only the points inside the triangle, \code{tri}.
#' That is,
#' if \code{in.tri.only=TRUE} edge density of the induced subgraph with the vertices
#' inside \code{tri} is computed, otherwise
#' otherwise edge density of the entire graph (i.e., graph with all the vertices) is computed.
#'
#' @return A \code{list} with the elements
#' \item{edge.dens}{Edge density of the underlying
#' or reflexivity graphs of the PE-PCD
#' whose vertices are the 2D numerical data set, \code{Xp};
#' PE proximity regions are defined
#' with respect to the triangle \code{tri} and \code{M}-vertex regions}
#' \item{std.edge.dens}{Edge density standardized
#' by the mean and asymptotic variance of the edge
#' density of the underlying or reflexivity graph
#' of the PE-PCD for uniform data in the triangle \code{tri}.
#' This will only be returned, if \code{M} is the center of mass.}
#'
#' @seealso \code{\link{ASedge.dens.tri}}, \code{\link{CSedge.dens.tri}},
#' and \code{\link[pcds]{PEarc.dens.tri}}
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
#' num.edgesPEtri(Xp,Tr,r=1.5,M)$num.edges
#' PEedge.dens.tri(Xp,Tr,r=1.5,M)
#' PEedge.dens.tri(Xp,Tr,r=1.5,M,in.tri.only = TRUE)
#'
#' #For the reflexivity graph
#' num.edgesPEtri(Xp,Tr,r=1.5,M,ugraph="r")$num.edges
#' PEedge.dens.tri(Xp,Tr,r=1.5,M,ugraph="r")
#' PEedge.dens.tri(Xp,Tr,r=1.5,M,in.tri.only = TRUE,ugraph="r")
#' #}
#'
#' @export PEedge.dens.tri
PEedge.dens.tri <- function(Xp,tri,r,M=c(1,1,1),
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

  nedges<-num.edgesPEtri(Xp,tri,r,M,ugraph)$num.edges

  mean.rho<-mu.undPE2D(r,ugraph)
  var.rho<-4*asy.var.undPE2D(r,ugraph)

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
  res=list(edge.dens=rho)

  CM=apply(tri,2,mean)
  if (isTRUE(all.equal(M,CM))){
    std.rho<-sqrt(n)*(rho-mean.rho)/sqrt(var.rho)
    res=list(
      edge.dens=rho, #edge density
      std.edge.dens=std.rho
    )}

  res
} #end of the function
#'

#################################################################

#' @title Number of edges of the underlying or reflexivity graph of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' multiple triangle case
#'
#' @description
#' An object of class \code{"NumEdges"}.
#' Returns the number of edges of
#' the underlying or reflexivity graph of
#' Proportional Edge Proximity Catch Digraph (PE-PCD)
#' and various other quantities and vectors such as
#' the vector of number of vertices (i.e., number of data points)
#' in the Delaunay triangles,
#' number of data points in the convex hull of \code{Yp} points,
#' indices of the Delaunay triangles for the data points, etc.
#'
#' PE proximity regions are defined with respect to the
#' Delaunay triangles based on \code{Yp} points
#' with expansion parameter \eqn{r \ge 1}
#' and vertex regions in each triangle
#' is based on the center \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of each
#' Delaunay triangle or based on circumcenter of each Delaunay triangle
#' (default for \eqn{M=(1,1,1)}
#' which is the center of mass of the triangle).
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
#' for more on PE-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds.ugraph})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graph of the PE-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 3D point in barycentric coordinates
#' which serves as a center in the interior of each Delaunay
#' triangle or circumcenter of each Delaunay triangle
#' (for this, argument should be set as \code{M="CC"}),
#' default for \eqn{M=(1,1,1)}
#' which is the center of mass of each triangle.
#' @param ugraph The type of the graph based on PE-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of edges
#' and related quantities for the induced subgraphs of the underlying
#' or reflexivity graphs (of PE-PCD) in the Delaunay triangles}
#' \item{und.graph}{Type of the graph as "Underlying"
#' or "Reflexivity" for the PE-PCD}
#' \item{num.edges}{Total number of edges in all triangles,
#' i.e., the number of edges for the entire underlying
#' or reflexivity graphs of the PE-PCD}
#' \item{num.in.conv.hull}{Number of \code{Xp} points
#' in the convex hull of \code{Yp} points}
#' \item{num.in.tris}{The vector of number of \code{Xp} points
#' in the Delaunay triangles based on \code{Yp} points}
#' \item{weight.vec}{The \code{vector} of the areas of
#' Delaunay triangles based on \code{Yp} points}
#' \item{tri.num.edges}{The \code{vector} of the number of edges
#' of the components of the PE-PCD in the
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
#' @seealso \code{\link{num.edgesPEtri}}, \code{\link{num.edgesAS}},
#' \code{\link{num.edgesCS}}, and \code{\link[pcds]{num.arcsPE}}
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
#' Nedges = num.edgesPE(Xp,Yp,r=1.5,M)
#' Nedges
#' summary(Nedges)
#' plot(Nedges)
#' #}
#'
#' @export num.edgesPE
num.edgesPE <- function(Xp,Yp,r,M=c(1,1,1),ugraph=c("underlying", "reflexivity"))
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

  if (!pcds::is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

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
           Tri<-Yp[Ytri[i,],])
    #vertices of ith triangle
    Wvec<-c(Wvec,pcds::area.polygon(Tri))
  }

  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  und.graph = ifelse(ugraph=="underlying",
                     "Underlying",
                     "Reflexivity")

  if (ny==3)
  { #tri<-pcds::as.basic.tri(Yp)$tri
    #NumEdges = num.edgesPEtri(Xp,tri,r,M,ugraph)
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
      NumEdges = num.edgesPEtri(Xdt,tri,r,M,ugraph) #for the vertices inside the triangle
      #Wvec<-pcds::area.polygon(tri)
      Tot.Edges<- NumEdges$num.edges
      #number of edges in the triangle Yp
      ni.vec = NumEdges$num.in.tri
      Tri.Ind = NumEdges$ind.in.tri #returns 1's if the points Xp[i,]'s are inside triangle based on Yp, NA otherwise
      data.tri.ind = rep(NA,nx)
      data.tri.ind[Tri.Ind] = 1
      edges = NumEdges$num.edges
      ind.in.CH = which(inCH) #which(!is.na(Tri.Ind))
    }

    Tot.Edges = Tot.Edges + sum(duplicated(Xp[!inCH,]))

    desc<-paste("Number of Edges of the ",und.graph,
                " Graphs of the PE-PCD with vertices Xp and the Related Quantities for the Induced Subdigraph for the Points in the Delaunay Triangle",sep="")

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
    {Tot.Edges<-0;
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
        #  data.del.tris=c(data.del.tris,list(Xpi))
        ifelse(ndt==1,
               Tri<-Yp[Ytri,],
               Tri<-Yp[Ytri[i,],])
        #vertices of ith triangle
        tri<-pcds::as.basic.tri(Tri)$tri
        #convert the triangle Tri into an nonscaled basic triangle, see as.basic.tri help page
        #  del.tris=rbind(del.tris,tri)
        ni.vec<-c(ni.vec,length(Xpi)/2)
        #number of points in ith Delaunay triangle

        ifelse(identical(M,"CC"),
               cent<-pcds::circumcenter.tri(tri),
               cent<-M)
        num.edges<-num.edgesPEtri(Xpi,tri,r,cent,ugraph)$num.edges
        #number of edges in ith triangle
        edges<-c(edges,num.edges)
        #number of edges in all triangles as A \code{vector}
      }

      Tot.Edges<-sum(edges)  #the total number of edges in all triangles
    }

    Tot.Edges = Tot.Edges + sum(duplicated(Xp[!inCH,]))

    desc<-paste("Number of Edges of the ",und.graph,
                " Graphs of the PE-PCD with vertices Xp and the Related Quantities for the Induced Subdigraphs for the Points in the Delaunay Triangles",sep="")

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
              # indices of the Delaunay triangles, each column corresponds to the vector of indices
              #of the vertices of one triangle
              data.tri.ind=data.tri.ind,
              #indices of the Delaunay triangles in which data points reside,
              #i.e., column number of del.tri.ind for each Xp point
              tess.points=Yp, #tessellation points
              vertices=Xp #vertices of the digraph
    )
  }
  class(res) <- "NumEdges"
  res$call <-match.call()

  res
} #end of the function
#'

################################################################# bura

#' @title A test of segregation/association based on edge density of underlying
#' or reflexivity graph of Proportional Edge Proximity Catch Digraph
#' (PE-PCD) for 2D data
#'
#' @description
#' An object of class \code{"htest"} (i.e., hypothesis test) function
#' which performs a hypothesis test of complete spatial
#' randomness (CSR) or uniformity of \code{Xp} points
#' in the convex hull of \code{Yp} points against the alternatives
#' of segregation (where \code{Xp} points cluster
#' away from \code{Yp} points) and association
#' (where \code{Xp} points cluster around
#' \code{Yp} points) based on the normal approximation
#' of the edge density of the underlying or reflexivity graph of
#' PE-PCD for uniform 2D data.
#'
#' The function yields the test statistic,
#' \eqn{p}-value for the corresponding \code{alternative},
#' the confidence interval,
#' estimate and null value for the parameter of interest
#' (which is the edge density),
#' and method and name of the data set used.
#'
#' Under the null hypothesis of uniformity of \code{Xp} points
#' in the convex hull of \code{Yp} points, edge density
#' of underlying or reflexivity graph of PE-PCD
#' whose vertices are \code{Xp} points equals
#' to its expected value under the uniform distribution and
#' \code{alternative} could be two-sided, or left-sided
#' (i.e., data is accumulated around the \code{Yp} points, or association)
#' or right-sided (i.e., data is accumulated
#' around the centers of the triangles,
#' or segregation).
#'
#' PE proximity region is constructed
#' with the expansion parameter \eqn{r \ge 1} and \eqn{CM}-vertex regions
#' (i.e., the test is not available for a general center \eqn{M}
#' at this version of the function).
#'
#' **Caveat:** This test is currently a conditional test,
#' where \code{Xp} points are assumed to be random,
#' while \code{Yp} points are
#' assumed to be fixed (i.e., the test is conditional on \code{Yp} points).
#' Furthermore,
#' the test is a large sample test when \code{Xp} points
#' are substantially larger than \code{Yp} points,
#' say at least 5 times more.
#' This test is more appropriate when supports of \code{Xp}
#' and \code{Yp} have a substantial overlap.
#' Currently, the \code{Xp} points
#' outside the convex hull of \code{Yp} points
#' are handled with a correction factor
#' which is derived under the assumption of
#' uniformity of \code{Xp} and \code{Yp} points in the study window,
#' (see the description below for the argument \code{ch.cor} and the function code.)
#' However, in the special case of no \code{Xp} points
#' in the convex hull of \code{Yp} points,
#' edge density is taken to be 1,
#' as this is clearly a case of segregation.
#' Removing the conditioning and extending it to
#' the case of non-concurring supports are
#' topics of ongoing research of the author of the package.
#'
#' \code{ch.cor} is for convex hull correction
#' (default is \code{"no convex hull correction"}, i.e., \code{ch.cor=FALSE})
#' which is recommended
#' when both \code{Xp} and \code{Yp} have the same rectangular support.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph})
#' for more on the test based on the edge density of
#' underlying or reflexivity graph of PE-PCDs.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graphs of the PE-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param ugraph The type of the graph based on PE-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#' @param ch.cor A logical argument for convex hull correction,
#' default \code{ch.cor=FALSE},
#' recommended when both \code{Xp} and \code{Yp}
#' have the same rectangular support.
#' @param alternative Type of the alternative hypothesis in the test,
#' one of \code{"two.sided"}, \code{"less"}, \code{"greater"}.
#' @param conf.level Level of the confidence interval,
#' default is \code{0.95}, for the edge density of
#' underlying or reflexivity graph of PE-PCD based on
#' the 2D data set \code{Xp}.
#'
#' @return A \code{list} with the elements
#' \item{statistic}{Test statistic}
#' \item{p.value}{The \eqn{p}-value for the hypothesis test
#' for the corresponding \code{alternative}}
#' \item{conf.int}{Confidence interval for the edge density
#' at the given confidence level \code{conf.level} and
#' depends on the type of \code{alternative}.}
#' \item{estimate}{Estimate of the parameter, i.e., edge density}
#' \item{null.value}{Hypothesized value for the parameter,
#' i.e., the null edge density, which is usually the
#' mean edge density under uniform distribution.}
#' \item{alternative}{Type of the alternative hypothesis in the test,
#' one of \code{"two.sided"}, \code{"less"}, \code{"greater"}}
#' \item{method}{Description of the hypothesis test}
#' \item{data.name}{Name of the data set}
#'
#' @seealso \code{\link{CSedge.dens.test}} and \code{\link[pcds]{PEarc.dens.test}}
#'
#' @references
#' \insertAllCited{}
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' \donttest{
#' #nx is number of X points (target) and ny is number of Y points (nontarget)
#' nx<-100; ny<-5;
#'
#' set.seed(1)
#' Xp<-cbind(runif(nx),runif(nx))
#' Yp<-cbind(runif(ny,0,.25),
#' runif(ny,0,.25))+cbind(c(0,0,0.5,1,1),c(0,1,.5,0,1))
#'
#' pcds::plotDelaunay.tri(Xp,Yp,xlab="",ylab="")
#'
#' PEedge.dens.test(Xp,Yp,r=1.25)
#' PEedge.dens.test(Xp,Yp,r=1.25,ch=TRUE)
#'
#' PEedge.dens.test(Xp,Yp,r=1.25,ugraph="r")
#' PEedge.dens.test(Xp,Yp,r=1.25,ugraph="r",ch=TRUE)
#' #since Y points are not uniform, convex hull correction is invalid here
#' }
#'
#' @export PEedge.dens.test
PEedge.dens.test <- function(Xp,Yp,r,ugraph=c("underlying", "reflexivity"),ch.cor=FALSE,
                             alternative = c("two.sided", "less", "greater"),
                             conf.level = 0.95)
{
  dname <-deparse(substitute(Xp))

  alternative <-match.arg(alternative)
  if (length(alternative) > 1 || is.na(alternative))
    stop("alternative must be one of \"greater\", \"less\", \"two.sided\"")

  if (!is.numeric(as.matrix(Xp)) || !is.numeric(as.matrix(Yp)))
  {stop('Xp and Yp must be numeric')}

  if (pcds::is.point(Xp))
  { Xp<-matrix(Xp,ncol=2)
  } else
  {Xp<-as.matrix(Xp)
  if (ncol(Xp)!=2 )
  {stop('Xp must be of dimension nx2')}
  }

  n<-nrow(Xp)  #number of X points
  if (n<=1)
  {stop('The graph is void or has only one vertex!
    So, there are not enough Xp points to compute the edge density!')}

  Yp<-as.matrix(Yp)
  if (ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('Yp must be of dimension kx2 with k>=3')}

  if (!pcds::is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) ||
        conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")

  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  Edges<-num.edgesPE(Xp,Yp,r,M=c(1,1,1),ugraph)
  #use the default, i.e., CM for the center M
  NinCH<-Edges$num.in.con

  num.edges<-Edges$num.edges #total number of edges in the PE-PCD
  num.edges.tris = Edges$tri.num.edges
  #vector of number of edges in the Delaunay triangles
  num.dat.tris = Edges$num.in.tris
  #vector of number of data points in the Delaunay triangles
  Wvec<-Edges$w
  LW<-Wvec/sum(Wvec)

  tri.ind = Edges$data.tri.ind
  ind.triCH =  t(Edges$del.tri)

  ind.Xp1 = which(num.dat.tris==1)
  if (length(ind.Xp1)>0)
  {
    for (i in ind.Xp1)
    {
      Xpi = Xp[which(tri.ind==i),]
      tri =  Yp[ind.triCH[i,],]
      npe = pcds::NPEtri(Xpi,tri,r)
      num.edges = num.edges+pcds::area.polygon(npe)/Wvec[i]
    }
  }
  asy.mean0<-mu.undPE2D(r,ugraph)  #asy mean value for the r value
  asy.mean<-asy.mean0*sum(LW^2)

  asy.var0<-4*asy.var.undPE2D(r,ugraph)  #asy variance value for the r value
  asy.var<-asy.var0*sum(LW^3)+4*asy.mean0^2*(sum(LW^3)-(sum(LW^2))^2)

  if (NinCH == 0) {
    warning('There is no Xp point in the convex hull of Yp points to compute edge density,
           but as this is clearly a segregation pattern, so edge density is taken to be 1!')
    edge.dens=1
    TS0<-sqrt(n)*(edge.dens-asy.mean)/sqrt(asy.var)  #standardized test stat
  } else
  {  edge.dens<-num.edges/(NinCH*(NinCH-1)/2)
  TS0<-sqrt(NinCH)*(edge.dens-asy.mean)/sqrt(asy.var)
  #standardized test stat}  #edge density
  }
  estimate1<-edge.dens; estimate2<-asy.mean

  method <- c("Large Sample z-Test Based on Edge Density of", ifelse(ugraph ==  "underlying", "underlying","reflexivity"), "graph of PE-PCD for Testing Uniformity of 2D Data ---")

  if (ch.cor==FALSE)
  {
    TS<-TS0
    method <-c(method, " without Convex Hull Correction")
  }
  else
  {
    m<-nrow(Yp)  #number of Y points
    NoutCH<-n-NinCH #number of points outside of the convex hull

    prop.out<-NoutCH/n #observed proportion of points outside convex hull
    exp.prop.out<-1.7932/m+1.2229/sqrt(m)
    #expected proportion of points outside convex hull

    TS<-TS0+abs(TS0)*sign(prop.out-exp.prop.out)*(prop.out-exp.prop.out)^2
    method <-c(method, " with Convex Hull Correction")
  }

  names(estimate1) <-c("edge density")
  null.dens<-asy.mean
  names(null.dens) <-"(expected) edge density"
  names(TS) <-"standardized edge density (i.e., Z)"

  if (alternative == "less") {
    pval <-pnorm(TS)
    cint <-edge.dens+c(-Inf, qnorm(conf.level))*sqrt(asy.var/NinCH)
  }
  else if (alternative == "greater") {
    pval <-pnorm(TS, lower.tail = FALSE)
    cint <-edge.dens+c(-qnorm(conf.level),Inf)*sqrt(asy.var/NinCH)
  }
  else {
    pval <-2 * pnorm(-abs(TS))
    alpha <-1 - conf.level
    cint <-qnorm(1 - alpha/2)
    cint <-edge.dens+c(-cint, cint)*sqrt(asy.var/NinCH)
  }

  attr(cint, "conf.level") <- conf.level

  rval <-list(
    statistic=TS,
    p.value=pval,
    conf.int = cint,
    estimate = estimate1,
    null.value = null.dens,
    alternative = alternative,
    method = method,
    data.name = dname
  )

  attr(rval, "class") <-"htest"
  return(rval)
} #end of the function
#'

#################################################################

#' @title Incidence matrix for the underlying or reflexivity graph of
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' one triangle case
#'
#' @description Returns the incidence matrix
#' for the underlying or reflexivity graph of the PE-PCD
#' whose vertices are the given 2D numerical data set, \code{Xp},
#' in the triangle \code{tri}\eqn{=T(v=1,v=2,v=3)}.
#'
#' PE proximity regions are constructed with respect to triangle \code{tri}
#' with expansion parameter \eqn{r \ge 1}
#' and vertex regions are based on the center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' Loops are allowed, so the diagonal entries are all equal to 1.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graph of the PE-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' @param ugraph The type of the graph based on PE-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return Incidence matrix for the underlying or reflexivity graph
#' of the PE-PCD with vertices
#' being 2D data set, \code{Xp}
#' in the triangle \code{tri} with vertex regions based on center \code{M}
#'
#' @seealso \code{\link{inci.mat.undPE}}, \code{\link{inci.mat.undAStri}},
#' \code{\link{inci.mat.undCStri}}, and \code{\link[pcds]{inci.matPEtri}}
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
#' (IM<-inci.mat.undPEtri(Xp,Tr,r=1.25,M))
#' pcds::dom.num.greedy(IM)
#' pcds::Idom.num.up.bnd(IM,3)
#'
#' (IM<-inci.mat.undPEtri(Xp,Tr,r=1.25,M,ugraph="r"))
#' pcds::dom.num.greedy(IM)
#' pcds::Idom.num.up.bnd(IM,3)
#' #}
#'
#' @export inci.mat.undPEtri
inci.mat.undPEtri <- function(Xp,tri,r,M=c(1,1,1),ugraph=c("underlying", "reflexivity"))
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

  if (!pcds::is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

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
  inc.mat<-matrix(0, nrow=n, ncol=n)
  # if (n>=1)
  # {
  for (i in 1:n)
  {p1<-Xp[i,]

  for (j in (i:n) )
  {p2<-Xp[j,]
  inc.mat[i,j]<-inc.mat[j,i]<-IedgePEtri(p1,p2,tri,r,M,ugraph)
  }
  }
  # }
  inc.mat
} #end of the function
#'

#################################################################

#' @title The edges of the underlying or reflexivity graph of
#' the Proportional Edge Proximity Catch Digraph
#' (PE-PCD) for 2D data - one triangle case
#'
#' @description
#' An object of class \code{"UndPCDs"}.
#' Returns edges of the underlying or reflexivity graph of PE-PCD
#' as left and right end points
#' and related parameters and the quantities of these graphs.
#' The vertices of these graphs are the data points in \code{Xp}
#' in the multiple triangle case.
#'
#' PE proximity regions are constructed
#' with respect to the triangle \code{tri} with expansion
#' parameter \eqn{r \ge 1}, i.e.,
#' edges may exist only for points inside \code{tri}.
#' It also provides various descriptions
#' and quantities about the edges of
#' the underlying or reflexivity graph of the PE-PCD
#' such as number of edges, edge density, etc.
#'
#' Vertex regions are based on center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of
#' the triangle \code{tri} or based on the circumcenter of \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' When the center is the circumcenter, \code{CC},
#' the vertex regions are constructed based on the
#' orthogonal projections to the edges,
#' while with any interior center \code{M},
#' the vertex regions are constructed using the extensions
#' of the lines combining vertices with \code{M}.
#' The different consideration of circumcenter vs
#' any other interior center of the triangle is because
#' the projections from circumcenter are orthogonal to the edges,
#' while projections of \code{M} on the edges are on the extensions
#' of the lines connecting \code{M} and the vertices.
#' \code{M}-vertex regions are recommended spatial inference,
#' due to geometry invariance property of the edge density
#' and domination number the PE-PCDs based on uniform data.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graph of the PE-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' @param ugraph The type of the graph based on PE-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the underlying
#' or reflexivity graph of the digraph}
#' \item{parameters}{Parameters of the underlying
#' or reflexivity graph of the digraph,
#' the center \code{M} used to
#' construct the vertex regions and the expansion parameter \code{r}.}
#' \item{tess.points}{Tessellation points, i.e., points on which
#' the tessellation of the study region
#' is performed, here, tessellation is the support triangle.}
#' \item{tess.name}{Name of the tessellation points \code{tess.points}}
#' \item{vertices}{Vertices of the underlying
#' or reflexivity graph of the digraph, \code{Xp} points}
#' \item{vert.name}{Name of the data set
#' which constitutes the vertices of
#' the underlying or reflexivity graph of the digraph}
#' \item{LE, RE}{Left and right end points of the edges of
#' the underlying or reflexivity graph of PE-PCD for 2D data set \code{Xp}
#' as vertices of the underlying or reflexivity graph of the digraph}
#' \item{mtitle}{Text for \code{"main"} title
#' in the plot of the underlying or reflexivity graph of the digraph}
#' \item{quant}{Various quantities for the underlying
#' or reflexivity graph of the digraph:
#' number of vertices, number of partition points,
#' number of intervals, number of edges, and edge density.}
#'
#' @seealso \code{\link{edgesPE}}, \code{\link{edgesAStri}}, \code{\link{edgesCStri}},
#' and \code{\link[pcds]{arcsPEtri}}
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
#' r<-1.5
#'
#' #for underlying graph
#' Edges<-edgesPEtri(Xp,Tr,r,M)
#' Edges
#' summary(Edges)
#' plot(Edges)
#'
#' #for reflexivity graph
#' Edges<-edgesPEtri(Xp,Tr,r,M,ugraph="r")
#' Edges
#' summary(Edges)
#' plot(Edges)
#'
#' #can add vertex regions
#' #but we first need to determine center is the circumcenter or not,
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
#' @export edgesPEtri
edgesPEtri <- function(Xp,tri,r,M=c(1,1,1),ugraph=c("underlying", "reflexivity"))
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

  if (!pcds::is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

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

  #the edges of the underlying or reflexivity graph of PE-PCDs
  lep<-rep<-NULL #left and right end points for the edges
  if (n2>1)
  {
    for (j in 1:(n2-1))
    {
      p1<-Xtri[j,];
      for (k in (j+1):n2)  #to avoid loops
      {
        p2<-Xtri[k,];
        if (IedgePEtri(p1,p2,tri,r,M,ugraph)==1)
        {
          lep <-rbind(lep,p1); rep <-rbind(rep,p2);
        }
      }
    }
  }

  param<-list(M,r)
  Mr<-round(M,2)

  if (identical(M,"CC") || isTRUE(all.equal(M,CC)))
  {
    cname <-"CC"
    names(param)<-c("circumcenter","expansion parameter")
    main.txt<-paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"), " Graph of PE-PCD with r = ",r," and Circumcenter",sep="")
    typ<-paste(ifelse(ugraph=="underlying","Underlying", "Reflexivity")," Graph of Proportional Edge Proximity Catch Digraph (PE-PCD) for 2D Points in the Triangle with Expansion Parameter r = ",r," and Circumcenter",sep="")
  } else
  {
    cname <-"M"
    names(param)<-c("center","expansion parameter")
    main.txt<-paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"), " Graph of PE-PCD with r = ",r," and Center ", cname," = (",Mr[1],",",Mr[2],")",sep="")
    typ<-paste(ifelse(ugraph=="underlying","Underlying", "Reflexivity")," Graph of Proportional Edge Proximity Catch Digraph (PE-PCD) for 2D Points in the Triangle with Expansion Parameter r = ",r," and Center ", cname," = (",Mr[1],",",Mr[2],")",sep="")
  }

  nvert<-n2; ny<-3; ntri<-1;
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
    tess.points=tri, tess.name=yname, #tessellation points
    vertices=Xp, vert.name=xname,
    #vertices of the underlying or reflexivity graph of the digraph
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
#' the Proportional Edge Proximity Catch Digraph
#' (PE-PCD) for 2D data - one triangle case
#'
#' @description Plots the edges of the underlying or reflexivity graph of
#' the Proportional Edge Proximity Catch Digraph
#' (PE-PCD) whose vertices are the data points, \code{Xp}
#' and the triangle \code{tri}.
#' PE proximity regions
#' are constructed with respect to the triangle \code{tri}
#' with expansion parameter \eqn{r \ge 1},
#' i.e., edges may exist only for \code{Xp} points
#' inside the triangle \code{tri}.
#'
#' Vertex regions are based on center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri}
#' or based on the circumcenter of \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' When the center is the circumcenter, \code{CC},
#' the vertex regions are constructed based on the
#' orthogonal projections to the edges,
#' while with any interior center \code{M},
#' the vertex regions are constructed using the extensions
#' of the lines combining vertices with \code{M}.
#' \code{M}-vertex regions are recommended spatial inference,
#' due to geometry invariance property of the edge density
#' and domination number the PE-PCDs based on uniform data.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graphs of the PE-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri}
#' or the circumcenter of \code{tri}
#' which may be entered as "CC" as well;
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' @param ugraph The type of the graph based on PE-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
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
#' or reflexivity graphs of the PE-PCD
#' whose vertices are the points in data set \code{Xp}
#' and the triangle \code{tri}
#'
#' @seealso \code{\link{plotPEedges}}, \code{\link{plotASedges.tri}},
#' \code{\link{plotCSedges.tri}}, and \code{\link[pcds]{plotPEarcs.tri}}
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
#' r<-1.5
#' plotPEedges.tri(Xp,Tr,r,M,vert.reg = TRUE,xlab="",ylab="")
#' plotPEedges.tri(Xp,Tr,r,M,ugraph="r",vert.reg = TRUE,xlab="",ylab="")
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
#' @export plotPEedges.tri
plotPEedges.tri <- function(Xp,tri,r,M=c(1,1,1),ugraph=c("underlying", "reflexivity"),
                            asp=NA,main=NULL,xlab=NULL,ylab=NULL,
                            xlim=NULL,ylim=NULL,vert.reg=FALSE,...)
{
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  if( any(duplicated(as.data.frame(Xp))) )
    #if there are duplicates for Xp values, only one is taken for each
  {Xp = unique(as.data.frame(Xp))
  warning("There were duplicate Xp values;
          only one value is kept for each duplicate Xp value (to avoid edges of zero length)!")}

  EdgesPE<-edgesPEtri(Xp,tri,r,M,ugraph)
  lep<-EdgesPE$LE; rep<-EdgesPE$RE
  #lep, rep are left and right end points of the edges of the graph
  cent = (EdgesPE$param)$c

  Xp<-matrix(Xp,ncol=2)
  if (is.null(xlim))
  {xlim<-range(tri[,1],Xp[,1],cent[1])}
  if (is.null(ylim))
  {ylim<-range(tri[,2],Xp[,2],cent[2])}

  if ( isTRUE(all.equal(cent,pcds::circumcenter.tri(tri))) )
  {M="CC"}

  if (is.null(main))
  {if (identical(M,"CC")){
    main=paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"), " Graph of PE-PCD\n with r = ",r," and Circumcenter",sep="")
  } else {Mr=round(cent,2)
  Mvec= paste(Mr, collapse=",")
  main=paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"), " Graph of PE-PCD\n with r = ",r," and M = (",Mvec,")",sep="")}
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
#' the Proportional Edge Proximity Catch Digraph
#' (PE-PCD) for 2D data - multiple triangle case
#'
#' @description
#' An object of class \code{"UndPCDs"}.
#' Returns edges of the underlying or reflexivity graph of PE-PCD
#' as left and right end points
#' and related parameters and the quantities of these graphs.
#' The vertices of these graphs are the data points in \code{Xp}
#' in the multiple triangle case.
#'
#' PE proximity regions are defined
#' with respect to the Delaunay triangles
#' based on \code{Yp} points with expansion parameter \eqn{r \ge 1} and
#' vertex regions in each triangle are
#' based on the center \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates
#' in the interior of each Delaunay triangle or
#' based on circumcenter of each Delaunay triangle
#' (default for \eqn{M=(1,1,1)}
#' which is the center of mass of the triangle).
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
#' for more on the PE-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds.ugraph})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graph of the PE-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 3D point in barycentric coordinates
#' which serves as a center in the interior of each Delaunay
#' triangle or circumcenter of each Delaunay triangle
#' (for this, argument should be set as \code{M="CC"}),
#' default for \eqn{M=(1,1,1)}
#' which is the center of mass of each triangle.
#' @param ugraph The type of the graph based on PE-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the underlying
#' or reflexivity graph of the digraph}
#' \item{parameters}{Parameters of the underlying
#' or reflexivity graph of the digraph,
#' the center \code{M} used to
#' construct the vertex regions and the expansion parameter \code{r}.}
#' \item{tess.points}{Tessellation points, i.e., points on which the tessellation
#' of the study region is performed, here, tessellation
#' is Delaunay triangulation based on \code{Yp} points.}
#' \item{tess.name}{Name of the tessellation points \code{tess.points}}
#' \item{vertices}{Vertices of the underlying
#' or reflexivity graph of the digraph, \code{Xp} points}
#' \item{vert.name}{Name of the data set
#' which constitute the vertices of the underlying
#' or reflexivity graph of the digraph}
#' \item{LE, RE}{Left and right end points of the edges of
#' the underlying or reflexivity graph of PE-PCD for 2D data set \code{Xp}
#' as vertices of the underlying or reflexivity graph of the digraph}
#' \item{mtitle}{Text for \code{"main"} title
#' in the plot of the underlying or reflexivity graph of the digraph}
#' \item{quant}{Various quantities for the underlying
#' or reflexivity graph of the digraph:
#' number of vertices, number of partition points,
#' number of intervals, number of edges, and edge density.}
#'
#' @seealso \code{\link{edgesPEtri}}, \code{\link{edgesAS}}, \code{\link{edgesCS}},
#' and \code{\link[pcds]{arcsPE}}
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
#' r<-1.5
#'
#' Edges<-edgesPE(Xp,Yp,r,M)
#' Edges
#' summary(Edges)
#' plot(Edges)
#' #}
#'
#' @export edgesPE
edgesPE <- function(Xp,Yp,r,M=c(1,1,1),ugraph=c("underlying", "reflexivity"))
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
  {stop('Xp must be of dimension nx2')}
  }

  Yp<-as.matrix(Yp)
  if (ncol(Yp)!=2 || nrow(Yp)<3)
  {stop('Yp must be of dimension kx2 with k>=3')}

  if (!pcds::is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  if (nrow(Yp)==3)
  {
    res<-edgesPEtri(Xp,Yp,r,M,ugraph)
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

    DTr<-matrix(interp::triangles(DTmesh)[,1:3],ncol=3)
    nt<-nrow(DTr)
    nx2<-nrow(Xch)

    lep<-rep<-NULL #left and right end points for the edges
    if (nx2>1)
    {
      i.tr<-rep(0,nx2)
      #the vector of indices for the triangles that contain the Xp points
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
            if (IedgePEtri(Xl[j,],Xl[k,],Yi.tri,r,cent,ugraph)==1 )
            {
              lep <-rbind(lep,Xl[j,]); rep <-rbind(rep,Xl[k,]);
            }
          }
        }
      }
    }

    cname <-ifelse(identical(M,"CC"),"CC","M")
    param<-list(M,r)
    names(param)<-c("center","expansion parameter")

    if (identical(M,"CC")){
      main.txt=paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"), " Graph of PE-PCD\n with r = ",r," and Circumcenter",sep="")
      typ<-paste(ifelse(ugraph=="underlying","Underlying", "Reflexivity")," Graph of Proportional Edge Proximity Catch Digraph (PE-PCD) for 2D points in Multiple Triangles with Expansion Parameter r = ",r," and Circumcenter",sep="")
    } else {Mvec= paste(M, collapse=",")
    main.txt=paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"), " Graph of PE-PCD\n with r = ",r," and Center ", cname," = (",Mvec,")",sep="")
    typ<-paste(ifelse(ugraph=="underlying","Underlying", "Reflexivity")," Graph of Proportional Edge Proximity Catch Digraph (PE-PCD) for 2D points in Multiple Triangles with Expansion parameter r = ",r," and Center ", cname," = (",Mvec,")",sep="")}

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
      vertices=Xp, vert.name=xname,
      #vertices of the underlying or reflexivity graph of the digraph
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
#' Proportional Edge Proximity Catch Digraphs (PE-PCDs) -
#' multiple triangle case
#'
#' @description Returns the incidence matrix
#' for the underlying or reflexivity graph of the PE-PCD
#' whose vertices are the data points in \code{Xp}
#' in the multiple triangle case.
#'
#' PE proximity regions are
#' defined with respect to the Delaunay triangles
#' based on \code{Yp} points with expansion parameter \eqn{r \ge 1} and
#' vertex regions in each triangle are
#' based on the center \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates
#' in the interior of each Delaunay triangle
#' or based on circumcenter of each Delaunay triangle
#' (default for \eqn{M=(1,1,1)}
#' which is the center of mass of the triangle).
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
#' for more on the PE-PCDs.
#' Also, see
#' (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds.ugraph})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graph of the PE-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 3D point in barycentric coordinates
#' which serves as a center in the interior of each Delaunay
#' triangle or circumcenter of each Delaunay triangle
#' (for this, argument should be set as \code{M="CC"}),
#' default for \eqn{M=(1,1,1)}
#' which is the center of mass of each triangle.
#' @param ugraph The type of the graph based on PE-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return Incidence matrix for the underlying or reflexivity graph
#' of the PE-PCD whose vertices are the 2D data set, \code{Xp}.
#' PE proximity regions are constructed
#' with respect to the Delaunay triangles and \code{M}-vertex regions.
#'
#' @seealso \code{\link{inci.mat.undPEtri}}, \code{\link{inci.mat.undAS}},
#' \code{\link{inci.mat.undCS}}, and \code{\link[pcds]{inci.matPE}}
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
#' r<-1.5
#'
#' IM<-inci.mat.undPE(Xp,Yp,r,M)
#' IM
#' pcds::dom.num.greedy(IM)
#' #}
#'
#' @export inci.mat.undPE
inci.mat.undPE <- function(Xp,Yp,r,M=c(1,1,1),ugraph=c("underlying", "reflexivity"))
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

  if (!pcds::is.point(r,1) || r<1)
  {stop('r must be a scalar >= 1')}

  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  if (nrow(Yp)==3)
  {
    inc.mat<-inci.mat.undPEtri(Xp,Yp,r,M,ugraph)
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
    inc.mat[i,j]<-inc.mat[j,i]<-IedgePEtri(p1,p2,Yi.tri,r,cent,ugraph)
    }
    }
  }
  inc.mat
} #end of the function
#'

#################################################################

#' @title The plot of the edges of the underlying or reflexivity graph of
#' the Proportional Edge Proximity Catch Digraph
#' (PE-PCD) for 2D data - multiple triangle case
#'
#' @description Plots the edges of the underlying or reflexivity graph of
#' the Proportional Edge Proximity Catch Digraph
#' (PE-PCD) whose vertices are the data
#' points in \code{Xp} in the multiple triangle case
#' and the Delaunay triangles based on \code{Yp} points.
#'
#' PE proximity regions are constructed
#' with respect to the Delaunay triangles based on \code{Yp} points, i.e.,
#' PE proximity regions are defined only for \code{Xp} points
#' inside the convex hull of \code{Yp} points.
#' That is, edges may exist for \code{Xp} points
#' only inside the convex hull of \code{Yp} points.
#'
#' Vertex regions in each triangle are
#' based on the center \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of each Delaunay triangle
#' or based on circumcenter of
#' each Delaunay triangle (default for \eqn{M=(1,1,1)}
#' which is the center of mass of the triangle).
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
#' for more on the PE-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds.ugraph})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graphs of the PE-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param r A positive real number
#' which serves as the expansion parameter in PE proximity region;
#' must be \eqn{\ge 1}.
#' @param M A 3D point in barycentric coordinates
#' which serves as a center in the interior of each Delaunay
#' triangle or circumcenter of each Delaunay triangle
#' (for this, argument should be set as \code{M="CC"}),
#' default for \eqn{M=(1,1,1)}
#' which is the center of mass of each triangle.
#' @param ugraph The type of the graph based on PE-PCDs,
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
#' or reflexivity graphs of the PE-PCD
#' whose vertices are the points in data set \code{Xp} and the Delaunay
#' triangles based on \code{Yp} points
#'
#' @seealso \code{\link{plotPEedges.tri}}, \code{\link{plotASedges}},
#' \code{\link{plotCSedges}}, and \code{\link[pcds]{plotPEarcs}}
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
#' r<-1.5
#'
#' plotPEedges(Xp,Yp,r,M,xlab="",ylab="")
#' plotPEedges(Xp,Yp,r,M,xlab="",ylab="",ugraph="r")
#' #}
#'
#' @export plotPEedges
plotPEedges <- function(Xp,Yp,r,M=c(1,1,1),ugraph=c("underlying", "reflexivity"),
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
    plotPEedges.tri(Xp,Yp,r,M,ugraph,asp,main,xlab,ylab,xlim,ylim)
  } else
  {
    EdgesPE<-edgesPE(Xp,Yp,r,M,ugraph)
    lep<-EdgesPE$LE; rep<-EdgesPE$RE

    DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
    Xch<-pcds::Xin.convex.hullY(Xp,Yp)

    if (is.null(main))
    {if (identical(M,"CC")){
      main=paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"), " Graph of PE-PCD\n with r = ",r," and Circumcenter",sep="")
    } else {Mvec= paste(M, collapse=",")
    main=paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"), " Graph of PE-PCD\n with r = ",r," and M = (",Mvec,")",sep="")}
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

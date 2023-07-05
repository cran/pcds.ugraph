#CentSimUG2D.R;
#Functions for Underlying or Reflexivity Graphs of CS-PCD in R^2
#################################################################

#' @title The indicator for the presence of an edge from a point to another
#' for the underlying or reflexivity graphs of
#' Central Similarity Proximity Catch Digraphs (CS-PCDs) -
#' standard basic triangle case
#'
#' @description Returns \eqn{I(}\code{p1p2} is an edge
#' in the underlying or reflexivity graph of CS-PCDs \eqn{)}
#' for points \code{p1} and \code{p2} in the standard basic triangle.
#'
#' More specifically, when the argument \code{ugraph="underlying"}, it returns
#' the edge indicator for the CS-PCD underlying graph,
#' that is, returns 1 if \code{p2} is
#' in \eqn{N_{CS}(p1,t)} or \code{p1} is in \eqn{N_{CS}(p2,t)},
#' returns 0 otherwise.
#' On the other hand,
#' when \code{ugraph="reflexivity"}, it returns
#' the edge indicator for the CS-PCD reflexivity graph,
#' that is, returns 1 if \code{p2} is
#' in \eqn{N_{CS}(p1,t)} and \code{p1} is in \eqn{N_{CS}(p2,t)},
#' returns 0 otherwise.
#'
#' In both cases \eqn{N_{CS}(x,t)} is the CS proximity region for point \eqn{x}
#' with expansion parameter \eqn{t > 0}.
#' CS proximity region is defined
#' with respect to the standard basic triangle \eqn{T_b=T((0,0),(1,0),(c_1,c_2))}
#' where \eqn{c_1} is
#' in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#'
#' Edge regions are based on the center, \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)} in
#' barycentric coordinates
#' in the interior of the standard basic triangle \eqn{T_b};
#' default is \eqn{M=(1,1,1)},
#' i.e., the center of mass of \eqn{T_b}.
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
#' @param p1 A 2D point whose CS proximity region is constructed.
#' @param p2 A 2D point. The function determines
#' whether there is an edge from \code{p1} to \code{p2} or not
#' in the underlying or reflexivity graphs of CS-PCDs.
#' @param t A positive real number
#' which serves as the expansion parameter
#' in CS proximity region; must be \eqn{> 0}
#' @param c1,c2 Positive real numbers
#' which constitute the vertex of the standard basic triangle
#' adjacent to the shorter edges;
#' \eqn{c_1} must be in \eqn{[0,1/2]}, \eqn{c_2>0} and \eqn{(1-c_1)^2+c_2^2 \le 1}.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the standard basic triangle;
#' default is \eqn{M=(1,1,1)} i.e., the center of mass of \eqn{T_b}.
#' @param ugraph The type of the graph based on CS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return Returns 1 if there is an edge between points \code{p1} and \code{p2}
#' in the underlying or reflexivity graph of CS-PCDs
#' in the standard basic triangle, and 0 otherwise.
#'
#' @seealso \code{\link{IedgeCStri}}, \code{\link{IedgeASbasic.tri}},
#' \code{\link{IedgePEbasic.tri}} and \code{\link[pcds]{IarcCSbasic.tri}}
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
#' t<-1.5
#'
#' P1<-as.numeric(pcds::runif.basic.tri(1,c1,c2)$g)
#' P2<-as.numeric(pcds::runif.basic.tri(1,c1,c2)$g)
#' IedgeCSbasic.tri(P1,P2,t,c1,c2,M)
#' IedgeCSbasic.tri(P1,P2,t,c1,c2,M,ugraph = "reflexivity")
#'
#' P1<-c(.4,.2)
#' P2<-c(.5,.26)
#' IedgeCSbasic.tri(P1,P2,t=2,c1,c2,M)
#' IedgeCSbasic.tri(P1,P2,t=2,c1,c2,M,ugraph="ref")
#' #}
#'
#' @export IedgeCSbasic.tri
IedgeCSbasic.tri <- function(p1,p2,t,c1,c2,M=c(1,1,1),
                             ugraph=c("underlying", "reflexivity"))
{
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  arc12 = pcds::IarcCSbasic.tri(p1,p2,t,c1,c2,M)
  arc21 = pcds::IarcCSbasic.tri(p2,p1,t,c1,c2,M)

  edge = ifelse(ugraph == "underlying",
                max(arc12,arc21),
                arc12*arc21)
  edge
} #end of the function
#'

#################################################################

#' @title The indicator for the presence of an edge from a point to another
#' for the underlying or reflexivity graphs of
#' Central Similarity Proximity Catch Digraphs (CS-PCDs) -
#' standard equilateral triangle case
#'
#' @description Returns \eqn{I(}\code{p1p2} is an edge
#' in the underlying or reflexivity graph of CS-PCDs \eqn{)}
#' for points \code{p1} and \code{p2} in the standard equilateral triangle.
#'
#' More specifically, when the argument \code{ugraph="underlying"}, it returns
#' the edge indicator for points \code{p1} and \code{p2}
#' in the standard equilateral triangle,
#' for the CS-PCD underlying graph,
#' that is, returns 1 if \code{p2} is
#' in \eqn{N_{CS}(p1,t)} or \code{p1} is in \eqn{N_{CS}(p2,t)},
#' returns 0 otherwise.
#' On the other hand,
#' when \code{ugraph="reflexivity"}, it returns
#' the edge indicator for points \code{p1} and \code{p2}
#' in the standard equilateral triangle,
#' for the CS-PCD reflexivity graph,
#' that is, returns 1 if \code{p2} is
#' in \eqn{N_{CS}(p1,t)} and \code{p1} is in \eqn{N_{CS}(p2,t)},
#' returns 0 otherwise.
#'
#' In both cases \eqn{N_{CS}(x,t)} is the CS proximity region
#' for point \eqn{x} with expansion parameter \eqn{t > 0}.
#' CS proximity region is defined
#' with respect to the standard equilateral triangle
#' \eqn{T_e=T(v=1,v=2,v=3)=T((0,0),(1,0),(1/2,\sqrt{3}/2))}
#' and edge regions are based on the center \eqn{M=(m_1,m_2)}
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
#' @param p1 A 2D point whose CS proximity region is constructed.
#' @param p2 A 2D point. The function determines
#' whether there is an edge from \code{p1} to \code{p2} or not
#' in the underlying or reflexivity graphs of CS-PCDs.
#' @param t A positive real number
#' which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center
#' in the interior of the standard equilateral triangle \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \eqn{T_e}.
#' @param ugraph The type of the graph based on CS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return Returns 1 if there is an edge between points \code{p1} and \code{p2}
#' in the underlying or reflexivity graph of CS-PCDs
#' in the standard equilateral triangle, and 0 otherwise.
#'
#' @seealso \code{\link{IedgeCSbasic.tri}}, \code{\link{IedgeCStri}},
#' and \code{\link[pcds]{IarcCSstd.tri}}
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
#' IedgeCSstd.tri(Xp[1,],Xp[3,],t=1.5,M)
#' IedgeCSstd.tri(Xp[1,],Xp[3,],t=1.5,M,ugraph="reflexivity")
#'
#' P1<-c(.4,.2)
#' P2<-c(.5,.26)
#' t<-2
#' IedgeCSstd.tri(P1,P2,t,M)
#' IedgeCSstd.tri(P1,P2,t,M,ugraph = "reflexivity")
#' #}
#'
#' @export IedgeCSstd.tri
IedgeCSstd.tri <- function(p1,p2,t,M=c(1,1,1),
                           ugraph=c("underlying", "reflexivity"))
{
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  arc12 = pcds::IarcCSstd.tri(p1,p2,t,M)
  arc21 = pcds::IarcCSstd.tri(p2,p1,t,M)

  edge = ifelse(ugraph == "underlying",max(arc12,arc21),arc12*arc21)
  edge
} #end of the function
#'

#################################################################

#' @title Number of edges in the underlying or reflexivity graphs of
#' Central Similarity Proximity Catch Digraphs (CS-PCDs) -
#' standard equilateral triangle case
#'
#' @description
#' An object of class \code{"NumEdges"}.
#' Returns the number of edges of
#' the underlying or reflexivity graphs of
#' Central Similarity Proximity Catch Digraphs (CS-PCDs)
#' whose vertices are the
#' given 2D numerical data set, \code{Xp}
#' in the standard equilateral triangle.
#' It also provides number of vertices
#' (i.e., number of data points inside the triangle)
#' and indices of the data points that reside in the triangle.
#'
#' CS proximity region \eqn{N_{CS}(x,t)} is defined
#' with respect to the standard equilateral triangle
#' \eqn{T_e=T(v=1,v=2,v=3)=T((0,0),(1,0),(1/2,\sqrt{3}/2))}
#' with expansion parameter \eqn{t > 0}
#' and edge regions are based on the center \eqn{M=(m_1,m_2)}
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
#' or reflexivity graphs based on the CS-PCD.
#' @param t A positive real number
#' which serves as the expansion parameter for CS proximity region.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center
#' in the interior of the standard equilateral triangle \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \eqn{T_e}.
#' @param ugraph The type of the graph based on CS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of edges
#' and quantities related to the standard equilateral triangle}
#' \item{und.graph}{Type of the graph as "Underlying" or "Reflexivity" for the CS-PCD}
#' \item{num.edges}{Number of edges of the underlying
#' or reflexivity graphs based on the CS-PCD
#' with vertices in the standard equilateral triangle \eqn{T_e}}
#' \item{num.in.tri}{Number of \code{Xp} points
#' in the standard equilateral triangle, \eqn{T_e}}
#' \item{ind.in.tri}{The vector of indices of the \code{Xp} points
#' that reside in \eqn{T_e}}
#' \item{tess.points}{Tessellation points,
#' i.e., points on which the tessellation of the study region is performed,
#' here, tessellation is the support triangle \eqn{T_e}.}
#' \item{vertices}{Vertices of the underlying or reflexivity graph, \code{Xp}.}
#'
#' @seealso \code{\link{num.edgesCStri}}, \code{\link{num.edgesCS}},
#' and \code{\link[pcds]{num.arcsCSstd.tri}}
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
#' Nedges = num.edgesCSstd.tri(Xp,t=1.5,M)
#' Nedges
#' summary(Nedges)
#' plot(Nedges)
#' #}
#'
#' @export num.edgesCSstd.tri
num.edgesCSstd.tri <- function(Xp,t,M=c(1,1,1),
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

  if (!pcds::is.point(t,1) || t<=0)
  {stop('t must be a scalar > 0')}

  if (!pcds::is.point(M) && !pcds::is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
        3D point for barycentric coordinates')}

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
        edges<-edges+IedgeCSstd.tri(p1,p2,t,M,ugraph)
      }
      }
    }
    }
  }

  NinTri = length(ind.in.tri)

  und.graph = ifelse(ugraph=="underlying","Underlying", "Reflexivity")
  desc<-paste("Number of Edges of the ",und.graph,
              " Graphs of the CS-PCD and the Related Quantities with vertices Xp in the Standard Equilateral Triangle",sep="")

  res<-list(desc=desc, #description of the output
            und.graph = und.graph, #"Underlying" or "Reflexivity"
            num.edges=edges-n, # -n is to avoid loops
            #number of edges for the underlying or reflexivity graph of the CS-PCD
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

#' @title Incidence matrix for the underlying or reflexivity graphs of
#' Central Similarity Proximity Catch Digraphs (CS-PCDs) -
#' standard equilateral triangle case
#'
#' @description Returns the incidence matrix
#' for the underlying or reflexivity graphs of the CS-PCD
#' whose vertices are the given 2D numerical data set, \code{Xp},
#' in the standard equilateral triangle
#' \eqn{T_e=T(v=1,v=2,v=3)=T((0,0),(1,0),(1/2,\sqrt{3}/2))}.
#'
#' CS proximity region is constructed
#' with respect to the standard equilateral triangle \eqn{T_e} with
#' expansion parameter \eqn{t > 0} and edge regions are based on
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
#' which constitute the vertices of the underlying or reflexivity graphs of the CS-PCD.
#' @param t A positive real number
#' which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center
#' in the interior of the standard equilateral triangle \eqn{T_e};
#' default is \eqn{M=(1,1,1)} i.e.
#' the center of mass of \eqn{T_e}.
#' @param ugraph The type of the graph based on CS-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return Incidence matrix for the underlying or reflexivity graphs
#' of the CS-PCD with vertices
#' being 2D data set, \code{Xp}
#' in the standard equilateral triangle where CS proximity
#' regions are defined with \code{M}-edge regions.
#'
#' @seealso \code{\link{inci.mat.undCStri}}, \code{\link{inci.mat.undCS}},
#' and \code{\link[pcds]{inci.matCSstd.tri}}
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
#' inc.mat<-inci.mat.undCSstd.tri(Xp,t=1.5,M)
#' inc.mat
#' (sum(inc.mat)-n)/2
#' num.edgesCSstd.tri(Xp,t=1.5,M)$num.edges
#'
#' pcds::dom.num.greedy(inc.mat)
#' pcds::Idom.num.up.bnd(inc.mat,2)
#' #}
#'
#' @export inci.mat.undCSstd.tri
inci.mat.undCSstd.tri <- function(Xp,t,M=c(1,1,1),
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

  if (!pcds::is.point(t,1) || t<=0)
  {stop('t must be a scalar > 0')}

  if (!pcds::is.point(M) && !pcds::is.point(M,3))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
        3D point for barycentric coordinates')}

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
  for (j in i:n)
  {p2<-Xp[j,]
  inc.mat[i,j]<-inc.mat[j,i]<-IedgeCSstd.tri(p1,p2,t,M,ugraph)
  }
  }
  inc.mat
} #end of the function
#'

#################################################################

# funsMuVarUndCS2D
#'
#' @title Returns the mean and (asymptotic) variance of edge density of
#' underlying or reflexivity graphs of
#' Central Similarity Proximity Catch Digraph (CS-PCD)
#' for 2D uniform data in one triangle
#'
#' @description
#' The mean and (asymptotic) variance functions
#' for the underlying or reflexivity graphs of
#' Central Similarity Proximity Catch Digraphs (CS-PCDs):
#' \code{muOrCS2D} and \code{asyvarOrCS2D} for the underlying graph
#' and
#' \code{muAndCS2D} and \code{asyvarAndCS2D} for the reflexivity graph.
#'
#' \code{muOrCS2D} and \code{muAndCS2D} return the mean of the (edge) density of
#' the underlying or reflexivity graphs of CS-PCDs, respectively,
#' for 2D uniform data in a triangle.
#' Similarly, \code{asyvarOrCS2D} and \code{asyvarAndCS2D} return the asymptotic variance
#' of the edge density of the underlying or reflexivity graphs of CS-PCDs, respectively,
#' for 2D uniform data in a triangle.
#'
#' CS proximity regions are defined with expansion parameter \eqn{t > 0}
#' with respect to the triangle in which the points reside and
#' edge regions are based on center of mass, \eqn{CM} of the triangle.
#'
#' See also (\insertCite{ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param t A positive real number which serves
#' as the expansion parameter in CS proximity region.
#' @param ugraph The type of the graph based on CS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return \code{mu.undCS2D} returns the mean
#' and \code{asyvarUndOrCS2D} returns the (asymptotic) variance of the
#' edge density of the underlying graph of the CS-PCD for uniform data in any triangle
#' if \code{ugraph="underlying"},
#' and those of the reflexivity graph if \code{ugraph="reflexivity"}.
#' The functions \code{muOrCS2D}, \code{muAndCS2D}, \code{asyvarOrCS2D},
#' and \code{asyvarAndCS2D} are the corresponding mean and asymptotic variance functions
#' for the edge density of the reflexivity graph of the CS-PCD, respectively,
#' for uniform data in any triangle.
#'
#' @name funsMuVarUndCS2D
NULL
#'
#' @seealso \code{\link{mu.undCS2D}}, \code{\link{asy.var.undCS2D}}
#' \code{\link[pcds]{muCS2D}}, and \code{\link[pcds]{asyvarCS2D}},
#'
#'
#' @rdname funsMuVarUndCS2D
#'
#' @author Elvan Ceyhan
#'
#' @examples
#' #\donttest{
#' mu.undCS2D(1.2)
#' mu.undCS2D(1.2,ugraph="r")
#'
#' tseq<-seq(0.01,10,by=.05)
#' ltseq<-length(tseq)
#'
#' muOR = muAND <- vector()
#' for (i in 1:ltseq)
#' {
#'   muOR<-c(muOR,mu.undCS2D(tseq[i]))
#'   muAND<-c(muAND,mu.undCS2D(tseq[i],ugraph="r"))
#' }
#'
#' plot(tseq, muOR,type="l",xlab="t",ylab=expression(mu(t)),lty=1,
#'      xlim=range(tseq),ylim=c(0,1))
#' lines(tseq,muAND,type="l",lty=2,col=2)
#' legend("bottomright", inset=.02,
#'        legend=c(expression(mu[or](t)),expression(mu[and](t))),
#'        lty=1:2,col=1:2)
#' #}
#'
#' @export muOrCS2D
muOrCS2D<-function(t)
{
  if (!pcds::is.point(t,1) || t<=0)
  {stop('The argument must be a scalar greater than 0')}

  mn<-0;
  if (t < 1/2)
  {
    mn<- 2/9*t^2*(3+t^2+t)/(t+3);
  } else {
    if (t < 1)
    {
      mn<- 1/3*t^2*(2*t^2+11*t+9)/(t+3)/(2*t+5);
    } else {
      mn<- (4*t^3+32*t^2+45*t-15)*t/(t+3)/(2*t+5)/(t+2)/(2*t+1);
    }}
  mn
} #end of the function
#'
#' @rdname funsMuVarUndCS2D
#'
#' @export muAndCS2D
muAndCS2D <- function(t)
{
  if (!pcds::is.point(t,1) || t<=0)
  {stop('The argument must be a scalar greater than 0')}

  mn<-0;
  if (t < 1/2)
  {
    mn<- -(2*t^2 - t - 3)*t^2/(9*(3 + t));
  } else {
    if (t < 1)
    {
      mn<- 2*t^2/((t + 3)*(2*t + 5));
    } else {
      mn<- 2*t^2/((t + 3)*(2*t + 5));
    }}
  mn
} #end of the function
#'
#' @rdname funsMuVarUndCS2D
#'
#' @export mu.undCS2D
mu.undCS2D <-function(t,ugraph=c("underlying", "reflexivity"))
{
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  mn = ifelse(ugraph == "underlying",muOrCS2D(t),muAndCS2D(t))
  mn
}
#' @references
#' \insertAllCited{}
#'
#' @rdname funsMuVarUndCS2D
#'
#' @examples
#' #\donttest{
#' asy.var.undCS2D(1.2)
#' asy.var.undCS2D(1.2,ugraph="r")
#'
#' asyvarOrCS2D(.2)
#'
#' tseq<-seq(.05,25,by=.05)
#' ltseq<-length(tseq)
#'
#' avarOR<-avarAND<-vector()
#' for (i in 1:ltseq)
#' {
#'   avarOR<-c(avarOR,asy.var.undCS2D(tseq[i]))
#'   avarAND<-c(avarAND,asy.var.undCS2D(tseq[i],ugraph="r"))
#' }
#'
#' oldpar <- par(mar=c(5,5,4,2))
#' plot(tseq, 4*avarAND,type="l",lty=2,col=2,xlab="t",
#'      ylab=expression(paste(sigma^2,"(t)")),xlim=range(tseq))
#' lines(tseq,4*avarOR,type="l")
#' legend(18,.1,
#'        legend=c(expression(paste(sigma["underlying"]^"2","(t)")),
#'                  expression(paste(sigma["reflexivity"]^"2","(t)")) ),
#'        lty=1:2,col=1:2)
#'
#' par(oldpar)
#' #}
#'
#' @export asyvarOrCS2D
asyvarOrCS2D<-function(t)
{
  if (!pcds::is.point(t,1) || t<=0)
  {stop('The argument must be a scalar greater than 0')}

  #internal OR underlying functions for covariance (i.e. asymptotic variance)
  CScov.or0 <-function(t)
  {2/1215*(192*t^19+1968*t^18+5048*t^17-13524*t^16-96750*t^15-153498*t^14+204921*t^13+1191711*t^12+1644184*t^11-
             782739*t^10-5793519*t^9-6962886*t^8+2566647*t^7+16522461*t^6+17799264*t^5+449793*t^4-21194217*t^3-24944922*t^2-
             11298042*t-1653372)/(t+3)^3/(2*t+3)/(2*t^2+6*t+3)/(2*t^2-t-3)/(t^2-t-3)/(t^2+t+3)/(t-3)/(t+2)/(2*t+1)*t^4}

  CScov.or1 <-function(t)
  {-1/1215*(20480*t^24+294912*t^23+853760*t^22-8420352*t^21-67953472*t^20-100448064*t^19+697856240*t^18+
              3143699712*t^17+1244152292*t^16-20840169012*t^15-48748490950*t^14+16404525522*t^13+228604700816*t^12+
              301025320353*t^11-197975773249*t^10-972536561376*t^9-940947719281*t^8+253203630786*t^7+1480868308641*t^6+
              1640646436698*t^5+980894429019*t^4+345696394563*t^3+67789274544*t^2+5712392970*t-12247200)/(4*t+5)/(2*t+3)/
      (t+3)^3/(4*t^2+11*t+4)/(2*t^2+6*t+3)/(2*t+5)^3/(2*t-5)/(t^2-t-3)/(t+2)/(t^2-4)/(2*t+1)/(t+1)*t^3}

  CScov.or2 <-function(t)
  {1/30*(111431680*t^20+2820980736*t^19+32710127872*t^18+229961057152*t^17+1093143530432*t^16+3707285176320*t^15+
           9212280477424*t^14+16915566759464*t^13+22729421834796*t^12+21484172948352*t^11+12608245109798*t^10+
           2015653755166*t^9-3733800929869*t^8-3725129049589*t^7-1551133260126*t^6-141206916516*t^5+170413896339*t^4+
           91611623187*t^3+21859209990*t^2+2646270000*t+131220000)/(t+3)^3/(4*t+3)^2/(2*t^2+5*t+1)^2/(t+1)/(4*t+5)/
      (4*t^2+11*t+4)/(2*t+3)/(2*t^2+6*t+3)/(2*t+5)^3/(t+2)^2/(2*t+1)^2}
####
  asyvar<-0;
  if (t < 1/2)
  {
    asyvar<-CScov.or0(t);
  } else {
    if (t < 1)
    {
      asyvar<- CScov.or1(t);
    } else {
      asyvar<- CScov.or2(t);
    }}
  asyvar #need to multiply this by 4 in the asymptotic approximation
} #end of the function
#'
#' @rdname funsMuVarUndCS2D
#'
#' @export asyvarAndCS2D
asyvarAndCS2D<-function(t)
{
  if (!pcds::is.point(t,1) || t<=0)
  {stop('The argument must be a scalar greater than 0')}

  #internal AND underlying functions for covariance (i.e. asymptotic variance)
  CScov.and0<-function(t)
  {1/1215*(224*t^13+2112*t^12+5872*t^11+1256*t^10-7854*t^9+22378*t^8+23973*t^7-109758*t^6-82737*t^5+
             186138*t^4+136323*t^3-132192*t^2-137781*t-30618)/(t+3)^3/(2*t^2+6*t+3)/(2*t+3)/(2*t^2-t-3)/(t+2)*t^4}

  CScov.and1<-function(t)
  {2/15*(16128*t^11+268864*t^10+1943688*t^9+7979884*t^8+20460770*t^7+33935927*t^6+36470074*t^5+
           24743401*t^4+9992165*t^3+2141319*t^2+184113*t+135)/(4*t+5)/(2*t+3)/(t+3)^3/(4*t^2+11*t+4)/
      (2*t^2+6*t+3)/(2*t+5)^3/(t+2)/(2*t+1)/(t+1)*t^3}

  CScov.and2<-function(t)
  {2/15*(16128*t^11+268864*t^10+1943688*t^9+7979884*t^8+20460770*t^7+33935927*t^6+36470074*t^5+
           24743401*t^4+9992165*t^3+2141319*t^2+184113*t+135)/(t+3)^3/(t+1)/(4*t+5)/(4*t^2+11*t+4)/(2*t+3)/
      (t+2)/(2*t^2+6*t+3)/(2*t+5)^3/(2*t+1)*t^3}
###
  asyvar<-0;
  if (t < 1/2)
  {
    asyvar<-CScov.and0(t);
  } else {
    if (t < 1)
    {
      asyvar<- CScov.and1(t);
    } else {
      asyvar<- CScov.and2(t);
    }}
  asyvar #need to multiply this by 4 in the asymptotic approximation
} #end of the function
#'
#' @rdname funsMuVarUndCS2D
#'
#' @export asy.var.undCS2D
asy.var.undCS2D <-function(t,ugraph=c("underlying", "reflexivity"))
{
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  asyvar = ifelse(ugraph == "underlying",asyvarOrCS2D(t),asyvarAndCS2D(t))
  asyvar #need to multiply this by 4 in the asymptotic approximation
}

#################################################################

#' @title The indicator for the presence of an edge from a point to another
#' for the underlying or reflexivity graphs of
#' Central Similarity Proximity Catch Digraphs (CS-PCDs) -
#' one triangle case
#'
#' @description Returns \eqn{I(}\code{p1p2} is an edge
#' in the underlying or reflexivity graph of CS-PCDs \eqn{)}
#' for points \code{p1} and \code{p2} in a given triangle.
#'
#' More specifically, when the argument \code{ugraph="underlying"}, it returns
#' the edge indicator for the CS-PCD underlying graph,
#' that is, returns 1 if \code{p2} is
#' in \eqn{N_{CS}(p1,t)} or \code{p1} is in \eqn{N_{CS}(p2,t)},
#' returns 0 otherwise.
#' On the other hand,
#' when \code{ugraph="reflexivity"}, it returns
#' the edge indicator for the CS-PCD reflexivity graph,
#' that is, returns 1 if \code{p2} is
#' in \eqn{N_{CS}(p1,t)} and \code{p1} is in \eqn{N_{CS}(p2,t)},
#' returns 0 otherwise.
#'
#' In both cases CS proximity region is constructed
#' with respect to the triangle \code{tri} and
#' edge regions are based on the center, \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of \code{tri};
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
#' @param p1 A 2D point whose CS proximity region is constructed.
#' @param p2 A 2D point. The function determines
#' whether there is an edge from \code{p1} to \code{p2} or not
#' in the underlying or reflexivity graphs of CS-PCDs.
#' @param t A positive real number
#' which serves as the expansion parameter in CS proximity region.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' @param ugraph The type of the graph based on CS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return Returns 1 if there is an edge between points \code{p1} and \code{p2}
#' in the underlying or reflexivity graph of CS-PCDs
#' in a given triangle \code{tri}, and 0 otherwise.
#'
#' @seealso \code{\link{IedgeCSbasic.tri}}, \code{\link{IedgeAStri}},
#' \code{\link{IedgePEtri}} and \code{\link[pcds]{IarcCStri}}
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
#' t<-1.5
#' n<-3
#' set.seed(1)
#' Xp<-pcds::runif.tri(n,Tr)$g
#'
#' IedgeCStri(Xp[1,],Xp[2,],Tr,t,M)
#' IedgeCStri(Xp[1,],Xp[2,],Tr,t,M,ugraph = "reflexivity")
#'
#' P1<-as.numeric(pcds::runif.tri(1,Tr)$g)
#' P2<-as.numeric(pcds::runif.tri(1,Tr)$g)
#' IedgeCStri(P1,P2,Tr,t,M)
#' IedgeCStri(P1,P2,Tr,t,M,ugraph="r")
#' #}
#'
#' @export IedgeCStri
IedgeCStri <- function(p1,p2,tri,t,M=c(1,1,1),ugraph=c("underlying", "reflexivity"))
{
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  arc12 = pcds::IarcCStri(p1,p2,tri,t,M)
  arc21 = pcds::IarcCStri(p2,p1,tri,t,M)

  edge = ifelse(ugraph == "underlying",
                max(arc12,arc21),
                arc12*arc21)
  edge
} #end of the function
#'

#################################################################

#' @title Number of edges in the underlying or reflexivity graphs of
#' Central Similarity Proximity Catch Digraphs (CS-PCDs) -
#' one triangle case
#'
#' @description
#' An object of class \code{"NumEdges"}.
#' Returns the number of edges of
#' the underlying or reflexivity graphs of
#' Central Similarity Proximity Catch Digraphs (CS-PCDs)
#' whose vertices are the
#' given 2D numerical data set, \code{Xp}
#' in a given triangle.
#' It also provides number of vertices
#' (i.e., number of data points inside the triangle)
#' and indices of the data points that reside in the triangle.
#'
#' CS proximity region \eqn{N_{CS}(x,t)} is defined
#' with respect to the triangle, \code{tri}
#' with expansion parameter \eqn{t > 0} and edge regions are
#' based on the center \eqn{M=(m_1,m_2)} in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' For the number of edges, loops are not allowed,
#' so edges are only possible for points
#' inside the triangle \code{tri} for this function.
#'
#' See also
#' (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-CS;textual}{pcds.ugraph}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of CS-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param t A positive real number
#' which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' @param ugraph The type of the graph based on CS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of edges
#' and quantities related to the triangle}
#' \item{und.graph}{Type of the graph as "Underlying" or "Reflexivity" for the CS-PCD}
#' \item{num.edges}{Number of edges of the underlying
#' or reflexivity graphs based on the CS-PCD
#' with vertices in the given triangle \code{tri}}
#' \item{num.in.tri}{Number of \code{Xp} points in the triangle, \code{tri}}
#' \item{ind.in.tri}{The vector of indices of the \code{Xp} points
#' that reside in the triangle}
#' \item{tess.points}{Tessellation points,
#' i.e., points on which the tessellation of the study region is performed,
#' here, tessellation is the support triangle.}
#' \item{vertices}{Vertices of the underlying or reflexivity graph, \code{Xp}.}
#'
#' @seealso \code{\link{num.edgesCS}}, \code{\link{num.edgesAStri}},
#' \code{\link{num.edgesPEtri}}, and \code{\link[pcds]{num.arcsCStri}}
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
#' Nedges = num.edgesCStri(Xp,Tr,t=1.5,M)
#' Nedges
#' summary(Nedges)
#' plot(Nedges)
#' #}
#'
#' @export num.edgesCStri
num.edgesCStri <- function(Xp,tri,t,M=c(1,1,1),
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

  if (!pcds::is.point(t,1) || t<=0)
  {stop('t must be a scalar > 0')}

  if (!(pcds::is.point(M) || pcds::is.point(M,3) ))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates ')}

  if (pcds::dimension(M)==3)
  {M<-pcds::bary2cart(M,tri)}

  if (!(pcds::in.triangle(M,tri,boundary=FALSE)$in.tri))
   {stop('M is not a center in the interior of the triangle')}

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
        edges.in.tri<-edges.in.tri+IedgeCStri(Xp[i,],Xp[k,],tri,t,M,ugraph=ugraph)
      }
      }

      for (j in (i:n)[-1]) #to avoid loops
      {
        tot.edges<-tot.edges+IedgeCStri(Xp[i,],Xp[j,],tri,t,M,ugraph=ugraph)
      }
    }
  }

  NinTri = length(ind.in.tri)

  und.graph = ifelse(ugraph=="underlying",
                     "Underlying",
                     "Reflexivity")
  desc<-paste("Number of Edges of the ",und.graph,
              " Graphs of the CS-PCD and the Related Quantities with vertices Xp in One Triangle",sep="")

  res<-list(desc=desc, #description of the output
            und.graph = und.graph, #"Underlying" or "Reflexivity"
            num.edges=tot.edges, #total number of edges in the entire underlying or reflexivity graph
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

#' @title Edge density of the underlying or reflexivity graphs of
#' Central Similarity Proximity Catch Digraphs (CS-PCDs) -
#' one triangle case
#'
#' @description Returns the edge density
#' of the underlying or reflexivity graphs of
#' Central Similarity Proximity Catch Digraphs (CS-PCDs)
#' whose vertex set is the given 2D numerical data set, \code{Xp},
#' (some of its members are) in the triangle \code{tri}.
#'
#' CS proximity regions is defined with respect to \code{tri} with
#' expansion parameter \eqn{t > 0} and edge regions are
#' based on center \eqn{M=(m_1,m_2)} in Cartesian coordinates or
#' \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri}; default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' The function also provides edge density standardized
#' by the mean and asymptotic variance of the edge density
#' of the underlying or reflexivity graphs of CS-PCD
#' for uniform data in the triangle \code{tri}
#' only when \code{M} is the center of mass.
#' For the number of edges, loops are not allowed.
#'
#' \code{in.tri.only} is a logical argument (default is \code{FALSE}) for considering only the points
#' inside the triangle or all the points as the vertices of the digraph.
#' if \code{in.tri.only=TRUE}, edge density is computed only for
#' the points inside the triangle (i.e., edge density of the subgraph of the underlying or reflexivity graph
#' induced by the vertices in the triangle is computed),
#' otherwise edge density of the entire graph
#' (i.e., graph with all the vertices) is computed.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying or reflexivity graphs of the CS-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param t A positive real number
#' which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' @param ugraph The type of the graph based on CS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
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
#' or reflexivity graphs based on the CS-PCD
#' whose vertices are the 2D numerical data set, \code{Xp};
#' CS proximity regions are defined
#' with respect to the triangle \code{tri} and \code{M}-edge regions}
#' \item{std.edge.dens}{Edge density standardized
#' by the mean and asymptotic variance of the edge
#' density of the underlying or reflexivity graphs
#' based on the CS-PCD for uniform data in the triangle \code{tri}.
#' This will only be returned, if \code{M} is the center of mass.}
#'
#' @seealso \code{\link{ASedge.dens.tri}}, \code{\link{PEedge.dens.tri}},
#' and \code{\link[pcds]{CSarc.dens.tri}}
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
#' num.edgesCStri(Xp,Tr,t=1.5,M)$num.edges
#' CSedge.dens.tri(Xp,Tr,t=1.5,M)
#' CSedge.dens.tri(Xp,Tr,t=1.5,M,in.tri.only = TRUE)
#'
#' #For the reflexivity graph
#' num.edgesCStri(Xp,Tr,t=1.5,M,ugraph="r")$num.edges
#' CSedge.dens.tri(Xp,Tr,t=1.5,M,ugraph="r")
#' CSedge.dens.tri(Xp,Tr,t=1.5,M,in.tri.only = TRUE,ugraph="r")
#' #}
#'
#' @export CSedge.dens.tri
CSedge.dens.tri <- function(Xp,tri,t,M=c(1,1,1),
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

  nedges<-num.edgesCStri(Xp,tri,t,M,ugraph)$num.edges

  mean.rho<-mu.undCS2D(t,ugraph)
  var.rho<-4*asy.var.undCS2D(t,ugraph)

  if (in.tri.only==TRUE)
  {
    ind.it<-c()
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

#' @title Number of edges of the underlying or reflexivity graphs of
#' Central Similarity Proximity Catch Digraphs (CS-PCDs) -
#' multiple triangle case
#'
#' @description
#' An object of class \code{"NumEdges"}.
#' Returns the number of edges of
#' the underlying or reflexivity graph of
#' Central Similarity Proximity Catch Digraph (CS-PCD)
#' and various other quantities and vectors such as
#' the vector of number of vertices (i.e., number of data points)
#' in the Delaunay triangles,
#' number of data points in the convex hull of \code{Yp} points,
#' indices of the Delaunay triangles for the data points, etc.
#'
#' CS proximity regions are defined with respect to the
#' Delaunay triangles based on \code{Yp} points
#' with expansion parameter \eqn{t > 0}
#' and edge regions in each triangle
#' is based on the center \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of each
#' Delaunay triangle (default for \eqn{M=(1,1,1)}
#' which is the center of mass of the triangle).
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
#' For the number of edges,
#' loops are not allowed so edges are only possible
#' for points inside the convex hull of \code{Yp} points.
#'
#' See (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph})
#' for more on CS-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds.ugraph})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying or reflexivity graphs of the CS-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param t A positive real number
#' which serves as the expansion parameter in CS proximity region.
#' @param M A 3D point in barycentric coordinates
#' which serves as a center in the interior of each Delaunay triangle,
#' default for \eqn{M=(1,1,1)}
#' which is the center of mass of each triangle.
#' @param ugraph The type of the graph based on CS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return A \code{list} with the elements
#' \item{desc}{A short description of the output: number of edges
#' and related quantities for the induced subgraphs of the underlying
#' or reflexivity graphs (of CS-PCD) in the Delaunay triangles}
#' \item{und.graph}{Type of the graph as "Underlying" or "Reflexivity" for the CS-PCD}
#' \item{num.edges}{Total number of edges in all triangles,
#' i.e., the number of edges for the entire underlying
#' or reflexivity graphs of the CS-PCD}
#' \item{num.in.conv.hull}{Number of \code{Xp} points
#' in the convex hull of \code{Yp} points}
#' \item{num.in.tris}{The vector of number of \code{Xp} points
#' in the Delaunay triangles based on \code{Yp} points}
#' \item{weight.vec}{The \code{vector} of the areas of
#' Delaunay triangles based on \code{Yp} points}
#' \item{tri.num.edges}{The \code{vector} of the number of edges
#' of the components of the CS-PCD in the
#' Delaunay triangles based on \code{Yp} points}
#' \item{del.tri.ind}{A matrix of indices of vertices of
#' the Delaunay triangles based on \code{Yp} points,
#' each column corresponds to the vector of
#' indices of the vertices of one triangle.}
#' \item{data.tri.ind}{A \code{vector} of indices of vertices of
#' the Delaunay triangles in which data points reside,
#' i.e., column number of \code{del.tri.ind} for each \code{Xp} point.}
#' \item{tess.points}{Tessellation points,
#' i.e., points on which the tessellation of the study region is performed,
#' here, tessellation is the Delaunay triangulation based on \code{Yp} points.}
#' \item{vertices}{Vertices of the underlying or reflexivity graph, \code{Xp}.}
#'
#' @seealso \code{\link{num.edgesCStri}}, \code{\link{num.edgesAS}},
#' \code{\link{num.edgesPE}}, and \code{\link[pcds]{num.arcsCS}}
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
#' Nedges = num.edgesCS(Xp,Yp,t=1.5,M)
#' Nedges
#' summary(Nedges)
#' plot(Nedges)
#' #}
#'
#' @export num.edgesCS
num.edgesCS <- function(Xp,Yp,t,M=c(1,1,1),ugraph=c("underlying", "reflexivity"))
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

  if (!pcds::is.point(t,1) || t<=0)
  {stop('t must be a scalar > 0')}

  if (!pcds::is.point(M,3) || !all(M>0))
  {stop('M must be a numeric 3D point with positive barycentric coordinates')}

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
  #NumEdges = num.edgesCStri(Xp,tri,t,M,ugraph)
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
    NumEdges = num.edgesCStri(Xdt,tri,t,M,ugraph) #for the vertices inside the triangle
  #  Wvec<-pcds::area.polygon(tri)
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
              " Graphs of the CS-PCD with vertices Xp and the Related Quantities for the Induced Subdigraph for the Points in the Delaunay Triangle",sep="")

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
      # data.del.tris = del.tris=list()
      data.tri.ind = rep(NA,nx)

      for (i in 1:ndt)
      {
        dt.ind=which(Tri.Ind==i)
        #which indices of data points residing in ith Delaunay triangle
        Xpi<-Xp[dt.ind,] #points in ith Delaunay triangle
        data.tri.ind[dt.ind] = i
        #assigning the index of the Delaunay triangle that contains the data point
        ifelse(ndt==1,
               Tri<-Yp[Ytri,],
               Tri<-Yp[Ytri[i,],])
        #vertices of ith triangle
        tri<-pcds::as.basic.tri(Tri)$tri
        #convert the triangle Tri into an nonscaled basic triangle, see as.basic.tri help page
        ni.vec<-c(ni.vec,length(Xpi)/2)
        #number of points in ith Delaunay triangle

        num.edges<-num.edgesCStri(Xpi,tri,t,M,ugraph)$num.edges
        #number of edges in ith triangle
        edges<-c(edges,num.edges)
        #number of edges in all triangles as A \code{vector}

      }

      Tot.Edges<-sum(edges)  #the total number of edges in all triangles
    }

    Tot.Edges = Tot.Edges + sum(duplicated(Xp[!inCH,]))

    desc<-paste("Number of Edges of the ",und.graph,
                " Graphs of the CS-PCD with vertices Xp and the Related Quantities for the Induced Subdigraphs for the Points in the Delaunay Triangles",sep="")

    res<-list(desc=desc, #description of the output
              und.graph = und.graph, #"Underlying" or "Reflexivity"
              num.edges=Tot.Edges, #number of edges for the entire CS-PCD
              tri.num.edges=edges,
              #vector of number of edges for the Delaunay triangles
              num.in.conv.hull=NinCH,
              # number of Xp points in CH of Yp points
              ind.in.conv.hull= ind.in.CH, #indices of Xp points in CH of Yp points
              num.in.tris=ni.vec,
              # vector of number of Xp points in the Delaunay triangles
              weight.vec=Wvec, #areas of Delaunay triangles
              del.tri.ind=t(Ytri),
              # indices of the Delaunay triangles, each column corresponds to the vector of indices of the vertices of one triangle
              data.tri.ind=data.tri.ind, #indices of the Delaunay triangles in which data points reside, i.e., column number of del.tri.ind for each Xp point
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

#' @title A test of segregation/association based on edge density
#' of underlying or reflexivity graph of Central Similarity Proximity Catch Digraph
#' (CS-PCD) for 2D data
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
#' CS-PCD for uniform 2D data.
#'
#' The function yields the test statistic,
#' \eqn{p}-value for the corresponding \code{alternative},
#' the confidence interval, estimate
#' and null value for the parameter of interest
#' (which is the edge density),
#' and method and name of the data set used.
#'
#' Under the null hypothesis of uniformity of \code{Xp} points
#' in the convex hull of \code{Yp} points, edge density
#' of underlying or reflexivity graph of CS-PCD
#' whose vertices are \code{Xp} points equals
#' to its expected value under the uniform distribution and
#' \code{alternative} could be two-sided, or left-sided
#' (i.e., data is accumulated around the \code{Yp} points, or association)
#' or right-sided (i.e., data is accumulated
#' around the centers of the triangles,
#' or segregation).
#'
#' CS proximity region is constructed
#' with the expansion parameter \eqn{t > 0} and \eqn{CM}-edge regions
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
#' the case of non-concurring supports is
#' an ongoing topic of research of the author of the package.
#'
#' \code{ch.cor} is for convex hull correction
#' (default is \code{"no convex hull correction"}, i.e., \code{ch.cor=FALSE})
#' which is recommended
#' when both \code{Xp} and \code{Yp} have the same rectangular support.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph})
#' for more on the test based on the edge density of
#' underlying or reflexivity graphs of CS-PCDs.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graphs of the CS-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param t A positive real number
#' which serves as the expansion parameter in CS proximity region.
#' @param ugraph The type of the graph based on CS-PCDs,
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
#' underlying or reflexivity graphs of CS-PCD based on
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
#' @seealso \code{\link{PEedge.dens.test}} and \code{\link[pcds]{CSarc.dens.test}}
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
#' CSedge.dens.test(Xp,Yp,t=1.5)
#' CSedge.dens.test(Xp,Yp,t=1.5,ch=TRUE)
#'
#' CSedge.dens.test(Xp,Yp,t=1.5,ugraph="r")
#' CSedge.dens.test(Xp,Yp,t=1.5,ugraph="r",ch=TRUE)
#' #since Y points are not uniform, convex hull correction is invalid here
#'
#' }
#'
#' @export CSedge.dens.test
CSedge.dens.test <- function(Xp,Yp,t,ugraph=c("underlying", "reflexivity"),ch.cor=FALSE,
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

  if (!pcds::is.point(t,1) || t<=0)
  {stop('t must be a scalar > 0')}

  if (!missing(conf.level))
    if (length(conf.level) != 1 || is.na(conf.level) ||
        conf.level < 0 || conf.level > 1)
      stop("conf.level must be a number between 0 and 1")

  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  Edges<-num.edgesCS(Xp,Yp,t,M=c(1,1,1),ugraph)
  #use the default, i.e., CM for the center M
  NinCH<-Edges$num.in.con

  num.edges<-Edges$num.edges #total number of edges in the CS-PCD
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
      npe = pcds::NCStri(Xpi,tri,t)
      num.edges = num.edges+pcds::area.polygon(npe)/Wvec[i]
    }
  }
  asy.mean0<-mu.undCS2D(t,ugraph)  #asy mean value for the r value
  asy.mean<-asy.mean0*sum(LW^2)

  asy.var0<-4*asy.var.undCS2D(t,ugraph)  #asy variance value for the r value
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

  method <- c("Large Sample z-Test Based on Edge Density of", ifelse(ugraph ==  "underlying", "underlying","reflexivity"), "graph of CS-PCD for Testing Uniformity of 2D Data ---")

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

#' @title Incidence matrix for the underlying or reflexivity graphs of
#' Central Similarity Proximity Catch Digraphs (CS-PCDs) -
#' one triangle case
#'
#' @description Returns the incidence matrix
#' for the underlying or reflexivity graphs of the CS-PCD
#' whose vertices are the given 2D numerical data set, \code{Xp},
#' in the triangle \code{tri}\eqn{=T(v=1,v=2,v=3)}.
#'
#' CS proximity regions are constructed with respect to triangle \code{tri}
#' with expansion parameter \eqn{t > 0}
#' and edge regions are based on the center \eqn{M=(m_1,m_2)}
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
#' or reflexivity graph of the CS-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param t A positive real number
#' which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' @param ugraph The type of the graph based on CS-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return Incidence matrix for the underlying or reflexivity graphs
#' of the CS-PCD with vertices
#' being 2D data set, \code{Xp}
#' in the triangle \code{tri} with edge regions based on center \code{M}
#'
#' @seealso \code{\link{inci.mat.undCS}}, \code{\link{inci.mat.undAStri}},
#' \code{\link{inci.mat.undPEtri}}, and \code{\link[pcds]{inci.matCStri}}
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
#' (IM<-inci.mat.undCStri(Xp,Tr,t=1.5,M))
#' pcds::dom.num.greedy(IM)
#' pcds::Idom.num.up.bnd(IM,3)
#'
#' (IM<-inci.mat.undCStri(Xp,Tr,t=1.5,M,ugraph="r"))
#' pcds::dom.num.greedy(IM)
#' pcds::Idom.num.up.bnd(IM,3)
#' #}
#'
#' @export inci.mat.undCStri
inci.mat.undCStri <- function(Xp,tri,t,M=c(1,1,1),
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

  if (!pcds::is.point(t,1) || t<=0)
  {stop('t must be a scalar > 0')}

  if (!(pcds::is.point(M) || pcds::is.point(M,3) ))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates ')}

  if (pcds::dimension(M)==3)
  {M<-pcds::bary2cart(M,tri)}

  if (!(pcds::in.triangle(M,tri,boundary=FALSE)$in.tri))
   {stop('M is not a center in the interior of the triangle')}

  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  n<-nrow(Xp)
  inc.mat<-matrix(0, nrow=n, ncol=n)
 # if (n>1)
 # {
    for (i in 1:n)
    {p1<-Xp[i,]

    for (j in (i:n) )
    {p2<-Xp[j,]
    inc.mat[i,j]<-inc.mat[j,i]<-IedgeCStri(p1,p2,tri,t,M,ugraph)
    }
    }
  #}
  inc.mat
} #end of the function
#'

#################################################################

#' @title The edges of the underlying or reflexivity graphs of
#' the Central Similarity Proximity Catch Digraph
#' (CS-PCD) for 2D data - one triangle case
#'
#' @description
#' An object of class \code{"UndPCDs"}.
#' Returns edges of the underlying or reflexivity graph of CS-PCD
#' as left and right end points
#' and related parameters and the quantities of these graphs.
#' The vertices of these graphs are the data points in \code{Xp}
#' in the multiple triangle case.
#'
#' CS proximity regions are constructed
#' with respect to the triangle \code{tri} with expansion
#' parameter \eqn{t > 0}, i.e.,
#' edges may exist only for points inside \code{tri}.
#' It also provides various descriptions
#' and quantities about the edges of
#' the underlying or reflexivity graphs of the CS-PCD
#' such as number of edges, edge density, etc.
#'
#' Edge regions are based on center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates or \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of
#' the triangle \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' With any interior center \code{M},
#' the edge regions are constructed using the extensions
#' of the lines combining vertices with \code{M}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying or reflexivity graphs of the CS-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param t A positive real number
#' which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' @param ugraph The type of the graph based on CS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the underlying
#' or reflexivity graph of the digraph}
#' \item{parameters}{Parameters of the underlying or reflexivity graph of the digraph,
#' the center \code{M} used to
#' construct the edge regions and the expansion parameter \code{t}.}
#' \item{tess.points}{Tessellation points,
#' i.e., points on which the tessellation of the study region
#' is performed, here, tessellation is the support triangle.}
#' \item{tess.name}{Name of the tessellation points \code{tess.points}}
#' \item{vertices}{Vertices of the underlying
#' or reflexivity graph of the digraph, \code{Xp} points}
#' \item{vert.name}{Name of the data set
#' which constitutes the vertices of
#' the underlying or reflexivity graph of the digraph}
#' \item{LE, RE}{Left and right end points of the edges of
#' the underlying or reflexivity graph of CS-PCD for 2D data set \code{Xp}
#' as vertices of the underlying or reflexivity graph of the digraph}
#' \item{mtitle}{Text for \code{"main"} title
#' in the plot of the underlying or reflexivity graph of the digraph}
#' \item{quant}{Various quantities for the underlying
#' or reflexivity graph of the digraph:
#' number of vertices, number of partition points,
#' number of intervals, number of edges, and edge density.}
#'
#' @seealso \code{\link{edgesCS}}, \code{\link{edgesAStri}}, \code{\link{edgesPEtri}},
#' and \code{\link[pcds]{arcsCStri}}
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
#' t<-1.5
#'
#' #for underlying graph
#' Edges<-edgesCStri(Xp,Tr,t,M)
#' Edges
#' summary(Edges)
#' plot(Edges)
#'
#' #for reflexivity graph
#' Edges<-edgesCStri(Xp,Tr,t,M,ugraph="r")
#' Edges
#' summary(Edges)
#' plot(Edges)
#'
#' #can add edge regions
#' cent<-M
#' cent.name<-"M"
#' Ds<-pcds::prj.cent2edges(Tr,M)
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
#' @export edgesCStri
edgesCStri <- function(Xp,tri,t,M=c(1,1,1),ugraph=c("underlying", "reflexivity"))
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

  if (!pcds::is.point(t,1) || t<=0)
  {stop('t must be a scalar > 0')}

  if (!(pcds::is.point(M) || pcds::is.point(M,3) ))
  {stop('M must be a numeric 2D point for Cartesian coordinates or
  3D point for barycentric coordinates ')}

  if (pcds::dimension(M)==3)
  {M<-pcds::bary2cart(M,tri)}

  if (!(pcds::in.triangle(M,tri,boundary=FALSE)$in.tri))
   {stop('M is not a center in the interior of the triangle')}

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

  #the edges of the underlying or reflexivity graphs of CS-PCDs
  lep<-rep<-NULL #left and right end points for the edges
  if (n2>1)
  {
    for (j in 1:(n2-1))
    {
      p1<-Xtri[j,];
      for (k in (j+1):n2)  #to avoid loops
      {
        p2<-Xtri[k,];
        if (IedgeCStri(p1,p2,tri,t,M,ugraph)==1)
        {
          lep <-rbind(lep,p1); rep <-rbind(rep,p2);
        }
      }
    }
  }

  param<-list(M,t)
  Mr<-round(M,2)
  cname <-"M"
  names(param)<-c("center","expansion parameter")
  main.txt<-paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"), " Graph of CS-PCD with t = ",t," and Center ", cname," = (",Mr[1],",",Mr[2],")",sep="")
  typ<-paste(ifelse(ugraph=="underlying","Underlying", "Reflexivity")," Graph of Central Similarity Proximity Catch Digraph (CS-PCD) for 2D Points in the Triangle with Expansion Parameter t = ",t," and Center ", cname," = (",Mr[1],",",Mr[2],")",sep="")

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

#' @title The plot of the edges of the underlying or reflexivity graphs of
#' the Central Similarity Proximity Catch Digraph
#' (CS-PCD) for 2D data - one triangle case
#'
#' @description Plots the edges of the underlying or reflexivity graphs of
#' the Central Similarity Proximity Catch Digraph
#' (CS-PCD) whose vertices are the data points, \code{Xp}
#' and the triangle \code{tri}.
#' CS proximity regions
#' are constructed with respect to the triangle \code{tri}
#' with expansion parameter \eqn{t > 0},
#' i.e., edges may exist only for \code{Xp} points inside the triangle \code{tri}.
#'
#' Edge regions are based on center \eqn{M=(m_1,m_2)}
#' in Cartesian coordinates
#' or \eqn{M=(\alpha,\beta,\gamma)} in barycentric coordinates
#' in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e.,
#' the center of mass of \code{tri}.
#' With any interior center \code{M},
#' the edge regions are constructed using the extensions
#' of the lines combining vertices with \code{M}.
#'
#' See also (\insertCite{ceyhan:Phd-thesis,ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graphs of the CS-PCD.
#' @param tri A \eqn{3 \times 2} matrix with each row
#' representing a vertex of the triangle.
#' @param t A positive real number
#' which serves as the expansion parameter in CS proximity region.
#' @param M A 2D point in Cartesian coordinates
#' or a 3D point in barycentric coordinates
#' which serves as a center in the interior of the triangle \code{tri};
#' default is \eqn{M=(1,1,1)}, i.e., the center of mass of \code{tri}.
#' @param ugraph The type of the graph based on CS-PCDs,
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
#' @param edge.reg A logical argument to add edge regions to the plot,
#' default is \code{edge.reg=FALSE}.
#' @param \dots	Additional \code{plot} parameters.
#'
#' @return A plot of the edges of the underlying
#' or reflexivity graphs of the CS-PCD
#' whose vertices are the points in data set \code{Xp}
#' and the triangle \code{tri}
#'
#' @seealso \code{\link{plotCSedges}}, \code{\link{plotASedges.tri}},
#' \code{\link{plotPEedges.tri}}, and \code{\link[pcds]{plotCSarcs.tri}}
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
#' t<-1.5
#' plotCSedges.tri(Xp,Tr,t,M,edge.reg = TRUE,xlab="",ylab="")
#' plotCSedges.tri(Xp,Tr,t,M,ugraph="r",edge.reg = TRUE,xlab="",ylab="")
#'
#' #can add vertex labels and text to the figure (with edge regions)
#' Ds<-pcds::prj.cent2edges(Tr,M); cent.name="M"
#'
#' txt<-rbind(Tr,M,Ds)
#' xc<-txt[,1]+c(-.02,.02,.02,.02,.04,-0.03,-.01)
#' yc<-txt[,2]+c(.02,.02,.02,.07,.02,.04,-.06)
#' txt.str<-c("A","B","C",cent.name,"D1","D2","D3")
#' text(xc,yc,txt.str)
#' #}
#'
#' @export plotCSedges.tri
plotCSedges.tri <- function(Xp,tri,t,M=c(1,1,1),ugraph=c("underlying", "reflexivity"),
                               asp=NA,main=NULL,xlab=NULL,ylab=NULL,
                               xlim=NULL,ylim=NULL,edge.reg=FALSE,...)
{
  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  if( any(duplicated(as.data.frame(Xp))) ) #if there are duplicates for Xp values, only one is taken for each
  {Xp = unique(as.data.frame(Xp))
  warning("There were duplicate Xp values;
          only one value is kept for each duplicate Xp value (to avoid edges of zero length)!")}

  EdgesCS<-edgesCStri(Xp,tri,t,M,ugraph)
  lep<-EdgesCS$LE; rep<-EdgesCS$RE #lep, rep are left and right end points of the edges of the graph
  cent = (EdgesCS$param)$c

  Xp<-matrix(Xp,ncol=2)
  if (is.null(xlim))
  {xlim<-range(tri[,1],Xp[,1],cent[1])}
  if (is.null(ylim))
  {ylim<-range(tri[,2],Xp[,2],cent[2])}

  if (is.null(main))
  {Mr=round(cent,2)
  Mvec= paste(Mr, collapse=",")
  main=paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"), " Graph of CS-PCD\n with t = ",t," and M = (",Mvec,")",sep="")
  }

  if (edge.reg)
  {main=c(main,"\n (edge regions added)")}

  plot(Xp,main=main,asp=asp, xlab=xlab, ylab=ylab,
       xlim=xlim,ylim=ylim,pch=".",cex=3,...)
  polygon(tri,...)
  if (!is.null(lep)) {segments(lep[,1], lep[,2], rep[,1], rep[,2], col= 4)}

  if (edge.reg){
    Ds<-pcds::prj.cent2edges(tri,M)
    L<-rbind(cent,cent,cent); R<-Ds
    segments(L[,1], L[,2], R[,1], R[,2], lty=2)
  }
} #end of the function
#'

#################################################################

#' @title The edges of the underlying or reflexivity graphs of
#' the Central Similarity Proximity Catch Digraph
#' (CS-PCD) for 2D data - multiple triangle case
#'
#' @description
#' An object of class \code{"UndPCDs"}.
#' Returns edges of the underlying or reflexivity graph of CS-PCD
#' as left and right end points
#' and related parameters and the quantities of these graphs.
#' The vertices of these graphs are the data points in \code{Xp}
#' in the multiple triangle case.
#'
#' CS proximity regions are
#' defined with respect to the Delaunay triangles
#' based on \code{Yp} points with expansion parameter \eqn{t > 0} and
#' edge regions in each triangle are
#' based on the center \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates
#' in the interior of each Delaunay triangle
#' (default for \eqn{M=(1,1,1)}
#' which is the center of mass of the triangle).
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
#' for more on the CS-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds.ugraph})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graphs of the CS-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param t A positive real number
#' which serves as the expansion parameter in CS proximity region.
#' @param M A 3D point in barycentric coordinates
#' which serves as a center in the interior of each Delaunay triangle,
#' default for \eqn{M=(1,1,1)}
#' which is the center of mass of each triangle.
#' @param ugraph The type of the graph based on CS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return A \code{list} with the elements
#' \item{type}{A description of the underlying
#' or reflexivity graph of the digraph}
#' \item{parameters}{Parameters of
#' the underlying or reflexivity graph of the digraph,
#' the center \code{M} used to
#' construct the edge regions and the expansion parameter \code{t}.}
#' \item{tess.points}{Tessellation points, i.e., points on which the tessellation
#' of the study region is performed, here, tessellation
#' is Delaunay triangulation based on \code{Yp} points.}
#' \item{tess.name}{Name of the tessellation points \code{tess.points}}
#' \item{vertices}{Vertices of the underlying and
#' reflexivity graph of the digraph, \code{Xp} points}
#' \item{vert.name}{Name of the data set
#' which constitute the vertices of
#' the underlying or reflexivity graph of the digraph}
#' \item{LE, RE}{Left and right end points of the edges of
#' the underlying or reflexivity graph of CS-PCD for 2D data set \code{Xp}
#' as vertices of the underlying or reflexivity graph of the digraph}
#' \item{mtitle}{Text for \code{"main"} title
#' in the plot of the underlying or reflexivity graph of the digraph}
#' \item{quant}{Various quantities for
#' the underlying or reflexivity graph of the digraph:
#' number of vertices, number of partition points,
#' number of intervals, number of edges, and edge density.}
#'
#' @seealso \code{\link{edgesCStri}}, \code{\link{edgesAS}}, \code{\link{edgesPE}},
#' and \code{\link[pcds]{arcsCS}}
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
#' t<-1.5
#'
#' Edges<-edgesCS(Xp,Yp,t,M)
#' Edges
#' summary(Edges)
#' plot(Edges)
#'
#' Edges<-edgesCS(Xp,Yp,t,M,ugraph="r")
#' Edges
#' summary(Edges)
#' plot(Edges)
#' #}
#'
#' @export edgesCS
edgesCS <- function(Xp,Yp,t,M=c(1,1,1),ugraph=c("underlying", "reflexivity"))
{
  xname <-deparse(substitute(Xp));
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

  if (!pcds::is.point(t,1) || t<=0)
  {stop('t must be a scalar > 0')}

  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  if (nrow(Yp)==3)
  {
    res<-edgesCStri(Xp,Yp,t,M,ugraph)
  } else
  {
    if ((!pcds::is.point(M,3) ) || !all(M>0))
    {stop('M must be a numeric 3D point with positive barycentric coordinates')}

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
          #convert the triangle Yi.Tri into an nonscaled basic triangle, see as.basic.tri help page

          nl<-nrow(Xl) ; cent<-M
          for (j in 1:(nl-1))
          { for (k in (j+1):nl)  #to avoid loops
            if (IedgeCStri(Xl[j,],Xl[k,],Yi.tri,t,cent,ugraph)==1 )
            {
              lep <-rbind(lep,Xl[j,]); rep <-rbind(rep,Xl[k,]);
            }
          }
        }
      }
    }

    cname <-"M"
    param<-list(M,t)
    names(param)<-c("center","expansion parameter")

    Mvec= paste(M, collapse=",")
    main.txt=paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"), " Graph of CS-PCD\n with t = ",t," and Center ", cname," = (",Mvec,")",sep="")
    typ<-paste(ifelse(ugraph=="underlying","Underlying", "Reflexivity")," Graph of Central Similarity Proximity Catch Digraph (CS-PCD) for 2D points in Multiple Triangles with Expansion parameter t = ",t," and Center ", cname," = (",Mvec,")",sep="")

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
      vertices=Xp, vert.name=xname, #vertices of the underlying or reflexivity graph of the digraph
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

#' @title Incidence matrix for the underlying or reflexivity graphs of
#' Central Similarity Proximity Catch Digraphs (CS-PCDs) -
#' multiple triangle case
#'
#' @description Returns the incidence matrix
#' for the underlying or reflexivity graphs of the CS-PCD
#' whose vertices are the data points in \code{Xp}
#' in the multiple triangle case.
#'
#' CS proximity regions are
#' defined with respect to the Delaunay triangles
#' based on \code{Yp} points with expansion parameter \eqn{t > 0} and
#' edge regions in each triangle are
#' based on the center \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates
#' in the interior of each Delaunay triangle
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
#' for more on the CS-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds.ugraph})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying or reflexivity graphs of the CS-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param t A positive real number
#' which serves as the expansion parameter in CS proximity region.
#' @param M A 3D point in barycentric coordinates
#' which serves as a center in the interior of each Delaunay triangle,
#' default for \eqn{M=(1,1,1)}
#' which is the center of mass of each triangle.
#' @param ugraph The type of the graph based on CS-PCDs,
#' \code{"underlying"} is for the underlying graph,
#' and \code{"reflexivity"} is for
#' the reflexivity graph (default is \code{"underlying"}).
#'
#' @return Incidence matrix for the underlying or reflexivity graphs
#' of the CS-PCD  whose vertices are the 2D data set, \code{Xp}.
#' CS proximity regions are constructed
#' with respect to the Delaunay triangles and \code{M}-edge regions.
#'
#' @seealso \code{\link{inci.mat.undCStri}}, \code{\link{inci.mat.undAS}},
#' \code{\link{inci.mat.undPE}}, and \code{\link[pcds]{inci.matCS}}
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
#' t<-1.5
#'
#' IM<-inci.mat.undCS(Xp,Yp,t,M)
#' IM
#' pcds::dom.num.greedy(IM)
#' #}
#'
#' @export inci.mat.undCS
inci.mat.undCS <- function(Xp,Yp,t,M=c(1,1,1),
                           ugraph=c("underlying", "reflexivity"))
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

  if (!pcds::is.point(t,1) || t<=0)
  {stop('t must be a scalar > 0')}

  ugraph <- match.arg(ugraph)
  if (length(ugraph) > 1 || is.na(ugraph))
    stop("ugraph must be one of \"underlying\" or \"reflexivity\"")

  if (nrow(Yp)==3)
  {
    inc.mat<-inci.mat.undCStri(Xp,Yp,t,M,ugraph)
  } else
  {
    if ((!pcds::is.point(M,3) ) || !all(M>0))
    {stop('M must be a numeric 3D point with positive barycentric coordinates')}

    DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")

    nx<-nrow(Xp)
    ch<-rep(0,nx)
    for (i in 1:nx)
      ch[i]<-interp::in.convex.hull(DTmesh,Xp[i,1],Xp[i,2],strict=FALSE)

    inc.mat<-matrix(0, nrow=nx, ncol=nx)

    DTr<-matrix(interp::triangles(DTmesh)[,1:3],ncol=3)
    nt<-nrow(DTr)  #number of Delaunay triangles

   # if (nx>1)
   # {
      i.tr<-rep(0,nx)  #the vector of indices for the triangles that contain the Xp points
      for (i in 1:nx)
        for (j in 1:nt)
        {
          tri<-Yp[DTr[j,],]
          if (pcds::in.triangle(Xp[i,],tri,boundary=TRUE)$in.tri )
            i.tr[i]<-j
        }

      cent<-M
      for (i in 1:nx)
      { p1<-Xp[i,]
      Yi.tri<-Yp[DTr[1,],]

      if (i.tr[i]!=0)
      {
        Yi.Tri<-Yp[DTr[i.tr[i],],] #vertices of the ith triangle
        Yi.tri<-pcds::as.basic.tri(Yi.Tri)$tri
        #convert the triangle Yi.Tri into an nonscaled basic triangle, see as.basic.tri help page
      }
      for (j in i:nx )
      {p2<-Xp[j,]
      inc.mat[i,j]<-inc.mat[j,i]<-IedgeCStri(p1,p2,Yi.tri,t,cent,ugraph)
      }
      # }
      }
      # }

    #diag(inc.mat)<-1
  }
  inc.mat
} #end of the function
#'

#################################################################

#' @title The plot of the edges of the underlying or reflexivity graphs of
#' the Central Similarity Proximity Catch Digraph
#' (CS-PCD) for 2D data - multiple triangle case
#'
#' @description Plots the edges of the underlying or reflexivity graphs of
#' the Central Similarity Proximity Catch Digraph
#' (CS-PCD) whose vertices are the data
#' points in \code{Xp} in the multiple triangle case
#' and the Delaunay triangles based on \code{Yp} points.
#'
#' CS proximity regions are constructed
#' with respect to the Delaunay triangles based on \code{Yp} points, i.e.,
#' CS proximity regions are defined only for \code{Xp} points
#' inside the convex hull of \code{Yp} points.
#' That is, edges may exist for \code{Xp} points
#' only inside the convex hull of \code{Yp} points.
#'
#' Edge regions in each triangle are
#' based on the center \eqn{M=(\alpha,\beta,\gamma)}
#' in barycentric coordinates in the interior of each Delaunay triangle
#' (default for \eqn{M=(1,1,1)}
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
#' for more on the CS-PCDs.
#' Also, see (\insertCite{okabe:2000,ceyhan:comp-geo-2010,sinclair:2016;textual}{pcds.ugraph})
#' for more on Delaunay triangulation and the corresponding algorithm.
#'
#' @param Xp A set of 2D points
#' which constitute the vertices of the underlying
#' or reflexivity graphs of the CS-PCD.
#' @param Yp A set of 2D points
#' which constitute the vertices of the Delaunay triangles.
#' @param t A positive real number
#' which serves as the expansion parameter in CS proximity region.
#' @param M A 3D point in barycentric coordinates
#' which serves as a center in the interior of each Delaunay triangle,
#' default for \eqn{M=(1,1,1)}
#' which is the center of mass of each triangle.
#' @param ugraph The type of the graph based on CS-PCDs,
#' \code{"underlying"} is for the underlying graph, and \code{"reflexivity"} is for
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
#' or reflexivity graphs of the CS-PCD
#' whose vertices are the points in data set \code{Xp} and the Delaunay
#' triangles based on \code{Yp} points
#'
#' @seealso \code{\link{plotCSedges.tri}}, \code{\link{plotASedges}},
#' \code{\link{plotPEedges}}, and \code{\link[pcds]{plotCSarcs}}
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
#' t<-1.5
#'
#' plotCSedges(Xp,Yp,t,M,xlab="",ylab="")
#' plotCSedges(Xp,Yp,t,M,xlab="",ylab="",ugraph="r")
#' #}
#'
#' @export plotCSedges
plotCSedges <- function(Xp,Yp,t,M=c(1,1,1),ugraph=c("underlying", "reflexivity"),
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
    plotCSedges.tri(Xp,Yp,t,M,ugraph,asp,main,xlab,ylab,xlim,ylim)
  } else
  {
    EdgesCS<-edgesCS(Xp,Yp,t,M,ugraph)
    lep<-EdgesCS$LE; rep<-EdgesCS$RE

    DTmesh<-interp::tri.mesh(Yp[,1],Yp[,2],duplicate="remove")
    Xch<-pcds::Xin.convex.hullY(Xp,Yp)

    if (is.null(main))
    {Mvec= paste(M, collapse=",")
    main=paste("Edges of ", ifelse(ugraph=="underlying","Underlying", "Reflexivity"), " Graph of CS-PCD\n with t = ",t," and M = (",Mvec,")",sep="")
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


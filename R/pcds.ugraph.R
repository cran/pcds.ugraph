#' pcds.ugraph: A package for the Underlying and Reflexivity Graphs of the Proximity Catch Digraphs and Their Applications
#'
#' \code{pcds.ugraph} is a package for construction and visualization of the underlying graphs based on
#' proximity catch digraphs and for computation of edge density of these graphs for testing spatial patterns.
#'
#' The PCD families considered are Arc-Slice PCDs,
#' Proportional-Edge PCDs and Central Similarity PCDs
#' (\insertCite{ceyhan:Phd-thesis,ceyhan:arc-density-PE,ceyhan:arc-density-CS;textual}{pcds.ugraph}).
#'
#' The graph invariant used in testing spatial point data are the edge density of
#' the underlying and reflexivity graphs of the PCDs
#' (see \insertCite{ceyhan:stamet2016;textual}{pcds.ugraph}).
#'
#' The package also contains visualization tools for these graphs for 1D-3D vertices.
#' The AS-PCD and CS-PCD related tools are provided for 1D and 2D data;
#' PE-PCD related tools are provided for 1D-3D data.
#'
#' @section The \code{pcds.ugraph} functions:
#' The \code{pcds.ugraph} functions can be grouped as AS-PCD Functions, PE-PCD Functions,
#' and CS-PCD Functions.
#'
#' @section Arc-Slice PCD Functions:
#' Contains the functions used in AS-PCD construction and computation of edge density of
#' the corresponding underlying and reflexivity graph.
#'
#' @section Proportional-Edge PCD Functions:
#' Contains the functions used in PE-PCD construction and computation of edge density of
#' the corresponding underlying and reflexivity graph.
#'
#' @section Central-Similarity PCD Functions:
#' Contains the functions used in CS-PCD construction and computation of edge density of
#' the corresponding underlying and reflexivity graph.
#'
#' @references
#' \insertAllCited{}
#'
#' @docType package
#' @name pcds.ugraph
NULL

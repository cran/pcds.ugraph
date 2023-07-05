## Startup functions ------------------------------------

#' .onAttach start message
#'
#' @param libname defunct
#' @param pkgname defunct
#'
#' @return invisible()
.onAttach <- function(libname, pkgname) {
  start_message <- c( "     pcds.ugraph: functions for the  \n\n",
                      " underlying and reflexivity graphs of proximity catch digraphs,  \n\n",
                      "  their visualisation and application in spatial data analysis\n\n",
                      "      by Dr. Elvan Ceyhan <elvanceyhan@gmail.com>\n\n"
  )
  packageStartupMessage(start_message)
  invisible()
}
#'

################################################################
#'
#' .onLoad getOption package settings
#'
#' @param libname defunct
#' @param pkgname defunct
#'
#' @return invisible()
#'
#' @examples
#' getOption("pcds.ugraph.name")
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.pcds.ugraph <- list(
    #pcds.path = "~/R-dev",
    pcds.ugraph.install.args  = "",
    pcds.ugraph.name          = "Elvan Ceyhan",
    pcds.ugraph.desc.author   = "Elvan Ceyhan <elvanceyhan@gmail.com> [aut, cre]",
    pcds.ugraph.desc.license  = "GPL-2",
    pcds.ugraph.desc.suggests = NULL,
    pcds.ugraph.desc          = list()
  )
  toset <- !(names(op.pcds.ugraph) %in% names(op))
  if (any(toset)) options(op.pcds.ugraph[toset])

  invisible()
}

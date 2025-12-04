#' My file reader
#'
#' Reads in data from a text file.
#'
#' @param file Path to the file.
#' @return A data frame.
#' @export
myreadx1 <- function(file) {
  read.table(file, header = TRUE)  # or whatever your class requires
}

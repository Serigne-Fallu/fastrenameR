#' Sanitise a character vector for FASTA headers
#'
#' Replaces spaces and punctuation with hyphens, collapses runs of hyphens,
#' removes leading/trailing hyphens, and sets empty strings to "UNKNOWN".
#'
#' @param x Character vector.
#' @param replace Character used as replacement (default "-").
#' @return A cleaned character vector.
#' @export
#'
#' @examples
#' clean_field(c("United Kingdom", "2023-01-15", ""))
clean_field <- function(x, replace = "-") {
  x <- as.character(x)
  # Replace forbidden characters (spaces, punctuation)
  x <- stringr::str_replace_all(x, "[\\s,;:()\\[\\]{}|/\\\\]+", replace)
  # Collapse multiple hyphens
  x <- stringr::str_replace_all(x, paste0(replace, "{2,}"), replace)
  # Trim leading/trailing hyphens
  x <- stringr::str_replace_all(x, paste0("^", replace, "|", replace, "$"), "")
  # Replace empty strings
  x[x == "" | is.na(x)] <- "UNKNOWN"
  x
}
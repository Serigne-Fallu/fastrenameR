#' Rename FASTA headers using a metadata mapping
#'
#' Reads a FASTA file (optionally .gz) and replaces each header with the
#' corresponding new name from the mapping. Sequences without a match are
#' kept unchanged and a warning is emitted.
#'
#' @param fasta_in Path to input FASTA (.gz allowed).
#' @param metadata_map Named character vector: names are old IDs, values new names.
#' @param fasta_out Output file; if NULL, the renamed FASTA is returned as a
#'   character vector (invisibly).
#' @param verbose If TRUE, print renaming progress.
#' @return If `fasta_out` is provided, the output path (invisibly);
#'   otherwise the renamed FASTA lines.
#' @export
rename_fasta <- function(fasta_in, metadata_map, fasta_out = NULL, verbose = TRUE) {
  con <- if (grepl("\\.gz$", fasta_in)) gzfile(fasta_in, "rt") else file(fasta_in, "r")
  on.exit(close(con))

  out_lines <- character(0)
  header_present <- FALSE
  first_seq <- TRUE

  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    if (grepl("^>", line)) {
      original <- sub("^>", "", line)
      original <- trimws(original)
      fasta_id <- strsplit(original, "\\s+")[[1]][1]
      if (fasta_id %in% names(metadata_map)) {
        new_header <- metadata_map[fasta_id]
        out_lines <- c(out_lines, paste0(">", new_header))
        if (verbose) message("Renamed ", fasta_id, " -> ", new_header)
      } else {
        warning("ID not in metadata, kept: ", fasta_id)
        out_lines <- c(out_lines, line)
      }
      header_present <- TRUE
      first_seq <- FALSE
    } else {
      if (!header_present && first_seq && nchar(trimws(line)) > 0) {
        warning("First record lacks '>', adding 'UNKNOWN' header.")
        out_lines <- c(out_lines, ">UNKNOWN")
        header_present <- TRUE
        first_seq <- FALSE
      }
      out_lines <- c(out_lines, line)
    }
  }

  if (!is.null(fasta_out)) {
    writeLines(out_lines, fasta_out)
    message("Renamed FASTA written to ", fasta_out)
    invisible(fasta_out)
  } else {
    invisible(out_lines)
  }
}
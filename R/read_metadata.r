#' Read metadata and build a renaming mapping
#'
#' Reads a TSV (or CSV with given delimiter), identifies an ID column and a set
#' of output columns, then creates a named vector that maps each sequence ID to
#' a new header string formed by joining the cleaned output column values with "_".
#'
#' @param path Path to the metadata file.
#' @param id_col ID column: character (name) or integer (1‑based index).
#' @param out_cols Output columns: vector of names or 1‑based indices.
#' @param delim Field delimiter (default: tab).
#' @return A list with elements `mapping` (named character vector),
#'   `raw_data` (tibble), `id_col_name`, and `out_col_names`.
#' @export
read_metadata_tsv <- function(path, id_col, out_cols, delim = "\t") {
  raw <- readr::read_delim(path, delim = delim, show_col_types = FALSE)
  cols <- names(raw)

  # ID column
  if (is.numeric(id_col)) {
    stopifnot(id_col >= 1, id_col <= ncol(raw))
    id_idx <- id_col
    id_col_name <- cols[id_idx]
  } else {
    if (!id_col %in% cols)
      stop("ID column '", id_col, "' not found. Available: ",
           paste(cols, collapse = ", "))
    id_idx <- match(id_col, cols)
    id_col_name <- id_col
  }

  # Output columns
  out_idx <- sapply(out_cols, function(oc) {
    if (is.numeric(oc)) {
      stopifnot(oc >= 1, oc <= ncol(raw))
      oc
    } else {
      idx <- match(oc, cols)
      if (is.na(idx)) stop("Output column '", oc, "' not found.")
      idx
    }
  })
  out_col_names <- cols[out_idx]

  # Build mapping
  ids <- raw[[id_idx]]
  if (anyDuplicated(ids)) warning("Duplicate IDs in metadata; last kept.")

  new_names <- apply(raw[, out_idx, drop = FALSE], 1, function(row) {
    paste(clean_field(as.character(row)), collapse = "_")
  })
  names(new_names) <- ids

  list(
    mapping      = new_names,
    raw_data     = raw,
    id_col_name  = id_col_name,
    out_col_names = out_col_names
  )
}
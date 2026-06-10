#' @importFrom lubridate leap_year
NULL

#' Fast FASTA header renaming using metadata
#'
#' Reads a metadata TSV, lets the user specify which columns to use for the new
#' FASTA header names based on the metadata header, and renames the input FASTA.
#'
#' @param standardize_dates Logical, whether to standardize date formats in output.
#' @param date_mid_month_day Day for month-year formats (default: 15).
#' @param date_mid_year_day Day for year-only formats (default: 182, July 1).
#' @param metadata Path to metadata TSV file.
#' @param fasta Path to input FASTA (.gz allowed).
#' @param output Output FASTA file. If NULL, uses input name + "_renamed.fasta".
#' @param id_col ID column: name or 1-based index from the metadata header.
#' @param out_cols Output columns: vector of names or indices. These will be
#'   joined with "_" to form the new FASTA header.
#' @param delim Metadata delimiter (default: tab).
#' @param verbose Print progress messages.
#' @param show_cols If TRUE, display available column names before renaming.
#' @param standardize_dates Logical, whether to standardize date formats in output.
#' @param date_mid_month_day Day for month-year formats (default: 15).
#' @param date_mid_year_day Day for year-only formats (default: 182, July 1).
#' @return The output file path (invisibly).
#' @export
#'
#' @examples
#' \dontrun{
#' # Show available columns first
#' fastrenameR::show_columns("metadata.tsv")
#'
#' # Rename with country + date
#' fastrenameR("metadata.tsv", "seqs.fasta",
#'             id_col = "accessionVersion",
#'             out_cols = c("geoLocCountry", "sampleCollectionDate"))
#'
#' # Rename with date standardization
#' fastrenameR("metadata.tsv", "seqs.fasta",
#'             id_col = "accessionVersion",
#'             out_cols = c("geoLocCountry", "sampleCollectionDate"),
#'             standardize_dates = TRUE)
#' }
fastrenameR <- function(metadata, fasta, output = NULL,
                        id_col = NULL, out_cols = NULL,
                        delim = "\t", verbose = TRUE,
                        show_cols = FALSE,
                        standardize_dates = FALSE,
                        date_mid_month_day = 15,
                        date_mid_year_day = 182) {
  
  # --- Lire l'en-tête pour connaître les colonnes disponibles ---
  header_line <- readLines(metadata, n = 1)
  available_cols <- strsplit(header_line, delim)[[1]]
  available_cols <- trimws(available_cols)
  
  # --- Afficher les colonnes si demandé ---
  if (show_cols) {
    show_columns(metadata, delim)
  }
  
  # --- Vérifier que les colonnes sont spécifiées ---
  if (is.null(id_col) || is.null(out_cols)) {
    stop(
      "You must specify id_col and out_cols.\n\n",
      "Examples:\n",
      "  fastrenameR('meta.tsv', 'seqs.fasta',\n",
      "          id_col = 'accessionVersion',\n",
      "          out_cols = c('geoLocCountry', 'sampleCollectionDate'))\n\n",
      "Available columns in metadata:\n  ",
      paste(available_cols, collapse = ", "),
      "\n\nUse show_columns('", metadata, "') to display them with indices."
    )
  }
  
  # --- Valider que les colonnes spécifiées existent ---
  if (is.character(id_col)) {
    if (!id_col %in% available_cols) {
      stop(
        "ID column '", id_col, "' not found.\n",
        "Available columns: ", paste(available_cols, collapse = ", ")
      )
    }
  }
  
  if (is.character(out_cols)) {
    missing_cols <- out_cols[!out_cols %in% available_cols]
    if (length(missing_cols) > 0) {
      stop(
        "Column(s) not found: ", paste(missing_cols, collapse = ", "), "\n",
        "Available columns: ", paste(available_cols, collapse = ", ")
      )
    }
  }
  
  # --- Construire le mapping ---
  if (verbose) cat("Reading metadata and building renaming map...\n")
  
  meta <- read_metadata_tsv(
    path      = metadata,
    id_col    = id_col,
    out_cols  = out_cols,
    delim     = delim
  )
  
  # --- Déterminer le fichier de sortie ---
  if (is.null(output)) {
    base <- tools::file_path_sans_ext(fasta)
    base <- sub("\\.gz$", "", base)
    output <- paste0(base, "_renamed.fasta")
  }
  
  # --- Renommer le FASTA ---
  if (verbose) cat("Renaming FASTA headers...\n")
  
  rename_fasta(
    fasta_in     = fasta,
    metadata_map = meta$mapping,
    fasta_out    = output,
    verbose      = verbose
  )
  
  # --- Standardize dates if requested ---
  if (standardize_dates) {
    if (verbose) cat("\nStandardizing date formats...\n")
    output_with_dates <- sub("\\.fasta$", "_ISOdates.fasta", output)
    output_with_dates <- sub("\\.gz$", "", output_with_dates)
    
    standardize_fasta_dates(
      fasta_in = output,
      fasta_out = output_with_dates,
      mid_month_day = date_mid_month_day,
      mid_year_day = date_mid_year_day,
      verbose = verbose
    )
    
    # Remove non-standardized version if different
    if (file.exists(output) && output != output_with_dates) {
      file.remove(output)
    }
    output <- output_with_dates
  }
  
  if (verbose) cat("\nDone! Renamed FASTA saved to:", output, "\n")
  
  invisible(output)
}


#' Show available columns in a metadata file
#'
#' Reads the header of a metadata TSV and displays column names with their
#' 1-based indices. Use these names or indices in `fastrenameR()`.
#'
#' @param metadata Path to metadata TSV file.
#' @param delim Field delimiter (default: tab).
#' @return Invisibly returns a character vector of column names.
#' @export
#'
#' @examples
#' \dontrun{
#' fastrenameR::show_columns("metadata.tsv")
#' }
show_columns <- function(metadata, delim = "\t") {
  header_line <- readLines(metadata, n = 1)
  cols <- strsplit(header_line, delim)[[1]]
  cols <- trimws(cols)
  
  cat("\n=== Available columns in", basename(metadata), "===\n")
  for (i in seq_along(cols)) {
    cat(sprintf("  %2d  %s\n", i, cols[i]))
  }
  cat("========================================\n\n")
  
  invisible(cols)
}
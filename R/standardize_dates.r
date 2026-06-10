#' Standardize dates in FASTA headers or metadata to ISO format
#'
#' Converts various date formats to ISO standard (YYYY-MM-DD) for use in
#' phylogenetic software like BEAST, BEAUti, or Nextstrain.
#'
#' Supported formats:
#' \itemize{
#'   \item Month-Day-Year: "5-15-2024" or "05-15-2024" -> "2024-05-15"
#'   \item Month-Year: "5-2024" -> "2024-05-15" (mid-month assumed)
#'   \item Decimal year: "2024.5" -> "2024-07-02" (day 183/365)
#'   \item Year only: "2024" -> "2024-07-01" (July 1st assumed)
#' }
#'
#' @param x Character vector of date strings or FASTA headers containing dates
#' @param pattern Regular expression pattern to find dates. Default detects
#'   common date formats (m-d-yyyy, m-yyyy, decimal year, yyyy)
#' @param replacement Logical, if TRUE replace dates within larger strings,
#'   if FALSE return only the converted date part
#' @param mid_month_day Day of month to use for month-year formats (default: 15)
#' @param mid_year_day Day of year to use for year-only formats (default: 182,
#'   July 1st)
#' @param return_na_if_invalid Return NA for invalid dates instead of original
#' @param verbose Print progress messages
#'
#' @return A character vector with standardized dates or modified strings
#' @export
#'
#' @examples
#' # Convert various date formats
#' standardize_dates(c("5-15-2024", "5-2024", "2024.5", "2024"))
#' # Returns: "2024-05-15", "2024-05-15", "2024-07-02", "2024-07-01"
#'
#' # Replace dates in FASTA headers
#' headers <- c(">USA_5-15-2024_sample1", ">France_2024.5_isolate")
#' standardize_dates(headers, replacement = TRUE)
#' # Returns: ">USA_2024-05-15_sample1", ">France_2024-07-02_isolate"
standardize_dates <- function(x,
                              pattern = NULL,
                              replacement = FALSE,
                              mid_month_day = 15,
                              mid_year_day = 182,
                              return_na_if_invalid = FALSE,
                              verbose = TRUE) {
  
  if (verbose) cat("Standardizing date formats...\n")
  
  # Default pattern to find dates in strings
  if (is.null(pattern)) {
    pattern <- "\\d{1,2}-\\d{1,2}-\\d{4}|\\d{1,2}-\\d{4}|\\d{4}\\.\\d+|\\d{4}"
  }
  
  # Core conversion function
  convert_to_iso <- function(date_str, mid_month_day, mid_year_day) {
    # Handle NA or empty
    if (is.na(date_str) || date_str == "") {
      return(if (return_na_if_invalid) NA_character_ else date_str)
    }
    
    # Month-Day-Year (e.g., 5-15-2024 or 05-15-2024)
    if (grepl("^\\d{1,2}-\\d{1,2}-\\d{4}$", date_str)) {
      parts <- strsplit(date_str, "-")[[1]]
      month <- sprintf("%02d", as.integer(parts[1]))
      day <- sprintf("%02d", as.integer(parts[2]))
      year <- parts[3]
      return(paste0(year, "-", month, "-", day))
    }
    
    # Month-Year (e.g., 5-2024)
    if (grepl("^\\d{1,2}-\\d{4}$", date_str)) {
      parts <- strsplit(date_str, "-")[[1]]
      month <- sprintf("%02d", as.integer(parts[1]))
      year <- parts[2]
      # Use mid-month day
      day <- sprintf("%02d", mid_month_day)
      return(paste0(year, "-", month, "-", day))
    }
    
    # Decimal year (e.g., 2024.5)
    if (grepl("^\\d{4}\\.\\d+$", date_str)) {
      year_val <- floor(as.numeric(date_str))
      frac <- as.numeric(date_str) - year_val
      
      days_in_year <- ifelse(lubridate::leap_year(year_val), 366, 365)
      day_of_year <- round(frac * days_in_year)
      
      # Ensure within bounds
      day_of_year <- max(1, min(day_of_year, days_in_year))
      
      date_obj <- as.Date(paste0(year_val, "-01-01")) + (day_of_year - 1)
      return(format(date_obj, "%Y-%m-%d"))
    }
    
    # Year only (e.g., 2024)
    if (grepl("^\\d{4}$", date_str)) {
      year <- date_str
      # Use mid-year day
      date_obj <- as.Date(paste0(year, "-01-01")) + (mid_year_day - 1)
      return(format(date_obj, "%Y-%m-%d"))
    }
    
    # No matching format
    if (return_na_if_invalid) {
      return(NA_character_)
    } else {
      return(date_str)
    }
  }
  
  # Apply conversion
  if (replacement) {
    # Replace dates within larger strings using base R gsub
    result <- sapply(x, function(string) {
      dates <- stringr::str_extract_all(string, pattern)[[1]]
      if (length(dates) == 0) return(string)
      
      for (date in dates) {
        iso_date <- convert_to_iso(date, mid_month_day, mid_year_day)
        if (!is.na(iso_date) && iso_date != date) {
          # Use fixed = TRUE for literal string replacement
          string <- gsub(date, iso_date, string, fixed = TRUE)
        }
      }
      string
    }, USE.NAMES = FALSE)
  } else {
    # Extract and convert just the date parts
    result <- sapply(x, function(string) {
      dates <- stringr::str_extract_all(string, pattern)[[1]]
      if (length(dates) == 0) {
        return(if (return_na_if_invalid) NA_character_ else string)
      }
      # Convert first date found (most common case)
      convert_to_iso(dates[1], mid_month_day, mid_year_day)
    }, USE.NAMES = FALSE)
  }
  
  if (verbose) {
    n_converted <- sum(result != x, na.rm = TRUE)
    cat(sprintf("Converted %d date(s) to ISO format\n", n_converted))
  }
  
  result
}

#' Apply date standardization to FASTA file
#'
#' Reads a FASTA file, standardizes all dates in headers, and writes a new file.
#' This is particularly useful for preparing sequences for BEAST/BEAUti analysis.
#'
#' @param fasta_in Path to input FASTA file (.gz allowed)
#' @param fasta_out Output FASTA file path. If NULL, adds "_ISOdates" suffix
#' @param pattern Regular expression pattern to find dates
#' @param mid_month_day Day of month for month-year formats (default: 15)
#' @param mid_year_day Day of year for year-only formats (default: 182, July 1)
#' @param verbose Print progress messages
#'
#' @return Output file path (invisibly)
#' @export
#'
#' @examples
#' \dontrun{
#' standardize_fasta_dates("sequences.fasta", "sequences_BEASTready.fasta")
#' }
standardize_fasta_dates <- function(fasta_in,
                                    fasta_out = NULL,
                                    pattern = NULL,
                                    mid_month_day = 15,
                                    mid_year_day = 182,
                                    verbose = TRUE) {
  
  if (verbose) cat("Reading FASTA file:", fasta_in, "\n")
  
  # Check if Biostrings is available
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package is required for this function. ",
         "Please install it with: BiocManager::install('Biostrings')")
  }
  
  # Read sequences
  if (grepl("\\.gz$", fasta_in)) {
    # Handle gzipped files
    temp_file <- tempfile(fileext = ".fasta")
    system(paste("gunzip -c", shQuote(fasta_in), ">", shQuote(temp_file)))
    sequences <- Biostrings::readDNAStringSet(temp_file)
    unlink(temp_file)
  } else {
    sequences <- Biostrings::readDNAStringSet(fasta_in)
  }
  
  if (verbose) cat(sprintf("Read %d sequences\n", length(sequences)))
  
  # Standardize headers
  old_names <- names(sequences)
  new_names <- standardize_dates(old_names,
                                 pattern = pattern,
                                 replacement = TRUE,
                                 mid_month_day = mid_month_day,
                                 mid_year_day = mid_year_day,
                                 verbose = verbose)
  
  names(sequences) <- new_names
  
  # Determine output file
  if (is.null(fasta_out)) {
    base <- tools::file_path_sans_ext(fasta_in)
    base <- sub("\\.gz$", "", base)
    fasta_out <- paste0(base, "_ISOdates.fasta")
  }
  
  # Write output
  Biostrings::writeXStringSet(sequences, fasta_out)
  
  if (verbose) {
    n_changed <- sum(old_names != new_names, na.rm = TRUE)
    cat(sprintf("Changed %d sequence headers\n", n_changed))
    cat("Output written to:", fasta_out, "\n")
  }
  
  invisible(fasta_out)
}

#' Enhance fastrenameR output with date standardization
#'
#' A wrapper that runs fastrenameR then standardizes dates in the output.
#' This combines metadata-based renaming with date standardization in one step.
#'
#' @param metadata Path to metadata TSV file.
#' @param fasta Path to input FASTA (.gz allowed).
#' @param output Output FASTA file. If NULL, uses input name + "_renamed.fasta".
#' @param id_col ID column: name or 1-based index from the metadata header.
#' @param out_cols Output columns: vector of names or indices. These will be
#'   joined with "_" to form the new FASTA header.
#' @param delim Metadata delimiter (default: tab).
#' @param verbose Print progress messages.
#' @param show_cols If TRUE, display available column names before renaming.
#' @param standardize_dates Logical, whether to standardize dates after renaming
#' @param date_pattern Pattern for finding dates (passed to standardize_dates)
#' @param mid_month_day Day for month-year formats (default: 15)
#' @param mid_year_day Day for year-only formats (default: 182, July 1)
#' @param ... Additional arguments passed to fastrenameR
#'
#' @return Output file path (invisibly)
#' @export
#'
#' @examples
#' \dontrun{
#' fastrenameR_with_dates("metadata.tsv", "seqs.fasta",
#'                        id_col = "accessionVersion",
#'                        out_cols = c("geoLocCountry", "sampleCollectionDate"))
#' }
fastrenameR_with_dates <- function(metadata,
                                   fasta,
                                   output = NULL,
                                   id_col = NULL,
                                   out_cols = NULL,
                                   delim = "\t",
                                   verbose = TRUE,
                                   show_cols = FALSE,
                                   standardize_dates = TRUE,
                                   date_pattern = NULL,
                                   mid_month_day = 15,
                                   mid_year_day = 182,
                                   ...) {
  
  # First run normal fastrenameR
  temp_output <- if (is.null(output)) {
    base <- tools::file_path_sans_ext(fasta)
    base <- sub("\\.gz$", "", base)
    paste0(base, "_temp_renamed.fasta")
  } else {
    output
  }
  
  # Call original fastrenameR
  fastrenameR(
    metadata = metadata,
    fasta = fasta,
    output = if (!standardize_dates && !is.null(output)) output else temp_output,
    id_col = id_col,
    out_cols = out_cols,
    delim = delim,
    verbose = verbose,
    show_cols = show_cols,
    standardize_dates = FALSE,
    ...
  )
  
  # Standardize dates if requested
  if (standardize_dates) {
    final_output <- if (is.null(output)) {
      base <- tools::file_path_sans_ext(fasta)
      base <- sub("\\.gz$", "", base)
      paste0(base, "_renamed_ISOdates.fasta")
    } else {
      output_parts <- strsplit(output, "\\.")[[1]]
      if (length(output_parts) > 1) {
        ext <- paste0(".", output_parts[length(output_parts)])
        base <- paste(output_parts[-length(output_parts)], collapse = ".")
        final_output <- paste0(base, "_ISOdates", ext)
      } else {
        final_output <- paste0(output, "_ISOdates.fasta")
      }
    }
    
    standardize_fasta_dates(
      fasta_in = temp_output,
      fasta_out = final_output,
      pattern = date_pattern,
      mid_month_day = mid_month_day,
      mid_year_day = mid_year_day,
      verbose = verbose
    )
    
    if (temp_output != final_output && file.exists(temp_output)) {
      unlink(temp_output)
    }
    
    if (verbose) cat("\nComplete! Final file with standardized dates:", final_output, "\n")
    invisible(final_output)
  } else {
    invisible(temp_output)
  }
}

#' Validate and report date formats in metadata or FASTA headers
#'
#' Scans dates and reports which are already ISO format and which need conversion.
#'
#' @param x Character vector of dates or FASTA headers
#' @param pattern Pattern to find dates (optional)
#' @param report_format Output format: "summary", "detailed", or "both"
#'
#' @return Depending on \code{report_format}: 
#'   \itemize{
#'     \item "summary": invisible summary statistics
#'     \item "detailed": data frame with detailed information
#'     \item "both": list containing both summary and details
#'   }
#' @export
#'
#' @examples
#' dates <- c("2024-05-15", "5-15-2024", "2024.5", "2024", "invalid")
#' validate_dates(dates, report_format = "summary")
validate_dates <- function(x,
                           pattern = NULL,
                           report_format = "summary") {
  
  if (is.null(pattern)) {
    pattern <- "\\d{1,2}-\\d{1,2}-\\d{4}|\\d{1,2}-\\d{4}|\\d{4}\\.\\d+|\\d{4}"
  }
  
  # Extract dates
  dates <- stringr::str_extract_all(x, pattern)
  dates <- unlist(lapply(dates, function(d) if (length(d) > 0) d[1] else NA))
  
  # Improved ISO detection
  is_iso <- grepl("^\\d{4}-\\d{2}-\\d{2}$", dates) & 
    !is.na(as.Date(dates, format = "%Y-%m-%d"))
  
  is_month_day_year <- grepl("^\\d{1,2}-\\d{1,2}-\\d{4}$", dates) & !is_iso
  is_month_year <- grepl("^\\d{1,2}-\\d{4}$", dates) & !is_iso & !is_month_day_year
  is_decimal <- grepl("^\\d{4}\\.\\d+$", dates) & !is_iso
  is_year_only <- grepl("^\\d{4}$", dates) & !is_iso & !is_decimal
  
  summary_stats <- list(
    total = length(dates),
    iso_format = sum(is_iso, na.rm = TRUE),
    month_day_year = sum(is_month_day_year, na.rm = TRUE),
    month_year = sum(is_month_year, na.rm = TRUE),
    decimal_year = sum(is_decimal, na.rm = TRUE),
    year_only = sum(is_year_only, na.rm = TRUE),
    missing_or_invalid = sum(is.na(dates) | (!is_iso & !is_month_day_year &
                                               !is_month_year & !is_decimal & !is_year_only))
  )
  
  summary_stats$needs_conversion <- summary_stats$total -
    summary_stats$iso_format -
    summary_stats$missing_or_invalid
  
  if (report_format == "summary") {
    cat("\n=== Date Format Validation Report ===\n")
    cat(sprintf("Total dates scanned: %d\n", summary_stats$total))
    cat(sprintf("Already ISO format (YYYY-MM-DD): %d\n", summary_stats$iso_format))
    cat(sprintf("Needs conversion: %d\n", summary_stats$needs_conversion))
    cat(sprintf("  - Month-Day-Year: %d\n", summary_stats$month_day_year))
    cat(sprintf("  - Month-Year: %d\n", summary_stats$month_year))
    cat(sprintf("  - Decimal year: %d\n", summary_stats$decimal_year))
    cat(sprintf("  - Year only: %d\n", summary_stats$year_only))
    cat(sprintf("Missing/Invalid: %d\n", summary_stats$missing_or_invalid))
    cat("=====================================\n\n")
    return(invisible(summary_stats))
    
  } else if (report_format == "detailed") {
    format_type <- ifelse(is_iso, "ISO (YYYY-MM-DD)",
                          ifelse(is_month_day_year, "Month-Day-Year",
                                 ifelse(is_month_year, "Month-Year",
                                        ifelse(is_decimal, "Decimal year",
                                               ifelse(is_year_only, "Year only",
                                                      "Missing/Invalid")))))
    
    details <- data.frame(
      original = dates,
      is_iso = is_iso,
      format_type = format_type,
      stringsAsFactors = FALSE
    )
    return(details)
    
  } else {
    format_type <- ifelse(is_iso, "ISO (YYYY-MM-DD)",
                          ifelse(is_month_day_year, "Month-Day-Year",
                                 ifelse(is_month_year, "Month-Year",
                                        ifelse(is_decimal, "Decimal year",
                                               ifelse(is_year_only, "Year only",
                                                      "Missing/Invalid")))))
    
    details <- data.frame(
      original = dates,
      format_type = format_type,
      stringsAsFactors = FALSE
    )
    return(list(
      summary = summary_stats,
      details = details
    ))
  }
}
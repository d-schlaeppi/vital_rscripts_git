
# Function to copy your contrast matrix to R 
# For regular tables just use the classic read_clip_tbl() function
# This function helps to import your excel contrast matrix and then
#     1. transforms it to the right format to be used for manual pairwise comparisons | contrast_matrix will be loaded in the environment
#     2. it prints the correct syntax for you to just copy form the console to the R script for futer use


copy_paste_contrast_matrix <- function() {
  if (!requireNamespace("clipr", quietly = TRUE)) {
    install.packages("clipr")
  }
  if (!requireNamespace("crayon", quietly = TRUE)) {
    install.packages("crayon")
  }
  library(crayon)
  library(clipr)
  tryCatch({
    contrast_matrix <- read_clip_tbl(header = FALSE, row.names = NULL)
    if (ncol(contrast_matrix) < 2) {
      stop("The clipboard data does not have enough columns.")
    }
    row_names <- contrast_matrix[, 1] # Extract first column as row names
    contrast_matrix <- as.matrix(sapply(contrast_matrix[, -1], as.numeric)) # Convert rest to a numeric matrix
    row.names(contrast_matrix) <- row_names # reassign row names and reset col names
    colnames(contrast_matrix) <- paste0("[,", seq_len(ncol(contrast_matrix)), "]")
    
    # Generate rbind syntax
    row_names_rbind_syntax <- rownames(contrast_matrix)
    formatted_rows <- vector("character", nrow(contrast_matrix))
    for (i in seq_len(nrow(contrast_matrix))) {
      formatted_rows[i] <- paste0("\"", row_names_rbind_syntax[i], "\"=c(", paste(contrast_matrix[i, ], collapse = ","), ")")
    }
    formatted_matrix <- paste(formatted_rows, collapse = ",\n")
    rbind_syntax <- paste("contrast_matrix <- rbind(\n", formatted_matrix, "\n)", sep = "")
    # output
    
    if (!is.null(contrast_matrix)) {
      assign("contrast_matrix", contrast_matrix, envir = .GlobalEnv)
    }
    
    cat(rbind_syntax, "\n")
    return(contrast_matrix)

  }, error = function(e) {
    cat(red("Error:\n"))
    cat("Unable to create the contrast matrix. Please ensure the clipboard data is in the correct format.\n")
    cat("Before running the function, make sure to: \n Highlight & copy your contrast matrix from Excel to the clipboard including row names but no column names.\n")
    cat("Details: ", e$message, "\n")
    return(NULL)
  })
}

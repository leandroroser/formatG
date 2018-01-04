
#' formatG class
#' @name formatG-class
#' @keywords internal
#' @aliases ecogen-class
#'
setClass( "formatG", representation( pointer = "externalptr" ) )

#' formatG initializer
#' @keywords internal

setMethod( "initialize", "formatG", function(.Object, data, ploidy, default_levels,
                                             obj_ncol, char_per_allele, has_rownames,
                                             has_colnames, colnames_vector) {

  mode(data) <- "character"
  .Object@pointer  <- .Call(formatG_method("new"), data, ploidy, default_levels,
                            obj_ncol, char_per_allele, has_rownames, has_colnames,
                            colnames_vector)
  .Object
  } )

#' formatG_method
#' @description Method constructor
#' @param name Name of the method
#' @keywords internal

formatG_method <- function (name) {
  paste ("_formatG_formatG", name, sep = "__")
}

#' $
#' @description Method accessor
#' @param x formatG object
#' @param name Name of the method
#' @keywords internal
#' @export

# syntactic sugar to allow object$method( ... )
setMethod("$", "formatG", function (x, name) {
  function (...) .Call (formatG_method(name) ,x@pointer, ...)
  } )


#' formatG
#' @description formatG constructor
#' @param input Matrix with data
#' @param ploidy data ploidy
#' @export

setClassUnion("matrixORmissing", c("matrix", "missing"))

setGeneric("formatG", function(input = matrix(0, 0, 0),
                               ploidy,
                               has_rownames = TRUE,
                               has_colnames = TRUE,
                               default_levels = character(0),
                               obj_ncol = OL,
                               char_per_allele = 0L,
                               colnames_vector = character(0)
                               ) standardGeneric("formatG"))

setMethod("formatG", "matrixORmissing", function(input = matrix(0, 0, 0),
                                                 ploidy,
                                                 has_rownames = TRUE,
                                        has_colnames = TRUE,
                                        default_levels = character(0),
                                        obj_ncol = 0L,
                                        char_per_allele = 0L,
                                        colnames_vector = character(0)) {


  if(ploidy <= 0L) {
    stop("ploidy must be > 0")
  }

  if(length(input) == 0L && (obj_ncol == 0L || char_per_allele == 0L)) {
    stop("length(input) == 0, requires obj_ncol and char_per_allele values > 0")
  }

  if(length(input) != 0) {
   if(has_rownames && is.null(rownames(input))){
    stop("Table without rown names. Set has_rownames = FALSE or provide a table with row names")
  } else if(has_colnames && is.null(colnames(input))){
    stop("Table without col names. Set has_colnames = FALSE or provide a table with col names")
  }
  out <- new("formatG", input, ploidy, default_levels, 0L, 0L, has_rownames, has_colnames, character(0))

  } else {
    if(obj_ncol == 0L || char_per_allele == 0L) {
      stop("obj_ncol and char_per_allele must be != 0")
    }
    out <- new("formatG", input, ploidy, default_levels, obj_ncol, char_per_allele, FALSE, FALSE, colnames_vector)
  }
  out
})

#' process chunk
#' @description Split locus in column into alleles in columns
#' @param input individuals x 1 StringMatrix, with locus in column
#' @export

setGeneric("process_chunk", function(obj, input_data, fun_name, has_rownames = TRUE,
                                     token = ".",
                                     NA_action = "NA", howmuch = 0,
                                     where = "last", what = "0",
                                     input_or_output = c("data_input", "dataS", "dataI"),
                                     set_ncol = 0) standardGeneric("process_chunk"))
setMethod("process_chunk", "formatG",
          function(obj, input_data,
                   fun_name = c("loci_to_alleles",
                           "loci_token_to_alleles",
                           "to_numeric",
                           "alleles_to_dummy",
                           " alleles_to_loci",
                           " dummy_to_alleles",
                           "add_chars",
                           "StringM_2_IntM"),
                   has_rownames,
                   token = ".", NA_action = "NA", howmuch = 0,
                   where = "last", what = "0",
                   input_or_output = c("data_input", "dataS", "dataI"),
                   set_ncol = 0
                   ){

  fun <- match.arg(fun_name, several.ok = TRUE)
  input_or_output <- match.arg(input_or_output)

  fun_args <-  c("loci_to_alleles",
                 "loci_token_to_alleles",
                 "to_numeric",
                 "alleles_to_dummy",
                 " alleles_to_loci",
                 " dummy_to_alleles",
                 "add_chars",
                 "StringM_2_IntM")

  fun_value <- match(fun, fun_args) - 1

  if(input_or_output == "data_input" || input_or_output == "dataS") {
    fake_matrix <- matrix(0,0,0)
    mode(fake_matrix) <- "integer"
    .Call(formatG_method("process_chunk"), obj@pointer, input_data,
          fake_matrix, fun_value, has_rownames, token, NA_action,
          howmuch, where, what, input_or_output, set_ncol)
  } else {
    fake_matrix <- matrix(0,0,0)
    mode(fake_matrix) <- "character"
    .Call(formatG_method("process_chunk"), obj@pointer, fake_matrix,
          input_data, fun_value, has_rownames, token, NA_action,
          howmuch, where, what, input_or_output, set_ncol)
  }
   })


#' loci_to_alleles
#' @description Split locus in column into alleles in columns
#' @param input individuals x 1 StringMatrix, with locus in column
#' @export

setGeneric("loci_to_alleles", function(obj) standardGeneric("loci_to_alleles"))
setMethod("loci_to_alleles", "formatG", function(obj) {
  .Call(formatG_method("loci_to_alleles"), obj@pointer)
})


#' to_numeric
#' @description Recodes a genetic matrix into numeric format
#' @param obj formatG object
#' @export

setGeneric("to_numeric", function(obj) standardGeneric("to_numeric"))
setMethod("to_numeric", "formatG", function(obj) {
  .Call(formatG_method("to_numeric"), obj@pointer)
})

#'  loci_token_to_alleles
#' @description Split locus in column with alleles separated by a character into alleles in columns
#' @param obj formatG object
#' @export

setGeneric("loci_token_to_alleles", function(obj, token = ".") standardGeneric("loci_token_to_alleles"))
setMethod("loci_token_to_alleles", "formatG", function(obj, token = ".") {
  .Call(formatG_method("loci_token_to_alleles"), obj@pointer, token)
})

#' alleles_to_dummy
#' @description Converts a matrix of alleles into a matrix of dummy factors
#' @param obj formatG object
#' @param NA_action The NA action ("0" or "NA") for recoding NAs.
#' @export

#' @export alleles_to_dummy"
setGeneric("alleles_to_dummy", function(obj, NA_action = "NA") standardGeneric("alleles_to_dummy"))
setMethod("alleles_to_dummy", "formatG", function(obj, NA_action = "NA") {
  .Call(formatG_method("alleles_to_dummy"), obj@pointer, NA_action)
})

 #' formatG constructor
 #' @export get_input
 setGeneric("get_input", function(obj) standardGeneric("get_input"))
 setMethod("get_input", "formatG", function(obj) {
   .Call(formatG_method("get_input"), obj@pointer)
 })

#' formatG constructor
#' @export get_output
setGeneric("get_output", function(obj) standardGeneric("get_output"))
setMethod("get_output", "formatG", function(obj) {
 # this_dim  <- .Call(formatG_method("get_dim"), obj@pointer)
#  structure(.Call(formatG_method("get_output_S"), obj@pointer),dim = this_dim)
  .Call(formatG_method("get_output"), obj@pointer)
})

#' add_chars
#' @description Add characters at the beggining/end of each cell of the matrix
#' referenced in a formatG object
#' @param obj formatG object
#' @param howmuch Number of characters
#' @param where "start" to add characters at the beginning of the string,
#' "end" to add characters at the end of the string
#' @param what Character to add
#' @export

#' formatG constructor
#' @export add_chars
setGeneric("add_chars", function(obj, howmuch = 0, where = 'fist', what = "0") standardGeneric("add_chars"))
setMethod("add_chars", "formatG", function(obj, howmuch = 0, where = 'fist', what = "0") {
  .Call(formatG_method("add_chars"), obj@pointer)
})

#' alleles_to_loci
#' @description Creates a individuals x 1 matrix,  with a loci in a column,
#' from a two-column matrix with the alleles
#' @param obj formatG object

#' @export alleles_to_loci
setGeneric("alleles_to_loci", function(obj, sep = "") standardGeneric("alleles_to_loci"))
setMethod("alleles_to_loci", "formatG", function(obj, sep = "") {
  .Call(formatG_method("alleles_to_loci"), obj@pointer, sep)
})


#' StringM_2_IntM
#' @description Converts StringMatrix into a IntegerMatrix
#' @param obj formatG object

#' @export StringM_2_IntM
setGeneric("StringM_2_IntM", function(obj) standardGeneric("StringM_2_IntM"))
setMethod("StringM_2_IntM", "formatG", function(obj) {
  .Call(formatG_method("StringM_2_IntM"), obj@pointer)
})

#' formatG constructor
#' @export dummy_to_alleles
#' @param obj formatG object

setGeneric("dummy_to_alleles", function(obj) standardGeneric("dummy_to_alleles"))
setMethod("dummy_to_alleles", "formatG", function(obj) {
  .Call(formatG_method("dummy_to_alleles_matrix"), obj@pointer)
})



#' delete_formatG
#' @description delete a formatG object and frees memory
#' @param obj formatG object
#' @export

setGeneric("delete_formatG", function(obj) standardGeneric("delete_formatG"))
setMethod("delete_formatG", "formatG", function(obj) {
  .Call(formatG_method("delete_formatG"), obj@pointer)
})


#' @export get_levels
setGeneric("get_levels", function(obj) standardGeneric("get_levels"))
setMethod("get_levels", "formatG", function(obj) {
  .Call(formatG_method("get_levels"), obj@pointer)
})

#' @export get_alleles_number
setGeneric("get_alleles_number", function(obj) standardGeneric("get_alleles_number"))
setMethod("get_alleles_number", "formatG", function(obj) {
  .Call(formatG_method("get_alleles_number"), obj@pointer)
})


#' formatG constructor
#' @export get_alleles_vector
setGeneric("get_alleles_vector", function(obj) standardGeneric("get_alleles_vector"))
setMethod("get_alleles_vector", "formatG", function(obj) {
  .Call(formatG_method("get_alleles_vector"), obj@pointer)
})


#' formatG constructor
#' @export get_total_alleles
setGeneric("get_total_alleles", function(obj) standardGeneric("get_total_alleles"))
setMethod("get_total_alleles", "formatG", function(obj) {
  .Call(formatG_method("get_total_alleles"), obj@pointer)
})



#' formatter_chunk
#' @export formatter_chunk
formatter_chunk <- function(dbR6_object, path, formatG_input = NULL, ploidy, format_fun, out_df,
                            has_rownames = TRUE, has_colnames = TRUE,
                            chunksize = 1000L,  token = ".", NA_action = "NA", sep_loc = " ", howmuch = 0,
                            where = "last", what = "0",
                            return_formatG = FALSE,
                            input_or_output = c("data_input", "dataS", "dataI"),
                            set_ncol = 0,
                            ...) {
  input_or_output <- match.arg(input_or_output)
  lines_processed <- 0
  if(dbR6_object$exist_table(out_df)) dbR6_object$remove_table(out_df) # remove table if exists

  init <- FALSE
  first_lines <- TRUE

  stream <- reader(path, sep = sep_loc, has_colnames = has_colnames , has_rownames = has_rownames, chunksize = chunksize)
  if(is.null(formatG_input)) {
    init <- TRUE
  } else {
    # clear formatG object
    .Call(formatG_method("reset_state"), formatG_input@pointer, TRUE, set_ncol)
  }

  while(next_chunk(stream)) {
    this_chunk <- get_data(stream)
    if(init) {
      formatG_input <- formatG(this_chunk[1, , drop = FALSE], ploidy =  ploidy,
                          has_rownames = has_rownames,  has_colnames = has_colnames)
      init <- FALSE
    }

    process_chunk(formatG_input , input_data = this_chunk,  fun = format_fun, has_rownames = has_rownames,
                  token = token, NA_action = NA_action, howmuch = howmuch,
                  where = where, what = what, input_or_output = input_or_output,
                  set_ncol = set_ncol)
    this_processed <- reader::matrix2df(get_output(formatG_input ))

    if(first_lines) {
      dbR6_object$add_table(out_df, new_df = this_processed, overwrite = TRUE)
      first_lines <- FALSE
    } else {
      dbR6_object$add_table(out_df, new_df = this_processed, append = TRUE)
    }

    lines_processed <- lines_processed  + nrow(this_chunk)
    cat("Processed ", lines_processed, " rows\n")
  }
  if(return_formatG){
    return(formatG_input)
  }

  rm(formatG_input)
  gc()
  invisible(NULL)
}

#' formatter_chunk_from_db
#' @export formatter_chunk
#'
formatter_chunk_from_db <- function(dbR6_object, table, formatG_input = NULL, ploidy, format_fun, out_df,
                                    has_rownames = TRUE, chunksize = 1000L,
                                    token = ".", NA_action = "NA", howmuch = 0,
                                    where = "last", what = "0", default_levels = character(0),
                                    return_formatG = FALSE,
                                    input_or_output = c("data_input", "dataS", "dataI"),
                                    set_ncol = 0,
                                    ...) {

  input_or_output <- match.arg(input_or_output)
  lines_processed <- 0
  init <- FALSE
  first_lines <- TRUE

  if(dbR6_object$exist_table(out_df)) dbR6_object$remove_table(out_df) # remove table if exists

  if(is.null(formatG_input)) {
    init <- TRUE
  }  else {
    # clear formatG object
    .Call(formatG_method("reset_state"), formatG_input@pointer, TRUE, set_ncol)
  }

  while(TRUE) {
    this_chunk <- dbR6_object$send_query(paste0("SELECT * FROM ", table, " LIMIT ", chunksize, " OFFSET ", lines_processed))

    if(colnames(this_chunk)[1] == "row_names") {
    this_chunk <- dbR6:::as_table_with_rownames(this_chunk)
    }

    this_chunk <- as.matrix(this_chunk)

    if(nrow(this_chunk) == 0) break
    if(init) {
      formatG_input <- formatG(this_chunk[1, , drop = FALSE], ploidy = ploidy, has_rownames = has_rownames, has_colnames = TRUE,
                          default_levels)
      init <- FALSE
    }
    process_chunk(formatG_input, input_data = this_chunk, fun = format_fun, has_rownames = has_rownames,
                  token = token, NA_action = NA_action, howmuch = howmuch,
                  where = where, what = what,  input_or_output = input_or_output, set_ncol = set_ncol)
    this_processed <- reader::matrix2df(get_output(formatG_input))
    if(first_lines) {
    dbR6_object$add_table(out_df, new_df = this_processed, overwrite = TRUE)
      first_lines <- FALSE
    } else {
    dbR6_object$add_table(out_df, new_df = this_processed, append = TRUE)
    }
    lines_processed <- lines_processed  + nrow(this_chunk)
    cat("Processed ", lines_processed, " rows\n")
  }
  if(return_formatG){
    return(formatG_input)
  }

  rm(formatG_input)
  gc()
  invisible(NULL)
}

#' @export get_formatG_parameters
setGeneric("get_formatG_parameters", function(obj) standardGeneric("get_formatG_parameters"))
setMethod("get_formatG_parameters", "formatG", function(obj) {
  .Call(formatG_method("get_formatG_parameters"), obj@pointer)
 })

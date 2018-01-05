
setClass("dbR6S4", slots = c(db = "dbR6"), contains = c("VIRTUAL"))
setClass("genetic_db", slots = c(ATTR = "list"), contains = c("dbR6S4"))

setMethod("initialize", "genetic_db",
          function(.Object, ...) {
            .Object@db<- dbR6$new()
            .Object@db$add_table("G", new_df = data.frame(0))
            .Object@db$add_table("A", new_df = data.frame(0))
            .Object@ATTR = list(names = character(0),   whereIs = new.env(emptyenv()), 
                                .call = call("."),  na_omit = logical(0), 
                                which_NA = integer(0), ploidy = integer(0),
                                char_per_allele = integer(0), n_loci = integer(0),
                                map = list(), alleles_number = integer(0),
                                alleles_vector = integer(0), total_alleles = integer(0),
                                user_levels = character(0), sep = character(0),
                                missing = character(0), type = character(0),
                                NA.char = character(0), rm.empty.ind = logical(0),
                                poly.level = integer(0))
            .Object
          }
)

setMethod("show", "genetic_db",
          function(object) {
            callNextMethod()
          }
)

#' Generate a genetic database

setGeneric("genetic_db", function(G = data.frame(), 
                  type = c("codominant", "dominant"),
                  ploidy = 2,
                  sep = "", 
                  n_loci = 2,
                  ncod = 1,
                  missing = c("0", "NA"),
                  NA.char = "NA", 
                  poly.level = 5,
                  default_levels = character(0),
                  chunksize = 1000,
                  has_rownames = TRUE,
                  has_colnames = TRUE,
                  char_per_allele = 1,
                  store_G = FALSE) {				
  
  
  # general configuration
  type <- tolower(type)
  type <- match.arg(type)
  missing <- toupper(as.character(missing))
  missing <- match.arg(missing)
  
  # creating a new ecogendb object
  object <- dbR6$new()
  output <- new("genetic_db")
  
  object$write_dataframe(G, "G", chunksize = chunksize, 
                         has_colnames = has_colnames, 
                         has_rownames = has_rownames)    
  G_colnames <- colnames(read.table(G, nrows = 1))
  
  ## coherence between data ploidy and ncod is checked for int.df2genind
  if(type == "codominant") {
    
    if(sep == "") {
      fun_alleles <- "loci_to_alleles"
    } else {
      fun_alleles <- "loci_token_to_alleles"
    }
    
    
    temp_formatG <- formatG(input = matrix(0, 0, 0), ploidy = ploidy, 
                            has_rownames = TRUE,  has_colnames = TRUE,
                            default_levels = character(0),
                            obj_ncol = n_loci, char_per_allele = ncod,
                            colnames_vector = G_colnames)
    
    if(length(default_levels) == 0) {
      
      formatter_chunk_from_db(object, table = "G",
                              formatG_input = temp_formatG,
                              format_fun = fun_alleles,
                              out_df = "temp", 
                              has_rownames = TRUE,  
                              has_colnames = TRUE,
                              ploidy = ploidy,
                              token = sep,
                              input_or_output = "data_input",
                              set_ncol = 0,
                              chunksize = chunksize)
      
      if(!store_G) {
        object$remove_table("G") 
      }
      
      formatter_chunk_from_db(object, table = "temp", 
                              formatG_input = temp_formatG,
                              format_fun = c("alleles_to_dummy"),
                              out_df = "A",
                              has_rownames = TRUE, 
                              has_colnames = TRUE,
                              ploidy = ploidy,
                              input_or_output = "dataS",
                              set_ncol = 0,
                              chunksize = chunksize)
      object$remove_table("temp")
      
    } else {
      formatter_chunk_from_db(object, table = "G",
                              formatG_input = temp_formatG,
                              format_fun = c(fun_alleles, "alleles_to_dummy"),
                              out_df = "A", 
                              has_rownames = TRUE,  
                              has_colnames = TRUE,
                              ploidy = ploidy, 
                              sep = sep,
                              default_levels = default_levels,
                              colnames_vector = G_colnames,
                              input_or_output = "data_input",
                              set_ncol = 0,
                              chunksize = chunksize)
      if(!store_G) {
        object$remove_table("G") 
      }
    }
    
  } else {
    tabnames <- colnames(object$send_query("SELECT * FROM G LIMIT 0"))
    if(tabnames[1] == "row_names") {
      cnames <- tabnames[-1]
      has_table_names = TRUE
    } else {
      has_table_names = FALSE
    }
    columns_binary <-  paste("AVG(", cnames, ")", sep = "", collapse = ",")
    if(missing == "NA") {
      where_clause <- paste(cnames, "!= 'NA'", collapse = " OR ")
      where_clause <- paste("WHERE", where_clause)
    } else {
      where_clause  <- ""
    }
    # compute polymorfism
    avg_query <- paste("SELECT ", columns_binary , " FROM G ", where_clause, collapse = "")
    avg_query_result <- object$send_query(avg_query)
    
    
    if(poly.level < 0 || poly.level > 100) {
      stop("poly.level must be a number between 0 and 100")
    }
    poly.level <- poly.level / 100
    which_retain <- (avg_query_result  >= poly.level) & (100 - avg_query_result  >= poly.level) 
    which_retain <- cnames[which_retain]
    if(has_table_names) {
      which_retain <- c(tabnames[1], which_retain)
    }
    poly_create <- paste("CREATE TABLE A AS SELECT", paste(which_retain, collapse = ", "), " FROM G")
    object$send_statement(poly_create)
    if(!store_G) {
      object$remove_table("G") 
    }
  }
 
  ## construct object
  output@db <- object
  parameters <- get_formatG_parameters(temp_formatG)
  
  output@ATTR$na_omit <- parameters$na_omit
  output@ATTR$which_NA <- parameters$which_NA
  output@ATTR$ploidy <- parameters$ploidy
  output@ATTR$char_per_allele <- parameters$char_per_allele
  output@ATTR$n_loci <- parameters$n_loci
  output@ATTR$map <- parameters$map
  output@ATTR$alleles_number <- parameters$alleles_number
  output@ATTR$alleles_vector <- parameters$alleles_vector
  output@ATTR$total_alleles <- parameters$total_alleles
  output@ATTR$user_levels <- parameters$user_levels
  output@ATTR$sep <- sep
  output@ATTR$missing <- missing
  output@ATTR$type <- type
  output@ATTR$NA.char <- NA.char
  #output@ATTR$rm.empty.ind <- rm.empty.ind
  output@ATTR$poly.level <- poly.level
  
  object
})

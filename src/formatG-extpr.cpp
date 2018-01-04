#include <Rcpp.h>
#include "formatG.h"
using namespace _formatG;
using namespace Rcpp;

//' @param x A single integer.
//' @export
// [[Rcpp::export]]
RcppExport SEXP formatG__new(SEXP obj_, SEXP ploidy_, SEXP default_levels_, 
                             SEXP obj_ncol_, SEXP char_per_allele_,
                             SEXP has_rownames_, SEXP has_colnames_,
                             SEXP colnames_vector_) {
  
  Rcpp::StringMatrix obj(obj_);
  int ploidy = Rcpp::as<int>(ploidy_);
  StringVector default_levels(default_levels_);
  int obj_ncol = Rcpp::as<int>(obj_ncol_); 
  
  if(obj.size() != 0) 
  {
    bool has_rownames = Rcpp::as<bool>(has_rownames_);
    bool has_colnames = Rcpp::as<bool>(has_colnames_);
    // dispatch constructor
    Rcpp::XPtr<formatG> ptr(new formatG(obj, ploidy, default_levels, has_rownames, has_colnames), true);
    return ptr;
  
  } else {
    int char_per_allele = Rcpp::as<int>(char_per_allele_); 
    CharacterVector colnames_vector = Rcpp::as<CharacterVector>(colnames_vector_);
    Rcpp::XPtr<formatG> ptr(new formatG(obj, ploidy, default_levels, obj_ncol, char_per_allele, colnames_vector), true);    
    return ptr;
  }
}


//' @param x A single integer.
//' @export
// [[Rcpp::export]]
RcppExport bool formatG__process_chunk(SEXP ptr, SEXP obj_, SEXP obj_I_,
                                       SEXP fun_name_, SEXP has_rownames_,
                                       SEXP token_, SEXP NA_action_,
                                       SEXP howmuch_, SEXP where_, SEXP what_,
                                       SEXP input_or_output_,
                                       SEXP set_ncol_) {
  
  Rcpp::StringMatrix obj(obj_);
  Rcpp::IntegerMatrix obj_I(obj_I_);
  IntegerVector fun_name(fun_name_);
  bool has_rownames =  Rcpp::as<bool>(has_rownames_);
  char token = Rcpp::as<char>(token_);
  std::string NA_action = Rcpp::as<std::string>(NA_action_);
  int howmuch = Rcpp::as<int>(howmuch_);
  std::string where = Rcpp::as<std::string>(where_);
  char what = Rcpp::as<char>(what_);
  std::string input_or_output = Rcpp::as<std::string>(input_or_output_);
  Rcpp::XPtr<formatG> data(ptr);
  size_t set_ncol = Rcpp::as<size_t>(set_ncol_);
  data->process_chunk(obj, obj_I, fun_name, has_rownames, token, NA_action,  howmuch,
                      where, what, input_or_output, set_ncol);
  return true;
}


// [[Rcpp::export]]
RcppExport bool formatG__delete_formatG(SEXP ext) {
	if (R_ExternalPtrAddr(ext) == NULL)
		return false;
	Rprintf("finalizing\n");
	formatG *ptr = (formatG *) R_ExternalPtrAddr(ext);
	ptr->clear();
	Free(ptr);
	R_ClearExternalPtr(ext);
	return true;
}

//' export attributes
// [[Rcpp::export]]
RcppExport SEXP formatG__get_formatG_parameters(SEXP ptr) {
  Rcpp::XPtr<formatG> data(ptr);
  return data->get_formatG_parameters();
}



//' @param x A single integer.
//' @export
// [[Rcpp::export]]
RcppExport bool formatG__reset_state(SEXP ptr, SEXP reset_lines_, SEXP set_ncol_) {
    bool reset_lines = Rcpp::as<bool>(reset_lines_);
    size_t set_ncol = Rcpp::as<size_t>(set_ncol_);
    Rcpp::XPtr<formatG> data(ptr);
    data->reset_state(reset_lines, set_ncol);
    return true;
  }

//' @param x A single integer.
//' @export
// [[Rcpp::export]]
RcppExport bool formatG__loci_to_alleles(SEXP ptr) {
	Rcpp::XPtr<formatG> data(ptr);
	data->loci_to_alleles_matrix();
	return true;
}

//' @param x A single integer.
//' @export
// [[Rcpp::export]]
RcppExport bool formatG__to_numeric(SEXP ptr) {
	Rcpp::XPtr<formatG> data(ptr);
	data->to_numeric_matrix();
	return true;
}

//' @param x A single integer.
//' @export
// [[Rcpp::export]]
RcppExport bool formatG__loci_token_to_alleles(SEXP ptr, SEXP token_) {
  char token = Rcpp::as<char>(token_);
	Rcpp::XPtr<formatG> data(ptr);
	data->loci_token_to_alleles_matrix(token);
	return true;
}

//' @param x A single integer.
//' @export
// [[Rcpp::export]]
RcppExport bool formatG__alleles_to_dummy(SEXP ptr, SEXP NA_policy_) {
  Rcpp::XPtr<formatG> data(ptr);
  std::string NA_policy = as<std::string>(NA_policy_);
  data->alleles_to_dummy_matrix(NA_policy);
  return true;
}

//' @param x A single integer.
//' @export
// [[Rcpp::export]]
RcppExport Rcpp::StringVector formatG__get_input(SEXP ptr) {
	Rcpp::XPtr<formatG> data(ptr);
	return data->get_input();
}

//' @param x A single integer.
//' @export
// [[Rcpp::export]]
RcppExport SEXP formatG__get_output(SEXP ptr) {
	Rcpp::XPtr<formatG> data(ptr);
	return data->get_output();
}


//' @param x A single integer.
//' @export
// [[Rcpp::export]]
RcppExport SEXP formatG__get_levels(SEXP ptr) {
	Rcpp::XPtr<formatG> data(ptr);
	return data->get_levels();
}


   //' @param x A single integer.
   //' @export
   // [[Rcpp::export]]
   RcppExport bool formatG__add_chars(SEXP ptr, int howmuch, std::string where,
                                      char what) {
     Rcpp::XPtr<formatG> data(ptr);
     data->add_chars(howmuch, where, what);
     return true;
   }

//' @param x A single integer.
//' @export
// [[Rcpp::export]]
RcppExport bool formatG__alleles_to_loci(SEXP ptr, char token) {
  Rcpp::XPtr<formatG> data(ptr);
  data->alleles_to_loci_matrix(token);
  return true;
}

//' @param x A single integer.
//' @export
// [[Rcpp::export]]
RcppExport SEXP formatG__get_alleles_number(SEXP ptr) {
	Rcpp::XPtr<formatG> data(ptr);
	return data->get_alleles_number();
}

//' @param x A single integer.
//' @export
// [[Rcpp::export]]
RcppExport SEXP formatG__get_alleles_vector(SEXP ptr) {
	Rcpp::XPtr<formatG> data(ptr);
	return data->get_alleles_vector();
}

//' @param x A single integer.
//' @export
// [[Rcpp::export]]
RcppExport int formatG__get_total_alleles(SEXP ptr) {
	Rcpp::XPtr<formatG> data(ptr);
	return data->get_total_alleles();
}

//' @param x A single integer.
//' @export
// [[Rcpp::export]]
RcppExport bool formatG__StringM_2_IntM(SEXP ptr) {
	Rcpp::XPtr<formatG> data(ptr);
	data->StringM_2_IntM();
	return true;
}

//' @param x A single integer.
//' @export
// [[Rcpp::export]]
RcppExport bool formatG__dummy_to_alleles_matrix(SEXP ptr) {
	Rcpp::XPtr<formatG> data(ptr);
	data->dummy_to_alleles_matrix();
	return true;
}


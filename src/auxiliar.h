#include <Rcpp.h>
#include "formatG.h"
using namespace Rcpp;

//' Empty objects definition

StringVector NULLstringv(0);
const StringVector empty_stringv(0);

IntegerVector NULLintv(0);
const IntegerVector empty_intv = NULLintv;

StringMatrix NULLstringm(0);
const StringMatrix empty_stringm = NULLstringm;

IntegerMatrix NULLintm(0);
const IntegerMatrix empty_intm = NULLintm;

List NULLlist(0);
const List empty_list = NULLlist;

// to use in switch
#define loci_to_alleles_ 0
#define loci_token_to_alleles_ 1
#define to_numeric_ 2
#define alleles_to_dummy_ 3
#define alleles_to_loci_ 4
#define dummy_to_alleles_ 5
#define add_chars_ 6
#define StringM_2_IntM_ 7

// in other params;
// 0: bool has_rownames, 
// 1: bool has_colnames, 
// 2: char token, 
// 3: std::string NA_action,
// 4: std::string sep,
// 5: int howmuch
// 6: std::string where, 
// 7: char what

namespace _formatG {

// using na_proxy under the hood when testing x[i]== NA
//StringVector NULLstringv(0);

//' which_is_NA_TMP
//' @description Template for NAs search in a Matrix
//' @param obj a Matrix object

template<int RTYPE>
IntegerVector which_is_NA_TMP(Matrix<RTYPE> obj) {

	int nx = obj.size();
	std::vector<int> y;
	y.reserve(nx);

	for (int i = 0; i < nx; i++) {
		if (obj[i] == NA)
			y.push_back(i);
	}

	IntegerVector out(y.size());
	out = y;
	return out;
}

//' which_is_NA
//' @description Dispatch the which_is_NA_TMP function
//' @param obj a Matrix object

IntegerVector which_is_NA(SEXP obj) {
	switch (TYPEOF(obj)) {
	case INTSXP:
		return which_is_NA_TMP<INTSXP>(obj);
		break;
	case STRSXP:
		return which_is_NA_TMP<STRSXP>(obj);
		break;
	}
}

//' set_NA_policy
//' @description formatG method that configures NA policy and a vector indicating which elements are NA.
//' @param obj input Matrix object
//' 
IntegerVector formatG::set_NA_policy(SEXP obj) {
	IntegerVector NA_present(which_is_NA(obj));
	if (NA_present.size() != 0) {
		na_omit = true;
	} else {
		na_omit = false;
	}
	return NA_present;
}

//' CountDigits
//' @description Count digits in a number
//' @param number A number

int countDigits(int number) {
	if (number < 10) {
		return 1;
	}
	int count = 0;
	while (number > 0) {
		number /= 10;
		count++;
	}
	return count;
}

} // end namespace

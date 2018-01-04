#ifndef formatG_H_
#define formatG_H_

#include <Rcpp.h>
#include <iostream>
#include <string>
using namespace Rcpp;

namespace _formatG {

//' formatG class

class formatG {

//---
public:
  // constructor storing a matrix
	formatG(StringMatrix obj, int ploidy, StringVector default_levels,
          bool has_rownames, bool has_colnames); 
  // empty constructor for a given matrix shape
  formatG(StringMatrix obj, int ploidy, StringVector default_levels, 
          int obj_ncol, int char_per_loci, CharacterVector colnames_vector);

	void process_chunk(StringMatrix data_chunk,  IntegerMatrix data_chunk_I, 
                    IntegerVector fun_name, bool has_rownames,
                    char token, std::string NA_action, 
                    int howmuch, std::string where, char what, 
                    std::string input_or_output,
                    size_t set_ncol);
	
	List get_formatG_parameters();
	//virtual ~formatG();

	StringMatrix loci_to_alleles(StringVector input);
	StringMatrix loci_token_to_alleles(StringVector input, char token);
	IntegerMatrix to_numeric(StringMatrix input, List this_map);
	IntegerMatrix alleles_to_dummy(StringMatrix input, List this_map,
			int n_alleles, std::string NA_action);
	StringMatrix alleles_to_loci(StringMatrix input,  char token);
	StringMatrix dummy_to_alleles(IntegerMatrix input, List this_map);

	void loci_to_alleles_matrix();
	void loci_token_to_alleles_matrix(char token);
	void to_numeric_matrix();
	void alleles_to_dummy_matrix(std::string NA_action);
	void alleles_to_loci_matrix(char token);
	void add_chars(int howmuch, std::string where, char what);
	void StringM_2_IntM();
	void dummy_to_alleles_matrix();
	void swap_data(StringMatrix x, int index_1_x, int index_1_y, int index_2_x,
			int index_2_y);
	void swap_data(IntegerMatrix x, int index_1_x, int index_1_y, int index_2_x,
			int index_2_y);
	void sort_matrix(SEXP x);

	bool is_unique_size(StringMatrix input);
	int size(StringMatrix what);
	bool is_unique_size(IntegerMatrix input);
	int size(IntegerMatrix what);
	int get_size();

	void set_levels(StringMatrix input);
	void set_levels_matrix();
	void set_alleles_vector();

	List get_levels();
	int get_sexptype();

	IntegerVector get_alleles_number();
	IntegerVector get_alleles_vector();
	int get_total_alleles();
	friend IntegerVector which_is_NA(StringMatrix obj);
	IntegerVector set_NA_policy(SEXP obj);
	void reset_state(bool reset_lines, size_t set_ncol);
	void clear();

	SEXP get_input();
	SEXP get_output();

	//StringVector
	StringVector set_header(std::string matrix_type);
	void set_generic_rownames(std::string what, int start_from);
	void set_generic_colnames(std::string what, int start_from);

	bool na_omit;
	IntegerVector which_NA;

//---
private:
	StringMatrix data_input;
	StringMatrix dataS;
	IntegerMatrix dataI;

	int ploidy;
	int char_per_allele;
	int nrow_data;
	int ncol_data;
	int n_loci;
	List map;
	std::string state;
	int sexptype;
	IntegerVector alleles_number;
	IntegerVector alleles_vector;
	int total_alleles;
	bool user_levels;

	IntegerVector levels_counter;
	int global_counter;

	// rownames and colnames must be character vectors to use 
	// the Rcpp colnames() and rownames() functions
	CharacterVector rownames_data;
	CharacterVector colnames_data;
	CharacterVector init_colnames;

	int lines_processed;

};

} /* namespace _formatG */

#endif /* formatG_H_ */

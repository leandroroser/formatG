#include <Rcpp.h>
#include <iostream>
#include <string>
#include <sstream>
#include "formatG.h"
#include "auxiliar.h"
using namespace Rcpp;

namespace _formatG {

//' formatG
//' @description formatG object constructor
//' @param obj StringMatrix input object  
//' @param ploidy Ploidy of the matrix 
//' @param na_omit NA omit policy (Logical)

// na_omit initialized with set_NA_policy
formatG::formatG(StringMatrix obj, int ploidy, StringVector default_levels,
		bool has_rownames, bool has_colnames) :
		data_input(obj), dataS(empty_stringm), dataI(empty_intm), ploidy(
				ploidy), char_per_allele(size(obj) / ploidy), map(
				[&] {List x(obj.ncol()); return x;}()), state("input"), nrow_data(
				obj.nrow()), ncol_data(obj.ncol()), n_loci(obj.ncol()), sexptype(
				TYPEOF(obj)), which_NA(set_NA_policy(obj)), alleles_number(
				[&] {IntegerVector out(obj.ncol()); return out;}()), alleles_vector(
				empty_intv), total_alleles(0), 
				user_levels([&default_levels] {if(default_levels.size() != 0) {return true;} else {return false;}}()), 
				levels_counter([&obj] {IntegerVector out(obj.ncol()); return out;}()), 
				global_counter(0), lines_processed(0)

{

	if (default_levels.size() != 0) {
		// set map & alleles_number
		int list_size = default_levels.size();
		List default_levels_list(list_size);
		int counter = 0;
		for (int i = 0; i < list_size; ++i) {
			default_levels_list[i] = counter++;
		}
		default_levels_list.names() = default_levels;

		for (int i = 0; i < n_loci; ++i) {
			map[i] = default_levels_list;
			alleles_number[i] = list_size;
		}

		// set total alleles & alleles_vector
		total_alleles = n_loci * list_size;
		alleles_vector = IntegerVector(total_alleles);

		counter = 0;
		for (int j = 0; j < n_loci; ++j) {
			for (int k = 0; k < list_size; ++k) {
				alleles_vector[counter++] = j;
			}
		}

	}

	if (has_colnames) {
	  colnames_data = colnames(obj);
	} else {
		set_generic_colnames("C", 0);
	}

	if (has_rownames) {
	  rownames_data = rownames(obj);
	} else {
		set_generic_rownames("R", 0);
	}
	
	init_colnames = colnames_data;

	Rcout << "New formatG object" << std::endl;

}


formatG::formatG(StringMatrix obj, int ploidy, StringVector default_levels,
                 int obj_ncol, int char_per_allele, CharacterVector colnames_vector) :
    data_input(empty_stringm), dataS(empty_stringm), dataI(empty_intm), 
    ploidy(ploidy), char_per_allele(char_per_allele), 
    map([&] {List x(obj_ncol); return x;}()), 
    state("input"), nrow_data(0), ncol_data(obj_ncol), n_loci(obj_ncol), 
    sexptype(STRSXP), which_NA(empty_intv), 
    alleles_number([&obj_ncol] {IntegerVector out(obj_ncol); return out;}()), 
    alleles_vector(empty_intv), total_alleles(0), 
    user_levels([&default_levels] {if(default_levels.size() != 0) {return true;} else {return false;}}()), 
    levels_counter([&obj_ncol] {IntegerVector out(obj_ncol); return out;}()), 
    global_counter(0), lines_processed(0), 
    rownames_data([&] {CharacterVector x(0); return x;}()),
    colnames_data([&obj_ncol] {CharacterVector x(obj_ncol); return x;}())

{
  
  if (default_levels.size() != 0) {
    // set map & alleles_number
    int list_size = default_levels.size();
    List default_levels_list(list_size);
    int counter = 0;
    for (int i = 0; i < list_size; ++i) {
      default_levels_list[i] = counter++;
    }
    default_levels_list.names() = default_levels;
    
    for (int i = 0; i < n_loci; ++i) {
      map[i] = default_levels_list;
      alleles_number[i] = list_size;
    }
    
    // set total alleles & alleles_vector
    total_alleles = n_loci * list_size;
    alleles_vector = IntegerVector(total_alleles);
    
    counter = 0;
    for (int j = 0; j < n_loci; ++j) {
      for (int k = 0; k < list_size; ++k) {
        alleles_vector[counter++] = j;
      }
    }
    
  }
  
  if (colnames_vector.size() != 0) {
    colnames_data = colnames_vector;
  } else {
    set_generic_colnames("C", 0);
  }
  init_colnames = colnames_vector;
  
  Rcout << "New empty formatG object" << std::endl;
  
}



//' process_chunk
//' @description process a chunk using one of the available functions
//' @param data_chunk StringMatrix input object  
//' @param fun_name Name of the function to call

// DEFINED IN auxiliar.cpp:
//#define loci_to_alleles_ 0
//#define loci_token_to_alleles_ 1
//#define to_numeric_ 2
//#define alleles_to_dummy_ 3
//#define alleles_to_loci_ 4
//#define dummy_to_alleles_ 5
//#define add_chars_ 6
//#define StringM_2_IntM_ 7


void formatG::process_chunk(StringMatrix data_chunk, IntegerMatrix data_chunk_I,
                            IntegerVector fun_name, bool has_rownames, 
                            char token, std::string NA_action, int howmuch,
                            std::string where, char what, std::string input_or_output,
                            size_t set_ncol) {

	// reset some slots of the object
  reset_state(false, set_ncol);
  
  if(input_or_output == "data_input") 
  {
  data_input = data_chunk;
    
  } else if(input_or_output == "dataS")
  {
    dataS = data_chunk;
    
  } else if(input_or_output == "dataI")
  {
    dataI = data_chunk_I;
  }
  
  which_NA = set_NA_policy(data_chunk);
  nrow_data = data_chunk.nrow();
  ncol_data = data_chunk.ncol();
  sexptype = TYPEOF(data_chunk);
	
	if (has_rownames) {
	  rownames_data = rownames(data_chunk);
	} else {
	  set_generic_rownames("R", lines_processed + 1);
	}
	
	lines_processed += nrow_data;

	for (int i = 0; i < fun_name.size(); ++i) {
	 
		// select function
		try {
			switch (fun_name[i]) {
			case loci_to_alleles_:
				loci_to_alleles_matrix();
				break;
			case loci_token_to_alleles_:
				loci_token_to_alleles_matrix(token);
				break;
			case to_numeric_:
				to_numeric_matrix();
				break;
			case alleles_to_dummy_:
				alleles_to_dummy_matrix(NA_action);
				break;
			case alleles_to_loci_:
				alleles_to_loci_matrix(token);
				break;
			case dummy_to_alleles_:
				dummy_to_alleles_matrix();
				break;
			case add_chars_:
				add_chars(howmuch, where, what);
				break;
			case StringM_2_IntM_:
				StringM_2_IntM();
				break;
			default:
				throw("invalid function selected");
				break;
			}
		} catch (std::string& e) {
			Rcout << e;
		}
	}
}


//' export_attributes

List formatG::get_formatG_parameters()
{
 List output;

output["na_omit"]= na_omit;
output["which_NA"] = which_NA;
output["ploidy"] =  ploidy;
output["char_per_allele"] = char_per_allele;
output["n_loci"] = n_loci;
output["map"] = map;
output["alleles_number"] = alleles_number;
output["alleles_vector"] = alleles_vector;
output["total_alleles"] = total_alleles;
output["user_levels"] = user_levels;
return output;
}

//' alleles_to_loci
//' @description Creates a individuals x 1 matrix,  with a loci in a column, 
//' from a two-column matrix with the alleles 
//' @param input StringMatrix input object  
//' @param sep character dividing alleles 

StringMatrix formatG::alleles_to_loci(StringMatrix input,  char token) {
	StringMatrix output(nrow_data, 1);
	std::ostringstream oss;

	for (int i = 0; i < input.nrow(); ++i) {
		for (int j = 0; j < input.ncol(); ++j) {
			if (j == 0) {
				oss << input(i, j);
			} else {
				oss << token << input(i, j);
			}
		}
		output[i] = oss.str();
		oss.str("");
	}
	return output;
	
}

//' alleles_to_loci_matrix
//' @description Wrapper that calls alleles_to_loci for a matrix with multiple loci
//' @param sep character dividing alleles 

void formatG::alleles_to_loci_matrix(char token) {

	StringMatrix output(nrow_data, ncol_data / ploidy);
	StringMatrix temp(nrow_data, 1);
	StringMatrix temp0(nrow_data);

	int counter = 0;
	for (int i = 0; i < ncol_data; i += ploidy) {
		temp0 = dataS(_, Range(i, i + ploidy - 1));

		temp = alleles_to_loci(temp0, token);
		
		output(_, counter++) = temp;
	}

	dataI = empty_intm;

	dataS = output;
	dataS.attr("dimnames") = List::create(rownames_data, colnames_data);

	state = "bind_string";
	ncol_data = ncol_data / ploidy;
	sexptype = TYPEOF(dataS);

}

//' add_chars
//' @description Add characters at the beggining/end of each cell of the matrix 
//' referenced in a formatG object
//' @param howmuch Number of characters
//' @param where "start" to add characters at the beginning of the string, 
//' "end" to add characters at the end of the string
//' @param what Character to add

void formatG::add_chars(int howmuch, std::string where, char what) {
	std::ostringstream oss;
	std::ostringstream to_paste;
	std::string zeros;
	StringMatrix output(nrow_data, ncol_data);

	try {
		if (where != "start" || where != "end") {
			throw("where must be 'start' or 'end' ");
		}
	} catch (std::string s) {
		Rcout << s;
	}

	for (int i = 0; i < howmuch; ++i) {
		to_paste << what;
	}
	zeros = to_paste.str();

	if (where == "end") {
		for (int i = 0; i < dataS.size(); ++i) {
			output[i] = as < std::string > (dataS[i]) + zeros;
		}

	} else if (where == "start") {
		for (int i = 0; i < dataS.size(); ++i) {
			output[i] = zeros + as < std::string > (dataS[i]);
		}
	}

	// dataI = IntegerVector::create();
	dataS = output;
	dataS.attr("dimnames") = List::create(rownames_data, colnames_data);

}

//' loci_to_alleles
//' @description Split locus in column into alleles in columns
//' @param input individuals x 1 StringMatrix, with locus in column

StringMatrix formatG::loci_to_alleles(StringVector input) {

	std::string buffer("");
	StringMatrix output(nrow_data, ploidy);
	std::string temp("");
	int counter = 0;

	for (int j = 0; j < nrow_data; ++j) {

		int counter = 0;

		if (na_omit) {
			if (input[j] == NA)
				continue;
		}

		temp = as < std::string > (input[j]);

		// cout << "texto completo: "; temp << endl;
		for (int i = 0; i < char_per_allele * ploidy; ++i) {
			buffer = buffer + temp[i];
			if (((i + 1) % char_per_allele) == 0) {
				output[counter * nrow_data + j] = buffer;
				buffer = "";
				counter++;
			}
		}

	}

	return output;

}

//'loci_to_alleles_matrix
//' @description Wrapper that calls loci_to_alleles, 
//' for a matrix with multiple loci, stored in a formatG object

void formatG::loci_to_alleles_matrix() {

  int ncol_output = ploidy * ncol_data;
	StringVector temp(nrow_data);
	StringMatrix splitted(nrow_data, ploidy);
	StringMatrix output(nrow_data, ncol_output);

	int counter = 0;
	for (int i = 0; i < ncol_data; ++i) {
		int index = 0;
		temp = data_input(_, i);
		splitted = loci_to_alleles(temp);
		while (index < ploidy) {
			output(_, counter++) = splitted(_, index++);
		}

		index = 0;
	}

	// NA policy
	if (na_omit) {

		for (int s = 0; s < which_NA.size(); ++s) {
			counter = 0;
			while (counter < ploidy) {
				output[counter * nrow_data + which_NA[s]] = NA_STRING;
				counter++;
			}
		}
	} // end NA policy

	dataS = output;
	dataS.attr("dimnames") = List::create(rownames_data, set_header("alleles"));
	ncol_data = ncol_output;

	dataI = empty_intm;
	state = "split_string";
	sexptype = TYPEOF(output);

	if (!user_levels) {
		set_levels_matrix();
		set_alleles_vector();
	}

}

//' loci_token_to_alleles
//' @description Split locus in column with alleles separated by a character into alleles in columns
//' @param input individuals x 1 StringMatrix, with locus in column
//' @param token  Character separating alleles

StringMatrix formatG::loci_token_to_alleles(StringVector input, char token) {

	StringMatrix output(nrow_data, ploidy);
	std::string buffer("");

	for (int i = 0; i < nrow_data; ++i) {

		if (na_omit) {
			if (input[i] == NA_STRING)
				continue;
		}

		std::istringstream iss(as<std::string> (input[i]));

		int counter = 0;
		while (std::getline(iss, (buffer), token)) {
			output[counter * nrow_data + i] = buffer;
			counter++;
		}
		iss.str();
		buffer = "";
	}

	return output;
}

//'loci_token_to_alleles_matrix
//' @description Wrapper that calls loci_token_to_alleles, 
//' for a matrix with multiple loci, stored in a formatG object
//' @param token Character separating alleles

void formatG::loci_token_to_alleles_matrix(char token) {

	int ncol_output = ploidy * ncol_data;
	StringVector temp(nrow_data);
	StringMatrix splitted(nrow_data, ploidy);
	StringMatrix output(nrow_data, ncol_output);

	int counter = 0;
	for (int i = 0; i < ncol_data; ++i) {

		int index = 0;
		temp = data_input(_, i);
		splitted = loci_token_to_alleles(temp, token);
		while (index < ploidy) {
			output(_, counter++) = splitted(_, index++);
		}
		index = 0;
	}

	dataS = output;
	dataS.attr("dimnames") = List::create(rownames_data, set_header("alleles"));
	ncol_data = ncol_output;

	dataI = empty_intm;
	state = "split_string";
	sexptype = TYPEOF(output);

	if (!user_levels) {
		set_levels_matrix();
		set_alleles_vector();
	}

}

//' to_numeric
//' @description Recodes a genetic matrix into numernc format
//' @param input individuals x 1 StringMatrix, with locus in column
//' @param this_map key-values list, generated with the formatG object (internal)

IntegerMatrix formatG::to_numeric(StringMatrix input, List this_map) {
	std::string this_string;

	IntegerMatrix output(nrow_data, ploidy);

	for (int i = 0; i < input.size(); ++i) {
		this_string = as < std::string > (input[i]);

		output[i] = this_map[this_string];
	}

	return output;

}

//' to_numeric_matrix
//' @description Wrapper that calls to_numeric, 
//' for a matrix with multiple loci, stored in a formatG object

void formatG::to_numeric_matrix() {

	StringMatrix temp(nrow_data, ploidy);
	IntegerMatrix results(nrow_data, ploidy);
	IntegerMatrix output(nrow_data, ncol_data);
	List this_map;

	int counter_numeric = 0;
	int counter_map = 0;

	for (int i = 0; i < ncol_data; i += ploidy) {

		temp = dataS(_, Range(i, i + ploidy - 1));
		this_map = map[counter_map++]; // dependo de map
		results = to_numeric(temp, this_map);

		for (int k = 0; k < ploidy; ++k) {
			output(_, counter_numeric++) = results(_, k);
		}
	}

	sexptype = TYPEOF(output);
	state = "split_int";
	dataI = output;
	dataI.attr("dimnames") = dataS.attr("dimnames");

	dataS = empty_stringm;
}

//' StringM_2_IntM
//' @description Converts StringMatrix into a IntegerMatrix 

void formatG::StringM_2_IntM() {

	IntegerMatrix out(nrow_data, ncol_data);
	if (dataS.size() != 0) {
		for (int i = 0; i < dataS.size(); ++i) {
			out[i] = std::stoi(as < std::string > ((dataS[i])));
		}
	} else {

		for (int i = 0; i < data_input.size(); ++i) {
			out[i] = std::stoi(as < std::string > (data_input[i]));
		}
	}

	sexptype == INTSXP;
	dataI = out;

	dataI.attr("dimnames") = dataS.attr("dimnames");

	dataS = empty_stringm;
}

//' alleles_to_dummy
//' @description Converts a matrix of alleles into a matrix of dummy factors
//' @param input individuals x 2 StringMatrix, with alleles in columns
//' @param NA_action The NA action ("0" or "NA") for recoding NAs.
//' @param this_map key-values list, generated with the formatG object (internal)
//' @param n_alleles Number of alleles present in the loci (internal)

IntegerMatrix formatG::alleles_to_dummy(StringMatrix input, List this_map,
		int n_alleles, std::string NA_action) {
	int this_value;
	// see that () initializes array to zero!
	IntegerMatrix output(nrow_data, n_alleles);

	for (int i = 0; i < nrow_data; ++i) {
		for (int j = 0; j < ploidy; ++j) {
			String temp = input(i, j);

			if (na_omit && (temp == NA_STRING)) {
				if (NA_action == "0") {
					continue;

				} else if (NA_action == "NA") {
					// set al alles as NA
					for (int k = 0; k < n_alleles; ++k) {
						output(i, k) = NA_INTEGER;
					}
					continue;
				}
			} // end NA policy

			this_value = this_map[temp];

			output(i, this_value)++;}

		}

	return output;
}

//' alleles_to_dummy_matrix
//' @description Wrapper that calls alleles_to_dummy,
//' for a matrix with multiple loci (with alleles in columns), stored in a formatG object

void formatG::alleles_to_dummy_matrix(std::string NA_action) {
	int counter_dummy = 0; // reset persistent counters
	int alleles_dummy = 0;
	
	IntegerMatrix output(nrow_data, total_alleles);
	// total_levels se genera en set_alleles_vector
	int alleles_in_loci = 0;

	for (int i = 0; i < ncol_data; i += ploidy) {

		// List this_map; this_map = map[allleles_dummy]; did not work
		List this_map(map[alleles_dummy]);

		alleles_in_loci = alleles_number[alleles_dummy++];
		IntegerMatrix temp0(nrow_data, alleles_in_loci);

		StringMatrix temp(nrow_data, ploidy);

		temp = dataS(_, Range(i, i + ploidy - 1));

		temp0 = alleles_to_dummy(temp, this_map, alleles_in_loci, NA_action);

		for (int k = 0; k < alleles_in_loci; ++k) {
			output(_, counter_dummy++) = temp0(_, k);
		}

	}

  dataI = output;
	dataI.attr("dimnames") = List::create(rownames_data, set_header("dummy"));
	ncol_data = total_alleles;

	dataS = empty_stringm;
	state = "dummy";
	sexptype = TYPEOF(dataI);
}

//' dummy_to_alleles
//' @description Inverse operation of alleles_to_dummy
//' @param input individuals x 2 StringMatrix, with alleles in columns
//' @param this_map key-values list, generated with the formatG object (internal)

StringMatrix formatG::dummy_to_alleles(IntegerMatrix input, List this_map) {

	StringVector this_names(this_map.names());
	StringMatrix output(nrow_data, ploidy);
	int this_value;

	for (int i = 0; i < nrow_data; ++i) {
		int counter = 0;
		for (int j = 0; j < input.ncol(); ++j) {
			this_value = input(i, j);
			int map_value = this_map[j];

			int k = 0;
			while (k < this_value) {
				output(i, counter++) = this_names[map_value];
				k++;
			}
		}
	}

	return output;
}

//' dummy_to_alleles_matrix
//' @description Wrapper that calls dummy_to_alleles 
//' for a matrix with dummy alleles of multiple loci 

void formatG::dummy_to_alleles_matrix() {

	StringMatrix output(nrow_data, n_loci * ploidy);
	// total_levels se genera en set_alleles_vector

	int counter_column = 0;
	int counter_loci = 0;
	int counter_alleles = 0;
	while (counter_alleles < total_alleles) {

		int alleles_in_loci = alleles_number[counter_loci];
		IntegerMatrix temp(nrow_data, alleles_in_loci);
		StringMatrix temp0(nrow_data, ploidy);

		temp = dataI(_,Range(counter_alleles, counter_alleles + alleles_in_loci - 1));
		
		List this_map = map[counter_loci];
    
		temp0 = dummy_to_alleles(temp, this_map);
    
		for (int k = 0; k < ploidy; ++k) {
			output(_, counter_column) = temp0(_, k);
			counter_column++;
		}
		counter_alleles += alleles_in_loci;
		counter_loci++;
	}

	if (na_omit) {

		for (int s = 0; s < which_NA.size(); ++s) {
			int counter = 0;
			while (counter < ploidy) {
				output[counter * nrow_data + which_NA[s]] = NA_STRING;
				counter++;
			}
		}
	} // end NA policy

	dataS = output;
	
	// first set n_col to create dimnames
	ncol_data = n_loci * ploidy;
  dataS.attr("dimnames") = List::create(rownames_data, set_header("alleles"));
  
	dataI = empty_intm;

	state = "column";
	sexptype = TYPEOF(dataS);
}

//----------

//' is_unique_size
//' @description Verify if all cells of a matrix have the same number of digits / characters
//' @param input input StringMatrix 
//' 
bool formatG::is_unique_size(StringMatrix input) {
	StringVector temp0(1);
	temp0 = input[0];
	StringVector temp(1);

	// check that first element is not NA, else pick the next
	if (temp0[0] == NA_STRING) {

		try {
			if (na_omit) {
				int counter = 0;
				while (counter < input.size()) {
					temp0[0] = input[++counter];
					if (temp0[0] != NA_STRING)
						break;
				}
			} else {
				throw 99;
			}

		} catch (int& e) {
			Rcout << "Error: NAs present, but na.omit is false" << std::endl;
		}
	}

	// still NA?
	if (temp0[0] == NA_STRING) {
		return NA_LOGICAL;
	}

	for (int i = 0; i < input.size(); ++i) {

		temp = input[i];
		if (temp[0] == NA_STRING) {

			try {
				if (na_omit) {
					continue;
				} else {
					throw 99;
				}
			} catch (int e) {
				Rcout << "Error: NAs present, but na.omit is false"
						<< std::endl;
			}
		}

		if (temp0[0].size() != temp[0].size()) {
			return false;
		} else {
			temp0 = temp;
		}
	}
	return true;
}

//' is_unique_size
//' @description Verify if all cells of a matrix have the same number of digits / characters
//' @param input input IntegerMatrix

bool formatG::is_unique_size(IntegerMatrix input) {
	StringVector temp0(1);
	temp0 = input[0];
	StringVector temp(1);

	// check that first element is not NA, else pick the next
	if (temp0[0] == NA_INTEGER) {

		try {
			if (na_omit) {
				int counter = 0;
				while (counter < input.size()) {
					temp0[0] = input[++counter];
					if (temp0[0] != NA_INTEGER)
						break;
				}
			} else {
				throw 99;
			}

		} catch (int& e) {
			Rcout << "Error: NAs present, but na.omit is false" << std::endl;
		}
	}

	// still NA?
	if (temp0[0] == NA_INTEGER) {
		return NA_LOGICAL;
	}

	for (int i = 0; i < input.size(); ++i) {

		temp = input[i];
		if (temp[0] == NA_INTEGER) {

			try {
				if (na_omit) {
					continue;
				} else {
					throw 99;
				}
			} catch (int e) {
				Rcout << "Error: NAs present, but na.omit is false"
						<< std::endl;
			}
		}

		if (temp0[0].size() != temp[0].size()) {
			return false;
		} else {
			temp0 = temp;
		}
	}
	return true;
}

//' size
//' @description Returns the number of characters / digits of a matrix
//' @param input input StingMatrix

int formatG::size(StringMatrix input) {

	// ver como ponerle una accion para na_omit
	// cuando el primer elemento es cero
	int charsize = 0;
	bool this_is_unique;

	try {
		this_is_unique = is_unique_size(input);

		// all na policy
		if (this_is_unique == NA_LOGICAL) {
			return NA_LOGICAL;
		}

		if (this_is_unique) {

			// na omit policy
			if (na_omit) {
				int k = 0;
				while (k < input.size()) {
					if (input[k] == NA_STRING) {
						k++;
						continue;
					} else {
						charsize = input[k].size();
						break;
					}
				}
				// end na policy
			} else {
				charsize = input[0].size();
			}

		} else {
			throw 1;
		}
	} catch (int& num) {
		Rcout << "non unique size in array" << std::endl;
	} // end try-catch

	return charsize;
}

//' size
//' @description Returns the number of characters / digits of a matrix
//' @param input input IntegerMatrix

int formatG::size(IntegerMatrix input) {

	// ver como ponerle una accion para na_omit
	// cuando el primer elemento es cero
	int charsize = 0;
	bool this_is_unique;

	try {
		this_is_unique = is_unique_size(input);

		// all na policy
		if (this_is_unique == NA_LOGICAL) {
			return NA_LOGICAL;
		}

		if (this_is_unique) {

			// na omit policy
			if (na_omit) {
				int k = 0;
				while (k < input.size()) {
					if (input[k] == NA_INTEGER) {
						k++;
						continue;
					} else {
						charsize = countDigits(input[k]);
						break;
					}
				}
				// end na policy
			} else {
				charsize = countDigits(input[0]);
			}

		} else {
			throw 1;
		}
	} catch (int& num) {
		Rcout << "non unique size in array" << std::endl;
	} // end try-catch

	return charsize;
}

//--------------
CharacterVector formatG::set_header(std::string matrix_type) {
  int header_elements;
  std::ostringstream os;
  if (matrix_type == "alleles") {
    CharacterVector  output(n_loci * ploidy);
    int i = 0;
    int j = 0;
    for (int i = 0; i < n_loci; ++i) {
      int k = 1;
      while (k <= ploidy) {
        os << colnames_data[i] << "_" << k;
        output[j++] = os.str();
        os.str("");
        os.clear();
        k++;
      }
    }
    return output;
    
  } else if (matrix_type == "dummy") {
    CharacterVector output(total_alleles);
    int k = 0;
    
    for (int i = 0; i < alleles_number.size(); ++i) {
      for (int j = 1; j <= alleles_number[i]; ++j) {
        os << colnames_data[i] << "_" << j;
        output[k++] = os.str();
        os.str("");
        os.clear();
      }
    }
    return output;
  }
}
//------------------

//' get_size
//' @description Accessor to the number of characters / digits of a matrix stored in a formatG object

int formatG::get_size() {
	return char_per_allele * ploidy;
}

//' set_levels
//' @description Set the levels of a loci (matrix with alleles in columns, individuals x ploidy)
//' @param input A StringMatrix

void formatG::set_levels(StringMatrix input) {

	std::string temp;
	std::string this_element_j("");

	bool add = true;
	List temp_map(map[global_counter]);

	//for(int i = 1;  i < nrow_data * ploidy ; ++i)
	for (int i = 0; i < input.size(); ++i) {
		temp = as < std::string > (input[i]);

		if (temp_map.size() != 0) {
			for (int j = 0; j < temp_map.size(); ++j) {
				if (na_omit) {
					if (input[i] == NA_STRING) {
						add = false;
						break;
					}
				}

				StringVector map_names(temp_map.names());

				this_element_j = map_names[j];

				if (this_element_j == temp) {
					add = false;
					break;
				}
			}

			if (add) {
				temp_map[temp] = levels_counter[global_counter];
				levels_counter[global_counter]++;
			}
			add = true;

		} else {
			temp_map[temp] = levels_counter[global_counter];
			levels_counter[global_counter]++;
		}
	}

	alleles_number[global_counter] = temp_map.size();
	map[global_counter] = temp_map;
}

//' set_levels_matrix
//' @description Wrapper for set_levels, to be used
//' @param input A StringMatrix

void formatG::set_levels_matrix() {

	StringMatrix temp(nrow_data, ploidy);

	for (int i = 0; i < ncol_data; i += ploidy) {
		temp = dataS(_, Range(i, i + ploidy - 1));
		set_levels(temp);
		global_counter++;
	}

	global_counter = 0; // reset becuse persists between calls!

}

//' set_alleles_vector
//' @description Set an internal vector with the number of alleles of each loci

void formatG::set_alleles_vector() {

	int n_alleles = alleles_number.size(); // check this

  total_alleles = 0;
  
	for (int i = 0; i < n_loci; ++i) {
		total_alleles += alleles_number[i];
	}

	IntegerVector output(total_alleles);

	int counter = 0;
	for (int j = 0; j < n_alleles; ++j) {
		for (int k = 0; k < alleles_number[j]; ++k) {
			output[counter++] = j;
		}

	}

	alleles_vector = output;
}

void formatG::set_generic_rownames(std::string what, int start_from) {
	std::ostringstream os;
	std::vector<std::string> output;
	output.reserve(nrow_data);
	for (int i = start_from; i < start_from + nrow_data; ++i) {
		os << what << "_" << i;
		output.push_back(os.str());
		os.str("");
		os.clear();
	}
	rownames_data = output;
}

void formatG::set_generic_colnames(std::string what,  int start_from) {
	std::ostringstream os;
	std::vector<std::string> output;
	output.reserve(ncol_data);
	for (int i = start_from; i < start_from + ncol_data; ++i) {
		os << what << "_" << i;
		output.push_back(os.str());
		os.str("");
		os.clear();
	}
	colnames_data = output;
}

//' get_input
//' @description Accessor to the input matrix stored in a formatG object

SEXP formatG::get_input() {
	return data_input;
}

//' get_output
//' @description Accessor to the output matrix stored in a formatG object

SEXP formatG::get_output() {
	if (dataI.size() == 0) {
		return dataS;

	} else {
		return dataI;
	}
}

//' get_levels
//' @description Accessor to internal List with value-keys pairs of a formatG object

List formatG::get_levels() {
	return map;
}

//' get_sexptype
//' @description Accessor to the SEXP type of the output matrix stored in a formatG object

int formatG::get_sexptype() {
	return sexptype;
}

//' get_alleles_number
//' @description Get the number of alleles stored in a formatG object

IntegerVector formatG::get_alleles_number() {
	return alleles_number;
}

//' get_alleles_vector
//' @description Accessor to internal vector with the number of alleles of each loci
//' 
IntegerVector formatG::get_alleles_vector() {
	return alleles_vector;
}

//' get_total_alleles
//' @description Get the total number of alleles of the output matrix stored in a formatG object 

int formatG::get_total_alleles() {
	return total_alleles;
}

//' reset_state
//' @description Reset state for chunk processing

void formatG::reset_state(bool reset_lines = false,  size_t set_ncol = 0) {
data_input = empty_stringm;
dataS = empty_stringm;
dataI = empty_intm;
state = "input";
nrow_data = 0;
if(set_ncol != 0)
{
ncol_data = set_ncol;
}
rownames_data = empty_stringv;
colnames_data = init_colnames;
if(reset_lines) 
  lines_processed = 0;
}

//' clear
//' @description Clear a formatG object

void formatG::clear() {
	data_input = empty_stringm;
	dataS = empty_stringm;
	dataI = empty_intm;
	map = empty_list;
	state = std::string("");
	alleles_number = empty_intv;
	which_NA = empty_intv;
	levels_counter = empty_intv;

	na_omit = false;
	ploidy = 0;
	char_per_allele = 0;
	nrow_data = 0;
	ncol_data = 0;
	n_loci = 0;
	sexptype = 0;
	total_alleles = 0;
	user_levels = false;
	global_counter = 0;
}

// it is required a validator to perform dispatch
// the following template does not work
/*
 template <typename T0, typename T1, typename T2>
 bool typed_valid( SEXP* args, int nargs ){
 return nargs == 3 && is<T0>(args[0]) && is<T1>(args[1]) && is<T2>(args[2]) ;
 }
 */

//' validate_string
//' @description Validate the input of a formatG object. Useful for dispatching
//' @param SEXP* args objects passed to the function
//' @param nargs Number of arguments
bool validate_C1(SEXP* args, int nargs) {
  if( nargs != 5 ) return false ;
	return TYPEOF(args[3]) == LGLSXP;
}

bool validate_C2(SEXP* args, int nargs) {
  if( nargs != 6 ) return false ;
  return TYPEOF(args[3]) == INTSXP;
}


//bool validate_int(SEXP* args, int nargs)
//{ return TYPEOF(args[0]) == INTSXP;}

//' formatG_finalizer
//' @description finalizer of a formatG object
//' @param formatG* ptr pointer to formatG object

static void formatG_finalizer(formatG* ptr) {
	ptr->clear();
	Rcout << "formatG finalizer called...\n";
}

//' formatG_module
//' @description formatG module

RCPP_MODULE(formatG_module)
{
	class_<formatG>("formatG")

	//.constructor<IntegerMatrix, int, bool>("int constructor", &validate_int)
 .constructor<StringMatrix, int, StringVector, bool, bool>("string constructor", &validate_C1)
 .constructor<StringMatrix, int, StringVector, int, int, CharacterVector>("empty constructor", &validate_C2)
	.finalizer( &formatG_finalizer )
	// .field( "output", &formatter::output)
	//  .field( "output_int", &formatter::output_int)
	// .field( "input", &formatter::input)

	.method("loci_to_alleles_matrix", &formatG::loci_to_alleles)
	.method("to_numeric_matrix", &formatG::to_numeric)
	.method("loci_token_to_alleles_matrix", &formatG::loci_token_to_alleles)
	.method("alleles_to_dummy_matrix", &formatG::alleles_to_dummy)
	.method("get_input", &formatG::get_input)
	.method("get_output", &formatG::get_output)
	.method("add_chars", &formatG::add_chars)
	.method("alleles_to_loci_matrix", &formatG::alleles_to_loci)
	.method("get_alleles_number", &formatG::get_alleles_number)
	.method("get_alleles_vector", &formatG::get_alleles_vector)
	.method("get_total_alleles", &formatG::get_total_alleles)
  .method("reset_state", &formatG::reset_state)
  .method("get_formatG_parameters", &formatG::get_formatG_parameters)
	;
}

} // end namespace

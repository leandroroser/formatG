## formatG

formatG is an R package for genetic data formatting, written in C++. 
The package allows to work with small and big data tables. For this latter purpose, 
the package uses the libraries "dbR6"" and "chunkR" (the latter also developed in C++) and available on this GitHub repository ([https://github.com/leandroroser/chunkR](https://github.com/leandroroser/chunkR), [https://github.com/leandroroser/dbR6](https://github.com/leandroroser/dbR6)).
The dbR6 package consists of an interface to SQLite based on a R6 class. The interaction with the formatG package allows to store, process and retrieve genetic data from a database in chunks, thus removing the need to store data in memory and allowing to process large files. The chunkR package allows to read large files from disk.


### Example

```diff

# let's create a table

a<-replicate(100, paste(sample(c("a", "t", "c", "g"), 2, replace=TRUE), collapse = ""))
a<-matrix(a, 50, 2)
colnames(a)<-c("L1", "L2")
rownames(a) <- 1:50

# create formatG object
my_g <-formatG(a, 2, default_levels = c("a", "t", "c", "g"), TRUE, TRUE)

# loci to alleles
loci_to_alleles(my_g)
get_output(my_g) ## get data

# alleles to dummy
alleles_to_dummy(my_g, NA_action = "NA")
get_output(my_g) ## get data

# dummy to alleles
dummy_to_alleles(my_g)
get_output(my_g) ## get data

# alleles to loci
alleles_to_loci(my_g, sep = "")
get_output(my_g) ## get data

# loci to alleles, for alleles separated by a character
## let's create a table
b<-replicate(100, paste(sample(c("a", "t", "c", "g"), 2, replace=TRUE), collapse = "."))
b<-matrix(b, 50, 2)
colnames(b)<-c("L1", "L2")
rownames(b) <- 1:50

# create formatG object
my_g2 <-formatG(b, 2, default_levels = c("a", "t", "c", "g"), TRUE, TRUE)

loci_token_to_alleles(my_g2, ".")
get_output(my_g2)

# to_numeric
to_numeric(my_g2)
get_output(my_g2)


# Integration with the dbR6 package

library(chunkR)
library(dbR6)

# let's create a big table
my_db <- dbR6$new("my_db", overwrite = TRUE)

a<-replicate(10000, paste(sample(c("a", "t", "c", "g"), 2, replace=TRUE), collapse = ""))
a<-matrix(a, 5000, 2)
write.table(a, "long_table.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

## processing from external file

                
## processing from table stored in db
my_db$write_matrix("long_table.txt", "long_table", has_colnames = FALSE, has_rownames = FALSE)
formatter_chunk_from_db(my_db, table = "long_table",
                format_fun = c("loci_to_alleles", "alleles_to_dummy"),
                out_df = "dummy_db", has_rownames = FALSE,  has_colnames = FALSE,
                ploidy = 2, default_levels = c("a", "t", "c", "g"))

```

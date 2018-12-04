vcf coding challenge

This project parses a vcf file and output a table as a result. The output table contains some fields that are interesting to the biologist.

Installation

This is an R script. You need to install R on your computer. The R script utilizes three third-party package : httr, jsonlite, and magrittr, which should be installed as well.



Usage

This script is intentionally designed for the coding challenge, so the input and output file names are hard-coded into the script.

    input_file <- "input_coding_challenge_final.vcf"
    output_file <- "output_vcf_useful_info.txt"

If the script will used for other same format input files, the user can uncomment the statements related to the command line arguments:

    args<-commandArgs(TRUE)
    if (length(args) < 2) {
        print("usage: Rscript parse_vcf.R <input_file> <output_file>")
        stop("incorrect arguments.")
    }
    input_file <- args[1]
    output_file <- args[2]



Information extracted from the vcf file

1. Variant type (e.g. insertion, deletion, etc.).
2. Variant effect (e.g. missense, synonymous, etc.). Note: If multiple variant types exist in the ExAC database, annotate with the most deleterious possibility. 
3. Read depth at the site of variation.
4. Number of reads supporting the variant.
5. Percentage of reads supporting the variant versus those supporting reference reads.
6. Allele frequency of variant.
7. Any other information from ExAC that you feel might be relevant.



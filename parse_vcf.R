rm(list=ls())
ptm <- proc.time()
options(stringsAsFactors = F)
# args<-commandArgs(TRUE)


# install.packages("httr")
# install.packages("jsonlite")
# install.packages("magrittr")


library("httr")
library("jsonlite")
library("magrittr")


# if (length(args) < 2) {
#     print("usage: Rscript parse_vcf.R <input_file> <output_file>")
#     stop("incorrect arguments.")
# }
# input_file <- args[1]
# output_file <- args[2]

# Please change the input and output accordingly.
input_file <- "input_coding_challenge_final.vcf"
output_file <- "output_vcf_useful_info.txt"


# Read the VCF file and tidy it in a table
raw_text <- readLines(input_file)
context_text <- raw_text[!grepl("^##", raw_text)]

header_text <- context_text[grepl("^#", context_text)]
header_string <- gsub("#", "", header_text)
header_vector <- strsplit(header_string, "\t")[[1]]

formal_text <- context_text[!grepl("^#", context_text)]

vcf_table <- read.table(header = F, text = formal_text)
colnames(vcf_table) <- header_vector



# if multiallelic, decompose into individual mutations.
# ALT and INFO's AF need to be split.
vcf_df <- vcf_table[0, ]
for (i in seq_len(nrow(vcf_table))) {
    multiallelic_condition <- grepl(",", vcf_table[i, "ALT"])

    if (multiallelic_condition == TRUE) {
        
        AF_vector <- vcf_table[i, "INFO"] %>%
            strsplit(., ";") %>%
            sapply(., function(x) x[grepl("^AF", x)]) %>%
            gsub("AF=", "", .) %>%
            strsplit(., ",") %>%
            .[[1]]
        
        multiallelic_vector <- strsplit(vcf_table[i, "ALT"], ",")[[1]]

        for (j in seq_along(multiallelic_vector)) {
            temporary_row <- vcf_table[i, ]
            temporary_row[, "ALT"] <- multiallelic_vector[j]

            AF_replacement <- paste(";AF=", AF_vector[j], ";", sep = "")
 
            temporary_row[, "INFO"] <- gsub(";AF=.*?;", AF_replacement, temporary_row[, "INFO"])
            vcf_df <- rbind.data.frame(vcf_df, temporary_row)
        }
    } else {
        vcf_df <- rbind.data.frame(vcf_df, vcf_table[i, ])
    }
}


# Fields are concatenated to be Variant Id: CHROMOSOME-POSITION-REFERENCE-VARIANT
variant_IDs <- paste(gsub("chr", "", vcf_df$CHROM), vcf_df$POS, vcf_df$REF, vcf_df$ALT, sep = "-")


# From now on, calculate the desired questions:

# 1.Variant type (e.g. insertion, deletion, etc.).
variant_types <- vcf_df$INFO %>%
    strsplit(., ";") %>%
    sapply(., function(x) x[grepl("TYPE", x)]) %>%
    gsub("TYPE=", "", .)


# 2.Variant effect (e.g. missense, synonymous, etc.). Note: If multiple variant types exist in the ExAC database, annotate with the most deleterious possibility.
exac_query_func <- function(a_variant_id) {
    
    base_url <- "http://exac.hms.harvard.edu"
    # rest_url <- "rest/awesome?query=CHD8&service=variants_in_gene"
    # rest_url <- "rest/variant/22-46615880-T-C"
    
    rest_url <- paste("rest/variant", a_variant_id, sep = "/")
    complete_url <- paste(base_url, rest_url, sep = "/")
    resource_json <- GET(complete_url)
    
    resource_text <- content(resource_json, as = "text", encoding = "UTF-8")
    resource_list <- fromJSON(resource_text)
    return(resource_list)
}

percent <- function(x, digits = 2, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

whole_exac_list <- list()
consequence_list <- list()
for (i in seq_along(variant_IDs)) {
    current_resource_list <- exac_query_func(variant_IDs[i])
    
    whole_exac_list[[variant_IDs[i]]] <- current_resource_list
    if (is.null(current_resource_list$consequence)) {
        consequence_list[[variant_IDs[i]]] <- NA
    } else {
        consequence_list[[variant_IDs[i]]] <- names(current_resource_list$consequence)
    }
    
    running_percent <- percent(i/length(variant_IDs))
    cat(running_percent, "has been queried.", sep = " ", "\n")

}

possible_variant_effects <- names(table(unname(unlist(consequence_list))))

deleterious_severity_vector <- c("stop_gained", "stop_lost", "stop_retained_variant", "initiator_codon_variant", "incomplete_terminal_codon_variant", "missense_variant", "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "mature_miRNA_variant", "intron_variant", "non_coding_transcript_exon_variant", "synonymous_variant")

find_the_most_deleterious_variant <- function(to_sort_vector, reference_vector = deleterious_severity_vector) {
    
    # How do I sort one vector based on values of another
    # https://stackoverflow.com/questions/1568511/how-do-i-sort-one-vector-based-on-values-of-another
    
    sorted_vector <- to_sort_vector[order(match(to_sort_vector,reference_vector))]
    return(sorted_vector[1])
}

# sort based on deleterious_intensity
variant_effects <- sapply(consequence_list, find_the_most_deleterious_variant)


# 3. Read depth at the site of variation.
DP_index <- vcf_df$FORMAT %>%
    strsplit(., ":") %>%
    sapply(., function(x) which(x == "DP")) %>%
    .[1]
read_depth_at_the_site_of_variations <- vcf_df$SL160967 %>%
    strsplit(., ":") %>%
    sapply(., function(x) x[DP_index]) %>%
    as.numeric


# 4. Number of reads supporting the variant.
AD_index <- vcf_df$FORMAT %>%
    strsplit(., ":") %>%
    sapply(., function(x) which(x == "AD")) %>%
    .[1]
number_of_reads_supporting_the_variants <- vcf_df$SL160967 %>%
    strsplit(., ":") %>%
    sapply(., function(x) x[AD_index]) %>%
    strsplit(., ",") %>%
    sapply(., function(x) x[2]) %>%
    as.numeric


# 5.Percentage of reads supporting the variant versus those supporting reference reads.
number_of_reads_supporting_the_references <- vcf_df$SL160967 %>%
    strsplit(., ":") %>%
    sapply(., function(x) x[AD_index]) %>%
    strsplit(., ",") %>%
    sapply(., function(x) x[1]) %>%
    as.numeric

percentage_of_reads_variants_vs_references <- number_of_reads_supporting_the_variants / number_of_reads_supporting_the_references


# 6. Allele frequency of variant
allele_frequency_of_variants <- vcf_df$INFO %>%
    strsplit(., ";") %>%
    sapply(., function(x) x[grepl("^AF", x)]) %>%
    gsub("AF=", "", .) %>%
    as.numeric


# 7. (Optional) Any other information from ExAC that you feel might be relevant 
get_exac_additional_info <- function(a_variant_id, all_exac_list = whole_exac_list, all_variant_effects = variant_effects) {
    
    if (is.null(all_exac_list[[a_variant_id]][["consequence"]])) {
        return(NA)
    }
    
    exac_info_df <- all_exac_list[[a_variant_id]] %>%
        .[["consequence"]] %>%
        .[[all_variant_effects[a_variant_id]]] %>%
        .[[1]]
    
    additional_information_from_exac <- paste(
        paste("Gene", paste(unique(exac_info_df$Gene), collapse = ","), sep = "="),
        paste("SYMBOL", paste(unique(exac_info_df$SYMBOL), collapse = ","), sep = "="),
        paste("HGNC_ID", paste(unique(exac_info_df$HGNC_ID), collapse = ","), sep = "="),
        paste("HGVSc", paste(unique(exac_info_df$HGVSc), collapse = ","), sep = "="),
        paste("HGVSp", paste(unique(exac_info_df$HGVSp), collapse = ","), sep = "="),
        paste("EXON", paste(unique(exac_info_df$EXON), collapse = ","), sep = "="),
        paste("Codons", paste(unique(exac_info_df$Codons), collapse = ","), sep = "="),
        paste("SIFT", paste(unique(exac_info_df$SIFT), collapse = ","), sep = "="),
        paste("PolyPhen", paste(unique(exac_info_df$PolyPhen), collapse = ","), sep = "="),
        sep = ";"
    )
    
    return(additional_information_from_exac)
}

extra_exac_informations <- sapply(variant_IDs, get_exac_additional_info)


# compile the extracted information:
# variant_IDs
# variant_types
# variant_effects
# read_depth_at_the_site_of_variations
# number_of_reads_supporting_the_variants
# percentage_of_reads_variants_vs_references
# allele_frequency_of_variants
# extra_exac_informations

summary_df <- data.frame(variant_id = variant_IDs,
                         variant_type = variant_types,
                         variant_effect = variant_effects,
                         read_depth = read_depth_at_the_site_of_variations,
                         read_supporting_variant = number_of_reads_supporting_the_variants,
                         read_ratio_var_vs_ref = percentage_of_reads_variants_vs_references,
                         variant_allele_freq = allele_frequency_of_variants,
                         ExAC_extra_info = extra_exac_informations)

write.table(summary_df, file = output_file, quote = F, sep = "\t", row.names = F)


proc.time() - ptm

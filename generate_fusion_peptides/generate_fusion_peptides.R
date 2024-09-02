"
Author Name: Helen Xie
Contact: helen_xie@brown.edu

Required packages:
1. data.table
2. dplyr
3. stringr
4. tidyr
5. Biostrings

Command line inputs: 
1. data_sv.txt - patient gene fusion data from MSK Impact
2. ucsc_refseq_hg19_mrna.txt - exon bed file for Chr37/hg19 
3. hg19.fa - hg19 fasta file (download from here: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/)

Output File:
data_fusions.rds - R data file containing column with the fusion peptide sequence 

data_fusions columns:
Sample_ID, Site1_Chromosome, Site2_Chromosome, Site1_Position, Site2_Position, 
Site1_Description, Site2_Description, Class, Event_Info, Connection_Type, 
Annotation, First_Transcript_ID, Second_Transcript_ID, First_Gene, Second_Gene, 
First_Exons, Second_Exons, First_Sequence, Second_Sequence, Fusion_Sequence, 
Peptide_Sequence, Peptide_Breakpoint, NetMHC_Peptide

"
# Load required packages
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)

# Set command line variables
data_sv <- args[1]
exon_bed_file <- args[2]
hg19_fasta <- args[3]

# STEP 1: Read and Parse the Fusion Data

#' Function to parse the Annotation column for transcript IDs
#' 
#' @param annotation: String containing the annotation with transcript IDs
#' 
#' @return: List of cleaned transcript IDs
parse_annotation <- function(annotation) {
  matches <- str_match_all(annotation, "\\(([^)]+)\\)")
  if (length(matches[[1]]) >= 2) {
    transcript_ids <- matches[[1]][1:2, 2]
    cleaned_transcript_ids <- sapply(list(transcript_ids), function(x) sub("\\..*$", "", x))
    return(list(cleaned_transcript_ids))
  } else {
    return(c(NA, NA))
  }
}

#' Function to determine which HugoSymbol comes first in the Annotation string
#' 
#' @param annotation: String containing the annotation
#' @param gene1: First gene name
#' @param gene2: Second gene name
#' 
#' @return: Gene name that appears first in the annotation
find_first_gene <- function(annotation, gene1, gene2) {
  pos_gene1 <- str_locate(annotation, gene1)[1]
  pos_gene2 <- str_locate(annotation, gene2)[1]
  
  if (is.na(pos_gene1)) pos_gene1 <- Inf
  if (is.na(pos_gene2)) pos_gene2 <- Inf
  
  if (pos_gene1 < pos_gene2) {
    return(gene1)
  } else {
    return(gene2)
  }
}

# Filter fusion data
data_fusions <- data_sv %>% 
  filter(grepl("Protein fusion: in frame", Event_Info)) %>% 
  filter(
    str_count(Annotation, "\\(N([^)]+)\\)") >= 2 |
      str_count(Comments, "\\(N([^)]+)\\)") >= 2
  ) %>%
  mutate(
    Annotation = ifelse(
      !str_detect(Annotation, "\\(N([^)]+)\\)"),  # Check if Annotation lacks the pattern
      Comments,  # Use Comments if pattern is not found
      Annotation  # Otherwise keep Annotation as is
    ),
    Annotation = ifelse(
      str_count(Comments, "exon") >= 2 & str_count(Comments, "\\(N([^)]+)\\)")  >= 2,  # Check if Comments has 2 instances of "exon"
      Comments,  # Replace Annotation with Comments if true
      Annotation  # Otherwise keep Annotation as is
    )
  )

# Add columns for parsed exon ranges and transcript IDs
data_fusions <- data_fusions %>%
  rowwise() %>%
  mutate(
    Transcript_IDs = parse_annotation(Annotation),
    Second_Transcript_ID = Transcript_IDs[2],
    First_Transcript_ID = Transcript_IDs[1],
    First_Gene = str_match(Event_Info, "[{(](.+?)[:-]")[,2],
    Second_Gene = str_match(Event_Info, "[{(][^:]+[:-](.+?)[})]")[,2]
  )

# STEP 2: Read and Parse the Bed File

# Rename columns for clarity
colnames(exon_bed_file) <- c("chr", "start", "end", "info", "score", "strand")

#' Function to parse transcript ID and exon number from bed file info column
#' 
#' @param info: String containing the BED file info
#' 
#' @return: Vector with transcript ID and exon number
parse_bed_info <- function(info) {
  parts <- strsplit(info, "_")[[1]]
  transcript_id <- paste(parts[1:2], collapse = "_")
  exon_number <- as.numeric(parts[4]) + 1 
  return(c(transcript_id, exon_number))
}

# Add columns for transcript_id and exon_number
exon_bed_file <- exon_bed_file %>%
  mutate(
    parsed_info = lapply(info, parse_bed_info),
    transcript_id = sapply(parsed_info, `[`, 1),
    exon_number = sapply(parsed_info, `[`, 2)
  ) %>%
  select(-parsed_info)

# Run bedtools getfasta on entire exon_bed_file 
system(paste("bedtools getfasta -fi", hg19_fasta, "-bed", exon_bed_file, "-fo", 'entire.fa')) 

entire_fasta <- readDNAStringSet("entire.fa")

# Convert the DNAStringSet to a data frame
fasta_df <- tibble(
  name = names(entire_fasta),
  sequence = as.character(entire_fasta)
)

# Merge the BED file data with the FASTA sequences
exon_bed_file <- exon_bed_file %>%
  mutate(sequence = fasta_df$sequence, name = fasta_df$name)

# STEP 3: Determine Upper Limit of Exons

#' Function to get the upper limit of exons for a given transcript ID
#' 
#' @param transcript: Transcript ID
#' @param chrom1: First chromosome
#' @param chrom2: Second chromosome
#' 
#' @return: Maximum exon number for the given transcript
get_upper_limit_exon <- function(transcript, chrom1, chrom2) {
  max_exon <- exon_bed_file %>%
    filter(transcript_id == transcript) %>% 
    filter(chr == paste0("chr", as.character(chrom1)) | chr == paste0("chr", as.character(chrom2))) %>%
    summarize(max_exon = max(as.numeric(exon_number))) %>%
    pull(max_exon)
  return(max_exon)
}

# Add upper limit of exons to fusion data
data_fusions <- data_fusions %>%
  rowwise() %>%
  mutate(
    First_Max_Exon = mapply(get_upper_limit_exon, First_Transcript_ID, Site1_Chromosome, Site2_Chromosome),
    Second_Max_Exon = mapply(get_upper_limit_exon, Second_Transcript_ID, Site1_Chromosome, Site2_Chromosome)
  )

# STEP 4: Parse Exon Ranges and Subset Bed File for Each Fusion

#' Function to extract exon range from a string
#' 
#' @param input_string: String containing exon information
#' @param instance: Instance of the match to extract
#' 
#' @return: Vector of exon numbers
extract_exon_range <- function(input_string, instance = 1) {
  pattern <- "(exon[s]?)\\s*(\\d+)(\\s*[-toand]+\\s*(\\d+))?"
  matches <- gregexpr(pattern, input_string, perl = TRUE)
  extracted_matches <- regmatches(input_string, matches)
  if (length(extracted_matches[[1]]) < instance) {
    return(NULL)
  }
  match <- extracted_matches[[1]][instance]
  numbers <- as.numeric(unlist(regmatches(match, gregexpr("\\d+", match, perl = TRUE))))
  if (length(numbers) == 1) {
    return(numbers:numbers)
  }
  if (length(numbers) == 2) {
    return(numbers[1]:numbers[2])
  }
  return(NULL)
}

#' Function to create exon ranges based on description and maximum exon number
#' 
#' @param description: Site description for first gene in the fusion
#' 
#' @return: Vector of exon numbers
create_exon_range_1 <- function(description) {
  parts <- strsplit(description, " ")[[1]]
  if (length(parts) < 7) {
    return(NULL)
  }
  before_or_after <- parts[5]
  exon_number <- as.numeric(parts[7])
  #gene_orientation <- substr(parts[3], nchar(parts[3]) - 2, nchar(parts[3]) - 2)
  
  if (before_or_after == "before"){
    return(1:(exon_number - 1))
  } else if (before_or_after == "after") {
    return(1:exon_number)
  } else {
    return(NULL)
  }
}



#' Function to create exon ranges based on description and maximum exon number
#' 
#' @param description Site description for second gene in the fusion
#' @param max_exon The maximum exon number for the transcript.
#' 
#' @return A vector containing the exon numbers in the range.
create_exon_range_2 <- function(description, max_exon) {
  parts <- strsplit(description, " ")[[1]]
  if (length(parts) < 7) {
    return(NULL)
  }
  before_or_after <- parts[5]
  exon_number <- as.numeric(parts[7])
  if (before_or_after == "before"){
    return(exon_number:max_exon)
  } else if (before_or_after == "after") {
    return((exon_number + 1):max_exon)
  } else {
    return(NULL)
  }
}

#' Function to parse exon ranges for the first gene
#' 
#' @param annotation A string containing the annotation data.
#' @param comments A string containing the comments data.
#' @param site_description A string containing the site description.
#' 
#' @return A vector containing the exon range for the first gene.
parse_exon_range_gene1 <- function(annotation, comments, site_description) {
  if (length(gregexpr("exon", annotation, ignore.case = TRUE)[[1]]) >= 2) {
    return(extract_exon_range(annotation, 1))
  } else if (length(gregexpr("exon", comments, ignore.case = TRUE)[[1]]) >= 2) {
    return(extract_exon_range(comments, 1))
  } else {
    return(create_exon_range_1(site_description))
  }
}

#' Function to parse exon ranges for the second gene
#' 
#' @param annotation A string containing the annotation data.
#' @param comments A string containing the comments data.
#' @param site_description A string containing the site description.
#' @param max_exon The maximum exon number for the transcript.
#' 
#' @return A vector containing the exon range for the second gene.
parse_exon_range_gene2 <- function(annotation, comments, site_description, max_exon) {
  if (length(gregexpr("exon", annotation, ignore.case = TRUE)[[1]]) >= 2) {
    return(extract_exon_range(annotation, 2))
  } else if (length(gregexpr("exon", comments, ignore.case = TRUE)[[1]]) >= 2) {
    return(extract_exon_range(comments, 2))
  } else {
    return(create_exon_range_2(site_description, max_exon))
  }
}

# Add columns for exon ranges for both genes
data_fusions <- data_fusions %>%
  mutate(
    First_Exons = ifelse(
      grepl(First_Gene, Site1_Description), 
      sapply(1:nrow(data_fusions), function(i) parse_exon_range_gene1(Annotation[i], Comments[i], Site1_Description[i])),
      sapply(1:nrow(data_fusions), function(i) parse_exon_range_gene1(Annotation[i], Comments[i], Site2_Description[i]))),
    Second_Exons = ifelse(
      grepl(Second_Gene, Site1_Description), 
      sapply(1:nrow(data_fusions), function(i) parse_exon_range_gene2(Annotation[i], Comments[i], Site1_Description[i], First_Max_Exon[i])),
      sapply(1:nrow(data_fusions), function(i) parse_exon_range_gene2(Annotation[i], Comments[i], Site2_Description[i], Second_Max_Exon[i]))
    )
  )

# STEP 5: Process fusion data from ARCHER 
# ARCHER data processing
extract_transcript_ids_archer <- function(annotation) {
  matches <- str_match_all(annotation, "\\(([^)]+)\\)")
  if (length(matches[[1]]) >= 2) {
    transcript_ids <- matches[[1]][1:2, 2]
    cleaned_transcript_ids <- sapply(list(transcript_ids), function(x) sub("\\..*$", "", x))
    return(list(cleaned_transcript_ids))
    return(list(transcript_ids))
  } else {
    return(c(NA, NA))
  }
}

extract_first_exons_archer <- function(annotation) {
  exons <- str_extract_all(annotation, "Exon\\d+")[[1]]
  exons <- as.numeric(sub("Exon", "", exons))
  first_exons <- 1:exons[1]
  return(first_exons)
}

extract_second_exons_archer <- function(annotation, max_exon) {
  exons <- str_extract_all(annotation, "Exon\\d+")[[1]]
  exons <- as.numeric(sub("Exon", "", exons))
  second_exons <- exons[2]:max_exon
  return(second_exons)
}

extract_genes_archer <- function(fusion_string) {
  parts <- strsplit(fusion_string, "-")[[1]]
  first_gene <- parts[1]
  second_gene <- strsplit(parts[2], " ")[[1]][1]
  return(list(first_gene, second_gene))
}

data_fusions_archer <- data_sv %>% 
  filter(grepl("Archer", Comments) & grepl("POSITIVE", Annotation) & grepl("in-frame fusion", Annotation)) %>%
  mutate(Transcript_IDs = extract_transcript_ids_archer(Annotation),
         First_Transcript_ID = Transcript_IDs[1],
         Second_Transcript_ID = Transcript_IDs[2],
         First_Gene = extract_genes_archer(Annotation)[1],
         Second_Gene = extract_genes_archer(Annotation)[2],
         First_Max_Exon = mapply(get_upper_limit_exon, First_Transcript_ID, Site1_Chromosome, Site2_Chromosome),
         Second_Max_Exon = mapply(get_upper_limit_exon, Second_Transcript_ID, Site1_Chromosome, Site2_Chromosome),
         First_Exons = extract_first_exons_archer(Annotation),
         Second_Exons = extract_second_exons_archer(Annotation, Second_Max_Exon))

data_fusions <- rbind(data_fusions, data_fusions_archer)

# STEP 6: Retrieve DNA and Peptide Sequences using Bedtools

#' Function to retrieve FASTA sequences from the bed file
#' 
#' @param df A dataframe containing the fusion data.
#' @param is_first A boolean indicating whether to retrieve the first or second sequence.
#' @param num The row number of the fusion data.
#' 
#' @return A string containing the concatenated DNA sequence for the exons.
get_fasta_sequences <- function(df, is_first, num) {
  print(num)
  if (is_first == TRUE) {
    subset <- exon_bed_file %>%
      rowwise()%>%
      mutate(exon_number = as.numeric(exon_number))%>%
      filter(transcript_id == df$First_Transcript_ID[num] & exon_number %in% df$First_Exons[[num]]) 
  } else{
    subset <- exon_bed_file %>%
      rowwise()%>%
      mutate(exon_number = as.numeric(exon_number))%>%
      filter(transcript_id == df$Second_Transcript_ID[num] & exon_number %in% df$Second_Exons[[num]])  
  }
  if (nrow(subset) == 0) {
    return("")
  }
  concatenated_sequence <- paste(subset$sequence, collapse = "")  
  if (any(subset$strand == "-")) {
    concatenated_sequence <- as.character(reverseComplement(DNAString(concatenated_sequence)))
  }
  return(concatenated_sequence)
}

#' Function to convert DNA sequences to peptide sequences
#' 
#' @param dna_string A string containing the DNA sequence.
#' 
#' @return A string containing the translated peptide sequence.
dna_to_peptide <- function(dna_string) {
  return(as.character(translate(DNAString(dna_string))))
}

#' Function to get sequences for the fusion data
#' 
#' @param fusion_data A dataframe containing the fusion data.
#' 
#' @return A dataframe containing the fusion data with added DNA and peptide sequences.
get_sequences <- function(fusion_data) {
  # Initialize empty columns
  fusion_data <- fusion_data %>%
    mutate(First_Sequence = "",
           Second_Sequence = "",
           Fusion_Sequence = "",
           Peptide_Sequence = "")
  
  for (i in seq_len(nrow(fusion_data))) {
    # Get sequences for each row
    first_sequence <- get_fasta_sequences(fusion_data, TRUE, i)
    second_sequence <- get_fasta_sequences(fusion_data, FALSE, i)
    fusion_sequence <- paste0(first_sequence, second_sequence)
    peptide_sequence <- dna_to_peptide(fusion_sequence)
    
    # Update the data frame with the results
    fusion_data$First_Sequence[i] <- first_sequence
    fusion_data$Second_Sequence[i] <- second_sequence
    fusion_data$Fusion_Sequence[i] <- fusion_sequence
    fusion_data$Peptide_Sequence[i] <- peptide_sequence
  }
  
  return(fusion_data)
}

data_fusions <- get_sequences(data_fusions)

# STEP 7: Get 20mer Peptide Sequence around Breakpoint for input around NetMHC

#' Function to calculate the peptide breakpoint position
#' 
#' @param dna_breakpoint_location An integer indicating the DNA breakpoint location.
#' 
#' @return An integer indicating the peptide breakpoint position.
calculate_peptide_breakpoint <- function(dna_breakpoint_location) {
  breakpoint_pos <- round(dna_breakpoint_location / 3)
  return(breakpoint_pos)
}

#' Function to get the surrounding peptides around the breakpoint
#' 
#' @param breakpoint_position An integer indicating the peptide breakpoint position.
#' @param peptide_sequence A string containing the peptide sequence.
#' 
#' @return A string containing the surrounding peptides.
get_surrounding_peptides <- function(breakpoint_position, peptide_sequence) {
  start_pos <- max(1, breakpoint_position - 10)
  end_pos <- min(nchar(peptide_sequence), breakpoint_position + 10)
  surrounding_peptides <- substr(peptide_sequence, start_pos, end_pos)
  return(surrounding_peptides)
}

# Add columns for: peptide breakpoint location, full peptide sequence with separator at breakpoint, and 20mer peptide sequence surrounding breakpoint
data_fusions <- data_fusions %>%
  rowwise() %>%
  mutate(Peptide_Breakpoint = calculate_peptide_breakpoint(nchar(First_Sequence))) %>%
  mutate(Peptide_Breakpoint_Sequence = paste0(substr(Peptide_Sequence, 1, Peptide_Breakpoint - 1), "|", substr(Peptide_Sequence, Peptide_Breakpoint, nchar(Peptide_Sequence)))) %>%
  mutate(NetMHC_Peptide = get_surrounding_peptides(Peptide_Breakpoint, Peptide_Sequence[[1]]))

# Save the resulting data to an RDS file
saveRDS(data_fusions, file = "data_fusions.rds")
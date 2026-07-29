#' Comprehensive Marker Gene Taxonomic Resolution Pipeline
#'
#' This function downloads genomes from NCBI based on targeted genera, dereplicates them, and uses .gff 
#' coordinates to extract genes from .fasta files. It aligns the 
#' sequences, performs in-silico PCR on gene regions based on primers inputed, and generates 
#' phylogenetic resolution metrics as well as a primer mismatch report to help with 
#' primer selection and taxonomic assignment analysis.
#'
#' @param target_genera A character vector of genera to analyze (e.g., c("Gilliamella", "Snodgrassella")) or input a .csv file with a list of genera and all will be processed with a log file to summarize resolution of primers for each genus.
#' @param output_dir String. Path to the folder where results should be saved. Defaults to current working directory. Each genera targeted generates a new folder in the output directory with all output files
#' @param db_dir String. Path to the directory where master databases (NCBI summary, LPSN) are stored. Defaults to current working directory. Setting this to a static folder saves your laptop's bandwidth by preventing the pipeline from re-downloading the databases for every new output directory.
#' @param mafft_path String. Path to MAFFT executable (default: "mafft"). If MAFFT is not found, the pipeline automatically falls back to DECIPHER. Set to "" to intentionally force DECIPHER and bypass the MAFFT system check.
#' @param target_gene String. The specific gene to extract from the GFF file (default: "16S"). RefSeq has a standardized annotation pipeline, so if having issues download a genomes GFF file and look through it to find the spelling for your gene of interest.
#' @param feature_type String. The GFF feature category to filter by (e.g., "rRNA", "CDS", "gene"). Set to "ANY" to bypass this filter (default: "rRNA"). Multiple genes have 16S in the name, but only one rRNA has 16S in the name, so use this to limit your search to the correct catagory.
#' @param custom_primers Optional. Path to a CSV file or an R list containing custom primer sets. 
#' @param only_reference Logical. If TRUE, only uses RefSeq/Representative genomes. Some RefSeq/Representative genomes are not validily published according to LPSN, so enable that setting as well if wanting only validly published reference species. TRUE is recommended as some genera have many poor quality genomes.
#' @param dereplicate_strains Logical. Keep only the best genome per strain.
#' @param remove_unclassified Logical. Remove unclassified strains like "sp." or "indicum".
#' @param enable_lpsn_check Logical. Validate names against LPSN database. HIGHLY RECOMMENDED* there are many synonym names in the refseq database which will give false negatives about a taxa being monophyletic
#' @param lpsn_db_path String. Path to the LPSN database CSV file. Can be downloaded at https://lpsn.dsmz.de/downloads
#' @param max_contigs Numeric. Maximum allowed contigs for draft genomes. RefSeq representative genomes are always kept regardless of contig #.
#' @param refseq_max_age Numeric. Maximum age (in days) of the local RefSeq summary file before a new one is downloaded. Set to Inf to always use the local file if it exists, or 0 to force a fresh download. Default is 30.
#' @param max_tax_level String. The highest taxonomic level to assess for clade resolution (e.g., "Genus", "Family", "Order", "Class", "Phylum"). Computation resources and time increase dramatically with increase taxa levels.
#' @param n_threats Numeric. The number of closest outgroup genera to pull full species data and align. Defaults to 2. If wanting to look at all species in a Family or higher taxonomic level, set to Inf and make sure max_scout_genera is also set to Inf
#' @param max_scout_genera Numeric. Maximum number of outgroup genera to fetch during the phylogenetic scout phase. Defaults to 50. Set to Inf for unlimited. A setting to tune as needed based on time and computational resources.
#' @param keep_genomes Logical. If TRUE (default), retains the downloaded .fasta and .gff files. If FALSE, deletes them after the run to save disk space. Unless you have a lot of extra space FALSE is recommended when going above the genus level.
#' 
#' @importFrom dplyr %>% mutate filter arrange desc select group_by slice ungroup case_when left_join bind_rows n n_distinct relocate summarize
#' @importFrom ggplot2 ggplot aes labs theme_bw theme geom_bar geom_text geom_area geom_rect scale_fill_manual scale_y_continuous coord_flip element_text ggsave
#' @importFrom stringr str_extract word
#' @importFrom Biostrings readDNAStringSet writeXStringSet DNAStringSet DNAString reverseComplement matchPattern matchLRPatterns subseq width
#' @importFrom IRanges subject
#' @importFrom DECIPHER RemoveGaps DistanceMatrix AlignSeqs
#' @importFrom ape njs
#' @importFrom ggtree ggtree geom_tiplab theme_tree2 hexpand
#' @importFrom curl curl_download
#' @importFrom digest digest
#' @importFrom forcats fct_relevel
#' @return Generates output folders containing alignments, trees, mismatch reports, and resolution plots.
#' @export
run_marker_pipeline <- function(target_genera = c("Commensalibacter", "Apilactobacillus"), 
                             output_dir = ".",
                             db_dir = ".",
                             mafft_path = "mafft",
                             target_gene = "16S",
                             feature_type = "rRNA",
                             custom_primers = NULL,
                             only_reference = TRUE,
                             dereplicate_strains = TRUE,
                             remove_unclassified = TRUE,
                             enable_lpsn_check = TRUE,
                             lpsn_db_path = "lpsn_gss.csv",
                             max_contigs = 100,
                             refseq_max_age = 30,
                             max_tax_level = "Genus",
                             n_threats = 2,
                             max_scout_genera = Inf,
                             keep_genomes = TRUE) {
  
  # Ensure the base output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # --- GENERA INPUT PARSING ---
  if (length(target_genera) == 1 && grepl("\\.csv$", target_genera, ignore.case = TRUE)) {
    if (!file.exists(target_genera)) stop("Error: The specified genera CSV file does not exist.")
    message(paste("  -> Loading target genera from file:", basename(target_genera)))
    genera_df <- read.csv(target_genera, stringsAsFactors = FALSE)
    target_genera <- trimws(genera_df[[1]]) 
    target_genera <- target_genera[target_genera != ""] 
  }
  
  # PRIMER CONFIGURATION
  default_primers <- list(
    "V1_V2_Lane_1991"              = list(Fwd="AGAGTTTGATCMTGGCTCAG", Rev="TGCTGCCTCCCGTAGGAGT", Min=280, Max=420),
    "V1_V3_Muyzer_1993"            = list(Fwd="AGAGTTTGATCMTGGCTCAG", Rev="ATTACCGCGGCTGCTGG", Min=450, Max=600),
    "V3_V4_Klindworth_2013"        = list(Fwd="CCTACGGGNGGCWGCAG", Rev="GACTACHVGGGTATCTAATCC", Min=400, Max=550),
    "V3_V4_Takahashi_2014"         = list(Fwd="CCTACGGGNGGCWGCAG", Rev="GGACTACNVGGGTWTCTAAT", Min=400, Max=550),
    "V3_V6_Huber_2007"             = list(Fwd="CCTACGGGAGGCAGCAG", Rev="ACGAGCTGACGACARCCATG", Min=620, Max=850),
    "V4_EMP_Parada_Apprill_2016"   = list(Fwd="GTGYCAGCMGCCGCGGTAA", Rev="GGACTACNVGGGTWTCTAAT", Min=230, Max=350),
    "V4_V5_Parada_2016"            = list(Fwd="GTGYCAGCMGCCGCGGTAA", Rev="CCGYCAATTYMTTTRAGTTT", Min=350, Max=480),
    "V5_V7_Engelbrektson_2010"     = list(Fwd="AACMGGATTAGATACCCKG", Rev="ACGTCATCCCCACCTTCC", Min=350, Max=500),
    "V6_V8_Engelbrektson_2010"     = list(Fwd="ACGCGHNRAACCTTACC", Rev="ACGGGCRGTGWGTRCAA", Min=350, Max=500),
    "Full_16S_Callahan_2019"       = list(Fwd="AGRGTTYGATYMTGGCTCAG", Rev="TACGGYTACCTTGTTACGACTT", Min=1300, Max=1600)
  )
  
  # ==============================================================================
  # --- Loading Custom Primers ---
  # ==============================================================================
  if (is.null(custom_primers)) {
    primer_sets <- default_primers
    message("  -> Using default built-in primer sets.")
    
  } else if (is.character(custom_primers)) {
    # Check if it's a file path
    if (file.exists(custom_primers)) {
      message(paste("  -> Loading custom primers from file:", basename(custom_primers)))
      custom_df <- read.csv(custom_primers, stringsAsFactors = FALSE)
      primer_sets <- list()
      for(i in 1:nrow(custom_df)) {
        p_name <- as.character(custom_df$Primer_Name[i])
        primer_sets[[p_name]] <- list(
          Fwd = as.character(custom_df$Fwd_Seq[i]), 
          Rev = as.character(custom_df$Rev_Seq[i]), 
          Min = as.numeric(custom_df$Min_Length[i]), 
          Max = as.numeric(custom_df$Max_Length[i])
        )
      }
      message(paste("  -> Successfully loaded", length(primer_sets), "primer sets:", paste(names(primer_sets), collapse=", ")))
    } else {
      stop(paste("Error: Custom primer file not found at:", custom_primers))
    }
    
  } else if (is.list(custom_primers)) {
    message("  -> Using custom R list for primer sets.")
    primer_sets <- custom_primers
    
  } else {
    stop("Error: 'custom_primers' must be NULL, a valid CSV file path, or an R list.")
  }
  
# General Marker Gene & Primer Safety Check
  is_default_16s <- grepl("16S", target_gene, ignore.case = TRUE)
  if(!is_default_16s && is.null(custom_primers)) {
    message(paste("  -> [WARNING] Custom target gene detected (", target_gene, ")."))
    message("  -> [WARNING] No custom primers provided. Default 16S V-region primers are disabled.")
    message("  -> [NOTE] Pipeline will proceed and calculate resolution for the Full-Length gene only.")
    primer_sets <- list() # Empties the list so the script safely skips the PCR loops!
  }

  # =======================
  # GLOBAL DATA LOAD 
  # =======================
  message("--- PRE-FLIGHT: Loading Master Databases ---")
  
  # A. FETCH NCBI SUMMARY
  dest_file <- file.path(db_dir, "assembly_summary_refseq.txt")
  should_download <- TRUE
  
  if(file.exists(dest_file)) {
    file_age <- as.numeric(difftime(Sys.time(), file.info(dest_file)$mtime, units = "days"))
    if(file_age <= refseq_max_age) {
      message(paste0("  -> Using existing NCBI summary (", round(file_age, 1), " days old)."))
      should_download <- FALSE
    } else {
      message(paste0("  -> Existing NCBI summary is too old (", round(file_age, 1), " days). Downloading fresh..."))
      unlink(dest_file)
    }
  }
  
  if(should_download) {
    message("  -> Downloading fresh RefSeq summary (this may take a minute)...")
    options(timeout = max(600, getOption("timeout"))) 
    url <- "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
    tryCatch({ download.file(url, destfile = dest_file, method = "libcurl", quiet = TRUE) }, 
             error = function(e) { stop(paste("\nNCBI Download Failed:\n", e$message)) })
  }
  
  raw_meta <- read.table(dest_file, sep = "\t", header = FALSE, comment.char = "#", quote = "", fill = TRUE, stringsAsFactors = FALSE, skip = 2)
  
  master_meta <- data.frame(
    assembly_accession = raw_meta[,1],
    refseq_category    = raw_meta[,5],
    organism_name      = raw_meta[,8],
    infraspecific_name = raw_meta[,9],
    assembly_level     = raw_meta[,12],
    seq_rel_date       = raw_meta[,15],
    ftp_path           = raw_meta[,20],
    contig_count       = as.numeric(raw_meta[,31]), 
    stringsAsFactors = FALSE
  )
  
 # B. LOAD LPSN
  master_lpsn_clean <- NULL
  if(enable_lpsn_check) {
    # Dynamically route to db_dir if only a filename is provided
    actual_lpsn_path <- if (basename(lpsn_db_path) == lpsn_db_path) file.path(db_dir, lpsn_db_path) else lpsn_db_path
    
    if(file.exists(actual_lpsn_path)) {
      message(paste("  -> Loading LPSN database:", basename(actual_lpsn_path)))
      lpsn <- tryCatch(read.csv(actual_lpsn_path, stringsAsFactors=FALSE, fill=TRUE), error=function(e) NULL)
      if(!is.null(lpsn)) {
        sp_col <- if("sp_epithet" %in% colnames(lpsn)) "sp_epithet" else "species_epithet"
        master_lpsn_clean <- lpsn %>%
          mutate(
            genus_clean = trimws(genus_name),
            species_clean = trimws(ifelse(is.na(.data[[sp_col]]), "", .data[[sp_col]])),
            Full_Name = trimws(gsub("\\s+", " ", paste(genus_clean, species_clean))),
            status_clean = tolower(status)
          ) %>%
          mutate(is_valid = grepl("validly published|correct name", status_clean)) %>%
          group_by(Full_Name) %>% arrange(desc(is_valid)) %>% dplyr::slice(1) %>% ungroup() %>% select(-is_valid)
      }
    } else {
      message(paste("  [WARNING] LPSN file not found at:", actual_lpsn_path, "- Check skipped."))
      enable_lpsn_check <- FALSE
    }
  }
  
  # =================
  # MAIN BATCH LOOP
  # =================
  master_resolution_log <- list()
  master_higher_taxa_log <- list()
  
  for (target_genus in target_genera) {
    message(paste("\n======================================================================"))
    message(paste("  INITIATING PIPELINE FOR:", toupper(target_genus)))
    message(paste("======================================================================\n"))
    
    OUT_DIR <- file.path(output_dir, paste0("analysis_refseq_", target_genus, "_", target_gene))
    if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR)
    
    # Establish Organized Subdirectories 
    dir_genomes    <- file.path(OUT_DIR, "01_Genomes")
    dir_annots     <- file.path(OUT_DIR, "02_Annotations")
    dir_alignments <- file.path(OUT_DIR, "03_Alignments")
    dir_trees      <- file.path(OUT_DIR, "04_Trees")
    dir_entropy    <- file.path(OUT_DIR, "05_Entropy_Maps")
    dir_logs       <- file.path(OUT_DIR, "06_Logs")
    dir_results    <- file.path(OUT_DIR, "07_Resolution_Results")
    
    for(d in c(dir_genomes, dir_annots, dir_alignments, dir_trees, dir_entropy, dir_logs, dir_results)) {
      if(!dir.exists(d)) dir.create(d)
    }
    
    # Metadata Filtering & QC
    message("--- Phase 3: Metadata & QC ---")
    meta <- master_meta %>% 
      mutate(organism_name = gsub("\\[|\\]", "", organism_name)) %>%
      filter(grepl(paste0("^", target_genus, "(\\s|$)"), organism_name, ignore.case = TRUE))
    
    if(nrow(meta) == 0) {
      message(paste("  -> [SKIP] No genomes found for", target_genus))
      next
    }
    
    if(only_reference) meta <- meta %>% filter(refseq_category %in% c("reference genome", "representative genome"))
    
    filtered_meta <- meta %>% 
      filter(assembly_level %in% c("Complete Genome", "Chromosome", "Scaffold", "Contig")) %>%
      filter(!is.na(ftp_path) & ftp_path != "na") %>%
      mutate(
        Strain_Label = gsub("strain=", "", ifelse(is.na(infraspecific_name)|infraspecific_name=="na", "Unknown", infraspecific_name)),
        Score_RefSeq = case_when(refseq_category=="reference genome"~3, refseq_category=="representative genome"~2, TRUE~1),
        Score_Level  = case_when(assembly_level=="Complete Genome"~4, assembly_level=="Chromosome"~3, assembly_level=="Contig"~2, TRUE~1),
        Is_VIP = refseq_category %in% c("reference genome", "representative genome")
      )
    
    if(remove_unclassified) filtered_meta <- filtered_meta %>% filter(!grepl(" sp\\.| sp$|Candidatus|uncultured", organism_name, ignore.case=TRUE))
    
####      LPSN CHECK
    if(enable_lpsn_check && !is.null(master_lpsn_clean)) {
      message("  -> Checking species validity against LPSN...")
      filtered_meta <- filtered_meta %>% mutate(Check_Name = gsub("\\[|\\]", "", word(organism_name, 1, 2)))
      validation_df <- filtered_meta %>% left_join(master_lpsn_clean %>% select(Full_Name, status_clean), by = c("Check_Name" = "Full_Name")) %>%
        mutate(
          Is_Not_Found = is.na(status_clean),
          Is_Not_Validly_Pub = grepl("not validly published", status_clean, fixed = TRUE),
          Is_Synonym   = grepl("synonym", status_clean, fixed = TRUE),
          Is_Correct   = grepl("correct name", status_clean, fixed = TRUE),
          Is_Valid_Pub = grepl("validly published", status_clean, fixed = TRUE),
          Is_Valid_LPSN = (Is_Correct | Is_Valid_Pub) & !Is_Synonym & !Is_Not_Validly_Pub,
          Failure_Reason = case_when(Is_Valid_LPSN ~ "Valid", Is_Not_Found ~ "Name_Not_Found_In_LPSN", Is_Not_Validly_Pub ~ "Not_Validly_Published", Is_Synonym ~ "Synonym_Old_Name", TRUE ~ paste0("Other_Status: ", status_clean))
        )

      filtered_meta$Is_Valid_LPSN <- validation_df$Is_Valid_LPSN
      filtered_meta$Is_Valid_LPSN[is.na(filtered_meta$Is_Valid_LPSN)] <- FALSE
      
      invalid <- filtered_meta %>% filter(!Is_Valid_LPSN)
      if(nrow(invalid) > 0) {
        message(paste("  -> [WARNING]", nrow(invalid), "genomes failed LPSN checks. VIP status revoked."))
        write.csv(validation_df %>% filter(!Is_Valid_LPSN) %>% select(assembly_accession, organism_name, Check_Name, Failure_Reason), file.path(dir_logs, "LPSN_Invalid_Species_Report.csv"), row.names=FALSE)
        
        # Revoke VIP status so they get caught by the contig filter!
        filtered_meta <- filtered_meta %>% mutate(Is_VIP = ifelse(!Is_Valid_LPSN, FALSE, Is_VIP))
      }
      filtered_meta$Check_Name <- NULL; filtered_meta$Is_Valid_LPSN <- NULL
    }

    # --- EXISTING QC LOGIC RESUMES HERE ---
    filtered_meta$contig_count[is.na(filtered_meta$contig_count)] <- 0
    filtered_meta$QC_Status <- ifelse((filtered_meta$contig_count <= max_contigs) | filtered_meta$Is_VIP, "Pass", "Drop_HighContigs")
    write.csv(filtered_meta %>% select(assembly_accession, organism_name, contig_count, QC_Status), file.path(dir_logs, "QC_Contig_Report.csv"), row.names = FALSE)
    filtered_meta <- filtered_meta %>% filter(QC_Status == "Pass")

    filtered_meta$contig_count[is.na(filtered_meta$contig_count)] <- 0
    filtered_meta$QC_Status <- ifelse((filtered_meta$contig_count <= max_contigs) | filtered_meta$Is_VIP, "Pass", "Drop_HighContigs")
    write.csv(filtered_meta %>% select(assembly_accession, organism_name, contig_count, QC_Status), file.path(dir_logs, "QC_Contig_Report.csv"), row.names = FALSE)
    filtered_meta <- filtered_meta %>% filter(QC_Status == "Pass")
    
    if(dereplicate_strains) {
      filtered_meta <- filtered_meta %>%
        mutate(Derep_Key = ifelse(Strain_Label %in% c("Unknown", "Strain_Unknown"), assembly_accession, Strain_Label)) %>%
        group_by(organism_name, Derep_Key) %>%
        arrange(desc(Score_RefSeq), desc(Score_Level), desc(seq_rel_date)) %>%
        dplyr::slice(1) %>% ungroup() %>% dplyr::select(-Derep_Key)
    }
    
    # Generate Clean Output Names
    filtered_meta <- filtered_meta %>%
      mutate(
        Short_Name = word(organism_name, 1, 2),
        Clean_Strain = ifelse(is.na(Strain_Label) | Strain_Label %in% c("Strain_Unknown", "Unknown", "na"), "", Strain_Label),
        Display_Name = paste(Short_Name, Clean_Strain, assembly_accession, sep="_"),
        Display_Name = gsub("_+", "_", Display_Name),
        Display_Name = gsub("^_|_$", "", Display_Name),
        Display_Name = gsub("[^a-zA-Z0-9_\\.]", "_", Display_Name) # FASTA Safe
      )
    
    filtered_meta$local_path <- NA; filtered_meta$local_gff <- NA
    message(paste("  -> Processing", nrow(filtered_meta), "genomes for", target_genus, "..."))
    
    for(i in 1:nrow(filtered_meta)) {
      acc <- filtered_meta$assembly_accession[i]
      ftp_base <- filtered_meta$ftp_path[i]
      folder_name <- basename(ftp_base)
      base_name <- filtered_meta$Display_Name[i]
      
      dest_genomic <- file.path(dir_genomes, paste0(base_name, ".fna.gz"))
      dest_gff     <- file.path(dir_annots, paste0(base_name, ".gff.gz"))
      
      if(!file.exists(dest_genomic)) tryCatch({ curl::curl_download(paste0(ftp_base, "/", folder_name, "_genomic.fna.gz"), dest_genomic, quiet=TRUE) }, error=function(e) {})
      if(!file.exists(dest_gff)) tryCatch({ curl::curl_download(paste0(ftp_base, "/", folder_name, "_genomic.gff.gz"), dest_gff, quiet=TRUE) }, error=function(e) {})
      
      if(file.exists(dest_genomic)) filtered_meta$local_path[i] <- dest_genomic
      if(file.exists(dest_gff)) filtered_meta$local_gff[i] <- dest_gff
      if(i %% 10 == 0) cat(".")
    }
    cat("\n")
    
    full_metadata <- filtered_meta %>% filter(!is.na(local_path)) %>% mutate(Isolation_Source="Unknown", Host="Unknown", Country="Unknown")
    
    # Fetch Full Taxonomy Lineage
    message(paste("  -> Fetching full NCBI taxonomy for", target_genus, "..."))
    tax_info <- tryCatch({
      # Search NCBI Taxonomy for the target genus
      search_res <- rentrez::entrez_search(db="taxonomy", term=paste0(target_genus, "[Scientific Name]"))
      
      if(length(search_res$ids) > 0) {
        # Fetch the XML lineage data
        tax_xml <- rentrez::entrez_fetch(db="taxonomy", id=search_res$ids[1], rettype="xml")
        xml_doc <- xml2::read_xml(tax_xml)
        
        # Extract the taxonomic ranks and names
        nodes <- xml2::xml_find_all(xml_doc, "//LineageEx/Taxon")
        ranks <- xml2::xml_text(xml2::xml_find_all(nodes, ".//Rank"))
        names <- xml2::xml_text(xml2::xml_find_all(nodes, ".//ScientificName"))
        dict <- setNames(names, ranks)
        
        list(
          Phylum = ifelse("phylum" %in% names(dict), dict["phylum"], NA),
          Class = ifelse("class" %in% names(dict), dict["class"], NA),
          Order = ifelse("order" %in% names(dict), dict["order"], NA),
          Family = ifelse("family" %in% names(dict), dict["family"], NA)
        )
      } else { 
        list(Phylum=NA, Class=NA, Order=NA, Family=NA) 
      }
    }, error = function(e) { 
      message("  -> [WARNING] Taxonomy fetch failed (Check internet or NCBI status).")
      list(Phylum=NA, Class=NA, Order=NA, Family=NA) 
    })

    # Enrich Metadata with Full Lineage
    full_metadata <- full_metadata %>%
      mutate(
        Phylum = tax_info$Phylum,
        Class = tax_info$Class,
        Order = tax_info$Order,
        Family = tax_info$Family,
        Genus = target_genus,
        Species = gsub("\\[|\\]", "", stringr::word(organism_name, 1, 2))
      ) %>%
      relocate(Phylum, Class, Order, Family, Genus, Species, .after = organism_name)

    write.csv(full_metadata, file.path(dir_logs, "genome_metadata_refseq_enriched.csv"), row.names=FALSE)
    
   # Dynamic Gene Extractor
    extract_gene_manual <- function(genome_path, gff_path, display_name, search_term, feat_type) {
      if(is.na(gff_path) || !file.exists(gff_path)) return(NULL)
      if(is.na(genome_path) || !file.exists(genome_path)) return(NULL)
      con <- if (grepl("\\.gz$", gff_path)) gzfile(gff_path, "rt") else file(gff_path, "rt")
      all_lines <- tryCatch(readLines(con, warn=FALSE), error = function(e) NULL)
      close(con)
      if(is.null(all_lines)) return(NULL)

      # 1. Broad dynamic grep to find candidate lines quickly
      candidate_indices <- grep(search_term, all_lines, ignore.case=TRUE)
      if(length(candidate_indices) == 0) return(NULL)

      genome <- tryCatch(readDNAStringSet(genome_path), error=function(e) NULL)
      if(is.null(genome)) return(NULL)
      names(genome) <- sapply(strsplit(names(genome), " "), `[`, 1)

      extracted_seqs <- DNAStringSet()
      for (idx in candidate_indices) {
        line <- all_lines[idx]; if(startsWith(line, "#")) next
        parts <- strsplit(line, "\t")[[1]]
        
        if (length(parts) >= 9) {
          # 2. Apply the Feature Type Filter (Column 3 of GFF)
          pass_feature <- (feat_type == "ANY" || parts[3] == feat_type)
          
          if (pass_feature) {
            tags <- strsplit(parts[9], ";")[[1]]
            is_valid_target <- FALSE
            for(tag in tags) {
              tag <- trimws(tag)
              
              # 3. Verify the product, gene, or Name matches the search term
              if(startsWith(tag, "product=") || startsWith(tag, "gene=") || startsWith(tag, "Name=")) { 
                tag_parts <- strsplit(tag, "=")[[1]]
                if(length(tag_parts) >= 2) {
                  val <- tag_parts[2]
                  if (grepl(search_term, val, ignore.case=TRUE) || (search_term == "16S" && val %in% c("16S ribosomal RNA", "16S rRNA"))) { 
                    is_valid_target <- TRUE; break 
                  }
                }
              }
            }
            if (is_valid_target) {
              contig <- parts[1]; start_pos <- as.numeric(parts[4]); end_pos <- as.numeric(parts[5]); strand <- parts[7]
              if (contig %in% names(genome)) {
                seq_frag <- subseq(genome[[contig]], start=start_pos, end=end_pos)
                if (strand == "-") seq_frag <- reverseComplement(seq_frag)
                
                # Minimum size dropped to 200bp for smaller genes
                if (length(seq_frag) > 200) extracted_seqs <- c(extracted_seqs, DNAStringSet(seq_frag))
              }
            }
          }
        }
      }
      
      # 4. Dereplicate identical copies
      if(length(extracted_seqs) > 0) { 
        unique_idx <- !duplicated(as.character(extracted_seqs))
        extracted_seqs <- extracted_seqs[unique_idx]
        names(extracted_seqs) <- paste0(display_name, "_copy", 1:length(extracted_seqs))
        return(extracted_seqs) 
      } else {
        return(NULL)
      }
    }
    
    all_gene_seqs <- DNAStringSet()
    status_df <- data.frame(Genome_Name = gsub("_", " ", full_metadata$Display_Name), Accession = full_metadata$assembly_accession, Status = "Pending", Count = 0, stringsAsFactors = FALSE)
    for(i in 1:nrow(full_metadata)) {
      res <- extract_gene_manual(full_metadata$local_path[i], full_metadata$local_gff[i], full_metadata$Display_Name[i], target_gene, feature_type)
      if(!is.null(res)) { all_gene_seqs <- c(all_gene_seqs, res); status_df$Status[i] <- "Found"; status_df$Count[i] <- length(res) } else status_df$Status[i] <- "Missing"
    }
    write.csv(status_df, file.path(dir_logs, "extraction_status_report.csv"), row.names=FALSE)
    writeXStringSet(all_gene_seqs, file.path(dir_alignments, paste0("all_genomes_", target_gene, "_general.fasta")))
    
    if(length(all_gene_seqs) == 0) { message(paste("  -> [SKIP] No", target_gene, "extracted. Skipping.")); next }
    
    #  Dereplication - making sure indentical sequences from a single genome are removed
    message("--- Phase 5: Dereplication ---")
    seq_data <- data.frame(Name = names(all_gene_seqs), Sequence = as.character(all_gene_seqs), Accession = str_extract(names(all_gene_seqs), "GCF_[0-9]+\\.[0-9]+"), stringsAsFactors = FALSE) %>% left_join(full_metadata %>% select(assembly_accession, refseq_category), by = c("Accession" = "assembly_accession"))
    seq_to_keep <- list()
    for(i in 1:nrow(seq_data)) {
      row <- seq_data[i, ]
      if(grepl("N{5,}", row$Sequence)) { if(row$refseq_category %in% c("reference genome", "representative genome")) seq_to_keep[[row$Name]] <- gsub("N{5,}", "", row$Sequence) } else seq_to_keep[[row$Name]] <- row$Sequence
    }
    general_gene <- DNAStringSet(unlist(seq_to_keep))
    
    seq_df <- data.frame(Seq_Name = names(general_gene), Sequence = as.character(general_gene), stringsAsFactors = FALSE) %>%
      mutate(Accession = str_extract(Seq_Name, "GCF_[0-9]+\\.[0-9]+"), Copy_Num = as.numeric(str_extract(Seq_Name, "(?<=_copy)[0-9]+")), Seq_Hash = sapply(Sequence, function(x) digest::digest(x, algo="md5"))) %>%
      group_by(Accession, Seq_Hash) %>% arrange(Copy_Num) %>% dplyr::slice(1) %>% ungroup()
    
    general_gene <- DNAStringSet(seq_df$Sequence); names(general_gene) <- seq_df$Seq_Name
    
    # PCR
    message("--- Phase 6: PCR & Alignment ---")
    smart_align <- function(seqs, mafft_exec = "mafft") {
      if (length(seqs) < 2) return(seqs)
      use_mafft <- FALSE
      if (!is.null(mafft_exec) && mafft_exec != "") {
        mafft_check <- suppressWarnings(try(system2(mafft_exec, "--version", stdout=TRUE, stderr=TRUE), silent=TRUE))
        if (!inherits(mafft_check, "try-error") && (is.null(attr(mafft_check, "status")) || attr(mafft_check, "status") == 0)) use_mafft <- TRUE
      }
      if (use_mafft) {
        tmp_in <- tempfile(fileext=".fasta"); tmp_out <- tempfile(fileext=".fasta")
        writeXStringSet(seqs, tmp_in); system2(mafft_exec, args = c("--auto", "--quiet", tmp_in), stdout = tmp_out, stderr = FALSE)
        aligned_seqs <- readDNAStringSet(tmp_out); unlink(c(tmp_in, tmp_out)); return(aligned_seqs)
      } else { return(AlignSeqs(seqs, verbose = FALSE)) }
    }
    
  # Custom function to generate IUPAC-aware visual alignments
    get_mismatch_info <- function(primer_seq, target_seq) {
      aln <- Biostrings::pairwiseAlignment(primer_seq, target_seq, type="global")
      
      # Use alignedPattern and alignedSubject to safely extract sequences and bypass the namespace firewall
      p_aln_str <- as.character(Biostrings::alignedPattern(aln))
      t_aln_str <- as.character(Biostrings::alignedSubject(aln))
      
      p_aln <- strsplit(p_aln_str, "")[[1]]
      t_aln <- strsplit(t_aln_str, "")[[1]]
      
      # SAFEGUARD: If the alignment is completely empty due to a truncated genome edge
      if(length(p_aln) == 0) return(list(Target_Seq = "", Visual = "", Misses = nchar(primer_seq)))
      
      visual <- character(length(p_aln))
      mismatches <- 0
      
      # Use seq_along instead of 1:length to prevent the 1:0 looping bug
      for(k in seq_along(p_aln)) {
        p <- toupper(p_aln[k]); t <- toupper(t_aln[k])
        if (p == t && p != "-") {
          visual[k] <- "."
        } else if (p == "-" || t == "-") {
          visual[k] <- t
          mismatches <- mismatches + 1
        } else {
          allowed <- Biostrings::IUPAC_CODE_MAP[p]
          if (!is.na(allowed) && grepl(t, allowed)) {
            visual[k] <- "~"
          } else {
            visual[k] <- t
            mismatches <- mismatches + 1
          }
        }
      }
      return(list(Target_Seq = t_aln_str, Visual = paste(visual, collapse=""), Misses = mismatches))
    }

    target_seqs <- smart_align(general_gene, mafft_exec = mafft_path)
    writeXStringSet(target_seqs, file.path(dir_alignments, paste0("Alignment_Full_Extracted_", target_gene, ".fasta")))
    primer_mismatch_report <- data.frame() 
    
    for (region_name in names(primer_sets)) {
      p_info <- primer_sets[[region_name]]
      fwd_primer <- DNAString(p_info$Fwd); rev_primer <- DNAString(p_info$Rev)
      current_amplicons <- DNAStringSet()
      
      for(i in seq_along(target_seqs)) {
        subj <- RemoveGaps(target_seqs[i])[[1]]; subj_rc <- reverseComplement(subj)
        
        # Explicitly setting Lfixed/Rfixed=FALSE allows matchLRPatterns to read IUPAC correctly
        # 1. Look for any characters that are NOT standard A, C, G, or T
        fwd_degen <- grepl("[^ACGT]", as.character(fwd_primer), ignore.case = TRUE)
        rev_degen <- grepl("[^ACGT]", as.character(rev_primer), ignore.case = TRUE)
        
        # 2. If EITHER primer has degenerate bases, turn indels OFF. Otherwise, leave them ON.
        allow_indels <- !(fwd_degen || rev_degen)

        # 3. Pass the dynamic variable to the function
        run_pcr <- function(s) matchLRPatterns(fwd_primer, reverseComplement(rev_primer), s, 
                                               max.gaplength=as.numeric(p_info$Max), max.Lmismatch=5, max.Rmismatch=5, 
                                               with.Lindels=allow_indels, with.Rindels=allow_indels, Lfixed=FALSE, Rfixed=FALSE)
        
        v1 <- run_pcr(subj); v2 <- run_pcr(subj_rc); hit <- NULL; orientation <- "Fwd"
        if(length(v1)>0) { hit <- v1[which.max(width(v1))] } else if(length(v2)>0) { hit <- v2[which.max(width(v2))]; orientation <- "Rev" }
        
        if(!is.null(hit)) {
          amp <- as(hit, "DNAStringSet")[[1]]
          
          # Safely read min/max from CSV without crashing on NAs
          safe_min <- suppressWarnings(as.numeric(p_info$Min)[1]); if(is.na(safe_min)) safe_min <- 0
          safe_max <- suppressWarnings(as.numeric(p_info$Max)[1]); if(is.na(safe_max)) safe_max <- 10000
          
          if(isTRUE(length(amp) >= safe_min && length(amp) <= safe_max)) {
            current_amplicons <- c(current_amplicons, DNAStringSet(amp))
            names(current_amplicons)[length(current_amplicons)] <- paste0(names(target_seqs)[i], "_", region_name)
            
            tryCatch({
              # Stop doing coordinate math on the full genome! 
              # Just slice the first and last N bases directly off the exact extracted amplicon (amp)
              
              fwd_bind <- as.character(subseq(amp, 1, min(length(amp), length(fwd_primer))))
              rev_bind <- as.character(subseq(amp, max(1, length(amp) - length(rev_primer) + 1), length(amp)))
              
              # Pad with gaps if an indel caused it to be slightly shorter than the primer
              if(nchar(fwd_bind) < length(fwd_primer)) fwd_bind <- paste0(fwd_bind, strrep("-", length(fwd_primer) - nchar(fwd_bind)))
              if(nchar(rev_bind) < length(rev_primer)) rev_bind <- paste0(strrep("-", length(rev_primer) - nchar(rev_bind)), rev_bind)
              
              # Run our custom visual alignment function (using your fixed manual iteration version!)
              f_res <- get_mismatch_info(as.character(fwd_primer), fwd_bind)
              r_res <- get_mismatch_info(as.character(reverseComplement(rev_primer)), rev_bind)
              
              total_miss <- f_res$Misses + r_res$Misses
              
              # Optional filter: If you don't want 35-mismatch garbage in your CSV, uncomment the 'if' statement
              # if(total_miss <= 15) {
                  primer_mismatch_report <- rbind(primer_mismatch_report, data.frame(
                    Bacterium = gsub("_", " ", sub("_copy.*", "", names(target_seqs)[i])), 
                    Region = region_name, 
                    Primer_Fwd = as.character(fwd_primer),
                    Target_Fwd_Seq = f_res$Target_Seq,
                    Visual_Fwd_Match = f_res$Visual,
                    Mismatches_Fwd = f_res$Misses,
                    Primer_Rev = as.character(reverseComplement(rev_primer)),
                    Target_Rev_Seq = r_res$Target_Seq,
                    Visual_Rev_Match = r_res$Visual,
                    Mismatches_Rev = r_res$Misses,
                    Total_Mismatches = total_miss,
                    stringsAsFactors = FALSE
                  ))
              # }
            }, error=function(e) { message(paste("  -> [WARNING] Mismatch mapping failed for a genome:", e$message)) })
          
          } else {
            message(paste("  -> [PCR INFO] Amplicon length outside bounds for:", names(target_seqs)[i], "using primer set", region_name))
          }
        } else {
          message(paste("  -> [PCR INFO] No amplicons generated for:", names(target_seqs)[i], "using primer set", region_name))
        }
      }
      
      if(length(current_amplicons) >= 1) { 
        aln <- smart_align(current_amplicons, mafft_exec = mafft_path) 
        writeXStringSet(aln, file.path(dir_alignments, paste0("Alignment_", region_name, ".fasta"))) 
      }
    }
    
    if(nrow(primer_mismatch_report) > 0) {
      write.csv(primer_mismatch_report, file.path(dir_logs, "primer_mismatch_report_summary.csv"), row.names=FALSE)
    }
    
    # Tree Generation
    message("--- Phase 7: Trees ---")
    plot_tree <- function(aln_file, suffix, filter_ref=FALSE) {
      aln <- readDNAStringSet(aln_file)
      aln_accs <- str_extract(names(aln), "GCF_[0-9]+\\.[0-9]+")
      if(filter_ref) aln <- aln[aln_accs %in% full_metadata$assembly_accession[full_metadata$Is_VIP]]
      if(length(aln) < 3) return(NULL)
      tree <- njs(DistanceMatrix(aln, verbose=FALSE)); num_tips <- length(tree$tip.label)
      
      clean_tips <- gsub("_copy", " (copy", tree$tip.label)
      tree$tip.label <- paste0(gsub("_", " ", clean_tips), ")")
      
      options(ignore.negative.edge=TRUE)
      p <- ggtree(tree) + geom_tiplab(size=2) + theme_tree2() + labs(title=paste0(basename(aln_file), suffix)) + hexpand(0.5)
      ggsave(file.path(dir_trees, paste0("Tree_", basename(aln_file), suffix, ".pdf")), p, width=15, height=max(10, num_tips*0.2), limitsize=FALSE)
    }

    files <- list.files(dir_alignments, pattern="^Alignment_.*\\.fasta$", full.names=TRUE)
    for(f in files) { 
      if(!only_reference) plot_tree(f, "_All")
      plot_tree(f, "_RefOnly", TRUE) 
    }
    
    # Entropy Maps
    message("--- Phase 8: Entropy Maps ---")
    aln_file <- file.path(dir_alignments, paste0("Alignment_Full_Extracted_", target_gene, ".fasta"))
    if(file.exists(aln_file)) {
      master_aln <- readDNAStringSet(aln_file)
      generate_entropy_analysis <- function(aln_set, subset_name, title_suffix) {
        if(length(aln_set) == 0) return(NULL)
        # --- DYNAMIC BACKBONE SELECTION ---
        # Simply pick the longest available sequence to serve as the reference backbone
        best_name <- (data.frame(Name = names(aln_set), Length = width(RemoveGaps(aln_set)), stringsAsFactors = FALSE) %>% 
                        arrange(desc(Length)))$Name[1]
        backbone_aligned <- aln_set[best_name]; backbone_clean <- RemoveGaps(backbone_aligned)
        entropy_scores <- apply(as.matrix(aln_set), 2, function(col) { bases <- col[!col %in% c("-","N",".")]; if(length(bases)==0) return(0); freqs <- table(bases)/length(bases); -sum(freqs*log2(freqs)) })
        entropy_df <- data.frame(Position=1:length(entropy_scores), Entropy=entropy_scores)
        # 1. Gather the backbone and all its corresponding amplicons into a single list
        seqs_to_align <- DNAStringSet(backbone_clean)
        names(seqs_to_align) <- "Backbone"
        primer_names <- c()
        
        for(f in list.files(dir_alignments, pattern="^Alignment_.*\\.fasta$", full.names=TRUE)[basename(list.files(dir_alignments, pattern="^Alignment_.*\\.fasta$", full.names=TRUE)) != paste0("Alignment_Full_Extracted_", target_gene, ".fasta")]) {
          primer_set_name <- sub("Alignment_", "", sub("\\.fasta$", "", basename(f)))
          r_aln <- readDNAStringSet(f)
          
          rel_idx <- which(sub(paste0("_", primer_set_name, "$"), "", names(r_aln)) %in% names(aln_set))
          if(length(rel_idx) > 0) {
            backbone_amp_name <- paste0(best_name, "_", primer_set_name)
            amp_to_map <- if(backbone_amp_name %in% names(r_aln)) RemoveGaps(r_aln[backbone_amp_name])[[1]] else RemoveGaps(r_aln[rel_idx[1]])[[1]]
            
            new_seq <- DNAStringSet(amp_to_map)
            names(new_seq) <- primer_set_name
            seqs_to_align <- c(seqs_to_align, new_seq)
            primer_names <- c(primer_names, primer_set_name)
          }
        }
        
        amplicon_map <- data.frame()
        if (length(seqs_to_align) > 1) {
          # 2. Align them all together (just like the Test_alignment.fasta you made)
          mini_aln <- smart_align(seqs_to_align, mafft_exec = mafft_path)
          
          bb_aln_str <- as.character(mini_aln["Backbone"][[1]])
          bb_is_base <- strsplit(bb_aln_str, "")[[1]] != "-"
          
          get_aln_pos <- function(pos, aln_str) { is_base <- strsplit(as.character(aln_str), "")[[1]] != "-"; match_col <- which(cumsum(is_base) == pos); if(length(match_col) > 0) return(match_col[1]) else return(pos) }
          
          # 3. Read the exact coordinates straight out of the alignment
          for (p_name in primer_names) {
            amp_aln_str <- as.character(mini_aln[p_name][[1]])
            amp_is_base <- strsplit(amp_aln_str, "")[[1]] != "-"
            
            if (any(amp_is_base)) {
              # Find first and last physical base in the aligned amplicon
              start_idx <- which(amp_is_base)[1]
              end_idx <- rev(which(amp_is_base))[1]
              
              # Count how many absolute bases exist in the backbone up to those points
              start_clean <- sum(bb_is_base[1:start_idx])
              end_clean <- sum(bb_is_base[1:end_idx])
              
              # Map those absolute coordinates onto the final PDF Master Alignment gaps
              amplicon_map <- rbind(amplicon_map, data.frame(
                Primer_Set = p_name, 
                Start_Aln = get_aln_pos(max(1, start_clean), backbone_aligned[[1]]), 
                End_Aln = get_aln_pos(max(1, end_clean), backbone_aligned[[1]])
              ))
            }
          }
        }
        
        # Save entropy map
        p_basic <- ggplot() + geom_area(data=entropy_df, aes(x=Position, y=Entropy), fill="gray90", color="gray70", linewidth=0.2) + labs(title=paste("Entropy:", target_genus, title_suffix), subtitle=paste("Backbone:", gsub("_", " ", best_name)), y="Entropy") + theme_bw()
        if(nrow(amplicon_map) > 0) {
          amplicon_map <- amplicon_map %>% arrange(Start_Aln); amplicon_map$Track <- NA; tracks <- c()
          for(i in 1:nrow(amplicon_map)) { s <- amplicon_map$Start_Aln[i]; e <- amplicon_map$End_Aln[i]; assigned <- FALSE; if(length(tracks)>0) for(t in 1:length(tracks)) if(tracks[t]+50 < s) { amplicon_map$Track[i] <- t; tracks[t] <- e; assigned <- TRUE; break }; if(!assigned) { nt <- length(tracks)+1; amplicon_map$Track[i] <- nt; tracks[nt] <- e } }
          p_valid <- p_basic + geom_rect(data=mutate(amplicon_map, Y = -0.2 - ((Track-1) * 0.1)), aes(xmin=Start_Aln, xmax=End_Aln, ymin=Y, ymax=Y+0.05, fill="Primer"), color="black") + geom_text(data=mutate(amplicon_map, Y = -0.2 - ((Track-1) * 0.1)), aes(x=(Start_Aln+End_Aln)/2, y=Y+0.025, label=Primer_Set), size=2)
          pdf(file.path(dir_entropy, paste0("Entropy_Map", subset_name, ".pdf")), width=12, height=8); print(p_valid); dev.off()
        } else {
          # Save the basic map even if there are no primers
          pdf(file.path(dir_entropy, paste0("Entropy_Map", subset_name, ".pdf")), width=12, height=8); print(p_basic); dev.off()
        }
      }
      if(!only_reference) generate_entropy_analysis(master_aln, "_All", "(All)")
      
      ref_aln <- master_aln[str_extract(names(master_aln), "GCF_[0-9]+\\.[0-9]+") %in% full_metadata$assembly_accession[full_metadata$Is_VIP]]
      if(length(ref_aln)>=1) generate_entropy_analysis(ref_aln, "_RefOnly", "(RefOnly)")
    }
    
    # Gene/Primer Resolution & Coverage
    message("--- Phase 9: Resolution & Coverage ---")
    ref_meta <- full_metadata %>% filter(Is_VIP) %>% mutate(Clean_Species = gsub("\\[|\\]", "", stringr::word(organism_name, 1, 2)))
    TOTAL_REF_SPECIES <- n_distinct(ref_meta$Clean_Species); if(TOTAL_REF_SPECIES==0) TOTAL_REF_SPECIES <- 1
    scores <- data.frame()
    full_aln_path <- file.path(dir_alignments, paste0("Alignment_Full_Extracted_", target_gene, ".fasta"))
    full_seqs <- if(file.exists(full_aln_path)) readDNAStringSet(full_aln_path) else NULL
    
    # Initialize the Species-Level Tracking Dataframe safely
    if (nrow(ref_meta) > 0) {
      species_res_df <- data.frame(Genus = target_genus, Species = unique(ref_meta$Clean_Species), stringsAsFactors = FALSE)
    } else {
      species_res_df <- data.frame(Genus = character(0), Species = character(0), stringsAsFactors = FALSE)
    }
    
    for(f in files) {
      aln <- readDNAStringSet(f)
      primer_name <- sub("Alignment_", "", sub("\\.fasta$", "", basename(f)))
      
      # Setup default status for this primer as missing (will be overwritten if found)
      current_primer_status <- rep("Missing_or_Failed_PCR", nrow(species_res_df))
      names(current_primer_status) <- species_res_df$Species
      
      aln_accs <- stringr::str_extract(names(aln), "GCF_[0-9]+\\.[0-9]+")
      ref_accs <- aln_accs[aln_accs %in% ref_meta$assembly_accession]
      count_resolved <- 0
      
      if(length(ref_accs) > 1) {
        ref_aln <- aln[aln_accs %in% ref_accs]
        
        if(length(ref_aln) >= 3) {
          d <- DECIPHER::DistanceMatrix(ref_aln, verbose=FALSE, includeTerminalGaps=TRUE)
          
          d[is.na(d)] <- 1
          d[is.nan(d)] <- 1
          
          ref_tree <- ape::njs(d)
          
          tip_accs <- stringr::str_extract(ref_tree$tip.label, "GCF_[0-9]+\\.[0-9]+")
          tip_sp <- gsub("\\[|\\]", "", stringr::word(full_metadata$organism_name[match(tip_accs, full_metadata$assembly_accession)], 1, 2))
          
          for(sp in unique(tip_sp)) { 
            if(!is.na(sp)) { 
              my_tips <- ref_tree$tip.label[tip_sp == sp]
              my_idx <- which(rownames(d) %in% my_tips)
              others <- which(!(rownames(d) %in% my_tips))
              
              min_dist_outgroup <- if(length(others) > 0 && length(my_idx) > 0) min(d[my_idx, others]) else 1
              
              if(min_dist_outgroup > 0) {
                 if(length(my_tips) == 1) {
                    current_primer_status[sp] <- "Monophyletic"
                    count_resolved <- count_resolved + 1
                 } else {
                    if(ape::is.monophyletic(ref_tree, my_tips)) {
                       current_primer_status[sp] <- "Monophyletic"
                       count_resolved <- count_resolved + 1
                    } else {
                       current_primer_status[sp] <- "Polyphyletic"
                    }
                 }
              } else {
                 current_primer_status[sp] <- "Polyphyletic (Identical to outgroup)"
              }
            } 
          }
        } else if (length(ref_aln) == 2) {
           d <- DECIPHER::DistanceMatrix(ref_aln, verbose=FALSE, includeTerminalGaps=TRUE)
           
           d[is.na(d)] <- 1
           d[is.nan(d)] <- 1
           
           tip_accs <- stringr::str_extract(names(ref_aln), "GCF_[0-9]+\\.[0-9]+")
           tip_sp <- gsub("\\[|\\]", "", stringr::word(full_metadata$organism_name[match(tip_accs, full_metadata$assembly_accession)], 1, 2))
           
           if(d[1,2] > 0) {
              count_resolved <- 2 
              current_primer_status[tip_sp] <- "Monophyletic"
           } else {
              current_primer_status[tip_sp] <- "Polyphyletic (Identical to outgroup)"
           }
        }
      }
      
      # Bind the status column to the dataframe
      species_res_df[[primer_name]] <- current_primer_status[species_res_df$Species]
      
      missing_accs <- setdiff(ref_meta$assembly_accession, ref_accs)
      dropout_str <- ""
      
      if(length(missing_accs) > 0) {
        if(primer_name %in% names(primer_sets) && !is.null(full_seqs)) {
          p_info <- primer_sets[[primer_name]]
          fwd_primer <- DNAString(p_info$Fwd); rev_primer_rc <- reverseComplement(DNAString(p_info$Rev))
          f_fail <- 0; r_fail <- 0; b_fail <- 0
          for(m_acc in missing_accs) {
            seq_idx <- grep(m_acc, names(full_seqs))
            if(length(seq_idx) > 0) {
              subj <- RemoveGaps(full_seqs[seq_idx[1]])[[1]]; subj_rc <- reverseComplement(subj)
              f_pass <- length(matchPattern(fwd_primer, subj, max.mismatch=4, fixed=FALSE)) > 0 || length(matchPattern(fwd_primer, subj_rc, max.mismatch=4, fixed=FALSE)) > 0
              r_pass <- length(matchPattern(rev_primer_rc, subj, max.mismatch=4, fixed=FALSE)) > 0 || length(matchPattern(rev_primer_rc, subj_rc, max.mismatch=4, fixed=FALSE)) > 0
              if(!f_pass && !r_pass) b_fail <- b_fail + 1 else if(!f_pass) f_fail <- f_fail + 1 else if(!r_pass) r_fail <- r_fail + 1
            }
          }
          parts <- c()
          if(f_fail > 0) parts <- c(parts, paste(f_fail, "Fwd Primer"))
          if(r_fail > 0) parts <- c(parts, paste(r_fail, "Rev Primer"))
          if(b_fail > 0) parts <- c(parts, paste(b_fail, "Both Primers"))
          if(length(parts) > 0) dropout_str <- paste(parts, collapse=", ")
        } else if (primer_name == paste0("Full_Extracted_", target_gene)) {
          dropout_str <- paste(length(missing_accs), paste("Missing", target_gene, "Gene"))
        }
      }
      scores <- rbind(scores, data.frame(Primer_Set = primer_name, Resolved_Count = count_resolved, Total_Species = TOTAL_REF_SPECIES, Percent_Resolved = (count_resolved / TOTAL_REF_SPECIES) * 100, Dropout_Label = dropout_str))
    }
    
    # Save the local genus-level species report
    write.csv(species_res_df, file.path(dir_results, paste0("Species_Resolution_Report_", target_genus, ".csv")), row.names=FALSE)
    
    # Store it in the master list so it can be combined at the very end
    if(!exists("master_species_log")) master_species_log <- list()
    master_species_log[[target_genus]] <- species_res_df
    
    if(nrow(scores) > 0) {
      scores$Count_Label <- paste0(scores$Resolved_Count, "/", scores$Total_Species)
      scores$Plot_Label <- ifelse(scores$Dropout_Label == "", scores$Primer_Set, paste0(scores$Primer_Set, "\n(Failed: ", scores$Dropout_Label, ")"))
      scores$Bar_Color <- ifelse(grepl("Full_Extracted", scores$Primer_Set), "Baseline", "Primer")
      scores <- scores %>% arrange(Percent_Resolved)
      scores$Plot_Label <- factor(scores$Plot_Label, levels = unique(scores$Plot_Label))
      baseline_label <- scores$Plot_Label[scores$Primer_Set == paste0("Full_Extracted_", target_gene)]
      scores$Plot_Label <- forcats::fct_relevel(scores$Plot_Label, as.character(baseline_label), after = Inf)
      
      p_strict <- ggplot(scores, aes(x=reorder(Plot_Label, Percent_Resolved), y=Percent_Resolved, fill=Bar_Color)) +
        geom_bar(stat="identity", color="black", width=0.7) + geom_text(aes(label=Count_Label), hjust=-0.2, size=3.5, fontface="bold") + scale_fill_manual(values = c("Baseline" = "#e67e22", "Primer" = "#3498db")) + coord_flip() + scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25)) +
        labs(title=paste("Resolution (Reference Strains Only) -", target_genus), subtitle="Excludes conflicts caused by database noise/draft genomes", x="Primer Set", y="Resolution %") + theme_bw() + theme(legend.position = "none", axis.text.y = element_text(size=9))
      pdf(file.path(dir_results, "Strict_Resolution_Comparison_RefOnly.pdf"), width=12, height=10); print(p_strict); dev.off()
      
      genus_log <- data.frame(Genus = target_genus, stringsAsFactors = FALSE)
      for (i in 1:nrow(scores)) genus_log[[scores$Primer_Set[i]]] <- ifelse(scores$Percent_Resolved[i] == 100, "Yes", "No") 
      master_resolution_log[[target_genus]] <- genus_log
    }

    # ==================================================
    #  Higher Taxonomic Resolution (The Scout & Swarm)
    # ==================================================
    tax_levels_ordered <- c("Genus", "Family", "Order", "Class", "Phylum")
    target_idx <- which(tax_levels_ordered == "Genus")
    max_idx <- which(tax_levels_ordered == max_tax_level)
    
    if (length(max_idx) == 1 && max_idx > target_idx) {
      current_idx <- target_idx + 1
      
     while(current_idx <= max_idx) {
        current_tax_level <- tax_levels_ordered[current_idx]
        
        # Safely pull the name of the current Family/Order/Class from our metadata
        target_clade_name <- as.character(full_metadata[[current_tax_level]][1])
        
        # 1. First, check if the clade actually exists!
        if(is.na(target_clade_name)) {
          message(paste("  -> [WARNING] NCBI did not provide a", current_tax_level, "name for", target_genus, "- Aborting climb."))
          break
        }
        
        # 2. Only create the subfolders if we didn't abort
        lvl_align_dir <- file.path(dir_alignments, current_tax_level)
        lvl_tree_dir  <- file.path(dir_trees, current_tax_level)
        
        if(!dir.exists(lvl_align_dir)) dir.create(lvl_align_dir)
        if(!dir.exists(lvl_tree_dir)) dir.create(lvl_tree_dir)
        
        message(paste("\n--- Upgrading to", current_tax_level, "Level:", target_clade_name, "---"))
        
        # ------------------------------------------------
        #  Higher taxonomic look: DYNAMIC SCOUT SWEEP
        # ------------------------------------------------
        message(paste("  -> Initiating Scout Sweep for", current_tax_level, ":", target_clade_name))
        
        # 1. Ping NCBI for every Genus inside this specific Clade
        query <- paste0(target_clade_name, "[Subtree] AND genus[Rank]")
        search_res <- tryCatch(rentrez::entrez_search(db="taxonomy", term=query, retmax=5000), error=function(e) NULL)
        
        outgroup_genera <- c()
        if(!is.null(search_res) && length(search_res$ids) > 0) {
          tax_xml <- tryCatch(rentrez::entrez_fetch(db="taxonomy", id=search_res$ids, rettype="xml"), error=function(e) NULL)
          if(!is.null(tax_xml)) {
            xml_doc <- xml2::read_xml(tax_xml)
            outgroup_genera <- xml2::xml_text(xml2::xml_find_all(xml_doc, "//ScientificName"))
            outgroup_genera <- outgroup_genera[outgroup_genera != target_genus]
          }
        }
        
        if(length(outgroup_genera) == 0) {
          message(paste("  -> [WARNING] No outgroup genera found for", target_clade_name, "- Aborting climb."))
          break
        }
        
        # 2. Filter the Master NCBI table first to see who actually has valid genomes
        scout_meta <- master_meta %>%
          # Add Check_Name so we can cross-reference LPSN exactly
          mutate(
            Temp_Genus = gsub("\\[|\\]", "", stringr::word(organism_name, 1)),
            Check_Name = gsub("\\[|\\]", "", stringr::word(organism_name, 1, 2))
          ) %>%
          filter(Temp_Genus %in% outgroup_genera) %>%
          filter(refseq_category %in% c("reference genome", "representative genome")) %>%
          filter(!grepl(" sp\\.| sp$|Candidatus|uncultured", organism_name, ignore.case=TRUE)) %>%
          filter(contig_count <= max_contigs)
          
        # LPSN Filter for Higher Taxa
        if(enable_lpsn_check && !is.null(master_lpsn_clean)) {
           scout_meta <- scout_meta %>% filter(Check_Name %in% master_lpsn_clean$Full_Name)
        }
        
        # Now safely group and grab the single best representative
        scout_meta <- scout_meta %>%
          group_by(Temp_Genus) %>%
          arrange(desc(assembly_level == "Complete Genome"), seq_rel_date) %>%
          dplyr::slice(1) %>%  
          ungroup()
          
        # 3. Apply User-Defined Limit AFTER confirming they have genomes
        if(!is.infinite(max_scout_genera) && nrow(scout_meta) > max_scout_genera) {
          set.seed(42) # Ensures reproducibility between runs
          scout_meta <- scout_meta %>% dplyr::sample_n(max_scout_genera)
          message(paste("  -> [NOTE] Scout sweep capped at", max_scout_genera, "valid genera out of", length(outgroup_genera), "taxonomy nodes."))
        } else {
          message(paste("  -> [NOTE] Uncapped sweep: Targeting all", nrow(scout_meta), "valid genera found in clade."))
        }
        
        message(paste("  -> Found", nrow(scout_meta), "representative outgroup genomes for", target_clade_name))
        
        # HELPER FUNCTION: Download & Extract target_gene for Outgroups (CACHED)
        process_outgroups <- function(outgroup_meta, run_type_name) {
          if(nrow(outgroup_meta) == 0) return(NULL)
          
          # 1. Clean Names
          outgroup_meta <- outgroup_meta %>%
            mutate(
              Short_Name = stringr::word(organism_name, 1, 2),
              Display_Name = paste(Short_Name, assembly_accession, sep="_"),
              Display_Name = gsub("[^a-zA-Z0-9_\\.]", "_", Display_Name)
            )
          
          # 2. Point to the Global Database Cache (Instead of the local OUT_DIR)
          dir_cache <- file.path(db_dir, "Genome_Cache")
          if(!dir.exists(dir_cache)) dir.create(dir_cache)
          
          # 3. Download Genomes & GFFs (Skips if already cached!)
          outgroup_meta$local_path <- NA; outgroup_meta$local_gff <- NA
          for(i in 1:nrow(outgroup_meta)) {
            ftp_base <- outgroup_meta$ftp_path[i]
            folder_name <- basename(ftp_base)
            dest_genomic <- file.path(dir_cache, paste0(outgroup_meta$Display_Name[i], ".fna.gz"))
            dest_gff <- file.path(dir_cache, paste0(outgroup_meta$Display_Name[i], ".gff.gz"))
            
            if(!file.exists(dest_genomic)) tryCatch({ curl::curl_download(paste0(ftp_base, "/", folder_name, "_genomic.fna.gz"), dest_genomic, quiet=TRUE) }, error=function(e) {})
            if(!file.exists(dest_gff)) tryCatch({ curl::curl_download(paste0(ftp_base, "/", folder_name, "_genomic.gff.gz"), dest_gff, quiet=TRUE) }, error=function(e) {})
            
            if(file.exists(dest_genomic)) outgroup_meta$local_path[i] <- dest_genomic
            if(file.exists(dest_gff)) outgroup_meta$local_gff[i] <- dest_gff
          }
          
        # 4. Extract terget_gene
          outgroup_seqs <- DNAStringSet()
          for(i in 1:nrow(outgroup_meta)) {
            if(!is.na(outgroup_meta$local_path[i]) && !is.na(outgroup_meta$local_gff[i])) {
              res <- extract_gene_manual(outgroup_meta$local_path[i], outgroup_meta$local_gff[i], paste0("OUTGROUP_", outgroup_meta$Display_Name[i]), target_gene, feature_type)
              
              if(!is.null(res)) {
                # Dereplicate sequences prior to moving on with analysis
                unique_idx <- !duplicated(as.character(res))
                res <- res[unique_idx]

                
                outgroup_seqs <- c(outgroup_seqs, res)
              }
            }
          }
          
          # Note: We completely removed the `unlink()` deletion command so the cache persists!
          if(length(outgroup_seqs) == 0) return(NULL)
          return(outgroup_seqs)
        }

        # 3. Execute Scout Extraction
        scout_seqs <- process_outgroups(scout_meta, "Scout")
        
        if(is.null(scout_seqs)) {
          message(paste("  -> [WARNING] Failed to extract", target_gene, "from Scout outgroups. Aborting climb."))
          break
        }
        
      # 4. Build the Scout Tree
        combined_scout_seqs <- c(general_gene, scout_seqs)
        scout_aln <- smart_align(combined_scout_seqs, mafft_exec = mafft_path)
        
        # Calculate the distance safely
        scout_dist <- DistanceMatrix(scout_aln, verbose=FALSE)
        scout_dist[is.na(scout_dist)] <- 1
        scout_dist[is.nan(scout_dist)] <- 1
        
        # Generate the absolute distance matrix directly (bypassing the need for a tree!)
        dist_matrix <- as.matrix(scout_dist)
        target_tips <- names(combined_scout_seqs)[!grepl("OUTGROUP_", names(combined_scout_seqs))]
        outgroup_tips <- names(combined_scout_seqs)[grepl("OUTGROUP_", names(combined_scout_seqs))]

        message("  -> Saving Scout Alignment...")
        writeXStringSet(scout_aln, file.path(lvl_align_dir, paste0("Alignment_Scout_", current_tax_level, "_", target_genus, ".fasta")))
        
        # Only build and save the Scout Tree if we have 3 or more sequences
        if(length(combined_scout_seqs) >= 3) {
            scout_tree <- njs(scout_dist)
            pdf(file.path(lvl_tree_dir, paste0("Tree_Scout_", current_tax_level, "_", target_genus, ".pdf")), width=15, height=max(10, length(scout_tree$tip.label)*0.2))
            p_scout <- ggtree(scout_tree) + geom_tiplab(size=2) + theme_tree2() + 
              labs(title=paste("Scout Tree:", target_genus, "within", current_tax_level, target_clade_name)) + hexpand(0.5)
            print(p_scout)
            dev.off()
        }

        # ---------------------------
        # PHASE 2: THREAT DETECTION
        # ---------------------------
        message("  -> Calculating distances to identify nearest threats...")
        
        # 1. Generate the absolute distance matrix from the tree
        dist_matrix <- ape::cophenetic.phylo(scout_tree)
        
        # 2. Separate Target tips from Outgroup tips
        target_tips <- scout_tree$tip.label[!grepl("OUTGROUP_", scout_tree$tip.label)]
        outgroup_tips <- scout_tree$tip.label[grepl("OUTGROUP_", scout_tree$tip.label)]
        
        # SAFEGUARD: Ensure we have both targets and outgroups before mathing
        if(length(target_tips) == 0 || length(outgroup_tips) == 0) {
           message("  -> [WARNING] Missing tips in scout tree. Aborting climb.")
           break
        }
        
        # 3. Calculate the average distance from the Target clade to each specific Outgroup tip
        target_distances <- dist_matrix[target_tips, outgroup_tips, drop=FALSE]
        avg_distances <- colMeans(target_distances)
        
        # 4. Map the tips back to their Genera and find the closest threats
        threat_df <- data.frame(
          Tip_Name = names(avg_distances),
          Distance = avg_distances,
          stringsAsFactors = FALSE
        ) %>%
          # Extract the Genus name by splitting at the first underscore!
          mutate(Threat_Genus = stringr::word(gsub("OUTGROUP_", "", Tip_Name), 1, sep="_")) %>%
          group_by(Threat_Genus) %>%
          summarize(Mean_Distance = mean(Distance)) %>%
          arrange(Mean_Distance)
        
        # Save the decision log so users can verify why these threats were chosen
        write.csv(threat_df, file.path(dir_logs, paste0(current_tax_level, "_", target_genus, "_Cophenetic_Threats.csv")), row.names=FALSE)
        
        # 5. Lock in the targets for Phase 3 (The Swarm)
        # We use min() to prevent crashing if n_threats is larger than the available outgroups
        actual_threats_to_pull <- min(n_threats, nrow(threat_df))
        top_threats <- head(threat_df$Threat_Genus, actual_threats_to_pull)
        message(paste("  -> Top threats identified:", paste(top_threats, collapse=", ")))
        
        # DEBUG LINE
        message(paste("  -> DEBUG: Inspecting exact threat names to be fetched:", paste(paste0("'", top_threats, "'"), collapse=", ")))

    
        # -----------------------
        # THE SWARM (DEEP FETCH) 
        # -----------------------
        message(paste("  -> Initiating Swarm Fetch for top", actual_threats_to_pull, "threats..."))
        
        swarm_meta <- master_meta %>%
          mutate(
            Temp_Genus = gsub("\\[|\\]", "", stringr::word(organism_name, 1)),
            Temp_Species = stringr::word(organism_name, 1, 2),
            Check_Name = gsub("\\[|\\]", "", stringr::word(organism_name, 1, 2))
          ) %>%
          filter(Temp_Genus %in% top_threats) %>%
          filter(refseq_category %in% c("reference genome", "representative genome")) %>%
          filter(!grepl(" sp\\.| sp$|Candidatus|uncultured", organism_name, ignore.case=TRUE)) %>%
          filter(contig_count <= max_contigs)

        # LPSN Filter for Threat Swarm
        if(enable_lpsn_check && !is.null(master_lpsn_clean)) {
           swarm_meta <- swarm_meta %>% filter(Check_Name %in% master_lpsn_clean$Full_Name)
        }

        swarm_meta <- swarm_meta %>%
          group_by(Temp_Species) %>%
          arrange(desc(assembly_level == "Complete Genome"), seq_rel_date) %>%
          dplyr::slice(1) %>%
          ungroup()
          
        message(paste("  -> Downloading", nrow(swarm_meta), "threat species genomes..."))
        
        # Use our trusty helper function again!
        swarm_seqs <- process_outgroups(swarm_meta, "Swarm")
        
        if(is.null(swarm_seqs)) {
          message(paste("  -> [WARNING] Failed to extract", target_gene, "from Swarm outgroups. Aborting climb."))
          break
        }
        
        # Combine target genus + full threat swarm
        combined_swarm_seqs <- c(general_gene, swarm_seqs)
        swarm_aln <- smart_align(combined_swarm_seqs, mafft_exec = mafft_path)

        message("  -> Saving Swarm Alignment...")
        writeXStringSet(swarm_aln, file.path(lvl_align_dir, paste0("Alignment_Swarm_", current_tax_level, "_", target_genus, ".fasta")))
        
        swarm_dist <- DistanceMatrix(swarm_aln, verbose=FALSE)
        swarm_dist[is.na(swarm_dist)] <- 1
        swarm_dist[is.nan(swarm_dist)] <- 1

        target_tips_final <- names(combined_swarm_seqs)[!grepl("OUTGROUP_", names(combined_swarm_seqs))]
        outgroup_tips_final <- names(combined_swarm_seqs)[grepl("OUTGROUP_", names(combined_swarm_seqs))]

        if (length(combined_swarm_seqs) >= 3) {
            swarm_tree <- njs(swarm_dist)
            pdf(file.path(lvl_tree_dir, paste0("Tree_Swarm_", current_tax_level, "_", target_genus, ".pdf")), width=15, height=max(10, length(swarm_tree$tip.label)*0.2))
            p_swarm <- ggtree(swarm_tree) + geom_tiplab(size=2) + theme_tree2() + labs(title=paste("Swarm Tree:", target_genus, "vs Top Threats at", current_tax_level, "Level")) + hexpand(0.5)
            print(p_swarm)
            dev.off()
            
            is_exclusive_full <- ape::is.monophyletic(swarm_tree, target_tips_final)
        } else {
            # Fallback for 2-sequence Swarm (1 Target, 1 Threat)
            is_exclusive_full <- as.matrix(swarm_dist)[target_tips_final[1], outgroup_tips_final[1]] > 0
        }

        # ---------------------------------
        #  THE MULTI-PRIMER MONOPHYLY TEST
        # ---------------------------------
        message(paste("  -> Testing Clade Exclusivity (Monophyly) at", current_tax_level, "Level for all primers..."))
        
        # 1. Evaluate the Baseline (Full target_gene)
        target_tips_final <- swarm_tree$tip.label[!grepl("OUTGROUP_", swarm_tree$tip.label)]
        is_exclusive_full <- ape::is.monophyletic(swarm_tree, target_tips_final)
        
        level_log <- data.frame(
          Primer_Set = paste0("Full_Extracted_", target_gene),
          Status = ifelse(is_exclusive_full, "Resolved", "Polyphyletic (Failed)"),
          stringsAsFactors = FALSE
        )
        
        # 2. Loop through all requested Primer Sets
        for (region_name in names(primer_sets)) {
          p_info <- primer_sets[[region_name]]
          fwd_primer <- DNAString(p_info$Fwd); rev_primer <- DNAString(p_info$Rev)
          
          current_amplicons <- DNAStringSet()
          
          # Perform in-silico PCR on the entire swarm!
          for(i in seq_along(combined_swarm_seqs)) {
             subj <- RemoveGaps(combined_swarm_seqs[i])[[1]]; subj_rc <- reverseComplement(subj)
            
            # 1. Look for any characters that are NOT standard A, C, G, or T
        fwd_degen <- grepl("[^ACGT]", as.character(fwd_primer), ignore.case = TRUE)
        rev_degen <- grepl("[^ACGT]", as.character(rev_primer), ignore.case = TRUE)
        
        # 2. If EITHER primer has degenerate bases, turn indels OFF. Otherwise, leave them ON.
        allow_indels <- !(fwd_degen || rev_degen)

        # 3. Pass the dynamic variable to the function
        run_pcr <- function(s) matchLRPatterns(fwd_primer, reverseComplement(rev_primer), s, 
                                               max.gaplength=as.numeric(p_info$Max), max.Lmismatch=10, max.Rmismatch=10, 
                                               with.Lindels=allow_indels, with.Rindels=allow_indels, Lfixed=FALSE, Rfixed=FALSE)

             v1 <- run_pcr(subj); v2 <- run_pcr(subj_rc); hit <- NULL
             if(length(v1)>0) { hit <- v1[which.max(width(v1))] } else if(length(v2)>0) { hit <- v2[which.max(width(v2))] }

             if(!is.null(hit)) {
               amp <- as(hit, "DNAStringSet")[[1]]
               
               safe_min <- suppressWarnings(as.numeric(p_info$Min)[1]); if(is.na(safe_min)) safe_min <- 0
               safe_max <- suppressWarnings(as.numeric(p_info$Max)[1]); if(is.na(safe_max)) safe_max <- 10000
               
               if(isTRUE(length(amp) >= safe_min && length(amp) <= safe_max)) {
                 current_amplicons <- c(current_amplicons, DNAStringSet(amp))
                 names(current_amplicons)[length(current_amplicons)] <- names(combined_swarm_seqs)[i]
               }
               # If it falls outside the length bounds, it safely does nothing and moves on
             } else {
               # If no hit is found at all, it safely bypasses everything and moves to the next sequence
             }
          }
          
          # Check if the PCR successfully captured both targets and outgroups
          has_targets <- any(!grepl("OUTGROUP_", names(current_amplicons)))
          has_outgroups <- any(grepl("OUTGROUP_", names(current_amplicons)))
          
          if(length(current_amplicons) >= 2 && has_targets && has_outgroups) {
             # Re-align the amplicons and build the V-region tree
             amp_aln <- smart_align(current_amplicons, mafft_exec = mafft_path)
             
             amp_dist <- DistanceMatrix(amp_aln, verbose=FALSE)
             amp_dist[is.na(amp_dist)] <- 1
             amp_dist[is.nan(amp_dist)] <- 1

             writeXStringSet(amp_aln, file.path(lvl_align_dir, paste0("Alignment_Swarm_", current_tax_level, "_", region_name, "_", target_genus, ".fasta")))
             
             t_tips <- names(current_amplicons)[!grepl("OUTGROUP_", names(current_amplicons))]
             o_tips <- names(current_amplicons)[grepl("OUTGROUP_", names(current_amplicons))]

             if(length(current_amplicons) >= 3) {
                 amp_tree <- njs(amp_dist)
                 pdf(file.path(lvl_tree_dir, paste0("Tree_Swarm_", current_tax_level, "_", region_name, "_", target_genus, ".pdf")), width=15, height=max(10, length(amp_tree$tip.label)*0.2))
                 p_amp <- ggtree(amp_tree) + geom_tiplab(size=2) + theme_tree2() + labs(title=paste(region_name, "Swarm Tree:", target_genus, "at", current_tax_level, "Level")) + hexpand(0.5)
                 print(p_amp)
                 dev.off()
                 
                 is_exc <- ape::is.monophyletic(amp_tree, t_tips)
             } else {
                 is_exc <- as.matrix(amp_dist)[t_tips[1], o_tips[1]] > 0
             }

             if(length(t_tips) > 0) {
                level_log <- rbind(level_log, data.frame(Primer_Set = region_name, Status = ifelse(is_exc, "Resolved", "Polyphyletic (Failed)")))
             } else {
                level_log <- rbind(level_log, data.frame(Primer_Set = region_name, Status = "PCR_Failure (Target Genus Dropped)"))
             }
          } else {
             level_log <- rbind(level_log, data.frame(Primer_Set = region_name, Status = "PCR_Failure (Missing Comparison Groups)"))
          }
        }
        
        # 3. Format and export the Level-Specific Resolution Log
        level_log$Target_Genus <- target_genus
        level_log$Taxonomic_Level <- current_tax_level
        level_log$Clade_Name <- target_clade_name
        level_log$Threats_Tested <- paste(top_threats, collapse=" | ")
        
        level_log <- level_log %>% relocate(Target_Genus, Taxonomic_Level, Clade_Name, Primer_Set, Status, Threats_Tested)
        
        write.csv(level_log, file.path(dir_results, paste0("Monophyly_Log_", current_tax_level, "_", target_genus, ".csv")), row.names=FALSE)
        
        master_higher_taxa_log[[length(master_higher_taxa_log) + 1]] <- level_log

        # Advance the loop to the next taxonomic tier!
        current_idx <- current_idx + 1
      }
    }

    # Space Saver Cleanup
    if(!keep_genomes) {
      message("  -> Space Saver Enabled: Deleting downloaded .fasta and .gff files...")
      unlink(dir_genomes, recursive = TRUE)
      unlink(dir_annots, recursive = TRUE)
    }
    
    message(paste(">>> ANALYSIS COMPLETE FOR:", target_genus, "<<<"))
    
  } # <--- Bracket 1: This closes the main 'for (target_genus in target_genera)' LOOP
  
  # Final Log generation
  if (length(master_resolution_log) > 0) {
    write.csv(dplyr::bind_rows(master_resolution_log), file.path(output_dir, "Master_Resolution_Summary.csv"), row.names = FALSE)
    message(paste(">>> PIPELINE COMPLETE. Master summary saved to:", file.path(output_dir, "Master_Resolution_Summary.csv")))
  } else {
    message(">>> PIPELINE COMPLETE. No resolution data was generated.")
  }

  if (length(master_higher_taxa_log) > 0) {
    write.csv(dplyr::bind_rows(master_higher_taxa_log), file.path(output_dir, "Master_Higher_Taxa_Summary.csv"), row.names = FALSE)
    message(paste(">>> PIPELINE COMPLETE. Higher Taxa summary saved to:", file.path(output_dir, "Master_Higher_Taxa_Summary.csv")))
  }

  if (exists("master_species_log") && length(master_species_log) > 0) {
    write.csv(dplyr::bind_rows(master_species_log), file.path(output_dir, "Master_Species_Resolution_Summary.csv"), row.names = FALSE)
    message(paste(">>> PIPELINE COMPLETE. Master Species Resolution summary saved to:", file.path(output_dir, "Master_Species_Resolution_Summary.csv")))
  }

} # <--- Bracket 2: This is the FINAL bracket that closes 'run_marker_pipeline'
# Squamate Genomics Pipeline: Overview

This pipeline automates the collection, processing, and evolutionary analysis of reptile genome size and GC content data.

## 1. Data Acquisition: 01_scrape_genome_size.py
Crawls and extracts C-value records for reptiles from the Animal Genome Size Database.
* **Arguments**: 
    - `-o, --outdir`: Directory for CSVs and cached HTML (default: genomes/records/...).
    - `--pages`: Total results pages to iterate through (default: 5).
    - `--sleep`: Rate-limiting delay between GET requests (default: 0.75s).
* **Inputs**: Live web data from genomesize.com.
* **Outputs**: 
    - reptile_c_values_detected_summary.csv: Cleaned summary of detected C-values.
    - reptile_species_pages_raw_table_rows.csv: All raw table data for manual validation.
    - /html: Local cache of every scraped page.

## 2. Data Processing: 02_summarize_scraped_data.py
Transforms raw scraped strings into research-ready tidy data and standardized genomic units.
* **Arguments**: 
    - `-i, --input`: Path to the summary CSV from the scraper.
    - `-o, --outdir`: Output directory for analysis results.
* **Inputs**: reptile_c_values_detected_summary.csv.
* **Outputs**: 
    - species_summary.csv: Aggregated statistics (mean, median, etc.) per species.
    - expanded_c_values.csv: Every measurement converted to Mb and Gb (1pg = 978Mb).
    - overall_summary.csv: Global dataset metrics.

## 3. Data Integration: 03_make_master_input_tsvs.py
Aggregates phenotypic data (genome size and body mass) for specific project focal species.
* **Arguments**: 
    - genomes_dir: Top-level directory for the project.
    - `--manifest`: CSV defining focal species (overrides default metadata).
    - `--genome-size`: Path to species_summary.csv.
* **Inputs**: Species summary data and ecological datasets (e.g., Meiri or Title mass records).
* **Outputs**: 
    - genome_sizes.tsv: Standardized table of mean C-values and calculated sizes.
    - species_mass.tsv: Unified mass records with hierarchical priority (Meiri > Title).

## 4. Evolutionary Analysis: 04_calculate_gc_divergence.py
Computes GC3 content and divergence rates (dij and dianc) from phylogenomic models.
* **Arguments**: 
    - genomes_dir: Base project path.
    - `--number-of-species`: REQUIRED for node logic (first internal node = n + 1).
    - `--nhphyml-dir`: Path to .f_anc and .tree output files.
    - `--node-numbers`: Mapping file for nhPhyML node assignments.
* **Inputs**: nhPhyML output files, node mapping text, and the TSVs from Step 03.
* **Outputs**: 
    - gc3_summary.tsv: Mean and median GC3 values per species.
    - dij_results.tsv / dianc_results.tsv: Evolutionary divergence rates.
    - master_gc_regression_input.tsv: Final dataset for statistical modeling.



Genome Size Data Acquisition (Animal Genome Size Database)
File: 01_scrape_genome_size.py

	This script automates the retrieval of C-value data for Reptiles from the Animal Genome Size Database. It is designed to facilitate the collection of phenotypic data for evolutionary genomics research, specifically focusing on the relationship between genome size and thermal adaptation or genome structure.

	Core Functionality
	Targeted Search: Initializes a session to target the "Reptiles" category via the database's search parameters.

	Multi-Level Scraping:

	Iterates through result pagination to compile a master list of species-specific URLs.

	Performs deep-scraping of individual species pages to capture full HTML records for reproducibility.

	Intelligent Parsing: * Uses a Flexible C-value Regex to detect 1C values and "pg" units within unstructured page text.

	Extracts all structured table rows into a raw CSV format to ensure no metadata (e.g., method of estimation, references) is lost during the scrape.

	Bioinformatics Integration: Includes a normalization utility to convert inconsistent organism names into a standardized Genus_species key, compatible with standard bioinformatics pipelines like pyfaidx or SQL-based genomic databases.

	Outputs
	reptile_c_values_detected_summary.csv: A cleaned summary of inferred species names and detected C-values.

	reptile_species_pages_raw_table_rows.csv: Every row from every table on the species pages, preserved for manual validation of outliers.

	/html: A local cache of every scraped page to minimize repeated server hits and provide an offline record of the data source.

Genomic Metric Normalization & SummaryFile: 
02_summarize_scraped_data.py

	This script functions as the primary data processing engine for the C-value records retrieved from the Animal Genome Size Database. It transforms semi-structured text data into standardized genomic units necessary for evolutionary analysis.C ore FunctionalityData Tidying: Parses semicolon-delimited C-value strings and "explodes" them into a tidy, long-format structure where each row represents a single measurement. Unit Conversion: Automatically calculates genome size in Megabases (Mb) and Gigabases (Gb) using the industry-standard conversion: 1 pg = 978 Mb. This ensures compatibility with sequence-based assembly data (e.g., RefSeq GCF accessions) used in your genomics pipelines.Statistical Aggregation: Computes essential descriptive statistics—including mean, median, and standard deviation—at the species level to account for variability in legacy records.Automated Reporting: Generates a terminal-based report of the top 10 largest genomes in the dataset, providing an immediate overview of outliers within the squamate lineage.
	
Natural History Data Integration
File: 03_make_master_input_tsvs.py

	This script is the final staging step before executing GC-metric pipelines. It creates standardized TSV files that link genomic accessions (RefSeq/GenBank) with phenotypic and ecological traits.

	Core Functionality
	Species Set Definition: Identifies the focal species set for the project, defaulting to compleasm/records/metadata.csv unless a custom project manifest is supplied.

	Standardized Taxonomic Keys: Applies the Genus_species normalization protocol to all incoming datasets, ensuring that disparate ecological databases (Meiri, Title, and Animal Genome Size Database) can be merged without naming conflicts.

	Mass Priority Logic: Implements a hierarchical selection for body mass data, prioritizing Meiri metrics while falling back to Title data when necessary to maximize coverage across the squamate phylogeny.

	Automated Data Discovery: Intelligently scans records/natural_history/ for relevant CSV/TSV files based on column headers (e.g., mean_c_value_pg) and filename keywords, reducing manual configuration requirements.

	Outputs
	genome_sizes.tsv: A standardized table containing mean C-values and calculated sizes in Mb/Gb for all focal species.

	species_mass.tsv: A unified mass record including raw values from both sources and a "preferred" column for simplified downstream modeling.

GC Divergence & Evolutionary Rate CalculationFile: 
04_calculate_gc_divergence.py

	This script performs the core evolutionary analysis of the pipeline, extracting ancestral and terminal GC3 states from nhPhyML models to calculate GC divergence across the phylogeny.Core FunctionalityNode Heuristics: Automatically identifies the root and internal nodes using the required number-of-species parameter, ensuring accurate divergence calculations even when analyzing different taxonomic subsets.Metric Computation:GC3 Summaries: Aggregates GC content at the third codon position across all genes.$d_{ij}$ and $d_{ianc}$: Calculates the divergence between sister taxa and between terminal taxa and their most recent common ancestor.Multi-Dataset Integration: Seamlessly joins evolutionary metrics with the "Natural History" data (body mass and genome size) generated by 03_make_master_input_tsvs.py.Robust Filtering: By default, the script filters the final output to include only species successfully modeled in the nhPhyML trees, preventing blank rows in subsequent regression analyses.Outputsgc3_summary.tsv: Mean and median GC3 values per species.dij_results.tsv / dianc_results.tsv: Calculated divergence values for each terminal node.master_gc_regression_input.tsv: The final, comprehensive dataset used for statistical modeling of thermal adaptation and genomic evolution.

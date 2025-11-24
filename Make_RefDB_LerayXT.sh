#!/bin/bash

# Functions

# Function to print yellow text with a border
print_yellow() {
  echo -e "\033[1;33m____________________________________________________________________________\033[0m"
  echo ""
  echo -e "\033[1;33m$1\033[0m"
  echo -e "\033[1;33m____________________________________________________________________________\033[0m"
  sleep 0.3
}

# Function to validate a FASTA file
validate_fasta() {
  local fasta_file="$1"
  # Check if the file exists and is not empty
  if [[ ! -f "$fasta_file" || ! -s "$fasta_file" ]]; then
    echo -e "\033[1;31mError: File $fasta_file is missing or empty.\033[0m"
    return 1
  fi
  # Check if the file starts with a valid FASTA header (">")
  if ! head -n 1 "$fasta_file" | grep -q "^>"; then
    echo -e "\033[1;31mError: File $fasta_file does not start with a valid FASTA header.\033[0m"
    return 1
  fi
  # Count the number of header lines (lines starting with ">")
  header_count=$(grep -c "^>" "$fasta_file")
  # Count the number of sequence lines (lines not starting with ">")
  sequence_count=$(grep -vc "^>" "$fasta_file")
  # Check if the number of headers matches the number of sequences
  if [[ "$header_count" -ne "$sequence_count" ]]; then
    echo -e "\033[1;31mError: File $fasta_file has mismatched headers and sequences.\033[0m"
    echo -e "\033[1;31mHeaders: $header_count, Sequences: $sequence_count\033[0m"
    return 1
  fi
  # If all checks pass
  echo ""
  echo -e "\033[1;32mFile $fasta_file passed integrity checks.\033[0m"
  return 0
}

# Check if the user is mrjensen and activate the conda environment
if [ "$USER" == "mrjensen" ]; then
  print_yellow "Activating conda environment for user mrjensen"
  module load Miniconda3/22.11.1-1
  source ${EBROOTMINICONDA3}/bin/activate
  conda activate /cluster/projects/nn10069k/conda/crabs4
else
  print_yellow "WARNING: Ensure that the appropriate environment is activated"
fi

# Suppress Python warnings
export PYTHONWARNINGS="ignore"

# Get the current date and time
current_date=$(date +"%Y-%m-%d") # Format: YYYY-MM-DD
current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")  # Format: YYYY-MM-DD_HH-MM-SS

# Directory to store the output
output_dir="crabs_${current_date}" # Include the date and time in the directory name
mkdir -p "$output_dir"

# Check if required files exist
if [[ -f "$output_dir/names.dmp" && -f "$output_dir/nodes.dmp" && -f "$output_dir/nucl_gb.accession2taxid" ]]; then
  print_yellow "Required taxonomy files already exist in $output_dir. Skipping download."
else
  print_yellow "Required taxonomy files not found. Downloading them..."
  crabs \
  --download-taxonomy \
  --output "$output_dir"
fi

# Get the total number of taxa
total_taxa=$(wc -l < taxonlist.txt)

# Initialize a counter
current_taxon=0

# Check if one of the final bold-files already exists
if [[ ! -f "$output_dir/bold_COI5P_noblacklisted.fasta" && ! -f "$output_dir/bold_COI5P.fasta" ]]; then
    # Loop through each taxon in the taxa list
    while read -r taxon; do
        # Increment the counter
        current_taxon=$((current_taxon + 1))
        
        # Print progress
        current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
        print_yellow "[$current_datetime] Downloading data for taxon $current_taxon: $taxon, $total_taxa taxa remain"
        
        # Check if the taxon is one of the large file taxa
        if [[ "$taxon" == "Cecidomyiidae" || "$taxon" == "Chironomidae" || "$taxon" == "Phoridae" || "$taxon" == "Sciaridae" || "$taxon" == "Ceratopogonidae" || "$taxon" == "Arachnida" || "$taxon" == "Chordata" || "$taxon" == "Collembola" || "$taxon" == "Mollusca" || "$taxon" == "Malacostraca" ]]; then
            print_yellow "[$current_datetime] Note that $taxon may take a little while to download given its file size"
        fi
        
        # Run the crabs command
        output_file="${output_dir}/bold_${taxon}.fasta"
        crabs --download-bold --taxon "$taxon" --output "$output_file"
        tail -n 2 "$output_file"
        print_yellow "[$current_datetime] Does the downloaded $taxon fasta file above seem uncompromised?"
        
        # Check if the file was downloaded properly
        validate_fasta "$output_file"
        if [[ $? -ne 0 ]]; then
            echo -e "\033[1;31mError: File $output_file failed integrity checks. Consider splitting $taxon into smaller taxonomic entities for download. Exiting script.\033[0m"
            exit 1
        fi
        
        calc1=$(wc -l < "$output_file")
        calc2=$(grep '^>' "$output_file" | wc -l)
        calc3=$((calc1 / calc2))
        echo "Lines in file: $calc1"
        echo "Headers in file: $calc2"
        echo "Line/header ratio: $calc3"
        
        # Remove the processed taxon from the taxonlist
        sed -i "/^$taxon$/d" taxonlist.txt
        
        # Update the total taxa count
        total_taxa=$(wc -l < taxonlist.txt)
    done < taxonlist.txt

    # Merging the bold files
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Merging all bold-files and removing bold-taxa found in the file blacklist.txt"
    # Extracting patterns from blacklist.txt
    cat blacklist.txt | grep 'BOLD' | grep 'COI' | cut -f1 > blacklist_bold.txt
    
    # Making a single bold-file containing all records downloaded
    cat "$output_dir"/bold_*.fasta > "$output_dir"/bold_COI5P.fasta
    
    # Extracting everything but blacklisted taxa
    print_yellow "[$current_datetime] Extracting everything but blacklisted taxa"
    patterns=$(paste -sd '|' blacklist_bold.txt)
    awk -v patterns="$patterns" '
        NR == FNR {
            blacklist[$1]
            next
        } 
        {
            split(patterns, patternArray, "|")
            for (p in patternArray) {
                if ($0 ~ patternArray[p]) {
                    getline
                    next
                }
            }
        } 
        1
    ' blacklist_bold.txt $output_dir/bold_COI5P.fasta > $output_dir/bold_COI5P_noblacklisted.fasta
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Cleaning up ..."
    find "$output_dir" -type f -name "*.fasta" ! -name "*noblacklisted.fasta" -exec rm {} \;
    rm blacklist_bold.txt
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Cleanup complete!"
else
    print_yellow "Final bold-files $output_dir/bold_COI5P_noblacklisted.fasta and/or $output_dir/bold_COI5P.fasta already exist. Check whether you have $output_dir/bold_COI5P_noblacklisted.fasta, otherwise consider making it manually or downloading from scratch. Skipping taxon processing."
fi

# Check if the final NCBI-file already exists
if [[ ! -f "$output_dir/ncbi_COI_hits.fasta" ]]; then
  # Download nt database matching specific eukaryote query
  current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
  print_yellow "[$current_datetime] Downloading NCBI-files ..."
  crabs \
  --download-ncbi \
  --query '"Eukaryota"[Organism] AND species[All Fields] AND (CO1[GENE] OR COI[GENE] OR COX1[GENE] OR COXI[GENE]) AND 2003[PDAT] : 2025[PDAT] AND (mitochondrion[filter] AND ("250"[SLEN] : "25000"[SLEN]))' \
  --output "$output_dir"/ncbi_COI_hits.fasta \
  --email mads.jensen@uit.no \
  --database nucleotide
else
  print_yellow "File $output_dir/ncbi_COI_hits.fasta already exists. Skipping download of NCBI files."
fi

# Extracting everything but blacklisted taxa from NCBI file
if [[ ! -f "$output_dir/ncbi_COI_noblacklisted.fasta" ]]; then
    # Extracting patterns from blacklist.txt
    cat blacklist.txt | grep 'NCBI' | grep -E 'COI|mtDNA' | cut -f1 > blacklist_ncbi.txt
    # Extracting everything but blacklisted taxa
    print_yellow "[$current_datetime] Extracting all NCBI taxa but the blacklisted ones.."
    patterns=$(paste -sd '|' blacklist_ncbi.txt)	
    awk -v patterns="$patterns" '
        BEGIN {
            skip = 0
        }
        NR == FNR {
            blacklist[$1]
            next
        }
        {
            split(patterns, patternArray, "|")
            if (skip == 1) {
                if ($0 ~ /^>/) {
                    skip = 0
                } else {
                    next
                }
            }
            for (p in patternArray) {
                if ($0 ~ patternArray[p]) {
                    skip = 1
                    next
                }
            }
        }
        skip == 0
    ' blacklist_ncbi.txt "$output_dir/ncbi_COI_hits.fasta" > "$output_dir/ncbi_COI_noblacklisted.fasta"
else
    print_yellow "File $output_dir/ncbi_COI_noblacklisted.fasta already exists. Skipping download of NCBI files."
fi

# Check if the output file already exists
if [ ! -f "$output_dir/crabs_bold_COI5P.txt" ]; then
    # Import downloaded data to CRABS format (bold)
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Importing bold-files to CRABS format..."
    crabs \
    --import \
    --import-format bold \
    --input "$output_dir"/bold_COI5P_noblacklisted.fasta \
    --names "$output_dir"/names.dmp \
    --nodes "$output_dir"/nodes.dmp \
    --acc2tax "$output_dir"/nucl_gb.accession2taxid \
    --output "$output_dir"/crabs_bold_COI5P.txt \
    --ranks 'domain;phylum;class;order;family;genus;species'
else
    print_yellow "File $output_dir/crabs_bold_COI5P.txt already exists. Skipping import."
fi

# Check if the output file already exists
if [ ! -f "$output_dir/crabs_ncbi.txt" ]; then
    # Import downloaded data to CRABS format (ncbi)
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Importing NCBI-files to CRABS format..."
    crabs \
    --import \
    --import-format ncbi \
    --input "$output_dir"/ncbi_COI_noblacklisted.fasta \
    --names "$output_dir"/names.dmp \
    --nodes "$output_dir"/nodes.dmp \
    --acc2tax "$output_dir"/nucl_gb.accession2taxid \
    --output "$output_dir"/crabs_ncbi.txt \
    --ranks 'domain;phylum;class;order;family;genus;species'
else
    print_yellow "File $output_dir/crabs_ncbi.txt already exists. Skipping import."
fi

# Merging (required when sequencing data from multiple repositories was used)
# Check if the merged output file already exists
if [ ! -f "$output_dir/merged.txt" ]; then
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Merging bold and NCBI files..."
    crabs \
    --merge \
    --input "$output_dir/crabs_ncbi.txt;$output_dir/crabs_bold_COI5P.txt" \
    --uniq \
    --output "$output_dir/merged.txt"
else
    print_yellow "File $output_dir/merged.txt already exists. Skipping merge."
fi

# Extract amplicon regions through in-silico PCR analysis
# Check if the in-silico PCR output file already exists
if [ ! -f "$output_dir/insilico.txt" ]; then
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Running relaxed in-silico PCR ..."
    print_yellow "[$current_datetime] Looking for Leray-XT priming sites ..."
    print_yellow "[$current_datetime] This step will remove non-COI references from the database ..."
    crabs \
    --in-silico-pcr \
    --input "$output_dir"/merged.txt \
    --output "$output_dir"/insilico.txt \
    --forward GGWACWRGWTGRACWITITAYCCYCC \
    --reverse TAIACYTCIGGRTGICCRAARAAYCA \
    --relaxed
else
    print_yellow "File $output_dir/insilico.txt already exists. Skipping in-silico PCR."
fi

# Filtering the local reference database
# Check if the filtered output file already exists
if [ ! -f "$output_dir/FilteredRefs_CRABS.txt" ]; then
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Filtering based on the following criteria:"
    print_yellow "[$current_datetime] Minimum length: 283"
    print_yellow "[$current_datetime] Maximum length: 343"
    print_yellow "[$current_datetime] Maximum Ns in sequence: 1"
    print_yellow "[$current_datetime] Remove environmental sequences"
    print_yellow "[$current_datetime] Maximum two out of seven NAs in taxonomic ranks"
	print_yellow "[$current_datetime] Removing a redundant entry of NCBI accession number EU148067"
    print_yellow "[$current_datetime] Removal of records without species ID is currently deactivated"
    crabs \
    --filter \
    --input "$output_dir"/insilico.txt \
    --output "$output_dir"/FilteredRefs_CRABS.txt \
    --minimum-length 283 \
    --maximum-length 343 \
    --maximum-n 1 \
    --environmental \
    --rank-na 3
	linetoremove=$(grep -n 'EU148067' "$output_dir"/FilteredRefs_CRABS.txt | cut -d ':' -f1 | tail -n 1)
	sed "${linetoremove}d" "$output_dir"/FilteredRefs_CRABS.txt > "$output_dir"/FilteredRefs2_CRABS.txt
	mv "$output_dir"/FilteredRefs2_CRABS.txt "$output_dir"/FilteredRefs_CRABS.txt
else
    print_yellow "File $output_dir/FilteredRefs_CRABS.txt already exists. Skipping filtering."
fi

# Exporting the database for use
# Check if the BLAST_TAX_COI output file already exists
if [ ! -f "$output_dir/BLAST_TAX_COI.nsq" ]; then
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Exporting filtered database to blast-tax format ..."
    crabs \
    --export \
    --input "$output_dir"/FilteredRefs_CRABS.txt \
    --output "$output_dir"/BLAST_TAX_COI \
    --export-format 'blast-tax'
else
    print_yellow "File $output_dir/BLAST_TAX_COI already exists. Skipping database export."
fi

# Cleaning up
current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
print_yellow "[$current_datetime] Clean-up of the Reference database folder is currently deactivated"
# Uncomment the following lines to enable cleanup
# rm "$output_dir"/insilico.txt "$output_dir"/merged.txt "$output_dir"/bold_COI5P.fasta "$output_dir"/nucl_gb.accession2taxid "$output_dir"/names.dmp "$output_dir"/nodes.dmp "$output_dir"/taxdb.btd "$output_dir"/taxdb.bti "$output_dir"/taxonomy4blast.sqlite3
# current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
# print_yellow "[$current_datetime] Cleanup done"
print_yellow "[$current_datetime] Enjoy your new reference database!"

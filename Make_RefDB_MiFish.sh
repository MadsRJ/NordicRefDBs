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

# Function to print green text with a border
print_green() {
  echo -e "\033[1;32m____________________________________________________________________________\033[0m"
  echo ""
  echo -e "\033[1;32m$1\033[0m"
  echo -e "\033[1;32m____________________________________________________________________________\033[0m"
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
output_dir="crabs_MiFish_${current_date}" # Include the date and time in the directory name
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

# Check if the final NCBI-file already exists
if [[ ! -f "$output_dir/ncbi_12S_hits.fasta" ]]; then
  # Download nt database matching specific eukaryote query
  current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
  print_yellow "[$current_datetime] Downloading NCBI-files ..."
  crabs \
  --download-ncbi \
  --query '("Vertebrata"[Organism] OR "Vertebrata"[All Fields]) AND (is_nuccore[filter] AND mitochondrion[filter] AND ("140"[SLEN] : "25000"[SLEN]))' \
  --output "$output_dir"/ncbi_12S_hits.fasta \
  --email mads.jensen@uit.no \
  --database nucleotide
else
  print_yellow "File $output_dir/ncbi_12S_hits.fasta already exists. Skipping download of NCBI files."
fi

# Check if the output file already exists
if [ ! -f "$output_dir/crabs_ncbi.txt" ]; then
    # Import downloaded data to CRABS format (ncbi)
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Importing NCBI-files to CRABS format..."
    crabs \
    --import \
    --import-format ncbi \
    --input "$output_dir"/ncbi_12S_hits.fasta \
    --names "$output_dir"/names.dmp \
    --nodes "$output_dir"/nodes.dmp \
    --acc2tax "$output_dir"/nucl_gb.accession2taxid \
    --output "$output_dir"/crabs_ncbi.txt \
    --ranks 'domain;phylum;class;order;family;genus;species'
else
    print_yellow "File $output_dir/crabs_ncbi.txt already exists. Skipping import."
fi

# Extract amplicon regions through in-silico PCR analysis
# Check if the in-silico PCR output file already exists
if [ ! -f "$output_dir/insilico.txt" ]; then
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Running relaxed in-silico PCR ..."
    print_yellow "[$current_datetime] Looking for MiFish (UiT-version) priming sites ..."
    print_yellow "[$current_datetime] This step will remove non-12S references from the database ..."
    crabs \
    --in-silico-pcr \
    --input "$output_dir"/crabs_ncbi.txt \
    --output "$output_dir"/insilico.txt \
    --forward GCCGGTAAAACTCGTGCCAGC \
    --reverse CATAGTGGGGTATCTAATCCCAGTTTG \
    --relaxed
else
    print_yellow "File $output_dir/insilico.txt already exists. Skipping in-silico PCR."
fi

# Retrieve amplicons without primer-binding regions
# Check if the in-silico PCR output file already exists
if [ ! -f "$output_dir/aligned.txt" ]; then
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Retrieving amplicons without primer-binding regions ..."
    print_yellow "[$current_datetime] Looking for MiFish (UiT-version) region ..."
    print_yellow "[$current_datetime] This step will fetch sequences that do not contain the priming sites ..."
    crabs \
    --pairwise-global-alignment \
    --input "$output_dir"/crabs_ncbi.txt \
	--amplicons "$output_dir"/insilico.txt \
    --output "$output_dir"/aligned.txt \
    --forward GCCGGTAAAACTCGTGCCAGC \
    --reverse CATAGTGGGGTATCTAATCCCAGTTTG \
	--size-select 10000 \
	--percent-identity 0.92 \
	--coverage 95
else
    print_yellow "File $output_dir/aligned.txt already exists. Skipping amplicon retrieval step"
fi

# Removing blacklisted accession numbers, adding missing amplicons and fixing synonyms and subspecies of Nordic species.
# Check if the in-silico PCR output file already exists
if [ ! -f "$output_dir/aligned_curated.txt" ]; then
    # Extracting patterns from blacklist.txt
    cat blacklist.txt | grep 'NCBI' | grep -E '12S|mtDNA' | cut -f1 > blacklist_ncbi.txt
    
    # Extracting everything but blacklisted taxa
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Extracting all accessions but the blacklisted ones.."
    patterns=$(paste -sd '|' blacklist_ncbi.txt)
    grep -v -f <(echo $patterns | tr '|' '\n') $output_dir/aligned.txt > $output_dir/aligned_tmp1.txt
    
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Adding some accessions missed by the CRABS amplicon retrieval step"
    cat $output_dir/aligned_tmp1.txt ToBeAdded_CRABS.txt > $output_dir/aligned_tmp2.txt
    
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Fixing some synonymy issues of fishes known to occur in Nordic countries"
    
    awk -F '\t' -v OFS='\t' '
    $2 == "Trigloporus lastoviza" {
        $2 = "Chelidonichthys lastoviza"; 
        $9 = "Chelidonichthys"; 
        $10 = "Chelidonichthys lastoviza"
    } 
    $2 == "Nannobrachium atrum" {
        $2 = "Lampanyctus ater"; 
        $9 = "Lampanyctus"; 
        $10 = "Lampanyctus ater"
    } 
    $2 == "Ulcina olrikii" {
        $2 = "Aspidophoroides olrikii"; 
        $9 = "Aspidophoroides"; 
        $10 = "Aspidophoroides olrikii"
    } 
    $2 == "Gonostoma elongatum" {
        $2 = "Sigmops elongatus"; 
        $9 = "Sigmops"; 
        $10 = "Sigmops elongatus"
    } 
    $2 == "Acipenser baerii" {
        $2 = "Huso baerii"; 
        $9 = "Huso"; 
        $10 = "Huso baerii"
    } 
    $2 == "Acipenser gueldenstaedtii" {
        $2 = "Huso gueldenstaedtii"; 
        $9 = "Huso"; 
        $10 = "Huso gueldenstaedtii"
    } 
    $2 == "Acipenser naccarii" {
        $2 = "Huso naccarii"; 
        $9 = "Huso"; 
        $10 = "Huso naccarii"
    } 
    $2 == "Acipenser ruthenus" {
        $2 = "Huso ruthenus"; 
        $9 = "Huso"; 
        $10 = "Acipenser ruthenus"
    } 
    $2 == "Acipenser stellatus" {
        $2 = "Huso stellatus"; 
        $9 = "Huso"; 
        $10 = "Huso stellatus"
    } 
    $2 == "Acipenser sturio" {
        $2 = "Huso sturio"; 
        $9 = "Huso"; 
        $10 = "Huso sturio"
    } 
    $2 == "Acipenser transmontanus" {
        $2 = "Sinosturio transmontanus"; 
        $9 = "Sinosturio"; 
        $10 = "Sinosturio transmontanus"
    } 
    $2 == "Abramis bjoerkna" {
        $2 = "Blicca bjoerkna"; 
        $9 = "Blicca"; 
        $10 = "Blicca bjoerkna"
    } 
    $2 == "Theragra finnmarchica" {
        $2 = "Gadus chalcogrammus"; 
        $3 = "1042646"; 
        $9 = "Gadus"; 
        $10 = "Gadus chalcogrammus"
    } 
    $2 == "Gadus finnmarchicus" {
        $2 = "Gadus chalcogrammus"; 
        $3 = "1042646"; 
        $9 = "Gadus"; 
        $10 = "Gadus chalcogrammus"
    } 
    $2 == "Gadus ogac" {
        $2 = "Gadus macrocephalus"; 
        $3 = "80720"; 
        $9 = "Gadus"; 
        $10 = "Gadus macrocephalus"
    } 
    $2 == "Lepidion eques" {
        $2 = "Lepidion lepidion"; 
        $3 = "1229969"; 
        $9 = "Lepidion"; 
        $10 = "Lepidion lepidion"
    } 
    $2 == "Gobiusculus flavescens" {
        $2 = "Pomatoschistus flavescens"; 
        $9 = "Pomatoschistus"; 
        $10 = "Pomatoschistus flavescens"
    } 
    $2 == "Liza aurata" {
        $2 = "Chelon auratus"; 
        $9 = "Chelon"; 
        $10 = "Chelon auratus"
    } 
    $2 == "Liza ramado" {
        $2 = "Chelon ramada"; 
        $9 = "Chelon"; 
        $10 = "Chelon ramada"
    } 
    $2 == "Cottus poecilopus" {
        $2 = "Alpinocottus poecilopus"; 
        $9 = "Alpinocottus"; 
        $10 = "Alpinocottus poecilopus"
    } 
    $2 == "Lycodes rossi" {
        $2 = "Lycodes reticulatus"; 
        $3 = "215418"; 
        $9 = "Lycodes"; 
        $10 = "Lycodes reticulatus"
    } 
    $2 == "Zeugopterus norvegicus" {
        $2 = "Phrynorhombus norvegicus"; 
        $9 = "Phrynorhombus"; 
        $10 = "Phrynorhombus norvegicus"
    } 
    $2 == "Psetta maxima" {
        $2 = "Scophthalmus maximus"; 
        $9 = "Scophthalmus"; 
        $10 = "Scophthalmus maximus"
    } 
    $2 == "Dipturus linteus" {
        $2 = "Rajella lintea"; 
        $9 = "Rajella"; 
        $10 = "Rajella lintea"
    } 
    $2 == "Salvelinus umbla" {
        $2 = "Salvelinus alpinus"; 
        $3 = "8036"; 
        $9 = "Salvelinus"; 
        $10 = "Salvelinus alpinus"
    } 
    $2 == "Centroscymnus crepidater" {
        $2 = "Centroselachus crepidater"; 
        $9 = "Centroselachus"; 
        $10 = "Centroselachus crepidater"
    } 
    $2 == "Cromileptes altivelis" {
        $2 = "Chromileptes altivelis"; 
        $9 = "Chromileptes"; 
        $10 = "Chromileptes altivelis"
    } 
    1' $output_dir/aligned_tmp2.txt > $output_dir/aligned_tmp3.txt

    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Changing subspecies taxIDs to species taxIDs for Nordic species with multiple taxIDs present"
    
    awk -F '\t' -v OFS='\t' '
    $2 == "Abramis brama" {
        $3 = "38527" 
    } 
    $2 == "Auxis rochei" {
        $3 = "217026" 
    } 
    $2 == "Dipturus batis" {
        $3 = "420460" 
    } 
    $2 == "Gadiculus argenteus" {
        $3 = "185737" 
    } 
    $2 == "Huso baerii" {
        $3 = "27689" 
    } 
    $2 == "Huso stellatus" {
        $3 = "7903" 
    } 
    $2 == "Phoxinus phoxinus" {
        $3 = "58324" 
    } 
    $2 == "Thunnus thynnus" {
        $3 = "8237" 
    } 
    $2 == "Trichiurus lepturus" {
        $3 = "13733" 
    } 
    $2 == "Acipenser oxyrinchus" {
        $2 = "Acipenser oxyrinchus"; 
        $3 = "36177"; 
        $10 = "Acipenser oxyrinchus"
    } 
    $2 == "Salvelinus alpinus" {
        $3 = "8036" 
    } 
    $2 == "Gasterosteus aculeatus" {
        $3 = "69293" 
    } 
    $2 == "Salmo trutta" {
        $2 = "Salmo trutta"; 
        $3 = "8032"; 
        $10 = "Salmo trutta"
    } 
    $2 == "Oncorhynchus clarkii" {
        $3 = "30962" 
    } 
    $2 == "Oncorhynchus mykiss" {
        $2 = "Oncorhynchus mykiss"; 
        $3 = "8022"; 
        $10 = "Oncorhynchus mykiss"
    } 
    $2 == "Cyprinus carpio" {
        $2 = "Cyprinus carpio"; 
        $3 = "7962"; 
        $10 = "Cyprinus carpio"
    } 
    $2 == "Carassius langsdorfii" {
        $2 = "Carassius langsdorfii"; 
        $3 = "138676"; 
        $10 = "Carassius langsdorfii"
    } 
    $2 == "Carassius auratus" {
        $2 = "Carassius auratus"; 
        $3 = "7957"; 
        $10 = "Carassius auratus"
    } 
    1' $output_dir/aligned_tmp3.txt > $output_dir/aligned_tmp4.txt

    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Changing genera to match Eschmeyer families for all species"
    
    awk -F '\t' -v OFS='\t' '
    $9 == "Anotopterus" {
        $8 = "Anotopteridae" 
    } 
    $9 == "Omosudis" {
        $8 = "Omosudidae" 
    } 
    $9 == "Carcharodon" {
        $8 = "Lamnidae" 
    } 
    $9 == "Cetorhinus" {
        $8 = "Cetorhinidae" 
    } 
    $9 == "Isurus" {
        $8 = "Lamnidae" 
    } 
    $9 == "Lamna" {
        $8 = "Lamnidae" 
    } 
    $9 == "Anoplogaster" {
        $8 = "Anoplogastridae" 
    } 
    $9 == "Atherion" {
        $8 = "Atherionidae" 
    } 
    $9 == "Cololabis" {
        $8 = "Scomberesocidae" 
    } 
    $9 == "Scomberesox" {
        $8 = "Scomberesocidae" 
    } 
    $9 == "Grammatobothus" {
        $8 = "Grammatobothidae" 
    } 
    $9 == "Monolene" {
        $8 = "Monolenidae" 
    } 
    $9 == "Taeniopsetta" {
        $8 = "Taeniopsettidae" 
    } 
    $9 == "Ogilbia" {
        $8 = "Dinematichthyidae" 
    } 
    $9 == "Typhlias" {
        $8 = "Dinematichthyidae" 
    } 
    $9 == "Antigonia" {
        $8 = "Antigoniidae" 
    } 
    $9 == "Galeocerdo" {
        $8 = "Galeocerdonidae" 
    } 
    $9 == "Sphyrna" {
        $8 = "Sphyrnidae" 
    } 
    $9 == "Alosa" {
        $8 = "Alosidae" 
    } 
    $9 == "Amblygaster" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Anodontostoma" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Brevoortia" {
        $8 = "Alosidae" 
    } 
    $9 == "Clupanodon" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Clupeichthys" {
        $8 = "Ehiravidae" 
    } 
    $9 == "Clupeoides" {
        $8 = "Ehiravidae" 
    } 
    $9 == "Clupeonella" {
        $8 = "Ehiravidae" 
    } 
    $9 == "Dorosoma" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Ehirava" {
        $8 = "Ehiravidae" 
    } 
    $9 == "Escualosa" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Ethmalosa" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Gilchristella" {
        $8 = "Ehiravidae" 
    } 
    $9 == "Gudusia" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Harengula" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Herklotsichthys" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Hilsa" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Jenkinsia" {
        $8 = "Spratelloididae" 
    } 
    $9 == "Konosirus" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Limnothrissa" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Microthrissa" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Nematalosa" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Odaxothrissa" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Opisthonema" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Pellonula" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Potamothrissa" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Sardina" {
        $8 = "Alosidae" 
    } 
    $9 == "Sardinella" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Sardinops" {
        $8 = "Alosidae" 
    } 
    $9 == "Spratelloides" {
        $8 = "Spratelloididae" 
    } 
    $9 == "Stolothrissa" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Sundasalanx" {
        $8 = "Ehiravidae" 
    } 
    $9 == "Tenualosa" {
        $8 = "Dorosomatidae" 
    } 
    $9 == "Aperioptus" {
        $8 = "Gonorynchidae" 
    } 
    $9 == "Alcichthys" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Argyrocottus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Artediellus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Artedius" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Ascelichthys" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Astrocottus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Bero" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Chitonotus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Clinocottus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Cottiusculus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Enophrys" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Furcina" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Gymnocanthus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Hemilepidotus" {
        $8 = "Agonidae" 
    } 
    $9 == "Icelinus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Icelus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Megalocottus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Micrenophrys" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Microcottus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Myoxocephalus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Ocynectes" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Oligocottus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Orthonopias" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Phasmatocottus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Porocottus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Pseudoblennius" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Radulinopsis" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Radulinus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Rastrinus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Ricuzenius" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Scorpaenichthys" {
        $8 = "Jordaniidae" 
    } 
    $9 == "Stlengis" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Taurocottus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Taurulus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Trichocottus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Triglops" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Vellitor" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Zesticelus" {
        $8 = "Psychrolutidae" 
    } 
    $9 == "Neocyema" {
        $8 = "Neocyematidae" 
    } 
    $9 == "Raniceps" {
        $8 = "Ranicipitidae" 
    } 
    $9 == "Acanthogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Amblychaeturichthys" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Amblyotrypauchen" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Apocryptodon" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Astrabe" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Awaous" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Boleophthalmus" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Brachygobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Buenia" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Caragobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Chaenogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Chaeturichthys" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Clariger" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Crystallogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Ctenogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Ctenotrypauchen" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Eucyclogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Eugnathogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Eutaeniichthys" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Evorthodus" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Gillichthys" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Gnatholepis" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Gobioides" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Gobionellus" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Gobiopterus" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Gymnogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Hemigobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Inu" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Knipowitschia" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Lentipes" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Lepidogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Lethops" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Leucopsarion" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Lophiogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Luciogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Mugilogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Odontamblyopus" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Oligolepis" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Oxuderces" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Oxyurichthys" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Pandaka" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Papuligobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Parapocryptes" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Paratrypauchen" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Periophthalmodon" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Periophthalmus" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Pomatoschistus" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Pseudaphya" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Pseudogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Pterogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Redigobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Rhinogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Sagamia" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Scartelaos" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Schismatogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Sicydium" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Sicyopterus" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Sicyopus" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Siphonogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Smilosicyopus" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Stenogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Stigmatogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Stiphodon" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Suruga" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Taenioides" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Tasmanogobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Tridentiger" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Trypauchen" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Trypauchenopsis" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Wuhanlinigobius" {
        $8 = "Oxudercidae" 
    } 
    $9 == "Ciliata" {
        $8 = "Gaidropsaridae" 
    } 
    $9 == "Bathygadus" {
        $8 = "Bathygadidae" 
    } 
    $9 == "Gadomus" {
        $8 = "Bathygadidae" 
    } 
    $9 == "Squalogadus" {
        $8 = "Trachyrincidae" 
    } 
    $9 == "Trachyrincus" {
        $8 = "Trachyrincidae" 
    } 
    $9 == "Aetobatus" {
        $8 = "Aetobatidae" 
    } 
    $9 == "Mobula" {
        $8 = "Mobulidae" 
    } 
    $9 == "Rhinoptera" {
        $8 = "Rhinopteridae" 
    } 
    $9 == "Brotula" {
        $8 = "Brotulidae" 
    } 
    $9 == "Petrotyx" {
        $8 = "Bythitidae" 
    } 
    $9 == "Xyelacyba" {
        $8 = "Acanthonidae" 
    } 
    $9 == "Mordacia" {
        $8 = "Mordaciidae" 
    } 
    $9 == "Stereolepis" {
        $8 = "Stereolepididae" 
    } 
    $9 == "Neocentropogon" {
        $8 = "Synanceiidae" 
    } 
    $9 == "Apristurus" {
        $8 = "Pentanchidae" 
    } 
    $9 == "Galeus" {
        $8 = "Triakidae" 
    } 
    $9 == "Halaelurus" {
        $8 = "Pentanchidae" 
    } 
    $9 == "Parmaturus" {
        $8 = "Pentanchidae" 
    } 
    $9 == "Helicolenus" {
        $8 = "Scorpaenidae" 
    } 
    $9 == "Hozukius" {
        $8 = "Scorpaenidae" 
    } 
    $9 == "Sebastes" {
        $8 = "Scorpaenidae" 
    } 
    $9 == "Sebastiscus" {
        $8 = "Scorpaenidae" 
    } 
    $9 == "Sebastolobus" {
        $8 = "Scorpaenidae" 
    } 
    $9 == "Acantholumpenus" {
        $8 = "Lumpenidae" 
    } 
    $9 == "Anisarchus" {
        $8 = "Lumpenidae" 
    } 
    $9 == "Askoldia" {
        $8 = "Opisthocentridae" 
    } 
    $9 == "Cebidichthys" {
        $8 = "Cebidichthyidae" 
    } 
    $9 == "Dictyosoma" {
        $8 = "Cebidichthyidae" 
    } 
    $9 == "Esselenichthys" {
        $8 = "Cebidichthyidae" 
    } 
    $9 == "Kasatkia" {
        $8 = "Opisthocentridae" 
    } 
    $9 == "Leptoclinus" {
        $8 = "Lumpenidae" 
    } 
    $9 == "Lumpenella" {
        $8 = "Lumpenidae" 
    } 
    $9 == "Lumpenus" {
        $8 = "Lumpenidae" 
    } 
    $9 == "Opisthocentrus" {
        $8 = "Opisthocentridae" 
    } 
    $9 == "Pholidapus" {
        $8 = "Opisthocentridae" 
    } 
    $9 == "Plectobranchus" {
        $8 = "Opisthocentridae" 
    } 
    $9 == "Poroclinus" {
        $8 = "Lumpenidae" 
    } 
    $9 == "Xenolumpenus" {
        $8 = "Lumpenidae" 
    } 
    $9 == "Dallia" {
        $8 = "Esocidae" 
    } 
    $9 == "Novumbra" {
        $8 = "Esocidae" 
    } 
    $9 == "Neozoarces" {
        $8 = "Neozoarcidae" 
    } 
    $9 == "Zoarchias" {
        $8 = "Neozoarcidae" 
    } 
    $9 == "Caraibops" {
        $8 = "Synagropidae" 
    } 
    $9 == "Parascombrops" {
        $8 = "Synagropidae" 
    } 
    $9 == "Synagrops" {
        $8 = "Synagropidae" 
    } 
    $9 == "Verilus" {
        $8 = "Malakichthyidae" 
    } 
    $9 == "Barathronus" {
        $8 = "Bythitidae" 
    } 
    $9 == "Apistus" {
        $8 = "Synanceiidae" 
    } 
    $9 == "Cocotropus" {
        $8 = "Synanceiidae" 
    } 
    $9 == "Erisphex" {
        $8 = "Synanceiidae" 
    } 
    $9 == "Paraploactis" {
        $8 = "Synanceiidae" 
    } 
    $9 == "Dolloidraco" {
        $8 = "Harpagiferidae" 
    } 
    $9 == "Neodraco" {
        $8 = "Harpagiferidae" 
    } 
    $9 == "Pogonophryne" {
        $8 = "Harpagiferidae" 
    } 
    $9 == "Aulichthys" {
        $8 = "Hypoptychidae" 
    } 
    $9 == "Horabagrus" {
        $8 = "Horabagridae" 
    } 
    $9 == "Lophiobagrus" {
        $8 = "Claroteidae" 
    } 
    $9 == "Rita" {
        $8 = "Ritidae" 
    } 
    $9 == "Brachionichthys" {
        $8 = "Antennariidae" 
    } 
    $9 == "Caracanthus" {
        $8 = "Scorpaenidae" 
    } 
    $9 == "Spicara" {
        $8 = "Sparidae" 
    } 
    $9 == "Lates" {
        $8 = "Latidae" 
    } 
    $9 == "Psammoperca" {
        $8 = "Latidae" 
    } 
    $9 == "Agoniates" {
        $8 = "Triportheidae" 
    } 
    $9 == "Cynodon" {
        $8 = "Cynodontidae" 
    } 
    $9 == "Elachocharax" {
        $8 = "Crenuchidae" 
    } 
    $9 == "Gnathocharax" {
        $8 = "Acestrorhynchidae" 
    } 
    $9 == "Heterocharax" {
        $8 = "Acestrorhynchidae" 
    } 
    $9 == "Hoplocharax" {
        $8 = "Acestrorhynchidae" 
    } 
    $9 == "Hydrolycus" {
        $8 = "Cynodontidae" 
    } 
    $9 == "Rhaphiodon" {
        $8 = "Cynodontidae" 
    } 
    $9 == "Thayeria" {
        $8 = "Acestrorhamphidae" 
    } 
    $9 == "Nemadactylus" {
        $8 = "Latridae" 
    } 
    $9 == "Caenotropus" {
        $8 = "Chilodidae" 
    } 
    $9 == "Chilodus" {
        $8 = "Chilodidae" 
    } 
    $9 == "Ichthyborus" {
        $8 = "Distichodontidae" 
    } 
    $9 == "Phago" {
        $8 = "Distichodontidae" 
    } 
    $9 == "Gilbertolus" {
        $8 = "Acestrorhynchidae" 
    } 
    $9 == "Roestes" {
        $8 = "Acestrorhynchidae" 
    } 
    $9 == "Datnioides" {
        $8 = "Lobotidae" 
    } 
    $9 == "Bostrychus" {
        $8 = "Butidae" 
    } 
    $9 == "Ophiocara" {
        $8 = "Butidae" 
    } 
    $9 == "Oxyeleotris" {
        $8 = "Butidae" 
    } 
    $9 == "Paloa" {
        $8 = "Butidae" 
    } 
    $9 == "Prionobutis" {
        $8 = "Butidae" 
    } 
    $9 == "Ereunias" {
        $8 = "Rhamphocottidae" 
    } 
    $9 == "Marukawichthys" {
        $8 = "Rhamphocottidae" 
    } 
    $9 == "Hapalogenys" {
        $8 = "Lobotidae" 
    } 
    $9 == "Phreatobius" {
        $8 = "Phreatobiidae" 
    } 
    $9 == "Heteropneustes" {
        $8 = "Clariidae" 
    } 
    $9 == "Phractolaemus" {
        $8 = "Phractolaemidae" 
    } 
    $9 == "Kraemeria" {
        $8 = "Gobiidae" 
    } 
    $9 == "Atypichthys" {
        $8 = "Microcanthidae" 
    } 
    $9 == "Labracoglossa" {
        $8 = "Scorpididae" 
    } 
    $9 == "Medialuna" {
        $8 = "Scorpididae" 
    } 
    $9 == "Microcanthus" {
        $8 = "Microcanthidae" 
    } 
    $9 == "Neatypus" {
        $8 = "Microcanthidae" 
    } 
    $9 == "Scorpis" {
        $8 = "Scorpididae" 
    } 
    $9 == "Tilodon" {
        $8 = "Microcanthidae" 
    } 
    $9 == "Neoclinus" {
        $8 = "Chaenopsidae" 
    } 
    $9 == "Leptochilichthys" {
        $8 = "Alepocephalidae" 
    } 
    $9 == "Lophichthys" {
        $8 = "Antennariidae" 
    } 
    $9 == "Branchiostegus" {
        $8 = "Latilidae" 
    } 
    $9 == "Caulolatilus" {
        $8 = "Latilidae" 
    } 
    $9 == "Lopholatilus" {
        $8 = "Latilidae" 
    } 
    $9 == "Gavialiceps" {
        $8 = "Congridae" 
    } 
    $9 == "Bembrops" {
        $8 = "Bembropidae" 
    } 
    $9 == "Lepidogalaxias" {
        $8 = "Lepidogalaxiidae" 
    } 
    $9 == "Paedocypris" {
        $8 = "Paedocyprididae" 
    } 
    $9 == "Haletta" {
        $8 = "Labridae" 
    } 
    $9 == "Heteroscarus" {
        $8 = "Labridae" 
    } 
    $9 == "Neoodax" {
        $8 = "Labridae" 
    } 
    $9 == "Odax" {
        $8 = "Labridae" 
    } 
    $9 == "Olisthops" {
        $8 = "Labridae" 
    } 
    $9 == "Siphonognathus" {
        $8 = "Labridae" 
    } 
    $9 == "Carcharias" {
        $8 = "Carchariidae" 
    } 
    $9 == "Arapaima" {
        $8 = "Arapaimidae" 
    } 
    $9 == "Heterotis" {
        $8 = "Arapaimidae" 
    } 
    $9 == "Parabembras" {
        $8 = "Bembridae" 
    } 
    $9 == "Percilia" {
        $8 = "Perciliidae" 
    } 
    $9 == "Acanthaphritis" {
        $8 = "Hemerocoetidae" 
    } 
    $9 == "Chrionema" {
        $8 = "Bembropidae" 
    } 
    $9 == "Osopsaron" {
        $8 = "Hemerocoetidae" 
    } 
    $9 == "Gargariscus" {
        $8 = "Triglidae" 
    } 
    $9 == "Paraheminodus" {
        $8 = "Triglidae" 
    } 
    $9 == "Peristedion" {
        $8 = "Triglidae" 
    } 
    $9 == "Satyrichthys" {
        $8 = "Triglidae" 
    } 
    $9 == "Scalicus" {
        $8 = "Triglidae" 
    } 
    $9 == "Ichthyococcus" {
        $8 = "Ichthyococcidae" 
    } 
    $9 == "Pollichthys" {
        $8 = "Vinciguerriidae" 
    } 
    $9 == "Polymetme" {
        $8 = "Yarrellidae" 
    } 
    $9 == "Vinciguerria" {
        $8 = "Vinciguerriidae" 
    } 
    $9 == "Yarrella" {
        $8 = "Yarrellidae" 
    } 
    $9 == "Eques" {
        $8 = "Sciaenidae" 
    } 
    $9 == "Rhamdia" {
        $8 = "Heptapteridae" 
    } 
    $9 == "Bleekeria" {
        $8 = "Ammodytidae" 
    } 
    $9 == "Zapteryx" {
        $8 = "Trygonorrhinidae" 
    } 
    $9 == "Neoachiropsetta" {
        $8 = "Achiropsettidae" 
    } 
    $9 == "Rhynchobatus" {
        $8 = "Rhinidae" 
    } 
    $9 == "Mallotus" {
        $8 = "Osmeridae" 
    } 
    $9 == "Eutropiichthys" {
        $8 = "Ailiidae" 
    } 
    $9 == "Schindleria" {
        $8 = "Gobiidae" 
    } 
    $9 == "Aethaloperca" {
        $8 = "Epinephelidae" 
    } 
    $9 == "Anyperodon" {
        $8 = "Epinephelidae" 
    } 
    $9 == "Aporops" {
        $8 = "Grammistidae" 
    } 
    $9 == "Aulacocephalus" {
        $8 = "Liopropomatidae" 
    } 
    $9 == "Bathyanthias" {
        $8 = "Liopropomatidae" 
    } 
    $9 == "Cephalopholis" {
        $8 = "Epinephelidae" 
    } 
    $9 == "Chromileptes" {
        $8 = "Epinephelidae" 
    } 
    $9 == "Dermatolepis" {
        $8 = "Epinephelidae" 
    } 
    $9 == "Diploprion" {
        $8 = "Liopropomatidae" 
    } 
    $9 == "Epinephelus" {
        $8 = "Epinephelidae" 
    } 
    $9 == "Grammistes" {
        $8 = "Grammistidae" 
    } 
    $9 == "Hyporthodus" {
        $8 = "Epinephelidae" 
    } 
    $9 == "Liopropoma" {
        $8 = "Liopropomatidae" 
    } 
    $9 == "Mycteroperca" {
        $8 = "Epinephelidae" 
    } 
    $9 == "Paranthias" {
        $8 = "Epinephelidae" 
    } 
    $9 == "Plectropomus" {
        $8 = "Epinephelidae" 
    } 
    $9 == "Pogonoperca" {
        $8 = "Grammistidae" 
    } 
    $9 == "Pseudogramma" {
        $8 = "Grammistidae" 
    } 
    $9 == "Rypticus" {
        $8 = "Grammistidae" 
    } 
    $9 == "Saloptia" {
        $8 = "Epinephelidae" 
    } 
    $9 == "Suttonia" {
        $8 = "Grammistidae" 
    } 
    $9 == "Triso" {
        $8 = "Epinephelidae" 
    } 
    $9 == "Variola" {
        $8 = "Epinephelidae" 
    } 
    $9 == "Ectreposebastes" {
        $8 = "Scorpaenidae" 
    } 
    $9 == "Lythrichthys" {
        $8 = "Scorpaenidae" 
    } 
    $9 == "Setarches" {
        $8 = "Scorpaenidae" 
    } 
    $9 == "Tetrabrachium" {
        $8 = "Antennariidae" 
    } 
    $9 == "Ablabys" {
        $8 = "Synanceiidae" 
    } 
    $9 == "Hypodytes" {
        $8 = "Synanceiidae" 
    } 
    $9 == "Neovespicula" {
        $8 = "Synanceiidae" 
    } 
    $9 == "Paracentropogon" {
        $8 = "Synanceiidae" 
    } 
    $9 == "Tetraroge" {
        $8 = "Synanceiidae" 
    } 
    $9 == "Lasiognathus" {
        $8 = "Oneirodidae" 
    } 
    $9 == "Zenion" {
        $8 = "Zeniontidae" 
    } 
    $9 == "Pareutropius" {
        $8 = "Schilbeidae" 
    } 
    $9 == "Schilbe" {
        $8 = "Schilbeidae" 
    } 
   1' $output_dir/aligned_tmp4.txt > $output_dir/aligned_tmp5.txt
   
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Changing families to match Eschmeyer orders for all species"

    awk -F '\t' -v OFS='\t' '
    $8 == "Antigoniidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Caproidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Chaetodontidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Leiognathidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Drepaneidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Ephippidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Gerreidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Lobotidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Haemulidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Lutjanidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Callanthiidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Emmelichthyidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Latilidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Malacanthidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Monodactylidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Moronidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Pomacanthidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Scatophagidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Sciaenidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Siganidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Sillaginidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Cepolidae (in: bony fishes)" {
        $8 = "Cepolidae";
        $7 = "Acanthuriformes"
    } 
    $8 == "Priacanthidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Lethrinidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Nemipteridae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Sparidae" {
        $7 = "Acanthuriformes"
    } 
    $8 == "Hemerocoetidae" {
        $7 = "Acropomatiformes"
    } 
    $8 == "Scombropidae" {
        $7 = "Acropomatiformes"
    } 
    $8 == "Holocentridae" {
        $7 = "Beryciformes"
    } 
    $8 == "Embiotocidae" {
        $7 = "Blenniiformes"
    } 
    $8 == "Grammatidae" {
        $7 = "Blenniiformes"
    } 
    $8 == "Opistognathidae" {
        $7 = "Blenniiformes"
    } 
    $8 == "Plesiopidae" {
        $7 = "Blenniiformes"
    } 
    $8 == "Pomacentridae" {
        $7 = "Blenniiformes"
    } 
    $8 == "Pseudochromidae" {
        $7 = "Blenniiformes"
    } 
    $8 == "Istiophoridae" {
        $7 = "Carangiformes"
    } 
    $8 == "Xiphiidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Centropomidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Lactariidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Latidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Leptobramidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Menidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Polynemidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Sphyraenidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Toxotidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Achiridae" {
        $7 = "Carangiformes"
    } 
    $8 == "Achiropsettidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Bothidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Citharidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Cyclopsettidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Cynoglossidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Grammatobothidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Monolenidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Paralichthyidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Pleuronectidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Poecilopsettidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Psettodidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Rhombosoleidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Samaridae" {
        $7 = "Carangiformes"
    } 
    $8 == "Scophthalmidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Soleidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Taeniopsettidae" {
        $7 = "Carangiformes"
    } 
    $8 == "Polycentridae" {
        $7 = "Cichliformes"
    } 
    $8 == "Stylephoridae" {
        $7 = "Gadiformes"
    } 
    $8 == "Apogonidae" {
        $7 = "Gobiiformes"
    } 
    $8 == "Kurtidae" {
        $7 = "Gobiiformes"
    } 
    $8 == "Gonorynchidae" {
        $7 = "Gonorynchiformes"
    } 
    $8 == "Ammodytidae" {
        $7 = "Labriformes"
    } 
    $8 == "Ammodytidae" {
        $7 = "Labriformes"
    } 
    $8 == "Cheimarrichthyidae" {
        $7 = "Labriformes"
    } 
    $8 == "Pinguipedidae" {
        $7 = "Labriformes"
    } 
    $8 == "Uranoscopidae" {
        $7 = "Labriformes"
    } 
    $8 == "Lepisosteidae" {
        $7 = "Lepisosteiformes"
    } 
    $8 == "Ambassidae" {
        $7 = "Mugiliformes"
    } 
    $8 == "Bembropidae" {
        $7 = "Perciformes"
    } 
    $8 == "Pristiophoridae" {
        $7 = "Pristiophoriformes"
    } 
    $8 == "Esocidae" {
        $7 = "Salmoniformes"
    } 
    $8 == "Umbridae" {
        $7 = "Salmoniformes"
    } 
    $8 == "Squatinidae" {
        $7 = "Squatiniformes"
    } 
    $8 == "Platyrhinidae" {
        $7 = "Torpediniformes"
    } 
    1' $output_dir/aligned_tmp5.txt > $output_dir/aligned_tmp6.txt

    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Changing orders to match Eschmeyer classes for all species"

    awk -F '\t' -v OFS='\t' '
    $7 == "Carcharhiniformes" {
        $6 = "Elasmobranchii"
    } 
    $7 == "Chimaeriformes" {
        $6 = "Holocephali"
    } 
    $7 == "Echinorhiniformes" {
        $6 = "Elasmobranchii"
    } 
    $7 == "Heterodontiformes" {
        $6 = "Elasmobranchii"
    } 
    $7 == "Hexanchiformes" {
        $6 = "Elasmobranchii"
    } 
    $7 == "Lamniformes" {
        $6 = "Elasmobranchii"
    } 
    $7 == "Myliobatiformes" {
        $6 = "Elasmobranchii"
    } 
    $7 == "Orectolobiformes" {
        $6 = "Elasmobranchii"
    } 
    $7 == "Pristiophoriformes" {
        $6 = "Elasmobranchii"
    } 
    $7 == "Rajiformes" {
        $6 = "Elasmobranchii"
    } 
    $7 == "Rhinopristiformes" {
        $6 = "Elasmobranchii"
    } 
    $7 == "Squaliformes" {
        $6 = "Elasmobranchii"
    } 
    $7 == "Squatiniformes" {
        $6 = "Elasmobranchii"
    } 
    $7 == "Torpediniformes" {
        $6 = "Elasmobranchii"
    } 
    $7 == "Polypteriformes" {
        $6 = "Cladistii"
    } 
    $7 == "Petromyzontiformes" {
        $6 = "Petromyzonti"
    } 
    1' $output_dir/aligned_tmp6.txt > $output_dir/aligned_curated.txt
	rm $output_dir/aligned_tmp*.txt
else
    print_yellow "File $output_dir/aligned_curated.txt already exists. Skipping blacklist-removal, additional amplicon fixing, synonymy, subspecies, genus-to family correction, family to order correction, and order to class correction step"
fi

# Filtering the reference database
# Check if the filtered output file already exists
if [ ! -f "$output_dir/FilteredRefs_CRABS.txt" ]; then
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Filtering based on the following criteria:"
    print_yellow "[$current_datetime] Minimum length: 150"
    print_yellow "[$current_datetime] Maximum length: 500"
    print_yellow "[$current_datetime] Maximum Ns in sequence: 1"
    print_yellow "[$current_datetime] Remove environmental sequences"
    print_yellow "[$current_datetime] Maximum two out of seven NAs in taxonomic ranks"
    print_yellow "[$current_datetime] Removal of records without species ID is currently deactivated"
    crabs \
    --filter \
    --input "$output_dir"/aligned_curated.txt \
    --output "$output_dir"/FilteredRefs_CRABS.txt \
    --minimum-length 150 \
    --maximum-length 500 \
    --maximum-n 1 \
    --environmental \
    --rank-na 3
	cat "$output_dir"/FilteredRefs_CRABS.txt | sort -k7 > "$output_dir"/tmp.txt
	mv "$output_dir"/tmp.txt "$output_dir"/FilteredRefs_CRABS.txt
	awk -F'\t' '{print $1 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' "$output_dir"/FilteredRefs_CRABS.txt | sort > "$output_dir"/refs_sorted.txt
else
    print_yellow "File $output_dir/FilteredRefs_CRABS.txt already exists. Skipping filtering."
fi

# Exporting the database for use
# Check if the BLAST_TAX_12S output file already exists
if [ ! -f "$output_dir/BLAST_TAX_12S/BLAST_TAX_12S.nsq" ]; then
    current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")
    print_yellow "[$current_datetime] Exporting filtered database to blast-tax format ..."
	mkdir -p "$output_dir"/BLAST_TAX_12S/
    crabs \
    --export \
    --input "$output_dir"/FilteredRefs_CRABS.txt \
    --output "$output_dir"/BLAST_TAX_12S/BLAST_TAX_12S \
    --export-format 'blast-tax'
	mv "$output_dir"/refs_sorted.txt "$output_dir"/BLAST_TAX_12S/refs_sorted.txt
	print_green "[$current_datetime] Success! Enjoy your new reference database."
else
    print_yellow "File $output_dir/BLAST_TAX_12S already exists. Skipping database export."
fi
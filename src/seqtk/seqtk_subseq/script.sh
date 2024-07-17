#!/bin/bash

## VIASH START
## VIASH END

# Function to check if a file is a valid BED file
check_bed_file() {
    local file="$1"
    num_columns=$(head -n 1 "$file" | cut -f 1- | tr '\t' '\n' | wc -l)

    # Check if the file exists
    if [[ ! -f "$file" ]]; then
        echo "Error: The specified file does not exist."
        return 1
    fi

    # Check if the file is non-empty
    if [[ ! -s "$file" ]]; then
        echo "Error: The specified file is empty."
        return 1
    fi

    # Check if the file is a valid BED file (minimum 3 tab-separated columns)
    if [[ $num_columns -lt 3 ]]; then
        echo "The specified file is not a valid BED file. It should have at least three tab-separated columns."
        return 1
    fi

    # Additional check: Ensure that the 6th column (if present) is either + or -
    if [[ $num_columns -eq 6 ]]; then
        while IFS=$'\t' read -r -a columns; do
            columns[5]=${columns[5]%$'\n'}
            if [[ ${columns[5]} == "+" || ${columns[5]} == "-" ]]; then
                return 0
            else 
                echo "Error: The 6th column of the specified BED file should be either + or -."
                echo "Offending line: ${columns[5]}"
                return 1
            fi
        done < "$file"
    fi

    return 0
}

# Function to check if a file is a valid list of FASTA file IDs
check_fasta_id_list() {
    local file="$1"

    # Check if the file exists
    if [[ ! -f "$file" ]]; then
        echo "Error: The specified file does not exist."
        return 1
    fi

    # Check if the file is non-empty
    if [[ ! -s "$file" ]]; then
        echo "Error: The specified file is empty."
        return 1
    fi

    # Additional check: Ensure that each line contains only one word (FASTA ID)
    if ! awk 'NF != 1 { exit 1 }' "$file"; then
        return 1
    fi

    return 0
}

# Check if the par_name_list is given and validate accordingly
if [[ -n "$par_name_list" ]]; then
    if check_fasta_id_list "$par_name_list"; then
        echo "The specified file is a valid list of FASTA IDs."
    elif check_bed_file "$par_name_list"; then
        echo "The specified file is a valid BED file."
    else
        echo "Error: The specified file is neither a valid BED file nor a valid list of FASTA IDs."
        exit 1
    fi
fi

[[ "$par_tab" == "false" ]] && unset par_tab
[[ "$par_strand_aware" == "false" ]] && unset par_strand_aware

seqtk subseq \
    ${par_tab:+-t} \
    ${par_strand_aware:+-s} \
    ${par_sequence_line_length:+-l "$par_sequence_line_length"} \
    "$par_input" \
    "$par_name_list" \
    > "$par_output"
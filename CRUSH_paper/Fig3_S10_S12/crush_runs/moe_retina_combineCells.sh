#!/bin/bash


check_directory() {
    if [ ! -d "$1" ]; then
        echo "Error: Directory $1 does not exist"
        exit 1
    fi
}

check_output_path() {
    output_dir=$(dirname "$1")
    if [ ! -d "$output_dir" ]; then
        echo "Error: Output directory $output_dir does not exist"
        exit 1
    fi
    if [ -f "$1" ]; then
        read -p "Output file exists. Overwrite? (y/n) " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            exit 1
        fi
    fi
}

validate_resolution() {
    if ! [[ "$1" =~ ^[0-9]+$ ]]; then
        echo "Error: Resolution must be a positive integer"
        exit 1
    fi
}

process_files() {
    local dir=$1
    local resolution=$2
    find "$dir" -name "GSM*" -type d -exec sh -c '
        resolution=$1
        shift
        for d; do
            if ls "$d"/*mergedCrush_"$resolution".bedgraph 1> /dev/null 2>&1; then
                cat "$d"/*mergedCrush_"$resolution".bedgraph
            else
                echo "Warning: No files with resolution $resolution found in $d" >&2
            fi
        done
    ' sh "$resolution" {} + \
    | awk '{if ($1 == 1) print $0}' \
    | sort -k 2bn,2b --stable
}

if [ $# -lt 3 ]; then
    echo "Usage: $0 output_file resolution directory1 directory2 [directory3 ...]"
    echo "Example: $0 output.txt 100000 /path/to/retina /path/to/moe"
    exit 1
fi

output_file="$1"
resolution="$2"
shift 2

check_output_path "$output_file"
validate_resolution "$resolution"

for dir in "$@"; do
    check_directory "$dir"
done

echo "Processing files with resolution ${resolution} and saving to $output_file..."
(for dir in "$@"; do
    process_files "$dir" "$resolution"
done | groupBy -i stdin -g 2 -c 4 -o collapse | sed 's/,/\t/g') > "$output_file"

if [ $? -eq 0 ]; then
    echo "Successfully saved results to $output_file"
    echo "Number of lines in output: $(wc -l < "$output_file")"
else
    echo "Error occurred while saving to $output_file"
    exit 1
fi

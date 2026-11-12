#!/bin/bash
hic_dir="moe/"
output_dir="crush"
crush_script="Crush_FV_OnePunch_2.bash"
eigenfile="EigenVCSQRT_100kb.bedgraph"

genome_sizes="mm10_num.sizes"
genes_bed="mm10genes.bed"
bbins_file="Bbins_mm10_1SD.bed"

if [[ ! -f "$eigenfile" ]]; then
    echo "Error: Constant eigenfile not found: $eigenfile"
    exit 1
fi

for hicfile in "$hic_dir"/*.hic
do
    if [[ ! -f "$hicfile" ]]; then
        echo "Error: Hi-C file not found: $hicfile"
        continue
    fi

    filename=$(basename -- "$hicfile")
    filename_no_ext="${filename%.*}"
    
    cd "$output_dir"
    output_folder="${filename_no_ext}"
    mkdir -p "$output_folder"

    echo "Running CRUSH for ${filename_no_ext}"
    cd "$output_folder" || exit 1

    (
        /bin/bash "$crush_script" \
            -i "$hicfile" \
            -r 25000 \
            -cpu 20 \
            -w 5 \
            -C 0 \
            -g "$genome_sizes" \
            -a "$genes_bed" \
            -b "$bbins_file" \
            -e "$eigenfile" \
    )

    echo "Completed processing for ${filename_no_ext}"
done

echo "All files have been processed."

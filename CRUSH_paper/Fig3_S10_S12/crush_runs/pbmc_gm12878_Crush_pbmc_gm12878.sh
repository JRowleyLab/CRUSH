#!/bin/bash
hic_dir="hicfiles/"

gm12878_output_dir="crush"
pbmc_output_dir="crush"

crush_script="Crush_FV_OnePunch_2.bash"

eigenfile_gm12878="EigenVCSQRT_100kb_gm12878.bedgraph"
eigenfile_pbmc="EigenVCSQRT_100kb_pbmc.bedgraph"

genome_sizes="hg19.sizes"
genes_bed="hg19genes.bed"
bbins_file="hg19_nochr.fa"

if [[ ! -f "$eigenfile_gm12878" ]]; then
    echo "Error: GM12878 eigenfile not found: $eigenfile_gm12878"
    exit 1
fi

if [[ ! -f "$eigenfile_pbmc" ]]; then
    echo "Error: PBMC eigenfile not found: $eigenfile_pbmc"
    exit 1
fi

mkdir -p "$gm12878_output_dir"
mkdir -p "$pbmc_output_dir"

for hicfile in "$hic_dir"/*.hic
do
    if [[ ! -f "$hicfile" ]]; then
        echo "Error: Hi-C file not found: $hicfile"
        continue
    fi

    filename=$(basename -- "$hicfile")
    filename_no_ext="${filename%.*}"

    if [[ "$filename_no_ext" == gm12878* ]]; then
        eigenfile="$eigenfile_gm12878"
        output_dir="$gm12878_output_dir"
    elif [[ "$filename_no_ext" == pbmc* ]]; then
        eigenfile="$eigenfile_pbmc"
        output_dir="$pbmc_output_dir"
    else
        echo "Skipping file $filename_no_ext: Not a GM12878 or PBMC file."
        continue
    fi

    output_folder="$output_dir/${filename_no_ext}"
    mkdir -p "$output_folder"

    echo "Running CRUSH for ${filename_no_ext} using eigenfile ${eigenfile}"
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
            -e "$eigenfile"
    )

    echo "Completed processing for ${filename_no_ext}"
done

echo "All files have been processed."

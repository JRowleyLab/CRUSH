#!/bin/bash
source_folder="crush"
metadata_file="MOE_Retina_metadataCells.tsv"
moe_output="moe_Totalcombined_25kb.tsv"
retina_output="retina_Totalcombined_25kb.tsv"


> "$moe_output"
> "$retina_output"

declare -A cell_mapping
while IFS=$'\t' read -r cell_number cell_type; do
    cell_mapping["$cell_number"]="$cell_type"
done < "$metadata_file"

for folder in "$source_folder"/*/; do
    folder_name=$(basename "$folder")
    cell_number=$(echo "$folder_name" | grep -o 'cell_[0-9]\+')

    if [[ -n "${cell_mapping[$cell_number]}" ]]; then
        cell_type="${cell_mapping[$cell_number]}"
        if [[ "$cell_type" == "MOE" ]]; then
            echo "Processing MOE folder: $folder"
            cat "$folder"*mergedCrush_25000.bedgraph | awk '{if ($1 == 1) print $0}' | \
            sort -k 2bn,2b --stable | groupBy -i stdin -g 2 -c 4 -o collapse | \
            sed 's/,/\t/g' >> "$moe_output"
        elif [[ "$cell_type" == "retina" ]]; then
            echo "Processing Retina folder: $folder"
            cat "$folder"*mergedCrush_25000.bedgraph | awk '{if ($1 == 1) print $0}' | \
            sort -k 2bn,2b --stable | groupBy -i stdin -g 2 -c 4 -o collapse | \
            sed 's/,/\t/g' >> "$retina_output"
        fi
    else
        echo "Unknown cell number or missing metadata for $folder"
    fi
done

echo "Processing completed. Outputs:"
echo "  MOE: $moe_output"
echo "  Retina: $retina_output"

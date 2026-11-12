chrom=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y)
sizes=(195471971 182113224 160039680 156508116 151834684 149736546 145441459 129401213 124595110 130694993 122082543 120129022 120421639 124902244 104043685 98207768 94987271 90702639 61431566 171031299 91744698)

moe_dir="moe"
retina_dir="retina"

output_dir="Possumm"

juicer_tool="juicer_tools_1.14.08.jar"
eigen_tool="mainEigen.exe"

for dir in "$moe_dir" "$retina_dir"; do
  for hicfile in "$dir"/*.hic; do
    sample_name=$(basename "$hicfile" .hic)
    output_path="$output_dir/$sample_name"
    
    mkdir -p "$output_path"
    cd "$output_path"

    echo "Processing $sample_name at 500kb resolution..."
    outie3="Eigen_${sample_name}_500kb.bedgraph"

    for j in "${!chrom[@]}"; do
      echo "Processing chromosome ${chrom[j]} with size ${sizes[j]}..."
      java -jar "$juicer_tool" dump oe NONE "$hicfile" "${chrom[j]}" "${chrom[j]}" BP 500000 "${chrom[j]}.out"
      "$eigen_tool" -s "${sizes[j]}" "${chrom[j]}.out" "${chrom[j]}.tmp" 500000
      cat "${chrom[j]}.tmp" | sed 's/ //g' | awk -v var="${chrom[j]}" '{print var"\t"(NR-1)*500000"\t"NR*500000"\t"$1}' | grep -v NaN > "${chrom[j]}.chrom_bedgraph"
      rm "${chrom[j]}.out" "${chrom[j]}.tmp"
    done

    cat *.chrom_bedgraph > "$outie3"
    rm *.chrom_bedgraph

    echo "Finished processing $sample_name."
  done
done

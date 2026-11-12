#!/bin/bash
chrom=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y)
sizes=(195471971 182113224 160039680 156508116 151834684 149736546 145441459 129401213 124595110 130694993 122082543 120129022 120421639 124902244 104043685 98207768 94987271 90702639 61431566 171031299 91744698)

hic_files=(
    "4DNFIIRMBE8Y_OSN_precursors.hic"
)

output_base="Possumm"

juicer_tools="juicer_tools_1.14.08.jar"
eigenvector_bin="mainEigen.exe"

for hicfile in "${hic_files[@]}"; do
    hic_name=$(basename "$hicfile" .hic)
    
    output_dir="$output_base/$hic_name"
    mkdir -p "$output_dir"
    
    cd "$output_dir" || exit

    echo "Processing Hi-C file: $hicfile"
    
    outie2="Eigen_${hic_name}_500kb.bedgraph"

    for j in "${!chrom[@]}"; do
        echo "Processing Chromosome: ${chrom[j]}, Size: ${sizes[j]}"
        
        java -jar "$juicer_tools" dump oe NONE "$hicfile" "${chrom[j]}" "${chrom[j]}" BP 500000 "${chrom[j]}.out"
        
        "$eigenvector_bin" -s "${sizes[j]}" "${chrom[j]}.out" "${chrom[j]}.tmp" 500000
        
        awk -v chrom=${chrom[j]} '{print chrom "\t" (NR-1)*500000 "\t" NR*500000 "\t" $1}' "${chrom[j]}.tmp" | grep -v NaN > "${chrom[j]}.chrom_bedgraph"
        
        rm "${chrom[j]}.out" "${chrom[j]}.tmp"
    done

    cat *.chrom_bedgraph > "$outie2"
    
    rm *.chrom_bedgraph

    echo "Compartment calling completed for: $hic_name"
done

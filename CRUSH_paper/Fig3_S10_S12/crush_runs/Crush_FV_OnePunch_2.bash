#!/bin/bash

version=1.0

hicpath=0
myres=0
help=0
genome=0
cpu=1
doNotMerge=0
adjustment=0
distance=0
upperlim=0
trackline=1
threshold=10
totreads=0
lowerthresh=0
window=0
verbose=0
switch=1
qthresh=0.05
whichchoice=`echo "o"`
cleanup=1
endZ=0
exclbed=0
crushdir=`echo "CRUSHtmp_""$RANDOM"`
outpre=""
norm=NONE
coarsestres=2500000
pcalculation=1
reshift=0
smoothing=0
eigenres=100000

function shutdown() {
    tput cnorm
}
trap shutdown EXIT

function cursorBack() {
    echo -en "\033[$1D"
}

function spinner() {
    local LC_CTYPE=C
    local pid=$1

    case $(($RANDOM % 17)) in
        0) local spin='⠁⠂⠄⡀⢀⠠⠐⠈'; local charwidth=3 ;;
        1) local spin='-\|/'; local charwidth=1 ;;
        2) local spin="▁▂▃▄▅▆▇█▇▆▅▄▃▂▁"; local charwidth=3 ;;
        3) local spin="▉▊▋▌▍▎▏▎▍▌▋▊▉"; local charwidth=3 ;;
        4) local spin='←↖↑↗→↘↓↙'; local charwidth=3 ;;
        5) local spin='▖▘▝▗'; local charwidth=3 ;;
        6) local spin='┤┘┴└├┌┬┐'; local charwidth=3 ;;
        7) local spin='◢◣◤◥'; local charwidth=3 ;;
        8) local spin='◰◳◲◱'; local charwidth=3 ;;
        9) local spin='◴◷◶◵'; local charwidth=3 ;;
        10) local spin='◐◓◑◒'; local charwidth=3 ;;
        11) local spin='⣾⣽⣻⢿⡿⣟⣯⣷'; local charwidth=3 ;;
        12) local spin='CCCRRRUUUSSSHHH'; local charwidth=1 ;;
        13) local spin='RRROOOWWWLLLEEEYYYLLLAAABBB'; local charwidth=1 ;;
        14) local spin='SSSAAAIIITTTAAAMMMAAA'; local charwidth=1 ;;
        15) local spin='MMMOOONNNKKKEEEYYYDDDLLLUUUFFFFFFYYY'; local charwidth=1 ;;
        16) local spin='RRROOORRROOONNNOOOAAAZZZOOORRROOO'; local charwidth=1 ;;
    esac

    local i=0
    tput civis
    while kill -0 $pid 2>/dev/null; do
        local i=$(((i + $charwidth) % ${#spin}))
        printf "%s" "${spin:$i:$charwidth}"
        cursorBack 1
        sleep .1
    done
    tput cnorm
    wait $pid
    return $?
}

function usage {
    echo -e "\n\nusage : crush -i HIC  -g SIZEFILE -a ABED -b BBED | FASTA -r FINERESOLUTION [-e EIGENVECTORBED] [-cpu CPU] [-w WINDOW] [-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo
    echo "CRUSH will create a temprory folder in your current directory, so make sure you have write access to your current directory."
    echo "See https://github.com/JRowleyLab/CRUSH/ for full documentation on CRUSH"
    echo "-----------------------------------------"
    echo "OPTIONS:"

    echo "-h|--help                :  Display this help menu"
    echo " "
    echo "--------------------------REQUIRED PARAMETERS------------------"
    echo "-i|--hic                 :  Input .hic/.cooler/.mcool file by specifying the path. (e.g., '/path/to/ file.hic or file.mcool or file.cooler)"
    echo "-g|--genomesize          :  Specify path to a chromosome size file with two columns corresponding to chromosome and size respectively."
    echo "-a|--initialA            :  Specify path to a bed file with the regions for initializing A. For example, gene annotations e.g. hg19genes.bed."
    echo "-b|--initialB            :  Specify path to either a fasta file or to a bed file for initializing B. If you specify a fasta file, we will calculate intialB from gc content." 
    echo "-e|--eigenfile           :  Specify path to an eigenfile to initialize A and B states. If you don't specify a file, CRUSH will calculate EigenVectors by default and pick the best PC possible." 
    echo "-r|--res                 :  Resolution desired."
    echo "---------------" 

    echo " "
    echo "--------------------------OPTIONAL PARAMETERS------------------"

    echo "-o|--outpre              :  Set this if you want to specify a prefix for the output files"
    echo "-c|--cpu                 :  Set the value for cpu number of threads to use. Default is 1."
    echo "-n|--no-merge            :  Set this option to 1 to keep each resolution as a separate output file. Default is to merge in a way that provides maximum resolution."
    echo "-A|--adjustment          :  Set this option to 1 to include a adjustment of CRUSH values at the end. This adjustment shifts values based on any internal skewing of the data. Do not set this if using CRUSH to compare between two Hi-C maps."
    echo "-d|--distance            :  Using this option will filter out the distance next to the diagonal. Default is 0 which considers everything."
    echo "-u|--upperlim            :  The upperlimit of the distance away from the diagonal to consider. Default is 0 so that it considers the whole chromosome."
    echo "-t|--trackline           :  Set to 0 if you want to disable printing a bedgraph trackline header."
    echo "-T|--threshold           :  Distance normalized threshold to filter out extreme outliers."
    echo "-l|--lowerthresh         :  Set the value for lowerthresh."
    echo "-w|--window              :  Set to perform a sliding window average of the scores at individual resolutions. Default is to calculate the appropriate window based on sequencing depth. Set to 1 to remove sliding window."
    echo "-v|--verbose             :  Set to 1 to enable verbose mode showing extensive messages."
    echo "-S|--switch              :  Set this option to 0 for bypassing re-initialization. Default is 1."  
    echo "-N|--norm                :  Set the normalization scheme to use. Options are NONE, VC, VC_SQRT, KR, SCALE. Default is NONE."
    echo "-C|--cleanup             :  Set the value for cleanup"
    echo "-E|--endZ                :  Set the value for endZ"
    echo "-x|--exclbed             :  Set the value for exclbed"
    echo "-p|--pcalculation        :  Set this option to 0 for bypassing pvalue calculation. Default is 1."
    echo "-q|--qvalue              :  Set the qvalue threshold. default 0.05. Set to 0 to not perform qvalue filtering. The qvalues will be reported as a separate track regardelss."
    echo "-u|--use                 :  Whether to use of overwrite existing GI tracks previously calculated at individual resolutions. Set this option to u to use previous calculations. Default is to recalculate. This option is useful for merging resolutions."     
    echo "-f|--tmpfolder           :  Set this if you want to name the temporary folder youself. Make sure it doesn't already exist in your current working directory. Default is to name it CRUSHtmp with a randomnumber."
    echo "-m|--maxres              :  Set this to the coarsest resolution you want to consider. Default is to check every resolution between 1000000 and your desired resolution to inform each other." 
    echo "-R|--reshift             :  Set this to option to 1 for turning on end-shifting. Default is 0." 
    echo "-s|--smoothing           :  Set this to option to 1 for smoothening the compartments. Default is 0." 
    echo "-Z|--eigenres            :  Set a resolution of your choice for eigenvector calculations. Default is 100kb (100000)."
}

while test $# -gt 0; do
    case "$1" in
        -h|--help)
            help
            exit 0
            ;;
        -i|--hic)
            hicpath=`readlink -e $2`
            ;;
        -o|--outpre)
            outpre=$2
            ;;
        -r|--res)
            res=$2
            ;;
        -g|--genomesize)
            sizefile=`readlink -e $2`
            ;;
        -a|--initialA)
            genesfile=`readlink -e $2`
            ;;
        -b|--initialB)
            fastafile=`readlink -e $2`
            ;;
        -e|--eigenfile) 
            eigenfile=`readlink -e $2`
            shift
            ;;
        -q|--qvalue)
            qthresh=$2
            ;;
        -x|--exclbed)
            exclbed=`readlink -e $2`
            ;;
        -c|--cpu)
            cpu=$2
            ;;
        -n|--no-merge)
            doNotMerge=$2
            ;;
        -A|--adjustment)
            adjustment=$2
            ;;
        -d|--distance)
            distance=$2
            ;;
        -u|--upperlim)
            upperlim=$2
            ;;
        -l|--lowerthresh)
            lowerthresh=$2
            ;;
        -t|--trackline)
            trackline=$2
            ;;
        -T|--threshold)
            threshold=$2
            ;;
        -w|--window)
            window=$2
            ;;
        -v|--verbose)
            verbose=$2
            ;;
        -S|--switch)
            switch=$2
            ;;
        -N|--norm)
            norm=$2
            ;;
        -use)
            whichchoice=$2
            ;;
        -f|--tmpfolder)
            crushdir=$2
            ;;
        -C|--cleanup)
            cleanup=$2
            ;;
        -E|--endZ)
            endZ=$2
            ;;
        -m|--maxres)
            coarsestres=$2
            ;;
        -p|--pcalculation)
            pcalculation=$2
            ;;
        -R|--reshift)
            reshift=$2
            ;;
        -s|--smoothing)
            smoothing=$2
            ;;
        -Z|--eigenres) 
            eigenres=$2
            ;;
    esac
    shift
done

if [ $hicpath == 0 ]; then
    echo "You must specify a .hic file!......"
    exit 0
fi

if [ $res == 0 ]; then
    echo "You must specify a resolution!......"
    exit 0
fi

if [ $switch == 0 ]; then
    echo "Warning: Re-evaluation & Re-iteration of A and B Bins has been deactivated. We recommend leaving this parameter alone to achieve the highest possible confidence and resolution."
fi

if [ -z "$sizefile" ] || [ -z "$genesfile" ] || [ -z "$fastafile" ]; then
    echo "Missing one or more required arguments: --genomesize, --initialA, or --initialB!..."
    help
    exit 1
fi

if [ $doNotMerge == 1 ]; then
    echo "Warning: doNotMerge set. We recommend merging to achieve the highest possible confidence and resolution. Only set doNotMerge if you plan to merge the individual resolutions later."
fi

if [ $adjustment == 1 ]; then
    echo "Warning: adjusting end compartmental values which may result in some loss of quantitative power. We recommend adjustment only for troubleshooting or when not comparing two samples."
fi

echo "CRUSH_v""$version"" --hic ""$hicpath"" --res ""$res"" --genomesize ""$sizefile"" --initialA ""$genesfile"" --initialB ""$fastafile"" --cpu ""$cpu"" --no-merge ""$doNotMerge"" --adjustment ""$adjustment"" --distance ""$distance"" --upperlim ""$upperlim"" --lowerthresh ""$lowerthresh"" --trackline ""$trackline"" --threshold ""$threshold"" --window ""$window"" --switch ""$switch"" --maxres ""$coarsestres"" --norm ""$norm"" --eigenres ""$eigenres"" --reshift ""$reshift"" --smoothing ""$smoothing" > "$outpre"CRUSHparamters.txt

juiceorcool=`echo "$hicpath" | sed 's/\./\t/g' | sed 's/\./\t/g' | awk '{if ($NF == "hic") print 0; else if ($NF == "mcool") print 1; else print 2}'`

if [ $juiceorcool == 2 ]; then
    echo "You must choose a Hi-C file that ends in either .hic (juicer) or .mcool (cooler)"
    exit 0
elif [ $juiceorcool == 0 ]; then
    echo "Identified .hic juicer format. We will use the pythonic hicstraw to extract the data and assume you have it installed, using pip install hic-straw."
else
    echo "Identified .mcool format. We will use cooler to extract the data and assume you have it installed in your path."
fi



run_dumper() {
    local readtype=$1
    local norm=$2
    local hicinnie=$3
    local chrom=$4
    local res=$5
    local dumpoutie=$6
    local dumperrors=$7

    echo "Dumping $readtype reads using $norm normalization at resolution $res."

    cat << EOF > run_dumper.py
import hicstraw
import sys

result = hicstraw.straw("$readtype", "$norm", "$hicinnie", "$chrom", "$chrom", "BP", int("$res"))
for record in result:
    print("{0}\t{1}\t{2}".format(record.binX, record.binY, record.counts))

EOF

    python run_dumper.py $readtype $norm $hicpath $chrom $res | grep -v WARN > $dumpoutie 2> $dumperrors
}

process_genes_Bbins() {
    local myres=$1
    local mychr=$2
    local genesblock=$3
    local sortedbed=$4
    local pseudoB=$5

    if [[ "$myres" == "$maxres" && "$genesblock" -eq 1 ]]; then

        echo "Using initial Genes and Bbins files for $maxres"

        maxres_processed=1

        cat $EVAstates | mawk -v myres=$maxres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$maxres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed

        cat $EVBstates |  mawk -v myres=$maxres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$maxres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $pseudoB

    elif [ "$switch" -lt 1 ]; then

        cat $EVAstates | mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed

        cat $EVBstates | mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed

    else  

        echo "Using combined bins for $myres"

        cat $combined_newAbins | mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed

        cat $combined_newBbins |  mawk -v myres=$myres -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$myres '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $pseudoB

    fi
}


run_EigenVector() {
    local hic_file="$1"
    local chrom="$2"
    local res="$3"
    local n_components=10
    local eigendumpoutie="Eigen_dump_${chrom}.txt"
    local eigendumperrors="Eigen_errors_${chrom}.log"

    echo "Generating EV's from 1 to $n_components for $chrom at resolution $res."

    run_dumper oe VC_SQRT $hic_file $chrom $res $eigendumpoutie $eigendumperrors

    cat << EOF > run_EigenVector.py
import numpy as np
import pandas as pd
from scipy.linalg import eigh
import sys

def save_eigenvector_to_file(eigenvector, bin_starts, res, chrom, output_file):
    with open(output_file, 'w') as f:
        for bin_start, value in zip(bin_starts, eigenvector):
            f.write(f"{chrom}\\t{bin_start}\\t{bin_start + res}\\t{value}\\n")

def process_hic_to_eigen(dump_file, res, chrom, n_components=10):
    dump_df = pd.read_csv(dump_file, delim_whitespace=True, header=None, names=['binX', 'binY', 'counts'])
    dump_df = dump_df.apply(pd.to_numeric, errors='coerce').dropna()

    binX, binY, counts = dump_df['binX'].values, dump_df['binY'].values, dump_df['counts'].values
    n_bins = max(max(binX), max(binY)) // res + 1
    matrix = np.zeros((n_bins, n_bins))

    for x, y, c in zip(binX, binY, counts):
        if x // res < n_bins and y // res < n_bins:
            matrix[x // res, y // res] = c
            if x != y:
                matrix[y // res, x // res] = c

    mask = np.sum(matrix, axis=0) > 0
    matrix = matrix[mask][:, mask]
    pearson_matrix = np.corrcoef(matrix)

    eigenvalues, eigenvectors = eigh(pearson_matrix)
    idx = eigenvalues.argsort()[::-1]
    eigenvectors = eigenvectors[:, idx]
    eigenvectors = eigenvectors[:, :n_components]

    for i in range(n_components):
        save_eigenvector_to_file(eigenvectors[:, i], np.arange(0, n_bins * res, res)[mask], res, chrom, f"PC{i+1}_{chrom}_{res}.bedgraph")

if __name__ == "__main__":
    dump_file = sys.argv[1]
    res = int(sys.argv[2])
    chrom = sys.argv[3]

    process_hic_to_eigen(dump_file, res, chrom)
EOF

    python run_EigenVector.py "$eigendumpoutie" "$res" "$chrom"
}

check_and_flip_sign() {
    local pc_file="$1"
    local genes_file="$2"
    local output_file="sign_check_output.txt"

    echo "Checking and potentially flipping signs for $pc_file based on gene content from $genes_file..."

    pos_neg=$(awk 'BEGIN {FS=OFS="\t"} {if ($4 > 0) print $0, "positive"; else if ($4 < 0) print $0, "negative"}' "$pc_file" |
        sed 's/chr//g' |
        bedtools intersect -u -a stdin -b "$genes_file" |
        sort -k1,1 -V -k5b,5b |
        bedtools groupby -i stdin -g 1,5 -c 1 -o count |
        bedtools groupby -i stdin -g 1 -c 3 -o collapse |
        sed 's/,/\t/g' |
        awk 'BEGIN {FS=OFS="\t"} {if ($2 > $3) print "-"; else print "+"}'
    )

    echo "Overlap check result for $pc_file: $pos_neg" >> "$output_file"

    if [ "$pos_neg" == "-" ]; then
        echo "Flipping signs for $pc_file..." >> "$output_file"
        awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3, -1 * $4}' "$pc_file" > temp && mv temp "$pc_file"
    else
        echo "No flipping needed for $pc_file." >> "$output_file"
    fi

    echo "Check and flip sign process completed for $pc_file. See $output_file for details."
}

process_all_pcs_for_chromosome() {
    local chrom="$1"
    local res="$2"
    local genes_file="$3"

    echo "Starting sign check and flip for all PCs for chromosome $chrom at resolution $res." > "sign_check_output.txt"

    for pc_file in PC*_${chrom}_${res}.bedgraph; do
        echo "Processing $pc_file..."
        check_and_flip_sign "$pc_file" "$genes_file"
    done

    echo "All PCs for chromosome $chrom at resolution $res have been processed. See sign_check_output.txt for details."
}

find_best_pc() {
    local chrom="$1"
    local res="$2"
    local genes_file="$3"
    local gc_file="$4"
    local chrom_size="$5"
    local best_pc=""
    local best_corr=-100
    local log_file="correlation_log_${chrom}_${res}.txt"

    echo "Calculating best PC for $chrom at resolution $res..." > "$log_file"

    local bin_size=$res
    local num_bins=$((chrom_size / bin_size))

    calculate_density() {
        local feature_file="$1"
        local chrom="$2"
        local bin_size="$3"
        local chrom_size="$4"

        declare -a density
        for ((i=0; i<num_bins; i++)); do
            density[$i]=0
        done

        awk -v bin_size="$bin_size" -v chrom="$chrom" '
        BEGIN {
            OFS = "\t"
        }
        $1 == chrom {
            bin_start = int($2 / bin_size)
            bin_end = int($3 / bin_size)
            for (i = bin_start; i <= bin_end; i++) {
                density[i]++
            }
        }
        END {
            for (i in density) print density[i]
        }
        ' "$feature_file"
    }

    gene_density=($(calculate_density "$genes_file" "$chrom" "$bin_size" "$chrom_size"))
    gc_density=($(calculate_density "$gc_file" "$chrom" "$bin_size" "$chrom_size"))

    for pc_file in PC*_${chrom}_${res}.bedgraph; do
        echo "Processing $pc_file..." >> "$log_file"

        local pc_values=($(awk -v bin_size="$bin_size" -v num_bins="$num_bins" '
        BEGIN {
            for (i = 0; i < num_bins; i++) pc[i] = 0
        }
        {
            bin = int($2 / bin_size)
            if (bin < num_bins) {
                pc[bin] = $4
            }
        }
        END {
            for (i = 0; i < num_bins; i++) print pc[i]
        }
        ' "$pc_file"))

        gene_density_str=$(IFS=,; echo "${gene_density[*]}")
        gc_density_str=$(IFS=,; echo "${gc_density[*]}")
        pc_values_str=$(IFS=,; echo "${pc_values[*]}")

        read -r gene_corr gc_corr <<< $(python3 - <<EOF
import numpy as np
from scipy.stats import pearsonr

pc_values = np.array([${pc_values_str}])
gene_density = np.array([${gene_density_str}])
gc_density = np.array([${gc_density_str}])

gene_corr, _ = pearsonr(pc_values, gene_density)
gc_corr, _ = pearsonr(pc_values, gc_density)

print(gene_corr, gc_corr)
EOF
)

        local total_corr=$(echo "$gene_corr - $gc_corr" | bc -l)

        echo "PC file: $pc_file, Gene correlation: $gene_corr, GC correlation: $gc_corr, Total correlation: $total_corr" >> "$log_file"

        if (( $(echo "$total_corr > $best_corr" | bc -l) )); then
            best_pc="$pc_file"
            best_corr="$total_corr"
            echo "New best PC: $best_pc with total correlation $best_corr" >> "$log_file"
        fi
    done

    echo "Best PC for $chrom at resolution $res: $best_pc with total correlation $best_corr" >> "$log_file"
    echo "Correlation calculation completed. See $log_file for details."

    cp "$best_pc" "best_Eigen_${chrom}_${res}.bedgraph"
}

generateEV() {
    local hic_file="$1"
    local chrom="$2"
    local res="$3"
    local genes_file="$4"
    local gc_file="$5"
    local chrom_size="$6"

    echo "Starting EV generation for chromosome $chrom at resolution $res..."
    
    run_EigenVector "$hic_file" "$chrom" "$eigenres"

    echo "PCs generated for $chrom at resolution $eigenres."

    echo "Flipping the signs for $chrom PC's at resolution $eigenres..."

    process_all_pcs_for_chromosome "$chrom" "$eigenres" "$genes_file"

    echo "Finding best PC for $chrom at resolution $eigenres..."

    find_best_pc "$chrom" "$eigenres" "$genes_file" "$gc_file" "$chrom_size"

    echo "Best PC identification completed for $chrom at resolution $eigenres."
}

concatenate_best_pcs() {
    echo "Concatenating the best eigenvectors for all chromosomes into a full genome file."
    cat best_Eigen_*_{$res}.bedgraph > EV_full_genome.bedgraph
    echo "Concatenation completed. Full genome eigenvector file: EV_full_genome.bedgraph"
}

run_smoothing() {
    local innieS=$1
    local outieS=$2

    cat "$innieS" | sort -k 1,1 -V -k 2bn,2b --stable > "$outieS"

    cat "$outieS" | tail -n +2 > tmp1_shiftedright

    cat "$outieS" | awk '{if (NR == 1) print $0"\n"$0; else print $0}' > tmp2_shiftedleft

    paste "$outieS" tmp1_shiftedright tmp2_shiftedleft | awk '{print $1"\t"$2"\t"$3"\t"(($4 + $8 + $12)/3)}' > shifted_outie

    mv shifted_outie "$innieS"
}

process_oppocheck_statement() {
    local innie_genes="$1"
    local innie_gc="$2"
    local outie_oppo="$3"
    local chr="$4"
    local res="$5"
    local input_dir="$6"
    local output_dir="$7"

    local zfile="${input_dir}/zerofile"
    local rowZ="${input_dir}/RowZs"

    if [ ! -f "$zfile" ] || [ ! -f "$rowZ" ]; then
        echo "Required file missing: $zfile or $rowZ"
        return 1
    fi

    mawk 'NR==FNR { c1[$1] = $2; next } { if ($1 in c1) print $2 "\t" $3; next } { print $0 }' "$innie_genes" "$rowZ" "$zfile" | awk '{ b[$1] += $2; c[$1]++; d[$1] += ($2**2) } END { for (i in b) { print i "\t" (b[i])*(c[i]) "\t" c[i] "\t" b[i]/c[$1] "\t" ((d[i]/(c[$1])-((b[i]/(c[$1]))**2)))**.5 } }' > "${output_dir}/A2removedfile" &

    mawk 'NR==FNR { c1[$1] = $2; next } { if ($1 in c1) print $2 "\t" $3; next } { print $0 }' "$innie_gc" "$rowZ" "$zfile" | awk '{ b[$1] += $2; c[$1]++; d[$1] += ($2**2) } END { for (i in b) { print i "\t" (b[i])*(c[i]) "\t" c[i] "\t" b[i]/c[$1] "\t" ((d[i]/(c[$1])-((b[i]/(c[$1]))**2)))**.5 } }' > "${output_dir}/B2removedfile" &

    wait

    myABsum=$(cat "${output_dir}/A2removedfile" "${output_dir}/B2removedfile" | mawk '{ sum += $3 } END { print sum }')

    wait

    if [ "$verbose" -gt 0 ]; then
        end=$(date +%s)
        seconds=$(echo "$end - $start" | bc)
        echo "It took $seconds seconds. Comparing A vs B scores for $chr chromosome."
        start=$(date +%s)
    fi

    mawk -v var="$res" -v myABsum="$myABsum" '{ if ($1 in b) b[$1] -= $2; else b[$1] = $2 } END { for (i in b) { if (b[i] != 0) print i "\t" (b[i])/myABsum } }' "${output_dir}/A2removedfile" "${output_dir}/B2removedfile" | sort -k 1bn,1b --stable | mawk -v mchr="$chr" -v myres="$res" -v window="$newwindow" -v OFS="\t" 'BEGIN { slide = 1 } { mod = NR % window; if (NR <= window) { print mchr, int($1), $2; count++ } else { sum -= array[mod] } sum += $2; array[mod] = $2;} (NR % slide) == 0 { print mchr, int($1-((myres*window)/2)), sum/count }' | mawk -v myres="$res" '{ print $1 "\t" $2-int(myres/2)+myres "\t" $2+(int(myres/2))+myres "\t" $3 }' | mawk '{ if ($2 > 0) print $0 }' > "$outie_oppo" 2> "${output_dir}/errorfile"

    if [ "$exclbed" != 0 ]; then
        exclbins="${output_dir}/exclbins"
        excloutie="${output_dir}/exclout_${res}_${chr}"

        cat "$exclbed" | mawk -v myres="$res" -v mychr="$chr" '{ if ($1 == mychr) print $1 "\t" int($2/myres)*myres "\t" int($3/myres)*myres }' | awk -v myres="$res" '{ for (i = $2; i <= $3; i += myres) b[i] += 1 } END { for (j in b) print j "\t" b[j] }' > "$exclbins"

        wait

        mawk 'NR==FNR { c1[$1] = $2; next } { if ($2 in c1); else print $0 }' "$exclbins" "$outie_oppo" > "$excloutie"

        wait

        mv "$excloutie" "$outie_oppo"
    fi

    if [[ "$oppocheck" -gt 0 && "$adjustment" -gt 0 ]]; then

        amed=`mawk 'NR==FNR { c1[$1] = $2; next} {if ($2 in c1) print $4}' $innie_genes $outie_oppo | mawk '{sum +=$1} END {print sum/NR}'`
        bmed=`mawk 'NR==FNR { c1[$1] = $2; next} {if ($2 in c1); else print $4}' $innie_genes $outie_oppo | mawk '{sum += $1} END {print sum/NR}'`

        wait

        cat $outie_oppo | mawk -v varA=$amed -v varB=$bmed '{print $1"\t"$2"\t"$3"\t"$4-(varA-(varB*-1))}' > $outie2

        wait

        mv $outie2 $outie_oppo

    fi

    if [ "$myres" -ge "$endZ" ]; then

        myendmean=`cat $outie_oppo | mawk '{sum += $4; cnum++} END {print sum/cnum}'`
        echo "$myendmean"
        tmpendZ="${output_dir}/tmpendnorm"
        cat $outie_oppo | mawk -v mymean=$myendmean '{print $1"\t"$2"\t"$3"\t"($4-mymean)}' > $tmpendZ

        wait

        mv $tmpendZ $outie_oppo
    fi
}

run_midshifter() {
    local finer_res=$1
    local coarser_res=$2
    local high_innie=$3
    local coarse_innie=$4
    local shifter_outie=$5

    echo "Running shifter with Finer Resolution: $finer_res and Coarser Resolution: $coarser_res"
    echo "High-res file: $high_innie, Coarse-res file: $coarse_innie"
    
    cat << EOF > hdrescompare.py
import numpy as np
from tqdm import tqdm

sizefile = "$sizefile"
res = int("$finer_res")
coarseres = int("$coarser_res")
infile1 = "$high_innie"
infile2 = "$coarse_innie"
outiewrite = "$shifter_outie"

chrnames = []
chrsizes = []
with open(sizefile, 'r') as innie:
    for line in innie:
        sline = line.strip()
        lines = sline.split("\t")
        chrnames.append(lines[0])
        chrsizes.append(int(lines[1]))

chrmat = [[] for _ in chrnames]
smat = [[] for _ in chrnames]
emat = [[] for _ in chrnames]
scoremat = [[] for _ in chrnames]
coscoremat = [[] for _ in chrnames]

for i in tqdm(range(len(chrnames))):
    for j in range(0, chrsizes[i] + coarseres + coarseres, res):
        chrmat[i].append(chrnames[i])
        smat[i].append(j)
        emat[i].append(j + res)
        scoremat[i].append(0)
        coscoremat[i].append(0)

with open(infile1, 'r') as innie:
    for line in innie:
        sline = line.strip()
        lines = sline.split("\t")
        try:
            myindex = chrnames.index(lines[0])
            mypos = int(float(lines[1])/res)
            scoremat[myindex][mypos] = float(lines[3])
        except (ValueError, IndexError):
            continue

window_size = 3 * (coarseres // res)

rolling_avg_scores = [[] for _ in chrnames]
for c in range(len(scoremat)):
    for i in range(len(scoremat[c])):
        window_start = max(0, i - window_size // 2)
        window_end = min(len(scoremat[c]), i + window_size // 2)
        window_scores = scoremat[c][window_start:window_end]
        rolling_avg_scores[c].append(np.mean(window_scores))

shifted_scores = [[0 for _ in range(len(scoremat[c]))] for c in range(len(scoremat))]

for c in range(len(scoremat)):
    for i in range(len(scoremat[c])):
        mymean = rolling_avg_scores[c][i]
        coarscore = coscoremat[c][i] if i < len(coscoremat[c]) else 0
        shiftval = mymean - coarscore
        newval = scoremat[c][i] - shiftval
        shifted_scores[c][i] = newval

with open(outiewrite, "w") as ofile:
    for w in range(len(shifted_scores)):
        for p in range(len(shifted_scores[w])):
            ofile.write(f"{chrmat[w][p]}\t{smat[w][p]}\t{emat[w][p]}\t{shifted_scores[w][p]}\n")
EOF

    python hdrescompare.py

    echo "After Shifter: Executed with Coarser Resolution: $coarser_res, Finer Resolution: $finer_res"
    echo "Shifter executed and output written to $shifter_outie"

}

run_endshifter() {
    local finer_res=$1
    local coarser_res=$2
    local high_innie=$3
    local coarse_innie=$4
    local shifter_outie=$5
    
    echo "Running shifter with Finer Resolution: $finer_res and Coarser Resolution: $coarser_res"
    echo "High-res file: $high_innie, Coarse-res file: $coarse_innie"

    awk -v var="$coarser_res" '
        {
            for (i=0; i<=$2; i+=var) {
                start = i-(var*2)
                end = i+(var*3)
                if (start < 0) start = 0
                if (end > $2) end = $2
                print $1 "\t" start "\t" end
            }
        }
    ' "$sizefile" > temp_intervals.bed
    
    intersectBed -wa -a temp_intervals.bed -wb -b $high_innie | groupBy -i stdin -g 1,2,3 -c 7 -o mean | intersectBed -wa -a $coarse_innie -wb -b stdin | awk '{print $1"\t"$2"\t"$3"\t"$8-$4}' | intersectBed -wa -a $high_innie -wb -b stdin | groupBy -i stdin -g 1,2,3,4 -c 8 -o mean | awk '{print $1"\t"$2"\t"$3"\t"$4-$5}' > $shifter_outie

    rm temp_intervals.bed

    echo "Shifter executed and output written to $shifter_outie"
}


separating_ABbins() {
    local separatinginnie=$1

AbinsR=`echo "AbinsfromGI_""$myres"".txt"`
BbinsR=`echo "BbinsfromGI_""$myres"".txt"`

newBbins=`echo "newBbinsGI_""$myres"".txt"`
newAbins=`echo "newAbinsGI_""$myres"".txt"`

non_newBbins=`echo "non_newBbinsGI_""$myres"".txt"`
non_newAbins=`echo "non_newAbinsGI_""$myres"".txt"`

combined_newBbins=`echo "combined_newBbinsGI_""$myres"".txt"`
combined_newAbins=`echo "combined_newAbinsGI_""$myres"".txt"`

cat $separatinginnie | awk '{if ($4 > 0) print $0}' > $AbinsR

cat $separatinginnie | awk '{if ($4 < 0) print $0}' > $BbinsR

if [ "$switch" -eq 1 ]; then

    echo "Re-initialization of Bins is switched ""ON"". Re-evaluating both A and B bins and re-initializing for next resolution."

    cat $Bbins | intersectBed -u -a stdin -b $BbinsR > $newBbins

    awk 'NR==FNR{a[$1]=$1;next}!a[$1]' $newBbins $Bbins > $non_newBbins

    cat $newBbins $non_newBbins > $combined_newBbins 

    cat $genesfile | cut -f 1-3 | intersectBed -u -a stdin -b $AbinsR > $newAbins

    awk 'NR==FNR{a[$1]=$1;next}!a[$1]' $newAbins $genesfile > $non_newAbins

    cat $newAbins $non_newAbins > $combined_newAbins 

else

    echo "Generating the output without re-initialization."

fi

}

BHcorrection() {
    local BHinnie=$1
    echo "Generating Benjamini-Hochberg Corrected Values"

cat << EOF > BHcorrection.py
from statsmodels.stats.multitest import multipletests
inputfile="$BHinnie"
outfile="bhcorrected"
mystarts=[]
mychr=[]
myends=[]
mypvals=[]
with open(inputfile, 'r') as innie:
        for line in innie:
                sline = line.strip()
                lines = sline.split("\t")
                mychr.append(str(lines[0]))
                mystarts.append(str(lines[1]))
                myends.append(str(lines[2]))
                mypvals.append(float(lines[3]))
hyptest, FDRlist, _, alpha_corrected =multipletests(mypvals, alpha=0.5, method='fdr_bh', is_sorted=False, returnsorted=False)
with open(outfile, 'w') as outie:
        for w in range(0,len(FDRlist)):
                outie.write(str(mychr[w]) + "\t" + str(mystarts[w]) + "\t" + str(myends[w]) + "\t" + str(FDRlist[w]) + "\n")

EOF

python BHcorrection.py

}

calculate_pvalues() {
    local output_dir=$1
    local chr=$2
    local res=$3
    local outie_pvalues=$4

    awk 'NR==FNR { c[$1] = $3; m[$1] = $4; sd[$1] = $5; next } { if (($1 in m) && ($3 != 1) && (c[$1] != 1)) print $1 "\t" (m[$1]-$4) / (((sd[$1]**2)/c[$1]) + (($5**2)/$3))**.5 "\t" (c[$1] + $3 -2) }' "${output_dir}/A2removedfile" "${output_dir}/B2removedfile" > "${output_dir}/ttable"

    zfileforpy="${output_dir}/ttable"
    pfileforpy="${output_dir}/ptable"

    cat << EOF > zlookup.py
import scipy.stats
inputfile = "$zfileforpy"
outfile = "$pfileforpy"
mybins = []
mypvals = []
with open(inputfile, 'r') as innie:
    for line in innie:
        sline = line.strip()
        lines = sline.split("\t")
        mybins.append(str(lines[0]))
        mytstat = float(lines[1])
        mydf = int(lines[2])
        mypvals.append(scipy.stats.norm.sf(abs(mytstat)))
with open(outfile, 'w') as outie:
    for w in range(len(mybins)):
        outie.write(str(mybins[w]) + "\t" + str(mypvals[w]) + "\n")
EOF
    python zlookup.py

    mawk -v mchr="$chr" -v myres="$res" '{ print mchr "\t" $1 "\t" $1+myres "\t" $2 }' "${output_dir}/ptable" > "$outie_pvalues"
    wait
}

declare -a maxres_processed=0

process_resolution() {
    local myres=$1
    local genesblock=$2

    if [ $window == 0 ]; then
        totbins=$(awk -v var=$myres '{print $1"\t"($2/var)**2}' $sizefile | awk '{sum += $2} END {print sum}')
        newwindow=1
    else
        newwindow=$(mawk -v win=$window -v res=$myres 'BEGIN{print int(win/(res/1000))+1}')
    fi

    outiefull="Crush_${myres}.bedgraph"

    echo -e "\nNow dumping reads and calculating initial CRUSH for all chromosomes at $myres"

    task() {
        mychr=$(awk -v var=$myiter 'NR==var {print $1}' $sizefile)
        mysize=$(awk -v var=$myiter 'NR==var {print $2}' $sizefile)
        rowbins=$(awk -v myres=$myres -v var=$myiter 'NR==var {print int($2/myres)+1}' $sizefile)
        
        outie="Crush_${myres}_${mychr}_tmp"
        outie2="Crush_${myres}_${mychr}_tmp2"
        tmpfiles="ABtmpfiles_${mychr}_${myres}"

        [ ! -d "$tmpfiles" ] && mkdir "$tmpfiles"

        mydist=$((distance + 0))
        myupper=$((upperlim + 0))

        if [ "$verbose" -gt 0 ]; then
            echo "Getting genic bins."
        fi
        
        sortedbed="${tmpfiles}/genebins"
        pseudoB="${tmpfiles}/pseudoB"
        dumped="${tmpfiles}/dumped_${myres}_${mychr}_tmp"

        process_genes_Bbins $myres $mychr $genesblock $sortedbed $pseudoB

        if [ "$verbose" -gt 0 ]; then
            echo "Dumping ${mychr} chromosomal reads."
        fi

        newofile="${tmpfiles}/observeds"
        newefile="${tmpfiles}/expecteds"
        dumperrors="dumperrors_${myres}"

        if [ $juiceorcool -eq 0 ]; then
            run_dumper 'observed' $norm $hicpath $mychr $myres $dumped $dumperrors
        else

            echo "Dumping using cooler..."

            cooler dump $hicpath::/resolutions/$myres -r $mychr --join | cut -f 2,5,7 > $dumped 2> $dumperrors
        fi

        chromcheck=$(grep -e KeyError -e name dumperrors_$myres | wc -l)
        if [ $chromcheck -gt 0 ]; then
            echo "WARNING: A chromosome in your size file was not found in the .mcool file. Make sure the names match up. For example, check for chr1 vs. 1."
        fi

        wait

        if [ "$upperlim" -gt 0 ]; then
            awk -v var=$myres -v distfilt=$mydist -v upper=$upperlim '{if (($2-$1 >= distfilt) && ($2-$1 <= upper)) print ($2-$1)/var"\t"$1"\t"$2"\t"$3}' $dumped | \
            awk -v var=$myres -v mychrsize=$mysize -v newofile=$newofile -v newefile=$newefile '{b[$1]+=$4; print $0 >> newofile } END { for (i in b) print i"\t"b[i]/((mychrsize/var)-i) >> newefile }' &
        else
            awk -v var=$myres -v distfilt=$mydist '{if ($2-$1 >= distfilt) print ($2-$1)/var"\t"$1"\t"$2"\t"$3}' $dumped | \
            awk -v var=$myres -v mychrsize=$mysize -v newofile=$newofile -v newefile=$newefile '{b[$1]+=$4; print $0 >> newofile } END { for (i in b) print i"\t"b[i]/((mychrsize/var)-i) >> newefile }' &
        fi

        wait 

        if [ "$verbose" -gt 0 ]; then
            echo "Done dumping."
            echo "Now getting distance normalized values for ${mychr} chromosome."
            start=$(date +%s)
        fi

        errorfile="errorlog"
        tmpmean="${tmpfiles}/tmpmean"
        symmfile="${tmpfiles}/symmonfile"

        mawk 'NR==FNR { c1[$1] = $2; next} {if ($2 != $3) print $2"\t"$3"\t"($4+1)/(c1[$1]+1)}' $newefile $newofile | mawk -v thresh=$threshold '{if ($3 < thresh) print $1"\t"$2"\t"$3"\n"$2"\t"$1"\t"$3; else print $1"\t"$2"\t"thresh"\n"$2"\t"$1"\t"thresh}' | awk -v symmfile=$symmfile -v tmpmean=$tmpmean -v rowbins=$rowbins '{m[$1]+=$3;d[$1]+=($3**2); print $0 >> symmfile} END {for (i in m) print i"\t"m[i]/rowbins"\t"((d[i]/rowbins-((m[i]/rowbins)**2))**.5) >> tmpmean }'  2> $tmpfiles/$errorfile

        wait

        if [ "$verbose" -gt 0 ]; then
            end=$(date +%s)
            seconds=$(echo "$end - $start" | bc)
            start=$(date +%s)
        fi

        if [ "$myres" -gt 50000 ]; then
            oppocheck=$(mawk 'NR==FNR { c1[$1] = $2; next} {if ($1 in c1) ; else if ($2 in c1); else print $0}' $sortedbed $symmfile | wc -l)
        else
            oppocheck=2
        fi

        if [ "$oppocheck" -lt 1 ] && [ "$verbose" -gt 0 ]; then
            echo "WARNING: All bins on chromosome ${mychr} overlap with your bed file at ${myres} resolution. This should not affect the scores if using multiple resolutions. If using a single resolution, then maybe not use this resolution?"
        fi

        if [ "$verbose" -gt 0 ]; then
            end=$(date +%s)
            seconds=$(echo "$end - $start" | bc)
            start=$(date +%s)
        fi

        RowZs=`echo "RowZs"`
        mawk 'NR==FNR { m[$1] = $2; d[$1] = $3; next} {if (d[$2] > 0) print $1 "\t" $2 "\t" ($3 - m[$2]) / (d[$2])}' $tmpmean $symmfile > $tmpfiles/$RowZs

        zfile="zerofile"
        mawk -v var=$mychr '{if ($1 == var) print $0}' $sizefile | mawk -v wlkr=$myres '{for (i=0;i<=$2;i+=wlkr) print i"\t"i"\t0.1"}' > $tmpfiles/$zfile

        process_oppocheck_statement $sortedbed $pseudoB $outie $mychr $myres "ABtmpfiles_${mychr}_${myres}" "ABtmpfiles_${mychr}_${myres}"

    }

    open_sem() {
        mkfifo pipe-$$
        exec 3<>pipe-$$
        rm pipe-$$
        local i=$1
        for((;i>0;i--)); do
            printf %s 000 >&3
        done
    }

    run_with_lock() {
        local x
        read -u 3 -n 3 x && ((0==x)) || exit $x
        (
        ( "$@"; )
        
        printf '%.3d' $? >&3
        )& spinner $!
    }

    open_sem $cpu 

    if [ "$whichchoice" == "o" ]; then
        for (( myiter=1; myiter<=$totchroms; myiter++ )); do
            run_with_lock task $myiter
        done
    fi

    wait

    echo "Resolution: $myres for Resolution output"
    echo -ne "Merging individual chromosome files\033[0K\r"

    cat Crush_${myres}_*_tmp | grep -v -i nan | \
    awk '
    {
        key = $1"\t"$2"\t"$3;
        sum[key] += $4;
        count[key]++;
    } 
    END {
        for (key in sum) {
            if (count[key] > 1) {
                avg = sum[key] / count[key];
                print key "\t" avg;
            } else {
                print key "\t" sum[key];
            }
        }
    }' | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $outiefull

    wait 

    rm Crush_${myres}_*_tmp


    if [ "$myres" == "$maxres" ]; then
        coarse_infile=`echo "GI_""$myres"".bedgraph"`
        cat $outiefull > $coarse_infile
        wait
    fi

    res_values=`echo "res_values"`
    echo "$myres" >> $res_values


    if [ "$countres" -gt 0 ]; then
        python_output=`echo "shifter.bedgraph"`

        if [ "$coarse_res" -eq "$maxres" ] && [ "$coarse_res" != "$high_res" ]; then

            high_res_infile=`echo "$outiefull"`
            coarse_res_infile=`echo "$coarse_infile"`

        else
            
            high_res_infile=`echo "$outiefull"` 
            coarse_res_infile=`echo "$python_output"`

        fi

        if [[ "$(( $coarse_res % $high_res ))" -eq 0 ]]; then

            crtmp=$coarse_res

        else

            if [[ "$(( $crtmp % $high_res ))"  -eq 0 ]]; then
                coarse_res=$crtmp
            fi

        fi

        run_midshifter $high_res $coarse_res $high_res_infile $coarse_res_infile $python_output

        echo "Resolution: "$myres" for shifter output"

        separating_ABbins $python_output

        coarse_res=$high_res

    else

        separating_ABbins $outiefull

    fi

    wait

    rm Crush*_tmp 2> smallerrors
    rm ACrush*_tmp 2> smallerrors
    rm BCrush*_tmp 2> smallerrors

    countres=$((countres+1))
}

reprocess_resolutions_with_shifter() {
    local current_res=$1
    local end_index=$2
    local reprocess=$3

    countres_reprocess=0
    coarse_res_reprocess=$maxres
    crtmp_reprocess=$maxres


    for (( j=0; j<=$end_index; j++ )); do

        local prev_res=${res_array[$j]}

        echo "Reprocessing for resolution: $prev_res"

        highres_reprocess=$prev_res

        local outiefull_reprocess=`echo "Crush_Re_""$prev_res"".bedgraph"`

        local chromosomes=($(cut -f1 $sizefile))

        for mychr in "${chromosomes[@]}"; do
            local sortedbed_reprocess="genebins_reprocess_${mychr}_${prev_res}"
            local pseudoB_reprocess="pseudoB_reprocess_${mychr}_${prev_res}"
            local outie_reprocess="Crush_reprocess_${prev_res}_${mychr}_tmp"
            
            reprocesstmpfiles=`echo "ABreprocesstmpfiles_""$mychr""_""$prev_res"`

            [ ! -d "$reprocesstmpfiles" ] && mkdir "$reprocesstmpfiles"

            echo "Using combined bins for $prev_res and the current value of resolution in this loop is $current_res"

            cat $combined_newAbins | mawk -v myres=$prev_res -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$prev_res '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $sortedbed_reprocess

            cat $combined_newBbins |  mawk -v myres=$prev_res -v mychr=$mychr '{if ($1 == mychr) print $1"\t"int($2/myres)*myres"\t"int($3/myres)*myres}' | awk -v myres=$prev_res '{for (i=$2;i<=$3;i+=myres) b[i]+=1}  END { for (j in b) print j"\t"b[j]} ' > $pseudoB_reprocess

            echo "Files: $sortedbed_reprocess, $pseudoB_reprocess"

            process_oppocheck_statement $sortedbed_reprocess $pseudoB_reprocess $outie_reprocess $mychr $prev_res "ABtmpfiles_${mychr}_${prev_res}" "ABreprocesstmpfiles_${mychr}_${prev_res}"

            wait 

            if [ $prev_res -eq $minres ] && [ $pcalculation -gt 0 ]; then
                pscores_reprocess="ttest_reprocess_${prev_res}_${mychr}_tmp"
                calculate_pvalues "ABreprocesstmpfiles_${mychr}_${prev_res}" $mychr $prev_res $pscores_reprocess
            fi            
        done    

        echo -ne "Merging individual chromosome files\033[0K\r"

        cat Crush_reprocess_${prev_res}_*_tmp | grep -v -i nan | \
        awk '
        {
            key = $1"\t"$2"\t"$3;
            sum[key] += $4;
            count[key]++;
        } 
        END {
            for (key in sum) {
                if (count[key] > 1) {
                    avg = sum[key] / count[key];
                    print key "\t" avg;
                } else {
                    print key "\t" sum[key];
                }
            }
        }' | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $outiefull_reprocess
        
        wait

        rm Crush_reprocess_${prev_res}_*_tmp

        if [ $prev_res -eq $minres ] && [ $pcalculation -gt 0 ]; then
        
            outiefullPval_reprocess="pvalues_Re_${prev_res}.bedgraph"
            echo "Resolution: $prev_res for Resolution output and minres is : $minres"
            cat ttest_reprocess_${prev_res}_*_tmp | grep -v -i nan | sort -k 1,1 -V -k 2bn,2b -k 3bn,3b --stable > $outiefullPval_reprocess
            wait 
        fi

        if [ "$prev_res" == "$maxres" ]; then

            coarse_infile_reprocess=`echo "GI_Re_""$prev_res"".bedgraph"`
            cat $outiefull_reprocess > $coarse_infile_reprocess

            wait
        fi

        prev_values=`echo "prev_values"`
        echo "$prev_res" >> $prev_values


        echo "$j"

        if [ "$j" -gt 0 ]; then

            if [ "$coarse_res_reprocess" -eq "$maxres" ] && [ "$coarse_res_reprocess" != "$highres_reprocess" ]; then
                highres_infile_reprocess=`echo "$outiefull_reprocess"`
                coarse_res_infile_reprocess=`echo "$coarse_infile_reprocess"`
            else
                highres_infile_reprocess=`echo "$outiefull_reprocess"` 
                coarse_res_infile_reprocess=`echo "$python_output"`

            fi

            if [[ "$(( $coarse_res_reprocess % $highres_reprocess ))" -eq 0 ]]; then
                crtmp_reprocess=$coarse_res_reprocess
            else
                if [[ "$(( $crtmp_reprocess % $highres_reprocess ))" -eq 0 ]]; then
                coarse_res_reprocess=$crtmp_reprocess
                fi
            fi

            if [ -f "$highres_infile_reprocess" ] && [ -f "$coarse_res_infile_reprocess" ]; then
                echo "Files exist: $highres_infile_reprocess, $coarse_res_infile_reprocess"
            else
                echo "One or both files do not exist: $highres_infile_reprocess, $coarse_res_infile_reprocess"
            fi

            echo "Running shifter for $prev_res"

            run_midshifter $highres_reprocess $coarse_res_reprocess $highres_infile_reprocess $coarse_res_infile_reprocess $python_output

            echo "Shifter run completed for $prev_res"

            echo "Resolution: "$prev_res" for shifter output"

            separating_ABbins $python_output

            coarse_res_reprocess=$highres_reprocess
        else
            separating_ABbins $outiefull_reprocess
        fi
    done
}


echo "Creating and moving into temporary directory"
mkdir $crushdir
cd $crushdir

if [ $juiceorcool -eq 0 ]
then
cat << EOF > listres.py
import hicstraw
hic = hicstraw.HiCFile("$hicpath")
totres=hic.getResolutions()
keepres=[]
for i in totres:
	if (int(i) >= int("$res")) and (int(i) <= int("$coarsestres")):
		keepres.append(str(i))
print(",".join(keepres))

EOF

res=`python listres.py`

else
reslist=`cooler ls $hicpath | sed 's/\//\t/g' | awk -v var=$res -v cres=$coarsestres '{if (($NF >= var) && ($NF <= cres)) print $NF}' | sort -k 1bnr,1b --stable | awk '{if (NR == 1) printf "%s", $1; else printf ",%s", $1}' | awk '{print $0}'`
res=$reslist
fi


isfasta=`head -1 $fastafile | awk '{if ($1 ~ /^>/) print 1; else print 0}'`

echo "initial-A: $genesfile"
echo "sizes-file: $sizefile"

numcolumns=`mawk '{if (NF != 2) print $0}' $sizefile | wc -l`
if [ "$numcolumns" -gt 0 ]
then
echo "size file is incorrect format. It should be two columns with chromosome name and size"
exit 0
fi

totchroms=`cat $sizefile | wc -l`
maxres=`echo "$res" | sed "s/,/ /g" | mawk '{max=$1;for(i=2;i<=NF;i++){if($i > max) max = $i} print max}'`
minres=`echo "$res" | sed "s/,/ /g" | mawk '{min=$1;for(i=2;i<=NF;i++){if($i < min) min = $i} print min}'`

if [ $isfasta == "1" ]; then
    echo "Fasta-file: $fastafile"
else
    echo "initial-B: $fastafile"
fi

resbins=`echo "size_bins"".bed"`
gcfile=`echo "gc_bins"".txt"`
gc_g_ga=`echo "gc_g_ga_bins"".bed"`
Bbins=`echo "Bbins"".bed"`

if [ $minres -ge 500 ]; then
    cat $sizefile | awk -v myres=500 '{for (i=0;i<=$2;i+=myres) print $1"\t"i"\t"i+myres}' > $resbins
    echo "$resbins calculated at 500bp resolution"
else
    cat $sizefile | awk -v myres=$minres '{for (i=0;i<=$2;i+=myres) print $1"\t"i"\t"i+myres}' > $resbins
    echo "$resbins calculated at ""$minres"" resolution"
fi

wait

if [ $isfasta == "1" ]; then
    echo "Now, calculating the gc content:"
    bedtools nuc -fi $fastafile -bed $resbins > $gcfile
    wait

    cat $gcfile | grep -v user | cut -f 1-5 | awk '{print $0"\t"($4+$5)}' | awk '{if ($6 > 0) print $0}' | awk '{print $0"\t"($5/$6)}' > $gc_g_ga
    wait

    gc_thresh=`cat $gc_g_ga | awk '{s+=$7; ss+=$7*$7; linecount+=1} END{print m=s/linecount, sqrt(ss/linecount-m^2)}' | awk '{print $0}' | awk '{print $1 - 1*$2}'`

    if [ $verbose -gt 0 ]; then
        echo "gc threshold is ""$gc_thresh"
    fi

    cat $gc_g_ga | awk -v thresh=$gc_thresh '{if ($7 < thresh) print $0}' > $Bbins
    wait
    echo "Finished generating Bbins file for all chromosomes "

else
    cat $fastafile > $Bbins
fi

if [ $endZ == 0 ]; then
    endZ=`echo "$minres"`
fi


EVAstates="EVAstates.bed"
EVBstates="EVBstates.bed"

if [ -n "$eigenfile" ]; then
    echo "Using user-provided Eigen Vector file: $eigenfile"
    
    awk '{if ($4 > 0) print $0}' "$eigenfile" > "$EVAstates"
    if [ $? -ne 0 ]; then
        echo "Error creating EVAstates from eigenfile. Please check the awk command."
        exit 1
    fi

    awk '{if ($4 < 0) print $0}' "$eigenfile" > "$EVBstates"
    if [ $? -ne 0 ]; then
        echo "Error creating EVBstates from eigenfile. Please check the awk command."
        exit 1
    fi

    if [ ! -s "$EVAstates" ]; then
        echo "Error: EVAstates is empty or not created."
        exit 1
    fi

    if [ ! -s "$EVBstates" ]; then
        echo "Error: EVBstates is empty or not created."
        exit 1
    fi

    echo "Eigen Vector files created successfully: $EVAstates, $EVBstates"
else

    echo "No user-provided Eigen Vector file. Generating Eigen Vectors using run_EigenVector script."

    for chrom in $(cut -f 1 $sizefile); do
        
        genesbinned="genesbinned.bed"
        gcbinned="gcbinned.bed"

        cat $genesfile | mawk -v myres=$eigenres -v mychr=$chrom '{if ($1 == mychr) print $1 "\t" int($2/myres)*myres "\t" int($3/myres)*myres}' > $genesbinned
        cat $Bbins | mawk -v myres=$eigenres -v mychr=$chrom '{if ($1 == mychr) print $1 "\t" int($2/myres)*myres "\t" int($3/myres)*myres}' > $gcbinned

        generateEV "$hicpath" "$chrom" "$eigenres" "$genesbinned" "$gcbinned" "$sizefile"
    
    done

    concatenate_best_pcs
fi

countres=0
coarse_res=$maxres
crtmp=$maxres
IFS=',' read -r -a res_array <<< "$res"
maxres_processed=0

for (( i=0; i<${#res_array[@]}; i++ )); do
    myres=${res_array[$i]}
    echo "Current resolution: $myres"
    echo "Array content: ${res_array[@]}"
    high_res=$myres

    genesblock=0
    if [ "$myres" == "$maxres" ] && [ "$maxres_processed" -eq 0 ]; then
        genesblock=1
    fi

    echo "Processing with myres = $myres, genesblock = $genesblock"

    process_resolution $myres $genesblock
 
    if [ $i -gt 0 ]; then
        reprocess_resolutions_with_shifter $myres $i 1


        finaloutie=`echo "$outpre""mergedCrush_""$myres"".bedgraph"`     
        Crush_todelete=`echo "Crush_todelete_""$myres""_reprocess"`

        finalpoutie=`echo "$outpre""mergedqvalue_""$myres"".bedgraph"`
        pval_todelete=`echo "pval_todelete_""$myres""_reprocess"`

        cat Crush_Re*.bedgraph | mawk -v minr=$myres '{print $1"\t"$2/minr"\t"$3/minr"\t"$4/int($3-$2)}' | mawk -v mr=$myres '{for (i=$2;i<$3;i++) print $1":"int(i)*mr":"(int(i)+1)*mr"\t"$4}' | mawk '{c1[$1] += $2; c2[$1]++} END {for (i in c1) print i"\t"c1[i]/c2[i]}' | sed 's/:/\t/g' > $finaloutie 2> smallerrors
        
        if [ $myres -eq $minres ]; then
            cat pvalues_Re*.bedgraph | mawk -v minr=$myres '{print $1"\t"$2/minr"\t"$3/minr"\t"$4}' | mawk -v mr=$myres '{for (i=$2;i<$3;i++) print $1":"int(i)*mr":"(int(i)+1)*mr"\t"$4}' | awk '{c1[$1] += log($2)/log(10); c2[$1]++} END {for (i in c1) print i"\t"10**(c1[i]/c2[i])}' | sed 's/:/\t/g' | sort -k 1,1 -V -k 2bn,2b --stable > $finalpoutie 2> smallerrors
        fi

        echo "Running BH correction"
        BHcorrection $finaloutie
        wait
        mv bhcorrected $finalpoutie
        wait

        GIave=`cat $finaloutie | mawk '{if ($4 < 0) sum+=($4*-1); else sum+=($4)} END {print sum/NR}' `

        if [ $(bc <<< "$qthresh > 0") -eq 1 ]; then

            finaloutiefilt=`echo "$outpre""mergedCrush_""$myres""_qfiltered_reprocess.bedgraph"`

            mawk -v fdr=$qthresh -v var=$GIave 'NR==FNR {a[$1":"$2":"$3] = $4; next} {if (a[$1":"$2":"$3] <= fdr) print $1"\t"$2"\t"$3"\t"$4/(var/100)}'  $finalpoutie $finaloutie> $finaloutiefilt

        fi

        if [ $trackline -eq 0 ]; then

            cat $finaloutie | mawk -v var=$GIave '{print $1"\t"$2"\t"$3"\t"$4/(var/100)}' > $Crush_todelete
            cat $finalpoutie | mawk '{print $1"\t"$2"\t"$3"\t"$4}' > $pval_todelete
        else
            cat $finaloutie | mawk -v var=$GIave '{print $1"\t"$2"\t"$3"\t"$4/(var/100)}' | mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=0,120,0 altColor=127,0,127 viewLimits=-20:20 autoScale off\n"$0; else print $0}' > $Crush_todelete
            cat $finalpoutie | mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=0,120,0 altColor=127,0,127 viewLimits=-20:20 autoScale off\n"$1"\t"$2"\t"$3"\t"$4; else print $1"\t"$2"\t"$3"\t"$4}' > $pval_todelete
            wait
        fi

        if [ $trackline -gt 0 ] && [ $(bc <<< "$qthresh > 0") -eq 1 ]; then
            filt_todelete=`echo "tmpcrushfiltered_""$myres""_reprocess"`
            cat $finaloutiefilt | mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=0,120,0 altColor=127,0,127 viewLimits=-20,20 autoScale off\n"$0; else print $0}' > $filt_todelete

            wait

            mv $filt_todelete $finaloutiefilt
            mv $Crush_todelete $finaloutie
            mv $pval_todelete $finalpoutie

            wait
        fi

        sizeBed=`echo "Chrsizes.bed"`
        finaloutie2=`echo "mergedCrush2_""$myres"".bedgraph"`

        cat $sizefile | awk '{print $1"\t""1""\t"$2}' > $sizeBed

        cat $finaloutie | grep -v track | intersectBed -wa -a stdin -wb -b $sizeBed | awk '{if ($3 <= $7) print $0}' | cut -f 1-4 | sort -k 1,1 -V -k 2bn,2b --stable | mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=204,0,0 altColor=0,0,0 viewLimits=-150:150 autoScale off\n"$0; else print $0}' > $finaloutie2
        mv $finaloutie2 $finaloutie
        
        if [ "$smoothing" -gt 0 ]; then
            echo "Smoothing Function is turned on"
            smoother="smoother.bedgraph"
            run_smoothing $newoutie $smoother
            rm $smoother tmp2_shiftedleft tmp1_shiftedright
        fi

        if [ $(bc <<< "$qthresh > 0") -eq 1 ] && [ $pcalculation -eq 1 ] && [ $myres -eq $minres ]; then 
            filt_todelete=`echo "tmpcrushfiltered_""$myres""_reprocess"`
            
            cat $finaloutiefilt | grep -v track | intersectBed -wa -a stdin -wb -b $sizeBed | \
                awk '{if ($3 <= $7) print $0}' | cut -f 1-4 | \
                mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=204,0,0 altColor=0,0,0 viewLimits=-150:150 autoScale off\n"$0; else print $0}' > $filt_todelete

            mv $filt_todelete $finaloutiefilt
            mv $finalpoutie ../
        fi

        mv $finaloutie ../
        
        if [ $(bc <<< "$qthresh > 0") -eq 1 ] && [ $pcalculation -eq 1 ] && [ $myres -eq $minres ]; then
            mv $finaloutiefilt ../
        fi
    fi
done

if [ $doNotMerge -eq 1 ]; then
    echo "Finished with each resolution listed, but not merging...."
    exit 0
fi

countres=0
for totres in $(echo $res | sed "s/,/\n/g" | sort -k 1bn,1b --stable | sed "s/\n/ /g" ); do
    echo "$totres" >> GIrestmpfile
    countres=$((countres+1))
done

if [ $countres -le 1 ]; then
    echo "Only one resolution listed, so nothing to merge. Check Crush_""$totres"".bedgraph for the single resolution GI score"

    rm AbinsfromGI_*.txt 2> smallerrors
    rm BbinsfromGI_*.txt 2> smallerrors
    rm newAbinsGI_*.txt 2> smallerrors
    rm newBbinsGI_*.txt 2> smallerrors
    rm non_newAbinsGI_*.txt 2> smallerrors
    rm non_newBbinsGI_*.txt 2> smallerrors
    rm combined_newAbinsGI_*.txt 2> smallerrors
    rm combined_newBbinsGI_*.txt 2> smallerrors
    rm smallerrors
    exit 0
fi

rm GIrestmpfile

wait

cd ../

if [ "$cleanup" -gt 0 ]; then
    echo "cleaning up"
    rm -r $crushdir
fi

if [ "$reshift" -gt 0 ]; then

    echo "Going into Final Shifter"

    for (( i=1; i<${#res_array[@]}; i++ )); do

        myres=${res_array[$i]}

        if [ $i -eq 1 ]; then
            echo "almost done..."
            oldres=`echo "$myres"`
            oldinnie=`echo "$outpre""mergedCrush_""$myres"".bedgraph"`
        else
            newres=`echo "$myres"`
            newinnie=`echo "$outpre""mergedCrush_""$myres"".bedgraph"`
            newoutie=`echo "mergedCrush2_""$myres"".bedgraph"`

            run_endshifter $newres $oldres $newinnie $oldinnie $newoutie

            if [ "$smoothing" -gt 0 ]; then
                echo "Smoothing Function is turned on"
                finalshifter="finalshifter.bedgraph"
                run_smoothing $newoutie $finalshifter
            fi 

            oldres=`echo "$myres"`
            oldinnie=`echo "mergedCrush2_""$myres"".bedgraph"`

            wait

            mawk '{if (NR == 1) print "track type=bedgraph visibility=full color=204,0,0 altColor=0,0,0 viewLimits=-150:150 autoScale off\n"$0; else print $0}' $newoutie > $newinnie    
            wait

        fi

    done
    
    if [ "$smoothing" -gt 0 ]; then
        rm $finalshifter
        rm tmp2_shiftedleft
        rm tmp1_shiftedright
    fi 

    rm mergedCrush2_*.bedgraph

fi

echo "Finished! Check the output."
echo -e "Note: Please keep in mind that we are using resolution walking.\nIf the coarsest resolution doesn't match the compartment pattern at that resolution, please consider re-running with the -m parameter set to start the walking with smaller bins."


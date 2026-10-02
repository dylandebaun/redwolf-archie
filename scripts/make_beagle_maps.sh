# Create a directory for the new Beagle-ready maps
mkdir -p beagle_maps
mkdir -p maps_sorted

for mapfile in canFam3.1/maps/*.txt; do
    #sort -k4,4n "$mapfile" | awk '($3 > prev) {print $0; prev=$3}' > maps_sorted/"$NAME"
    # Extract the chromosome name from the first line (skipping header)
    # or from the filename. If filename is chr10_average..., we grab 'chr10'
    RAW_CHR=$(echo "$mapfile" | grep -o 'chr[0-9]\+')
    # 2. Extract just the number, pad it to 2 digits, and prefix with 'chr'
    CHR=$(printf "chr%02d" "${RAW_CHR#chr}")
    echo "Processing $CHR..."
    sort -k2,2n "$mapfile" | \
awk 'NR>1 && ($4 > prev) {print chr "\t.\t" $4 "\t" $2; prev=$4}' chr="$CHR" \
> beagle_maps/"$CHR".beagle.map
    # Reformat: 
    # 1. Skip header (NR>1)
    # 2. Output: Chromosome, ".", MAP(cM), POS
done

mv beagle_maps/chr00.beagle.map beagle_maps/chrX.beagle.map
sed -i 's/chr00/chrX/g' beagle_maps/chrX.beagle.map 




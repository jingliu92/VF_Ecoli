# VF_Ecoli
## Find Virulence Genes in the whole dataset
```
conda activate virulencefinder_env
```

Run the pipeline 
```
for f in /home/jing/E.coli_test/ecoli_all/*.fasta; do
    sample=$(basename $f .fasta)
    echo "Running $sample..."

    mkdir -p /home/jing/E.coli_test/vf_out/$sample   # 👈 ADD THIS

    virulencefinder.py \
      -i $f \
      -o /home/jing/E.coli_test/vf_out/$sample \
      -p /home/jing/E.coli_test/virulencefinder_db \
      -x
done
```
Classify EHEC/STECT/NA based on the resuls
```
echo -e "Sample\tClass\tstx_subtype\teae" > classification_with_subtype.tsv

for f in /home/jing/E.coli_test/vf_out/*/results_tab.tsv; do
    sample=$(basename "$(dirname "$f")")

    stx_types=$(grep -oE 'stx[12][a-z]?' "$f" | sort -u | paste -sd "," -)
    has_stx=$(grep -E 'stx1|stx2' "$f")
    has_eae=$(grep -w 'eae' "$f")

    if [ -n "$has_stx" ] && [ -n "$has_eae" ]; then
        class="EHEC"
        eae_status="eae"
    elif [ -n "$has_stx" ]; then
        class="STEC"
        eae_status="NA"
    else
        class="Non-STEC"
        eae_status="NA"
    fi

    echo -e "$sample\t$class\t${stx_types:-NA}\t$eae_status"
done >> classification_with_subtype.tsv
```
Among 1,353 E.coli, 247 of them are classified as STEC, and 847 of them are EHEC and 254 of them are others.

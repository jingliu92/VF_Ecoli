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
## Classify EHEC/STECT/NA based on the resuls
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
Moving forward, we will extract stx2 from both EHEC and STEC sequences.

## Extract STEC + EHEC:
```
awk 'NR>1 && ($2=="STEC" || $2=="EHEC") {print $1}' classification_with_subtype.tsv > stec_samples.txt
```
## Extract stx2 sequences only from these
```
mkdir -p stx2_analysis

while read sample; do
    f="/home/jing/E.coli_test/vf_out/$sample/Hit_in_genome_seq.fsa"

    awk -v s="$sample" '
        /^>/ {
            keep = ($0 ~ /stx2/)
            if (keep) print ">" s "|" substr($0,2)
            next
        }
        keep { print }
    ' "$f"

done < stec_samples.txt > stx2_analysis/stx2_all.fasta
```
## Filter out the alignment that less than 300bp
```
awk '
/^>/ {
    if (seq && length(seq) >= 300) {
        print header; print seq
    }
    header=$0; seq=""; next
}
{ seq=seq $0 }
END {
    if (seq && length(seq) >= 300) {
        print header; print seq
    }
}
' stx2_analysis/stx2_all.fasta > stx2_analysis/stx2_filtered.fasta
```
No change. Meaning all stx2 genes extract are larger than 300bp. 2,063 hits are detected within 1,094 samples, means:
1) Multiple prophages (most common)

The genes for stx2 live on lambdoid bacteriophages. A single E. coli isolate can carry more than one phage, each with its own stx2.

Example in one genome:
stx2a (phage 1)
stx2c (phage 2)
2) Different subtypes in the same strain

You might see:

stx2a
stx2c

in one sample.
3) Assembly/contig effects

Because phages can be repetitive:

assemblies may split them across contigs
VirulenceFinder can report multiple hits for the same gene region

👉 You may see:

same subtype repeated
slightly different coordinates
## 

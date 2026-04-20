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
4,126 sequences were left. (Orginal: 45387). More than one stx2 genes detected in some of the isolates, and some are at eact locus, this could because database redundancy: stx2 is highly variable, so the database include ref sequences of slightly different sequences, from different strains,from different studies. Database = many known versions of the same gene
<img width="785" height="285" alt="image" src="https://github.com/user-attachments/assets/ed712348-5ede-4f89-8445-bb3dc02ff2e7" />


## deduplicate by locus
```
awk '
BEGIN {RS=">"; FS="\n"}
NR>1 {
    header=$1
    seq=""

    for(i=2;i<=NF;i++) seq=seq $i

    # extract unique locus key
    match(header, /(NODE_[^:]+:[0-9]+\.\.[0-9]+)/, c)
    locus=c[1]

    split(header, arr, "|")
    sample=arr[1]

    key = sample "|" locus

    if (!(key in seen)) {
        seen[key]=1
        print ">" header
        print seq
    }
}
' stx2_analysis/stx2_filtered.fasta > stx2_analysis/stx2_unique.fasta
```


## Extract useful info from original header
```
awk '
BEGIN {RS=">"; FS="\n"}
NR>1 {
    header=$1
    seq=""
    for(i=2;i<=NF;i++) seq=seq $i

    # ---- keep reference unchanged ----
    if (header ~ /^REF_stx2A/) {
        print ">" header
        print seq
        next
    }

    # ---- parse sample from header ----
    split(header, arr, "|")
    sample=arr[1]

    # ---- subtype (stx2, stx2a, stx2c, ...) ----
    subtype="stx2"
    if (match(header, /(stx2[a-z]?)/, m)) {
        subtype=m[1]
    }

    # ---- locus (contig + coordinates) ----
    # try NODE_... first; fallback to any contig:pos..pos pattern
    coord="NA"
    if (match(header, /(NODE_[^:]+:[0-9]+\.\.[0-9]+)/, c)) {
        coord=c[1]
    } else if (match(header, /([^|]+:[0-9]+\.\.[0-9]+)/, c2)) {
        coord=c2[1]
    }

    print ">" sample "|" subtype "|" coord
    print seq
}
' stx2_analysis/stx2_with_ref_aligned.fasta > stx2_analysis/stx2_with_ref_labeled.fasta

```
## Mafft Alignment
```
mafft --auto stx2_analysis/stx2_with_ref_labeled.fasta > stx2_analysis/stx2_with_ref_aligned.fasta
```
## Add STECT/EHEC in the alignment Header
```
awk '
BEGIN {
    RS=">"; FS="\n"

    # load sample → class mapping
    while ((getline < "sample_info.tsv") > 0) {
        class_map[$1] = $2
    }
}
NR>1 {
    header=$1
    seq=""

    for(i=2;i<=NF;i++) seq=seq $i

    # ---- keep reference unchanged ----
    if (header ~ /^REF_stx2A/) {
        print ">" header
        print seq
        next
    }

    # ---- extract sample ----
    split(header, arr, "|")
    sample=arr[1]

    # ---- get STEC/EHEC ----
    class = (sample in class_map) ? class_map[sample] : "NA"

    # ---- rebuild header ----
    # keep original info, just insert class after sample
    rest = substr(header, length(sample)+2)

    print ">" sample "|" class "|" rest
    print seq
}
' stx2_analysis/stx2_with_ref_aligned.fasta > stx2_analysis/stx2_final_annotated.fasta
```

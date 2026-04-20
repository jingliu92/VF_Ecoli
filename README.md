# Find Deletions in variable regions of stx2 A subunit
The goal of current project is to find deletions of stx2 in subunit A in STECs in our in-house sequenced E. coli database (1,353 isolates). 

## Use virulencefinder to find all the virulence genes in the whole dataset
```
conda activate virulencefinder_env
```

Run the pipeline 
```
for f in /home/jing/E.coli_test/ecoli_all/*.fasta; do
    sample=$(basename $f .fasta)
    echo "Running $sample..."

    mkdir -p /home/jing/E.coli_test/vf_out/$sample 
    virulencefinder.py \
      -i $f \
      -o /home/jing/E.coli_test/vf_out/$sample \
      -p /home/jing/E.coli_test/virulencefinder_db \
      -x
done
```
Note, 5 samples failed virulennce gene detection due to the short length of contigs. Left 1,348 isolates for downstream analysis.

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
2063 sequences were left. (Orginal: 2063). More than one stx2 genes detected in some of the isolates, and some are at eact locus, this could because database redundancy: stx2 is highly variable, so the database include ref sequences of slightly different sequences, from different strains,from different studies. Database = many known versions of the same gene
<img width="785" height="285" alt="image" src="https://github.com/user-attachments/assets/ed712348-5ede-4f89-8445-bb3dc02ff2e7" />


## Deduplicate by locus
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
```
 grep -c ">" stx2_unique.fasta
```
875 sequences left. * In the in-house dataset 838 isolates have stx2 *. So some of the isolates have more than 1 hits.

## Include stx2 A reference sequence
```
cat stx2A_ref.fasta stx2_unique.fasta > stx2_with_ref.fasta
```
## Extract useful info from original header Add STECT/EHEC in the Header
```
awk '
# ---- load mapping FIRST ----
FNR==NR {
    if (NR==1) next
    class_map[$1] = $2
    next
}

# ---- switch to FASTA mode AFTER loading ----
FNR==1 {
    RS=">"; FS="\n"
}

{
    if ($0 == "") next

    header=$1
    seq=""

    for(i=2;i<=NF;i++) seq=seq $i

    # ---- handle reference ----
    if (header ~ /^stx2A/ || header ~ /^REF/) {
        print ">REF_stx2A"
        print seq
        next
    }

    split(header, arr, "|")
    sample=arr[1]

    class = (sample in class_map) ? class_map[sample] : "NA"

    # subtype
    subtype="stx2"
    if (match(header, /(stx2[a-z]?)/, m)) {
        subtype=m[1]
    }

    # locus
    locus="NA"
    if (match(header, /(NODE_[^:]+:[0-9]+\.\.[0-9]+)/, c)) {
        locus=c[1]
    }

    print ">" sample "|" class "|" subtype "|" locus
    print seq
}
' sample_info.tsv stx2_analysis/stx2_with_ref.fasta > stx2_analysis/stx2_annotated.fasta

```

```
sed 's/\.contigs//g' stx2_analysis/stx2_annotated.fasta > stx2_analysis/stx2_final.fasta
```
## Seperate analysis for EHEC and STEC
```
awk '
BEGIN {RS=">"; FS="\n"}
NR>1 {
    header=$1
    if (header ~ /\|EHEC\|/) {
        print ">" $0
    }
}
' stx2_analysis/stx2_final.fasta > stx2_analysis/EHEC.fasta
```

```
 awk '
BEGIN {RS=">"; FS="\n"}
NR>1 {
    header=$1
    if (header ~ /\|STEC\|/) {
        print ">" $0
    }
}
' stx2_analysis/stx2_final.fasta > stx2_analysis/STEC.fasta
```
🧬 ✅ STEP 2 — Add reference to BOTH groups

Add reference
```
cat stx2A_ref.fasta stx2_analysis/EHEC.fasta > stx2_analysis/EHEC_with_ref.fasta
cat stx2A_ref.fasta stx2_analysis/STEC.fasta > stx2_analysis/STEC_with_ref.fasta
```
🧬 ✅ STEP 3 — Re-align separately

Using MAFFT:
```
mafft --auto stx2_analysis/EHEC_with_ref.fasta > stx2_analysis/EHEC_aligned.fasta
mafft --auto stx2_analysis/STEC_with_ref.fasta > stx2_analysis/STEC_aligned.fasta
```


## Mafft Alignment for both STEC and EHEC
```
mafft --auto stx2_analysis/stx2_final.fasta > stx2_analysis/stx2_aligned.fasta
```


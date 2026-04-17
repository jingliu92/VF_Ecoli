# VF_Ecoli

conda activate virulencefinder_env

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

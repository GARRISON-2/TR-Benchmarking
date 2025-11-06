##run nextflow
sbatch -J tr_bench -p amilan --qos=normal --time=08:00:00 --mem=2G --output=logs/%x_%j.out --error=logs/%x_%j.err --wrap="bash run_benchmarking_slurm.sh run_benchmarking.nf --list ont_bam.list --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta"


sbatch -J tr_bench -p amilan --qos=normal --time=08:00:00 --mem=2G --output=logs/%x_%j.out --error=logs/%x_%j.err --wrap="bash run_benchmarking_slurm.sh run_benchmarking_akshay.nf --list ont_bam.list --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta"

##Run_Atarva 
cd /pl/active/dashnowlab/projects/TR-benchmarking/

while IFS=$'\t' read -r CRAM KARYOTYPE; do   [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue;   SAMPLE="$(basename "$CRAM" .haplotagged.cram)";   SAMPLE="${SAMPLE%.cram}";   sbatch -J "atarva_${SAMPLE}" run_atarva.sh "$CRAM" "$KARYOTYPE"; done < test_bam.list

##Run LongTR
cd /pl/active/dashnowlab/projects/TR-benchmarking/

while IFS=$'\t' read -r CRAM KARYOTYPE; do   [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue;   SAMPLE="$(basename "$CRAM" .haplotagged.cram)";   SAMPLE="${SAMPLE%.cram}";   sbatch -J "longtr_${SAMPLE}" run_longtr.sh "$CRAM" "$KARYOTYPE"; done < test_bam.list

##Run Straglr
cd /pl/active/dashnowlab/projects/TR-benchmarking/

while IFS=$'\t' read -r CRAM KARYOTYPE; do   [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue;   SAMPLE="$(basename "$CRAM" .haplotagged.cram)";   SAMPLE="${SAMPLE%.cram}";   sbatch -J "straglr_${SAMPLE}" run_straglr.sh "$CRAM" "$KARYOTYPE"; done < test_bam.list

##Run TRSV
cd /pl/active/dashnowlab/projects/TR-benchmarking/

while IFS=$'\t' read -r CRAM KARYOTYPE; do   [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue;   SAMPLE="$(basename "$CRAM" .haplotagged.cram)";   SAMPLE="${SAMPLE%.cram}";   sbatch -J "trsv_${SAMPLE}" run_trsv.sh "$CRAM" "$KARYOTYPE"; done < test_bam.list
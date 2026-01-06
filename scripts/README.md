##Header fixer
sh header_fixer.sh /pl/active/dashnowlab/projects/TR-benchmarking/ont_bam.list   /pl/active/dashnowlab/projects/TR-benchmarking/fixed_header_bam/


##Catalog builder
https://www.twistbioscience.com/sites/default/files/resources/2023-02/Alliance_Dark_Genes_LR_Pnl_TargetsCaptured_hg38_ann.bed
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/SegmentalDuplications/GRCh38_segdups.bed.gz
https://raw.githubusercontent.com/dashnowlab/STRchive/refs/heads/main/data/catalogs/STRchive-disease-loci.hg38.general.bed


bash -x build_catalogs.sh "TR_catalog_for_vamos.shard_*_of_05.*.tsv.gz" /pl/active/dashnowlab/projects/TR-benchmarking/catalogs/hell_catalog/bed_files_filtering/test/ --include /pl/active/dashnowlab/projects/TR-benchmarking/catalogs/hell_catalog/bed_files_filtering/Alliance_Dark_Genes_LR_Pnl_TargetsCaptured_hg38_ann.bed /pl/active/dashnowlab/projects/TR-benchmarking/catalogs/hell_catalog/bed_files_filtering/GRCh38_segdups.bed


# Separate Shell Scripts

## Run_Atarva
```bash
cd /pl/active/dashnowlab/projects/TR-benchmarking/

while IFS=$'\t' read -r CRAM KARYOTYPE; do
  [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue
  SAMPLE="$(basename "$CRAM" .haplotagged.cram)"
  SAMPLE="${SAMPLE%.cram}"
  sbatch -J "atarva_${SAMPLE}" run_atarva.sh "$CRAM" "$KARYOTYPE"
done < test_bam.list
```

## Run LongTR
```bash
cd /pl/active/dashnowlab/projects/TR-benchmarking/

while IFS=$'\t' read -r CRAM KARYOTYPE; do
  [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue
  SAMPLE="$(basename "$CRAM" .haplotagged.cram)"
  SAMPLE="${SAMPLE%.cram}"
  sbatch -J "longtr_${SAMPLE}" run_longtr.sh "$CRAM" "$KARYOTYPE"
done < test_bam.list
```

## Run Straglr
```bash
cd /pl/active/dashnowlab/projects/TR-benchmarking/

while IFS=$'\t' read -r CRAM KARYOTYPE; do
  [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue
  SAMPLE="$(basename "$CRAM" .haplotagged.cram)"
  SAMPLE="${SAMPLE%.cram}"
  sbatch -J "straglr_${SAMPLE}" run_straglr.sh "$CRAM" "$KARYOTYPE"
done < test_bam.list
```

## Run TRSV
```bash
cd /pl/active/dashnowlab/projects/TR-benchmarking/

while IFS=$'\t' read -r CRAM KARYOTYPE; do
  [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue
  SAMPLE="$(basename "$CRAM" .haplotagged.cram)"
  SAMPLE="${SAMPLE%.cram}"
  sbatch -J "trsv_${SAMPLE}" run_trsv.sh "$CRAM" "$KARYOTYPE"
done < test_bam.list

```
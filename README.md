# CITUS

# Installation
```
curl -s https://get.nextflow.io | bash
```

# Usage

```
nextflow run bgsi-id/citus
    -r              dev \
    --input         sample.csv \
    --fasta         genome.fasta \
    --bwa_index     genome.bwa/ \
    --region        calling.bed \
    --known_site    indels.vcf.gz \
    --outdir        outdir/ \
    --gpu           single \
    -work-dir       work/ \
    -with-tower -resume
```


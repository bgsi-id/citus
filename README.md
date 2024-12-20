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
    --svd_prefix    's3://bgsi-data-dev/JY/1000g.phase3.100k.b38.vcf.gz.dat' \
    --outdir        outdir/ \
    --gpu           single \
    -work-dir       work/ \
    -with-tower -resume
```


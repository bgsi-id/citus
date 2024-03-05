process PB_GERMLINE {
  tag { meta.id }

  input:
  tuple val( meta ), path(reads_1), path(reads_2)
  path(fasta)
  path(bwa)
  path(known_site)

  output:
  tuple val( meta ), path('*.bam'), emit: bam
  tuple val( meta ), path('*.vcf'), emit: vcf
  tuple val( meta ), path('*chrs.txt'), emit: chrs
  tuple val( meta ), path('*recal.txt'), emit: recal
  path  "versions.yml", emit: versions

  script:
  """
  # nvidia-smi
  mv ${bwa}/* .
  pbrun germline \
  --ref ${fasta} \
  --in-fq ${reads_1} ${reads_2} \
  --knownSites ${known_site} \
  --out-bam ${ meta.id }.bam \
  --out-variants ${ meta.id }.vcf \
  --out-recal-file ${ meta.id }_recal.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      fastqc: \$( pbrun --version | sed -n 's/pbrun: \\([0-9.]*\\)/\\1/p' )
  END_VERSIONS
  """
}
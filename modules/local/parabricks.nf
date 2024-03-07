process PB_GERMLINE {
  tag { meta.id }

  input:
  tuple val( meta ), path(reads)
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
  def prefix = task.ext.prefix ?: "${meta.id}"
  def old_new_pairs = reads instanceof Path || reads.size() == 1 ? [[ reads, "${prefix}.${reads.extension}" ]] : reads.withIndex().collect { entry, index -> [ entry, "${prefix}_${index + 1}.${entry.extension}" ] }
  def rename_to = old_new_pairs*.join(' ').join(' ')
  def renamed_files = old_new_pairs.collect{ old_name, new_name -> new_name }.join(' ')
  """
  printf "%s %s\\n" $rename_to | while read old_name new_name; do
      [ -f "\${new_name}" ] || ln -s \$old_name \$new_name
  done

  nvidia-smi
  
  mv ${bwa}/* .
  pbrun germline \
  --ref ${fasta} \
  --in-fq ${renamed_files} \
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
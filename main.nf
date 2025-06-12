def errorMessage() {
  log.info(
    """
    ====================
    please provide input parameters --SAMPLEFILE or --BARCODEFILE
    ====================
    """.stripIndent()
  )
  exit(1)
}


process make_ref {
  publishDir "${params.outdir}/", mode: 'copy'

  input:
  path gtf

  output:
  path 'segments.csv'
  path 'gtf.rds'
  path 'segments.sajr'
  path 'functional_annotation'

  shell:
  '''
  # Convert GTF to SAJR reference format
  java -Xmx10G -jar !{projectDir}/bin/sajr.ss.jar \
    gff2sajr !{projectDir}/bin/sajr.config \
      -ann_foreign=!{gtf} \
      -ann_out=segments.sajr
  
  # Convert SAJR objects to segments.csv, gtf.rds, and functional_annotation
  Rshell !{projectDir}/bin/prepare_reference.R !{gtf} segments.sajr
  '''
}


process get_data {
  input:
  tuple val(id), val(bam_path)

  output:
  tuple val(id), path('*bam'), path('*bam.bai')

  shell:
  '''
  if "!{params.bam_on_irods}"; then
    iget -f -v -K "!{bam_path}" "!{id}.bam"
    iget -f -v -K "!{bam_path}.bai" "!{id}.bam.bai"
  else
    ln -s !{bam_path} !{id}.bam
    ln -s !{bam_path}.bai !{id}.bam.bai
  fi
  '''
}


process make_chr_list {
  input:
  path ref

  output:
  path ('chrs.txt'), emit: rds

  shell:
  '''
  # larger chrs first to be used for strand determination
  cut -f 1 !{ref}/segments.sajr | grep -v '#' | sort | uniq -c | sort -nr | sed  's/^[[:blank:]]*//' | cut -d ' ' -f2 > chrs.txt
  '''
}


process determine_strand {
  input:
  tuple val(id), path(bam_path), path(bami_path), path(ref), val(chr)

  output:
  tuple val(id), path(bam_path), path(bami_path), path(ref), stdout

  shell:
  '''
  mkdir p m
  
  !{projectDir}/bin/run_sajr.sh !{bam_path} p/p !{ref}/segments.sajr  1 !{chr} 1000000 !{params.use_bam_dedupl} > p/log &
  !{projectDir}/bin/run_sajr.sh !{bam_path} m/p !{ref}/segments.sajr -1 !{chr} 1000000 !{params.use_bam_dedupl} > m/log &
  wait
  
  m=`grep 'exon records' m/log | cut -d ' ' -f4`
  p=`grep 'exon records' p/log | cut -d ' ' -f4`
  
  # Output 1 if plus > minus (forward-stranded), -1 if plus <= minus (reverse-stranded)
  if [ $p -gt $m ]; then
    echo -n '1'
  else
    echo -n '-1'
  fi
  '''
}


process run_sajr_per_chr {
  input:
  tuple val(id), path(bam_path), path(bami_path), path(ref), val(strand), val(chr)

  output:
  tuple val(id), val(chr), path(id), path(bam_path), path(bami_path), val(strand)

  shell:
  '''
  mkdir !{id}
  !{projectDir}/bin/run_sajr.sh !{bam_path} !{id}/!{chr} !{ref}/segments.sajr !{strand} !{chr} -1 !{params.use_bam_dedupl} > !{id}/!{chr}.log
  '''
}


process combine_sajr_output {
  label "pseudobulk"

  input:
  path samples
  path barcodes
  path ref
  val n_samples

  output:
  path ('rds'), emit: rds

  shell:
  '''
  Rshell !{projectDir}/bin/combine_sajr_output.R !{samples} !{barcodes} !{ref} !{projectDir}/bin !{params.ncores}
  '''
}


process remake_pseudobulk {
  label "pseudobulk"

  input:
  path samples
  path barcodes
  path ref
  val n_samples

  output:
  path ('rds'), emit: rds

  shell:
  '''
  Rshell !{projectDir}/bin/remake_pseudobulk.R !{samples} !{barcodes} !{params.preprocessed_rds} !{ref} !{projectDir}/bin !{params.ncores}
  '''
}


process postprocess {
  publishDir "${params.outdir}/", mode: 'copy'

  input:
  path rds
  path samples
  path barcodes
  path ref
  val n_samples

  output:
  path ('rds'), emit: rds
  path 'examples'

  shell:
  '''
  Rshell !{projectDir}/bin/postprocess.R !{rds} !{params.mincells} !{params.minsamples} !{samples} !{barcodes} !{ref} !{projectDir}/bin !{params.ncores}
  '''
}


process generate_summary {
  publishDir "${params.outdir}/", mode: 'copy'

  input:
  path rds

  output:
  path 'summary.html'

  shell:
  '''
  cp !{projectDir}/bin/summary.Rmd .
  cp !{projectDir}/bin/sajr_utils.R .
  Rscript -e "wd=getwd();rmarkdown::render('summary.Rmd',
                     output_file = 'summary.html',
                     clean = TRUE)"
 '''
}

workflow reference {
  // Generate SAJR reference files from GTF
  make_ref(params.gtf)
}


workflow repseudobulk {
  ch_barcodes = Channel.fromPath(params.BARCODEFILE)
  ch_ref = Channel.fromPath(params.ref)
  ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
  ch_sample_list | flatMap { it.readLines() } | map { it -> [it.split()[0], it.split()[1]] } | get_data | set { ch_data }
  bam_path_file = ch_data.collectFile { item -> ["bam_paths.txt", item[0].toString() + " " + item[1] + "\n"] }

  remake_pseudobulk(ch_sample_list, Channel.fromPath(params.BARCODEFILE), ch_ref, ch_sample_list.countLines())
  postprocess(remake_pseudobulk.out.rds, bam_path_file, ch_barcodes, ch_ref, ch_sample_list.countLines())
  generate_summary(postprocess.out.rds)
}


workflow {
  // Parameter parsing & Channel setup
  ch_barcodes = Channel.fromPath(params.BARCODEFILE)
  ch_ref = Channel.fromPath(params.ref)
  ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
  ch_sample_list | flatMap { it.readLines() } | map { it -> [it.split()[0], it.split()[1]] } | get_data | set { ch_data }

  // Reference Preparation (future: integrate make_ref when --gtf is provided instead of --ref?)
  ch_chrs = make_chr_list(ch_ref).splitText().map { it.trim() }

  // Chromosome list & Strand determination
  ch_data = ch_data.combine(ch_ref).combine(ch_chrs.first())
  ch_data | determine_strand | set { ch_data }

  // Run SAJR per chromosome
  ch_data = ch_data.combine(ch_chrs)
  ch_data | run_sajr_per_chr | set { sajr_outs }
  sajrout_path_file = sajr_outs.collectFile { item -> ["sajr_outs.txt", item[0] + " " + item[1] + " " + item[2] + " " + item[3] + " " + item[5] + "\n"] }

  // Pseudobulk aggregation
  combine_sajr_output(sajrout_path_file, ch_barcodes, ch_ref, ch_sample_list.countLines())
  bam_path_file = ch_data.collectFile { item -> ["bam_paths.txt", item[0].toString() + " " + item[1] + "\n"] }

  // Post-processing & filtering
  postprocess(combine_sajr_output.out.rds, bam_path_file, ch_barcodes, ch_ref, ch_sample_list.countLines())

  // Report generation
  generate_summary(postprocess.out.rds)
}

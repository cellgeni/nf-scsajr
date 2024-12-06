def errorMessage() {
    log.info"""
    ====================
    please provide input parameters --SAMPLEFILE or --BARCODEFILE
    ====================
    """.stripIndent()
    exit 1
}


process make_ref {
  publishDir "${params.outdir}/", mode: 'copy'
  
  input:
  path gtf 
  
  output:
  path('segments.csv')
  path('gtf.rds')
  path('segments.sajr')
  path('functional_annotation')
  
  shell:
  '''
  java -Xmx10G -jar !{projectDir}/bin/sajr.ss.jar \
  gff2sajr !{projectDir}/bin/sajr.config \
  -ann_foreign=!{gtf} \
  -ann_out=segments.sajr
  
  Rscript !{projectDir}/bin/prepare_reference.R !{gtf} segments.sajr
  '''
}


process get_data {
  input:
  tuple val(id), val(bam_path)

  output:
  tuple val(id),path('*bam'),path('*bam.bai')

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


process run_sajr {
  
  memory = {10.GB + bam_path.size() / 10 * 11.B + 10.GB * task.attempt}
  
  input:
  tuple val(id), path(bam_path), path(bam_path_id), val(strand), path(ref)
  
  output:
  tuple path(id), val(id), path(bam_path), path(bam_path_id), val(strand)
  
  shell:
  '''
  mkdir !{id}
  
  java -Xmx60G -jar !{projectDir}/bin/sajr.ss.jar \
    count_reads !{projectDir}/bin/sajr.config \
    -batch_in=!{bam_path}\
    -batch_out=!{id}/!{id} \
    -ann_in=!{ref}/segments.sajr \
    -stranded=!{strand} \
    -count_only_border_reads=true
  '''
}


process combine_sajr_output {
 label "pseudobulk"
 
 input:
 path(samples)
 path(barcodes)
 path(ref)
 val(n_samples)
 
 output:
 path('rds'), emit: rds
 
 shell:
 '''
 Rscript !{projectDir}/bin/combine_sajr_output.R !{samples} !{barcodes} !{ref} !{projectDir}/bin  !{params.ncores}
 '''
}

process remake_pseudobulk {
 label "pseudobulk"
 
 input:
 path(samples)
 path(barcodes)
 path(ref)
 val(n_samples)
 
 output:
 path('rds'), emit: rds
 
 shell:
 '''
 Rscript !{projectDir}/bin/remake_pseudobulk.R !{samples} !{barcodes} !{params.preprocessed_rds} !{ref} !{projectDir}/bin !{params.ncores}
 '''
}

process postprocess {
 publishDir "${params.outdir}/", mode: 'copy'
 
 memory = {30.GB + 20.GB * task.attempt}
 
 input:
 path(rds)
 path(samples)
 path(barcodes)
 path(ref)
 
 output:
 path('rds'), emit: rds
 path('examples.pdf')
 
 shell:
 '''
 Rscript !{projectDir}/bin/postprocess.R !{rds} !{params.mincells} !{params.minsamples} !{samples} !{barcodes} !{ref} !{projectDir}/bin !{params.ncores}
 '''
}


process generate_summary {
 publishDir "${params.outdir}/", mode: 'copy'
 
 input:
 path(rds)
 
 output:
 path('summary.html')
 
 shell:
 '''
 cp !{projectDir}/bin/summary.Rmd .
 cp !{projectDir}/bin/sajr_utils.R .
 Rscript -e "wd=getwd();rmarkdown::render('summary.Rmd',
                    output_file = 'summary.html',
                    clean = TRUE)"
 '''
}

workflow  {
  ch_barcodes = Channel.fromPath(params.BARCODEFILE)
  ch_ref = Channel.fromPath(params.ref)
  ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
  ch_sample_list | flatMap{ it.readLines() } | map { it -> [ it.split()[0], it.split()[1] ] } | get_data | set { ch_data }
  ch_data = ch_data.combine(Channel.of(-1,1))
  ch_data = ch_data.combine(ch_ref)
  ch_data | run_sajr | set {sajr_outs}
  sajrout_path_file = sajr_outs.collectFile{item -> ["sajr_outs.txt", item[0].toString()  + " " + item[1] + " " + item[2] + " " + item[3] + " " + item[4] + "\n"]}
  combine_sajr_output(sajrout_path_file,ch_barcodes,ch_ref,sajr_outs.count())
  
  bam_path_file = ch_data.collectFile{item -> ["bam_paths.txt", item[0].toString()  + " " + item[1] + "\n"]}
  postprocess(combine_sajr_output.out.rds,bam_path_file,ch_barcodes,ch_ref)
  generate_summary(postprocess.out.rds)
}

workflow repseudobulk {
  ch_barcodes = Channel.fromPath(params.BARCODEFILE)
  ch_ref = Channel.fromPath(params.ref)
  ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
  ch_sample_list | flatMap{ it.readLines() } | map { it -> [ it.split()[0], it.split()[1] ] } | get_data | set { ch_data }
  bam_path_file = ch_data.collectFile{item -> ["bam_paths.txt", item[0].toString()  + " " + item[1] + "\n"]}
  
  remake_pseudobulk(ch_sample_list,Channel.fromPath(params.BARCODEFILE),ch_ref,ch_sample_list.countLines())
  postprocess(remake_pseudobulk.out.rds,bam_path_file,ch_barcodes,ch_ref)
  generate_summary(postprocess.out.rds)
}


workflow reference {
  make_ref(params.gtf)
}

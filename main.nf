def errorMessage() {
    log.info"""
    ====================
    please provide input parameters --SAMPLEFILE or --BARCODEFILE
    ====================
    """.stripIndent()
    exit 1
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
    #cp "!{bam_path}" "!{id}.bam"
    #cp "!{bam_path}.bai" "!{id}.bam.bai"
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
    -ann_in=!{ref}\
    -stranded=!{strand} \
    -count_only_border_reads=true
  '''
}


process combine_sajr_output {
 publishDir "${params.outdir}/", mode: 'copy'
 
 memory = {10.GB + n_samples * 300.MB  + 10.GB * task.attempt}
 
 input:
 path(samples)
 path(barcodes)
 val(n_samples)
 
 output:
 path('rds'), emit: rds
 path('examples.pdf')
 
 shell:
 '''
 Rscript !{projectDir}/bin/combine_sajr_output.R !{samples} !{barcodes} !{projectDir}/data  !{params.ref}
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
 Rscript -e "wd=getwd();rmarkdown::render('summary.Rmd',
                    output_file = 'summary.html',
                    clean = TRUE)"
 '''
}

workflow {
  ch_sample_list = params.SAMPLEFILE != null ? Channel.fromPath(params.SAMPLEFILE) : errorMessage()
  ch_sample_list | flatMap{ it.readLines() } | map { it -> [ it.split()[0], it.split()[1] ] } | get_data | set { ch_data }
  ch_data = ch_data.combine(Channel.of(-1,1))
  ch_data = ch_data.combine(Channel.fromPath(params.ref))
  ch_data | run_sajr | set {sajr_outs}
  sajrout_path_file = sajr_outs.collectFile{item -> ["samples.txt", item[0].toString()  + " " + item[1] + " " + item[2] + " " + item[3] + " " + item[4] + "\n"]}
  combine_sajr_output(sajrout_path_file,Channel.fromPath(params.BARCODEFILE),sajr_outs.count())
  generate_summary(combine_sajr_output.out.rds)
}

#!/usr/bin/env nextflow

//take input as --reports and--metadata

/*
 * A Python script which parses the output of the previous script
 */
process make_pcoa_script {

    container 'community.wave.seqera.io/library/pip_mako_matplotlib_natsort_pruned:55548172648d6af5'

    input:
    path reports
    path metadata
    path site_key

    output:
    path "text_files"

    script:
    """
    make_pcoa_files.py \
        --reports ${reports.join(' ')} \
        --metadata ${metadata} \
        --site_key ${site_key} \
        --r_dir text_files
    """
}

/*
 * An R script which produces output for shannon's diversity graphs
 */

process get_pcoa_plot {

    container 'community.wave.seqera.io/library/r-ade4_r-ecodist_r-vegan:5d7d5c3d409fc0d1'
    input:
    path text_files
    output:
    path "plots"

    script:
    """
    make_site_pcoa.R ${text_files}/ plots/
    """
}

/*
 * A python script which produces output for making the html report
 */

process make_site_report {

    container 'community.wave.seqera.io/library/pip_mako_matplotlib_natsort_pruned:55548172648d6af5'

    input:
    path reports
    path metadata
    path site_key
    path template
    path reference
    path plots

    output:
    path "site_reports/"

    publishDir "${params.outdir}/", mode: 'copy' // Publish final report to local directory specified in params.config

    script:
    """
    make_site_report.py \
      --reports ${reports.join(' ')} \
      --metadata ${metadata} \
      --site_key ${site_key} \
      --plots_dir ${plots}/ \
      --final_report site_reports/ \
      --template ${template} \
      --reference ${reference}
    """
}


workflow evaluate_by_site {
    report_list = params.reports?.split('\n') as List
    Channel
        .fromPath(report_list)
        .flatten()
        .collect()
        .set { reports }
    //reports.view()


    metadata_file = file(params.metadata, type: "file", checkIfExists:true)
    Channel
        .fromPath(metadata_file)
        .set { metadata }

    //hcids.view()
    //rename(hcids)
    

    site_key = file(params.site_key, type: "file", checkIfExists:true)

    template = file("$baseDir/bin/site_report_template.html")
    reference = file("$baseDir/bin/contaminant_literature.xlsx")

    make_pcoa_script(reports, metadata, site_key)
    get_pcoa_plot(make_pcoa_script.out)
    make_site_report(reports, metadata, site_key, template, reference, get_pcoa_plot.out)
    println "Report will be generated in ${params.outdir}"
}

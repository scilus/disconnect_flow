#!/usr/bin/env nextflow

import groovy.json.*

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["output_dir":"$params.output_dir",
                "no_pruning":"$params.no_pruning",
                "no_remove_loops":"$params.no_remove_loops",
                "no_remove_outliers":"$params.no_remove_outliers",
                "min_length":"$params.min_length",
                "max_length":"$params.max_length",
                "loop_max_angle":"$params.loop_max_angle",
                "outlier_threshold":"$params.outlier_threshold",
                "nbr_subjects_for_avg_connections":"$params.nbr_subjects_for_avg_connections",
                "processes_register":"$params.processes_register",
                "processes_connectivity":"$params.processes_connectivity",
                "cpu_count":"$cpu_count"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info "Disconnectome pipeline"
log.info "==================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}


if (params.input){
    log.info "Input: $params.input"
    root = file(params.input)
    in_lesions = Channel
        .fromFilePairs("$root/**/*{cavity*.nii.gz}",
                       size: 1,
                       maxDepth:1,
                       flat: true) {it.parent.name}
}

if (params.tractograms){
    log.info "Input: $params.tractograms"
    tractograms = file(params.tractograms)
    in_tractograms = Channel
      .fromFilePairs("$tractograms/**/*trk",
                     size: 1,
                     maxDepth: 1,
                     flat: true) {it.parent.name}
}

if (params.atlases){
    log.info "Input: $params.atlases"
    atlases = file(params.atlases)
    in_atlases = Channel
      .fromFilePairs("$atlases/**/{atlas_*nii.gz,labels_*.txt}",
                     size: 2,
                     maxDepth: 1,
                     flat: true) {it.parent.name}
}

in_tractograms.combine(in_atlases)
    .set{trk_atlases_for_decompose_connectivity}

process README {
    cpus 1
    publishDir = params.Readme_Publish_Dir
    tag = "README"

    output:
    file "readme.txt"

    script:
    String list_options = new String();
    for (String item : params) {
        list_options += item + "\n"
    }

    """
    echo "Disconnectome pipeline\n" >> readme.txt
    echo "Start time: $workflow.start\n" >> readme.txt
    echo "[Command-line]\n$workflow.commandLine\n" >> readme.txt
    echo "[Options]\n" >> readme.txt
    echo "$list_options" >> readme.txt
    """
}

process Decompose_Connectivity {
    cpus 1
    memory { 6 * trackings.size() }

    input:
    set sid, file(trackings), atlas_name, file(atlas), file(atlas_label) from trk_atlases_for_decompose_connectivity

    output:
    set sid, atlas_name, file(atlas), file(atlas_label), "${sid}_${atlas_name}__decompose.h5" into h5_for_combine_with_lesion

    script:
    no_pruning_arg = ""
    if (params.no_pruning) {
        no_pruning_arg = "--no_pruning"
    }
    no_remove_loops_arg = ""
    if (params.no_remove_loops) {
        no_remove_loops_arg = "--no_remove_loops"
    }
    no_remove_outliers_arg = ""
    if (params.no_pruning) {
        no_remove_outliers_arg = "--no_pruning"
    }
    no_remove_outliers_arg = ""
    if (params.no_remove_outliers) {
        no_remove_outliers_arg = "--no_remove_outliers"
    }
    """
    if [ `echo $trackings | wc -w` -gt 1 ]; then
        scil_streamlines_math.py lazy_concatenate $trackings tracking_concat.trk
    else
        mv $trackings tracking_concat.trk
    fi

    scil_decompose_connectivity.py tracking_concat.trk $atlas "${sid}_${atlas_name}__decompose.h5" --no_remove_curv_dev \
        $no_pruning_arg $no_remove_loops_arg $no_remove_outliers_arg --min_length $params.min_length \
        --max_length $params.max_length --loop_max_angle $params.loop_max_angle \
        --outlier_threshold $params.outlier_threshold
    """
}


h5_for_combine_with_lesion.combine(in_lesions)
  .into{h5_labels_lesion_for_compute_connectivity;toto}

process Compute_Connectivity_without_similiarity {
    cpus params.processes_connectivity
    publishDir = {"${params.output_dir}/$lesion_id/$sid/Compute_Connectivity"}

    input:
    set sid, atlas_name, file(atlas),  file(atlas_label), file(h5), lesion_id, file(lesion) from h5_labels_lesion_for_compute_connectivity

    output:
    set sid, "$atlas_name/*.npy" into matrices_for_visualize

    script:
    """
    mkdir $atlas_name

    scil_compute_connectivity.py $h5 $atlas --force_labels_list $atlas_label \
        --volume vol.npy --streamline_count sc.npy \
        --length len.npy --density_weighting \
        --no_self_connection --include_dps ./ --lesion_load $lesion $atlas_name \
        --processes $params.processes_connectivity

    rm rd_fixel.npy -f
    scil_normalize_connectivity.py sc.npy sc_edge_normalized.npy \
        --parcel_volume $atlas $atlas_label
    scil_normalize_connectivity.py vol.npy sc_vol_normalized.npy \
        --parcel_volume $atlas $atlas_label
    """
}

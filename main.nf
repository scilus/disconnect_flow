#!/usr/bin/env nextflow

import groovy.json.*

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["output_dir":"$params.output_dir",
                "run_bet":"$params.run_bet",
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
    log.info "Input lesions: $params.input"
    root = file(params.input)
    in_lesions = Channel.fromFilePairs("$root/**/*{$params.lesion_name*.nii.gz}",
                                    size: 1,
                                    maxDepth:1,
                                    flat: true) {it.parent.name}

    Channel.fromPath("$root/**/*t1.nii.gz",
                     maxDepth:1)
                     .map{[it.parent.name, it]}
                     .into{t1s_for_register; check_t1s}
}

if (params.tractograms){
    log.info "Input tractograms: $params.tractograms"
    tractograms = file(params.tractograms)
    in_tractograms = Channel.fromFilePairs("$tractograms/**/*trk",
                                           size: 1,
                                           maxDepth: 1,
                                           flat: true) {it.parent.name}
}

if (params.atlas){
    log.info "Input atlas: $params.atlas"
    atlas = file(params.atlas)
    in_atlas = Channel
      .fromFilePairs("$atlas/{atlas_labels.nii.gz,atlas_t1.nii.gz,atlas_list.txt}",
                     size: 3,
                     maxDepth: 0,
                     flat: true) {it.parent.name}
}

in_atlas.into{atlas_for_combine; atlas_for_registration; check_atlas}

in_lesions.into{lesions; lesions_for_registration; check_lesions}

in_tractograms.into{tractograms_for_combine; check_trks}

// Check Lesions
check_lesions.count().into{check_lesions_number; check_lesions_number_compare_t1}
check_t1s.count().into{check_t1s_number; check_t1s_number_for_registration}

check_lesions_number
.subscribe{a -> if (a == 0)
    error "Error ~ No subjects found. Please check the naming convention, your --input path."}

check_lesions_number_compare_t1
  .concat(check_t1s_number)
  .toList()
  .subscribe{a, b -> if (a != b && b > 0)
  error "Error ~ Some subjects have a T1w and others don't.\n" +
        "Please be sure to have the same acquisitions for all subjects."}

// Check TRKs
check_trks
.subscribe{a -> if (a == 0)
    error "Error ~ No tractograms found. Please check the naming convention, your --tractograms path."
}

tractograms_for_combine.combine(atlas_for_combine)
    .set{trk_atlases_for_decompose_connectivity}


if (params.linear_registration & params.nonlinear_registration){
  error "Error ~ You cannot select both linear_registration and nonlinear_registration."
}
else if(params.linear_registration){
  registration_strategy="-t a"
}
else if(params.nonlinear_registration){
  registration_strategy="-t s"
}

check_t1s_number_for_registration.subscribe{a -> if (a > 0 &&  !params.linear_registration && !params.nonlinear_registration)
  error "You provided t1 in order to register the lesion into the atlas space without specifying the registration strategy. Please add --linear_registration or --nonlinear_registration."
}

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


process Register_T1 {
    cpus params.processes_bet_register_t1

    input:
    set sid, file(t1) from t1s_for_register
    set atlas_name, file(atlas_labels), file(atlas_txt), file(atlas_t1) from atlas_for_registration

    output:
    set sid, atlas_name, "${sid}__output0GenericAffine.mat", "${sid}__t1_${atlas_name}_space.nii.gz" into transformation_for_registration_lesions
    file "${sid}__t1_bet_mask.nii.gz" optional true
    file "${sid}__t1_bet.nii.gz" optional true

    script:
    if (params.run_bet){
    """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export ANTS_RANDOM_SEED=1234

        antsBrainExtraction.sh -d 3 -a $t1 -e $params.template_t1/t1_template.nii.gz\
            -o bet/ -m $params.template_t1/t1_brain_probability_map.nii.gz -u 0
        scil_image_math.py convert bet/BrainExtractionMask.nii.gz ${sid}__t1_bet_mask.nii.gz --data_type uint8
        scil_image_math.py multiplication $t1 ${sid}__t1_bet_mask.nii.gz ${sid}__t1_bet.nii.gz

        antsRegistrationSyN.sh -d 3 -m ${sid}__t1_bet.nii.gz -f ${atlas_t1} -n ${task.cpus} -o "${sid}__output" ${registration_strategy}
        mv ${sid}__outputWarped.nii.gz ${sid}__t1_${atlas_name}_space.nii.gz
    """
    }
    else{
    """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export ANTS_RANDOM_SEED=1234

        antsRegistrationSyN.sh -d 3 -m ${t1} -f ${atlas_t1} -n ${task.cpus} -o "${sid}__output" ${registration_strategy}
        mv ${sid}__outputWarped.nii.gz ${sid}__t1_${atlas_name}_space.nii.gz
    """
    }
}

lesions_for_registration
    .join(transformation_for_registration_lesions)
    .set{lesion_mat_for_transformation}


process Transform_Lesions {
    cpus 1

    input:
    set sid, file(lesion), atlas_name, file(mat), file(t1_ref) from lesion_mat_for_transformation

    output:
    set sid, "${sid}__${params.lesion_name}_${atlas_name}_space_int16.nii.gz" into transformed_lesions
    file "${sid}__${params.lesion_name}_${atlas_name}_space.nii.gz"

    script:
    """
    antsApplyTransforms -d 3 -i $lesion -r $t1_ref -o ${sid}__${params.lesion_name}_${atlas_name}_space.nii.gz -t $mat -n NearestNeighbor
    scil_image_math.py convert ${sid}__${params.lesion_name}_${atlas_name}_space.nii.gz "${sid}__${params.lesion_name}_${atlas_name}_space_int16.nii.gz" --data_type int16
    """
}


process Decompose_Connectivity {
    cpus 1
    memory { 6 * trackings.size() }

    input:
    set sid, file(trackings), atlas_name, file(atlas), file(atlas_label), file(atlas_t1) from trk_atlases_for_decompose_connectivity

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
        scil_streamlines_math.py lazy_concatenate $trackings tracking_concat.trk --ignore_invalid
    else
        mv $trackings tracking_concat.trk
    fi

    scil_decompose_connectivity.py tracking_concat.trk $atlas "${sid}_${atlas_name}__decompose.h5" --no_remove_curv_dev \
        $no_pruning_arg $no_remove_loops_arg $no_remove_outliers_arg --min_length $params.min_length \
        --max_length $params.max_length --loop_max_angle $params.loop_max_angle \
        --outlier_threshold $params.outlier_threshold
    """
}

h5_for_combine_with_lesion.combine(transformed_lesions)
  .set{h5_labels_lesion_for_compute_connectivity}

process Compute_Connectivity_Lesion_without_similiarity {
    cpus params.processes_connectivity
    publishDir = {"${params.output_dir}/$lesion_id/$sid/Compute_Connectivity"}

    input:
    set sid, atlas_name, file(atlas),  file(atlas_label), file(h5), lesion_id, file(lesion) from h5_labels_lesion_for_compute_connectivity

    output:
    set sid, lesion_id, "*.npy", "connectivity_w_lesion/*.npy" into matrices_for_connectivity_in_csv

    script:
    """
    mkdir connectivity_w_lesion

    scil_compute_connectivity.py $h5 $atlas --force_labels_list $atlas_label \
        --volume atlas_vol.npy --streamline_count atlas_sc.npy \
        --length atlas_len.npy \
        --include_dps ./ --lesion_load $lesion connectivity_w_lesion/ \
        --processes $params.processes_connectivity

    rm rd_fixel.npy -f
    scil_normalize_connectivity.py atlas_sc.npy atlas_sc_edge_normalized.npy \
        --parcel_volume $atlas $atlas_label
    scil_normalize_connectivity.py atlas_vol.npy atlas_sc_vol_normalized.npy \
        --parcel_volume $atlas $atlas_label
    """
}


process Connectivity_in_csv {
    cpus 1
    publishDir = {"${params.output_dir}/$lesion_id/$sid/Compute_Connectivity"}

    input:
    set sid, lesion_id, file(atlas_matrices), file(matrices_w_lesion) from matrices_for_connectivity_in_csv

    output:
    set sid, "*csv", "connectivity_w_lesion/*.csv"

    script:
    String matrices_list = atlas_matrices.join("\",\"")
    String matrices_w_lesion = matrices_w_lesion.join("\",\"")
    """
    #!/usr/bin/env python3
    import numpy as np
    import os, sys

    os.mkdir("connectivity_w_lesion")

    for data in ["$matrices_list","$matrices_w_lesion"]:
      fmt='%1.8f'
      if 'sc' in data:
        fmt='%i'

      curr_data = np.load(data)
      if "lesion" in data:
        np.savetxt(os.path.join("connectivity_w_lesion/", data.replace(".npy", ".csv")), curr_data, delimiter=",", fmt=fmt)
      else:
        np.savetxt(data.replace(".npy", ".csv"), curr_data, delimiter=",", fmt=fmt)
    """
}

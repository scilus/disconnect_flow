#!/usr/bin/env nextflow

import groovy.json.*

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["registration_script":"$params.registration_script",
                "linear_registration": "$params.linear_registration",
                "slow_registration": "$params.slow_registration",
                "output_dir":"$params.output_dir",
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


log.info "Inputs"
log.info "============================="
log.info "Input lesions: $params.input"
log.info "Input tractograms: $params.tractograms"
log.info "Input atlas: $params.atlas"
log.info ""


log.info "Options - registration"
log.info "============================="
log.info "Registration strategy: $params.registration_strategy"
log.info "Registration script: $params.registration_script"
log.info ""


workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}


if (params.input){
    root = file(params.input)
    Channel.fromPath("$root/**/cavity.nii.gz",
                     maxDepth:1)
          .map{[it.parent.name, it]}
          .into{lesions; lesions_for_registration; check_lesions}


    Channel.fromPath("$root/**/*t1.nii.gz",
                     maxDepth:1)
          .map{[it.parent.name, it]}
          .into{t1s_for_register; check_t1s}
}

if (params.tractograms){
    tractograms = file(params.tractograms)
    Channel.fromFilePairs("$tractograms/**/*trk",
                          size: 1,
                          maxDepth: 1,
                          flat: true) {it.parent.name}
          .into{tractograms_for_combine; check_trks}
}

if (params.atlas){
    atlas = file(params.atlas)
    Channel.fromFilePairs("$atlas/{atlas_labels.nii.gz,atlas_labels.txt,atlas_t1.nii.gz,atlas_list.txt}",
                          size: 4,
                          maxDepth: 0,
                          flat: true) {it.parent.name}
      .into{atlas_for_combine; atlas_for_registration; check_atlas}

}

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
    set atlas_name, file(atlas), file(atlas_labels), file(atlas_list), file(atlas_t1) from atlas_for_registration

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

        ${params.registration_script} -d 3 -m ${sid}__t1_bet.nii.gz -f ${atlas_t1} -n ${task.cpus} -o "${sid}__output" -t ${params.registration_strategy}
        mv ${sid}__outputWarped.nii.gz ${sid}__t1_${atlas_name}_space.nii.gz
    """
    }
    else{
    """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export ANTS_RANDOM_SEED=1234

        ${params.registration_script} -d 3 -m ${t1} -f ${atlas_t1} -n ${task.cpus} -o "${sid}__output" ${params.registration_strategy}
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
    scil_image_math.py convert ${sid}__${params.lesion_name}_${atlas_name}_space.nii.gz ${sid}__${params.lesion_name}_${atlas_name}_space_int16.nii.gz --data_type int16
    """
}


process Decompose_Connectivity {
    cpus 1
    memory { 6 * trackings.size() }

    input:
    set sid, file(trackings), atlas_name, file(atlas), file(atlas_labels), file(atlas_list), file(atlas_t1) from trk_atlases_for_decompose_connectivity

    output:
    set sid, atlas_name, file(atlas), file(atlas_labels), file(atlas_list), "${sid}_${atlas_name}__decompose.h5" into h5_for_combine_with_lesion, h5_for_combine_with_lesion_echo

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

transformed_lesions.ifEmpty(lesions).set{lesion_for_connectivity}

h5_for_combine_with_lesion.combine(lesion_for_connectivity)
  .set{h5_labels_lesion_for_compute_connectivity}

process Compute_Connectivity_Lesion_without_similiarity {
    cpus params.processes_connectivity
    publishDir = {"${params.output_dir}/$lesion_id/$sid/Compute_Connectivity"}

    input:
    set sid, atlas_name, file(atlas), file(atlas_labels), file(atlas_list), file(h5), lesion_id, file(lesion) from h5_labels_lesion_for_compute_connectivity

    output:
    set sid, lesion_id, "*.npy", "Connectivity_w_lesion/*.npy" into matrices_for_connectivity_in_csv
    set sid, lesion_id, "$atlas_labels", "$atlas_list", "Connectivity_w_lesion/lesion_sc.npy" into lesion_sc_for_visualisation

    script:
    """
    mkdir Connectivity_w_lesion

    scil_compute_connectivity.py $h5 $atlas --force_labels_list $atlas_list \
        --volume atlas_vol.npy --streamline_count atlas_sc.npy \
        --length atlas_len.npy \
        --include_dps ./ --lesion_load $lesion Connectivity_w_lesion/ \
        --processes $params.processes_connectivity

    rm rd_fixel.npy -f
    scil_normalize_connectivity.py atlas_sc.npy atlas_sc_edge_normalized.npy \
        --parcel_volume $atlas $atlas_list
    scil_normalize_connectivity.py atlas_vol.npy atlas_sc_vol_normalized.npy \
        --parcel_volume $atlas $atlas_list
    """
}


process Connectivity_in_csv {
    cpus 1
    publishDir = {"${params.output_dir}/$lesion_id/$sid/Compute_Connectivity"}

    input:
    set sid, lesion_id, file(atlas_matrices), file(matrices_w_lesion) from matrices_for_connectivity_in_csv

    output:
    set sid, "*csv", "Connectivity_w_lesion/*.csv"

    script:
    String matrices_list = atlas_matrices.join("\",\"")
    String matrices_w_lesion = matrices_w_lesion.join("\",\"")
    """
    #!/usr/bin/env python3
    import numpy as np
    import os, sys

    os.mkdir("Connectivity_w_lesion")

    for data in ["$matrices_list","$matrices_w_lesion"]:
      fmt='%1.8f'
      if 'sc' in data:
        fmt='%i'

      curr_data = np.load(data)
      if "lesion" in data:
        np.savetxt(os.path.join("Connectivity_w_lesion/", data.replace(".npy", ".csv")), curr_data, delimiter=",", fmt=fmt)
      else:
        np.savetxt(data.replace(".npy", ".csv"), curr_data, delimiter=",", fmt=fmt)
    """
}


process Visualize_Connectivity {
    cpus 1
    publishDir = {"${params.output_dir}/$lesion_id/$sid/Compute_Connectivity/Connectivity_w_lesion"}

    input:
    set sid, lesion_id, file(atlas_labels), file(atlas_list), file(matrices) from lesion_sc_for_visualisation

    output:
    set sid, "*.png"

    script:
    String matrices_list = matrices.join(", ").replace(',', '')
    """
    for matrix in "$matrices_list"; do
        scil_visualize_connectivity.py \$matrix \${matrix/.npy/_matrix.png} --labels_list $atlas_list --name_axis \
            --display_legend --lookup_table $atlas_labels --histogram \${matrix/.npy/_hist.png} --nb_bins 50 --exclude_zeros --axis_text_size 5 5
    done
    """
}

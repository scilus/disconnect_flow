#!/usr/bin/env nextflow

import groovy.json.*

params.root = false
params.tractograms = false
params.help = false

if(params.help) {
    usage = file("$baseDir/USAGE")
    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make()

    bindings = ["tractogramAtlas":"$params.tractogramAtlas"]

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

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

if (params.root){
    log.info "Input: $params.root"
    root = file(params.root)
    in_data = Channel
        .fromFilePairs("$root/**/*{cavity*.nii.gz}",
                       size: 1,
                       maxDepth:1,
                       flat: true) {it.simpleName}
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

in_data.into{in_data_cc;
             in_data_cs;
             in_data_ct}

in_tractograms.into{in_tractograms_cc;
                    in_tractograms_cs;
                    in_tractograms_ct}

indices_cc = params.indices_cc?.tokenize(',')
indices_cs = params.indices_cs?.tokenize(',')
indices_ct = params.indices_ct?.tokenize(',')
sides = params.sides?.tokenize(',')

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
    echo "[Git Info]\n" >> readme.txt
    echo "$workflow.repository - $workflow.revision [$workflow.commitId]\n" >> readme.txt
    echo "[Options]\n" >> readme.txt
    echo "$list_options" >> readme.txt
    """
}

process filterCorticoCortical{
  cpus 1
  publishDir = params.FilterCorticoThalamic_Publish_Dir  
  tag = "Filter Cortico-Cortical"

  input:
    set tid, file(tractogram) from in_tractograms_cc
    each idx from indices_cc

  output:
    set tid, idx, "${tid}_CorticoCortical_${idx}.trk" into corticocortical

  when:
    params.tractogramAtlas==false

  script:
    """
        scil_remove_invalid_streamlines.py ${tractogram} valid_${tid}_${idx}.trk
        scil_filter_tractogram.py valid_${tid}_${idx}.trk \
          ${tid}_CorticoCortical_${idx}.trk \
          --drawn_roi ${params.atlasFolder}cortical17/${idx}.nii.gz both_ends include
    """
}

process filterCorticoStriatal{
  cpus 1
  publishDir = params.FilterCorticoStriatal_Publish_Dir
  tag = "Filter Cortical-Striatal"

  input:
    set tid, file(tractogram) from in_tractograms_cs
    each idx from indices_cs

  output:
    set tid, idx, "${tid}_CorticoStriatal_${idx}.trk" into corticostriatal

  when:
    params.tractogramAtlas==false

  script:
    """
        scil_remove_invalid_streamlines.py ${tractogram} valid_${tid}_${idx}.trk
        scil_filter_tractogram.py valid_${tid}_${idx}.trk \
          ${tid}_CorticoStriatal_${idx}.trk \
          --drawn_roi ${params.atlasFolder}cortical17/${idx}.nii.gz either_end include \
          --drawn_roi ${params.atlasFolder}striatal17/${idx}.nii.gz either_end include
    """
}

process filterCorticoThalamic{
  cpus 1
  publishDir = params.FilterCorticoThalamic_Publish_Dir
  tag = "Filter Cortical-Thalamic"


  input:
    set tid, file(tractogram) from in_tractograms_ct
    each idx from indices_ct

  output:
    set tid, idx, "${tid}_CorticoThalamic_${idx}.trk" into corticothalamic

  when:
    params.tractogramAtlas==false

  script:
    """
        scil_remove_invalid_streamlines.py ${tractogram} valid_${tid}_${idx}.trk
        scil_filter_tractogram.py valid_${tid}_${idx}.trk \
          ${tid}_CorticoThalamic_${idx}.trk \
          --drawn_roi ${params.atlasFolder}cortical17/${idx}.nii.gz either_end include \
          --drawn_roi ${params.atlasFolder}thalamic17/${idx}.nii.gz either_end include
    """
}

in_data_cc.combine(corticocortical).set{in_data_corticocortical}
in_data_cs.combine(corticostriatal).set{in_data_corticostriatal}
in_data_ct.combine(corticothalamic).set{in_data_corticothalamic}

process filterLesionCorticoCortical{
  input:
    set sid, file(lesion), tid, idx, file(tractogram) from in_data_corticocortical
    each side from sides

  output:
    set sid, tid, idx, side, "${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}.txt", "${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}_${params.minL}.txt" into corticocorticalLesion
    set sid, idx into corticocorticalIndices
    file "${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}.trk"
    file "${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}_${params.minL}.trk"

  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}.trk \
                --bdo ${params.atlasFolder}${side}HemisphereMNI.bdo both_ends include \
                --drawn_roi ${lesion} any include \
                -f --display_counts > ${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}.txt

  scil_filter_streamlines_by_length.py ${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}.trk \
              ${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}_${params.minL}.trk \
              --minL ${params.minL} -f --display_counts > ${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}_${params.minL}.txt
  """
}

process filterLesionCorticoStriatal{
  input:
    set sid, file(lesion), tid, idx, file(tractogram) from in_data_corticostriatal
    each side from sides

  output:
    set sid, tid, idx, side, "${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}.txt", "${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}_${params.minL}.txt" into corticostriatalLesion
    file "${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}.trk"
    file "${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}_${params.minL}.trk"

  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}.trk \
                --bdo ${params.atlasFolder}${side}HemisphereMNI.bdo both_ends include \
                --drawn_roi ${lesion} any include \
                -f --display_counts > ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}.txt

  scil_filter_streamlines_by_length.py ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}.trk \
              ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}_${params.minL}.trk \
              --minL ${params.minL} -f --display_counts > ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}_${params.minL}.txt
  """
}

process filterLesionCorticoThalamic{
  input:
    set sid, file(lesion), tid, idx, file(tractogram) from in_data_corticothalamic
    each side from sides

  output:
    set sid, tid, idx, side, "${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}.txt", "${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}_${params.minL}.txt" into corticothalamicLesion
    file "${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}.trk"
    file "${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}_${params.minL}.trk"

  script:
  """
  scil_filter_tractogram.py ${tractogram} ${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}.trk \
                --bdo ${params.atlasFolder}${side}HemisphereMNI.bdo both_ends include \
                --drawn_roi ${lesion} any include \
                -f --display_counts > ${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}.txt

  scil_filter_streamlines_by_length.py ${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}.trk \
              ${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}_${params.minL}.trk \
              --minL ${params.minL} -f --display_counts > ${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}_${params.minL}.txt
  """
}

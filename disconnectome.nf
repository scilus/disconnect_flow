#!/usr/bin/env nextflow

import groovy.json.*

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

if (params.filteredTractograms){
    log.info "Input: $params.filteredTractograms"
    root = file(params.filteredTractograms)
    filteredTractograms = Channel
      .fromFilePairs("$root/**/*trk",
                      size: 1,
                      maxDepth:1,
                      flat: true) {it.simpleName.split('_')}

    filteredTractograms.map{[it].flatten()}.into{find_cc; find_cs; find_ct}
    find_cc.filter{ it[1] == 'CorticoCortical'}.set{corticocortical_precomputed}
    find_cs.filter{ it[1] == 'CorticoStriatal'}.set{corticostriatal_precomputed}
    find_ct.filter{ it[1] == 'CorticoThalamic'}.set{corticothalamic_precomputed}

    corticocortical_precomputed.map{[it[0], it[2], it[3]]}.set{corticocortical_precomputed}
    corticostriatal_precomputed.map{[it[0], it[2], it[3]]}.set{corticostriatal_precomputed}
    corticothalamic_precomputed.map{[it[0], it[2], it[3]]}.set{corticothalamic_precomputed}
}

in_data.into{in_data_cc;
             in_data_ccc;
             in_data_cs;
             in_data_csc;
             in_data_ct;
             in_data_ctc}

if (params.tractograms) {
  log.info "Im HERE"
  in_tractograms.into{in_tractograms_cc;
                      in_tractograms_cs;
                      in_tractograms_ct}
}

indices_cc = params.indices_cc?.tokenize(',')
indices_cs = params.indices_cs?.tokenize(',')
indices_ct = params.indices_ct?.tokenize(',')
sides = params.sides?.tokenize(',')

atlas = Channel.from(file(params.atlasFolder))
atlas.into{in_atlas_1;
           in_atlas_2;
           in_atlas_3;
           in_atlas_4;
           in_atlas_5;
           in_atlas_6;
           in_atlas_7;
           in_atlas_8}

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

if (params.tractograms){
    process filterCorticoCortical{
      cpus 1
      publishDir = params.FilterCorticoCortical_Publish_Dir
      tag = "Filter Cortico-Cortical"

      input:
        set tid, file(tractogram) from in_tractograms_cc
        each idx from indices_cc
        file(atlas_f) from in_atlas_1.first()

      output:
        set tid, idx, "${tid}_CorticoCortical_${idx}.trk" into corticocortical, corticocorticalcommisural optional true

      script:
        """
            scil_remove_invalid_streamlines.py ${tractogram} valid_${tid}_${idx}.trk
            scil_filter_tractogram.py valid_${tid}_${idx}.trk \
              ${tid}_CorticoCortical_${idx}.trk \
              --drawn_roi ${atlas_f}/cortical17/${idx}.nii.gz both_ends include --no_empty
        """
    }
}

if (params.tractograms){
      process filterCorticoStriatal{
        cpus 1
        publishDir = params.FilterCorticoStriatal_Publish_Dir
        tag = "Filter Cortical-Striatal"

        input:
          set tid, file(tractogram) from in_tractograms_cs
          each idx from indices_cs
          file(atlas_f) from in_atlas_2.first()

        output:
          set tid, idx, "${tid}_CorticoStriatal_${idx}.trk" into corticostriatal, corticostriatalcommisural optional true

        script:
          """
              scil_remove_invalid_streamlines.py ${tractogram} valid_${tid}_${idx}.trk
              scil_filter_tractogram.py valid_${tid}_${idx}.trk \
                ${tid}_CorticoStriatal_${idx}.trk \
                --drawn_roi ${atlas_f}/cortical17/${idx}.nii.gz either_end include \
                --drawn_roi ${atlas_f}/striatal17/${idx}.nii.gz either_end include --no_empty
          """
      }
}

if (params.tractograms){
    process filterCorticoThalamic{
      cpus 1
      publishDir = params.FilterCorticoThalamic_Publish_Dir
      tag = "Filter Cortical-Thalamic"


      input:
        set tid, file(tractogram) from in_tractograms_ct
        each idx from indices_ct
        file(atlas_f) from in_atlas_3.first()

      output:
        set tid, idx, "${tid}_CorticoThalamic_${idx}.trk" into corticothalamic, corticothalamiccommisural optional true

      script:
        """
            scil_remove_invalid_streamlines.py ${tractogram} valid_${tid}_${idx}.trk
            scil_filter_tractogram.py valid_${tid}_${idx}.trk \
              ${tid}_CorticoThalamic_${idx}.trk \
              --drawn_roi ${atlas_f}/cortical17/${idx}.nii.gz either_end include \
              --drawn_roi ${atlas_f}/thalamic17/${idx}.nii.gz either_end include --no_empty
        """
    }
}

if (params.filteredTractograms) {
    corticocortical_precomputed.into{corticocortical; corticocorticalcommisural}
    corticostriatal_precomputed.into{corticostriatal; corticostriatalcommisural}
    corticothalamic_precomputed.into{corticothalamic; corticothalamiccommisural}
}

in_data_cc.combine(corticocortical).set{in_data_corticocortical}
in_data_ccc.combine(corticocorticalcommisural).set{in_data_corticocorticalcommisural}
in_data_cs.combine(corticostriatal).set{in_data_corticostriatal}
in_data_csc.combine(corticostriatalcommisural).set{in_data_corticostriatalcommisural}
in_data_ct.combine(corticothalamic).set{in_data_corticothalamic}
in_data_ctc.combine(corticothalamiccommisural).set{in_data_corticothalamiccommisural}

process filterLesionCorticoCortical{
  input:
    set sid, file(lesion), tid, idx, file(tractogram) from in_data_corticocortical
    each side from sides
    file(atlas_f) from in_atlas_4.first()

  output:
    set sid, tid, idx, side, "${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}.txt", "${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}_${params.minL}.txt" into corticocorticalLesion optional true
    set sid, idx into corticocorticalIndices
    file "${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}.trk" optional true
    file "${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}_${params.minL}.trk" optional true
    file "${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}_filter.txt"

  script:
  """
  echo "bdo ${atlas_f}/${side}HemisphereMNI.bdo both_ends include" >> ${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}_filter.txt
  echo "drawn_roi ${lesion} any include" >> ${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}_filter.txt

  scil_filter_tractogram.py ${tractogram} ${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}.trk \
                --filtering_list ${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}_filter.txt \
                -f --no_empty --display_counts > ${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}.txt --no_empty

  if ${params.filterLength} && [ -f "${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}.trk" ]
  then
  scil_filter_streamlines_by_length.py ${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}.trk \
              ${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}_${params.minL}.trk \
              --minL ${params.minL} -f --no_empty --display_counts > ${sid}_${tid}_CorticoCortical_lesion_${idx}_${side}_${params.minL}.txt
  fi
  """
}

process filterLesionCorticoCorticalCommissural{
  input:
    set sid, file(lesion), tid, idx, file(tractogram) from in_data_corticocorticalcommisural
    file(atlas_f) from in_atlas_5.first()

  output:
    set sid, tid, idx, "${sid}_${tid}_CorticoCortical_lesion_${idx}.txt", "${sid}_${tid}_CorticoCortical_lesion_${idx}_${params.minL}.txt" into corticocorticalcommisuralLesion optional true
    set sid, idx into corticocorticalcommisuralIndices
    file "${sid}_${tid}_CorticoCortical_lesion_${idx}.trk" optional true
    file "${sid}_${tid}_CorticoCortical_lesion_${idx}_${params.minL}.trk" optional true
    file "${sid}_${tid}_CorticoCortical_lesion_${idx}_filter.txt"

  script:
  """
  echo "bdo ${atlas_f}/LHemisphereMNI.bdo either_end include" >> ${sid}_${tid}_CorticoCortical_lesion_${idx}_filter.txt
  echo "bdo ${atlas_f}/RHemisphereMNI.bdo either_end include" >> ${sid}_${tid}_CorticoCortical_lesion_${idx}_filter.txt
  echo "drawn_roi ${lesion} any include" >> ${sid}_${tid}_CorticoCortical_lesion_${idx}_filter.txt

  scil_filter_tractogram.py ${tractogram} ${sid}_${tid}_CorticoCortical_lesion_${idx}.trk \
                --filtering_list ${sid}_${tid}_CorticoCortical_lesion_${idx}_filter.txt \
                -f --no_empty --display_counts > ${sid}_${tid}_CorticoCortical_lesion_${idx}.txt

  if ${params.filterLength} && [ -f "${sid}_${tid}_CorticoCortical_lesion_${idx}.trk" ]
  then
  scil_filter_streamlines_by_length.py ${sid}_${tid}_CorticoCortical_lesion_${idx}.trk \
             ${sid}_${tid}_CorticoCortical_lesion_${idx}_${params.minL}.trk \
             --minL ${params.minL} -f --no_empty --display_counts > ${sid}_${tid}_CorticoCortical_lesion_${idx}_${params.minL}.txt
  fi
  """
}

process filterLesionCorticoStriatal{
  input:
    set sid, file(lesion), tid, idx, file(tractogram) from in_data_corticostriatal
    each side from sides
    file(atlas_f) from in_atlas_6.first()

  output:
    set sid, tid, idx, side, "${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}.txt", "${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}_${params.minL}.txt" into corticostriatalLesion optional true
    file "${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}.trk" optional true
    file "${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}_${params.minL}.trk" optional true
    file "${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}_filter.txt"

  script:
  """
  echo "bdo ${atlas_f}/${side}HemisphereMNI.bdo both_ends include" >> ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}_filter.txt
  echo "drawn_roi ${lesion} any include" >> ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}_filter.txt

  scil_filter_tractogram.py ${tractogram} ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}.trk \
                --filtering_list ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}_filter.txt \
                -f --no_empty --display_counts > ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}.txt

  if ${params.filterLength} && [ -f "${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}.trk" ]
  then
    scil_filter_streamlines_by_length.py ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}.trk \
              ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}_${params.minL}.trk \
              --minL ${params.minL} -f --no_empty --display_counts > ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}_${params.minL}.txt
  fi
  """
}

process filterLesionCorticoStriatalCommisural{
  input:
    set sid, file(lesion), tid, idx, file(tractogram) from in_data_corticostriatalcommisural
    each side from sides
    file(atlas_f) from in_atlas_7.first()

  output:
    set sid, tid, idx, "${sid}_${tid}_CorticoStriatal_lesion_${idx}_LR.txt", "${sid}_${tid}_CorticoStriatal_lesion_${idx}_LR_${params.minL}.txt" into corticostriatalLRcommisuralLesion optional true
    set sid, tid, idx, "${sid}_${tid}_CorticoStriatal_lesion_${idx}_RL.txt", "${sid}_${tid}_CorticoStriatal_lesion_${idx}_RL_${params.minL}.txt" into corticostriatalRLcommisuralLesion optional true
    file "${sid}_${tid}_CorticoStriatal_lesion_${idx}_LR.trk" optional true
    file "${sid}_${tid}_CorticoStriatal_lesion_${idx}_RL.trk" optional true
    file "${sid}_${tid}_CorticoStriatal_lesion_${idx}_LR_${params.minL}.trk" optional true
    file "${sid}_${tid}_CorticoStriatal_lesion_${idx}_RL_${params.minL}.trk" optional true
    file "${sid}_${tid}_CorticoStriatal_lesion_${idx}_LR_filter.txt" optional true
    file "${sid}_${tid}_CorticoStriatal_lesion_${idx}_RL_filter.txt" optional true

  script:
  """
  if [ "${side}" == "L" ]
  then
    opside="R"
  else
    opside="L"
  fi

  echo "drawn_roi ${atlas_f}/cortical17/${idx}_${side}.nii.gz either_end include" >> ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}\${opside}_filter.txt
  echo "drawn_roi ${atlas_f}/striatal17/${idx}_\${opside}.nii.gz either_end include" >> ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}\${opside}_filter.txt
  echo "drawn_roi ${lesion} any include" >>  ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}\${opside}_filter.txt

  if [ -f "${atlas_f}/striatal17/${idx}_\${opside}.nii.gz" ]
  then
  scil_filter_tractogram.py ${tractogram} ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}\${opside}.trk \
              --filtering_list ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}\${opside}_filter.txt \
              -f --no_empty --display_counts > ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}\${opside}.txt
  fi

  if ${params.filterLength} && [ -f "${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}\${opside}.trk" ]
  then
  scil_filter_streamlines_by_length.py ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}\${opside}.trk \
            ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}\${opside}_${params.minL}.trk \
           --minL ${params.minL} -f --no_empty --display_counts > ${sid}_${tid}_CorticoStriatal_lesion_${idx}_${side}\${opside}_${params.minL}.txt
  fi
  """
}

process filterLesionCorticoThalamic{
  input:
    set sid, file(lesion), tid, idx, file(tractogram) from in_data_corticothalamic
    each side from sides
    file(atlas_f) from in_atlas_8.first()

  output:
    set sid, tid, idx, side, "${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}.txt", "${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}_${params.minL}.txt" into corticothalamicLesion optional true
    file "${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}.trk" optional true
    file "${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}_${params.minL}.trk" optional true
    file "${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}_filter.txt"

  script:
  """
  echo "bdo ${atlas_f}/${side}HemisphereMNI.bdo both_ends include" >> ${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}_filter.txt
  echo "drawn_roi ${lesion} any include" >> ${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}_filter.txt

  scil_filter_tractogram.py ${tractogram} ${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}.trk \
                --filtering_list ${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}_filter.txt \
                -f --no_empty --display_counts > ${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}.txt

  if ${params.filterLength} && [ -f "${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}.trk" ]
  then
  scil_filter_streamlines_by_length.py ${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}.trk \
             ${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}_${params.minL}.trk \
             --minL ${params.minL} -f --no_empty --display_counts > ${sid}_${tid}_CorticoThalamic_lesion_${idx}_${side}_${params.minL}.txt
  fi
  """
}

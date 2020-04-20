DisconetsFlow:

Example command line (Build filtered tractograms really long):
```
nextflow run disconnectome.nf \
        --root disconects_dataset/ \
        --atlasFolder [FullPathTo]/disconect_atlas/ \
        --tractograms [FullPathTo]/disconect_tractograms/ -resume
```

or

```
nextflow run disconnectome.nf \
        --root disconects_dataset/ \
        --filteredTractograms [FullPathTo]/disconets_filtered_tractograms/ \
        --atlasFolder [FullPathTo]/disconect_atlas_short_version/ \
        --output_dir [FullPathTo]/output/ \
        --resume
``


```
disconects_dataset
.
├── sub-ROUSSEAU
│   ├── cavity_ROUSSEAU_ANTS.nii.gz
│   └── cavity_ROUSSEAU_FSL.nii.gz
└── sub-ROUSSELIE
    ├── cavity_ROUSSELIE_ANTS.nii.gz
    └── cavity_ROUSSELIE_FSL.nii.gz
 ```

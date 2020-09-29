Disconnectoflow:

Here are some [data](https://drive.google.com/file/d/1VX1dUgGBYr6JJ8_osHGmnsVqO_CJVl7y/view?usp=sharing)

Example command line (Build filtered tractograms really long):
```
nextflow run disconnectome.nf \
        --root dataset/ \
        --atlasFolder [FullPathTo]/atlas/ \
        --tractograms [FullPathTo]/tractograms/ \
        -with-singularity scilpy-1.0.0-rc1.sif \
        resume
```

or

```
nextflow run disconnectome.nf \
        --root dataset/ \
        --filteredTractograms [FullPathTo]/filtered_tractograms/ \
        --atlasFolder [FullPathTo]/atlas/ \
        --output_dir [FullPathTo]/output/ \
        -with-singularity scilpy-1.0.0-rc1.sif \
        --resume
```


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

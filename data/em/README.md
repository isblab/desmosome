# EM map


### Source
The main cryo-EM density map used as a restraint is [EMDB-1703](https://www.ebi.ac.uk/pdbe/entry/emdb/EMD-1703) and the depositing paper being [Al-Amoudi et al, PNAS 2011](doi.org/10.1073/pnas.1019469108)\[1\]. 

### Deposition Details
The deposition contains the main density, `maps/emd_1703.map` and two additional masked densities `maps/emd_1703_msk.map` and `maps/emd_1703_msk_1.map`; which respectively correspond to a symmetrized tomograph reconstruction and the denoised cylindrical section of the original (asymmetric) density. For all the modeling, we are using the `emd_1703_msk_1.map`. The header xml files establish the `pixelspacing` of the three maps to be `6A` and the recommended contour level is `0.006`. 

### Processing 
The EM map `voxelsize` is set to `6A` before all downstream processing in chimera and this *resized* density is stored separately as `maps/em_1703_msk_1_voxel6.mrc` (command to change the voxelsize: `volume #model-id voxelsize 6`). [Segger](https://github.com/gregdp/segger) is the default [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/) plugin to perform watershed-filtering based segmentation\[2\]. We segmented the three major "layers" (PG-DP, PKP and the plasma membrane) and this segmentation is saved as `segments/layers.seg`. The isolated PKP layer density (masked by the pkp layer segmentation) is saved as `maps/pkp_separate.mrc` and a corresponding segmentation of the the pkp-only density to identify the individual density segments (blobs) is saved as `segments/pkp_individual.seg`. All cross-correlation calculations and fitting of the PKP molecules for preliminary analysis can be done using `fit to segments` toolbox in Chimera after generating a `molmap` from the PKP molecule at `32A` resolution (the reported overall resolution of the deposition) and a `gridSpacing` of `6` with a contour level of `0.03`.




\[1\]: Al-Amoudi, A., Castaño-Diez, D., Devos, D. P., Russell, R. B., Johnson, G. T., & Frangakis, A. S. (2011). The three-dimensional molecular structure of the desmosomal plaque. Proceedings of the National Academy of Sciences, 108(16), 6480-6485

\[2\]: Pintilie, G. D., Zhang, J., Goddard, T. D., Chiu, W., & Gossard, D. C. (2010). Quantitative analysis of cryo-EM density map segmentation by watershed and scale-space filtering, and fitting of structures by alignment to regions. Journal of structural biology, 170(3), 427-438

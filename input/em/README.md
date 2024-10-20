# cryo-tomogram/GMM files


### Source
The main cryo-tomogram used as a restraint is [EMDB-1703](https://www.ebi.ac.uk/pdbe/entry/emdb/EMD-1703) and the depositing paper being [Al-Amoudi et al, PNAS 2011](doi.org/10.1073/pnas.1019469108)\[1\]. 

### Deposition Details
The deposition contains the main density, `maps/emd_1703.map` and two additional masked densities `maps/emd_1703_msk.map` and `maps/emd_1703_msk_1.map`; which respectively correspond to a symmetrized tomograph reconstruction and the denoised cylindrical section of the original (asymmetric) density. For all the modeling, we are using the `emd_1703_msk_1.map` after processing (see below). The header xml files establish the `pixelspacing` of the three maps to be `6A` and the author-recommended contour level to be `0.006`.

### Restraint summary
The restraint requires a GMM approximation of the tomogram and a correspoding model-generated GMM (forward model) the overlap of which scores the model (i.e. determines the probability of the model during MCMC sampling).

### Processing 
1. The EM map `voxelsize` is set to `6A` before all downstream processing in chimera and this *resized* full unmodified-density is stored separately as `maps/em_1703_msk_1_voxel6.mrc` (command to change the voxelsize: `volume #model-id voxelsize 6`). [Segger](https://github.com/gregdp/segger) is the default [UCSF Chimera](https://www.cgl.ucsf.edu/chimera/) plugin to perform watershed-filtering based segmentation\[2\]. We segmented the three major "layers" (PG-DP, PKP and the plasma membrane) and this segmentation is saved as `segments/layers.seg`. A plasma-membrane-removed version is saved as `maps/full_PMremoved.mrc` (by masking the map with the non-PM layer segmentations). Furthermore, the isolated PKP layer density (masked by the pkp layer segmentation) is saved as `maps/pkp_separate.mrc` and a corresponding segmentation of the the pkp-only density to identify the individual density segments (blobs) is saved as `segments/pkp_individual.seg` (used for calculating the centers of the blobs below). Likewise, the isolated PG-layer density is saved as `maps/pgdp_separate.mrc`.
2. For the GMM calculations and for obtaining the corresponding `gmm` text files for input into IMP, the script `gmm_selection.py` runs the IMP script `IMP.isd.create_gmm` for multiple values of the `n_centers` parameter (number of GMM gaussians used) and saves all the corresponding `gmm.txt` (to be found in `gmm/`) and the `.mrc` files (not uploaded due to size); see `gmm_selection.bat`. A chimera script (`compare_gmm.com`) compares all the generated `.mrc` files against the original to get the cross-correlation. The final selected `gmm` file (`gmm/selection/full_PMremoved_gmm_30.txt`) is the one with the minimum `n_centers` giving a cross-correlation around mean greater than 0.9 (alternatively, the cross-correlation around zero greater than 0.95). (see `gmm/chimera_reply_log.txt` for the `.mrc` comparisons)
3. The `bild_files` directory contains [Chimera BILD files](https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/bild.html) for easy visualisation of the axes/various cutoffs/bounding boxes used in relation to the EM density. See the in-file comments for more details.
4. The `measure_centers_and_fit` chimera command file calculates the centers of all the separate PKP blobs (using the segmentation corresponding to each blob stored separately in a subfolder `sub_pkp_segments` in `segments`. Each of these segements is the same as in `segments/pkp_individual.seg`) and the fit of the PKP PDB (`1XM9`) to each of the blobs. The partial reply log with all the relevant output is given in the file `measure_centers_and_fit_replylog.txt`
5. `docking.com` is a chimera command file to find the best positions/orientations of `1XM9` in 3 of the 7 PKP blobs starting from a set of initial manual placements (`docking_start_x_1xm9.pdb`) and the final placements as computed and saved by the command file are saved as `docking_x_gpkp.pdb` for use downstream

**NOTE**: Step 2 is done for the full density with the plasma membrane segment removed (`full_PMremoved`), the isolated PKP-layer density (`pkp_separate`) and the isolated PG/PG-DP-layer density (`pgdp_separate`). The `.bat` and `.com` scripts include commands for all the three steps sequentially.

**NOTE**: The field linking to the path of the original map in the header of the `gmm.txt` may show broken paths depending on the platform due to the script used. Manually fixing the single field if needed (modified for the ones in `gmm/` according to downstream use) should avoid any subsequent errors.

\[1\]: Al-Amoudi, A., Castaño-Diez, D., Devos, D. P., Russell, R. B., Johnson, G. T., & Frangakis, A. S. (2011). The three-dimensional molecular structure of the desmosomal plaque. Proceedings of the National Academy of Sciences, 108(16), 6480-6485

\[2\]: Pintilie, G. D., Zhang, J., Goddard, T. D., Chiu, W., & Gossard, D. C. (2010). Quantitative analysis of cryo-EM density map segmentation by watershed and scale-space filtering, and fitting of structures by alignment to regions. Journal of structural biology, 170(3), 427-438

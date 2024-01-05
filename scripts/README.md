# Modeling, Analysis and Helper Scripts

### Layout

This contains the scripts of the main modeling setup (stratified epithelial ODP with PKP1) for which the results are presented in the paper. It also contains the scripts for several auxilary runs including other tissue-specific ODPs, Stoichiometry Runs, Autocorrelation Runs and running a Multimeric version of Alphafold-2.

1. **modeling**: Modeling scripts for all the runs
2. **analysis**: Analysis and plotting scripts (including figure papers), Alphafold scripts
3. **utils_helpers_archive**: Miscellaneous

----

### Custom Module

`custom_imp_module`: Contains the custom IMP module, "Desmosome", that contains C++ implementations of certain restraints (Cylinder and Immuno-EM). These files are located in `src` as `CylinderLocalizationRestraint` and `SingleAxisMinGaussianRestraint` (the other files are either used in auxiliary runs or not used for modeling). The `tests` contains some unit-tests for the same implementations. Compiling IMP with this module is necessary to run the subsequent scripts. See the instructions for compilation see the [official guide](https://integrativemodeling.org/nightly/doc/manual/installation.html). Also see `auto_compile.sh` and `install_all.sh` in `utils_helpers_archive`.

----

### Modeling ODPs

Modeling each ODP involved parallelizing 45 runs across several servers using [GNU Parallel](link) (for SSH based parallelizing) and MPI for parallelizing the multiple replicas (8) in each run. Each ODP-specific folder (`main_run`, `other_runs/basal_pkp3_dc3`, `other_runs/upper_pkp3`) contains 3 files:
0. `all_run.sh`: To be run from the main server; uses SSH to run `server_run.sh` from each of the other servers. The data from each server is collected manually and collated at the end of sampling for analysis.
1. `representation_sampling.py`: The main modeling script that `server_run.sh` runs for each run for each replica. All the specific parameters relevant to the modeling are defined here.
2. `server_run.sh`: This is called by `all_run.sh` in each server to run one individual run (containing 8 replicas)

**Note**: `representation_sampling.py` expects the input data to be in a different directory structure than the `input` directory in this repository. 

The expected structure is as follows: all the required fasta files and the PDB files inside a directory called `data` in the working directory (from which the script is called), with a subdirectory `data/gmm` that contains the `.mrc` files and the GMM `.txt` files corresponding to the tomogram GMMs (for the specific densities see README in `data/em`). An example setup is given under `inputs/data` which corresponds to the `data` folder for the main run. Minor modifications of the header in `FASTA` files or the GMM files might be needed for other runs but is already done for the example `data` directory (see the note [here](https://github.com/isblab/desmosome/blob/main/input/em/README.md))

The results are in the `results` folder and Zenodo (see the [results README](https://github.com/isblab/desmosome/blob/main/results/README.md)).

----

### Analysis of ODPs
Analysis involved running sampling exhaustiveness, plotting, contact-map generation and all the processing after sampling the models. Each ODP-specific folder contains several common files:
1. `extract_all.sh`: The set of commands to reproduce the analysis results. This requires all the "output" folders, corresponding to the 45 runs sampled above, to be present in the directory this script is run from. The output directory, `analysis_output` has to be created manually. Ensure that the paths used in the script (for example, the path to the IMP-installation setup script, or to the `.mrc` files) are correct.
2. `analysis.py`: This implements HDBSCAN with `n_skip=2` (based on autocorrelation decay, see below) to create a set of `selected_models`
3. `variable_filter.py`: Filters the `selected_models` of the largest cluster identified above to create a subset, `gsm` (good scoring models), by filtering all models with restraint-specific and the total score worse than `mean + c * SD` where `c` is calculated to maximize computational efficiency
4. `extract_models.py`: Extracts the filtered models into RMF files
5. `extract_sampcon.py`: This is called after running sampcon (see [imp-sampcon](https://github.com/salilab/imp-sampcon)) to create an RMF file of the largest cluster identified during clustering in sampcon.
6. Helper scripts needed for the above: `ambiguity.txt` defines the ambiguity of protein copies, `density.txt` defines the density ranges for sampcon, `mrc_parser.py` has helper functions for cross-correlation of densities (only in the main run directory), `cm.py` contains helper functions for contact map computation (only in the main run directory)
7. `housekeeping.py`: (Output not used in the paper) Generates several run-specific housekeeping plots (to see the equilibration of restraint and total score) as well as the general distributions of the restraint/total scores. Also checks for the MCMC acceptance ratios being in a target range.

Additionally, the main run directory contains some extra files:
1. `sampcon_fit_to_data.py`: Generates the contact maps, plots for fit to data used in modeling.
2. `validation_fit_to_data.py`: Same as above for the fit to data not used in modeling.
3. `cuttoff_compute.py`: Generates the table of cutoffs for different thresholds that can be used in creating the contact maps. See also `utils_helpers_archive/rectangle_overlap.py` to see how the regions of a lower threshold (>0.2 in the paper) are separated from those at a higher threshold (>0.25 in the paper)
4. `compare_sample_A_B.py`: Compares the Sample A/B Localization Densities using an implementation of cross-correlation as part of statistical tests for sampling exhaustiveness
5. `chimera_density.py`: The Chimera script that can load the densities and apply the visualization thresholds used in the paper. **The same thresholds are used for all the three runs.**. 
6. `sikora_dsc_dsg.py`: Generates one of the subfigures and computes fit to the data from Sikora et al.
7. `align_pdb_to_rmf.py`: Aligns and outputs PDBs to the corresponding rigid body locations in a given RMF (eg: cluster center model)

----

### Stoichiometry and Autocorrelation Runs

Additionally, separate runs were done to calculate the number of protein copies of PG and DP to fit the tomogram (stoichiometry runs) as well as to determine the rate of autocorrelation decay to filter the models during analysis for the main run (autocorrelation runs). The corresponding modeling files and analysis files (including the ones needed to generate Fig S1) are found in the `modeling/other_runs/{name}` directory. The scripts in the autocorrelation directory may have some code overlap with other files.

----

### Alphafold

We employed [AF2-Multimer](https://github.com/deepmind/alphafold) on all protein-pairs of interest. The best predicted PDB for each protein-pair as well as the input sequences are found in protein-pair specific directories in `analysis/alphafold` (without the `.pkl` files due to space constraints). The script `alphafold_analyze.py` generates the pLDDT/PAE figures (not used in paper, needs the `.pkl` files), contact maps used in the paper and the interface lists. The script accepts a folder containing the best performing model (`ranked_0.pdb`), the `ranking_debug.json` file to identify the corresponding `.pkl` file, and the corresponding `.pkl` file. Generated PAE/PLDDT figures and the contact lists are uploaded (`results/alphafold_analysis_output`). Note that the PAE matrix returned by AF2-Multimer is not strictly symmetric, hence, the contacts between chain A -> B and B -> A can differ. For the paper, we use the former, even though the overall disagreement between the two is minimal.

----

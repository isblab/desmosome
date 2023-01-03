import sys
sys.path = ['/home/satwik/pmi_analysis/pyext/src'] + sys.path
import os
from analysis_trajectories import AnalysisTrajectories as AT

# sys.argv -> path, number of cores, number of runs, n-skip (default 1)
n_skip = 1
if len(sys.argv) == 5:
    n_skip = int(sys.argv[4])
n_runs = int(sys.argv[3])
n_processes = int(sys.argv[2])
if not os.path.isdir(f'{sys.argv[1]}/analysis_output'):
    os.mkdir(f'{sys.argv[1]}/analysis_output')

at_obj = AT([f'{sys.argv[1]}/output{i}' for i in range(n_runs)], dir_name='output', nproc=n_processes,
            analysis_dir=f'{sys.argv[1]}/analysis_output', burn_in_fraction=0.1, nskip=n_skip, detect_equilibration=True)
at_obj.set_analyze_Connectivity_restraint()
at_obj.set_analyze_Excluded_volume_restraint()
at_obj.set_analyze_score_only_restraint(handle='GaussianEMRestraint', short_name='EM', do_sum=True)
at_obj.set_analyze_score_only_restraint(handle='SingleAxisMinGaussianRestraint', short_name='SAMGR', do_sum=True)
at_obj.set_analyze_score_only_restraint(handle='MinimumPairDistanceBindingRestraint', short_name='MPDBR', do_sum=True)

at_obj.read_stat_files()
at_obj.write_models_info()
at_obj.hdbscan_clustering(['EV_sum', 'CR_sum', 'EM_sum', 'SAMGR_sum', 'MPDBR_sum'], min_cluster_size=100, min_samples=5,
                          skip=n_skip)


for ((j = 1; j < 31; j++))
do
  mkdir output$j
done
echo Made all

../../imp-custom/build/setup_environment.sh python representation_sampling.py --topology_only -DP_NUMBER=7 -PG_NUMBER=7
../../imp-custom/build/setup_environment.sh python -u batch_run.py run_stoichiometry_2_3_4.txt 0
../../imp-custom/build/setup_environment.sh python -u batch_run.py run_stoichiometry_5_6_7.txt 15
python plot_all_results.py
python plot_all_results_2.py
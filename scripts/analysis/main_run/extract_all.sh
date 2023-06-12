#!/bin/bash
python -u housekeeping.py . plots 45 1> housekeeping_report_1.txt 2> housekeeping_report_2.txt
~/imp-custom/build/setup_environment.sh python -u analysis.py . 45 45 2 1>analysis_report_1.txt 2>analysis_report_2.txt
~/imp-custom/build/setup_environment.sh python variable_filter.py -c 3 -n 25000 -ss 0.01 -lc -3 -hc 4 -g ./analysis_output/
~/imp-custom/build/setup_environment.sh python extract_models.py . 45 45 3
cd analysis_output
mkdir sampcon_output
cd sampcon_output
~/imp-custom/build/setup_environment.sh python ~/imp-custom/imp/modules/sampcon/pyext/src/exhaust.py --sysname desmosome -ra ../A_gsm_clust3.rmf3 -rb ../B_gsm_clust3.rmf3 -sa ../A_gsm_clust3.txt -sb ../B_gsm_clust3.txt -gp -g 1.0 -m cpu_omp -pr -c 5 -amb ../../ambiguity.txt -d ../../density.txt 1>sampcon_report_1.txt
~/imp-custom/build/setup_environment.sh python ../../extract_sampcon.py sampcon_0_extracted.rmf3 ../A_gsm_clust3.rmf3 cluster.0.sample_A.txt ../B_gsm_clust3.rmf3 cluster.0.sample_B.txt
mkdir plots
~/imp-custom/build/setup_environment.sh python ../../sampcon_fit_to_data.py sampcon_0_extracted.rmf3 plots ../../original_mrc_file.mrc . ../../original_mrc_file_pkp.mrc
python ../../cutoff_compute.py contact_maps contact_maps
~/imp-custom/build/setup_environment.sh python ../../validation_fit_to_data.py plots
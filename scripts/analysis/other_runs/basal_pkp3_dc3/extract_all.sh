#!/bin/bash
python -u housekeeping.py . plots 45 1> housekeeping_report_1.txt 2> housekeeping_report_2.txt
~/imp-custom/build/setup_environment.sh python -u analysis.py . 45 45 2 1>analysis_report_1.txt 2>analysis_report_2.txt
~/imp-custom/build/setup_environment.sh python variable_filter.py -c 1 -n 25000 -ss 0.01 -lc -3 -hc 4 -g ./analysis_output/
~/imp-custom/build/setup_environment.sh python extract_models.py . 45 45 1
rm -r output*
cd analysis_output
mkdir sampcon_output
cd sampcon_output
~/imp-custom/build/setup_environment.sh python ~/imp-custom/imp/modules/sampcon/pyext/src/exhaust.py --sysname desmosome -ra ../A_gsm_clust1.rmf3 -rb ../B_gsm_clust1.rmf3 -sa ../A_gsm_clust1.txt -sb ../B_gsm_clust1.txt -gp -g 1.0 -m cpu_omp -pr -c 5 -amb ../../ambiguity.txt -d ../../density.txt 1>sampcon_report_1.txt
~/imp-custom/build/setup_environment.sh python ../../extract_sampcon.py sampcon_0_extracted.rmf3 ../A_gsm_clust1.rmf3 cluster.0.sample_A.txt ../B_gsm_clust1.rmf3 cluster.0.sample_B.txt
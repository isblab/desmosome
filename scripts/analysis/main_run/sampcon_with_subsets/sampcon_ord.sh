#!/bin/bash
~/imp-custom/build/setup_environment.sh python ~/imp-custom/imp/modules/sampcon/pyext/src/exhaust.py --sysname desmosome -ra A_gsm_clust3.rmf3 -rb B_gsm_clust3.rmf3 -sa A_gsm_clust3.txt -sb B_gsm_clust3.txt -gp -g 1.0 -m cpu_omp -pr -c 5 -amb ambiguity_ord.txt -sn selection_ord.txt -d selection_ord.txt 1>sampcon_report_1.txt



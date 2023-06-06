module load mpi/openmpi-x86_64
seq 0 9 | parallel --progress -P 10 '~/imp-custom/build/setup_environment.sh python -u representation_sampling.py all_outputs/output{} 1> logs/run{}_stream2.txt'
~/imp-custom/build/setup_environment.sh python autocorrel_and_cm.py all_outputs 0:500:5 plots 10 20000
cp -R all_outputs temporary
~/imp-custom/build/setup_environment.sh python analysis.py temporary 10 10
~/imp-custom/build/setup_environment.sh python extract_models.py temporary 10 10 0
~/imp-custom/build/setup_environment.sh python extract_models.py temporary 10 10 1
mv temporary/analysis_output ./analysis_output
rm -r temporary
~/imp-custom/build/setup_environment.sh python ~/imp-custom/imp/modules/sampcon/pyext/src/exhaust.py desmosome -ra analysis_output/A_models_clust0.rmf3 -rb analysis_output/B_models_clust0.rmf3 -sa analysis_output/A_models_clust0.txt -sb analysis_output/B_models_clust0.txt -gp -g 1.0 -m cpu_omp -c 10 -amb ambiguity.txt -d density.txt --align

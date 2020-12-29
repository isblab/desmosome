ECHO This script will output log files everywhere possible for debugging
ECHO =========================================================================
:: Getting residue number information from the PDB files
python scripts/pdb_info_extract.py input_data/pdb_files/1xm9.pdb output_data 1>> script_debug_run.log 2>&1
python scripts/pdb_info_extract.py input_data/pdb_files/3ifq.pdb output_data 1>> script_debug_run.log 2>&1
python scripts/pdb_info_extract.py input_data/pdb_files/3r6n.pdb output_data 1>> script_debug_run.log 2>&1
:: Getting MODELLER ready alignments (and residue number details) from .hhr files
python scripts/hhr_extract_alignment.py input_data/hhr_files/pkp2.hhr 1XM9_A input_data/pdb_files/1xm9.pdb A output_data/modeller_alignment_files pkp2 1>> script_debug_run.log 2>&1
python scripts/hhr_extract_alignment.py input_data/hhr_files/pkp3.hhr 1XM9_A input_data/pdb_files/1xm9.pdb A output_data/modeller_alignment_files pkp3 1>> script_debug_run.log 2>&1
python scripts/hhr_extract_alignment.py input_data/hhr_files/dsg3.hhr 3IFQ_D input_data/pdb_files/3ifq.pdb D output_data/modeller_alignment_files dsg3 1>> script_debug_run.log 2>&1
python scripts/hhr_extract_alignment.py input_data/hhr_files/dsg1.hhr 3IFQ_D input_data/pdb_files/3ifq.pdb D output_data/modeller_alignment_files dsg1 1>> script_debug_run.log 2>&1
python scripts/hhr_extract_alignment.py input_data/hhr_files/dsc3.hhr 3IFQ_D input_data/pdb_files/3ifq.pdb D output_data/modeller_alignment_files dsc3 1>> script_debug_run.log 2>&1
python scripts/hhr_extract_alignment.py input_data/hhr_files/dsc1.hhr 3IFQ_D input_data/pdb_files/3ifq.pdb D output_data/modeller_alignment_files dsc1 1>> script_debug_run.log 2>&1
python scripts/hhr_extract_alignment.py input_data/hhr_files/dsg3.hhr 3IFQ_D input_data/pdb_files/3ifq.pdb C output_data/modeller_alignment_files/alternative dsg3 1>> script_debug_run.log 2>&1
python scripts/hhr_extract_alignment.py input_data/hhr_files/dsg1.hhr 3IFQ_D input_data/pdb_files/3ifq.pdb C output_data/modeller_alignment_files/alternative dsg1 1>> script_debug_run.log 2>&1
python scripts/hhr_extract_alignment.py input_data/hhr_files/dsc3.hhr 3IFQ_D input_data/pdb_files/3ifq.pdb C output_data/modeller_alignment_files/alternative dsc3 1>> script_debug_run.log 2>&1
python scripts/hhr_extract_alignment.py input_data/hhr_files/dsc1.hhr 3IFQ_D input_data/pdb_files/3ifq.pdb C output_data/modeller_alignment_files/alternative dsc1 1>> script_debug_run.log 2>&1
python scripts/fasta_extract_alignment.py input_data/pdb_files/1xm9.pdb A output_data/modeller_alignment_files pkp1 1>> script_debug_run.log 2>&1
:: Combine alignments to get final alignment files for multi-chain models
python scripts/collate_alignments.py PDB:input_data/pdb_files/3ifq.pdb:A,PDB:input_data/pdb_files/3ifq.pdb:B,output_data/modeller_alignment_files/alternative/alignment_MODELLER_dsg3-3ifq.ali,output_data/modeller_alignment_files/alignment_MODELLER_dsg3-3ifq.ali 3ifq,pg-dsg3-pg-dsg3,D output_data/modeller_alignment_files/combined_3ifq_dsg3.ali 1>> script_debug_run.log 2>&1
python scripts/collate_alignments.py PDB:input_data/pdb_files/3ifq.pdb:A,PDB:input_data/pdb_files/3ifq.pdb:B,output_data/modeller_alignment_files/alternative/alignment_MODELLER_dsg1-3ifq.ali,output_data/modeller_alignment_files/alignment_MODELLER_dsg1-3ifq.ali 3ifq,pg-dsg1-pg-dsg1,D output_data/modeller_alignment_files/combined_3ifq_dsg1.ali 1>> script_debug_run.log 2>&1
python scripts/collate_alignments.py PDB:input_data/pdb_files/3ifq.pdb:A,PDB:input_data/pdb_files/3ifq.pdb:B,output_data/modeller_alignment_files/alternative/alignment_MODELLER_dsc3-3ifq.ali,output_data/modeller_alignment_files/alignment_MODELLER_dsc3-3ifq.ali 3ifq,pg-dsc3-pg-dsc3,D output_data/modeller_alignment_files/combined_3ifq_dsc3.ali 1>> script_debug_run.log 2>&1
python scripts/collate_alignments.py PDB:input_data/pdb_files/3ifq.pdb:A,PDB:input_data/pdb_files/3ifq.pdb:B,output_data/modeller_alignment_files/alternative/alignment_MODELLER_dsc1-3ifq.ali,output_data/modeller_alignment_files/alignment_MODELLER_dsc1-3ifq.ali 3ifq,pg-dsc1-pg-dsc1,D output_data/modeller_alignment_files/combined_3ifq_dsc1.ali 1>> script_debug_run.log 2>&1
:: Run MODELLER
PAUSE
python scripts/run_modeller_automodel.py output_data/modeller_alignment_files/alignment_MODELLER_pkp1-1xm9.ali 1xm9 pkp1 output_data/modeller_run_logs/pkp1_1xm9 input_data/pdb_files 1> output_data/modeller_run_logs/pkp1_1xm9/run.log 2>&1
python scripts/run_modeller_automodel.py output_data/modeller_alignment_files/alignment_MODELLER_pkp3-1xm9.ali 1xm9 pkp3 output_data/modeller_run_logs/pkp3_1xm9 input_data/pdb_files 1> output_data/modeller_run_logs/pkp3_1xm9/run.log 2>&1
python scripts/run_modeller_loopmodel.py output_data/modeller_alignment_files/alignment_MODELLER_pkp2-1xm9.ali 1xm9 pkp2 output_data/modeller_run_logs/pkp2_1xm9 input_data/pdb_files 1> output_data/modeller_run_logs/pkp2_1xm9/run.log 2>&1
python scripts/run_modeller_automodel.py output_data/modeller_alignment_files/combined_3ifq_dsg3.ali 3ifq pg-dsg3-pg-dsg3 output_data/modeller_run_logs/combined_3ifq_dsg3 input_data/pdb_files 1> output_data/modeller_run_logs/combined_3ifq_dsg3/run.log 2>&1
python scripts/run_modeller_automodel.py output_data/modeller_alignment_files/combined_3ifq_dsg1.ali 3ifq pg-dsg1-pg-dsg1 output_data/modeller_run_logs/combined_3ifq_dsg1 input_data/pdb_files 1> output_data/modeller_run_logs/combined_3ifq_dsg1/run.log 2>&1
python scripts/run_modeller_automodel.py output_data/modeller_alignment_files/combined_3ifq_dsc3.ali 3ifq pg-dsc3-pg-dsc3 output_data/modeller_run_logs/combined_3ifq_dsc3 input_data/pdb_files 1> output_data/modeller_run_logs/combined_3ifq_dsc3/run.log 2>&1
python scripts/run_modeller_automodel.py output_data/modeller_alignment_files/combined_3ifq_dsc1.ali 3ifq pg-dsc1-pg-dsc1 output_data/modeller_run_logs/combined_3ifq_dsc1 input_data/pdb_files 1> output_data/modeller_run_logs/combined_3ifq_dsc1/run.log 2>&1
:: Superimpose stuff (if necessary)
:: Collate PDB files (if necessary)
:: Done!
PAUSE
import sys
sys.path = ['/home/satwik/pmi_analysis/pyext/src'] + sys.path
import IMP
import IMP.rmf
import RMF
import tqdm
import pandas as pd
import numpy as np

# sys.argv -> path, number of cores, number of runs, cluster number
n_runs = int(sys.argv[3])
nproc = int(sys.argv[2])
c = int(sys.argv[4])
analys_dir = f'{sys.argv[1]}/analysis_output'

HA = pd.read_csv(f'{sys.argv[1]}/analysis_output/good_scoring_models_A_cluster' + str(c) + '_detailed.csv')
HB = pd.read_csv(f'{sys.argv[1]}/analysis_output/good_scoring_models_B_cluster' + str(c) + '_detailed.csv')

rmf_file_out_A = 'A_gsm_clust' + str(c) + '.rmf3'
rmf_file_out_B = 'B_gsm_clust' + str(c) + '.rmf3'

for df, rmf_out in zip([HA, HB], [rmf_file_out_A, rmf_file_out_B]):
    scores = []
    scorefile = rmf_out.split('.')[0]
    scorefile = f'{sys.argv[1]}/analysis_output/{scorefile}.txt'
    row1 = df.iloc[0]
    rmf_file = '/'.join([x for x in (row1.rmf3_file).split('/')[1:] if x])  # remove the leading directory name
    rmf_file = f'{sys.argv[1]}/{rmf_file}'
    m = IMP.Model()
    f = RMF.open_rmf_file_read_only(rmf_file)
    h0 = IMP.rmf.create_hierarchies(f, m)[0]
    fh_out = RMF.create_rmf_file(f'{analys_dir}/{rmf_out}')
    IMP.rmf.add_hierarchy(fh_out, h0)
    del f

    all_rmfs = df['rmf3_file'].unique()
    all_rmfs2 = ['/'.join([y for y in x.split('/')[1:] if y]) for x in all_rmfs]
    
    print(df['rmf3_file'].value_counts())
    tq = tqdm.tqdm(total=len(df), desc='total', smoothing=0)
    for i in range(len(all_rmfs)):
        rep = f'{sys.argv[1]}/{all_rmfs2[i]}'
        rep_df = df[df['rmf3_file'] == all_rmfs[i]]
        f = RMF.open_rmf_file_read_only(rep)
        IMP.rmf.link_hierarchies(f, [h0])
        tqq = tqdm.tqdm(total=len(rep_df), desc='specific_rmf', smoothing=0, leave=False)
        for row_id, row in rep_df.iterrows():
            fr_rmf = int(row.rmf_frame_index)
            IMP.rmf.load_frame(f, RMF.FrameID(fr_rmf))
            IMP.rmf.save_frame(fh_out, str(i))
            scores.append(row.Total_Score)
            tqq.update()
            tq.update()
        tqq.close()
        del f
    tq.close()
    del fh_out
    np.savetxt(scorefile, np.array(scores))



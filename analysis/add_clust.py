import numpy as np
from Bio.PDB import PDBParser
from sklearn.cluster import DBSCAN
import os
METAL_IONS = {
    'ZN', 'CA', 'MG', 'MN', 'FE', 'CU', 'K', 'NA', 'CD', 'CO', 'NI', 'HG', 'PB', 'SR', 'BA'
}
tmp_file_path = '/mnt/data/public/AF3_benchmark/tmp'
method_list = ['ProteniX']
"""
peptide_type_list = [
    'SAA_seq',
    'USAA_ptm',
    'SAA_smiles',
    'USAA_smiles',
    'CYC_smiles',
    'CYC_ptm',
    'CYC_ptm_bond'
]
"""
peptide_type_list = [

    'CYC_ptm_bond'
]

def get_central(pdb_file):
    all_atom_coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                conf_b = line[16]
                if conf_b != ' ' and conf_b!='A':
                    continue
                
                atom_index = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                if atom_name in METAL_IONS:
                    continue
                
                # 1. 去除氢原子
                if atom_name.startswith('H'):
                    continue
                x, y, z = map(float, (line[30:38], line[38:46], line[46:54]))
                all_atom_coords.append([x, y, z])
    assert(len(all_atom_coords)!= 0)
    return all_atom_coords

def cal_clu(method, peptide_type, key, eps=10.0):
    centers = []

    for i in range(5):
        pdb_file = f'{tmp_file_path}/{method}/{peptide_type}/{key}/{i}/pred_pep.pdb'
        coords = get_central(pdb_file)
        if len(coords) == 0:
            raise ValueError(f"No valid atoms found in {pdb_file}")
        center = np.mean(coords, axis=0)
        centers.append(center)

    centers = np.array(centers)

    db = DBSCAN(eps=eps, min_samples=1).fit(centers)
    labels = db.labels_
    n_clusters = len(set(labels))

    return n_clusters


ori_dir = '/mnt/data/public/AF3_benchmark/all_scores_v1'
aim_dir = '/mnt/data/public/AF3_benchmark/all_scores_v2'

for method in method_list:
    for peptide_type in peptide_type_list:
        ori_score_file = f'{ori_dir}/{method}/{peptide_type}_all.csv'
        aim_score_file = f'{aim_dir}/{method}/{peptide_type}_all.csv'
        
        if not os.path.isfile(ori_score_file):
            continue
        if not os.path.isdir(f'{aim_dir}/{method}/'):
            os.makedirs(f'{aim_dir}/{method}/')
        wf = open(aim_score_file,'w')
        for line in open(ori_score_file,'r'):
            key = line.split(',')[0]
            line = line.replace('\n','')
            n_clusters = cal_clu(method, peptide_type, key)
            print(f'{method},{peptide_type},{key},{n_clusters}')
            wf.write(f'{line},{n_clusters}\n')
        wf.close()
from utils import cal_all_score
import os


peptide_type_list = [
    'SAA_seq',
    'USAA_ptm',
    'SAA_smiles',
    'USAA_smiles',
    'CYC_smiles',
    'CYC_ptm',
    'CYC_ptm_bond'
]
root_dir = '/mnt/data/public/AF3_benchmark/'
method = 'ProteniX'
score_file = 'scores.txt'
gt_dir = '/mnt/data/public/AF3_benchmark/gt'
temp_dir = '/mnt/data/public/AF3_benchmark/tmp'
gt_prot_chain = 'A'
gt_pep_chain = 'B'

checkpoint = ''
miss_list = []
no_repeat = []

bond = False
repeat = False
given_repeat_dir = ''

if __name__ == '__main__':
    
    for peptide_type in peptide_type_list:
        count = 0
        result_dir = f'{root_dir}/{method}/{peptide_type}'
        save_score_file = f'/mnt/data/public/AF3_benchmark/all_scores/{method}/{peptide_type}_all.csv'
        parent_folder = os.path.dirname(save_score_file)
        if not os.path.isdir(parent_folder):
            os.makedirs(parent_folder)
        
        wf = open(save_score_file,'w')
        for item in sorted(os.listdir(result_dir)):
            count += 1
            
            item_dir = f'{result_dir}/{item}'
            pdb_id = item[:4]
            if pdb_id < checkpoint:
                continue
            if pdb_id in miss_list:
                continue
            
            item_score_file = f'{item_dir}/{score_file}'
            ranking_scores = []
            for line in open(item_score_file,'r'):
                ranking_scores.append(float(line))
            assert(len(ranking_scores)==5)
            write_line = f'{pdb_id}'
            for i in range(5):
                print(f'Handling {method}/{peptide_type}/{pdb_id}/{i}   {count}/{len(os.listdir(result_dir))}')
                predicted_file = f'{item_dir}/rank_{i}.cif'
                gt_file = f'{gt_dir}/{pdb_id}.pdb'
                item_temp_dir = f'{temp_dir}/{method}/{peptide_type}/{pdb_id}/{i}'
                repeat_dir = ''
                if repeat and pdb_id not in no_repeat:
                    repeat_dir = given_repeat_dir
                protein_rmsd, pep_dist, pep_rmsd, pep_rmsd_alignd_prot = cal_all_score(predicted_file, gt_file,gt_prot_chain,gt_pep_chain,item_temp_dir,bond=bond,repeat=repeat_dir)
                write_line = f'{write_line},{ranking_scores[i]},{protein_rmsd},{pep_dist},{pep_rmsd},{pep_rmsd_alignd_prot}'
            wf.write(write_line+'\n')
        wf.close()
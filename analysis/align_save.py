from utils import  just_align
import os


peptide_type_list = [
    'SAA_seq',
    #'USAA_ptm',
    #'SAA_smiles',
    #'USAA_smiles',
    #'CYC_smiles',
    #'CYC_ptm',
    #'CYC_ptm_bond'
]
root_dir = '/mnt/data/public/AF3_benchmark/'
method = 'AlphaFold3'
score_file = 'score.txt'
gt_dir = '/mnt/data/public/AF3_benchmark/gt'
temp_dir = '/mnt/data/public/AF3_benchmark/tmp'
gt_prot_chain = 'A'
gt_pep_chain = 'B'

checkpoint = ''
bond = False

if __name__ == '__main__':
    for peptide_type in peptide_type_list:
        result_dir = f'{root_dir}/{method}/{peptide_type}'
        
        for item in sorted(os.listdir(result_dir)):
            
            item_dir = f'{result_dir}/{item}'
            pdb_id = item[:4]
            if pdb_id < checkpoint:
                continue
            
            item_score_file = f'{item_dir}/{score_file}'
            # each method will be different
            ranking_scores = []
            for line in open(item_score_file,'r'):
                ranking_scores.append(float(line))
            assert(len(ranking_scores)==5)
            write_line = f'{pdb_id}'
            for i in range(5):
                print(f'Handling {method}/{peptide_type}/{pdb_id}/{i}')
                predicted_file = f'{item_dir}/rank_{i}.cif'
                gt_file = f'{gt_dir}/{pdb_id}.pdb'
                item_temp_dir = f'{temp_dir}/{method}/{peptide_type}/{pdb_id}/{i}'
                just_align(predicted_file, gt_file,gt_prot_chain,gt_pep_chain,item_temp_dir,bond=bond,repeat='')
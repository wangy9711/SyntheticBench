# calculate center distance
from pymol import cmd  
from collections import defaultdict
import math
import numpy as np
import os
import networkx as nx


def modify_connect(file_a, file_b, key='CONECT'):
    # Read file A and remove lines starting with the key
    with open(file_a, 'r') as f:
        lines = f.readlines()
    new_lines = [line for line in lines if not line.startswith(key)]

    # Read file B and extract lines starting with the key
    with open(file_b, 'r') as f:
        key_lines = [line for line in f if line.startswith(key)]

    # Write modified content back to file A
    with open(file_a, 'w') as f:
        f.writelines(new_lines)
        f.writelines(key_lines)


def parse_pdb_to_graph(pdb_file, aim_chain='C'):
    atom_coords = {}
    graph = nx.Graph()
    node_set = set()
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21:22]
                if chain_id != aim_chain:
                    continue

                conf_b = line[16]
                if conf_b != ' ' and conf_b != 'A':
                    continue
                
                atom_index = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                
                # 1. Remove hydrogen atoms
                if atom_name.startswith('H'):
                    continue
                
                # 2. Use first valid character as atom type
                atom_type = atom_name[0]
                
                x, y, z = map(float, (line[30:38], line[38:46], line[46:54]))
                atom_coords[atom_index] = (x, y, z)
                
                # Force node ID as int and type as str to avoid hashing issues
                graph.add_node(int(atom_index), type=str(atom_type))
                node_set.add(int(atom_index))
            
            elif line.startswith("CONECT"):
                fields = list(map(int, line.split()[1:]))
                if fields[0] not in node_set:
                    continue
                for i in range(1, len(fields)):
                    if fields[i] not in node_set:
                        continue
                    if not graph.has_edge(fields[0], fields[i]):
                        graph.add_edge(int(fields[0]), int(fields[i]))
    
    return graph, atom_coords


def match_graphs(pdb_file_a, pdb_file_b, chain_a='A', chain_b='B'):
    graph_a, coords_a = parse_pdb_to_graph(pdb_file_a, chain_a)
    graph_b, coords_b = parse_pdb_to_graph(pdb_file_b, chain_b)
    
    # Define node matching rule (match by atom type)
    def node_match(n1, n2):
        return n1["type"] == n2["type"]
    
    # Use graph isomorphism to match
    if len(graph_a) < len(graph_b):
        graph_a, graph_b = graph_b, graph_a
        coords_a, coords_b = coords_b, coords_a

    matcher = nx.algorithms.isomorphism.GraphMatcher(graph_a, graph_b, node_match=node_match)
    
    matched_coords_a = []
    matched_coords_b = []
    
    for match in matcher.subgraph_isomorphisms_iter():
        matched_coords_a = [coords_a[node_a] for node_a in match.keys()]
        matched_coords_b = [coords_b[node_b] for node_b in match.values()]
        break
    if len(graph_a) < len(graph_b):
        return matched_coords_b, matched_coords_a
    return matched_coords_a, matched_coords_b


def kabsch_algorithm(P, Q):
    """
    Kabsch algorithm to compute the optimal rotation matrix and translation vector

    Parameters:
        P: Reference coordinates (n x 3 numpy array)
        Q: Coordinates to align (n x 3 numpy array)

    Returns:
        R: Rotation matrix (3 x 3)
        t: Translation vector (3,)
    """
    assert P.shape == Q.shape
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)

    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q
    H = P_centered.T @ Q_centered
    U, S, Vt = np.linalg.svd(H)
    R = U @ Vt

    if np.linalg.det(R) < 0:
        print("Warning: Chirality inversion detected!")
        Vt[-1, :] *= -1
        R = U @ Vt

    t = centroid_P - R @ centroid_Q
    return R, t


def apply_transformation(coords, R, t):
    """
    Apply rotation and translation to coordinates

    Parameters:
        coords: Coordinates (n x 3)
        R: Rotation matrix (3 x 3)
        t: Translation vector (3,)

    Returns:
        Transformed coordinates
    """
    return (R @ coords.T).T + t


def calculate_rmsd(P, Q):
    """
    Calculate RMSD between two sets of coordinates
    """
    diff = P - Q
    return np.sqrt(np.mean(np.sum(diff * diff, axis=1)))


def align_and_rmsd(coords_a, coords_b):
    R, t = kabsch_algorithm(coords_a, coords_b)
    aligned_coords_b = apply_transformation(coords_b, R, t)
    return calculate_rmsd(coords_a, aligned_coords_b)


def save_pdb_with_bonds(save_file_name, bonds, obj_name):
    with open(save_file_name, "w") as f:
        f.write(cmd.get_pdbstr(obj_name))
        # 追加 CONECT 记录
        for bond in bonds:
            atom1, atom2, _ = bond
            f.write(f"CONECT{atom1+1:5}{atom2+1:5}\n")

def just_align(prediction_file, gt_file, gt_prot_chain, gt_pep_chain, save_dir,bond=False,repeat=''):
    # load
    cmd.delete("all")
    cmd.load(prediction_file, f"pred")
    cmd.load(gt_file, "gt")
    
    
    def find_shorter_and_longer_chain(object_name):  
        """Identify the shorter chain in the given object."""  
        chain_lengths = defaultdict(int)  
        model = cmd.get_model(object_name)  
        for atom in model.atom:  
            chain_lengths[atom.chain] += 1  
        shorter_chain = min(chain_lengths, key=chain_lengths.get)  
        longer_chain = max(chain_lengths, key=chain_lengths.get)
        return longer_chain, shorter_chain
    
    # find prot and pep chain in prediction
    pred_prot_chain, pred_pep_chain = find_shorter_and_longer_chain('pred')
    if bond:
        pred_pep_chain = 'B'

    # align proteinm and get score1
    protein_rmsd = cmd.align(f"pred and chain {pred_prot_chain}", f"gt and chain {gt_prot_chain}")

    protein_rmsd = protein_rmsd[0]
    # get pep center distance(gt)
    gt_model = cmd.get_model(f"gt and chain {gt_pep_chain}")
    gt_pep_coords = np.array([[atom.coord[0], atom.coord[1], atom.coord[2]] for atom in gt_model.atom])
    gt_center = gt_pep_coords.mean(axis=0)

    # get pep center distance(pre)
    pred_model = cmd.get_model(f"pred and chain {pred_pep_chain}")        
    pred_pep_coords = np.array([[atom.coord[0], atom.coord[1], atom.coord[2]] for atom in pred_model.atom])
    pred_center = pred_pep_coords.mean(axis=0)

    # score2
    pep_dist = np.linalg.norm(gt_center - pred_center)

    # save file to align pep
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)

    gt_pep_file_name = f'{save_dir}/gt_pep.pdb'
    pred_pep_file_name = f'{save_dir}/pred_pep.pdb'
    cmd.set("connect_mode", 2)
    gt_bonds = cmd.get_bonds(f"gt and chain {gt_pep_chain}")
    pred_bonds = cmd.get_bonds(f"pred and chain {pred_pep_chain}")
    if not os.path.isfile(gt_pep_file_name):
        if repeat!='':
            ref_gt_file = f'{repeat}/0/gt_pep.pdb'
            os.system(f'cp {ref_gt_file} {gt_pep_file_name}')
            print('use repeat gt')
        elif save_dir[-1]!='0':
            ref_gt_file = f'{save_dir[:-1]}0/gt_pep.pdb'
            os.system(f'cp {ref_gt_file} {gt_pep_file_name}')
        else:
            save_pdb_with_bonds(gt_pep_file_name,gt_bonds,f"gt and chain {gt_pep_chain}")
    if not os.path.isfile(pred_pep_file_name):
        save_pdb_with_bonds(pred_pep_file_name,pred_bonds,f"pred and chain {pred_pep_chain}")
        if repeat!='':
            ref_connect_file = f'{repeat}/0/pred_pep.pdb'
            modify_connect(pred_pep_file_name,ref_connect_file)
            print('use repeat connect')
        elif save_dir[-1]!='0':
            ref_connect_file = f'{save_dir[:-1]}0/pred_pep.pdb'
            modify_connect(pred_pep_file_name,ref_connect_file)


def cal_all_score(prediction_file, gt_file, gt_prot_chain, gt_pep_chain, save_dir,bond=False,repeat=''):
    # load
    cmd.delete("all")
    cmd.load(prediction_file, f"pred")
    cmd.load(gt_file, "gt")
    
    
    def find_shorter_and_longer_chain(object_name):  
        """Identify the shorter chain in the given object."""  
        chain_lengths = defaultdict(int)  
        model = cmd.get_model(object_name)  
        for atom in model.atom:  
            chain_lengths[atom.chain] += 1  
        shorter_chain = min(chain_lengths, key=chain_lengths.get)  
        longer_chain = max(chain_lengths, key=chain_lengths.get)
        return longer_chain, shorter_chain
    
    # find prot and pep chain in prediction
    pred_prot_chain, pred_pep_chain = find_shorter_and_longer_chain('pred')
    if bond:
        pred_pep_chain = 'B'

    # align proteinm and get score1
    protein_rmsd = cmd.align(f"pred and chain {pred_prot_chain}", f"gt and chain {gt_prot_chain}")

    protein_rmsd = protein_rmsd[0]
    # get pep center distance(gt)
    gt_model = cmd.get_model(f"gt and chain {gt_pep_chain}")
    gt_pep_coords = np.array([[atom.coord[0], atom.coord[1], atom.coord[2]] for atom in gt_model.atom])
    gt_center = gt_pep_coords.mean(axis=0)

    # get pep center distance(pre)
    pred_model = cmd.get_model(f"pred and chain {pred_pep_chain}")        
    pred_pep_coords = np.array([[atom.coord[0], atom.coord[1], atom.coord[2]] for atom in pred_model.atom])
    pred_center = pred_pep_coords.mean(axis=0)

    # score2
    pep_dist = np.linalg.norm(gt_center - pred_center)

    # save file to align pep
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)

    gt_pep_file_name = f'{save_dir}/gt_pep.pdb'
    pred_pep_file_name = f'{save_dir}/pred_pep.pdb'
    cmd.set("connect_mode", 2)
    gt_bonds = cmd.get_bonds(f"gt and chain {gt_pep_chain}")
    pred_bonds = cmd.get_bonds(f"pred and chain {pred_pep_chain}")
    if not os.path.isfile(gt_pep_file_name):
        if repeat!='':
            ref_gt_file = f'{repeat}/0/gt_pep.pdb'
            os.system(f'cp {ref_gt_file} {gt_pep_file_name}')
            print('use repeat gt')
        elif save_dir[-1]!='0':
            ref_gt_file = f'{save_dir[:-1]}0/gt_pep.pdb'
            os.system(f'cp {ref_gt_file} {gt_pep_file_name}')
        else:
            save_pdb_with_bonds(gt_pep_file_name,gt_bonds,f"gt and chain {gt_pep_chain}")
    if not os.path.isfile(pred_pep_file_name):
        save_pdb_with_bonds(pred_pep_file_name,pred_bonds,f"pred and chain {pred_pep_chain}")
        if repeat!='':
            ref_connect_file = f'{repeat}/0/pred_pep.pdb'
            modify_connect(pred_pep_file_name,ref_connect_file)
            print('use repeat connect')
        elif save_dir[-1]!='0':
            ref_connect_file = f'{save_dir[:-1]}0/pred_pep.pdb'
            modify_connect(pred_pep_file_name,ref_connect_file)


    # read and build graph to align peptide
    # 
    if len(pred_pep_chain)==2:
        pred_pep_chain = pred_pep_chain[0]
    gt_mapped_coords, pred_mapped_coords = match_graphs(gt_pep_file_name, pred_pep_file_name, gt_pep_chain, pred_pep_chain)

    assert(np.array(gt_mapped_coords).shape == np.array(pred_mapped_coords).shape )

    gt_mapped_coords = np.array(gt_mapped_coords)
    pred_mapped_coords = np.array(pred_mapped_coords)

    # score3 
    pep_aligned_rmsd = align_and_rmsd(gt_mapped_coords, pred_mapped_coords)

    # score4
    pep_rmsd_prot_align = calculate_rmsd(gt_mapped_coords, pred_mapped_coords)

    return protein_rmsd, pep_dist, pep_aligned_rmsd, pep_rmsd_prot_align
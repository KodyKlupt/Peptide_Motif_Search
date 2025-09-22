
import argparse
import json
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch
from Bio.PDB.DSSP import DSSP
from Bio.PDB.vectors import calc_angle, calc_dihedral
from itertools import product
import os
import glob
import csv

dss_warning_printed = False #remove annoying warnings from dssp

def parse_motif_file(motif_path): ##grab the json file with all the definitions added
    with open(motif_path, 'r') as f:
        return json.load(f)

def get_atoms_from_structure(model, component): ##parse through the pdb and look for atoms from the motif definition
    matching_atoms = [] 
    for chain in model:
        for residue in chain:
            res_name = residue.get_resname().strip()
            if res_name == component['residue_type']:
                try:
                    atom_dict = {key: residue[atom_name] for key, atom_name in component['atom_selectors'].items()}
                    atom_dict['residue_info'] = {
                        'res_name': res_name,
                        'chain_id': chain.id,
                        'res_id': residue.id[1],
                        'biopython_residue_object': residue
                    }
                    matching_atoms.append(atom_dict) 
                except KeyError:
                    continue
    return matching_atoms

def check_primary_constraints(atom_combination, constraints, component_map, dssp_data):
    global dss_warning_printed
    for constraint in constraints:
        c_type = constraint['type']
        if c_type == 'distance': ##distance vector equation between two atoms
            id1, name1 = constraint['atoms'][0].split('.')
            id2, name2 = constraint['atoms'][1].split('.')
            dist = np.linalg.norm(atom_combination[component_map[id1]][name1].get_coord() - atom_combination[component_map[id2]][name2].get_coord())
            if not (constraint['value'] - constraint['tolerance'] <= dist <= constraint['value'] + constraint['tolerance']):
                return False
        elif c_type == 'angle': ##angle
            v1 = atom_combination[component_map[constraint['atoms'][0].split('.')[0]]][constraint['atoms'][0].split('.')[1]].get_vector()
            v2 = atom_combination[component_map[constraint['atoms'][1].split('.')[0]]][constraint['atoms'][1].split('.')[1]].get_vector()
            v3 = atom_combination[component_map[constraint['atoms'][2].split('.')[0]]][constraint['atoms'][2].split('.')[1]].get_vector()
            angle = np.rad2deg(calc_angle(v1, v2, v3))
            if not (constraint['value'] - constraint['tolerance'] <= angle <= constraint['value'] + constraint['tolerance']):
                return False
        elif c_type == 'dihedral': ##four atom vector dihedral angle
            v1 = atom_combination[component_map[constraint['atoms'][0].split('.')[0]]][constraint['atoms'][0].split('.')[1]].get_vector()
            v2 = atom_combination[component_map[constraint['atoms'][1].split('.')[0]]][constraint['atoms'][1].split('.')[1]].get_vector()
            v3 = atom_combination[component_map[constraint['atoms'][2].split('.')[0]]][constraint['atoms'][2].split('.')[1]].get_vector()
            v4 = atom_combination[component_map[constraint['atoms'][3].split('.')[0]]][constraint['atoms'][3].split('.')[1]].get_vector()
            dihedral = np.rad2deg(calc_dihedral(v1, v2, v3, v4))
            if not (constraint['value'] - constraint['tolerance'] <= dihedral <= constraint['value'] + constraint['tolerance']):
                return False
        elif c_type in ['secondary_structure', 'accessibility']:
            if not dssp_data:
                if not dss_warning_printed: ##sometimes cant find dssp so just skip it
                    print("WARNING: DSSP not found. Skipping secondary_structure and accessibility constraints.")
                    dss_warning_printed = True
                continue
            res_info = atom_combination[component_map[constraint['component_id']]]['residue_info']
            res_key = (res_info['chain_id'], (' ', res_info['res_id'], ' '))
            if res_key not in dssp_data:
                return False
            if c_type == 'secondary_structure' and dssp_data[res_key][2] not in constraint['value']: ##key 2 is the secondary structure assignment
                return False
            elif c_type == 'accessibility':
                rsa_value = dssp_data[res_key][3] ##this gives the relative solvent accessibility, from key 3 in the dssp
                if constraint['comparison'] == 'greater_than' and not rsa_value > constraint['value']:
                    return False
                if constraint['comparison'] == 'less_than' and not rsa_value < constraint['value']:
                    return False
    return True

def check_exclusion_constraints(atom_combination, constraints, component_map, ns):
    for constraint in constraints: 
        if constraint['type'] == 'exclusion_sphere': ##radius around an atom you dont want any other atoms to be in
            center_atom = atom_combination[component_map[constraint['component_id']]][constraint['atom_selector']]
            
            #there will probably be other atoms from the same residue so just ignore those
            parent_residue = center_atom.get_parent()
            atoms_to_ignore = {atom for atom in parent_residue.get_atoms()}

            neighbor_atoms = ns.search(center_atom.get_coord(), constraint['radius'], level='A')
            for neighbor in neighbor_atoms:
                if neighbor not in atoms_to_ignore: 
                    return False #found an atom in the exclusion sphere (not false is basically true)
    return True

def search_single_file(structure_path, motif_def):
    try:
        parser = MMCIFParser(QUIET=True) if structure_path.lower().endswith('.cif') else PDBParser(QUIET=True)
        model = parser.get_structure(os.path.basename(structure_path), structure_path)[0] ##convert from CIF to PDB for parsing
    except Exception as e:
        print(f"Error parsing {structure_path}: {e}")
        return []

    dssp_data = None
    try:
        dssp_data = DSSP(model, structure_path, dssp='dssp') 
    except Exception: pass

    all_atoms = [a for a in model.get_atoms() if a.get_parent().id[0] == ' ']
    ns = NeighborSearch(all_atoms)

    primary_constraints = [c for c in motif_def['constraints'] if c['type'] != 'exclusion_sphere']
    exclusion_constraints = [c for c in motif_def['constraints'] if c['type'] == 'exclusion_sphere']

    component_map = {comp['id']: i for i, comp in enumerate(motif_def['components'])} 
    candidate_atoms_per_component = [get_atoms_from_structure(model, c) for c in motif_def['components']] 
    if not all(candidate_atoms_per_component):
        return []

    found_motifs = []
    for atom_combination in product(*candidate_atoms_per_component): 
        if len(set(res['residue_info']['res_id'] for res in atom_combination)) < len(motif_def['components']):
            continue 

        if check_primary_constraints(atom_combination, primary_constraints, component_map, dssp_data):
            if check_exclusion_constraints(atom_combination, exclusion_constraints, component_map, ns):
                pruned_residue_info = []
                for res in atom_combination:
                    info = res['residue_info']
                    pruned_residue_info.append({k: v for k, v in info.items() if k != 'biopython_residue_object'})
                
                match = {"residues": sorted(pruned_residue_info, key=lambda x: x['res_id'])} 
                if match not in found_motifs:
                    found_motifs.append(match)
    return found_motifs

def main():
    parser = argparse.ArgumentParser(description="Batch search for 3D structural motifs.")
    parser.add_argument("-i", "--input_folder", required=True, help="Folder with PDB/CIF files.")
    parser.add_argument("-m", "--motif_file", required=True, help="JSON motif definition file.")
    parser.add_argument("-o", "--output_folder", required=True, help="Folder to save JSON results.")
    parser.add_argument("-s", "--summary_csv", required=True, help="Final summary CSV file.")
    args = parser.parse_args()

    motif_def = parse_motif_file(args.motif_file) 
    os.makedirs(args.output_folder, exist_ok=True)
    print(f"Loaded motif '{motif_def['motif_name']}'. Searching in folder: {args.input_folder}")

    all_files = sorted(glob.glob(os.path.join(args.input_folder, '*.pdb')) + glob.glob(os.path.join(args.input_folder, '*.cif')))
    if not all_files:
        print("No .pdb or .cif files found.")
        return

    csv_header = ['source_file', 'motif_id'] + [f'residue_{i+1}' for i in range(len(motif_def['components']))]
    csv_summary_data = [] 

    for structure_file in all_files:
        base_name = os.path.basename(structure_file) 
        print(f"-- Processing {base_name}...")
        found_motifs = search_single_file(structure_file, motif_def) 

        with open(os.path.join(args.output_folder, base_name.rsplit('.', 1)[0] + '.json'), 'w') as f:
            json.dump({"source_file": base_name, "motifs_found": len(found_motifs), "matches": found_motifs}, f, indent=2)

        if found_motifs:
            print(f"  Found {len(found_motifs)} motif(s).")
            for i, motif in enumerate(found_motifs):
                row = {'source_file': base_name, 'motif_id': f"{base_name.rsplit('.', 1)[0]}_motif_{i+1}"}
                for j, res in enumerate(motif['residues']):
                    row[f'residue_{j+1}'] = f"{res['res_name']}-{res['chain_id']}-{res['res_id']}"
                csv_summary_data.append(row)
        else:
            row = {key: '' for key in csv_header}
            row['source_file'] = base_name
            row['motif_id'] = 'NO_MOTIF_FOUND'
            csv_summary_data.append(row)

    print(f"Writing summary for {len(all_files)} files to {args.summary_csv}")
    with open(args.summary_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=csv_header)
        writer.writeheader()
        writer.writerows(csv_summary_data)
    
    print("Processing complete.")

if __name__ == '__main__':
    main()

import collections
import json
import os
import csv
import pandas
import time
import re


AMINO_ACID_NOMENCLATURE = {
    'A': ['A'], 'C': ['C'], 'D': ['D'], 'E': ['E'], 'F': ['F'], 'G': ['G'], 'H': ['H'], 'I': ['I'], 'K': ['K'], 'L': ['L'], 'M': ['M'], 'N': ['N'], 'P': ['P'], 'Q': ['Q'], 'R': ['R'], 'S': ['S'], 'T': ['T'], 'V': ['V'], 'W': ['W'], 'Y': ['Y'],
    '%': ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'], #hydrophobic
    '@': ['F', 'Y', 'W', 'H'], #aromatic....if you include histidine, its not super strict here, delete if needed
    '&': ['R', 'N', 'D', 'Q', 'E', 'K', 'H', 'S', 'T', 'Y'], #polar
    '[+]': ['K', 'R', 'H'], #positively charged
    '[-]': ['D', 'E'], #negatively charged
    '#': ['A', 'V', 'L', 'I'], #aliphatic
    '~': ['A', 'C', 'D', 'G', 'N', 'P', 'S', 'T', 'V'], #small
    'x': ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'], #any amino acid
}
PTM_NOMENCLATURE = {
    '[Y:po]': ['Y'], '[Y:su]': ['Y'], '[S:gl]': ['S'], '[N:gl]': ['N'], '[R:me]': ['R'], '[R:sme]': ['R'], '[R:ame]': ['R'], '[K:ac]': ['K'], '[P:hy]': ['P']
}
AMINO_ACID_NOMENCLATURE.update(PTM_NOMENCLATURE)

def expand_motif(motif):
    """
    This will write the motifs into concrete sequences and then you can search through the protein sequence 
    """
    if not motif:
        return [""]

    token = None
    token_end_index = 0

    if motif.startswith('['):
        token_end_index = motif.find(']') + 1
        token = motif[:token_end_index]
    elif motif.startswith('{'):
        token_end_index = motif.find('}') + 1
        token = motif[:token_end_index]
    else:
        token_end_index = 1
        token = motif[0]

    options = []
    if token in AMINO_ACID_NOMENCLATURE:
        options = AMINO_ACID_NOMENCLATURE[token]
    elif token.startswith('[') and token.endswith(']'):
        options = list(token[1:-1])
    elif token.startswith('{') and token.endswith('}'):
        excluded_chars = set(list(token[1:-1]))
        all_amino_acids = AMINO_ACID_NOMENCLATURE['x']
        options = [aa for aa in all_amino_acids if aa not in excluded_chars]
    else:
        options = [token]

    remaining_motif_expansions = expand_motif(motif[token_end_index:])
    return [aa + r for aa in options for r in remaining_motif_expansions]


class AhoCorasick:
    """
    aho corasick for motif searching
    """
    def __init__(self, keywords):
        self._node_count = 1
        self._trie = collections.defaultdict(dict)
        self._outputs = collections.defaultdict(list)
        self._failure = {}
        self._build_trie(keywords)
        self._build_failure_links()

    def _build_trie(self, keywords):
        for keyword in keywords:
            node = 0
            for char in keyword:
                if char not in self._trie[node]:
                    self._trie[node][char] = self._node_count
                    self._node_count += 1
                node = self._trie[node][char]
            self._outputs[node].append(keyword)
    
    def _build_failure_links(self):
        queue = collections.deque()
        for char, next_node in self._trie[0].items():
            queue.append(next_node)
            self._failure[next_node] = 0
            
        while queue:
            current_node = queue.popleft()
            for char, next_node in self._trie[current_node].items():
                queue.append(next_node)
                fail_node = self._failure[current_node]
                
                while char not in self._trie[fail_node] and fail_node != 0:
                    fail_node = self._failure[fail_node]
                
                if char in self._trie[fail_node]:
                    self._failure[next_node] = self._trie[fail_node][char]
                else:
                    self._failure[next_node] = 0
                
                fail_output_node = self._failure[next_node]
                self._outputs[next_node].extend(self._outputs[fail_output_node])

    def search(self, text):
        results = collections.defaultdict(list)
        node = 0
        for i, char in enumerate(text):
            while char not in self._trie[node] and node != 0:
                node = self._failure[node]
            
            if char in self._trie[node]:
                node = self._trie[node][char]
            
            if self._outputs[node]:
                for keyword in self._outputs[node]:
                    results[keyword].append(i)
        return results

#run motif search

def sanitize_filename(name):
    """Sanitizes a string to be a valid filename."""
    #clean the filename
    s = re.sub(r'[<>:"/\\|?*()]', '', name)
    s = s.replace(' ', '_')
    # Truncate to avoid "File name too long" errors
    return s[:100]

def run_motif_search(motifs_file, motif_column, sequences_file, name_column, sequence_column, output_file):
    """
    motif searching
    """
    motifs_to_find = []
    if os.path.exists(motifs_file):
        with open(motifs_file, 'r', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            motifs_to_find = [row[motif_column] for row in reader]
    else:
        print(f"Warning: Motif file '{motifs_file}' not found. No motifs will be searched.")

    protein_sequences = {}
    if os.path.exists(sequences_file):
        with open(sequences_file, 'r', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            protein_sequences = {row[name_column]: row[sequence_column] for row in reader}
    else:
        print(f"Error: Sequences file '{sequences_file}' not found.")
        return

    if not motifs_to_find:
        print("No motifs to find. Please populate your motifs file.")
        return

    # timestamped json files for each sequence
    output_dir_jsons = "Outputs/" + time.strftime("%Y%m%d_%H%M%S") + "_motif_search_results_jsons"
    os.makedirs(output_dir_jsons, exist_ok=True)
    print(f"Individual JSON results will be saved in: '{output_dir_jsons}/'")
    
    df_sequences_motifs = pandas.DataFrame(columns=['name', 'sequence', 'motifs'])

    for name, sequence in protein_sequences.items():
        print(f"\nProcessing sequence '{name}' (length: {len(sequence)})")
        
        expansion_map = {}
        all_expanded_motifs = []
        for m in motifs_to_find:
            processed_motif = m.replace('-', '')
            expanded_list = expand_motif(processed_motif)
            all_expanded_motifs.extend(expanded_list)
            for concrete_motif in expanded_list:
                expansion_map[concrete_motif] = m
        
        unique_expanded_motifs = list(set(all_expanded_motifs))
        automaton = AhoCorasick(unique_expanded_motifs)
        found_concrete_motifs = automaton.search(sequence)
        
        final_results = collections.defaultdict(list)
        for concrete_motif, positions in found_concrete_motifs.items():
            original_motif = expansion_map[concrete_motif]
            for end_pos in positions:
                final_results[original_motif].append((end_pos, concrete_motif))
        ##replace the name with a sanitized version, cause UniProt sometimes has very long names causing issues
        name = sanitize_filename(name)
        
        # saving the json
        json_output_data = {
            "name": name,
            "sequence": sequence,
            "results": dict(final_results)
        }
        json_output_path = os.path.join(output_dir_jsons, f"motif_search_results_{name}.json")
        with open(json_output_path, 'w', encoding='utf-8') as f:
            json.dump(json_output_data, f, indent=4)
        
        # adding sequences to dataframe with motifs found, writing a column from the json
        new_row_data = {
            'name': name,
            'sequence': sequence,
            'motifs': json.dumps(dict(final_results))
        }
        df_sequences_motifs = pandas.concat([df_sequences_motifs, pandas.DataFrame([new_row_data])], ignore_index=True)
        
        if not final_results:
            print(f"  -> No motifs found in '{name}'.")
        else:
            print(f"  -> Found {len(final_results)} unique motifs in '{name}'.")
        print(f"  -> Detailed results saved to {json_output_path}")

    # export results to one csv
    df_sequences_motifs.to_csv(output_dir_jsons + "_all_results.csv", index=False)
    output_file = output_dir_jsons + "_all_results.csv"
    print(f"\nAll results aggregated and saved to {output_file}")
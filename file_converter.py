import os
import csv

# Biopython for PDB parsing
from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.SeqUtils import seq1

# Biotite for CIF and FASTA processing
import biotite.structure.io.pdbx as pdbx
import biotite.sequence.io.fasta as fasta
from biotite.sequence import NucleotideSequence, AlphabetError

def translate_dna_orf_search(dna_sequence):
    protein_sequences, _ = dna_sequence.translate(complete=False)
    return "".join([str(p) for p in protein_sequences])

def translate_dna_forced(dna_sequence):
    """
    Translates the entire DNA sequence, ignoring start/stop codons.
    This assumes the input is a clean coding sequence.
    """
    protein_sequence = dna_sequence.translate(complete=True)
    return str(protein_sequence)

def process_protein_files(folder_path, output_csv):
    """
    Processes a folder of PDB, CIF, or FASTA files, converts all chains and
    translated DNA to one-letter protein sequences, concatenates them,
    and saves to a CSV file.

    Args:
        folder_path (str): The path to the folder containing the protein files.
        output_csv (str): The path to the output CSV file.
    """
    pdb_parser = PDBParser(QUIET=True)

    with open(output_csv, 'w', newline='', encoding='utf-8') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['name', 'sequence'])

        for filename in os.listdir(folder_path):
            input_path = os.path.join(folder_path, filename)
            concatenated_sequence = ""

            try:
                if filename.endswith(".pdb"):
                    structure = pdb_parser.get_structure(filename, input_path)
                    for model in structure:
                        for chain in model:
                            # Use a list comprehension to get the sequence of 3-letter codes
                            residues = [res.get_resname() for res in chain if res.id[0] == ' ']
                            # Convert the 3-letter code sequence to a 1-letter sequence
                            seq = seq1("".join(residues), custom_map={"MSE": "M"})
                            concatenated_sequence += seq

                elif filename.endswith(".cif"):
                    cif_file = pdbx.CIFFile.read(input_path)
                    sequences = pdbx.get_sequence(cif_file)
                    for seq in sequences.values():
                         concatenated_sequence += str(seq)

                elif filename.endswith(".fasta") or filename.endswith(".fa"):
                    fasta_file = fasta.FastaFile.read(input_path)
                    for header, sequence_str in fasta_file.items():
                        try:
                            # Attempt to treat as DNA first
                            dna_seq = NucleotideSequence(sequence_str)
                            concatenated_sequence += translate_dna_forced(dna_seq)

                        except AlphabetError:
                            # If not DNA as the input then will concatenate protein sequences
                            concatenated_sequence += sequence_str

                if concatenated_sequence:
                    csv_writer.writerow([filename, concatenated_sequence])
                    print(f"Processed {filename}")

            except Exception as e:
                print(f"Could not process {filename}: {e}")
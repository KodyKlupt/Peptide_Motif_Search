import os
import pandas as pd
import csv

##for every motif e.g. ADADA, write it as A-D-A-D-A
def convert_motif_to_regex(motif):
    motif = motif.replace(' ', '')
    motif = motif.replace('-', '')
    motif = motif.replace(',', '')
    motif = motif.replace(';', '')
    regex_motif = '-'.join(list(motif))
    return regex_motif
def process_kinase_motif_file(input_file, output_file):
    """
    Convert kinase motif data to a standardized format.
    """
    df = pd.read_csv(input_file)
    # Assuming the input file has columns 'Kinase' and 'Motif'
    df['motifs'] = df['Centralized AA Sequence'].apply(convert_motif_to_regex)
    df.to_csv(output_file, index=False)
    print(f"Converted motifs saved to {output_file}")
    
if __name__ == "__main__":
    input_file = "motif_libraries/Sugiyama_Kinase_Motifs.csv"  # Input file path
    output_file = "motif_libraries/Sugiyama_converted_kinase_motifs.csv"  # Output file path
    process_kinase_motif_file(input_file, output_file)
    
    ##clean the output file to only have kinase and motifs columns
    df = pd.read_csv(output_file)
    df = df[['Kinase', 'motifs']]
    df.to_csv(output_file, index=False)
    print(f"Cleaned motifs saved to {output_file}")
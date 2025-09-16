import requests
import pandas as pd
import sys
from io import StringIO

API_URL = "https://rest.uniprot.org/uniprotkb/search"

def search_uniprot(query, limit=500):
    """
    Searches UniProt with a given query for sequences
    """
    params = {
        "query": query,
        "format": "tsv",
        "fields": "accession,protein_name,sequence,organism_name,ec",
        "size": limit
    }
    try:
        response = requests.get(API_URL, params=params)
        response.raise_for_status()
        return response.text
    except requests.exceptions.RequestException as e:
        print(f"Error during UniProt API request: {e}", file=sys.stderr)
        sys.exit(1)

def to_csv(data, output_file):
    """
    Converts the UniProt data in TSV format to a CSV file.
    """
    if not data:
        print("No data received from UniProt.", file=sys.stderr)
        return

    lines = data.strip().split('\n')
    if len(lines) <= 1:
        print("No results found for the given query.", file=sys.stderr)
        return

    #create dataframe and store the results from uniprot
    df = pd.read_csv(StringIO(data), sep='\t')
    df.rename(columns={'Protein name': 'protein_name', 'Sequence': 'sequence', 'Organism': 'organism_name', 'EC number': 'ec'}, inplace=True)
    
    df.rename(columns={'Protein name': 'name'}, inplace=True)
    df.to_csv(output_file, index=False)
    print(f"Data saved to {output_file}")

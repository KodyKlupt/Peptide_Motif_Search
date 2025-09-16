import argparse
from argparse import _SubParsersAction
import os
import sys
from file_converter import process_protein_files
from motif_searcher import run_motif_search
from uniprot_api import search_uniprot, to_csv

def main():
    """
    Main function to parse command-line arguments and orchestrate the
    conversion of protein files and/or the search for motifs.
    """
    parser = argparse.ArgumentParser(
        description='A tool for protein motif searching. Can either convert protein files '
                    'to a sequence CSV or search for motifs in a sequence CSV.'
    )
    
    subparsers = parser.add_subparsers(dest='command', required=True, help='Available commands')

    # --- Subparser for the 'convert' command ---
    parser_convert = subparsers.add_parser('convert', help='Convert protein files (.pdb, .cif, .fasta) into a CSV file.')
    parser_convert.add_argument(
        '--input_folder', 
        type=str, 
        default='protein_files',
        help='Path to the folder containing protein files. Default: protein_files'
    )
    parser_convert.add_argument(
        '--output_csv', 
        type=str, 
        default='sequences.csv',
        help='Path to the output CSV file. Default: sequences.csv'
    )

    # --- Subparser for the 'search' command ---
    parser_search = subparsers.add_parser('search', help='Search for motifs in a CSV file of protein sequences.')
    parser_search.add_argument(
        '--motifs', 
        type=str, 
        default='motifs.csv', 
        help='Path to the CSV file containing motifs. Default: motifs.csv'
    )
    parser_search.add_argument(
        '-mc', '--motif_column', 
        dest='motif_column', 
        type=str, 
        default='motifs',
        help='Name of the column containing motifs in the motif file.'
    )
    parser_search.add_argument(
        '--sequences', 
        type=str, 
        default='sequences.csv', 
        help='Path to the CSV file containing protein sequences. Default: sequences.csv'
    )
    parser_search.add_argument(
        '-sc', '--sequence_column', 
        dest='sequence_column', 
        type=str, 
        default='sequence',
        help='Name of the column containing sequences in the sequences file.'
    )
    parser_search.add_argument(
        '--output', 
        type=str, 
        default='motif_search_results.csv', 
        help='Path to the output CSV file for results. Default: motif_search_results.csv'
    )
    parser_search.add_argument(
        '--name_column', 
        dest='name_column',
        type=str, 
        default='name',
        help='Name of the column containing sequence names/IDs')

    # --- Subparser for the 'uniprot' command ---
    parser_uniprot = subparsers.add_parser('uniprot', help='Fetch protein data from UniProt.')
    parser_uniprot.add_argument("--query", type=str, help="UniProt query string.")
    parser_uniprot.add_argument("--organism", type=str, help="Filter by organism.")
    parser_uniprot.add_argument("--enzyme_family", type=str, help="Filter by enzyme family.")
    parser_uniprot.add_argument("--protein_name", type=str, help="Filter by protein name.")
    parser_uniprot.add_argument("--accession", type=str, help="Filter by accession number.")
    parser_uniprot.add_argument("--output_csv", type=str, default="sequences.csv", help="Output CSV file.")
    parser_uniprot.add_argument("--limit", type=int, default=500, help="Limit number of results.")

    args = parser.parse_args()

    # runs either the convert or the search commands here
    if args.command == 'convert':
        if not os.path.isdir(args.input_folder):
            print(f"Error: The input folder '{args.input_folder}' does not exist.", file=sys.stderr)
            sys.exit(1)
        
        print(f"Starting file conversion from '{args.input_folder}'...")
        process_protein_files(args.input_folder, args.output_csv)
        print(f"\nConversion complete. Sequences saved to '{args.output_csv}'.")

    elif args.command == 'search':
        if not os.path.exists(args.sequences):
            print(f"Error: The sequences file '{args.sequences}' was not found.", file=sys.stderr)
            print("Please generate it first using the 'convert' command or provide a valid file path.", file=sys.stderr)
            sys.exit(1)
            
        print("Starting motif search...")
        run_motif_search(
            motifs_file=args.motifs,
            motif_column=args.motif_column,
            sequences_file=args.sequences,
            sequence_column=args.sequence_column,
            output_file=args.output,
            name_column=args.name_column
        )
        print("\nMotif search completed.")

    elif args.command == 'uniprot':
        query_parts = []
        if args.query:
            query_parts.append(f"({args.query})")
        if args.organism:
            query_parts.append(f"(organism_name:\"{args.organism}\")")
        if args.enzyme_family:
            query_parts.append(f"(family:\"{args.enzyme_family}\")")
        if args.protein_name:
            query_parts.append(f"(protein_name:\"{args.protein_name}\")")
        if args.accession:
            query_parts.append(f"(accession:{args.accession})")

        if not query_parts:
            parser.error("At least one query parameter must be specified for the 'uniprot' command.")

        query = " AND ".join(query_parts)
        print(f"Executing UniProt query: {query}")

        data = search_uniprot(query, args.limit)
        to_csv(data, args.output_csv)

if __name__ == '__main__':
    main()

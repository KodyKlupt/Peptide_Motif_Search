# Protein Motif Search Tool Documentation

This document provides a guide on how to use the different features of the Protein Motif Search Tool, including file conversion, motif searching, and fetching data from UniProt.

## Project Structure

For the tool to function correctly, it is recommended to have the following directory structure:

```
.
├── main.py                 # Main script to run the tool
├── file_converter.py       # Logic for file conversion
├── motif_searcher.py       # Logic for motif searching
├── uniprot_api.py          # Logic for UniProt API integration
|
├── protein_files/          # FOLDER: Place your input .pdb, .cif, .fasta files here
|
├── motif_libraries/        # FOLDER: Store your motif definitions here
|
└── sequences.csv           # OUTPUT: Default output for the 'convert' command
```

## How to Use

The tool is controlled through `main.py` and has three primary commands: `convert`, `search`, and `uniprot`.

### 1. `convert`: Convert Protein Files

If your protein sequences are in `.pdb`, `.cif`, or `.fasta` files, you first need to convert them into a single `sequences.csv` file.

**Usage:**
```bash
python main.py convert [OPTIONS]
```

**Arguments:**

| Argument          | Default Value     | Description                                             |
| ----------------- | ----------------- | ------------------------------------------------------- |
| `--input_folder`  | `protein_files`   | Path to the folder containing your protein source files.  |
| `--output_csv`    | `sequences.csv`   | Path for the generated CSV file of sequences.           |

**Example:**

This command reads all files from the `protein_files` directory and creates `sequences.csv` with the extracted protein sequences.

```bash
python main.py convert --input_folder protein_files --output_csv sequences.csv
```

### 2. `search`: Search for Motifs

Once you have a `sequences.csv` file, you can search for motifs. You will also need a separate CSV file containing the motifs you want to find.

**Usage:**
```bash
python main.py search [OPTIONS]
```

**Arguments:**

| Argument            | Flag | Required | Default Value              | Description                                                          |
| ------------------- | ---- | -------- | -------------------------- | -------------------------------------------------------------------- |
| `--motifs`          |      | **Yes**  | `motifs.csv`               | Path to the CSV file containing the motifs to search for.            |
| `--motif_column`    | `-mc`| **Yes**  |                            | Name of the column in the motif file that contains the motifs.         |
| `--sequences`       |      | **Yes**  | `sequences.csv`            | Path to the CSV file containing the protein sequences.               |
| `--sequence_column` | `-sc`| **Yes**  |                            | Name of the column in the sequences file that contains the sequences.  |
| `--output`          |      | No       | `motif_search_results.csv` | Path for the output CSV file that will store the aggregate results.  |
| `--name_column`     |      | No       | `name`                     | Name of the column containing sequence names/IDs.                    |

**Example:**

This command searches for motifs defined in `examples/motif_libraries/glycosylation_motifs.csv` within the sequences from `examples/sequences/example_sequences.csv` and saves the results to `final_results.csv`.

```bash
python main.py search \
    --motifs examples/motif_libraries/glycosylation_motifs.csv \
    -mc "motifs" \
    --sequences examples/sequences/example_sequences.csv \
    -sc "sequence" \
    --output "final_results.csv" \
    --name_column "name"
```

### 3. `uniprot`: Fetch Sequences from UniProt

You can use the UniProt API to retrieve sequences based on various criteria.

**Usage:**
```bash
python main.py uniprot [OPTIONS]
```

**Arguments:**

| Argument | Required | Default Value | Description |
| --- | --- | --- | --- |
| `--query` | No | | A general UniProt query string. |
| `--organism` | No | | Filter by organism name. |
| `--enzyme_family` | No | | Filter by enzyme family. |
| `--protein_name` | No | | Filter by protein name. |
| `--accession` | No | | Filter by accession number. |
| `--output_csv` | No | `sequences.csv` | The name of the output file. |
| `--limit` | No | 500 | The maximum number of sequences to retrieve. |

**Example:**

This command fetches up to 10 sequences for "ice-binding protein" from the UniProt database and saves them to `uniprot_sequences.csv`.

```bash
python main.py uniprot \
    --protein_name "ice-binding protein" \
    --output_csv "uniprot_sequences.csv" \
    --limit 10
```

## Motif Nomenclature

The tool supports a flexible motif syntax. Here are the special characters you can use in your motif definitions:

| Symbol | Example      | Description                                                               |
| :----: | ------------ | ------------------------------------------------------------------------- |
| `x`    | `AxC`        | Any of the 20 standard amino acids.                                       |
| `[ ]`  | `[AGV]xP`    | **Custom Set**: Matches any single character inside the brackets (A, G, or V). |
| `{ }`  | `{PG}x[DE]`  | **Exclusion**: Matches any amino acid *except* those in the curly braces.   |
| `%`    | `%xL`        | **Hydrophobic**: `A, V, I, L, M, F, Y, W`                                 |
| `@`    | `@S`         | **Aromatic**: `F, Y, W, H`                                                |
| `&`    | `&T`         | **Polar**: `R, N, D, Q, E, K, H, S, T, Y`                                 |
| `[+]`  | `[+]G`       | **Positively Charged**: `K, R, H`                                         |
| `[-]`  | `[-]G`       | **Negatively Charged**: `D, E`                                            |
| `#`    | `#L`         | **Aliphatic**: `A, V, L, I`                                               |
| `~`    | `~P`         | **Small**: `A, C, D, G, N, P, S, T, V`                                    |
```
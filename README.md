# Protein Motif Search Tool

This command-line tool provides a powerful and flexible way to search for protein motifs within sequences. It uses the efficient Aho-Corasick algorithm for fast pattern matching and supports a rich nomenclature for defining complex motifs.

The tool is organized into two main functions:
1.  **Convert**: Processes various biological file formats (`.pdb`, `.cif`, `.fasta`) into a standardized CSV format containing protein sequences.
2.  **Search**: Scans the sequence CSV file to find all occurrences of specified motifs.

## Features

*   **Efficient Searching**: Utilizes the Aho-Corasick algorithm to find multiple motifs in a single pass.
*   **Flexible Motif Definition**: Supports a rich nomenclature including wildcards, custom character sets, exclusions, and biochemical property groups.
*   **Multiple Input Formats**: Converts `.pdb`, `.cif`, and `.fasta` files (including DNA sequences, which are translated into proteins) into a uniform format.
*   **Modular and Extensible**: The logic is separated into distinct scripts for file conversion and motif searching, orchestrated by a central command-line interface.
*   **Detailed Outputs**: Generates an aggregate CSV file with all results, plus individual JSON files for each sequence with detailed match information.

## Project Structure

For the tool to function correctly, please organize your files in the following directory structure:

```
.
├── main.py                 # Main script to run the tool
├── file_converter.py       # Logic for file conversion
├── motif_searcher.py       # Logic for motif searching
|
├── protein_files/          # FOLDER: Place your input .pdb, .cif, .fasta files here
│   ├── protein1.pdb
│   └── sequence2.fasta
|
├── motif_libraries/        # FOLDER: Store your motif definitions here
│   └── my_motifs.csv
|
└── sequences.csv           # OUTPUT: Default output for the 'convert' command
```

## Installation

This project requires a few external Python libraries. You can install them using pip:

```bash
pip install pandas biopython biotite
```
or:
```bash
pip install -r requirements.txt
```


## How to Use

The tool is controlled through `main.py` and has two primary commands: `convert` and `search`. You must choose one of them.

### Step 1: `convert` (Optional)

If your protein sequences are in `.pdb`, `.cif`, or `.fasta` files, you first need to convert them into a single `sequences.csv` file. If you already have a CSV file with protein sequences, you can skip this step.

**Usage:**
```bash
python main.py convert [OPTIONS]
```

**Arguments for `convert`:**

| Argument          | Default Value     | Description                                             |
| ----------------- | ----------------- | ------------------------------------------------------- |
| `--input_folder`  | `protein_files`   | Path to the folder containing your protein source files.  |
| `--output_csv`    | `sequences.csv`   | Path for the generated CSV file of sequences.           |

**Example:**
```bash
# This command reads all files from the 'protein_files' directory
# and creates 'sequences.csv' with the extracted protein sequences.
python main.py convert --input_folder protein_files --output_csv sequences.csv
```

### Step 2: `search`

Once you have a `sequences.csv` file, you can search for motifs. You will also need a separate CSV file containing the motifs you want to find.

**Usage:**
```bash
python main.py search [OPTIONS]
```

**Arguments for `search`:**

| Argument            | Flag | Required | Default Value              | Description                                                          |
| ------------------- | ---- | -------- | -------------------------- | -------------------------------------------------------------------- |
| `--motifs`          |      | **Yes**  | `motifs.csv`               | Path to the CSV file containing the motifs to search for.            |
| `--motif_column`    | `-mc`| **Yes**  |                            | Name of the column in the motif file that contains the motifs.         |
| `--sequences`       |      | **Yes**  | `sequences.csv`            | Path to the CSV file containing the protein sequences.               |
| `--sequence_column` | `-sc`| **Yes**  |                            | Name of the column in the sequences file that contains the sequences.  |
| `--output`          |      | No       | `motif_search_results.csv` | Path for the output CSV file that will store the aggregate results.  |

**Example:**
```bash
# This command searches for motifs defined in 'my_motifs.csv' within the
# sequences from 'sequences.csv' and saves the results to 'final_results.csv'.
python main.py search \
    --motifs motif_libraries/my_motifs.csv \
    -mc "motif_pattern" \
    --sequences sequences.csv \
    -sc "sequence" \
    --output "final_results.csv"
```

## Input File Formats

#### Sequences File (`--sequences`)

The sequences file must be a CSV with at least two columns: one for a unique identifier and one for the protein sequence.

*Example (`sequences.csv`):*
```csv
name,sequence
ProteinA,MDSGSEYGPLVHEFKNADLSLDKFN...
ProteinB,MAAVVGGASFGGHJKLMNPQRS...
```

#### Motifs File (`--motifs`)

The motifs file is a CSV that lists the patterns to search for.

*Example (`motif_libraries/my_motifs.csv`):*
```csv
motif_id,motif_pattern,description
M001,R[ST]xP,"A common kinase motif"
M002,{P}G,"Proline-Glysine exclusion"
M003,#x[+],"Aliphatic followed by any aa and a positive charge"```

## Motif Nomenclature

The power of this tool comes from its flexible motif syntax. The following special characters can be used in the `--motif_column`:

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

*Note: Post-translational modifications (e.g., `[Y:po]`) are also supported and map to their base amino acid.*

## Output Format

The `search` command produces two types of output:

1.  **Aggregate CSV File (`--output`)**: A single CSV file summarizing all findings.
    *   `name`: The name of the sequence from the input file.
    *   `sequence`: The full protein sequence.
    *   `motifs`: A JSON string containing a dictionary where keys are the original motifs and values are a list of all matches found. Each match includes the end position and the specific concrete motif that matched.

2.  **Individual JSON Files**: For each sequence processed, a detailed JSON file is saved in a timestamped directory (e.g., `20250725_163000_motif_search_results_jsons/`).
    *   `name`: The name of the sequence.
    *   `sequence`: The full protein sequence.
    *   `results`: A dictionary where keys are the original motifs found. The values are lists of `[end_position, "concrete_motif"]` pairs.

*Example JSON (`motif_search_results_ProteinA.json`):*
```json
{
    "name": "ProteinA",
    "sequence": "MDSGSEYGPLVHEFKRSTP...",
    "results": {
        "R[ST]xP": [
            [
                18,
                "RSTP"
            ]
        ]
    }
}
```
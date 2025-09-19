# Protein Motif Search Tools

This command-line tool provides a powerful and flexible way to search for protein motifs within sequences and structures. It is organized into two main tools:

1.  **Sequence Motif Search**: Searches for motifs in protein sequences using a rich nomenclature for defining complex patterns.
2.  **Structure Motif Search**: Searches for 3D structural motifs in protein structures using a flexible JSON-based definition format.

---

## Sequence Motif Search

This tool uses the efficient Aho-Corasick algorithm for fast pattern matching and supports a rich nomenclature for defining complex motifs.

The tool is organized into two main functions:
1.  **Convert**: Processes various biological file formats (`.pdb`, `.cif`, `.fasta`) into a standardized CSV format containing protein sequences.
2.  **Search**: Scans the sequence CSV file to find all occurrences of specified motifs.

### Features

*   **Efficient Searching**: Utilizes the Aho-Corasick algorithm to find multiple motifs in a single pass.
*   **Flexible Motif Definition**: Supports a rich nomenclature including wildcards, custom character sets, exclusions, and biochemical property groups.
*   **Multiple Input Formats**: Converts `.pdb`, `.cif`, and `.fasta` files (including DNA sequences, which are translated into proteins) into a uniform format.
*   **Modular and Extensible**: The logic is separated into distinct scripts for file conversion and motif searching, orchestrated by a central command-line interface.
*   **Detailed Outputs**: Generates an aggregate CSV file with all results, plus individual JSON files for each sequence with detailed match information.
*   **Uniprot API**: Integration with the UniProt API for fast retrieval of protein sequences with search features.

### Project Structure

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

### Installation

This project requires a few external Python libraries. You can install them using pip:

```bash
pip install pandas biopython biotite requests numpy
```
or:
```bash
pip install -r requirements.txt
```


### How to Use

The tool is controlled through `main.py` and has two primary commands: `convert` and `search`. You must choose one of them.

#### `convert` (Optional)

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

#### `search`

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
| `--motif_name_column` | `-mnc` | No       | `motif_name`               | Name of the column in the motif file that contains the motif names.    |
| `--sequences`       |      | **Yes**  | `sequences.csv`            | Path to the CSV file containing the protein sequences.               |
| `--sequence_column` | `-sc`| **Yes**  |                            | Name of the column in the sequences file that contains the sequences.  |
| `--output`          |      | No       | `motif_search_results.csv` | Path for the output CSV file that will store the aggregate results.  |

#### `uniprot` (optional)

You can use the UniProt API to retrieve sequences. 

**Usage:**
```bash
python main.py uniprot [OPTIONS]
```

**Arguments for `uniprot`:**

| Argument | Required | Default Value | Description |
| --- | --- | --- | --- |
| `--query` | No | | A general UniProt query string. |
| `--organism` | No | | Filter by organism name. |
| `--enzyme_family` | No | | Filter by enzyme family. |
| `--protein_name` | No | | Filter by protein name. |
| `--accession` | No | | Filter by accession number. |
| `--output_csv` | No | `sequences.csv` | The name of the output file. |
| `--limit` | No | 500 | The maximum number of sequences to retrieve. |


### Input File Formats

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

### Motif Nomenclature

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

### Output Format

The `search` command produces two types of output:

1.  **Aggregate CSV File (`--output`)**: A single CSV file summarizing all findings.
    *   `name`: The name of the sequence from the input file.
    *   `sequence`: The full protein sequence.
    *   `motifs`: A JSON string containing a dictionary where keys are the original motifs. Each value is a dictionary containing the `motif_name` and a list of `matches`. Each match includes the end position and the specific concrete motif that matched.

2.  **Individual JSON Files**: For each sequence processed, a detailed JSON file is saved in a timestamped directory (e.g., `20250725_163000_motif_search_results_jsons/`).
    *   `name`: The name of the sequence.
    *   `sequence`: The full protein sequence.
    *   `results`: A dictionary where keys are the original motifs found. The values are dictionaries containing the `motif_name` and a list of `matches`, where each match is a `[end_position, "concrete_motif"]` pair.
```

## Structural Motif Search

This tool searches for 3D structural motifs in protein structures (`.pdb`, `.cif`).

### How to Use

**Usage:**
```bash
python structure_motif/search_3d_motif.py [OPTIONS]
```

**Arguments:**

| Argument | Flag | Required | Description |
| --- | --- | --- | --- |
| `--input_folder` | `-i` | **Yes** | Folder with PDB/CIF files. |
| `--motif_file` | `-m` | **Yes** | JSON motif definition file. |
| `--output_folder` | `-o` | **Yes** | Folder to save JSON results. |
| `--summary_csv` | `-s` | **Yes** | Final summary CSV file. |

### Example Usage

To run the structural motif search, you need to provide the input folder containing your protein structure files, the motif definition file, the output folder for the detailed JSON results, and the summary CSV file.

Here is an example command:

```bash
python structure_motif/search_3d_motif.py \
    --input_folder protein_files \
    --motif_file structure_motif/motifs/catalytic_triad.json \
    --output_folder outputs \
    --summary_csv summary.csv
```

This command will:
- Search for the catalytic triad motif defined in `structure_motif/motifs/catalytic_triad.json`.
- Look for matching structures in all `.pdb` and `.cif` files within the `protein_files` directory.
- Save a detailed JSON file for each input structure in the `outputs` directory.
- Generate a `summary.csv` file containing a summary of all found motifs.



### Structural Motif Definition

Structural motifs are defined in a JSON file. The format is described in detail in `structure_motif/motifs/motif_format_documentation.md`.

Here is an example of a catalytic triad motif:

```json
{
  "motif_name": "Catalytic Triad",
  "description": "A classic Ser-His-Asp catalytic triad.",
  "components": [
    {
      "id": "ser",
      "residue_type": "SER",
      "atom_selectors": {
        "hydroxyl_oxygen": "OG",
        "beta_carbon": "CB"
      }
    },
    {
      "id": "his",
      "residue_type": "HIS",
      "atom_selectors": {
        "imidazole_nitrogen_delta": "ND1"
      }
    },
    {
      "id": "asp",
      "residue_type": "ASP",
      "atom_selectors": {
        "carboxyl_oxygen_delta1": "OD1"
      }
    }
  ],
  "constraints": [
    {
      "type": "distance",
      "atoms": ["ser.hydroxyl_oxygen", "his.imidazole_nitrogen_delta"],
      "value": 3.0,
      "tolerance": 0.5
    },
    {
      "type": "distance",
      "atoms": ["his.imidazole_nitrogen_delta", "asp.carboxyl_oxygen_delta1"],
      "value": 3.0,
      "tolerance": 0.5
    }
  ]
}
```

### Output Format

The structural motif search produces two types of output:

1.  **Summary CSV File (`--summary_csv`)**: A single CSV file summarizing all findings.
    *   `source_file`: The name of the PDB/CIF file.
    *   `motif_id`: A unique identifier for the found motif.
    *   `residue_1`, `residue_2`, ...: The residues that form the motif, in the format `RES-CHAIN-RESID` (e.g., `SER-A-123`).

2.  **Individual JSON Files (`--output_folder`)**: For each input structure, a detailed JSON file is saved.
    *   `source_file`: The name of the PDB/CIF file.
    *   `motifs_found`: The number of motifs found in the file.
    *   `matches`: A list of found motifs, where each motif is a list of residues with their details (name, chain, and residue ID).

*Example JSON output (`1AQ7.json`):*
```json
{
  "source_file": "1AQ7.cif",
  "motifs_found": 1,
  "matches": [
    {
      "residues": [
        {
          "res_name": "SER",
          "chain_id": "A",
          "res_id": 195
        },
        {
          "res_name": "HIS",
          "chain_id": "A",
          "res_id": 57
        },
        {
          "res_name": "ASP",
          "chain_id": "A",
          "res_id": 102
        }
      ]
    }
  ]
}
```
### Motif Format Documentation
For additional info on parsing motifs structurally using constraints such as bond angles, dihedrals, radius exclusion spheres, distance, and secondary structure, see ```motif_format_documentation.md``` for help with JSON motif definition formatting.

## Benchmarking
Due to its efficient use of an Aho-Corasick search algorithm to locate sequence motifs, the sequence motif search looked for 151,394 unique kinase consensus motifs on 100 substrate proteins of varying lengths (UniProt query: human kinase substrate protein) in ~57 seconds on M3 Macbook Air.
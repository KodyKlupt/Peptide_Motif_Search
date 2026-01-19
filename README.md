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
pip install pandas biopython biotite requests numpy ipykernel py3Dmol
```
or:
```bash
pip install -r requirements.txt
```

#### Conda Environment Setup (Recommended)

To create an isolated conda environment and register it as a Jupyter kernel:

```bash
# Create a new conda environment
conda create -n motifsearch python=3.11 -y

# Activate the environment
conda activate motifsearch

# Install dependencies from requirements.txt
pip install -r requirements.txt

# Register the environment as a Jupyter kernel
python -m ipykernel install --user --name motifsearch --display-name "Motif Search"
```

After this, "Motif Search" will appear as a kernel option in Jupyter notebooks.


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

---

## Command Line Examples

This section provides extensive examples demonstrating the full capabilities of both the sequence and structural motif search tools.

### File Conversion Examples

#### Basic Conversion
```bash
# Convert all protein files in the default folder to sequences.csv
python main.py convert

# Convert files from a custom folder
python main.py convert --input_folder my_protein_data

# Convert to a custom output file
python main.py convert --input_folder pdb_structures --output_csv my_sequences.csv

# Convert a folder of FASTA files
python main.py convert --input_folder fasta_files --output_csv fasta_sequences.csv

# Convert mixed file types (PDB, CIF, FASTA)
python main.py convert --input_folder mixed_structures --output_csv all_sequences.csv
```

---

### UniProt Query Examples

#### Basic Queries
```bash
# Search for human proteins
python main.py uniprot --organism "Homo sapiens" --output_csv human_proteins.csv

# Search for kinases
python main.py uniprot --enzyme_family kinase --output_csv kinases.csv

# Search for a specific protein by name
python main.py uniprot --protein_name "insulin receptor" --output_csv insulin_receptors.csv

# Search by accession number
python main.py uniprot --accession P00533 --output_csv egfr.csv
```

#### Combined Filters
```bash
# Human kinases
python main.py uniprot --organism "Homo sapiens" --enzyme_family kinase --output_csv human_kinases.csv

# Mouse tyrosine kinases
python main.py uniprot --organism "Mus musculus" --protein_name "tyrosine kinase" --output_csv mouse_tyr_kinases.csv

# E. coli proteases
python main.py uniprot --organism "Escherichia coli" --enzyme_family protease --output_csv ecoli_proteases.csv

# Human dehydrogenases with limit
python main.py uniprot --organism "Homo sapiens" --enzyme_family dehydrogenase --limit 100 --output_csv human_dehydrogenases.csv
```

#### Advanced Queries
```bash
# Custom query string for complex searches
python main.py uniprot --query "(organism:human) AND (keyword:kinase) AND (length:[100 TO 500])" --output_csv custom_query.csv

# Retrieve specific enzyme classes
python main.py uniprot --query "ec:2.7.1.*" --limit 1000 --output_csv transferases.csv

# Search for reviewed entries only
python main.py uniprot --query "(reviewed:true) AND (organism:9606)" --output_csv reviewed_human.csv

# Search membrane proteins
python main.py uniprot --query "(keyword:membrane) AND (organism:human)" --limit 250 --output_csv membrane_proteins.csv
```

---

### Sequence Motif Search Examples

#### Basic Searches
```bash
# Search using default files
python main.py search --motifs motifs.csv --motif_column pattern --sequences sequences.csv --sequence_column sequence

# Search with explicit column names
python main.py search \
    --motifs motif_libraries/kinase_motifs.csv \
    --motif_column motif_pattern \
    --motif_name_column motif_id \
    --sequences sequences.csv \
    --sequence_column sequence \
    --output kinase_motif_results.csv
```

#### Searching Specific Motif Types
```bash
# Search for phosphorylation sites (simple)
python main.py search \
    --motifs motif_libraries/phospho_sites.csv \
    --motif_column motif \
    --sequences human_kinases.csv \
    --sequence_column sequence \
    --output phospho_results.csv

# Search for SH2 domain binding motifs
python main.py search \
    --motifs motif_libraries/sh2_motifs.csv \
    --motif_column consensus \
    --sequences sequences.csv \
    --sequence_column sequence \
    --output sh2_binding_results.csv

# Search for SH3 domain binding motifs
python main.py search \
    --motifs motif_libraries/sh3_motifs.csv \
    --motif_column pattern \
    --sequences my_proteins.csv \
    --sequence_column seq \
    --output sh3_results.csv
```

#### Using Different Motif Nomenclatures
```bash
# Wildcard motifs (x = any amino acid)
# Example motif: "RxxS" matches R-any-any-S
python main.py search \
    --motifs motif_libraries/simple_motifs.csv \
    --motif_column wildcard_motif \
    --sequences sequences.csv \
    --sequence_column sequence \
    --output wildcard_results.csv

# Custom character sets (brackets)
# Example motif: "R[ST]xP" matches R-S/T-any-P
python main.py search \
    --motifs motif_libraries/kinase_consensus.csv \
    --motif_column consensus \
    --sequences substrates.csv \
    --sequence_column sequence \
    --output kinase_sites.csv

# Exclusion patterns (curly braces)
# Example motif: "{P}G" matches any-except-P followed by G
python main.py search \
    --motifs motif_libraries/exclusion_motifs.csv \
    --motif_column pattern \
    --sequences sequences.csv \
    --sequence_column sequence \
    --output exclusion_results.csv

# Hydrophobic residues (%)
# Example motif: "%xL" matches hydrophobic-any-L
python main.py search \
    --motifs motif_libraries/hydrophobic_motifs.csv \
    --motif_column motif \
    --sequences membrane_proteins.csv \
    --sequence_column sequence \
    --output hydrophobic_results.csv

# Aromatic residues (@)
# Example motif: "@S" matches aromatic-S
python main.py search \
    --motifs motif_libraries/aromatic_motifs.csv \
    --motif_column pattern \
    --sequences sequences.csv \
    --sequence_column sequence \
    --output aromatic_results.csv

# Polar residues (&)
# Example motif: "&xT" matches polar-any-T
python main.py search \
    --motifs motif_libraries/polar_motifs.csv \
    --motif_column motif \
    --sequences sequences.csv \
    --sequence_column sequence \
    --output polar_results.csv

# Positively charged residues ([+])
# Example motif: "[+]GxE" matches K/R/H-G-any-E
python main.py search \
    --motifs motif_libraries/charged_motifs.csv \
    --motif_column pattern \
    --sequences sequences.csv \
    --sequence_column sequence \
    --output positive_charge_results.csv

# Negatively charged residues ([-])
# Example motif: "[-]xL" matches D/E-any-L
python main.py search \
    --motifs motif_libraries/negative_motifs.csv \
    --motif_column pattern \
    --sequences sequences.csv \
    --sequence_column sequence \
    --output negative_charge_results.csv

# Aliphatic residues (#)
# Example motif: "#xP" matches A/V/L/I-any-P
python main.py search \
    --motifs motif_libraries/aliphatic_motifs.csv \
    --motif_column motif \
    --sequences sequences.csv \
    --sequence_column sequence \
    --output aliphatic_results.csv

# Small residues (~)
# Example motif: "~G~" matches small-G-small
python main.py search \
    --motifs motif_libraries/small_motifs.csv \
    --motif_column pattern \
    --sequences sequences.csv \
    --sequence_column sequence \
    --output small_residue_results.csv
```

#### Complex Motif Patterns
```bash
# Combined nomenclature: hydrophobic + custom set + any
# Motif: "%[ST]x[+]" matches hydrophobic-S/T-any-positive
python main.py search \
    --motifs motif_libraries/complex_motifs.csv \
    --motif_column complex_pattern \
    --sequences kinase_substrates.csv \
    --sequence_column sequence \
    --output complex_results.csv

# Multiple exclusions with wildcards
# Motif: "{PG}xx{DE}K" matches not-P/G, any, any, not-D/E, K
python main.py search \
    --motifs motif_libraries/exclusion_complex.csv \
    --motif_column motif \
    --sequences sequences.csv \
    --sequence_column sequence \
    --output multi_exclusion_results.csv

# Long patterns with mixed nomenclature
# Motif: "Rx[ST]%x#[+]" 
python main.py search \
    --motifs motif_libraries/long_motifs.csv \
    --motif_column extended_pattern \
    --sequences proteome.csv \
    --sequence_column sequence \
    --output long_pattern_results.csv
```

#### Batch Processing Examples
```bash
# Search multiple organisms
for organism in "human" "mouse" "zebrafish"; do
    python main.py uniprot --organism "$organism" --output_csv ${organism}_seqs.csv
    python main.py search \
        --motifs kinase_motifs.csv \
        --motif_column pattern \
        --sequences ${organism}_seqs.csv \
        --sequence_column sequence \
        --output ${organism}_kinase_results.csv
done

# Search with multiple motif libraries
for lib in phospho sumo ubiquitin glycosyl; do
    python main.py search \
        --motifs motif_libraries/${lib}_motifs.csv \
        --motif_column pattern \
        --sequences human_proteome.csv \
        --sequence_column sequence \
        --output ${lib}_sites.csv
done
```

---

### Structural Motif Search Examples

#### Basic Structure Searches
```bash
# Search for catalytic triad
python structure_motif/search_3d_motif.py \
    --input_folder protein_files \
    --motif_file structure_motif/motifs/catalytic_triad.json \
    --output_folder outputs/catalytic_triad \
    --summary_csv catalytic_triad_summary.csv

# Search for zinc finger motif
python structure_motif/search_3d_motif.py \
    --input_folder pdb_structures \
    --motif_file structure_motif/motifs/zinc_finger.json \
    --output_folder outputs/zinc_finger \
    --summary_csv zinc_finger_summary.csv

# Search for EF-hand calcium binding motif
python structure_motif/search_3d_motif.py \
    --input_folder ca_binding_proteins \
    --motif_file structure_motif/motifs/ef_hand.json \
    --output_folder outputs/ef_hand \
    --summary_csv ef_hand_summary.csv
```

#### Using Short Flags
```bash
# Compact notation for catalytic triad search
python structure_motif/search_3d_motif.py \
    -i protein_files \
    -m structure_motif/motifs/catalytic_triad.json \
    -o outputs \
    -s summary.csv

# Compact notation for custom motif
python structure_motif/search_3d_motif.py \
    -i my_structures \
    -m my_motif.json \
    -o results \
    -s results.csv
```

#### Different Constraint Types
```bash
# Distance-only constraints
python structure_motif/search_3d_motif.py \
    -i protein_files \
    -m structure_motif/motifs/distance_only_motif.json \
    -o outputs/distance_search \
    -s distance_results.csv

# Angle constraints (e.g., for metal coordination)
python structure_motif/search_3d_motif.py \
    -i metalloenzymes \
    -m structure_motif/motifs/metal_coordination.json \
    -o outputs/metal_coord \
    -s metal_coordination_results.csv

# Dihedral constraints (e.g., for backbone conformations)
python structure_motif/search_3d_motif.py \
    -i loop_structures \
    -m structure_motif/motifs/loop_conformation.json \
    -o outputs/loops \
    -s loop_results.csv

# Secondary structure constraints
python structure_motif/search_3d_motif.py \
    -i all_structures \
    -m structure_motif/motifs/helix_turn_helix.json \
    -o outputs/hth \
    -s helix_turn_helix_results.csv

# Accessibility constraints (surface exposed residues)
python structure_motif/search_3d_motif.py \
    -i enzyme_structures \
    -m structure_motif/motifs/surface_active_site.json \
    -o outputs/surface_sites \
    -s surface_site_results.csv

# Exclusion sphere constraints (no atoms within radius)
python structure_motif/search_3d_motif.py \
    -i protein_files \
    -m structure_motif/motifs/buried_residue.json \
    -o outputs/buried \
    -s buried_residue_results.csv
```

#### Enzyme Active Site Searches
```bash
# Serine proteases
python structure_motif/search_3d_motif.py \
    -i protease_structures \
    -m structure_motif/motifs/serine_protease.json \
    -o outputs/serine_protease \
    -s serine_protease_results.csv

# Cysteine proteases
python structure_motif/search_3d_motif.py \
    -i protease_structures \
    -m structure_motif/motifs/cysteine_protease.json \
    -o outputs/cysteine_protease \
    -s cysteine_protease_results.csv

# Aspartyl proteases
python structure_motif/search_3d_motif.py \
    -i protease_structures \
    -m structure_motif/motifs/aspartyl_protease.json \
    -o outputs/aspartyl_protease \
    -s aspartyl_protease_results.csv

# Metalloprotease active sites
python structure_motif/search_3d_motif.py \
    -i metalloenzyme_structures \
    -m structure_motif/motifs/metalloprotease.json \
    -o outputs/metalloprotease \
    -s metalloprotease_results.csv

# Kinase active sites
python structure_motif/search_3d_motif.py \
    -i kinase_structures \
    -m structure_motif/motifs/kinase_active_site.json \
    -o outputs/kinase_sites \
    -s kinase_active_site_results.csv

# Phosphatase active sites
python structure_motif/search_3d_motif.py \
    -i phosphatase_structures \
    -m structure_motif/motifs/phosphatase.json \
    -o outputs/phosphatase \
    -s phosphatase_results.csv
```

#### Binding Site Searches
```bash
# ATP binding sites
python structure_motif/search_3d_motif.py \
    -i atp_binding_proteins \
    -m structure_motif/motifs/atp_binding_site.json \
    -o outputs/atp_binding \
    -s atp_binding_results.csv

# NAD binding (Rossmann fold)
python structure_motif/search_3d_motif.py \
    -i oxidoreductases \
    -m structure_motif/motifs/rossmann_fold.json \
    -o outputs/rossmann \
    -s rossmann_fold_results.csv

# DNA binding motifs
python structure_motif/search_3d_motif.py \
    -i transcription_factors \
    -m structure_motif/motifs/dna_binding.json \
    -o outputs/dna_binding \
    -s dna_binding_results.csv

# Ligand binding pocket
python structure_motif/search_3d_motif.py \
    -i receptor_structures \
    -m structure_motif/motifs/ligand_pocket.json \
    -o outputs/ligand_pockets \
    -s ligand_pocket_results.csv
```

#### Structural Feature Searches
```bash
# Disulfide bridges
python structure_motif/search_3d_motif.py \
    -i protein_files \
    -m structure_motif/motifs/disulfide_bond.json \
    -o outputs/disulfide \
    -s disulfide_results.csv

# Salt bridges
python structure_motif/search_3d_motif.py \
    -i protein_files \
    -m structure_motif/motifs/salt_bridge.json \
    -o outputs/salt_bridges \
    -s salt_bridge_results.csv

# Hydrogen bond networks
python structure_motif/search_3d_motif.py \
    -i protein_files \
    -m structure_motif/motifs/hbond_network.json \
    -o outputs/hbond_networks \
    -s hbond_network_results.csv

# Buried ionizable residues
python structure_motif/search_3d_motif.py \
    -i protein_files \
    -m structure_motif/motifs/buried_ionizable.json \
    -o outputs/buried_ionizable \
    -s buried_ionizable_results.csv

# Aromatic clusters
python structure_motif/search_3d_motif.py \
    -i protein_files \
    -m structure_motif/motifs/aromatic_cluster.json \
    -o outputs/aromatic_clusters \
    -s aromatic_cluster_results.csv
```

#### Batch Structure Processing
```bash
# Search multiple motifs on same structures
for motif in catalytic_triad zinc_finger ef_hand disulfide_bond; do
    python structure_motif/search_3d_motif.py \
        -i protein_files \
        -m structure_motif/motifs/${motif}.json \
        -o outputs/${motif} \
        -s ${motif}_summary.csv
done

# Search same motif on different structure sets
for pdb_set in kinases proteases phosphatases; do
    python structure_motif/search_3d_motif.py \
        -i structures/${pdb_set} \
        -m structure_motif/motifs/catalytic_triad.json \
        -o outputs/${pdb_set}_catalytic \
        -s ${pdb_set}_catalytic_results.csv
done

# Process entire PDB mirror
python structure_motif/search_3d_motif.py \
    -i /data/pdb_mirror \
    -m structure_motif/motifs/catalytic_triad.json \
    -o /data/results/catalytic_triad \
    -s /data/results/catalytic_triad_pdb_wide.csv
```

---

### Combined Workflow Examples

#### Complete Analysis Pipeline
```bash
# Step 1: Get human kinases from UniProt
python main.py uniprot \
    --organism "Homo sapiens" \
    --enzyme_family kinase \
    --limit 1000 \
    --output_csv human_kinases.csv

# Step 2: Search for phosphorylation consensus motifs
python main.py search \
    --motifs motif_libraries/phospho_consensus.csv \
    --motif_column pattern \
    --sequences human_kinases.csv \
    --sequence_column sequence \
    --output kinase_phospho_sites.csv

# Step 3: If you have structures, search for structural motifs
python structure_motif/search_3d_motif.py \
    -i kinase_structures \
    -m structure_motif/motifs/kinase_active_site.json \
    -o outputs/kinase_active_sites \
    -s kinase_active_site_summary.csv
```

#### Drug Target Analysis
```bash
# Get GPCR sequences
python main.py uniprot \
    --query "(family:GPCR) AND (organism:human)" \
    --limit 500 \
    --output_csv gpcrs.csv

# Search for conserved motifs
python main.py search \
    --motifs motif_libraries/gpcr_motifs.csv \
    --motif_column consensus \
    --sequences gpcrs.csv \
    --sequence_column sequence \
    --output gpcr_motif_analysis.csv

# Structural analysis of binding pockets
python structure_motif/search_3d_motif.py \
    -i gpcr_structures \
    -m structure_motif/motifs/gpcr_binding_pocket.json \
    -o outputs/gpcr_pockets \
    -s gpcr_pocket_analysis.csv
```
---

### Quick Reference

| Task | Command |
|------|---------|
| Convert PDB files | `python main.py convert --input_folder pdbs --output_csv seqs.csv` |
| Get UniProt sequences | `python main.py uniprot --organism human --output_csv human.csv` |
| Search sequence motifs | `python main.py search --motifs m.csv -mc pattern --sequences s.csv -sc seq` |
| Search structural motifs | `python structure_motif/search_3d_motif.py -i pdbs -m motif.json -o out -s sum.csv` |


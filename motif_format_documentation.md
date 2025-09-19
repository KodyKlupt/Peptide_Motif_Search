# Structural Motif JSON Format Documentation

This document details the JSON format used to define structural motifs for searching in protein structures.

## Root Level Fields

The JSON file has the following top-level fields:

- `motif_name` (string, required): A descriptive name for the motif.
- `description` (string, required): A brief description of the motif, its function, or where it's found.
- `components` (array of objects, required): Defines the individual residues that constitute the motif.
- `constraints` (array of objects, required): Defines the geometric and structural relationships between the components.

## Components

Each object in the `components` array defines a single residue of the motif and has the following fields:

- `id` (string, required): A unique identifier for the residue within the motif (e.g., "ser", "his1"). This is used to reference the component in the `constraints` section.
- `residue_type` (string, required): The three-letter code for the amino acid residue (e.g., "SER", "HIS", "ASP").
- `atom_selectors` (object, required): A dictionary that maps user-defined names to standard PDB atom names. This provides a convenient way to refer to specific atoms in the `constraints`.
  - **key** (string): A user-defined name for the atom (e.g., "hydroxyl_oxygen").
  - **value** (string): The standard PDB atom name (e.g., "OG", "CA", "ND1").

### Example Component

```json
{
  "id": "ser",
  "residue_type": "SER",
  "atom_selectors": {
    "hydroxyl_oxygen": "OG",
    "alpha_carbon": "CA"
  }
}
```

## Constraints

The `constraints` array defines the spatial and structural properties of the motif. Each constraint is an object with a `type` field that determines the other required fields.

### Constraint Types

#### 1. `distance`

Measures the distance between two atoms.

- `type`: "distance"
- `atoms` (array of strings, required): An array containing two strings, each identifying an atom in the format `component_id.atom_selector`.
- `value` (float, required): The ideal distance in angstroms.
- `tolerance` (float, required): The allowed deviation from `value` in angstroms.

**Example:**

```json
{
  "type": "distance",
  "atoms": ["ser.hydroxyl_oxygen", "his.imidazole_nitrogen_delta"],
  "value": 4.8,
  "tolerance": 0.6
}
```

#### 2. `angle`

Measures the angle formed by three atoms.

- `type`: "angle"
- `atoms` (array of strings, required): An array of three strings identifying the atoms that form the angle. The angle is measured at the second atom.
- `value` (float, required): The ideal angle in degrees.
- `tolerance` (float, required): The allowed deviation from `value` in degrees.

**Example:**

```json
{
  "type": "angle",
  "atoms": ["his.imidazole_nitrogen_delta", "ser.hydroxyl_oxygen", "ser.beta_carbon"],
  "value": 89.0,
  "tolerance": 15.0
}
```

#### 3. `dihedral`

Measures the dihedral (torsion) angle formed by four atoms.

- `type`: "dihedral"
- `atoms` (array of strings, required): An array of four strings identifying the atoms that form the dihedral angle.
- `value` (float, required): The ideal dihedral angle in degrees.
- `tolerance` (float, required): The allowed deviation from `value` in degrees.

**Example:**

```json
{
  "type": "dihedral",
  "atoms": ["his.alpha_carbon", "his.beta_carbon", "his.gamma_carbon", "his.imidazole_nitrogen_delta"],
  "value": -109.0,
  "tolerance": 20.0
}
```

#### 4. `secondary_structure`

Constrains the secondary structure of a specific residue.

- `type`: "secondary_structure"
- `component_id` (string, required): The `id` of the component to constrain.
- `value` (array of strings, required): A list of allowed secondary structure types. Common types include:
  - "H": Alpha-helix
  - "B": Beta-bridge
  - "E": Strand
  - "G": 3-10 helix
  - "I": Pi-helix
  - "T": Turn
  - "S": Bend
  - "-": Other/None

**Example:**

```json
{
  "type": "secondary_structure",
  "component_id": "ser",
  "value": ["T", "S", "-"]
}
```

#### 5. `accessibility`

Constrains the solvent accessibility of a residue.

- `type`: "accessibility"
- `component_id` (string, required): The `id` of the component to constrain.
- `comparison` (string, required): The comparison operator. Can be "greater_than" or "less_than".
- `value` (float, required): The relative solvent accessibility (RSA) value to compare against (a value between 0.0 and 1.0).

**Example:**

```json
{
  "type": "accessibility",
  "component_id": "his",
  "comparison": "greater_than",
  "value": 0.2
}
```

#### 6. `exclusion_sphere`

Defines a sphere around a specific atom that must not contain any other atoms from the protein (other than those belonging to the same residue). This is useful for ensuring a binding pocket or active site is open.

- `type`: "exclusion_sphere"
- `component_id` (string, required): The `id` of the component containing the central atom of the sphere.
- `atom_selector` (string, required): The user-defined name of the atom at the center of the sphere.
- `radius` (float, required): The radius of the exclusion sphere in angstroms.

**Example:**

```json
{
  "type": "exclusion_sphere",
  "component_id": "ser",
  "atom_selector": "hydroxyl_oxygen",
  "radius": 2.5
}
```

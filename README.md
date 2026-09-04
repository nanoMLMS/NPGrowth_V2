# NPGrowth V2

NPGrowth V2 provides a workflow for generating and running nanoparticle deposition simulations with LAMMPS.

The project combines:

* a modified implementation of the LAMMPS `fix deposit` command;
* a Python tool for generating deposition instructions from JSON configuration files;
* LAMMPS input files for running nanoparticle growth simulations;
* nanoparticle seed structures and interatomic potentials.

The modified deposition command allows particles to be generated inside a **source region** and directed toward one or more user-defined **target regions**.

---

## Repository structure

```text
NPGrowth_V2/
├── README.md
├── INSTALL.md
│
├── deposition_generator/
│   ├── generate_depo_input.py
│   ├── configs/
│   │   └── 1source_1targets.json
│   └── generated_inputs/
│       └── 1source_1targets.in
│
├── lammps_run/
│   ├── runs/
│   │   └── 1source_1target/
│   │       ├── lammps_input.in
│   │       └── run.sh
│   ├── seeds/
│   └── potentials/
│
└── lammps_source_modified/
    └── lammps_deposition_modified/
        ├── fix_deposit.cpp
        ├── fix_deposit.h
        └── cmake_command
```

Installation and compilation of the modified LAMMPS version are described in [`INSTALL.md`](INSTALL.md).

---

# Workflow

The basic workflow is:

```text
JSON deposition configuration
            │
            ▼
generate_depo_input.py
            │
            ▼
generated LAMMPS deposition input
            │
            ▼
main LAMMPS input file
            │
            ▼
modified LAMMPS executable
            │
            ▼
simulation
```

The deposition geometry is therefore kept separate from the rest of the LAMMPS simulation setup.

---

# 1. Deposition configuration

Deposition configurations are stored in:

```text
deposition_generator/configs/
```

A minimal configuration has the following structure:

```json
{
  "n_source_regions": 1,
  "sources": [
    {
      "source_region": [
        30.659863,
        34.659863,
        -2.0,
        2.0,
        21.094011,
        25.094011
      ],

      "n_target_regions": 1,

      "targets": [
        {
          "size": [4, 4, 4],
          "position": [0, -11.04, 12.88]
        }
      ],

      "type": "Cu",
      "velocity_range": [5.0, 8.0],
      "near": 4.0,
      "n_depo": 100,
      "every": 5000,
      "random_seed_LAMMPS": 42
    }
  ]
}
```

## Source regions

`n_source_regions` specifies the number of independent deposition sources.

Each source is defined inside the `sources` list.

The source geometry is specified by:

```json
"source_region": [xlo, xhi, ylo, yhi, zlo, zhi]
```

and corresponds to a LAMMPS block region:

```text
region source_box_i block xlo xhi ylo yhi zlo zhi
```

The deposition position is randomly selected inside this region.

Multiple source regions can be defined in the same configuration. Each source has its own target geometry and deposition parameters.

---

# 2. Target regions

Each source contains one or more target regions.

The number of targets is specified with:

```json
"n_target_regions": 1
```

and the corresponding geometries are listed in:

```json
"targets": [...]
```

Each target is an axis-aligned cuboid defined by:

```json
{
  "size": [sx, sy, sz],
  "position": [x, y, z]
}
```

`size` contains the dimensions of the cuboid.

`position` contains the coordinates of its **center**.

For example:

```json
{
  "size": [4, 4, 4],
  "position": [0, 0, 0]
}
```

generates a region extending from:

```text
x = -2 ... 2
y = -2 ... 2
z = -2 ... 2
```

If several targets are specified for the same source, the generator creates the individual LAMMPS regions and combines them into a `union` region.

The union therefore represents the complete target geometry associated with that source.

---

# 3. Deposition parameters

Each source has its own deposition parameters.

### `type`

```json
"type": "Cu"
```

Atomic type or label passed to the LAMMPS deposition command.

### `velocity_range`

```json
"velocity_range": [5.0, 8.0]
```

Defines the range used to generate the initial deposition velocity.

The generated input currently passes this range through the `vx` option of `fix deposit`.

When a target region is specified, the modified `fix deposit` implementation redirects the generated velocity toward a randomly selected point inside the target region while preserving the magnitude of the original velocity vector.

### `near`

```json
"near": 4.0
```

Minimum allowed distance used by `fix deposit` when attempting an insertion.

### `n_depo`

```json
"n_depo": 100
```

Number of deposition events requested for that source.

### `every`

```json
"every": 5000
```

Number of timesteps between deposition attempts.

### `random_seed_LAMMPS`

```json
"random_seed_LAMMPS": 42
```

Random-number seed passed to LAMMPS.

Different independent sources should normally use different seeds.

---

# 4. Generate the LAMMPS deposition input

From the repository root, run:

```bash
python deposition_generator/generate_depo_input.py \
    --config deposition_generator/configs/1source_1targets.json \
    --output deposition_generator/generated_inputs/1source_1targets.in
```

The generator first validates the configuration and then creates the corresponding LAMMPS regions and deposition fixes.

A successful execution prints messages similar to:

```text
Configuration OK: 1 source regions found.
LAMMPS deposition input written to: deposition_generator/generated_inputs/1source_1targets.in
```

For a single source and target, the generated file has the general form:

```text
region source_box_0 block ...
region target_0_0 block ...

fix deposit_0 added deposit ... region source_box_0 ... target target_0_0
```

For multiple targets, an additional union is generated:

```text
region target_box_0 union N target_0_0 target_0_1 ...
```

---

# 5. Prepare a LAMMPS run

Simulation cases are stored under:

```text
lammps_run/runs/
```

For example:

```text
lammps_run/runs/1source_1target/
├── lammps_input.in
└── run.sh
```

The main `lammps_input.in` defines the complete simulation, including:

* simulation box;
* initial nanoparticle;
* atom groups;
* interatomic potential;
* thermostat and integration settings;
* trajectory output;
* deposition input;
* total simulation time.

The generated deposition file is included from the main LAMMPS input with:

```text
include /path/to/deposition_generator/generated_inputs/1source_1targets.in
```

Before running a simulation, update paths in `lammps_input.in` to match your local repository location.

In particular, check the paths used for:

```text
read_dump
pair_coeff
include
```

---

# 6. Seeds and potentials

Initial nanoparticle structures are stored in:

```text
lammps_run/seeds/
```

The example simulation reads an XYZ nanoparticle structure using `read_dump`.

Interatomic potentials are stored in:

```text
lammps_run/potentials/
```

The current Cu example uses:

```text
Cu_u3.eam
```

with:

```text
pair_style eam
```

Users should modify the potential, atom types, masses and related LAMMPS settings when simulating different materials.

---

# 7. Run the simulation

After compiling the modified version of LAMMPS as described in `INSTALL.md`, enter a run directory:

```bash
cd lammps_run/runs/1source_1target
```

and execute:

```bash
bash run.sh
```

The current example uses MPI with four processes:

```bash
mpirun -np 4 /path/to/lammps/build/lmp -in lammps_input.in
```

Change the executable path and number of MPI processes according to your installation and computational environment.

---

# Modified `fix deposit`

This project uses a modified version of the LAMMPS `fix deposit` implementation.

The relevant source files are:

```text
lammps_source_modified/lammps_deposition_modified/
├── fix_deposit.cpp
└── fix_deposit.h
```

The modification introduces the use of a target region in the deposition process.

Conceptually:

```text
source region
     │
     │ insertion position
     ▼
     ●
      \
       \ velocity
        \
         ▼
   target region
```

For each deposition event, an insertion position is selected in the source region and a target point is selected inside the specified target region.

The particle velocity is then directed from the insertion position toward that target point while retaining the magnitude of the initially generated velocity.

See [`INSTALL.md`](INSTALL.md) for instructions on compiling LAMMPS with the modified source files.

---

# Notes

* Source and target regions must lie inside the LAMMPS simulation box.
* `n_source_regions` must match the number of entries in `sources`.
* For each source, `n_target_regions` must match the number of entries in `targets`.
* Each `source_region` must contain exactly six values.
* Each target `size` and `position` must contain exactly three values.
* Target sizes must be positive.
* The generated deposition input is intended to be included in a complete LAMMPS input file rather than used as a complete simulation by itself.


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
│   │   ├── 1source_1targets.json
│   │   └── 4source_tetra_sintering.json
│   └── generated_inputs/
│       ├── 1source_1targets.in
│       └── 4source_tetra_zerovel.in
│
├── lammps_run/
│   ├── runs/
│   │   ├── 1source_1target/
│   │   │   ├── lammps_input.in
│   │   │   └── run.sh
│   │   └── sintering_tetra/
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

Two configuration modes are supported by the same generator:

* the standard **source/target mode** (sections 1-3), where you hand-specify the source and target regions and the deposition proceeds incrementally (nanoparticle growth);
* the **geometry/sintering mode** (section 4), where you give a set of equilibrated seed clusters and a geometry, and the generator derives the source placements automatically so the clusters are inserted at rest as a sintering-start configuration.

The two modes can be mixed freely: they are only two alternative input schemas for `generate_depo_input.py`.

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

# Alternative usage: initial configurations for sintering experiments

Instead of hand-specifying each source region, you can let the generator place a set of **equilibrated seed clusters** at a chosen geometry and separation. This is a convenient way to build the initial configuration for soft-landing sintering simulations: N identical (or distinct) clusters are inserted simultaneously **at rest** (zero velocity), with no substrate, and left to evolve.

The generator synthesizes one source per cluster: a tiny `source_region` at the cluster's center of mass, the cluster as a molecule template, a single target at the central position, and a `fix deposit` with `n_depo 1`, `every 1` and `vx 0.0 0.0`. Because the velocity is zero, each cluster is placed at rest.

## Configuration

A sintering configuration uses a `geometry` block at the top level (instead of a `sources` list). The example below is `configs/4source_tetra_sintering.json`:

```json
{
  "n_source_regions": 4,
  "geometry": {
    "element": "Cu",
    "seeds": [
      "../../lammps_run/seeds/Cu147_Ih.xyz",
      "../../lammps_run/seeds/Cu147_Ih.xyz",
      "../../lammps_run/seeds/Cu147_Ih.xyz",
      "../../lammps_run/seeds/Cu147_Ih.xyz"
    ],
    "geometry": "tetra",
    "displacement": 15.0,
    "max_tries": 50000,
    "type": "Cu",
    "velocity_range": [0.0, 0.0],
    "delta": 0.01,
    "target_size": [0.1, 0.1, 0.1],
    "n_depo": 1,
    "every": 1,
    "random_seed": 42,
    "temperature": 300.0
  }
}
```

`n_source_regions` must equal the number of `seeds`.

### Fields

* `element` / `type` — atomic element label passed to the deposition fixes.
* `seeds` — one XYZ cluster file per source (length must equal `n_source_regions`). Relative paths are resolved against the config file's directory.
* `geometry` — one of the named geometries `line`, `tri`, `square`, `tetra`, `y`, or a user-provided list of displacement directions (see below).
* `displacement` — **mandatory** in geometry mode. It sets the inter-cluster separation. The geometry vectors are **not normalized**: for `line`, consecutive clusters are spaced by exactly `displacement`; for `tri`/`square`/`tetra`, `displacement` is the radial scale from the common center.
* `max_tries` — number of random orientations/placements attempted by the overlap fit (default 50000).
* `overlap_radius` — optional; half the minimum inter-atomic distance used by the overlap check. If omitted, it is derived automatically from the seed's smallest internal pair distance / 2.
* `delta` — source-region half-width (placement precision), default `0.01`.
* `near` — LAMMPS runtime minimum-distance check (`near 0` disables it). Default `2.56`.
* `target_size`, `target_position` — size and center of the (single) target region per source.
* `velocity_range` — deposition velocity range; use `[0.0, 0.0]` to place clusters at rest.
* `n_depo`, `every` — use `1` and `1` so each cluster is inserted once at the first timestep.
* `random_seed`, `temperature` — LAMMPS seed and thermostat target temperature.

### Selecting a geometry

The named geometries are chosen with the `geometry` string. Which geometries are available depends on the number of seeds (`n_source_regions`):

* **2 seeds** — the separation is unique, so the two clusters are placed on a line; the `geometry` value is ignored.
* **3 seeds** — `tri` (equilateral triangle) or `line` (collinear).
* **4 seeds** — `tetra` (tetrahedron), `square`, `line` (collinear), or `y`.

For example, four seeds in a tetrahedral arrangement (as in the example config) use `"geometry": "tetra"`.

Every positions is computed as `pos_i = offset + direction_i * displacement`. For pure directions of unit length (the `tri`/`tetra`/`square` cases) `displacement` is the radial distance from the shared center; for `line`, the directions are multiples of the unit step so the spacing between consecutive clusters equals `displacement`.

### Custom directions

For any number of seeds (for example two, or more than four), or to define a geometry with no dedicated name, provide an explicit list of displacement vectors, one per seed:

```json
"geometry": {
  "geometry": "custom",
  "displacement": 1.0,
  "directions": [
    [15.0, 0.0, 0.0],
    [-15.0, 0.0, 0.0]
  ]
}
```

Each cluster's center of mass is placed at `offset + directions[i] * displacement`. Because the vectors are scaled by `displacement`, set `"displacement": 1.0` if you want the vectors to be used directly as the absolute separations.

## Overlap fit

The generator places each cluster at the derived center of mass and performs a random rotation about its centroid, then checks the minimum inter-cluster distance against `2 * overlap_radius`. It counts how many of the `max_tries` arrangements are fully non-overlapping.

If a fully non-overlapping arrangement is found the fit stops immediately. Otherwise, after `max_tries` attempts it keeps the best arrangement found and prints a warning. Note that a physically impossible or too-compact geometry will spend the full `max_tries` loop, so for such cases lower `max_tries` or increase `displacement`.

## What the generator prints

For a sintering configuration the generator prints the derived center-of-mass positions, the overlap radius, `delta`, `near`, and a semi-automatic type-setup reminder (see section on preparing a run):

```text
Configuration OK: geometry mode with 4 source regions.
OK: found a non-overlapping arrangement after 1 tries (0 overlapping, 1 non-overlapping).
COM positions (box coords):
  cluster 0: (49.000, 49.000, 49.000)
  ...
```

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

## Generating a sintering configuration

A geometry/sintering configuration is generated exactly the same way, by pointing the generator at the sintering config:

```bash
python deposition_generator/generate_depo_input.py \
    --config deposition_generator/configs/4source_tetra_sintering.json \
    --output deposition_generator/generated_inputs/4source_tetra_zerovel.in
```

The generated deposition file has the general form (note that all `molecule` templates are grouped at the top):

```text
molecule molecule_0 <tmpfile>
molecule molecule_1 <tmpfile>
...

region source_box_0 block ...
region target_0_0 block ...
fix deposit_0 added deposit ... region source_box_0 ... target target_0_0 mol molecule_0

...
```

### Ordering of the `molecule` commands

The generator deliberately writes **all** `molecule` commands together at the top of the file, before any region or `fix deposit`. This is required for correctness: `fix deposit ... mol` resolves its molecule template to a pointer inside the LAMMPS `atom->molecules[]` array when the fix is constructed. If additional `molecule` commands are parsed afterwards, that array can grow (realloc), leaving the earlier fixes holding a dangling pointer, which segfaults at the first insertion. Grouping the molecules up front avoids that reallocation.

As a consequence, the generated sintering file references temporary molecule files created during generation. Keep the same `.in` and its companion temp files around (i.e. re-run the generator before launching LAMMPS), and reference the `.in` via `include` in the main input.

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

## Preparing a sintering run

A sintering case (no substrate) is stored, for example, in:

```text
lammps_run/runs/sintering_tetra/
├── lammps_input.in
└── run.sh
```

The main input differs from the growth workflow in a few important ways:

* **No substrate** — there is no `read_dump` of a substrate, and no `coordination` / `attached` group logic. All atoms are the placed clusters.
* **`group added empty`** — the deposit fixes insert into the group `added`, so it must exist (it is created empty).
* **Two atom types** — each deposited molecule atom gets type `base + internal_type`, so with a monoatomic Cu cluster you need `create_box 2`, `mass 1 63.546` and `mass 2 63.546`. Only type 1 is given a `labelmap atom 1 Cu`; **do not** add `labelmap atom 2 Cu` (LAMMPS rejects a duplicate element label). To show both types as Cu in the trajectory use `dump_modify dump1 element Cu Cu`.
* **Thermostat on `all`** — the standard input thermostats the `initial` (substrate) group. Since there is no substrate, use `fix lgv_initial all langevin <T> <T> 0.1 <seed>` so the placed clusters are thermostatted.
* **Tight simulation box** — the box only needs to contain the generated source/target regions plus the cluster radius.
* **Include** — `include` the sintering-generated file (e.g. `4source_tetra_zerovel.in`).

The generator prints the exact number of atom types required and ready-to-paste `labelmap` / `mass` stubs; use those when preparing the main input.

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

See [`INSTALLATION.md`](INSTALLATION.md) for instructions on compiling LAMMPS with the modified source files.

---

# Notes

* Source and target regions must lie inside the LAMMPS simulation box.
* `n_source_regions` must match the number of entries in `sources`.
* For each source, `n_target_regions` must match the number of entries in `targets`.
* Each `source_region` must contain exactly six values.
* Each target `size` and `position` must contain exactly three values.
* Target sizes must be positive.
* The generated deposition input is intended to be included in a complete LAMMPS input file rather than used as a complete simulation by itself.
* Choose the **source/target mode** to grow a nanoparticle by depositing atoms/clusters onto a substrate over time; choose the **geometry/sintering mode** when you want to start from several equilibrated clusters placed at rest and let them sinter.
* In geometry mode a `geometry` block with `seeds` replaces the `sources` list, and `n_source_regions` must equal the number of seeds.
* `overlap_radius` is optional in geometry mode and, if omitted, is derived from the seed's internal minimum pair distance.
* `delta` is the source-region half-width and hence controls the placement precision of each cluster.
* The source/target mode is dependency-free; the geometry mode additionally requires `numpy`, `ase` and `scipy`.

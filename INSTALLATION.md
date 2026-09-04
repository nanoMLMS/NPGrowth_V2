# Installation

This guide describes how to install and compile the modified version of LAMMPS required by NPGrowth V2.

NPGrowth V2 modifies the standard LAMMPS `fix deposit` implementation in order to support deposition toward user-defined target regions.

Because this modification is implemented directly in the LAMMPS source code, LAMMPS must be compiled with the modified `fix_deposit.cpp` and `fix_deposit.h` files provided by this repository.

---

## Requirements

The following software is required:

* Git
* CMake
* a C++ compiler with C++17 support
* MPI
* Python 3

You can check whether these tools are available with:

```bash
git --version
cmake --version
mpirun --version
python3 --version
```

The exact procedure for installing these dependencies depends on the operating system and computing environment.

---

## LAMMPS version used

NPGrowth V2 was developed and tested with:

```text
LAMMPS version: 22 Jul 2025 - Update 4
LAMMPS branch: stable
LAMMPS commit: 611ca3b7f8525ba802373b04f9e3d632d515f7e3
```

The local development build reports:

```text
stable_22Jul2025_update4-modified
```

because the standard LAMMPS source files:

```text
src/fix_deposit.cpp
src/fix_deposit.h
```

were modified for the NPGrowth deposition workflow.

Using the same LAMMPS commit is strongly recommended to ensure reproducibility and compatibility with the modified `fix deposit` implementation.

Other LAMMPS versions may also work, but they have not necessarily been tested. Changes in the internal LAMMPS API may require adapting the modified source files.

---

## 1. Clone NPGrowth V2

Clone the NPGrowth V2 repository:

```bash
git clone <NPGROWTH_REPOSITORY_URL>
cd NPGrowth_V2
```

The modified LAMMPS deposition files are located in:

```text
lammps_source_modified/
└── lammps_deposition_modified/
    ├── fix_deposit.cpp
    └── fix_deposit.h
```

These files replace the corresponding files in the standard LAMMPS source tree.

---

## 2. Clone LAMMPS

From the directory where you want to install LAMMPS:

```bash
git clone https://github.com/lammps/lammps.git
cd lammps
```

Check out the exact commit used by NPGrowth V2:

```bash
git checkout 611ca3b7f8525ba802373b04f9e3d632d515f7e3
```

You can verify the current commit with:

```bash
git rev-parse HEAD
```

The output should be:

```text
611ca3b7f8525ba802373b04f9e3d632d515f7e3
```

At this point, the LAMMPS source tree should contain:

```text
lammps/
├── cmake/
├── src/
├── examples/
├── potentials/
└── ...
```

The standard deposition implementation is located in:

```text
lammps/src/fix_deposit.cpp
lammps/src/fix_deposit.h
```

---

## 3. Replace `fix_deposit`

Copy the modified files from NPGrowth V2 into the LAMMPS `src` directory.

Assuming the two repositories are located side by side:

```text
parent_directory/
├── NPGrowth_V2/
└── lammps/
```

from inside the `lammps` directory run:

```bash
cp ../NPGrowth_V2/lammps_source_modified/lammps_deposition_modified/fix_deposit.cpp \
   src/fix_deposit.cpp

cp ../NPGrowth_V2/lammps_source_modified/lammps_deposition_modified/fix_deposit.h \
   src/fix_deposit.h
```

You can verify that only these two LAMMPS source files have been modified with:

```bash
git status
```

The expected result is:

```text
modified: src/fix_deposit.cpp
modified: src/fix_deposit.h
```

This confirms that the standard `fix deposit` implementation has been replaced by the NPGrowth V2 version.

---

## 4. Create the build directory

From the LAMMPS source directory:

```bash
mkdir -p build
cd build
```

LAMMPS is built out-of-source, so compiled files are kept separate from the source tree.

---

## 5. Configure LAMMPS with CMake

Configure the build with:

```bash
cmake \
    -D BUILD_MPI=yes \
    -D PKG_MANYBODY=yes \
    ../cmake
```

The NPGrowth V2 development build uses:

```text
BUILD_MPI=yes
PKG_MANYBODY=yes
```

`BUILD_MPI=yes` enables MPI support.

`PKG_MANYBODY=yes` enables the LAMMPS MANYBODY package, which includes interaction styles used by many many-body potentials.

The reference NPGrowth V2 build includes the following LAMMPS packages:

```text
KSPACE
MANYBODY
MOLECULE
```

Depending on the local LAMMPS configuration, some packages may already be enabled automatically.

If CMake completes without errors, LAMMPS is ready to be compiled.

---

## 6. Compile LAMMPS

From inside:

```text
lammps/build/
```

compile with:

```bash
cmake --build . -j 4
```

The value after `-j` specifies the number of parallel compilation jobs.

For example:

```bash
cmake --build . -j 8
```

uses eight parallel compilation jobs.

After a successful build, the executable should be available at:

```text
lammps/build/lmp
```

---

## 7. Verify the LAMMPS installation

Run:

```bash
./lmp -help
```

The beginning of the output should report:

```text
Large-scale Atomic/Molecular Massively Parallel Simulator - 22 Jul 2025 - Update 4
```

The build should also include the `deposit` fix style.

In the `Fix styles` section of:

```bash
./lmp -help
```

you should find:

```text
deposit
```

At this point, the modified LAMMPS executable is ready to be used with NPGrowth V2.

---

## 8. Generate a deposition input

Return to the NPGrowth V2 repository:

```bash
cd /path/to/NPGrowth_V2
```

Generate a deposition input from a JSON configuration:

```bash
python3 deposition_generator/generate_depo_input.py \
    --config deposition_generator/configs/1source_1targets.json \
    --output deposition_generator/generated_inputs/1source_1targets.in
```

If the configuration is valid, the generator prints messages similar to:

```text
Configuration OK: 1 source regions found.
LAMMPS deposition input written to: deposition_generator/generated_inputs/1source_1targets.in
```

The generated file contains:

* source regions;
* target regions;
* target unions when multiple targets are used;
* the modified `fix deposit` commands.

---

## 9. Prepare a simulation

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

The main `lammps_input.in` defines the complete LAMMPS simulation and includes the deposition input generated in the previous step.

Before running the simulation, verify the paths used by:

```text
read_dump
pair_coeff
include
```

These paths must point to:

* the nanoparticle seed;
* the interatomic potential;
* the generated deposition input.

Nanoparticle structures are stored in:

```text
lammps_run/seeds/
```

Interatomic potentials are stored in:

```text
lammps_run/potentials/
```

---

## 10. Configure the run script

The example run script launches LAMMPS through MPI.

The general syntax is:

```bash
mpirun -np <N_PROCESSES> /path/to/lammps/build/lmp -in lammps_input.in
```

For example, using four MPI processes:

```bash
mpirun -np 4 /path/to/lammps/build/lmp -in lammps_input.in
```

Update the path to `lmp` according to the location of your local LAMMPS installation.

---

## 11. Run the simulation

Enter the simulation directory:

```bash
cd lammps_run/runs/1source_1target
```

Run:

```bash
bash run.sh
```

Alternatively, execute LAMMPS directly:

```bash
mpirun -np 4 /path/to/lammps/build/lmp -in lammps_input.in
```

LAMMPS reads the main `lammps_input.in`, which in turn includes the deposition instructions generated by `generate_depo_input.py`.

---

## Installation workflow

The complete installation and execution workflow is:

```text
Clone NPGrowth V2
        │
        ▼
Clone LAMMPS
        │
        ▼
Checkout commit
611ca3b7f8525ba802373b04f9e3d632d515f7e3
        │
        ▼
Replace fix_deposit.cpp
and fix_deposit.h
        │
        ▼
Configure LAMMPS with CMake
        │
        ▼
Compile LAMMPS
        │
        ▼
Generate deposition input
        │
        ▼
Prepare lammps_input.in
        │
        ▼
Run with MPI
```

---

## Updating LAMMPS

NPGrowth V2 modifies internal LAMMPS source files.

The modified `fix_deposit.cpp` and `fix_deposit.h` files use internal LAMMPS classes and APIs. Therefore, compatibility with newer or older LAMMPS versions is not guaranteed.

For reproducible simulations, use:

```text
LAMMPS 22 Jul 2025 - Update 4
commit 611ca3b7f8525ba802373b04f9e3d632d515f7e3
```

If a different LAMMPS version is used, the modified deposition implementation should be tested carefully before production simulations.

---

## Troubleshooting

### `mpirun: command not found`

MPI is either not installed or is not available in the current shell environment.

Install or load an MPI implementation before configuring LAMMPS with:

```text
BUILD_MPI=yes
```

---

### CMake cannot find MPI

Make sure that the MPI installation is available in the current environment before running CMake.

On HPC systems, MPI may need to be loaded through the module system before configuring LAMMPS.

---

### `lmp` was compiled before replacing `fix_deposit`

If `fix_deposit.cpp` and `fix_deposit.h` were replaced after LAMMPS had already been compiled, rebuild the executable:

```bash
cd lammps/build
cmake --build . -j 4
```

---

### LAMMPS does not recognize the expected target-region deposition behavior

Check that the modified:

```text
fix_deposit.cpp
fix_deposit.h
```

were copied into:

```text
lammps/src/
```

before compilation.

You can check with:

```bash
git status
```

The two files should appear as modified.

---

### LAMMPS cannot open the seed, potential, or generated input

Check the paths in `lammps_input.in`, especially those used by:

```text
read_dump
pair_coeff
include
```

---

### Python reports a configuration error

The deposition generator validates the JSON configuration before creating the LAMMPS input.

Check that:

* `n_source_regions` matches the number of entries in `sources`;
* each `source_region` contains exactly six values;
* `n_target_regions` matches the number of entries in `targets`;
* each target contains `size` and `position`;
* `size` contains three positive values;
* `position` contains exactly three values;
* all required deposition parameters are present.

---

## Reference development environment

The reference development build was compiled on macOS ARM64 with:

```text
OS: Darwin arm64
Compiler: Apple Clang
C++ standard: C++17
MPI: Open MPI
```

These exact versions are not required, but they document the environment in which the current implementation was tested.

---

## LAMMPS license

LAMMPS is developed by Sandia National Laboratories and the LAMMPS development team and is distributed under the GNU General Public License.

The modified `fix_deposit.cpp` and `fix_deposit.h` files provided by NPGrowth V2 are derived from the corresponding LAMMPS source files and retain the original LAMMPS copyright and license notices.

Refer to the official LAMMPS project for complete licensing, documentation, and citation information.


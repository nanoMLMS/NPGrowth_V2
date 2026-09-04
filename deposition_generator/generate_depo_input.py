import json
import tempfile


def format_block(region):
    return "{:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f}".format(*region)


def load_json(path):
    with open(path, "r") as f:
        return json.load(f)


def write_region(output_file, name, region):
    output_file.write(
        f"region {name} block {format_block(region)}\n"
    )


def region_from_target(target):
    """
    Build a LAMMPS block region from:
        size     = [sx, sy, sz]
        position = [x, y, z]

    position is the CENTER of the cuboid.
    """

    sx, sy, sz = target["size"]
    x, y, z = target["position"]

    dx = sx / 2.0
    dy = sy / 2.0
    dz = sz / 2.0

    return [
        x - dx,
        x + dx,
        y - dy,
        y + dy,
        z - dz,
        z + dz
    ]


def parse_xyz(path):
    with open(path, "r") as f:
        natoms = int(f.readline())

        # Comment line
        _ = f.readline()

        atoms = []

        for line in f:
            parts = line.strip().split()

            if len(parts) >= 4:
                atom_type, x, y, z = parts[:4]

                atoms.append(
                    (
                        atom_type,
                        float(x),
                        float(y),
                        float(z)
                    )
                )

    return natoms, atoms


def write_molecule(molecule_file):
    """
    Convert an XYZ molecule file to the temporary LAMMPS molecule format.
    """

    lammps_file = tempfile.NamedTemporaryFile(
        "w",
        delete=False
    )

    natoms, atoms = parse_xyz(molecule_file)

    with open(lammps_file.name, "w") as f:
        f.write(f"# Nanoparticle {molecule_file}\n")
        f.write(f"{natoms} atoms\n\n")

        f.write("Coords\n\n")

        for i, (_, x, y, z) in enumerate(atoms, start=1):
            f.write(f"{i} {x} {y} {z}\n")

        f.write("\nTypes\n\n")

        for i, (atom_type, _, _, _) in enumerate(atoms, start=1):
            f.write(f"{i} {atom_type}\n")

        f.write("\n")

    return lammps_file.name


def validate_input(data):
    """
    Validate the JSON before generating the LAMMPS input.
    """

    # ==========================================================
    # GLOBAL SOURCE CHECK
    # ==========================================================

    if "n_source_regions" not in data:
        raise ValueError(
            "Missing 'n_source_regions' in the configuration file."
        )

    if "sources" not in data:
        raise ValueError(
            "Missing 'sources' in the configuration file."
        )

    n_source_regions = data["n_source_regions"]
    sources = data["sources"]

    if not isinstance(n_source_regions, int):
        raise ValueError(
            "'n_source_regions' must be an integer."
        )

    if n_source_regions <= 0:
        raise ValueError(
            "'n_source_regions' must be greater than zero."
        )

    if not isinstance(sources, list):
        raise ValueError(
            "'sources' must be a list."
        )

    if n_source_regions != len(sources):
        raise ValueError(
            "\nSOURCE REGION ERROR\n"
            f"Declared n_source_regions = {n_source_regions}\n"
            f"But {len(sources)} source regions were provided "
            "inside 'sources'."
        )

    # ==========================================================
    # CHECK EACH SOURCE
    # ==========================================================

    for i, source in enumerate(sources):

        # ------------------------------------------------------
        # Source region
        # ------------------------------------------------------

        if "source_region" not in source:
            raise ValueError(
                f"Source {i}: missing 'source_region'."
            )

        source_region = source["source_region"]

        if not isinstance(source_region, list):
            raise ValueError(
                f"Source {i}: 'source_region' must be a list."
            )

        if len(source_region) != 6:
            raise ValueError(
                f"Source {i}: 'source_region' must contain "
                "exactly 6 values:\n"
                "[xlo, xhi, ylo, yhi, zlo, zhi]"
            )

        xlo, xhi, ylo, yhi, zlo, zhi = source_region

        if xlo >= xhi:
            raise ValueError(
                f"Source {i}: xlo must be smaller than xhi."
            )

        if ylo >= yhi:
            raise ValueError(
                f"Source {i}: ylo must be smaller than yhi."
            )

        if zlo >= zhi:
            raise ValueError(
                f"Source {i}: zlo must be smaller than zhi."
            )

        # ------------------------------------------------------
        # Number of targets
        # ------------------------------------------------------

        if "n_target_regions" not in source:
            raise ValueError(
                f"Source {i}: missing 'n_target_regions'."
            )

        if "targets" not in source:
            raise ValueError(
                f"Source {i}: missing 'targets'."
            )

        n_target_regions = source["n_target_regions"]
        targets = source["targets"]

        if not isinstance(n_target_regions, int):
            raise ValueError(
                f"Source {i}: 'n_target_regions' must be an integer."
            )

        if n_target_regions <= 0:
            raise ValueError(
                f"Source {i}: 'n_target_regions' "
                "must be greater than zero."
            )

        if not isinstance(targets, list):
            raise ValueError(
                f"Source {i}: 'targets' must be a list."
            )

        if n_target_regions != len(targets):
            raise ValueError(
                "\nTARGET REGION ERROR\n"
                f"Source {i}:\n"
                f"Declared n_target_regions = {n_target_regions}\n"
                f"But {len(targets)} targets were provided."
            )

        # ------------------------------------------------------
        # Check each target
        # ------------------------------------------------------

        for j, target in enumerate(targets):

            if "size" not in target:
                raise ValueError(
                    f"Source {i}, target {j}: missing 'size'."
                )

            if "position" not in target:
                raise ValueError(
                    f"Source {i}, target {j}: missing 'position'."
                )

#            if "rotation" not in target:
#                raise ValueError(
#                    f"Source {i}, target {j}: missing 'rotation'."
#                )

            size = target["size"]
            position = target["position"]
#            rotation = target["rotation"]

            if len(size) != 3:
                raise ValueError(
                    f"Source {i}, target {j}: "
                    "'size' must contain 3 values [sx, sy, sz]."
                )

            if len(position) != 3:
                raise ValueError(
                    f"Source {i}, target {j}: "
                    "'position' must contain 3 values [x, y, z]."
                )

#            if len(rotation) != 3:
#                raise ValueError(
#                    f"Source {i}, target {j}: "
#                    "'rotation' must contain 3 values [rx, ry, rz]."
#                )

            if any(s <= 0 for s in size):
                raise ValueError(
                    f"Source {i}, target {j}: "
                    "all values in 'size' must be greater than zero."
                )

        # ------------------------------------------------------
        # Other required source parameters
        # ------------------------------------------------------

        required_keys = [
            "type",
            "velocity_range",
            "near",
            "n_depo",
            "every",
            "random_seed_LAMMPS"
        ]

        for key in required_keys:
            if key not in source:
                raise ValueError(
                    f"Source {i}: missing '{key}'."
                )

        if len(source["velocity_range"]) != 2:
            raise ValueError(
                f"Source {i}: 'velocity_range' must contain exactly 2 values."
            )


def generate_script(config_path, output_path):

    # ==========================================================
    # LOAD + VALIDATE
    # ==========================================================

    data = load_json(config_path)

    validate_input(data)

    configs = data["sources"]

    print(
        f"Configuration OK: "
        f"{len(configs)} source regions found."
    )

    # ==========================================================
    # WRITE LAMMPS INPUT
    # ==========================================================

    with open(output_path, "w") as output_file:

        for i, config in enumerate(configs):

            # ==================================================
            # SOURCE REGION
            # ==================================================

            source_region_name = f"source_box_{i}"

            write_region(
                output_file,
                source_region_name,
                config["source_region"]
            )

            # ==================================================
            # OPTIONAL MOLECULE
            # ==================================================

            molecule_file = config.get("molecule_xyz")

            molecule_name = f"molecule_{i}"

            mol_clause = ""

            if molecule_file:

                lammps_molecule_file = write_molecule(
                    molecule_file
                )

                output_file.write(
                    f"molecule {molecule_name} "
                    f"{lammps_molecule_file}\n"
                )

                mol_clause = f"mol {molecule_name}"

            # ==================================================
            # TARGET REGIONS
            # ==================================================

            target_region_names = []

            for j, target in enumerate(config["targets"]):

                target_region_name = f"target_{i}_{j}"

                region = region_from_target(target)

                write_region(
                    output_file,
                    target_region_name,
                    region
                )

                target_region_names.append(
                    target_region_name
                )

            # ==================================================
            # TARGET UNION
            # ==================================================

            if len(target_region_names) > 1:

                target_union_name = f"target_box_{i}"

                output_file.write(
                    f"region {target_union_name} union "
                    f"{len(target_region_names)} "
                    f"{' '.join(target_region_names)}\n"
                )

            else:

                target_union_name = target_region_names[0]

            # ==================================================
            # FIX DEPOSIT
            # ==================================================

            output_file.write(
                f"fix deposit_{i} added deposit "
                f"{config['n_depo']} "
                f"{config['type']} "
                f"{config['every']} "
                f"{config['random_seed_LAMMPS']} "
                f"region {source_region_name} "
                f"near {config['near']} "
                f"vx {config['velocity_range'][0]} {config['velocity_range'][1]} "
                f"target {target_union_name}"
            )

            if mol_clause:
                output_file.write(
                    f" {mol_clause}"
                )

            output_file.write("\n")

            # Blank line between sources
            output_file.write("\n")


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        required=True,
        help="Path to deposition config JSON"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Path to output LAMMPS file"
    )

    args = parser.parse_args()

    generate_script(
        args.config,
        args.output
    )

    print(
        f"LAMMPS deposition input written to: "
        f"{args.output}"
    )

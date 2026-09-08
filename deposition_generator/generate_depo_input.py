import json
import os
import math
import tempfile
import numpy as np
from scipy.spatial.distance import cdist, pdist
from ase.io import read as ase_read


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


# ==============================================================
# Geometry-driven source placement (sintering-style configs)
#
# When the top-level JSON contains a "geometry" block, the source
# regions are not written by hand: this code places N equilibrated
# clusters (seeds) at the COM positions of a chosen geometry at a
# controlled separation, reusing the same logic as the V1
# make-sintering-config.py (self-contained copy, no cross-repo
# import).  Each cluster becomes one source deposition cannon with
# zero insertion velocity (vx 0 0), i.e. the cluster is placed at
# rest at its source region.
#
# Geometry semantics (identical to V1): direction vectors are used
# UN-normalised (direction * displacement).  For 'line' this makes
# `displacement` the consecutive COM separation; for the other
# geometries it sets the radial scale (actual pair separation is
# displacement * the un-normalised pair distance).
# ==============================================================

COS30 = math.sqrt(3.0) / 2.0

GEOMETRIES_3 = {
    "tri":  np.array([[0., 1., 0.], [-COS30, -0.5, 0.], [COS30, -0.5, 0.]]),
    "line": np.array([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]]),
}

GEOMETRIES_4 = {
    "line":   np.array([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.], [3., 0., 0.]]),
    "square": np.array([[1., 1., 0.], [1., -1., 0.], [-1., -1., 0.], [-1., 1., 0.]]),
    "tetra":  np.array([[1., 1., 1.], [1., -1., -1.], [-1., 1., -1.], [-1., -1., 1.]]),
    "y":      np.array([[-0.5, COS30, 0.], [0.5, COS30, 0.], [0., 0., 0.], [0., -1., 0.]]),
}

GEOMETRIES = {3: GEOMETRIES_3, 4: GEOMETRIES_4}


def random_rotation_matrix():
    """Uniformly random 3D rotation via QR of a random normal matrix."""
    A = np.random.randn(3, 3)
    Q, R = np.linalg.qr(A)
    Q *= np.sign(np.linalg.det(Q))
    return Q


def rotate_cluster(atoms, R):
    """Rotate an ASE Atoms object in-place about its centroid (COM pinned)."""
    pos = atoms.get_positions()
    centroid = pos.mean(axis=0)
    atoms.set_positions((R @ (pos - centroid).T).T + centroid)
    return atoms


def randomly_rotate_all(clusters):
    """Apply a different random rotation to each cluster in-place."""
    for c in clusters:
        rotate_cluster(c, random_rotation_matrix())
    return clusters


def overlap_pair_count(clusters, radius):
    """Number of atom pairs across different clusters closer than 2*radius."""
    n = 0
    threshold = 2 * radius
    for i in range(len(clusters)):
        for j in range(i + 1, len(clusters)):
            a = clusters[i].get_positions()
            b = clusters[j].get_positions()
            n += int(np.sum(cdist(a, b) < threshold))
    return n


def resolve_directions_geometry(n, geometry, directions):
    """Return the (n, 3) un-normalised direction array for the geometry."""
    if directions is not None:
        dirs = np.asarray(directions, dtype=float)
        if dirs.shape != (n, 3):
            raise ValueError(
                f"directions must be ({n}, 3), got {dirs.shape}"
            )
        return dirs

    if n == 2:
        return np.array([[0., 0., 0.], [1., 0., 0.]])

    if n not in (3, 4):
        raise ValueError(f"n_clusters must be 2, 3, or 4, got {n}")

    geom_dict = GEOMETRIES.get(n, {})
    if geometry not in geom_dict:
        raise ValueError(
            f"Unknown geometry '{geometry}' for {n} clusters. "
            f"Available: {sorted(geom_dict.keys())}"
        )
    return geom_dict[geometry].copy()


def infer_overlap_radius(clusters, explicit):
    """Auto-compute overlap_radius from the seed internal min bond distance."""
    if explicit is not None:
        return float(explicit)
    mins = []
    for c in clusters:
        pos = c.get_positions()
        d = pdist(pos)
        mins.append(float(np.min(d)))
    if not mins or min(mins) <= 0.0:
        return 1.28
    return min(mins) / 2.0


def load_clusters(seed_paths):
    """Load N seed clusters as ASE Atoms."""
    return [ase_read(path) for path in seed_paths]


def compute_source_regions(clusters, cfg):
    """
    Compute the N cluster COM positions (absolute box coords, positive):
        pos_i = offset + direction_i * displacement
    Also pin each cluster's COM to its box position (rotation about the
    centroid keeps it fixed).
    Returns (positions, offset).
    """
    displacement = cfg["displacement"]
    directions = cfg["_directions"]

    radii = []
    for c in clusters:
        com = np.mean(c.positions, axis=0)
        radii.append(np.max(np.linalg.norm(c.positions - com, axis=1)))
    rmax = max(radii) if radii else 0.0

    raw = directions * displacement

    offset = cfg.get("offset")
    if offset is None:
        lo = np.min(raw, axis=0) - rmax
        margin = rmax + cfg["box_margin"]
        offset = (margin - lo).astype(int)
    offset = np.asarray(offset, dtype=float)

    positions = offset + raw

    # pin each cluster COM to its box position
    for c, pos in zip(clusters, positions):
        c.positions += (pos - np.mean(c.positions, axis=0))

    return positions, offset


def fit_no_overlap(clusters, cfg):
    """
    Randomly rotate all clusters until no inter-cluster overlap.
    Tally tries / non-overlapping successes, track best arrangement,
    warn + proceed on exhaustion.
    Returns (best_clusters, summary).
    """
    overlap_radius = cfg["overlap_radius"]
    max_tries = cfg["max_tries"]
    np.random.seed(int(cfg["random_seed"]))

    best_clusters = [c.copy() for c in clusters]
    best_violations = None
    n_overlap = n_ok = 0

    for _ in range(max_tries):
        randomly_rotate_all(clusters)
        violations = overlap_pair_count(clusters, overlap_radius)
        if violations > 0:
            n_overlap += 1
            if best_violations is None or violations < best_violations:
                best_violations = violations
                best_clusters = [c.copy() for c in clusters]
        else:
            n_ok += 1
            best_clusters = [c.copy() for c in clusters]
            return best_clusters, {
                "found": True, "tries": n_overlap + n_ok,
                "n_overlapping": n_overlap, "n_nonoverlapping": n_ok,
                "best_violations": 0,
            }

    return best_clusters, {
        "found": False, "tries": max_tries,
        "n_overlapping": n_overlap, "n_nonoverlapping": n_ok,
        "best_violations": best_violations,
    }


def max_molecule_type(clusters):
    """
    Max internal atom type index across the seed templates.

    In LAMMPS molecule templates the per-atom types are the species
    indices (1..n_species of the template), not the element atomic
    numbers.  For a monoatomic cluster this is 1.
    """
    m = 0
    for c in clusters:
        symbols = c.get_chemical_symbols()
        n_species = len(set(symbols))
        m = max(m, n_species)
    return m


def expand_geometry(data):
    """
    Turn a top-level 'geometry' block into (sources, meta) where
    * sources is a list of source dicts (one per cluster) ready for the
      standard generation path,
    * meta carries geometry summary, type/thermostat info printed to user.
    """
    g = data["geometry"]
    n = data["n_source_regions"]

    seed_paths = g["seeds"]
    if len(seed_paths) != n:
        raise ValueError(
            f"geometry: len(seeds)={len(seed_paths)} != n_source_regions={n}"
        )

    config_dir = os.path.dirname(os.path.abspath(data.get("_config_path", ".")))
    seed_paths = [
        os.path.abspath(os.path.join(config_dir, s)) if not os.path.isabs(s) else s
        for s in seed_paths
    ]

    clusters = load_clusters(seed_paths)

    # ---- displacement is mandatory in geometry mode ----
    if "displacement" not in g:
        raise ValueError(
            "geometry: 'displacement' is required in geometry mode. "
            "It sets the inter-cluster separation."
        )

    # ---- defaults for the geometry block ----
    g.setdefault("box_margin", 5.0)
    g.setdefault("max_tries", 50000)
    g.setdefault("delta", 0.01)
    g.setdefault("target_size", [0.1, 0.1, 0.1])
    g.setdefault("near", 2.56)
    g.setdefault("n_depo", 1)
    g.setdefault("every", 1)
    g.setdefault("random_seed", 42)
    g.setdefault("velocity_range", [0.0, 0.0])
    g.setdefault("temperature", 300.0)
    g.setdefault("type", g.get("element"))

    directions = g.get("directions")
    g["_directions"] = resolve_directions_geometry(n, g.get("geometry"), directions)

    # ---- overlap radius: explicit or auto-computed from seeds ----
    g["overlap_radius"] = infer_overlap_radius(
        clusters, g.get("overlap_radius")
    )

    positions, offset = compute_source_regions(clusters, g)

    best_clusters, summary = fit_no_overlap(clusters, g)

    # ---- target position: default mean COM (central "seed" aim point) ----
    tpos = g.get("target_position")
    if tpos is None:
        tpos = np.mean(positions, axis=0)
    tpos = [float(v) for v in tpos]

    delta = g["delta"]
    sources = []
    for i, pos in enumerate(positions):
        sources.append({
            "source_region": [
                float(pos[0] - delta), float(pos[0] + delta),
                float(pos[1] - delta), float(pos[1] + delta),
                float(pos[2] - delta), float(pos[2] + delta),
            ],
            "molecule_xyz": seed_paths[i],
            "n_target_regions": 1,
            "targets": [
                {"size": list(g["target_size"]), "position": tpos}
            ],
            "type": g["type"],
            "velocity_range": [float(v) for v in g["velocity_range"]],
            "near": g["near"],
            "n_depo": g["n_depo"],
            "every": g["every"],
            "random_seed_LAMMPS": int(g["random_seed"]) + i,
        })

    base_type_index = g.get("type_index", 1)
    max_type = base_type_index + max_molecule_type(best_clusters)

    meta = {
        "summary": summary,
        "positions": positions,
        "offset": offset,
        "delta": delta,
        "overlap_radius": g["overlap_radius"],
        "near": g["near"],
        "temperature": g["temperature"],
        "thermo_seed": int(g["random_seed"]),
        "base_type_index": base_type_index,
        "max_type": max_type,
        "element": g.get("element", g.get("type")),
        "seed_paths": seed_paths,
    }
    return sources, meta


def print_geometry_meta(meta):
    """Print the geometry summary, type stubs and thermostat reminder."""
    bar = "=" * 70
    print(bar)
    s = meta["summary"]
    if s["found"]:
        print(
            f"OK: found a non-overlapping arrangement after {s['tries']} tries "
            f"({s['n_overlapping']} overlapping, "
            f"{s['n_nonoverlapping']} non-overlapping)."
        )
    else:
        print(
            f"WARNING: no fully non-overlapping arrangement in {s['tries']} tries "
            f"({s['n_overlapping']} overlapping, "
            f"{s['n_nonoverlapping']} non-overlapping). "
            f"Writing best arrangement (min overlap pairs = "
            f"{s['best_violations']}). Reduce overlap_radius / increase "
            f"displacement / max_tries."
        )

    print("COM positions (box coords):")
    for i, pos in enumerate(meta["positions"]):
        print(f"  cluster {i}: ({pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f})")
    print(f"offset: {meta['offset']}")
    print(f"overlap_radius (used by fit): {meta['overlap_radius']}")
    print(f"source half-width (delta): {meta['delta']}")
    print(
        f"near (LAMMPS runtime min distance): {meta['near']} "
        f"(format: rsq < near^2 rejects an insertion)"
    )
    print(bar)

    # molecule template has M internal types; base type index -> total required
    print("SEMI-AUTOMATIC TYPE SETUP (cut & paste into the main input):")
    print(
        f"  # geometry mode target: base type '{meta['element']}' index "
        f"{meta['base_type_index']}, molecule template max type reached = "
        f"{meta['max_type']}"
    )
    print("  # declare this many atom types in the box and labelmap them:")
    for t in range(meta["base_type_index"], meta["max_type"] + 1):
        print(f"  labelmap atom {t} {meta['element']}")
        print(f"  mass {t} 63.546")
    print("  # ensure pair_coeff * * <potential> covers all these types")

    print(bar)
    print("THERMOSTAT REMINDER (geometry mode has NO substrate):")
    print("  # all atoms are the cannon-placed clusters, so thermostats must act on ALL atoms:")
    print(f"  fix my_nve all nve")
    print(
        f"  fix lgv_initial all langevin {meta['temperature']} "
        f"{meta['temperature']} 0.1 {meta['thermo_seed']}"
    )
    print("  # (do NOT use the stock 'initial'-group thermostat here)")
    print(bar)


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

    # record config file location so relative seed paths resolve correctly
    data["_config_path"] = config_path

    meta = None

    if "geometry" in data:
        # geometry-driven mode: synthesize the sources automatically
        configs, meta = expand_geometry(data)
        print(
            f"Configuration OK: geometry mode with "
            f"{len(configs)} source regions."
        )
        print_geometry_meta(meta)
    else:
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

        # ----------------------------------------------------------
        # All molecule templates go first, contiguously, BEFORE any
        # region/fix that will be emitted below.
        #
        # Why: fix deposit resolves its molecule template to a
        # pointer into the atom->molecules[] array at construction
        # time. If more molecule commands are parsed afterwards the
        # array can grow (realloc), leaving those saved pointers
        # dangling -> segfault at the first insertion. Defining all
        # molecules up front avoids that reallocation entirely.
        # ----------------------------------------------------------

        mol_clauses = {}

        for i, config in enumerate(configs):

            molecule_file = config.get("molecule_xyz")

            molecule_name = f"molecule_{i}"

            if molecule_file:

                lammps_molecule_file = write_molecule(
                    molecule_file
                )

                output_file.write(
                    f"molecule {molecule_name} "
                    f"{lammps_molecule_file}\n"
                )

                mol_clauses[i] = f"mol {molecule_name}"

        if mol_clauses:
            output_file.write("\n")

        for i, config in enumerate(configs):

            mol_clause = mol_clauses.get(i, "")

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

    try:
        generate_script(
            args.config,
            args.output
        )
    except ValueError as e:
        print(f"Configuration error: {e}")
        raise SystemExit(1)

    print(
        f"LAMMPS deposition input written to: "
        f"{args.output}"
    )

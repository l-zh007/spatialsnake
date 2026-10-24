import sys
import types
from importlib import util
from importlib import metadata


def apply_numpy_compat():
    import numpy as np

    # PySCENIC 0.12.1 still references NumPy aliases removed in NumPy 2.x.
    compat_aliases = {
        "object": object,
        "float": float,
        "int": int,
        "bool": bool,
        "complex": complex,
        "str": str,
    }
    for alias, target in compat_aliases.items():
        if alias not in np.__dict__:
            setattr(np, alias, target)  # type: ignore[attr-defined]


def apply_pkg_resources_compat():
    if util.find_spec("pkg_resources") is not None:
        return

    compat_module = types.ModuleType("pkg_resources")

    class DistributionNotFound(Exception):
        pass

    def get_distribution(dist_name):
        try:
            version = metadata.version(dist_name)
        except metadata.PackageNotFoundError as exc:
            raise DistributionNotFound(str(exc)) from exc

        return types.SimpleNamespace(version=version)

    compat_module.DistributionNotFound = DistributionNotFound
    compat_module.get_distribution = get_distribution
    sys.modules["pkg_resources"] = compat_module


def apply_scipy_compat():
    from scipy.sparse import spmatrix

    # Arboreto still uses the sparse .A alias removed in SciPy 1.14.
    if not hasattr(spmatrix, "A"):
        spmatrix.A = property(lambda matrix: matrix.toarray())


def main():
    if len(sys.argv) < 2:
        raise SystemExit("Usage: pyscenic_numpy_compat.py [arboreto|pyscenic] ...")

    command = sys.argv[1]
    argv = sys.argv[2:]

    if command == "arboreto":
        from pyscenic.cli.arboreto_with_multiprocessing import main as arboreto_main

        sys.argv = [sys.argv[0]] + argv
        arboreto_main()
        if "--output" in argv:
            output = argv[argv.index("--output") + 1]
            with open(output) as handle:
                next(handle, None)
                if not any(line.strip() for line in handle):
                    raise ValueError("GRN inference produced no regulatory links; check the inference warnings above.")
        return

    if command == "pyscenic":
        from pyscenic.cli.pyscenic import main as pyscenic_main

        pyscenic_main(argv)
        return

    raise SystemExit(f"Unsupported command: {command}")


# Also initialize compatibility when multiprocessing uses spawn.
apply_numpy_compat()
apply_pkg_resources_compat()
apply_scipy_compat()


if __name__ == "__main__":
    main()

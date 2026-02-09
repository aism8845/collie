#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

def load_dofmap_json(path: Path):
    """
    Supports MOOSE DOFMapOutput JSON schema:
      {"ndof": <int>, "vars": [{"name":..., "subdomains":[{"id":..., "kernels":[...], "dofs":[...]}]}]}
    Also keeps a fallback for older/other record-list schemas if present.
    """
    data = json.loads(path.read_text())

    # --- Schema A: DOFMapOutput (variable -> dof list) ---
    if isinstance(data, dict) and "vars" in data and "ndof" in data:
        mapping = {}
        for v in data.get("vars", []):
            vname = v.get("name", None)
            for sd in v.get("subdomains", []):
                sd_id = sd.get("id", None)
                kernels = []
                for k in sd.get("kernels", []):
                    if isinstance(k, dict) and "name" in k:
                        kernels.append(k["name"])
                for dof in sd.get("dofs", []):
                    try:
                        did = int(dof)
                    except Exception:
                        continue
                    mapping[did] = {
                        "variable": vname,
                        "subdomain": sd_id,
                        "kernels": kernels,
                    }
        if not mapping:
            raise RuntimeError(f"Parsed DOFMapOutput JSON but found no DOFs in {path}")
        return mapping

    # --- Schema B (fallback): list/dict of per-DOF records (best-effort) ---
    def _find_records(obj):
        if isinstance(obj, list) and obj and isinstance(obj[0], dict):
            return obj
        if isinstance(obj, dict):
            for k in ["dofmap", "dof_map", "dofs", "dof_records", "records", "data"]:
                if k in obj:
                    rec = _find_records(obj[k])
                    if rec:
                        return rec
            for v in obj.values():
                rec = _find_records(v)
                if rec:
                    return rec
        return None

    recs = _find_records(data)
    if not recs:
        raise RuntimeError(f"Unrecognized dofmap JSON schema in {path}")

    mapping = {}
    for r in recs:
        dof = None
        for k in ["dof", "dof_id", "dofID", "id", "row", "global_dof"]:
            if k in r:
                dof = r[k]
                break
        if dof is None:
            continue
        try:
            did = int(dof)
        except Exception:
            continue
        mapping[did] = {
            "variable": r.get("variable") or r.get("var") or r.get("name"),
            "subdomain": r.get("subdomain"),
            "kernels": None,
        }
    if not mapping:
        raise RuntimeError(f"Found records but extracted zero integer DOFs from {path}")
    return mapping

def load_petsc_mat_binary(path: Path):
    try:
        from petsc4py import PETSc
    except ImportError as e:
        raise RuntimeError(
            "petsc4py is required to read PETSc binary matrices.\n"
            "In your (moose) conda env, try:\n"
            "  conda install -c conda-forge petsc4py\n"
        ) from e

    viewer = PETSc.Viewer().createBinary(str(path), "r")
    A = PETSc.Mat().load(viewer)
    return A

def find_empty_rows(A, treat_near_zero_as_empty=False, eps=0.0):
    empty = []
    rstart, rend = A.getOwnershipRange()

    for i in range(rstart, rend):
        cols, vals = A.getRow(i)
        nnz = len(cols)
        if nnz == 0:
            empty.append(i)
        elif treat_near_zero_as_empty and eps > 0.0:
            if all(abs(v) <= eps for v in vals):
                empty.append(i)
        A.restoreRow(i, cols, vals)

    # MPI gather (serial: just returns)
    try:
        comm = A.comm
        rank = comm.getRank()
        all_empty = comm.gather(empty, root=0)
        if rank == 0:
            merged = sorted({i for part in all_empty for i in part})
            return merged
        return None
    except Exception:
        return sorted(empty)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mat", required=True, help="PETSc binary Jacobian file, e.g. J.petsc")
    ap.add_argument("--dofmap", required=True, help="DOFMapOutput JSON file, e.g. RZ1_out.json")
    ap.add_argument("--out", default="empty_rows_report.txt")
    ap.add_argument("--near-zero", action="store_true",
                    help="Also flag rows whose entries are all |v|<=eps (requires --eps).")
    ap.add_argument("--eps", type=float, default=0.0, help="Threshold for --near-zero.")
    args = ap.parse_args()

    mat_path = Path(args.mat)
    dof_path = Path(args.dofmap)
    out_path = Path(args.out)

    A = load_petsc_mat_binary(mat_path)
    dofmap = load_dofmap_json(dof_path)

    empty = find_empty_rows(A, treat_near_zero_as_empty=args.near_zero, eps=args.eps)
    if empty is None:  # MPI non-root
        return

    lines = []
    lines.append(f"Matrix: {mat_path}  size = {A.getSize()}\n")
    lines.append(f"DOF map: {dof_path}\n")
    lines.append(f"Found {len(empty)} empty rows (nnz==0).\n\n")

    for row in empty:
        meta = dofmap.get(row)
        if meta is None:
            lines.append(f"row {row}: (no dofmap entry found)\n")
            continue
        lines.append(
            f"row {row}: var={meta.get('variable')} subdomain={meta.get('subdomain')} "
            f"kernels={meta.get('kernels')}\n"
        )

    out_path.write_text("".join(lines))
    print("".join(lines[: min(len(lines), 60)]))  # preview
    print(f"\nWrote full report to: {out_path.resolve()}")

if __name__ == "__main__":
    main()

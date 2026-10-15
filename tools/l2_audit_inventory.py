"""Compact per-comparison inventory of the L2 groups (reference kind, bound, stored margin, run).

Read-only view of ``docs/validation/rows``; it assigns no verdict. Usage:
    python tools/l2_audit_inventory.py [--groups 4 6 ...] [--csv out.csv]
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
L2_GROUPS = ["4", "6", "7", "8", "9", "17", "18", "21", "24", "26", "36", "43", "55"]


def inventory(groups: list[str]) -> list[dict]:
    out = []
    for g in groups:
        for path in (ROOT / "docs" / "validation" / "rows").glob(f"{g}-*.yaml"):
            for row in yaml.safe_load(path.read_text())["rows"]:
                for cmp_ in row.get("comparisons", []):
                    ref, tol, meas = cmp_.get("reference") or {}, cmp_.get("tolerance") or {}, cmp_.get("measured") or {}
                    out.append(
                        {
                            "group": g,
                            "row": row["id"],
                            "ref_kind": ref.get("kind", ""),
                            "ref": (ref.get("label") or "")[:60],
                            "tol": f"{tol.get('kind', '')} {tol.get('value', '')}",
                            "source": tol.get("source", ""),
                            "status": meas.get("status", ""),
                            "margin_pct": meas.get("margin_pct", ""),
                            "run": meas.get("run", ""),
                            "date": str(meas.get("date", "")),
                        }
                    )
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--groups", nargs="*", default=L2_GROUPS)
    ap.add_argument("--csv")
    args = ap.parse_args()
    rows = inventory(args.groups)
    w = csv.DictWriter(open(args.csv, "w", newline="") if args.csv else sys.stdout, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)


if __name__ == "__main__":
    main()

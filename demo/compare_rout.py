#!/usr/bin/env python3
"""
Compare serial vs MPI .Rout outputs for numerical equivalence (ignoring timing).

Usage:
  python compare_rout.py serial n_2
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

TIME_PATTERNS = (
    "Estimation Time",
    "Elapsed time",
    "user  system elapsed",
    "Timing stopped at",
)

FLOAT_RE = re.compile(r"[-+]?(?:\\d+\\.\\d+|\\d+\\.|\\.\\d+|\\d+)(?:[eE][-+]?\\d+)?")


def read_filtered_numbers(path: Path) -> list[float]:
    nums: list[float] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if any(pat in line for pat in TIME_PATTERNS):
                continue
            for m in FLOAT_RE.findall(line):
                try:
                    nums.append(float(m))
                except ValueError:
                    continue
    return nums


def compare_files(serial_path: Path, mpi_path: Path, atol: float = 1e-8, rtol: float = 1e-6):
    s_nums = read_filtered_numbers(serial_path)
    m_nums = read_filtered_numbers(mpi_path)

    if len(s_nums) != len(m_nums):
        return {
            "status": "count_mismatch",
            "serial_count": len(s_nums),
            "mpi_count": len(m_nums),
            "max_abs": None,
            "max_rel": None,
        }

    max_abs = 0.0
    max_rel = 0.0
    for a, b in zip(s_nums, m_nums):
        diff = abs(a - b)
        max_abs = max(max_abs, diff)
        denom = max(abs(a), abs(b), 1.0)
        max_rel = max(max_rel, diff / denom)

    status = "ok" if (max_abs <= atol or max_rel <= rtol) else "diff"
    return {
        "status": status,
        "serial_count": len(s_nums),
        "mpi_count": len(m_nums),
        "max_abs": max_abs,
        "max_rel": max_rel,
    }


def main():
    if len(sys.argv) != 3:
        print("Usage: compare_rout.py <serial_dir> <mpi_dir>")
        sys.exit(2)

    serial_dir = Path(sys.argv[1])
    mpi_dir = Path(sys.argv[2])

    serial_files = sorted(serial_dir.glob("*.Rout"))
    if not serial_files:
        print(f"No .Rout files found in {serial_dir}")
        sys.exit(1)

    rows = []
    for s in serial_files:
        mpi = mpi_dir / s.name.replace("_serial", "_npRmpi")
        if not mpi.exists():
            rows.append((s.name, "missing", "-", "-", "-"))
            continue
        res = compare_files(s, mpi)
        rows.append(
            (
                s.name,
                res["status"],
                res["serial_count"],
                res["mpi_count"],
                res["max_abs"],
                res["max_rel"],
            )
        )

    header = ("file", "status", "serial_n", "mpi_n", "max_abs", "max_rel")
    print("{:<28} {:<14} {:>8} {:>8} {:>12} {:>12}".format(*header))
    for row in rows:
        print(
            "{:<28} {:<14} {:>8} {:>8} {:>12} {:>12}".format(
                row[0],
                row[1],
                row[2] if row[2] is not None else "-",
                row[3] if row[3] is not None else "-",
                f\"{row[4]:.3g}\" if isinstance(row[4], float) else row[4],
                f\"{row[5]:.3g}\" if isinstance(row[5], float) else row[5],
            )
        )


if __name__ == "__main__":
    main()

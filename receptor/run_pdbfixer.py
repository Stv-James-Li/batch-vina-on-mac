#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys


def fix_pdb(input_pdb, output_pdb, ph=7.0, keep_water=False):
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError as exc:
        raise RuntimeError(
            "缺少 pdbfixer/openmm，请先在 autodock 环境安装："
            "mamba install -c conda-forge pdbfixer openmm"
        ) from exc

    fixer = PDBFixer(filename=input_pdb)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=keep_water)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)

    os.makedirs(os.path.dirname(output_pdb) or ".", exist_ok=True)
    with open(output_pdb, "w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)


def main():
    parser = argparse.ArgumentParser(description="Run PDBFixer before receptor preparation.")
    parser.add_argument("input_pdb", help="原始受体 PDB")
    parser.add_argument("output_pdb", help="PDBFixer 修复后的 PDB")
    parser.add_argument("--ph", type=float, default=7.0, help="加氢 pH，默认 7.0")
    parser.add_argument("--keep-water", action="store_true", help="保留水分子")
    args = parser.parse_args()

    try:
        fix_pdb(args.input_pdb, args.output_pdb, ph=args.ph, keep_water=args.keep_water)
    except Exception as exc:
        print(f"PDBFixer failed for {args.input_pdb}: {exc}", file=sys.stderr)
        return 1

    print(f"PDBFixer fixed: {args.input_pdb} → {args.output_pdb}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

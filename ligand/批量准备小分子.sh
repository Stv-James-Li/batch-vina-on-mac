#!/bin/bash
# 2D SDF/PDB → 3D SDF → PDBQT

source "$(dirname "$0")/../activate_env.sh"

# PDB → 3D SDF
for f in *.pdb; do
    [ -f "$f" ] || continue
    name="${f%.pdb}"
    echo "Generating 3D: $f → ${name}_3d.sdf"
    obabel "$f" -O "${name}_3d.sdf" --gen3d -h
done

# 2D SDF → 3D SDF（重命名输出避免覆盖原文件）
for f in *.sdf; do
    [ -f "$f" ] || continue
    # 跳过已经生成的 _3d.sdf
    [[ "$f" == *_3d.sdf ]] && continue
    name="${f%.sdf}"
    echo "Generating 3D: $f → ${name}_3d.sdf"
    obabel "$f" -O "${name}_3d.sdf" --gen3d -h
done

# 3D SDF → PDBQT
for f in *_3d.sdf; do
    [ -f "$f" ] || continue
    name="${f%_3d.sdf}"
    echo "Preparing ${name}.pdbqt"
    mk_prepare_ligand.py -i "$f" -o "${name}.pdbqt"
done

echo "小分子准备完成。"
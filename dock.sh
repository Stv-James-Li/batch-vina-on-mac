#!/bin/bash
# 批量分子对接主脚本
# 在项目根目录下运行：bash dock.sh

source "$(dirname "$0")/activate_env.sh"

mkdir -p results

for receptor in receptor/*.pdbqt; do
    [ -f "$receptor" ] || continue
    receptor_name=$(basename "$receptor" .pdbqt)

    for ligand in ligand/*.pdbqt; do
        [ -f "$ligand" ] || continue
        ligand_name=$(basename "$ligand" .pdbqt)

        echo "Running Vina: ${receptor_name}  +  ${ligand_name}"
        vina \
            --receptor "$receptor" \
            --ligand "$ligand" \
            --exhaustiveness 25 \
            --config "receptor/${receptor_name}.txt" \
            --out "results/${receptor_name}_${ligand_name}.pdbqt" \
            >> "results/${receptor_name}_${ligand_name}.txt"

        sleep 1
    done
done

echo "所有对接任务完成。"

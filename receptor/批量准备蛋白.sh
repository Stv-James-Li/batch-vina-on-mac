#!/bin/bash
# 批量准备受体蛋白（PDB → PDBFixer → PDBQT）
# 在 receptor/ 目录下运行：bash 批量准备蛋白.sh

source "$(dirname "$0")/../activate_env.sh"

fixed_dir="pdbfixer_fixed"
mkdir -p "$fixed_dir"

for f in *.pdb; do
    [ -f "$f" ] || continue
    name="${f%.pdb}"
    fixed_pdb="${fixed_dir}/${name}_fixed.pdb"
    prepared_prefix="${name}_prepared"
    prepared_pdbqt="${prepared_prefix}.pdbqt"
    final_pdbqt="${name}.pdbqt"

    echo "Running PDBFixer: $f → $fixed_pdb"
    python3 run_pdbfixer.py "$f" "$fixed_pdb"
    if [ $? -ne 0 ]; then
        echo "错误：PDBFixer 处理失败，跳过 $f"
        continue
    fi

    echo "Preparing receptor: $fixed_pdb → $final_pdbqt"
    rm -f "$prepared_pdbqt"
    mk_prepare_receptor.py -i "$fixed_pdb" -o "$prepared_prefix" --default_altloc A -p -a
    if [ ! -f "$prepared_pdbqt" ]; then
        echo "错误：未生成 $prepared_pdbqt，跳过 $f"
        continue
    fi
    mv "$prepared_pdbqt" "$final_pdbqt"
    sleep 1
done

python3 box.py

echo "受体准备完成。"
echo "提示：如果目录里还有 *.pdbqt.pdbqt，是旧版脚本生成的异常文件名；确认不需要后可手动删除。"

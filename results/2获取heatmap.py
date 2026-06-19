#!/usr/bin/env python3
# 在 results/ 目录下运行：python3 2获取heatmap.py

import os
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(script_dir, "results.txt")

data = []
with open(file_path, "r", encoding="utf-8") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        try:
            match = re.match(r"([^,]+),\s*([-\d.]+)", line)
            if match:
                receptor_ligand = match.group(1).strip()
                affinity = float(match.group(2).strip())
                parts = receptor_ligand.split("_")
                if len(parts) == 2:
                    receptor, ligand = parts
                    data.append([receptor.strip(), ligand.strip(), affinity])
                else:
                    print(f"无法分离受体和配体: {receptor_ligand}")
            else:
                print(f"行格式不正确: {line}")
        except Exception as e:
            print(f"处理行时出错: {line}, 错误: {e}")

print(f"成功解析数据点: {len(data)}")
unique_receptors = set(row[0] for row in data)
unique_ligands   = set(row[1] for row in data)
print(f"受体: {len(unique_receptors)}  配体: {len(unique_ligands)}")

df = pd.DataFrame(data, columns=["receptor", "ligand", "affinity"])

df_pivot = pd.pivot_table(
    df,
    values="affinity",
    index="ligand",
    columns="receptor",
    aggfunc="min"
)

vmin = df["affinity"].min()
vmax = df["affinity"].max()

plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({"font.family": "Arial", "font.size": 10})

fig, ax = plt.subplots(figsize=(16, 14))

cmap = sns.diverging_palette(10, 240, as_cmap=True)
sns.heatmap(
    df_pivot,
    cmap=cmap,
    annot=True,
    fmt=".2f",
    linewidths=0.5,
    cbar_kws={"label": "Affinity (kcal/mol)", "shrink": 0.8},
    vmin=vmin,
    vmax=vmax,
    annot_kws={"size": 8},
    ax=ax
)

plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(
    np.arange(len(df_pivot.index)) + 0.5,
    df_pivot.index,
    fontsize=9
)
ax.yaxis.set_major_locator(plt.FixedLocator(np.arange(len(df_pivot.index)) + 0.5))
ax.yaxis.set_major_formatter(plt.FixedFormatter(df_pivot.index))

plt.subplots_adjust(bottom=0.3, left=0.2)
plt.title("Vina Affinity Score", fontsize=16, pad=20)
plt.xlabel("Receptor", fontsize=12, labelpad=10)
plt.ylabel("Ligand",   fontsize=12, labelpad=10)
plt.tight_layout()

output_path = os.path.join(script_dir, "heatmap.png")
plt.savefig(output_path, dpi=300, bbox_inches="tight")
plt.show()
print(f"热图已保存至 {output_path}")

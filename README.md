# AutoDock Batch Docking Pipeline for macOS Apple Silicon

[中文](#中文说明) | [English](#english)

A lightweight batch molecular docking workflow for macOS Apple Silicon. The pipeline prepares receptors and ligands, runs AutoDock Vina in batch mode, extracts docking scores, generates heatmaps, builds receptor-ligand complex PDB files, and optionally runs PLIP interaction analysis with PyMOL session output.

---

## English

### Features

- Batch receptor preparation from `.pdb` to `.pdbqt`
- PDBFixer preprocessing before receptor preparation
- Batch ligand preparation from `.sdf` or `.pdb` to `.pdbqt`
- AutoDock Vina batch docking for all receptor-ligand combinations
- Docking score extraction and heatmap generation
- Receptor-ligand complex PDB reconstruction from Vina output
- Optional PLIP interaction analysis for each complex
- Optional PyMOL `.pse` output from PLIP when supported by the local PLIP installation

### Tested Platform

- macOS on Apple Silicon
- Miniforge / Mambaforge
- Python 3.11
- AutoDock Vina 1.2.7 macOS ARM64 binary

Linux should also work with minor changes to the Vina installation step.

### Repository Structure

```text
.
├── activate_env.sh                 # Shared conda environment activation helper
├── dock.sh                         # Batch AutoDock Vina docking script
├── combine.py                      # Build complex PDB files and run PLIP
├── receptor/
│   ├── 批量准备蛋白.sh              # Batch receptor preparation script
│   ├── run_pdbfixer.py             # PDBFixer preprocessing helper
│   ├── box.py                      # Generate Vina config files from receptor PDBQT files
│   └── *.pdb                       # Input receptor PDB files
├── ligand/
│   ├── 批量准备小分子.sh            # Batch ligand preparation script
│   └── *.sdf / *.pdb               # Input ligand structures
└── results/
    ├── 1提取rank1自由能.sh          # Extract rank-1 docking scores
    └── 2获取heatmap.py              # Generate docking score heatmap
```

Generated files are intentionally not required as source files. Common outputs include:

```text
receptor/pdbfixer_fixed/*_fixed.pdb
receptor/*.pdbqt
receptor/*.txt
ligand/*_3d.sdf
ligand/*.pdbqt
results/*.pdbqt
results/*.txt
results/results.txt
results/heatmap.png
complex_pdbs_plip/*_complex.pdb
complex_pdbs_plip/plip/<complex_name>/
```

### Installation

#### 1. Install Miniforge

```bash
curl -LO https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3-MacOSX-arm64.sh
```

Restart the terminal after installation if `conda` is not immediately available.

#### 2. Create the conda environment

```bash
mamba create -n autodock -c conda-forge python=3.11
mamba activate autodock
mamba install -c conda-forge openbabel pymol-open-source numpy pandas seaborn matplotlib pdbfixer openmm
pip install meeko plip
```

#### 3. Install AutoDock Vina

```bash
curl -LO https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_mac_arm64
chmod +x vina_1.2.7_mac_arm64
sudo mv vina_1.2.7_mac_arm64 /usr/local/bin/vina
xattr -dr com.apple.quarantine /usr/local/bin/vina
```

Check the installation:

```bash
vina --help
```

### Usage

#### One-command workflow

After placing receptor `.pdb` files in `receptor/` and ligand `.sdf` or `.pdb` files in `ligand/`, copy and paste the following block into the terminal:

```bash
cd ~/Desktop/autodock-final && \
cd receptor && bash 批量准备蛋白.sh && \
cd ../ligand && bash 批量准备小分子.sh && \
cd .. && bash dock.sh && \
cd results && bash 1提取rank1自由能.sh && python3 2获取heatmap.py && \
cd .. && python3 combine.py
```

If PLIP is not installed and you only want to generate complex PDB files, replace the last line with:

```bash
cd .. && python3 combine.py --no-plip
```

#### Step-by-step workflow

#### Step 1. Prepare receptors

Place receptor `.pdb` files in `receptor/`, then run:

```bash
cd receptor
bash 批量准备蛋白.sh
```

This step:

1. runs PDBFixer on each receptor PDB;
2. writes fixed structures to `receptor/pdbfixer_fixed/`;
3. prepares receptor `.pdbqt` files;
4. generates Vina box config files through `box.py`.

#### Step 2. Prepare ligands

Place ligand `.sdf` or `.pdb` files in `ligand/`, then run:

```bash
cd ../ligand
bash 批量准备小分子.sh
```

This step generates 3D ligand structures and prepares ligand `.pdbqt` files.

#### Step 3. Run batch docking

From the repository root:

```bash
cd ..
bash dock.sh
```

Docking outputs are written to `results/` as:

```text
results/<receptor>_<ligand>.pdbqt
results/<receptor>_<ligand>.txt
```

#### Step 4. Extract docking scores

```bash
cd results
bash 1提取rank1自由能.sh
```

The extracted score table is written to `results/results.txt`.

#### Step 5. Generate heatmap

```bash
python3 2获取heatmap.py
```

The heatmap is written to `results/heatmap.png`.

#### Step 6. Build complexes and run PLIP

From the repository root:

```bash
cd ..
python3 combine.py
```

This step:

1. takes each docked ligand pose from `results/*.pdbqt`;
2. combines it with the corresponding receptor `.pdbqt`;
3. writes standard receptor-ligand complex PDB files to `complex_pdbs_plip/`;
4. runs PLIP for each complex if the `plip` command is available.

To build complex PDB files without running PLIP:

```bash
python3 combine.py --no-plip
```

Custom folders can be specified:

```bash
python3 combine.py \
  --receptor-folder receptor \
  --results-folder results \
  --output-folder complex_pdbs_plip \
  --plip-folder complex_pdbs_plip/plip
```

### Notes on Docking Box Generation

The current `receptor/box.py` script calculates the docking box from the full receptor PDBQT coordinate range. This is useful for blind docking, but the search box may be large. For site-specific docking, edit the generated `receptor/<receptor>.txt` file manually and set a focused center and box size.

Example Vina config:

```text
center_x=166.885
center_y=151.703
center_z=183.600
size_x=25
size_y=25
size_z=25
num_modes=9
```

### Troubleshooting

#### `pdbfixer` is missing

Install it inside the `autodock` environment:

```bash
mamba install -c conda-forge pdbfixer openmm
```

#### `plip` command is not found

Install PLIP:

```bash
pip install plip
```

Then check:

```bash
plip -h
```

#### Receptor output becomes `*.pdbqt.pdbqt`

This happens when `mk_prepare_receptor.py` is given an output prefix that already includes `.pdbqt`. The included receptor preparation script avoids this by using a temporary prefix and then renaming the final file to `<receptor>.pdbqt`.

#### Complex PDB appears misaligned in PyMOL

Regenerate complexes with the current `combine.py`:

```bash
python3 combine.py --no-plip
```

The script writes standard fixed-width PDB `ATOM` and `HETATM` lines so PyMOL and PLIP can parse coordinates correctly.

### Recommended `.gitignore`

If you plan to publish only the workflow scripts and examples, consider ignoring generated outputs:

```gitignore
.DS_Store
__pycache__/
*.pyc

receptor/pdbfixer_fixed/
receptor/*.pdbqt
receptor/*.txt
ligand/*_3d.sdf
ligand/*.pdbqt
results/*.pdbqt
results/*.txt
results/*.png
complex_pdbs_plip/
```

Keep small example inputs only if redistribution is allowed.

---

## 中文说明

这是一个面向 macOS Apple Silicon 的轻量级 AutoDock Vina 批量分子对接流程。流程支持受体和配体批量准备、Vina 批量对接、自由能提取、热图生成、受体-配体复合物 PDB 重建，以及可选的 PLIP 相互作用分析和 PyMOL `.pse` 输出。

### 功能特点

- 批量将受体 `.pdb` 准备为 `.pdbqt`
- 在受体准备前自动运行 PDBFixer
- 批量将配体 `.sdf` 或 `.pdb` 准备为 `.pdbqt`
- 对所有受体-配体组合运行 AutoDock Vina
- 提取 rank-1 docking score 并生成热图
- 根据 Vina 输出重建标准 PDB 格式的复合物
- 可选：对每个复合物运行 PLIP 相互作用分析
- 可选：在 PLIP 支持时输出 PyMOL `.pse` 会话文件

### 测试平台

- Apple Silicon Mac
- Miniforge / Mambaforge
- Python 3.11
- AutoDock Vina 1.2.7 macOS ARM64 版本

Linux 理论上也可以使用，只需要根据系统调整 Vina 安装方式。

### 目录结构

```text
.
├── activate_env.sh                 # 自动激活 autodock conda 环境
├── dock.sh                         # 批量 AutoDock Vina 对接脚本
├── combine.py                      # 生成复合物 PDB，并运行 PLIP
├── receptor/
│   ├── 批量准备蛋白.sh              # 批量准备受体
│   ├── run_pdbfixer.py             # PDBFixer 预处理脚本
│   ├── box.py                      # 根据受体 PDBQT 生成 Vina config
│   └── *.pdb                       # 输入受体 PDB 文件
├── ligand/
│   ├── 批量准备小分子.sh            # 批量准备小分子
│   └── *.sdf / *.pdb               # 输入配体结构
└── results/
    ├── 1提取rank1自由能.sh          # 提取 rank-1 自由能
    └── 2获取heatmap.py              # 生成 docking score 热图
```

常见输出文件包括：

```text
receptor/pdbfixer_fixed/*_fixed.pdb
receptor/*.pdbqt
receptor/*.txt
ligand/*_3d.sdf
ligand/*.pdbqt
results/*.pdbqt
results/*.txt
results/results.txt
results/heatmap.png
complex_pdbs_plip/*_complex.pdb
complex_pdbs_plip/plip/<complex_name>/
```

### 安装

#### 1. 安装 Miniforge

```bash
curl -LO https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3-MacOSX-arm64.sh
```

如果安装后终端暂时找不到 `conda`，请重启终端。

#### 2. 创建 conda 环境

```bash
mamba create -n autodock -c conda-forge python=3.11
mamba activate autodock
mamba install -c conda-forge openbabel pymol-open-source numpy pandas seaborn matplotlib pdbfixer openmm
pip install meeko plip
```

#### 3. 安装 AutoDock Vina

```bash
curl -LO https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_mac_aarch64
chmod +x vina_1.2.7_mac_aarch64
sudo mv vina_1.2.7_mac_aarch64 /usr/local/bin/vina
xattr -dr com.apple.quarantine /usr/local/bin/vina
```

检查是否安装成功：

```bash
vina --help
```

### 使用流程

#### 一键运行命令

将受体 `.pdb` 文件放入 `receptor/`，将配体 `.sdf` 或 `.pdb` 文件放入 `ligand/` 后，可以直接复制下面整段命令到终端运行：

```bash
cd ~/Desktop/autodock-final && \
cd receptor && bash 批量准备蛋白.sh && \
cd ../ligand && bash 批量准备小分子.sh && \
cd .. && bash dock.sh && \
cd results && bash 1提取rank1自由能.sh && python3 2获取heatmap.py && \
cd .. && python3 combine.py
```

如果没有安装 PLIP，只想生成复合物 PDB，可以把最后一行替换为：

```bash
cd .. && python3 combine.py --no-plip
```

#### 分步运行流程

#### 第 1 步：准备受体

将受体 `.pdb` 文件放入 `receptor/`，然后运行：

```bash
cd receptor
bash 批量准备蛋白.sh
```

该步骤会：

1. 对每个受体 PDB 运行 PDBFixer；
2. 将修复后的结构写入 `receptor/pdbfixer_fixed/`；
3. 生成受体 `.pdbqt`；
4. 通过 `box.py` 生成 Vina config 文件。

#### 第 2 步：准备配体

将配体 `.sdf` 或 `.pdb` 文件放入 `ligand/`，然后运行：

```bash
cd ../ligand
bash 批量准备小分子.sh
```

该步骤会生成 3D 配体结构，并准备配体 `.pdbqt` 文件。

#### 第 3 步：批量对接

回到项目根目录：

```bash
cd ..
bash dock.sh
```

对接输出位于 `results/`：

```text
results/<receptor>_<ligand>.pdbqt
results/<receptor>_<ligand>.txt
```

#### 第 4 步：提取 docking score

```bash
cd results
bash 1提取rank1自由能.sh
```

结果表格会写入 `results/results.txt`。

#### 第 5 步：生成热图

```bash
python3 2获取heatmap.py
```

热图会写入 `results/heatmap.png`。

#### 第 6 步：生成复合物并运行 PLIP

回到项目根目录：

```bash
cd ..
python3 combine.py
```

该步骤会：

1. 读取 `results/*.pdbqt` 中的 Vina 最优构象；
2. 与对应的受体 `.pdbqt` 合并；
3. 输出标准 PDB 格式的复合物到 `complex_pdbs_plip/`；
4. 如果系统中存在 `plip` 命令，则对每个复合物运行 PLIP。

如果只想生成复合物 PDB，不运行 PLIP：

```bash
python3 combine.py --no-plip
```

也可以自定义目录：

```bash
python3 combine.py \
  --receptor-folder receptor \
  --results-folder results \
  --output-folder complex_pdbs_plip \
  --plip-folder complex_pdbs_plip/plip
```

### 关于 docking box

当前 `receptor/box.py` 会根据整个受体 PDBQT 的坐标范围生成 docking box，适合 blind docking，但 box 可能较大。如果需要特定位点对接，建议手动修改生成的 `receptor/<receptor>.txt`，设置更聚焦的中心和范围。

示例：

```text
center_x=166.885
center_y=151.703
center_z=183.600
size_x=25
size_y=25
size_z=25
num_modes=9
```

### 常见问题

#### 缺少 `pdbfixer`

在 `autodock` 环境中安装：

```bash
mamba install -c conda-forge pdbfixer openmm
```

#### 找不到 `plip` 命令

安装 PLIP：

```bash
pip install plip
```

检查：

```bash
plip -h
```

#### 受体输出变成 `*.pdbqt.pdbqt`

这是因为 `mk_prepare_receptor.py` 的 `-o` 参数使用了已经带 `.pdbqt` 后缀的输出名前缀。本项目中的受体准备脚本会使用临时前缀，然后统一重命名为 `<receptor>.pdbqt`，以避免该问题。

#### PyMOL 中复合物看起来错位

请用当前版本的 `combine.py` 重新生成复合物：

```bash
python3 combine.py --no-plip
```

当前脚本会输出标准定宽 PDB `ATOM` 和 `HETATM` 行，便于 PyMOL 和 PLIP 正确解析坐标。

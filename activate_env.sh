#!/bin/bash
# 公共：激活 autodock conda 环境
# 被其他脚本 source 调用，不要直接运行

CONDA_BASE=$(conda info --base 2>/dev/null)
if [ -z "$CONDA_BASE" ]; then
    echo "错误：找不到 conda，请先安装 Miniforge"
    exit 1
fi

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate autodock

if [ $? -ne 0 ]; then
    echo "错误：无法激活 autodock 环境，请先运行："
    echo "  mamba create -n autodock -c conda-forge python=3.11"
    exit 1
fi

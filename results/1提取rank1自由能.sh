#!/bin/bash
# 提取每个对接结果中 rank 1 的结合自由能
# 在 results/ 目录下运行：bash 1提取rank1自由能.sh

source "$(dirname "$0")/../activate_env.sh"

> results.txt

for f in *.txt; do
    [ -f "$f" ] || continue
    [[ "$f" == "results.txt" ]] && continue

    name="${f%.txt}"
    affinity=$(awk 'NR>4 && $1=="1" {print $2; exit}' "$f")

    if [ -n "$affinity" ]; then
        echo "${name},${affinity}" >> results.txt
    fi
done

echo "完成，结果已保存至 results.txt"

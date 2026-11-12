# SEM02b 独立重建程序与实际运行记录

本程序包完成指定模型 D0 块 0–64、D2 块 65–268 的完整定义级重建。
所有 941 个变量与 269 个块全部通过；最终比较 1,380,389 个非零系数键。
完整结论、数学推导、有效性证明、独立性披露和实际耗时见 `SEM02b_report.md`。

## 完整复现

环境：Linux、Python 3.10+ 标准库、支持 `__int128` 的 g++。本次实跑为
Python 3.12.14、g++ 14.2.0。没有网络依赖，不需要安装第三方 Python 包。
进入本 README 所在目录后运行：

```bash
python3 run_all.py --workers 4
```

默认从数学定义重新计算，先变量表、再 D0、再 D2，不使用旧矩阵系数。
可以用 `--workers 1` 串行计算 D0。结果写到 `sem02b/results/`。
如需保留交付时的原始日志与计时记录，请先保留一份解压副本。

只回读、逐项比较已交付的系数文件：

```bash
python3 run_all.py --check-recorded
```

这个选项不会重建系数，不能把它作为定义级重建运行。

## 单独重建

```bash
python3 sem02b/orbits.py --model inputs/A18_6_6_model.json.gz --out sem02b/results
python3 sem02b/check_labels.py --model inputs/A18_6_6_model.json.gz --out sem02b/results/labels.json
python3 sem02b/theory_selftest.py
python3 sem02b/d0.py --model inputs/A18_6_6_model.json.gz --blocks all --workers 4
python3 sem02b/check_omitted_shapes.py
python3 sem02b/d2.py --model inputs/A18_6_6_model.json.gz --lookup sem02b/results/orbit_lookup.tsv --output sem02b/results
python3 sem02b/audit_results.py --model inputs/A18_6_6_model.json.gz
```

D0 可用 `--blocks 57,60,64` 指定小阶先行块；D2 可用
`--blocks 103,179,214,230,268` 指定小阶先行块。小批次运行不等于完整核验。
本包中的 D0 `summary.json` 和 D2 `d2_summary.json` 都来自完整运行；先行结果另保留。

## 关键数据文件

- `sem02b/results/orbits.json`：全部轨道和变量比较。
- `sem02b/results/canonical_variables.json`：独立生成的规范形。
- `sem02b/results/orbit_lookup.tsv`：含重复词折叠的完整四词直方图查找表。
- `sem02b/results/d0/block_*.raw.tsv`：选列前、约分前全部非零系数。
- `sem02b/results/d0/block_*.rebuilt.json.gz`：完整 D0 重建条目、g、K。
- `sem02b/results/d0/block_*.result.json`：完整 Gram、精确消元主元和比较记录。
- `sem02b/results/d2_reconstructed.json.gz`：全部 D2 重建条目。
- `sem02b/results/d2_summary.json`：D2 全部 Gram、主元、g、K、比较与计时。
- `sem02b/results/block_index.json` / `.csv`：全部 269 块的紧凑索引。
- `sem02b/results/delivery_audit.json`：交付文件回读、严格整数类型及全项比较。
- `sem02b/results/omitted_shapes/`：五个额外候选形状的零投影检查。

所有数学计算都是精确整数或有理数。D0 使用带溢出检查的 signed 128-bit
整数，并在报告中给出严格小于 2^111 的中间贡献界；Gram 消元使用任意精度
Fraction。D2 使用 Python 任意精度整数和 Fraction。浮点只用于计时。

`MANIFEST_SHA256.txt` 记录输入、源程序、报告和实际结果文件的摘要。
它用于交付完整性核对；矩阵通过判据是全部键和值直接相等，绝非只比哈希。
压缩包未包含已编译的本机可执行文件；复现入口会从所附 C++ 源码编译。

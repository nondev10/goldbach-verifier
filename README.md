# Goldbach Verifier

大规模验证强哥德巴赫猜想（每个 n ≥ 4 的偶数都是两个素数之和）的计算验证器。

**状态：未证明哥德巴赫猜想。** 本仓库交付的是到 n = 10^13 为止的零反例计算验证、
Hardy–Littlewood 渐近结构的精确验证，以及一份刻画现有证明方法边界的笔记。

## 已验证区间

| 区间 | 偶数个数 | 耗时 | 速度 | 反例 |
|------|---------|------|------|------|
| 4 ~ 10^10 | 4,999,999,999 | 121.5s | 41.1M/s | 0 |
| 4 ~ 10^11 | 49,999,999,999 | 1773.5s | 28.2M/s | 0 |
| 10^11 ~ 10^12 | 450,000,000,000 | 3057.7s | 147.2M/s | 0 |
| 10^12 ~ 10^13 | 4,500,000,000,000 | 31944.6s | 140.9M/s | 0 |

每段都做 Miller–Rabin 独立交叉验证（10^13 段抽样 2,145,772/2,145,772 通过）。

## 核心算法

位并行分段筛。对固定素数 p，64 个连续偶数的补数 q_j = n0 - p + 2j 在奇数索引空间
是连续的，因此一次 64 位移位加载即可同时判定 64 个偶数。

## Hardy–Littlewood 渐近验证

`src/hl_ntt.c` 用 NTT（模 2013265921 = 15·2²⁷+1）精确计算每个偶数的表示数 R(n)，
与 HL 理论预测逐 n 比对（n ≈ 5×10⁶，19842 样本）：

```
mean R/pred = 0.999910    stddev = 0.004564
geometric mean = 0.999899
evens with zero representations = 0
```

主项必须用积分形式 `2·C₂·S(n)·∫₂^{n/2} dp/(ln p · ln(n−p))`；
教科书的 `n/ln²n` 粗略式有 O(1/ln n) 相对误差，n≈5×10⁶ 处系统偏低约 15%。

## 构建与运行

```bash
gcc -O3 -march=native -fopenmp -o goldbach4 src/goldbach4.c -lm
gcc -O3 -march=native -fopenmp -o hl_ntt     src/hl_ntt.c    -lm

./goldbach4 4 1000000000000 2 2000     # <start> <end> <threads> <PMAX>
./hl_ntt 5000000 2 20000               # HL 渐近校验
bash src/test.sh                        # 回归测试
```

## 复现性

- `src/independent_check.py`：完全独立的 Python 实现，用于交叉校验
- `src/test.sh`：6 组回归测试（C vs Python 逐位校验和、多线程不变性、NTT vs 暴力、
  增量区间一致性等）
- `logs/run_1eNN.json`：机器可读结果（区间 / 计数 / 校验和 / 耗时）

## 证明尝试笔记

见 `PROOF_ATTEMPT.md`：圆周方法（卡在二元相位控制）、筛法（卡在 parity problem）、
随机模型（反而预测反例存在）。

**要真正证明 GC，需要的不是更快的计算机，而是 30 年来没人想到的新数学。**

## 人类验证记录

4×10¹⁸（Oliveira e Silva, Herzog, Pardi, 2014；585 核 / 32TB 内存）。

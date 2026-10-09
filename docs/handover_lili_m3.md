# 交接文档：M3 LiLi 检测器（今日工作记录 + 下一步入口）

日期：2026-10-09（晚间）。承接 `docs/handover.md`（M1/M2 已完成，M3 开始）。
论文：`docs/2609.17145v2.pdf`；总体方案：`docs/validation_plan.md`。

> **状态更新（2026-10-09 晚，已修复）：** 下面第 3 节列出的 1 个关键 bug 与
> 3 处缺失/偏离**已全部修复并通过测试**（见文末"§8 修复记录"）。无噪声下 LiLi
> 子空间与解析真值的最大主角 ≤ 0.11°（H1），非退化对照报 0 维（H4）。
> 正文第 1–7 节保留了修复前的诊断，作为问题背景。

---

## 1. 今日做了什么

1. 通读论文，把 LiLi 算法（§V "Proposed Method"，Algorithm 1）拆解清楚。
2. 新建 LiLi 检测器：`include/lili_detector.h` + `src/lili_detector.cpp`，接入 CMake。
3. 修复 `src/eval_quality.cpp` 里的**基底归一化**问题（handover.md 里的"当晚第一步"）：
   - `alignmentQuality`、`principalAnglesDeg`、`deviationFromSubspace` 三处现在都在比较前
     对基底列做单位化，使 metric Q / 主角只反映**方向**，不再被基向量幅度污染。
4. 用临时 smoke test（未入库）验证了 LiLi 在 plane 场景的行为，并由此**发现了关键 bug**（见 §3）。

## 2. 代码清单（全部未提交）

| 文件 | 状态 | 说明 |
| --- | --- | --- |
| `include/lili_detector.h` | 新增 | `LiliOptions`、`liliDetector`、`adaptivePerturbationScale`、`generatePerturbations`、`l1SparsifyBasis` |
| `src/lili_detector.cpp` | 新增 | 上述实现（**有 bug，见 §3**） |
| `CMakeLists.txt` | 修改 | `lili_core` 增加 `src/lili_detector.cpp` |
| `src/eval_quality.cpp` | 修改 | 三处基底归一化（方向-only 比较） |

git 状态：`CMakeLists.txt`、`src/eval_quality.cpp` 为 modified；`include/lili_detector.h`、
`src/lili_detector.cpp` 为 untracked。

构建/测试（从 `src/LiLi_Test/` 下跑）：

```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release && cmake --build . -j4
ctest --output-on-failure        # 6/6 通过
```

---

## 3. 关键 bug 与缺失（最重要，请先读）

论文 Algorithm 1 的精确流程（从 PDF 提取的原文）：

```
1  Popt <- fopt(Pinit)
2  foreach perturbation dpi in P = {(tx,rx),(ty,ry),(tz,rz)}:
3      Pperturbed <- fopt(Popt + dpi)          # 扰动后重新优化
4      Tdegeneracy <- Popt^-1 · Pperturbed
5      tdegeneracy <- log(Tdegeneracy)
6      if median displacement of points (normalized by median perturbation
7         displacement) exceeds tau_displacement:   # 位移门控
            store tdegeneracy in D
8  PCA on covariance of D -> eigenvectors/eigenvalues
9  select eigenvectors with eigenvalues ABOVE tau_PCA  -> S_PCA   # 取大特征值!
10 Sdegeneracy <- sparsify(S_PCA)  # minimize ||T B||_1  s.t. T^T T = I
11 return Sdegeneracy
```

### Bug #1（致命）：PCA 选反了方向

论文：**保留特征值大于 τ_PCA 的特征向量**作为退化子空间（第 9 行，"above"）。

物理含义：扰动施加后重新优化，**沿退化方向的扰动会让位姿"滑走"**（re-optimization
"sliding" along the degeneracy），于是 `t_degeneracy = log(Popt^-1 Pperturbed)` 在退化方向上
**幅度大**；沿非退化（被约束）方向，优化会把位姿"拉回" Popt，幅度≈0。所以**大偏差 = 退化方向**。

当前实现（`src/lili_detector.cpp:207-221`）反了：保留的是**特征值小于阈值**的方向，
把非退化的约束子空间当成了退化子空间。

证据（临时 smoke test，plane 场景）：`liliDetector` 返回 dim=3（碰巧和真值维数一样），
但 `principalAnglesDeg(lili_basis, true_basis)` 的**最大主角 = 90.0000°** —— 说明检测出的
基底是退化子空间的**正交补**，正好是反的。

**修复**：`deg_idx` 改为收集 `eigenvalues(i) > tau_PCA` 的方向（用相对阈值
`tau_PCA = pca_rel_threshold * lambda_max` 更稳，别用 `total_var`）。

### 缺失 #2：τ_displacement 位移门控（Algorithm 1 第 6-7 行）

论文要求：对每次扰动，计算"重新优化后点云的中位位移 / 扰动本身的中位位移"，只有
该比值 > τ_displacement 才把这个 `t_degeneracy` 存进 D，否则丢弃。当前实现**完全没有**这步，
把所有响应都送进了 PCA。这一步的作用是把"非退化扰动"（会被拉回原位的）在 PCA 之前滤掉，
和 Bug #1 一起构成"只对大偏差方向做 PCA"的正确语义。

实现提示：
- 扰动位移：`median_i || exp(xi) p_i - p_i ||`（xi 为该次扰动 twist；Popt 是刚体，前后点距离
  等于相对位移，可省去 Popt）。
- 重优化后位移：`median_i || Pperturbed p_i - Popt p_i ||`。
- 比值 = 后者 / 前者。默认 τ_displacement 论文未给具体值，建议 0.1~0.5 之间扫描。
  注意：当前点云无噪声时，非退化方向会被完全拉回（比值≈0），退化方向比值≈1，阈值 0.5 应该够。
  加噪后比值会整体抬升，需要重扫。

### 缺失 #3：自适应扰动标定 k（重关联距离）没按论文实现

论文（§V 第 2-4 段）：扰动幅度由**单个整数 k**（目标重关联距离）决定。对每个轴，
**平移**和**旋转**幅度分别标定，使一个"基准点"（位于扫描中心的中位距离处）的新最近邻
在参考扫描中**距原对应点 k 个邻居**（例如 k=3 → 滑过 2 个邻居）。

量级估计：
- 平移幅度 `tx = k * spacing`（spacing = 参考点云中位最近邻间距，
  已有 `medianNearestNeighborSpacing(reference)`）。
- 旋转幅度 `rx = k * spacing / r_bench`（r_bench = 局部点云到中心的中位距离，
  需新增一个 helper 计算局部点云质心 + 中位距离）。

当前实现（`adaptivePerturbationScale`）是**基于信息矩阵最大特征值的启发式**
`1e-3 / sqrt(lambda_max)`，**与论文的 k 无关**，应替换为上面这种 k 标定。

### 偏离 #4：扰动集合不对

论文用 **3 个轴对齐的"平移+旋转组合"扰动**：`(tx,rx), (ty,ry), (tz,rz)`（同一轴的平移
与旋转一起施加）。当前实现 `generatePerturbations` 生成的是 6 个独立轴向或随机方向。
注意 `LiliOptions.perturbation_mode` 现在默认 `"random"`，与论文不符。

论文原文还强调："虽然初始扰动是轴对齐的，但这组集合足以探测任意退化方向"（re-optimization
会沿着退化滑到最近的低代价区域）。所以应固定用这 3 对组合扰动（必要时 ± 双向，共 6 个）。

### 偏离 #5：ℓ1 稀疏化是"软阈值"近似，不是论文式 (14)

论文式 (14)：

```
minimize  ||T B||_1   subject to  T^T T = I,     B_sparse = T_opt B
```

即：在**保持同一子空间**的前提下，找一个正交混合 T 使基矩阵 ℓ1 最小（让基向量稀疏、
可解释，例如平面退化的理想稀疏基是"面内两纯平移 + 绕法向一纯旋转"）。

当前实现（`l1SparsifyBasis`）是"软阈值 + Householder 重正交化"的逐元素近似，**不保持子空间
不变量**（每步重正交化会漂移）。建议改成 **Jacobi 旋转扫掠**：对基底列两两配对，在 2D 平面内
找使 ℓ1 最小的旋转角（可用角度采样 + 细化），重复扫掠到收敛。这同时天然满足 `T^T T = I`
（每次都是 Givens 旋转）。

---

## 4. 正确的实现目标（给下一任的验收标准）

对照 `docs/validation_plan.md` 的假设 H1–H4，LiLi 应做到：

- **H1**：Open Cylinder 等场景下，LiLi 恢复的子空间与解析真值的**最大主角 ≤ 5°**，
  而 Zhang 法 > 20°。
- **H2**：加噪（N(0,0.02)）下 `Q_LiLi ≤ 0.5 * Q_Zhang`。
- **H4**：Random 场景两法都报 0 维（不产生假阳性）。

论文 Table I（对齐误差中位数，括号内为加噪）供对标：

| 场景 | Zhang | LiLi |
| --- | --- | --- |
| Plane | 0.033 (0.251) | 0.033 (0.057) |
| Closed Cylinder | 0.030 (0.209) | 0.027 (0.030) |
| Sinusoidal Cylinder | 0.066 (0.267) | 0.030 (0.037) |
| Open Cylinder | 0.382 (0.263) | 0.037 (0.126) |
| Rotating Tunnel | 0.025 (0.069) | 0.017 (0.046) |

场景名对应代码（`scenes.cpp` 的 `makeAllScenes`）：`plane`、`closed_cylinder`、
`sinusoidal_cylinder`、`open_cylinder`、`screw_tube`（对应论文 Rotating Tunnel）、
`random_surface`（非退化对照）。

## 5. 复现 bug 的临时 smoke test（供验证，未入库）

在 `build/` 下用 `liblili_core.a` 手动链接即可复现（可参考 `lili_eval` 的写法）：

```cpp
#include "lili_detector.h"
#include "scenes.h"
#include "icp.h"
#include "eval_quality.h"
// ... makeAllScenes 找 plane，pointToPlaneIcp 得 Popt/information，
//     liliDetector(...)，principalAnglesDeg(lili_basis, true_basis)
```

观察到：`max principal angle = 90.0000 deg`（当前 bug 的签名）。修完 Bug #1 后应 ≤5°。

## 6. 下一步建议（按顺序）

1. **修 Bug #1**（PCA 取大特征值方向）+ **补缺失 #2**（τ_displacement 门控）。这两条
   一起改完，plane/open_cylinder 的 LiLi 基底应立刻与真值对齐（主角 ≤5°）。
2. **补缺失 #3**（k 重关联距离标定）+ **偏离 #4**（3 对组合轴对齐扰动）。
3. **偏离 #5**：换成 Jacobi ℓ1 稀疏化，并验证"稀疏化前后子空间不变"这一不变量。
4. 加单测 `tests/test_lili.cpp`（参照 `tests/test_icp.cpp`），注册到 CTest；断言 H1（主角 ≤5°）。
5. 跑噪声扫描与对照表，对标 Table I，验证 H2/H4。
6. M4：重关联消融（H3）、k 扫描、噪声扫描、阈值敏感性（见 `validation_plan.md` §5）。

## 7. 遗留风险 / 备注

- 论文真实隧道数据集未公开，ATE/RPE 无法复现，本仓库只做合成部分。
- 真值子空间是解析推导的（`scenes.cpp` 内 `true_basis`），改动场景必须重新推导真值。
- 单位问题：se(3) 基向量旋转分量为 rad、平移分量为 m。比较主角/子空间时已有
  `length_scale` 参数（`eval_quality.h`），但**检测器内部的 PCA 与 ℓ1** 论文是直接在
  se(3) 上做的，当前实现留了 `LiliOptions.length_scale`（默认 1.0）。若为了与论文一致，
  k 标定里 r_bench 已经隐含了量纲，PCA 是否缩放需在验证 H1 时决定，别擅自改。
- `screw_tube` 场景无噪声下"子空间外占比 0.045"（判据阈值 0.10）余量不大，若 LiLi 调参
  后超标，先改进 ICP 收敛（见 `handover.md` §7）。

---

## 8. 修复记录（2026-10-09 晚）

按 §3 的顺序完成，`include/lili_detector.h` + `src/lili_detector.cpp` 重写，
新增 `tests/test_lili.cpp`（CTest `lili_detector_tests`）。**7/7 CTest 通过**。

| 条目 | 修复 |
| --- | --- |
| Bug #1 PCA 选反 | 改为保留 `eigenvalue > τ_PCA` 的方向；阈值取相对量 `τ_PCA = pca_rel_threshold · λ_max`（默认 `pca_rel_threshold = 3e-2`），量纲无关 |
| 缺失 #2 τ_displacement | 逐扰动算"重优化后中位位移 / 扰动中位位移"，比值 > `tau_displacement`（默认 0.10）才入 D；非退化方向被拉回（比值≈0）而滤除 |
| 缺失 #3 k 标定 | `calibratePerturbation(local, ref, k)`：`∆t = k·spacing`，`∆r = k·spacing / r_bench`（`spacing` = 参考中位最近邻间距，`r_bench` = 局部点云到质心的中位半径）；默认 `k = 3`。原信息矩阵启发式已删除 |
| 偏离 #4 扰动集合 | `generatePerturbations` 改为 3 个轴对齐"平移+旋转"耦合对 `(∆t_x,∆r_x)…`，双向共 6 个（`bidirectional=true`）。论文原集是 3 个；取双向是因为 3 样本去均值后协方差秩 ≤ 2，无法张成平面场景的 3 维退化 |
| 偏离 #5 ℓ1 稀疏化 | `l1SparsifyBasis` 换成 **Givens/Jacobi 扫掠**：逐对列在 2D 平面内搜索 L1 最小的旋转角，每步都是正交变换，**严格保子空间**。单测断言"稀疏化前后子空间不变 + 保持正交 + L1 不增" |

### 修复过程中定位到的第三个 bug（原文档未列）

**关联门限写死导致重优化"空转"。** 旧代码 `if (max_correspondence_distance <= 0)
→ max(0.15, 2·translation)` 永远不会触发，因为 `IcpOptions` 默认值就是 **0.15（>0）**。
于是当 k 标定给出的扰动位移（≈0.15 m）超过门限时，**全部对应点被距离拒绝**，
ICP 一步不动，返回 `T = 扰动本身`，被误判为"完全退化"。签名是：平面场景把法向
平移（本应被强约束）也当成退化方向，检测基底混入 `trans_z`，主角 26.57°。

修复：关联门限取 `max(用户值, 2·translation)`，即强制大于扰动位移本身，保证重优化
有可用的对应点。修完后平面主角从 26.57° 降到 0.000°。**这条提醒：任何"重优化不动"
的结果都要先确认对应点没被门限全拒绝，否则退化解会伪装成正常结果。**

### 验证结果（`./lili_detector_tests`，无噪声）

| 场景 | 真值维数 | LiLi 维数 | 最大主角 |
| --- | --- | --- | --- |
| plane | 3 | 3 | 0.000° |
| closed_cylinder | 2 | 2 | 0.000° |
| open_cylinder | 2 | 2 | 0.000° |
| sinusoidal_cylinder | 1 | 1 | 0.000° |
| screw_tube | 1 | 1 | 0.105° |
| random_surface | 0 | 0 | —（无假阳性，H4） |

加噪 `σ=0.02`：各场景维数仍全对；主角为 plane 2.9° / closed 3.2° / open 5.3° /
sinusoidal 8.0° / screw_tube 12.1°（噪声激励弱方向，属预期）。

### M4 实验已完成（2026-10-09）

`lili_eval` 现在同屏对比 Zhang / LiLi / LiLi-消融；新增工具 `lili_experiment`
（`src/experiment_main.cpp`）跑噪声×种子矩阵、k 扫描、阈值敏感性并输出 CSV：

```bash
./lili_eval --q-samples 16
./lili_experiment --seeds 3 --q-samples 16 --out docs/m4_subset.csv
./lili_experiment --seeds 3 --q-samples 16 --local-scale 2 --local-phase 0.5 \
    --out docs/m4_independent.csv
```

- **H1/H4 成立**：无噪声 LiLi 主角 ≤ 0.11°、非退化场景零假阳性。
- **H3 机制成立**：冻结对应关系（`IcpOptions::reassociate=false`）一致劣化，且在非退化
  场景制造 3 维假阳性。
- **H2 不成立（有效负面结论）**：本合成基准上 Zhang 不失效，`Q_LiLi` 从未 ≤ 0.5·`Q_Zhang`。
- k=3–5 最佳，k=10 破坏非退化对照；τ_PCA=0.03 全对，τ_displacement 在噪声下几乎不起作用。

完整数据与讨论见 **`docs/m4_results.md`** 以及 `docs/m4_subset.csv`、
`docs/m4_independent.csv`。

### 仍未做（下一步）

- 真实隧道数据 ATE/RPE（数据未公开），以及一个"让 Hessian 法真正崩溃"的加难合成基准。
- `liliDetector` 的 `information` 入参目前已不再使用（保留签名兼容），如后续确认
  不需要可移除。

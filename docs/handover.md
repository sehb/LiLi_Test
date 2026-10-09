# 交接文档：LiLi 论文算法有效性验证

日期：2026-10-09。本文件记录当日进度、关键技术决策、以及当晚继续的入口。
总方案见 `docs/validation_plan.md`，论文为 `docs/2609.17145v2.pdf`。

## 0. 一句话状态

M1、M2 完成并全部接入 CTest（6/6 通过，约 15 秒）。**M3 尚未开始**，卡在一个
前置门槛上：当前设置无法让 Zhang 基线真正失效（详见第 5 节）。今天最后定位到
一个方法论问题——$Q$ 的比较没有归一化基向量幅度，导致结论不可用。这是当晚的
第一个待办。

## 1. 仓库与构建

注意目录层次容易踩坑：工作目录是 `LiLi_Test_ws/`，**Git 仓库在
`LiLi_Test_ws/src/LiLi_Test/`**，所有命令都要在这个子目录里跑。

```bash
cd src/LiLi_Test
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j4
ctest --output-on-failure          # 6 个用例，约 15 秒
```

工具：

| 工具 | 用途 |
| --- | --- |
| `lili_scenes [--noise s] [--spacing h]` | M1 场景验收：解析真值基是否保持 Q 在表面内 |
| `lili_eval [--noise s] [--corr d] [--q-samples n] [--local-scale k] [--local-phase p]` | M2：共享 ICP + Zhang 基线 vs 解析真值 |
| `lili_cli` / `lili_synthetic` / `lili_tests` | 原有工具，未改动行为 |

## 2. 代码清单（全部未提交）

新增库文件：`include/se3.h` + `src/se3.cpp`、`include/point_cloud.h` +
`src/point_cloud.cpp`、`include/scenes.h` + `src/scenes.cpp`、`include/icp.h` +
`src/icp.cpp`、`include/detectors.h` + `src/zhang_detector.cpp`、
`include/eval_quality.h` + `src/eval_quality.cpp`。

新增工具与测试：`src/scene_report.cpp`（`lili_scenes`）、`src/icp_eval.cpp`
（`lili_eval`）、`tests/test_se3.cpp`、`tests/test_quality_metric.cpp`、
`tests/test_icp.cpp`；`CMakeLists.txt` 与 `README.md` 有改动。

**一个历史遗留修复**：`src/cli.cpp` 原本缺少 `#include <random>`，`lili_cli`
根本编译不过。已补上，否则整个构建都起不来。

`git status`：`docs/` 与上述新文件均未跟踪；`CMakeLists.txt`、`README.md`、
`src/cli.cpp` 为已修改。**没有任何提交**，需要时请先 `git add`。

## 3. M1 已完成

SE(3) 的 `exp/log/adjoint/twistAtPoint`（含小角度级数与近 π 分支）、六个合成
场景、均匀格最近邻索引、Eq. 8–10 的 $Q$ 指标与子空间主角。

验收（`lili_scenes`）：无噪声与 $\sigma=0.02$ 两种条件下全部 PASS。无噪声时
$Q(P_{gt})$ 精确为 0，真值基的 $Q$ 远小于对照方向。

两处真值修正：**正弦圆柱的精确对称只有绕轴旋转**（$r(z)$ 与 $\theta$ 无关，
轴向平移不精确）；**矩形扭转隧道换成了螺旋对称的非直纹管**（见第 4 节）。

## 4. 必须保留的关键技术决策

这几条都是踩过坑之后定的，推翻前请先确认理由。

1. **body 系右乘约定**。ICP 增量 $T \leftarrow T\exp(\delta^\wedge)$，
   $J=[\,-n^\top R[p]_\times,\ n^\top R\,]$，于是 $H=\sum J^\top J$ 直接落在
   body 系，与 LiLi 基（论文 Eq. 7/9）同系，**不需要伴随变换**。遗留的
   `InformationMatrixCalculator::addPointPlane` 用的是 world 系平移增量，两者
   不可互换。
2. **ICP 必须用 LM 阻尼**。裸 Gauss-Newton 在退化场景下沿零空间产生无界步长：
   实测平面场景旋转 76°、圆柱沿轴滑 3.2 m。
3. **验收判据不能是"恢复到 $P_{gt}$"**。退化场景解不唯一。正确判据是误差 twist
   $b=\log(P_{gt}^{-1}P_{opt})$ 落在已知退化子空间内。无噪声实测子空间外占比：
   plane $1.6\times10^{-15}$、两圆柱 $1.3\times10^{-13}$、sinusoidal
   $3.9\times10^{-8}$、screw_tube $0.045$。噪声下该性质不成立（噪声激励所有
   弱约束方向），只能要求 $|b|$ 有界。
4. **局部视图默认是参考点云的精确子集**（`local_spacing_scale=1`,
   `local_phase=0`），因此无噪声时 $Q(P_{gt})=0$。独立采样通过
   `--local-scale` / `--local-phase` 开启。抽稀必须在派生两个视图**之前**做，
   否则子集关系被破坏（曾因此让 210/420 个局部点失去精确对应）。
5. **$Q$ 的搜索半径有上限**（默认 8 倍中位点间距）并在超限处删失。原因是方向
   严重错误时最近邻搜索会膨胀到极端环数（一度把一次评估拖到 10 分钟以上）。
6. **格索引 cell 取搜索半径的一半**。此前用"接近点间距"的 cell 让测试套件耗时
   266 秒；改后 15 秒（约 19 倍）。

## 5. 今天最后定位的问题（当晚第一件要做的事）

### 观察到的事实

把局部视图改成独立采样（`--local-scale 2 --local-phase 0.5`，$\sigma=0.02$）
后，Zhang 的**方向**误差确实变大了，但 $Q$ 几乎不变：

| 场景 | 真值维数 | Zhang 维数 | 最大主角 | $Q_{true}$ | $Q_{zhang}$ |
| --- | --- | --- | --- | --- | --- |
| sinusoidal_cylinder | 1 | 1 | 20.19° | 0.03153 | 0.03033 |
| open_cylinder | 2 | 2 | 3.26° | 0.04010 | 0.03920 |
| screw_tube | 1 | 2 | 10.35° | — | — |

正弦圆柱是 1 维、比较是适定的（2 维场景内部基的选取任意，不能逐列比较）。它的
真值基与 Zhang 基分别为：

```
truth b = (0, 0, 1, 0, 0.55, 0)                       |b| = 1.1413
zhang b = (-0.0327, -0.0626, -0.855, -0.0583, -0.4538, -0.2331)  |b| = 1.0000
```

### 根因

把真值基人为旋转 5°/20°/45° 后 $Q$ 为 0.0324 / 0.0380 / 0.0488（真值 0.0315），
说明 $Q$ 对角度**是**单调敏感的。但 Zhang 实际基给出 0.0303，比真值还低——
因为它被评为 1.0 的范数、而真值基是 1.1413：**幅度小了 14%，扩展位移整体变小，
把方向错误带来的惩罚抵消掉了**。

结论：**$Q$ 同时混合了"方向是否正确"与"基向量幅度"两个因素**，不同方法产生的
基范数不同时，比较不是同一把尺子。论文 Eq. (9) 直接使用基向量乘系数，因此这个
陷阱在论文本身也存在。

### 当晚第一步

在计算 $Q$ 之前把两个基都归一化（在按特征长度 $L$ 缩放旋转行的同一空间里做单位
化），或改为固定扩展位移幅度后再比较。重跑后**再判断 Zhang 是否真的失效**。
如果归一化后 Zhang 仍然只是 20° 量级的方向误差而 $Q$ 差异很小，说明需要换更有
区分度的指标（例如直接报告子空间主角 + 扩展扫描的 $Q$，两者并列，不指望单一
标量）。

## 6. 后续步骤

0. **先跑一次 `ctest`。** 今天最后一次完整 `ctest`（6/6 通过）是在 $Q$ 删失与
   `local_spacing_scale` / `local_phase` 改动**之前**跑的；此后只单独跑过
   `lili_eval` 与临时调试程序，`lili_scenes` 与 `lili_quality_tests` 未复跑。
   预计仍会通过（所有"在表面附近"的距离都远小于 8 倍点间距的删失上限），但请
   先验证再继续。
1. 归一化 $Q$ 比较（见第 5 节），重新评估 Zhang 在哪些场景/参数下失效。
2. 若仍未失效，候选杠杆（按优先级）：局部覆盖参考的一小部分；更大的 $P_{init}$
   偏移；点云密度差更大；噪声更大。
3. 门槛通过后再实现 LiLi：$k$ 自适应扰动标定、逐扰动重优化、
   $T=P_{opt}^{-1}P_{perturbed}$ 与 $\log(T)$、$\tau_{displacement}$ 筛选、
   PCA 子空间、$\ell_1$ 稀疏化。
4. M4：重关联消融（对应假设 H3，价值最高）、$k$ 扫描、噪声扫描、阈值敏感性。

## 7. 未决问题与风险

- **论文的真实隧道数据集未公开**，Trial 1/2 的 ATE/RPE 无法复现，本验证仅覆盖
  合成部分。
- 当前所有结论都建立在"解析退化子空间"正确之上；正弦圆柱与螺旋管的真值已推导
  并数值确认，但如果后续改场景，真值必须重新推导。
- `Q` 的搜索半径上限（8 倍点间距）会删失大的距离。当前所有"在表面附近"的数值
  不受影响，但若将来需要与论文的绝对 $Q$ 数值对比，需要重新确认该上限是否足够。
- `screw_tube` 场景在无噪声下的子空间外占比为 $0.045$（判据阈值 0.10），余量
  不大，来自曲面点面残差的一阶近似误差；若调参后超标，需要先改进 ICP 收敛。

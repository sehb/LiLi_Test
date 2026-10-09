# LiLi 算法有效性验证方案

验证对象：`docs/2609.17145v2.pdf`（*LiLi: Lie Theory Based 3D LiDAR Scan
Alignment Degeneracy Detection*）提出的扰动式退化检测方法。

范围：在合成数据上独立复现论文的核心结论（Table I）。真实隧道数据的
ATE/RPE 复现暂不纳入，原因见"范围与风险"。

## 1. 验证目标：把论文主张拆成可证伪假设

| 编号 | 论文主张 | 可证伪形式 |
| --- | --- | --- |
| H1 | 扰动法能识别完整退化子空间（含耦合的旋转 + 平移） | Open Cylinder 等场景下，LiLi 恢复的子空间与解析真值的最大主角 ≤ 5°，而 Hessian 法 > 20° |
| H2 | 噪声下优势显著，误差降低约 50% | 加噪条件下 `Q_LiLi ≤ 0.5 × Q_Zhang`（论文 Table I 开放圆柱为 0.037 vs 0.382） |
| H3 | 优势来源是显式数据重关联 | 冻结对应关系的消融使 LiLi 的优势至少减半 |
| H4 | 不产生假阳性 | Random3D 等非退化场景，两法都报 0 维退化子空间 |

H3 直接检验论文的机制解释，成本低、价值最高，必须实现。

## 2. 实验设计

场景沿用论文 Fig. 2 的五种几何：Plane、Closed Cylinder、Open Cylinder、
Sinusoidal Cylinder、Rotating Tunnel，并刻意让几何轴与坐标轴不重合。

每种场景两种条件：无噪声；每点每轴加 $N(0, 0.02)$ 高斯噪声。点云平均
最近邻间距取 0.05，与论文一致。

公平性要求（关键）：

- 两种方法共用同一个配准器和同一个最终位姿。
- 扩展扫描（Eq. 8）使用同一组采样系数、相同种子、相同自由度 $k$、相同
  样本数 $n$。

满足以上两点后，$Q$ 的差异才只来自检测到的子空间本身。

规模：每格至少 10 个随机种子，报告均值 ± 标准差，最终产出
5 场景 × 2 噪声条件 × 2 方法的对照表，直接对标论文 Table I。

## 3. 评价指标

1. **$Q$（论文 Eq. 10）**：扩展扫描点到参考扫描最近邻距离的中位数。
   与论文可比的主指标。
2. **相对改善**：`(Q_Zhang - Q_LiLi) / Q_Zhang`，用于量化 H2。
3. **子空间主角**：检测基正交化后对 $B_{true}^\top B_{est}$ 做 SVD，
   $\theta_i = \arccos \sigma_i$，报告最大角与维度误差 $|d_{est} - d_{true}|$。
   论文只报 $Q$，加入主角才能分离"方向对不对"与"扩展扫描质量好不好"。

真值子空间需解析推导：平面为面内两平移 + 绕法向旋转（3 维）；圆柱为轴向
平移 + 绕轴旋转（2 维，论文明确给出）。正弦圆柱与旋转隧道的真值需推导并
用高密度无噪声数据数值确认，不可照搬。

注意量纲：$se(3)$ 基向量中旋转分量为 rad、平移分量为 m，直接计算主角会
因量纲差异失真。引入特征长度 $L$（扫描点到中心的中位距离），把旋转分量
乘以 $L$ 后再比较，并在报告中说明该选择。

## 4. 需要新建的模块

| 文件 | 职责 |
| --- | --- |
| `include/se3.h`, `src/se3.cpp` | `hat`、`exp`、`log`、左右雅可比；当前仓库无李群工具 |
| `include/point_cloud.h`, `src/point_cloud.cpp` | 点云/法向结构、最近邻搜索（点数小，暴力即可） |
| `include/icp.h`, `src/icp.cpp` | 点到面 Gauss-Newton ICP，每次迭代重关联；所有验证的前置条件 |
| `src/scenes.cpp` | 五个场景生成器（含轴错位与噪声注入），替代现有 `synthetic.cpp` |
| `include/degeneracy_detector.h` | 统一接口 `detect(...) -> Basis` |
| `src/zhang_detector.cpp` | 复用现有 `analyzeDegeneracy` |
| `src/lili_detector.cpp` | $k$ 自适应扰动标定、逐扰动重优化、$T = P_{opt}^{-1} P_{perturbed}$、`log(T)`、$\tau_{displacement}$ 筛选、PCA、$\ell_1$ 稀疏化 |
| `src/eval_quality.cpp` | 质量指标，Eq. 8–10 |
| `src/experiment_main.cpp` | 运行实验矩阵并输出 CSV |
| `tests/test_se3.cpp`, `tests/test_icp.cpp`, `tests/test_quality_metric.cpp`, `tests/test_subspace.cpp` | 单元测试 |

$\ell_1$ 稀疏化（论文 Eq. 14）是非凸问题，建议用迭代加权或交替优化，并同时
报告稀疏化前/后的结果。稀疏化理论上不改变子空间，这一不变量本身可被验证。

## 5. 消融与鲁棒性扫描

- **重关联消融（对应 H3）**：保持相同扰动与重优化，只锁定初始对应关系。
- 稀疏化开关。
- 扰动尺度 $k$ 扫描：1、3、5、10。
- 噪声扫描 $\sigma \in [0, 0.05]$，绘制 $Q$–噪声曲线，验证 Hessian 法是否
  如论文所述在复杂场景下崩溃。
- 阈值 $\tau_{PCA}$、$\tau_{displacement}$ 敏感性分析。使用留出种子，避免在
  测试集上调参。

## 6. 判定标准

H1、H2（复现论文量级）、H3（优势随消融减半以上）、H4（零假阳性）四条同时
成立，即判定论文算法有效性得到独立验证。若 H2 只在部分场景成立，需明确报告
哪些场景不成立——这是有效结论，且比笼统宣称"有效"更有价值。

## 7. 实施顺序

1. **M1**：SE(3) 工具 + 场景生成 + 最近邻 + 质量指标 $Q$；用手工构造的
   已知基验证 $Q$。
2. **M2**：共享 ICP 配准器 + Zhang 检测器，先复现 baseline 的 Table I 行。
3. **M3**：LiLi 检测器（扰动标定、重优化、PCA、稀疏化），产出主对照表。
4. **M4**：消融与鲁棒性扫描，输出 CSV 与曲线。

M1 是其余一切的前提，未通过 M1 不应进入 M3。

### 实施状态

**M1 已完成并接入 CTest**。新增文件：`include/se3.h` + `src/se3.cpp`（exp/log/
adjoint/twistAtPoint）、`include/point_cloud.h` + `src/point_cloud.cpp`（均匀
格最近邻索引、最近邻间距诊断）、`include/scenes.h` + `src/scenes.cpp`（五个
场景 + 解析真值基）、`include/eval_quality.h` + `src/eval_quality.cpp`（Eq. 8–10
的 $Q$ 与主角）。验收工具为 `lili_scenes`（`tests/` 之外，注册为 CTest 用例
`lili_scene_acceptance`）。

关键实现约定（与论文的对应关系）：

- 扰动为左乘 $P_{perturbed} = \exp(\Delta^\wedge) P_{opt}$，与 `info_matrix.h`
  中 $J = [\,n^\top(-R[p]_\times),\ n^\top\,]$ 的切空间次序 $[\varphi;\rho]$
  一致。
- 退化基表达在**传感器（body）坐标系**，扩展扫描按式 (9)
  $T_i = P_{opt}\exp(\sum_j c_{ij}b_j^\wedge)$ 右乘，因此
  $t_{degeneracy} = \log(P_{opt}^{-1}P_{perturbed})$ 给出的正是 body 系 twist。
- 系数 $c_{ij}$ 取自确定性 Halton 序列（而非随机数），保证两种方法在基维数
  不同时仍收到完全相同的系数，比较是公平的。
- 局部视图是参考点集的**子集**（按形状参数窗口选取，且抽稀在派生之前完成），
  因此无噪声时 $Q(P_{gt})$ 精确为 0。

验收结果（`./lili_scenes`，两项判据：真实基不使 $Q$ 超出
$\max(Q_{plain}, h/2)$ 的 1.5 倍；对照方向使 $Q$ 至少变差 1.5 倍）：

| 场景 | 真值维数 | 无噪声 | 噪声 $\sigma=0.02$ |
| --- | --- | --- | --- |
| plane | 3 | PASS | PASS |
| closed_cylinder | 2 | PASS | PASS |
| open_cylinder | 2 | PASS | PASS |
| sinusoidal_cylinder | 1 | PASS | PASS |
| rotating_tunnel | 1 | PASS | PASS |

无噪声典型值：$Q_{plain} \approx 0$，$Q_{true} = 0.013\text{–}0.027$，
$Q_{control} = 0.059\text{–}0.111$。参考点云平均最近邻间距约 0.047–0.05，
与论文设定一致。

真值推导的两处修正（相对最初方案）：

- 正弦圆柱的精确对称只有**绕轴旋转**（$r(z)$ 与 $\theta$ 无关），轴向平移
  并不精确，故为 1 维而非 2 维。
- 旋转隧道用**非圆截面**（矩形）沿轴扭转才真正构成螺旋对称；圆截面扭转
  等价于直圆柱，不产生耦合退化。真值为单条螺旋 twist
  $b=(0,0,\alpha,0,0,1)$。

**M2 已完成并接入 CTest**。新增 `include/icp.h` + `src/icp.cpp`（点到面 ICP）、
`include/detectors.h` + `src/zhang_detector.cpp`（Zhang 基线）、
`src/icp_eval.cpp`（驱动工具 `lili_eval`，CTest 用例 `lili_m2_baseline`）。

ICP 的设计要点（都影响结论，必须记录）：

- 增量在 body 系右乘 $T \leftarrow T\exp(\delta^\wedge)$，Jacobian 为
  $J=[\,-n^\top R[p]_\times,\ n^\top R\,]$，因此 $H=\sum J^\top J$ 直接落在
  body 系，与 LiLi 基的约定一致，**无需伴随变换**。注意这与遗留的
  `InformationMatrixCalculator::addPointPlane`（平移块为 $n^\top$，world 系平移
  增量）不是同一约定，二者不可互换。
- 用 Levenberg–Marquardt 阻尼而非裸 Gauss–Newton。退化场景下 $H$ 奇异，裸 GN
  会沿零空间产生无界步长（实测平面场景旋转 76°、圆柱平移 3.2 m）。

验收判据（与 M1 的区别很重要）：退化场景下配准解**本来就不唯一**，所以不能要求
$P_{opt} \approx P_{gt}$。正确判据是误差 twist
$b=\log(P_{gt}^{-1}P_{opt})$ 是否落在已知退化子空间内：

- 无噪声：退化场景要求子空间外分量占比 < 10%；非退化对照要求 $|b| < 0.02$。
- 有噪声：噪声会激励所有弱约束方向，$b$ 不再局限于解析子空间，只能要求
  $|b| < 0.25$ 有界。

实测（无噪声）子空间外占比：plane $1.6\times10^{-15}$、两个圆柱
$1.3\times10^{-13}$、sinusoidal $3.9\times10^{-8}$、screw_tube $0.045$、
random_surface 满足 $|b|<7\times10^{-15}$。screw_tube 的 $0.045$ 来自曲面点面
残差的一阶近似误差，不是子空间错误。

### M2 基线结果（Zhang）

| 场景 | 真值维数 | 无噪声 Zhang 维数 | 有噪声维数 | 最大主角（无噪声 / $\sigma=0.02$） |
| --- | --- | --- | --- | --- |
| plane | 3 | 3 | 3 | 0.00° / 0.30° |
| closed_cylinder | 2 | 2 | 2 | 0.00° / 1.03° |
| open_cylinder | 2 | 2 | 2 | 0.00° / 1.98° |
| sinusoidal_cylinder | 1 | 1 | 1 | 0.00° / 6.94° |
| screw_tube | 1 | 2（多报一个弱方向） | 2 | 0.78° / 8.08° |
| random_surface | 0 | 0 | 0 | 0.00° / 0.00° |

### 关键负面发现（必须在 M3 前解决）

**当前设置没有复现论文所述的 Zhang 失效。** 论文中 Zhang 在 Open Cylinder 加噪后
完全误判退化空间，而这里它的方向误差仅 1.98°，$Q$ 与真值几乎相同
（0.0393 vs 0.0395）。因此不能说"已复现论文 Table I"。

差异的可能来源（M3 前需要逐一排查）：

- 论文的局部扫描是真实扫描——**不是参考点云的子集**，重关联的歧义远大于本设置
  （这里 $P_{gt}$ 使每个局部点都有精确对应）。
- 论文用 $k$ 自适应标定的扰动幅值；当前 $P_{init}$ 偏移很小（$\approx 0.04$ m），
  配准误差的 basin 很窄。
- 论文的噪声是每点独立 $N(0,0.02)$ 且**点云更稀疏**，最近邻对应更容易跳变。

在解决这一点之前，M3 的 LiLi 与 Zhang 对比可能给出"两者都正确"的平凡结论，无法
检验论文的核心主张。建议 M3 的第一步是构造一个让 Zhang 真正失效的设置（例如局部
扫描只覆盖参考的一小部分、或用不同的采样密度生成局部与参考），并把它作为 M3 的
前置门槛。

**M3 前置工作进行中**，日常进度与入口见 `docs/handover.md`。当前状态：为了让
Zhang 真正失效，已加入"局部视图独立采样"选项（`SceneOptions::local_spacing_scale`
与 `local_phase`，对应 `lili_eval --local-scale/--local-phase`）。独立采样确实
使 Zhang 的方向误差增大（正弦圆柱从 6.9° 增到 20.2°），但 $Q$ 几乎没有变化。

**这暴露了一个指标缺陷：$Q$ 混合了方向正确性与基向量幅度。** 正弦圆柱的 Zhang
基范数为 1.0、真值基为 1.1413，幅度小 14% 抵消了方向错误的惩罚。把真值基人为
旋转 5°/20°/45° 得到 $Q$ = 0.0324/0.0380/0.0488（真值 0.0315），说明 $Q$ 本身
对角度敏感，问题出在幅度未归一化。论文 Eq. (9) 直接使用基向量乘系数，同样存在
这个陷阱。

因此在 M3 之前必须先归一化基向量（或在按特征长度缩放旋转行的同一空间内单位
化）再比较 $Q$，然后重新判断 Zhang 是否失效。`Q` 的搜索半径现已有上限（默认
8 倍中位点间距）并在超限处删失，否则方向严重错误时最近邻搜索会膨胀到极端环数。

尚未完成：M3（LiLi 检测器）、M4（消融与扫描）。

## 8. 范围与风险

论文的 Tunnel 数据集因双盲评审未公开，Trial 1/2 的 ATE/RPE 无法复现。若需
实车验证，替代方案是自采长走廊/隧道数据，或使用公开的长走廊序列（如 DARPA
SubT 类数据集，需先确认可用性）。建议先只做合成部分，把真实场景明确标为
未覆盖。

主要风险：

- 真值子空间推导错误会使全部结论失效，需独立复核。
- ICP 陷入局部极小会污染 $Q$，需固定初值并检查收敛。
- $Q$ 对点云密度与覆盖范围敏感，两法必须使用完全相同的参考点云与采样。

以上三点需在实现前用单元测试钉死。

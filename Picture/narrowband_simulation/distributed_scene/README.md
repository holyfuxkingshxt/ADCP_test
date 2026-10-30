# 连续距离分布散射：先生成整段回波，再只处理一个 gate

**历史版本／不作为默认验收输入（2026-10-30）。** 本目录的 `distributed_narrowband_simulation.mat` 是连续散射群的单次随机场景，保留用于分析随机误差。后续算法开发与单次精度验收使用 [`../algorithm_validation/narrowband_validation.mat`](../algorithm_validation/narrowband_validation.mat)，运行同目录的 `run_validation.m`；连续散射的统计验收使用 `continuous_ensemble_validation.mat`。

**依赖保留：** 当前验证数据生成器仍引用本目录的 `distributed_demo.py`、`add_tgc.py`、`distributed_signals.npz` 与 `distributed_results.json`，因此不要将本目录整体删除或移动。单次结果偏离真值不能直接判定算法错误；该场景需要按统计目标评估。

这是原单点示例的扩展。所有距离区的实数 RF 回波先叠加，再统一做 I/Q 解调、低通与抽取。测速只使用一个3 ms接收窗，但生成接收信号时没有把其他距离区裁掉。

## 当前场景与文件

- 单声束、固定 ADCP、声速1500 m/s、载频300 kHz，远离为正。
- 发射无编码的4 ms长正弦脉冲；相关lag仍为0.5 ms。
- 初始距离3至21 m，分为3–6、6–9、9–12、12–15、15–18、18–21 m六区，每区240个随机散射体，共1440个。
- 默认所有散射体以+0.2 m/s运动。另有把12–15 m区速度改为+0.6 m/s的对照。
- 完整接收记录0至34 ms：RF为81600点，I/Q为2040点。
- 只处理14.5至17.5 ms这个gate：180个I/Q点，lag30点，150对乘积。
- 接收幅度采用随机复数散射幅度，并加入假设的双程点目标几何衰减和吸收。各区每米点数相同，不是完整三维声束体积或绝对声压模型。

主数据：`distributed_narrowband_simulation.mat`。图册：`distributed_scene_figures.pdf`。运行 `inspect_distributed.m` 可在 MATLAB 画完整接收记录、各区在所选 gate 内的贡献，并重新计算该 gate 的速度。

已加入接收时间增益补偿（TGC，也常称TVG），对比图为 `tgc_comparison.pdf`。MAT中 `rx_rf_noisy` 仍表示未补偿的完整含噪回波，`rx_rf_noisy_tgc` 表示补偿后的完整含噪回波；两者都是81600×1。查看脚本同时绘出两条处理路径。

## 时间增益补偿如何计算？

当前点散射体的声压衰减假设为 `(10/r)^2 exp(-2 alpha r)`，其中距离r以米计，单程声压吸收系数alpha=0.005 Np/m。补偿取该幅度衰减的反比，以3 m为增益1的参考点。

由于4 ms发射脉冲中的每个时刻都会产生回波，一次接收时刻不能唯一对应一个距离。这里用发射脉冲中心作为代表：

$$
r_{\mathrm{rep}}(t)=\operatorname{clip}\!\left[\frac{c(t-T_{\mathrm{tx}}/2)}2,\ 3\ \mathrm m,\ 21\ \mathrm m\right].
$$

时间t从发射开始计算，单位秒；c=1500 m/s，Ttx=0.004 s。clip表示把距离限制在给定范围内。代表距离达到3 m的时刻为6 ms，达到21 m的时刻为30 ms。

$$
g(t)=\min\!\left\{10^{40/20},\ 
\left[\frac{r_{\mathrm{rep}}(t)}{3\ \mathrm m}\right]^2
\exp\!\left[2\alpha\bigl(r_{\mathrm{rep}}(t)-3\ \mathrm m\bigr)\right]\right\}.
$$

这是正实数幅度增益，增益分贝为20 log10(g)。默认最高允许40 dB；本场景实际到21 m时为58.6637倍、35.3674 dB，之后保持这个值，避免在记录末尾继续无限放大。此曲线是与当前简化衰减模型对应的教学选择，不是某型号ADCP的厂家默认曲线。

处理顺序为：

```text
全部距离区的真实RF叠加 → 加接收噪声 → 逐样本乘TGC增益
                                         ↓
                              I/Q混频 → FIR低通 → 抽取
                                         ↓
                              取同一个180点gate → 150对相关
```

噪声与回波一起增益：`rx_rf_noisy_tgc = tgc_gain_rf .* rx_rf_noisy`。这是采样后的数字增益演示，没有模拟模拟接收放大器、ADC饱和或量化。它可以补偿幅度趋势，但不能恢复已被噪声掩盖的信息；每个RF样本的信号和噪声之比不因乘同一增益而改善。长脉冲混合多个距离，随机相干叠加也仍存在，因此补偿后包络不会成为平直直线。

正实数增益不直接旋转信号相位，但不同时间的增益不同，会改变滤波和各对相关乘积的权重。本次固定散射配置的无噪声单次速度从0.027047变为0.078752 m/s，含噪声时从0.027079变为0.078696 m/s；真实速度仍为0.2 m/s。这个变化不能解释为已经消除单次随机误差，也不能用来断言精度提高。

检查已通过：RF除以增益后恢复原数据，最大误差2.8e-17；每个RF样本的信号与噪声增益之比与1的最大差4.5e-16；补偿后各区I/Q相加与总I/Q的最大差2.9e-16；40 dB上限、增益单调性及全部MAT数组读回一致性也已检查。MATLAB查看脚本尚未在MATLAB中执行。

`distributed_demo.py` 生成完整模拟，`export_distributed_mat.py` 将结果导出MAT并重新加载、逐数组验证。`distributed_results.json` 保存统计和检查结果，`range_zone_contributions.csv` 列出各区进入gate的散射体数与独立区场能量。

## 为什么一个gate含有不同距离的回波？

暂忽略散射体在传播期间的微小移动。声波在发射脉冲内的时刻s出发，从距离r返回的时刻是：

$$
t=s+\frac{2r}{c},\qquad 0\leq s\leq T_{\mathrm{tx}}.
$$

接收门为T1至T2。只要一个散射体的回波与门有重叠，就能贡献信号。因此，门内可能出现的直达散射距离为：

$$
\frac{c(T_1-T_{\mathrm{tx}})}2\leq r\leq\frac{cT_2}2.
$$

代入本例T1=14.5 ms、T2=17.5 ms、Ttx=4 ms，得到约7.875至13.125 m。用代码中的精确匀速传播关系，以发射开始时的初始距离计，得到7.87315至13.12325 m。

因此，9–12 m只是为比较而划出的物理区，不能等同于这个gate的全部空间采样体积。6–9 m区的较远部分、12–15 m区的较近部分也会贡献回波。

FIR低通会额外扩展时间响应。本例半长度为0.333 ms，按距离约为0.25 m；图和MAT里的几何gate/pair权重是未计滤波展宽的值。代码离线补偿了对称FIR的群延迟，故其时间戳对应的有效采样区间在两侧都有短暂扩展。这不是来自任意远处的回波。

## 幅度为什么既能升高也能降低？

不同区域回波必须先带相位相加：

$$
x_{\mathrm{total}}(t)=\sum_zx_z(t).
$$

总功率含有交叉项：

$$
|x_{\mathrm{total}}|^2=\sum_z|x_z|^2+
2\sum_{z<k}\operatorname{Re}\{x_zx_k^*\}.
$$

交叉项可正可负。因此附近区域能在某些时刻增强幅值，在另一些时刻相消。图中的“incoherent RSS”仅为去掉交叉项后的比较量，不是实际接收波形。

本次固定随机种子的实现中，与移除12–15 m区相比，加入它后约1/3的gate采样点幅值升高，其他点幅值降低。这个比例只是一次散射配置，不是通用规律。CSV中各区单独场的能量相加，也不必等于相干总场的能量。

## 更远区为何没有抬高这个gate？

18–21 m区的最早直达回波约24 ms才到，而所选门在17.5 ms结束，即使考虑短FIR展宽也不相交。代码将该区的回波幅值放大10倍，重新对完整RF解调后，所选gate内的最大I/Q差约3.6e-17，处于浮点误差级别。

所以不是所有更深处回波都会干扰前面的时间窗。实际多途、前一ping的迟到回波和边界强回波等需要其他模型；本例没有加入这些效应。

## 相同真实流速，为何单次速度仍不等于0.2？

本次均匀流速的单次估计为约0.02705 m/s。它是完整散射群的一个随机实现：有限脉冲扫过时，各颗粒复数回波的加入和退出会导致相位起伏。所选窗中180点并不是180个独立信息源；相干系数约0.907也不能保证该次估计准确。

为检查程序而不是挑选一个恰好接近真值的随机种子，另外做了两项验证：

1. 构造从散射复幅度到该gate I/Q的精确线性映射，包含同一个RF采样规则、FIR和实数信号的残余镜像。与逐散射体生成整段实数RF再处理的结果，最大差约6.0e-14。
2. 用该精确映射计算独立散射配置的统计期望：滤波后相关期望对应0.199963 m/s。1000组独立散射配置平均复数相关后得到约0.20336 m/s；单次速度估计标准差约0.23862 m/s。

这些配置用来检查随机误差和统计期望，没有模拟真实ping之间的水体演化。不能把1000组全部看作特定ping频率下可独立平均的数据。

## MAT数组与MATLAB索引

| 变量 | MATLAB维度 | 含义 |
|---|---|---|
| `t_rf_s` | 81600×1 | RF接收时间，秒，从发射开始计 |
| `tx_rf` | 81600×1 | 完整发射波形 |
| `rx_rf_clean`、`rx_rf_noisy` | 81600×1 | 所有距离区叠加后的接收波形 |
| `rx_rf_zones` | 81600×6 | 六个物理距离区各自的实数回波，按列排列 |
| `tgc_gain_rf`、`tgc_gain_db_rf` | 81600×1 | 接收时间增益，分别为幅度倍数与dB |
| `tgc_range_rf_m` | 81600×1 | 用于TGC曲线的代表距离，米，已限制在3–21 m |
| `rx_rf_clean_tgc`、`rx_rf_noisy_tgc`、`noise_rf_tgc` | 81600×1 | 完整记录的补偿后回波与噪声 |
| `iq_total`、`iq_noisy` | 2040×1 | 完整复数基带 |
| `iq_total_tgc`、`iq_noisy_tgc`、`iq_noise_tgc` | 2040×1 | 从补偿后的RF重新解调、滤波和抽取得到的基带 |
| `iq_zones_tgc`、`iq_gate_zones_tgc` | 2040×6、180×6 | 各物理区经过同一TGC与解调后的贡献 |
| `iq_gate_tgc`、`iq_gate_noisy_tgc` | 180×1 | 同一个gate内补偿后的无噪声和含噪基带 |
| `lag_products_tgc`、`lag_products_noisy_tgc` | 150×1 | 补偿后无噪声和含噪gate的滞后共轭乘积 |
| `iq_zones` | 2040×6 | 各距离区单独解调后的基带，可验证线性叠加 |
| `iq_gate_zones` | 180×6 | 所选gate中每区的贡献 |
| `iq_gate` | 180×1 | 所选gate的总复数信号 |
| `lag_products` | 150×1 | 相隔30点的共轭乘积 |
| `iq_far_boosted` | 2040×1 | 18–21 m区幅值×10后的基带 |
| `rx_rf_shear`、`iq_shear` | 81600×1、2040×1 | 将12–15 m区速度改为+0.6 m/s的对照 |
| `scatterer_*` | 1440×1 | 每个散射体的初始距离、区号、复幅度、速度和回波起止时间 |
| `ensemble_velocity_mps` | 1000×1 | 独立散射配置的单次速度估计 |
| `params`、`diagnostics`、`notes` | struct | 参数、结果与字段说明 |

`gate_indices_matlab`可直接索引MATLAB数组，为871至1050；`gate_indices_python`保留从0开始的870至1049。所有时间单位为秒；正速度为远离；区号从1开始。

MATLAB中查看叠加关系：

```matlab
S = load('distributed_narrowband_simulation.mat');
g = S.gate_indices_matlab;
total = sum(S.iq_zones,2);
max(abs(total-S.iq_total))
plot(S.t_gate_s*1000,abs(S.iq_zones(g,[2 3 4])));
hold on; plot(S.t_gate_s*1000,abs(S.iq_total(g)),'k','LineWidth',2);
legend('6-9 m','9-12 m','12-15 m','Total');
```

## 重新生成

```powershell
& 'C:\Users\pc\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' 'C:\Users\pc\Desktop\ADCP\Picture\narrowband_simulation\distributed_scene\distributed_demo.py'
& 'C:\Users\pc\AppData\Local\Programs\Python\Python314\python.exe' 'C:\Users\pc\Desktop\ADCP\Picture\narrowband_simulation\distributed_scene\export_distributed_mat.py'
```

第二个Python含NumPy和SciPy，负责导出MAT并逐数组重读验证。MATLAB查看脚本尚未在MATLAB中实际运行；保存数据已通过Python读回检查。

完整生成器已接入TGC。若已有旧版完整场景数据，只需运行下面两步即可增加或更新TGC，无须重新生成散射体或统计集合：

```powershell
& 'C:\Users\pc\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' 'C:\Users\pc\Desktop\ADCP\Picture\narrowband_simulation\distributed_scene\add_tgc.py'
& 'C:\Users\pc\AppData\Local\Programs\Python\Python314\python.exe' 'C:\Users\pc\Desktop\ADCP\Picture\narrowband_simulation\distributed_scene\export_distributed_mat.py'
```

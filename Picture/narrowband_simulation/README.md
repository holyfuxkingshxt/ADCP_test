# 窄带 ADCP 单声束模拟：从发射声波到相关测速

**当前唯一默认算法验收入口：[`algorithm_validation/narrowband_validation.mat`](algorithm_validation/narrowband_validation.mat)。** 查看和重算方式见 [`algorithm_validation/README.md`](algorithm_validation/README.md)，完整验证脚本为同目录的 `run_validation.m`。

**旧版状态（2026-10-30）：顶层 `narrowband_simulation.mat`、`signals.npz` 及 `distributed_scene/distributed_narrowband_simulation.mat` 已废弃为默认算法验收输入，仅保留为历史教学或随机误差分析资料。** 后续开发、数值对照均以 `algorithm_validation/` 为准。该目录中的边界测试和连续散射集合测试仍然有效。

当前主数据包含完整多层RF、噪声、TGC及解析真值。2026-10-30从保存的完整含噪TGC RF独立重算六层速度，最大绝对误差为0.000405807519 m/s，小于0.001 m/s的既定容差。它是一维束向、受控散射薄层基准，未设置声束相对垂直方向的倾角；真值不是直接给定的水平流速。完整验证记录见 `algorithm_validation/validation_report.json`。

**依赖保留：** 新版生成器仍导入顶层 `narrowband_demo.py`，以及 `distributed_scene/distributed_demo.py`、`add_tgc.py`，并读取 `distributed_signals.npz`、`distributed_results.json`。因此保留旧目录结构与文件，不移动或删除这些依赖；旧数据不作为默认验收输入，并不表示其随机散射现象是错误的。

本例已经实际运行。先学习单声束、一个匀速散射体、无噪声的最小模型，再看散射体群、接收噪声与相位模糊。它用于核对物理关系和处理算法，不代表任何具体厂商的固件实现。

**历史连续散射场**保留在 `distributed_scene/`：覆盖3–21 m的1440个散射体，先合成整段实数RF，再只取一个gate测速。该单次随机场景用于分析随机误差；其统计平均验证使用 `algorithm_validation/continuous_ensemble_validation.mat`。

## 文件与运行

- `narrowband_demo.py`：完整代码。依赖 NumPy 与 ReportLab，不依赖 SciPy、Matplotlib。
- `narrowband_figures.pdf`：三页可放大的波形图。
- `figure-1.png`、`figure-2.png`、`figure-3.png`：对应的页面预览。
- `iq_gate.csv`：180 个用于测速的复数样本，列出 I、Q、时间和相位。
- `lag_products.csv`：150 对样本的复数乘积，可逐行核算。
- `signals.npz`：发射 RF、接收 RF、完整 I/Q、测速窗口和乘积数组。
- `narrowband_simulation.mat`：同一组信号及参数的 MATLAB 文件，所有一维信号均保存为列向量。
- `results.json`：参数、估计结果、验证和蒙特卡洛统计。
- `cloud_lag_comparison.csv`：固定接收窗口下，散射体群的 lag 对比。

在当前电脑的 PowerShell 中运行：

```powershell
& 'C:\Users\pc\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' 'C:\Users\pc\Desktop\ADCP\Picture\narrowband_simulation\narrowband_demo.py'
```

脚本会更新本目录中的模拟输出；其他笔记和文件不会被修改。预览 PNG 是用 Poppler 从图册导出的，脚本本身只生成 PDF、CSV、JSON 和 NPZ。

MAT 文件由 `export_mat.py` 从 NPZ 与 JSON 导出，并用 SciPy 重新读取逐数组验证。重新运行模拟后，可使用安装了 NumPy、SciPy 的 Python 再运行导出器：

```powershell
& 'C:\Users\pc\AppData\Local\Programs\Python\Python314\python.exe' 'C:\Users\pc\Desktop\ADCP\Picture\narrowband_simulation\export_mat.py'
```

MATLAB 读取方法：

```matlab
S = load('C:\Users\pc\Desktop\ADCP\Picture\narrowband_simulation\narrowband_simulation.mat');
plot(S.t_gate_s * 1000, real(S.iq_gate));
xlabel('Time (ms)'); ylabel('I');
q = S.iq_gate(31:180) .* conj(S.iq_gate(1:150));
phi = angle(mean(q));
fd = phi / (2*pi*S.estimate.lag_seconds);
v = -S.params.c * fd / (2*S.params.f0);
```

主要变量：`tx_rf`/`rx_rf` 为无噪声实数发射/接收波形，`rx_rf_noisy` 为含噪回波，`noise_rf` 为噪声本身，`iq_full` 为完整复数基带，`iq_gate` 为180点测速窗，`lag_products` 为150个复数乘积。`params`、`truth`、`estimate` 分别保存参数、理论值和估计值。`gate_indices_matlab`、`pair_indices_matlab` 已转换为 MATLAB 从1开始的索引；`estimate` 中原有的索引字段保留 Python 规则并在 `notes` 中注明。MAT 文件采用压缩 Level-5 格式，可由常见 MATLAB 版本读取；未实际启动 MATLAB。

## 1. 固定坐标、正负号与参数

ADCP 固定，散射体只沿一条声束运动。定义正速度为远离 ADCP，因此远离时多普勒频移为负。这与某些仪器输出的束向正方向可能不同，换方向时必须同步改符号。

| 量 | 符号/代码 | 示例设定 |
|---|---|---:|
| 声速 | `c` | 1500 m/s |
| 载频 | `f0` | 300000 Hz |
| 发射长脉冲持续时间 | `pulse_duration` | 0.004 s = 4 ms |
| 发射开始时的距离 | `initial_range` | 3 m |
| 散射体远离速度 | `velocity_away` | +0.2 m/s |
| 实数 RF 采样率 | `fs_rf` | 2400000 Hz |
| 反射幅度系数 | `amplitude` | 0.4，无量纲 |
| 散射初相位 | `scattering_phase` | 0.7 rad |
| 低通截止频率 | `lowpass_hz` | 6000 Hz |
| 抽取倍数 | `decimation` | 40 |
| 抽取后 I/Q 采样率 | `fs_rf / decimation` | 60000 Hz |
| 测速窗口点数 | `gate_samples` | 180 点，即名义 3 ms |
| 相关 lag | `lag_seconds` | 0.0005 s = 0.5 ms |

这些数字都是教学假设。忽略传播损耗、换能器指向性、硬件振铃、多途和流场变化。一个点散射体不能体现真实水体散射群的全部随机误差，因此脚本另外加入了散射群实验。

## 2. 发射什么波形？

设 `s` 为声波离开换能器的时刻。定义矩形包络 `p(s)`：当 `0 <= s < Tp` 时为 1，其余为 0。

$$
s_{\mathrm{tx}}(s)=p(s)\cos(2\pi f_0s).
$$

本例是无编码的单频长脉冲。4 ms 内有 1200 个载波周期。实数采样率为载频的 8 倍，每个载波周期有 8 个样本。代码的 `generate_rf()` 负责生成它。

## 3. 从往返传播推导接收回波，而非直接指定频移

设 `h` 是声波遇到颗粒的时刻，`t` 是回到 ADCP 的时刻。颗粒位置为：

$$
R(h)=R_0+vh.
$$

去程距离等于声速乘去程时间：

$$
c(h-s)=R_0+vh,
\qquad h=\frac{R_0+cs}{c-v}.
$$

回程为 `t = h + R(h)/c`，所以：

$$
t=\frac{2R_0}{c-v}+\frac{c+v}{c-v}s.
$$

它表示：第一个波峰经过往返传播后晚到，后续波峰的间隔还被乘上 `(c+v)/(c-v)`。这正是经典声学的时间拉伸。

把关系反解，得到接收时刻 `t` 对应的发射时刻：

$$
s=\alpha t-\beta,
\qquad \alpha=\frac{c-v}{c+v},
\qquad \beta=\frac{2R_0}{c+v}.
$$

于是实数接收信号是：

$$
r(t)=A\,p(\alpha t-\beta)
\cos\bigl(2\pi f_0(\alpha t-\beta)+\phi_s\bigr).
$$

在代码中，先计算 `emit_time = alpha*t - beta`，再将它代入发射波形。这一步同时产生了到达延迟、包络拉伸和载频变化。

本例精确理论值为：

- 回波开始时间约 4.000533 ms。
- 回波持续约 4.001067 ms，比发射长约 1.067 us。
- 接收载频约 299920.010665 Hz。
- 多普勒频移约 -79.989335 Hz。

频移的精确式与低速近似分别是：

$$
f_d=f_0(\alpha-1)=\frac{-2f_0v}{c+v}
\approx-\frac{2f_0v}{c},\qquad |v|\ll c.
$$

去程和回程都受到散射体运动影响，故低速近似中有系数 2。

## 4. 为什么还要 I/Q 解调？

RF 回波仍以约 300 kHz 振荡，要估计的是其中约 -80 Hz 的小差值。用已知发射频率将它搬到基带：

$$
z(t)=2r(t)e^{-j2\pi f_0t}.
$$

将余弦写成正、负频率两项后，这一步产生：

- 一项以 `fd` 为频率缓慢转动，携带所需多普勒信息；
- 另一项在约 `-(f_receive + f0)`、即 -600 kHz 附近。

低通滤波去掉后一项，得到：

$$
x(t)=I(t)+jQ(t)\approx A e^{j(2\pi f_dt+\phi_c)}.
$$

本例使用 `I = LPF(2*r*cos(...))`，`Q = LPF(-2*r*sin(...))` 的约定。因此远离时，复数相位随时间减小。

滤波器为 1601 点 Blackman 窗 sinc FIR，截止 6 kHz。对称 FIR 的群延迟是 800 个 RF 样本，即约 0.333 ms。代码为离线教学补偿了这个已知延迟；实时实现应保留因果延迟并修正时间戳。滤波卷积先补零，不是循环卷积。

低通之后再每 40 点取 1 点，得到 60 kHz 的复数 I/Q。这样可避免抽取时把高频镜像折叠到基带。所选最大多普勒频移还必须落在该低通通带内。

## 5. 精确说清楚：本例取多少数据？

回波边缘存在滤波暂态。本例在回波开始后至少 0.5 ms 再选窗，截取连续 180 个 I/Q 点；末端也远离包络边缘。

实际选取完整 I/Q 数组中的 `iq[271:451]`，Python 右端不包含：

- 第一个样本在 4.516667 ms。
- 最后一个样本在 7.500000 ms。
- 按采样窗计算长度是 `180 / 60000 = 3 ms`；首末采样时刻之差是 `179 / 60000`。

这个窗口仅用于隔离单个点目标的稳定回波内部，不能直接当成真实分布散射水体某个厂商的深度单元配置。真实 range gate 的权重还受有限发射脉冲影响。

## 6. 相关 lag 转成索引

I/Q 采样间隔是 `1/60000 s`，所以 lag 为：

$$
L=T_Lf_{s,IQ}=0.0005\times60000=30\ \text{个样本}.
$$

设局部窗口叫 `x`，长度为 180。本例只保留早、晚两个样本都在窗口内的乘积：

```python
early = x[:-30]   # x[0] ... x[149]，共150点
late  = x[30:]    # x[30] ... x[179]，共150点
products = late * np.conj(early)
corr = np.mean(products)
phase = np.angle(corr)
```

因此配对是 `(0,30), (1,31), ..., (149,179)`，一共 150 对。名义重叠窗长为 `3 ms - 0.5 ms = 2.5 ms`。

一般地，若固定一个包含 N 个点的有限窗口，并要求配对的两个点都在窗内，可用乘积数是 `N-L`。有的实现按样本对中心分配 gate，规则会不同；这份代码明确采用上述有限窗口规则。

## 7. 为什么复数乘积能消掉未知初相位？

在单点、匀速、无噪声、避开包络边缘的条件下：

$$
x[n]\approx A e^{j(2\pi f_dn/f_{s,IQ}+\phi_c)}.
$$

两个相隔 L 点的样本做共轭乘法：

$$
q[n]=x[n+L]x^*[n]
\approx A^2e^{j2\pi f_dL/f_{s,IQ}}.
$$

样本编号 n 与未知初相位都抵消。所有乘积近似指向同一角度，先对它们做复数平均，再取角度：

$$
\widehat R=\frac1{150}\sum_{n=0}^{149}q[n],
\qquad \widehat\phi=\arg\widehat R.
$$

实际运行得到：

$$
\widehat R\approx0.1549731745-j\,0.0397848135,
\qquad \widehat\phi\approx-0.2512939064\ \mathrm{rad}.
$$

这相当于 -14.398080 度。`lag_products.csv` 能逐对查看；`iq_gate.csv` 可以核算任意一对的乘积。150 对并不等于 150 次统计独立的测量。

## 8. 相位、频移和速度如何换算？

$$
\widehat f_d=\frac{\widehat\phi}{2\pi T_L}
\approx-79.98933475\ \mathrm{Hz}.
$$

由于规定远离为正速度，低速公式带负号：

$$
\widehat v\approx-\frac{c\widehat f_d}{2f_0}
=0.1999733369\ \mathrm{m/s}.
$$

真实输入为 +0.2 m/s。这里的微小差别主要是小速度近似，不是噪声。用精确传播式反解：

$$
\widehat v_{\mathrm{exact}}
=-\frac{c\widehat f_d}{2f_0+\widehat f_d}
\approx0.199999999996\ \mathrm{m/s}.
$$

精确式在这里用于核对模型；若相位已经混叠，直接换成精确式也无法自动恢复正确速度。

## 9. 已运行的检查

| 真实束向速度 m/s | 低速公式估计 m/s | 说明 |
|---:|---:|---|
| 0 | 约0 | 静止时频移近似零 |
| +0.2 | +0.199973 | 远离、频移为负 |
| -0.2 | -0.200027 | 靠近、频移为正 |
| +3.0 | -2.005988 | 相位绕回，故意展示失败 |

脚本检查了测得频移与理论绕回频移的一致性，并检查小速度下用精确式恢复真值。通过检查不表示已经验证了真实仪器。

本例的小速度近似无模糊范围为正负 `c/(4*f0*TL) = 2.5 m/s`。+3 m/s 案例超出范围后被错误映射到负速度；需要其他 lag 或可靠先验确定整周期数。

## 10. 接收噪声和散射群为何要分开？

脚本先给 RF 接收信号加入 20 dB 输入 SNR 的实噪声，重新走同样的解调滤波流程，结果见 JSON 的 `noisy_rf`。

另外以相同干净 I/Q 为基础，在基带加入 10 dB SNR 的复噪声，做300次独立试验。运行结果：均值约0.200311 m/s，标准差约0.010866 m/s。这里的10 dB定义在滤波后的 I/Q 窗口，不能与 RF 输入的20 dB直接混为一谈。

单点目标始终贡献同一组相位信息，因此理想条件下结果非常准确。真实水体包含许多散射体：有限脉冲扫过时，不同颗粒的复数回波不断加入和离开，总回波幅值、相位会起伏。

脚本增加了600个随机点散射体，全部以相同真实速度+0.2 m/s运动。为了减少蒙特卡洛计算量，直接计算同一个精确模型对应的复数基带；并未随机指定一个“速度噪声”。接收窗口固定为9至11 ms，共120点。200组独立散射配置的结果列在 `cloud_lag_comparison.csv`。

该固定窗口实验中，lag 从0.1 ms增至1.0 ms时，可用乘积从114对降至60对，速度标准差约从0.192增至0.320 m/s。它说明“其他条件相同、相位误差不变时，速度误差正比于1/TL”不能不加条件地应用。相位误差、有效独立样本数和模糊概率也会随lag变化。

这里更换的是独立散射配置，用来估计该模型的统计误差；它不是一个真实ping间隔和水体演化模型，也没有复现任何具体型号的单ping规格。

## 11. 图7的下降相关曲线与本例的相干系数为何不同？

本例输出的 `coherence` 定义为：

$$
\gamma=\frac{|\sum x[n+L]x^*[n]|}
{\sqrt{\sum|x[n+L]|^2\sum|x[n]|^2}}.
$$

求和只覆盖有效样本对。纯净单正弦在这个归一化下接近1，即使窗口缩短也如此；但有效数据减少会影响有噪声时的估计方差。

FST001图7的理想有限脉冲自相关包络，是相对于零lag的整段能量归一化，包含了重叠长度减少的效应，因而呈三角下降。这两种归一化不是同一个量。

散射群图中长lag处的估计相干系数可能再次上升：有效窗口变短、样本高度相关，有限样本的相干估计会受偏置影响。不能据此判定测速更准确。这里的0至1系数也不是厂商输出的相关计数值，没有为真实ADCP设定通用质量阈值。

## 12. 图册怎样读？

1. 第一页：A看发射/回波包络位置；B放大看真实载波；C显示脉冲被拉长约1.067 us；D显示混频后仍须滤除的高频项。B中的起始相位差还包含散射相位，不能仅凭这幅短局部图直观估计80 Hz小频移。
2. 第二页：A看I/Q；B看相位线性变化；C看150个相隔30点的乘积近似具有同一相位；D看相位差与lag的关系。图中极小的剩余数值误差来自有限滤波和浮点计算。
3. 第三页：把单点目标换成散射群后，I/Q不再是平稳单正弦，单次相位估计出现随机误差。后两图是在固定2 ms窗口下比较不同lag。

## 13. 资料对应关系

- FST001，第4.1节，印刷页3至4，Eq.1至3：相位、lag、频移和束向速度关系。
- FST001，印刷页6：同一深度单元内的相位测量并不都独立。
- BBPRIME，第3节：传播时间拉伸与多普勒频移。
- BBPRIME，第5节：接收时间门与实际水体单元的加权关系。

本例的精确匀速传播模型、I/Q滤波器参数、索引规则和随机场配置均为教学选择，文献并未公开指定这些具体实现参数。

原版链接：

- https://www.teledynemarine.com/en-us/support/SiteAssets/RDI/Technical%20Notes/FST%20Documents/FST001.PDF
- https://www.teledynemarine.com/en-us/support/SiteAssets/RDI/Manuals%20and%20Guides/General%20Interest/BBPRIME.pdf

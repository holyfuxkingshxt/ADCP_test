# 可验收的窄带算法验证数据

**当前默认版本（2026-10-30确认）：`narrowband_validation.mat`。** 后续编写和核对算法以此为主数据；相邻顶层和 `distributed_scene/` 中的旧MAT已标记为历史教学数据，不再作为默认单次精度验收输入。本目录的边界测试和连续散射集合测试继续保留使用。

2026-10-30从此MAT保存的完整 `rx_rf_noisy_tgc` 重新解调、滤波、抽取并计算六个gate的速度，最大绝对误差0.000405807519 m/s，满足既定的0.001 m/s容差；默认第三层真值0.2 m/s，精确反演结果0.199958941802 m/s。该次复核为只读计算，未重新生成或修改数据与既有验证报告。

**几何范围：** 当前为沿声束距离轴的一维模型，没有显式设置声束相对垂直方向的倾角。`truth.velocity_mps` 是束向速度；将它解释为水平速度或三维速度前，需要另加明确的投影几何或多束模型。

旧的连续散射单次回波具有很大的固有随机误差，不适合作为“单次计算应接近真值”的基础验收样例。新版受控薄层的通过结果也不代表均匀连续水体的单次测量都具备同样精度。

本目录提供两种明确的验证目标：受控散射层验证确定性的RF处理链路；连续散射群验证统计平均。所有接收信号先按声波往返传播生成，没有根据真值修改估计结果，也没有挑选能得到好结果的随机种子。

## 先运行

在MATLAB打开并运行本目录的 `run_validation.m`，会从三个测试文件的完整RF数据重新计算I/Q和速度，并对连续散射群进行集合验证。脚本使用MATLAB核心函数，不依赖额外工具箱；每个验证函数包含误差断言。`inspect_validation.m` 用于画主场景。

已在本机安装的MATLAB中实际执行 `run_validation.m`，进程正常退出，72项受控场景数值对比和8192组连续散射配置的集合平均断言全部通过。执行日志为 `matlab_validation.log`，`validation_report.json` 已记录这次成功执行。

仅处理默认的0.2 m/s区域，可以直接使用：

```matlab
folder = 'C:/Users/pc/Desktop/ADCP/Picture/narrowband_simulation/algorithm_validation';
S = load(fullfile(folder,'narrowband_validation.mat'));
p = S.params;

% 输入：完整含噪RF，经TGC，仍包含所有距离层。
rf = S.rx_rf_noisy_tgc;
mixed = 2*rf.*exp(-1j*2*pi*p.f0*S.t_rf_s);
baseband = conv(mixed,S.fir_coefficients,'same');
iq = baseband(1:p.decimation:end);
x = iq(S.gate_indices_matlab); % 180点；MATLAB索引871至1050。

L = round(p.lag_seconds*p.fs_rf/p.decimation); % 30点，0.5 ms。
R = mean(x(1+L:end).*conj(x(1:end-L)));        % 150对。
fd = angle(R)/(2*pi*p.lag_seconds);
vSmall = -p.c*fd/(2*p.f0);                  % ADCP常用的小速度近似。
vExact = -p.c*fd/(2*p.f0+fd);               % 与本仿真精确传播模型对应。
fprintf('truth %.8f, small-v %.8f, exact %.8f m/s\n', ...
    S.truth.selected_velocity_mps,vSmall,vExact);
assert(abs(vExact-S.truth.selected_velocity_mps)<p.noisy_tolerance_mps);
```

默认区域的结果为：物理真值0.2 m/s；无噪声、无TGC的精确反演约0.20000000 m/s；含噪声、TGC后的精确反演约0.19995894 m/s，小速度近似约0.19993229 m/s。

## 主数据的物理场景与适用范围

- 单声束、仪器固定；正速度表示散射体远离仪器，输出为束向速度，单位m/s。
- 声速1500 m/s，发射载频300 kHz，无编码4 ms正弦脉冲，RF采样率2.4 MHz，记录34 ms，共81600点。
- 六个散射薄层中心为4.5、7.5、10.5、13.5、16.5、19.5 m。每层厚0.16 m、240个随机复幅度散射体，共1440个。每个散射体按自己的初始距离和运动速度生成实数RF，再叠加。
- 六层速度依次为0、+0.1、+0.2、-0.2、+0.4、-0.4 m/s。同一薄层中的散射体速度相同。
- TGC保持上一版的衰减反补偿模型：以3 m为增益1，最高允许40 dB，21 m处实际约35.37 dB。噪声与回波一起增益，之后才进行解调。
- 噪声为固定种子的白噪声；按六个稳定接收窗中最弱窗的原始RF功率设为20 dB SNR，其余窗SNR更高。FIR会进一步限制噪声带宽。
- FIR为1601点Blackman窗低通，截止6 kHz；抽取40倍，I/Q采样率60 kHz。奇数长度FIR通过`same`卷积离线补偿800个RF点的群延迟。
- 每层各有一个3 ms gate，180个I/Q点，lag30点、150个共轭乘积。默认处理第三层。

这是物理上受控的散射薄层基准，不把它宣称为均匀连续水体的单次测速精度。通过几何约束，窗内及FIR整个支持区间中，本层所有散射体的回波保持存在，其他层的直接回波不重叠。因此复包络相位保持稳定，可以对算法做严格验收。完整接收记录仍包含所有六层及窗外回波边缘；生成过程没有先截取某一个gate。

程序逐层检查上述几何条件，并用闭式复幅度叠加与FIR频率响应独立计算窗内参考I/Q，检查其与完整RF处理结果一致。

## 验收标准与独立复核

事先固定的验收误差为：无噪声1e-5 m/s；含噪声1e-3 m/s。都针对精确反演，或针对小速度近似本身应得到的参考值，不混淆两种关系。

物理真值为v时，本仿真采用精确匀速传播：

$$
f_D=-\frac{2f_0v}{c+v},\qquad
\phi=\operatorname{wrap}(2\pi f_DT_L).
$$

无相位模糊时：

$$
\hat f_D=\frac{\arg\hat R}{2\pi T_L},\qquad
\hat v_{\mathrm{exact}}=-\frac{c\hat f_D}{2f_0+\hat f_D}.
$$

通常的小速度近似为 `vSmall=-c*fd/(2*f0)`。因此，代码另外保存 `truth.velocity_small_v_expected_mps`，方便使用常规ADCP公式的人逐项对照。例如精确真值0.2 m/s时，小速度近似理论上为0.19997334 m/s，这个小差异是模型近似，不是程序错误。

`verify_and_export.py`不调用生成程序的解调器或测速器，用SciPy重新设计同参数FIR、独立卷积，并从完整RF重新算出相关与速度。三个场景、四条路径、六个gate，共72项数值检查已通过；所有MAT数组均重新加载并逐项检查。每次重新导出时，报告中的MATLAB执行状态重置为false，实际MATLAB执行记录见 `matlab_validation.log`。

`validation_results.csv`保存每项真值、输出、误差和容差；`validation_report.json`为机器可读报告。图册为 `validation_figures.pdf`。

## 文件用途

| 文件 | 用途 |
|---|---|
| `narrowband_validation.mat` | 默认：完整多层RF，零速度、正负速度、噪声、TGC及已通过的真值对照 |
| `narrowband_boundary_valid.mat` | 第三层改为+2.4 m/s，验证接近模糊边界但仍可恢复的速度 |
| `narrowband_boundary_alias.mat` | 第三层改为+3 m/s，明确验证单lag折叠结果，状态为EXPECTED_ALIAS |
| `continuous_ensemble_validation.mat` | 连续散射群8192组独立配置，验证复相关集合平均 |
| `run_validation.m` | 在MATLAB依次运行全部断言 |
| `verify_narrowband_validation.m` | 从完整实数RF独立复算，可单独指定任一个受控场景MAT |
| `verify_continuous_ensemble.m` | 从保存的I/Q独立复算集合平均 |
| `inspect_validation.m` | 查看默认场景的完整RF、相位、速度及误差 |

单lag只能得到主值相位。+3 m/s案例的折叠频移约+802.395 Hz，对应精确公式的主值速度约-2.00331 m/s，小速度公式约-2.00599 m/s。EXPECTED_ALIAS只表示程序得到应有的折叠结果，不能解读为恢复了真实+3 m/s，或单lag算法能自行识别模糊区间。

## 默认MAT变量

| 变量 | 维度 | 含义 |
|---|---|---|
| `t_rf_s`, `tx_rf` | 81600×1 | 时间，秒；发射实数波形 |
| `rx_rf_clean`, `rx_rf_noisy` | 81600×1 | 全部六层叠加的补偿前RF |
| `rx_rf_zones` | 81600×6 | 每层各自的实数RF贡献 |
| `rx_rf_clean_tgc`, `rx_rf_noisy_tgc` | 81600×1 | 补偿后的完整RF |
| `noise_rf`, `noise_rf_tgc` | 81600×1 | 补偿前后的噪声 |
| `tgc_gain_rf`, `tgc_gain_db_rf` | 81600×1 | 增益倍数和幅度增益dB |
| `fir_coefficients` | 1601×1 | 低通FIR系数 |
| `t_iq_s`, `iq_total`, `iq_noisy`, `iq_total_tgc`, `iq_noisy_tgc` | 2040×1 | 抽取后时间与完整I/Q |
| `gate_indices_matlab`, `t_gate_s` | 180×1 | 默认第三层gate索引与时间 |
| `gate_indices_all_matlab`, `t_gates_s` | 180×6 | 六个gate索引与时间，每列一层 |
| `iq_gate`, `iq_gate_noisy`, `iq_gate_tgc`, `iq_gate_noisy_tgc` | 180×1 | 默认第三层gate的四条路径 |
| `iq_gates`, `iq_gates_noisy`, `iq_gates_tgc`, `iq_gates_noisy_tgc` | 180×6 | 六层gate的四条路径 |
| `oracle_iq_gates`, `oracle_iq_gates_tgc` | 180×6 | 闭式受控层参考I/Q，用于核对解调和滤波 |
| `lag_products`, `lag_products_noisy_tgc` | 150×1 | 默认第三层相关乘积 |
| `scatterer_*` | 1440×1 | 各散射体的初始距离、幅度、速度及回波起止时刻 |
| `params`, `truth`, `reference`, `verification` | struct | 参数、解析真值、生成器计算结果、独立复核结果 |

`truth.velocity_mps`、`truth.fd_true_hz`、`truth.phase_wrapped_expected_rad`均为六层的列向量；`truth.selected_velocity_mps`为默认第三层真值。所有MATLAB索引从1开始。先用`rx_rf_noisy_tgc`验证完整链路，也可用`iq_gate_noisy_tgc`单独验证相关模块。

## 连续散射群的统计验收

该文件沿用原来3–21 m连续散射场、1440个散射体、4 ms窄带脉冲、相同FIR和TGC。真实速度为0.2 m/s。输入 `iq_gates_tgc` 为180×8192，每列是一组独立散射配置；逐列做相关，再平均复数相关值，最后求相位和速度。固定验收误差为0.01 m/s；本次结果0.19798851 m/s，绝对误差0.00201149 m/s，通过。

此处不要求每一列的速度接近真值。8192组配置由固定种子独立抽取复散射幅度，不模拟真实相邻ping的时间演化。使用与完整RF处理一致的精确线性映射生成gate，映射已与保存的完整RF例子逐项核对；没有保存8192条完整RF。因此，该文件验证I/Q相关与统计平均，完整RF链路由前面的受控层文件验证。

## 重现

```powershell
& 'C:\Users\pc\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' 'C:\Users\pc\Desktop\ADCP\Picture\narrowband_simulation\algorithm_validation\generate_validation.py'
& 'C:\Users\pc\AppData\Local\Programs\Python\Python314\python.exe' 'C:\Users\pc\Desktop\ADCP\Picture\narrowband_simulation\algorithm_validation\verify_and_export.py'
& 'C:\Users\pc\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' 'C:\Users\pc\Desktop\ADCP\Picture\narrowband_simulation\algorithm_validation\plot_validation.py'
```

连续散射统计数据读取相邻目录 `distributed_scene` 的原场景数据。受控层数据生成使用相同的传播模型和基础接收链路；Python独立复核及MATLAB脚本分别使用自己的滤波、索引和估计实现。

"""Narrowband ADCP teaching example: one ideal constant-velocity echo.

Model: fixed instrument, positive velocity away from the instrument,
first-order Doppler, unity amplitude, ideal I/Q demodulation, no noise.
The I/Q sequence is generated analytically; ADC filtering is not modeled.
Run with Python, NumPy and Pillow to reproduce the timing diagram.
"""

from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


C = 1500.0                 # Sound speed, m/s
F0 = 300_000.0             # Transmit carrier, Hz
TP = 1e-3                 # Transmit duration, s
TL = 0.2e-3               # Correlation lag within the echo, s
F_ADC = 2_400_000.0        # Passband waveform sampling rate, samples/s
F_IQ = 50_000.0           # Ideal complex-baseband sampling rate, samples/s
V_TRUE = 1.0              # Beam velocity, m/s, positive away
V_MAX_TARGET = 2.0        # Intended test range: +/- this value, m/s

N_ADC = round(F_ADC * TP)
N_IQ = round(F_IQ * TP)
LAG = round(F_IQ * TL)
TL_ACTUAL = LAG / F_IQ

# Rectangularly windowed, single-frequency transmit pulse.
t_tx = np.arange(round(F_ADC * 1.2 * TP)) / F_ADC
tx = np.where(t_tx < TP, np.cos(2 * np.pi * F0 * t_tx), 0.0)

# Ideal I/Q echo within a single analysis record. The carrier is removed.
t_iq = np.arange(N_IQ) / F_IQ
f_d = -2 * F0 * V_TRUE / C
x = np.exp(1j * 2 * np.pi * f_d * t_iq)
r_lag = np.mean(x[LAG:] * np.conj(x[:-LAG]))
phase = np.angle(r_lag)
v_est = -C * phase / (4 * np.pi * F0 * TL_ACTUAL)
v_ambiguity = C / (4 * F0 * TL_ACTUAL)
range_scale = C * TP / 2

print(f"Transmit carrier cycles: {F0 * TP:.0f}")
print(f"ADC samples in transmit pulse: {N_ADC}")
print(f"I/Q record samples: {N_IQ}; lag samples: {LAG}")
print(f"Pulse-limited along-beam range scale: {range_scale:.3f} m")
print(f"Doppler shift: {f_d:.3f} Hz")
print(f"Correlation phase: {phase:.6f} rad")
print(f"Beam velocity: true={V_TRUE:.6f}, estimated={v_est:.6f} m/s")
print(f"Unambiguous velocity interval: +/-{v_ambiguity:.3f} m/s")
print(f"Lag bound for +/-{V_MAX_TARGET:g} m/s: "
      f"{1e3 * C / (4 * F0 * V_MAX_TARGET):.3f} ms")

def font(size, bold=False):
    p = Path(r"C:\Windows\Fonts") / ("msyhbd.ttc" if bold else "msyh.ttc")
    return ImageFont.truetype(str(p), size) if p.exists() else ImageFont.load_default()


img = Image.new("RGB", (1600, 980), "white")
draw = ImageDraw.Draw(img)
ink, blue, green, orange = "#203447", "#1c6aa6", "#27836a", "#bc621a"
left, right = 140, 1450


def label(x, y, text, size=28, fill=ink, bold=False, center=False):
    draw.text((x, y), text, font=font(size, bold), fill=fill, anchor="mt" if center else "lt")


def arrow(x0, x1, y, color=orange):
    draw.line((x0, y, x1, y), fill=color, width=3)
    draw.polygon([(x0, y), (x0 + 13, y - 7), (x0 + 13, y + 7)], fill=color)
    draw.polygon([(x1, y), (x1 - 13, y - 7), (x1 - 13, y + 7)], fill=color)


label(80, 45, "单束窄带仿真：T_P 与 T_L 的位置", 42, bold=True)
label(80, 110, "载频 300 kHz｜脉宽 1 ms｜相关间隔 0.2 ms", 28)

label(80, 200, "发射：一次持续 T_P 的单频脉冲", 30, bold=True)
draw.rounded_rectangle((left, 310, right, 420), radius=12, fill="#eaf3fa", outline=blue, width=2)
points = [(left + (right-left)*q, 365 - 35*np.sin(2*np.pi*12*q)) for q in np.linspace(0, 1, 1200)]
draw.line(points, fill=blue, width=3)
arrow(left, right, 280)
label((left+right)/2, 235, "T_P = 1 ms", 28, orange, center=True)
label(left, 432, "0", 24)
label(right, 432, "1 ms", 24, center=True)
label(80, 485, "实际包含 300 个载波周期；图中正弦周期仅作示意。", 25)

label(80, 575, "接收：在同一段回波中取相隔 T_L 的样本", 30, bold=True)
draw.rounded_rectangle((left, 695, right, 805), radius=12, fill="#eaf5ef", outline=green, width=2)
label((left+right)/2, 735, "同一次发射对应的回波记录", 28, green, center=True)
x0, x1 = left + 0.2*(right-left), left + 0.4*(right-left)
for x_mark in [x0, x1]:
    draw.line((x_mark, 680, x_mark, 815), fill=orange, width=3)
arrow(x0, x1, 655)
label((x0+x1)/2, 613, "T_L = 0.2 ms", 25, orange, center=True)
label(x0, 825, "x[n]", 25, orange, center=True)
label(x1, 825, "x[n+10]", 25, orange, center=True)
label(80, 900, "I/Q 采样率 50 kHz：每隔 20 微秒一个样本，0.2 ms 对应 10 个采样间隔。", 25)

fig_path = Path(__file__).resolve().parent / "narrowband_tp_tl.png"
img.save(fig_path)
print(f"Figure: {fig_path}")

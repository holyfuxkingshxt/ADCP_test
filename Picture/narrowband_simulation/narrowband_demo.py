"""Single-beam narrowband Doppler teaching simulation.

Requires numpy and reportlab. All outputs stay beside this script.
Positive radial velocity means moving AWAY from the fixed transducer.
The real RF echo is generated from the exact constant-speed propagation map.
The correlator operates on filtered complex I/Q, not on the transmit waveform.
"""
from dataclasses import dataclass, asdict, replace
from pathlib import Path
import csv
import json
import math

import numpy as np
from reportlab.pdfgen import canvas
from reportlab.lib import colors
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics.shapes import Drawing


@dataclass(frozen=True)
class Parameters:
    c: float = 1500.0
    f0: float = 300_000.0
    fs_rf: float = 2_400_000.0
    pulse_duration: float = 0.004
    initial_range: float = 3.0
    velocity_away: float = 0.2
    amplitude: float = 0.4
    scattering_phase: float = 0.7
    record_duration: float = 0.010
    fir_taps: int = 1601
    lowpass_hz: float = 6000.0
    decimation: int = 40
    gate_guard: float = 0.0005
    gate_samples: int = 180
    lag_seconds: float = 0.0005
    seed: int = 20261029


P = Parameters()
OUT = Path(__file__).resolve().parent


def propagation(p):
    # t_emit = alpha * t_receive - beta; R(t) = R0 + v*t.
    alpha = (p.c - p.velocity_away) / (p.c + p.velocity_away)
    beta = 2.0 * p.initial_range / (p.c + p.velocity_away)
    onset = beta / alpha
    duration = p.pulse_duration / alpha
    fd = p.f0 * (alpha - 1.0)
    return alpha, beta, onset, duration, fd


def generate_rf(p):
    t = np.arange(round(p.record_duration * p.fs_rf)) / p.fs_rf
    tx_gate = (t >= 0) & (t < p.pulse_duration)
    tx = tx_gate.astype(float) * np.cos(2 * np.pi * p.f0 * t)
    alpha, beta, onset, duration, fd = propagation(p)
    emit_time = alpha * t - beta
    rx_gate = (emit_time >= 0) & (emit_time < p.pulse_duration)
    rx = p.amplitude * rx_gate * np.cos(
        2 * np.pi * p.f0 * emit_time + p.scattering_phase
    )
    return t, tx, rx


def lowpass_kernel(p):
    k = np.arange(p.fir_taps) - (p.fir_taps - 1) / 2
    h = (2 * p.lowpass_hz / p.fs_rf
         * np.sinc(2 * p.lowpass_hz * k / p.fs_rf)
         * np.blackman(p.fir_taps))
    return h / h.sum()


def fir_aligned(x, h):
    # Linear convolution using zero padding, then compensate known FIR delay.
    # Offline time alignment uses future samples; a real-time implementation
    # would retain the causal delay and compensate its timestamps.
    size = len(x) + len(h) - 1
    nfft = 1 << (size - 1).bit_length()
    full = np.fft.ifft(np.fft.fft(x, nfft) * np.fft.fft(h, nfft))[:size]
    delay = (len(h) - 1) // 2
    return full[delay:delay + len(x)]


def demodulate(t, rx, p):
    mixed = 2.0 * rx * np.exp(-2j * np.pi * p.f0 * t)
    filtered = fir_aligned(mixed, lowpass_kernel(p))
    return t[::p.decimation], filtered[::p.decimation]


def select_gate(t_bb, iq, p):
    _, _, onset, _, _ = propagation(p)
    first = int(np.searchsorted(t_bb, onset + p.gate_guard))
    last = first + p.gate_samples  # Python end-exclusive indexing.
    return first, last, t_bb[first:last], iq[first:last]


def estimate_velocity(x, fs_bb, lag_seconds, p):
    lag = int(round(lag_seconds * fs_bb))
    assert 0 < lag < len(x)
    actual_lag = lag / fs_bb
    early, late = x[:-lag], x[lag:]
    products = late * np.conj(early)
    corr = np.mean(products)
    phase = float(np.angle(corr))
    fd = phase / (2 * np.pi * actual_lag)
    speed_small_v = -p.c * fd / (2 * p.f0)
    # Exact inversion is valid only after Doppler ambiguity has been resolved.
    speed_exact = -p.c * fd / (2 * p.f0 + fd)
    denom = np.sqrt(np.sum(abs(early) ** 2) * np.sum(abs(late) ** 2))
    coherence = float(np.clip(abs(np.sum(products)) / denom, 0, 1)) if denom else 0.0
    result = dict(samples=len(x), lag_samples=lag, pairs=len(products),
                  lag_seconds=actual_lag, correlation_real=float(corr.real),
                  correlation_imag=float(corr.imag), phase_rad=phase,
                  phase_degrees=float(np.degrees(phase)), doppler_hz=float(fd),
                  velocity_small_v=float(speed_small_v),
                  velocity_exact=float(speed_exact), coherence=coherence)
    return result, products


def single_target(p):
    t, tx, rx = generate_rf(p)
    tb, iq = demodulate(t, rx, p)
    first, last, tg, x = select_gate(tb, iq, p)
    result, q = estimate_velocity(x, p.fs_rf / p.decimation, p.lag_seconds, p)
    result.update(gate_first_index=first, gate_last_exclusive=last,
                  gate_first_time=float(tg[0]), gate_window_seconds=len(x)/(p.fs_rf/p.decimation))
    return result, (t, tx, rx, tb, iq, tg, x, q)


def cloud_iq(p, seed, fs_bb=60_000.0, count=600):
    """Distributed point-scatterer echo, evaluated directly in complex baseband.

    Each point has an exact Doppler-scaled rectangular return. Cumulative
    interval sums avoid synthesizing the RF carrier in Monte Carlo trials.
    No decorrelation or random Doppler is inserted artificially.
    """
    rng = np.random.default_rng(seed)
    t = np.arange(round(0.021 * fs_bb)) / fs_bb
    ranges = rng.uniform(0.5, 12.0, count)
    scatter = (rng.normal(size=count) + 1j * rng.normal(size=count)) / np.sqrt(2 * count)
    v = p.velocity_away
    alpha = (p.c - v) / (p.c + v)
    beta = 2 * ranges / (p.c + v)
    start = beta / alpha
    end = start + p.pulse_duration / alpha
    i0 = np.clip(np.searchsorted(t, start), 0, len(t))
    i1 = np.clip(np.searchsorted(t, end), 0, len(t))
    amp = scatter * np.exp(-2j * np.pi * p.f0 * beta)
    change = np.zeros(len(t) + 1, dtype=complex)
    np.add.at(change, i0, amp)
    np.add.at(change, i1, -amp)
    envelope = np.cumsum(change[:-1])
    fd = p.f0 * (alpha - 1)
    iq = envelope * np.exp(2j * np.pi * fd * t)
    k = np.arange(41) - 20
    h = (2 * p.lowpass_hz / fs_bb * np.sinc(2 * p.lowpass_hz / fs_bb * k)
         * np.blackman(41))
    iq = fir_aligned(iq, h / h.sum())
    gate = (t >= 0.009) & (t < 0.011)
    return t[gate], iq[gate]


def monte_carlo(p, baseline_gate):
    fs_bb = p.fs_rf / p.decimation
    rng = np.random.default_rng(p.seed + 1000)
    snr_db = 10.0  # Added AFTER I/Q processing; not raw-ADC SNR.
    power = float(np.mean(abs(baseline_gate) ** 2))
    sigma = np.sqrt(power / 10 ** (snr_db / 10) / 2)
    speeds, coherences = [], []
    for _ in range(300):
        noise = sigma * (rng.normal(size=len(baseline_gate)) + 1j * rng.normal(size=len(baseline_gate)))
        est, _ = estimate_velocity(baseline_gate + noise, fs_bb, p.lag_seconds, p)
        speeds.append(est['velocity_small_v'])
        coherences.append(est['coherence'])
    noise_result = dict(trials=300, iq_snr_db=snr_db,
                        mean_velocity=float(np.mean(speeds)),
                        std_velocity=float(np.std(speeds, ddof=1)),
                        mean_coherence=float(np.mean(coherences)))
    lags = np.array([0.0001, 0.00025, 0.0005, 0.001, 0.0015])
    cloud_speeds = np.zeros((200, len(lags)))
    cloud_coherence = np.zeros_like(cloud_speeds)
    for trial in range(200):
        _, x = cloud_iq(p, p.seed + 2000 + trial, fs_bb)
        for j, lag in enumerate(lags):
            est, _ = estimate_velocity(x, fs_bb, lag, p)
            cloud_speeds[trial, j] = est['velocity_small_v']
            cloud_coherence[trial, j] = est['coherence']
    rows = [dict(lag_ms=float(lag * 1000), pairs=120-int(round(lag*fs_bb)),
                 mean_velocity=float(cloud_speeds[:, j].mean()),
                 std_velocity=float(cloud_speeds[:, j].std(ddof=1)),
                 mean_coherence=float(cloud_coherence[:, j].mean()))
            for j, lag in enumerate(lags)]
    return noise_result, rows, speeds, cloud_speeds


PALETTE = [colors.HexColor(s) for s in ['#125A92','#D46A30','#31856C','#8A4F9E']]


def panel(pdf, box, series, title, xlabel, ylabel, bounds=None):
    """ReportLab standard LinePlot; data stay numerical and export as vectors."""
    left, bottom, width, height = box
    pdf.setFillColor(colors.HexColor('#152E40'))
    pdf.setFont('Helvetica-Bold', 11)
    pdf.drawString(left + 51, bottom + height - 16, title)
    drawing = Drawing(width, height)
    chart = LinePlot()
    chart.x, chart.y = 54, 43
    chart.width, chart.height = width-72, height-77
    converted = []
    for sx, sy, label in series:
        sx, sy = np.asarray(sx), np.asarray(sy)
        keep = np.linspace(0, len(sx)-1, min(1800,len(sx))).astype(int)
        converted.append([(float(sx[k]), float(sy[k])) for k in keep])
    chart.data = converted
    for j in range(len(series)):
        chart.lines[j].strokeColor = PALETTE[j % len(PALETTE)]
        chart.lines[j].strokeWidth = 1.05
    for axis in [chart.xValueAxis, chart.yValueAxis]:
        axis.labels.fontName = 'Helvetica'
        axis.labels.fontSize = 8
        axis.strokeColor = colors.HexColor('#69808F')
        axis.labelTextFormat = '%g'
    if bounds:
        chart.xValueAxis.valueMin, chart.xValueAxis.valueMax = bounds[:2]
        if len(bounds) == 4:
            chart.yValueAxis.valueMin, chart.yValueAxis.valueMax = bounds[2:]
    drawing.add(chart)
    drawing.drawOn(pdf, left, bottom)
    pdf.setFont('Helvetica', 9)
    pdf.drawCentredString(left+width/2, bottom+12, xlabel)
    pdf.saveState()
    pdf.translate(left+10,bottom+height/2)
    pdf.rotate(90)
    pdf.drawCentredString(0,0,ylabel)
    pdf.restoreState()
    legend_x = left + 58
    for j, (_, _, label) in enumerate(series):
        pdf.setStrokeColor(PALETTE[j % len(PALETTE)])
        pdf.line(legend_x,bottom+height-30,legend_x+15,bottom+height-30)
        pdf.setFillColor(colors.HexColor('#253D4A'))
        pdf.setFont('Helvetica',8)
        pdf.drawString(legend_x+19,bottom+height-33,label)
        legend_x += 30 + len(label)*4.5


def page_title(pdf, title, subtitle, number):
    pdf.setFillColor(colors.HexColor('#102C41'))
    pdf.setFont('Helvetica-Bold',20)
    pdf.drawString(35,752,title)
    pdf.setFont('Helvetica',10)
    pdf.drawString(35,730,subtitle)
    pdf.setFillColor(colors.HexColor('#526876'))
    pdf.setFont('Helvetica',8)
    pdf.drawString(35,18,'Narrowband teaching simulation | Positive radial velocity = away | Offline FIR delay compensated')
    pdf.drawRightString(1005,18,str(number))


def write_plots(p, baseline, mc_rows):
    t, tx, rx, tb, iq, tg, x, q = baseline
    pdf = canvas.Canvas(str(OUT/'narrowband_figures.pdf'),pagesize=(1040,780))
    pdf.setTitle('Narrowband Doppler: waveform, I/Q and correlation')
    boxes=[(25,390,490,320),(525,390,490,320),(25,50,490,320),(525,50,490,320)]
    _, _, onset, duration, fd = propagation(p)
    page_title(pdf,'1 / Transmit and receive waveforms','Single point scatterer: c=1500 m/s, f0=300 kHz, Tp=4 ms, R0=3 m, v=+0.2 m/s',1)
    tx_envelope=(t<p.pulse_duration).astype(float)
    rx_envelope=p.amplitude*((t>=onset)&(t<onset+duration))
    panel(pdf,boxes[0],[(t*1000,tx_envelope,'TX envelope'),(t*1000,rx_envelope,'RX envelope')],
          'A. Pulse envelopes locate TX and RX','Time from transmit start (ms)','Envelope amplitude',(0,10,0,1.1))
    mask=(t>=0)&(t<0.00004)
    rx_local=(t>=onset)&(t<onset+0.00004)
    panel(pdf,boxes[1],[(t[mask]*1e6,tx[mask],'TX'),((t[rx_local]-onset)*1e6,rx[rx_local]/p.amplitude,'RX / A')],
          'B. First 40 us: many carrier cycles','Local time after each onset (us)','Normalized amplitude',(0,40,-1.1,1.1))
    emit=np.linspace(0,p.pulse_duration,200)
    alpha,beta,*_=propagation(p)
    delay_extra=emit/alpha-emit
    panel(pdf,boxes[2],[(emit*1000,delay_extra*1e6,'Extra receive spacing')],
          'C. Time stretching (onset removed)','Transmit-time separation (ms)','Extra receive separation (us)')
    mask=(t>=onset)&(t<onset+0.00004)
    mixed=2*rx[mask]*np.exp(-2j*np.pi*p.f0*t[mask])
    panel(pdf,boxes[3],[((t[mask]-onset)*1e6,mixed.real,'Mixed real'),((t[mask]-onset)*1e6,mixed.imag,'Mixed imag')],
          'D. Mixer output before low-pass','Time after echo onset (us)','Amplitude')
    pdf.showPage()
    page_title(pdf,'2 / I/Q, phase and lag correlation','I/Q sample rate = 60 kHz | Gate = 180 samples (3 ms) | Lag = 30 samples (0.5 ms) | 150 products',2)
    local=(tg-tg[0])*1000
    panel(pdf,boxes[0],[(local,x.real,'I'),(local,x.imag,'Q')],
          'A. Filtered I/Q inside selected echo gate','Time within gate (ms)','Amplitude')
    phase=np.unwrap(np.angle(x)); phase-=phase[0]
    theory=2*np.pi*fd*(tg-tg[0])
    panel(pdf,boxes[1],[(local,phase,'Measured'),(local,theory,'Exact model')],
          'B. Phase accumulates at the Doppler rate','Time within gate (ms)','Relative unwrapped phase (rad)')
    panel(pdf,boxes[2],[(np.arange(len(q)),np.degrees(np.angle(q)),'Each lag product')],
          'C. Every product compares samples 30 apart','Product index n','Phase difference (degrees)',
          (0,len(q)-1,-14.5,-14.3))
    lags=np.arange(1,91)
    ph=[np.angle(np.mean(x[k:]*np.conj(x[:-k]))) for k in lags]
    panel(pdf,boxes[3],[(lags/60,ph,'Measured'),(lags/60,2*np.pi*fd*lags/60000,'Exact model')],
          'D. Different lag, same constant Doppler','Lag (ms)','Correlation phase (rad)')
    pdf.showPage()
    tc,xc=cloud_iq(p,p.seed+2000)
    page_title(pdf,'3 / A cloud of scatterers changes the error','600 moving scatterers | Same constant true v=+0.2 m/s | 2 ms receive gate | No added receiver noise',3)
    panel(pdf,boxes[0],[((tc-tc[0])*1000,xc.real,'I'),((tc-tc[0])*1000,xc.imag,'Q')],
          'A. Echo is a changing sum of many phasors','Time within range gate (ms)','Amplitude')
    phase=np.unwrap(np.angle(xc)); phase-=phase[0]
    panel(pdf,boxes[1],[((tc-tc[0])*1000,phase,'Cloud echo'),((tc-tc[0])*1000,2*np.pi*fd*(tc-tc[0]),'Single-target slope')],
          'B. Scatterer turnover adds phase fluctuations','Time within range gate (ms)','Relative phase (rad)')
    lm=[r['lag_ms'] for r in mc_rows]
    panel(pdf,boxes[2],[(lm,[r['std_velocity'] for r in mc_rows],'200 cloud realizations')],
          'C. Measured single-ping random spread','Lag (ms)','Velocity standard deviation (m/s)')
    panel(pdf,boxes[3],[(lm,[r['mean_coherence'] for r in mc_rows],'Cloud mean coherence')],
          'D. Fixed window: coherence depends on lag','Lag (ms)','Pair-normalized coherence',(0,1.6,0,1.05))
    pdf.save()


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    baseline_result, baseline = single_target(P)
    validation=[]
    for v in [0.0,0.2,-0.2,3.0]:
        p=replace(P,velocity_away=v)
        est,_=single_target(p)
        *_,fd=propagation(p)
        lag=est['lag_seconds']
        expected_wrapped=float(np.angle(np.exp(2j*np.pi*fd*lag))/(2*np.pi*lag))
        assert abs(est['doppler_hz']-expected_wrapped)<0.02
        if abs(v)<1:
            assert abs(est['velocity_exact']-v)<5e-5
        validation.append(dict(true_velocity=v,true_doppler_hz=fd,
                               measured_doppler_hz=est['doppler_hz'],
                               measured_velocity=est['velocity_small_v'],
                               ambiguity_expected=abs(fd*lag)>0.5))
    noise,cloud_rows,speeds,cloud_speeds=monte_carlo(P,baseline[6])
    rng=np.random.default_rng(P.seed)
    rf_noise_std=np.sqrt(P.amplitude**2/2/10**(20/10))
    noisy_rf=baseline[2]+rng.normal(0,rf_noise_std,len(baseline[2]))
    nt,niq=demodulate(baseline[0],noisy_rf,P)
    _,_,_,ng=select_gate(nt,niq,P)
    noisy_rf_result,_=estimate_velocity(ng,P.fs_rf/P.decimation,P.lag_seconds,P)
    noisy_rf_result['input_rf_snr_db']=20.0
    first=baseline_result['gate_first_index']; last=baseline_result['gate_last_exclusive']
    *_,tg,x,q=baseline
    with (OUT/'iq_gate.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.writer(f); writer.writerow(['local_index','full_iq_index','time_s','I','Q','magnitude','phase_rad'])
        for n in range(len(x)):
            writer.writerow([n,first+n,tg[n],x[n].real,x[n].imag,abs(x[n]),np.angle(x[n])])
    with (OUT/'lag_products.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.writer(f); writer.writerow(['pair_index','early_local_index','late_local_index','product_real','product_imag','product_phase_rad'])
        L=baseline_result['lag_samples']
        for n,z in enumerate(q):
            writer.writerow([n,n,n+L,z.real,z.imag,np.angle(z)])
    alpha,beta,onset,duration,fd=propagation(P)
    results=dict(parameters=asdict(P),truth=dict(time_scale=1/alpha,receive_onset_s=onset,
                   receive_duration_s=duration,exact_doppler_hz=fd,
                   small_v_doppler_hz=-2*P.f0*P.velocity_away/P.c),
                 baseline=baseline_result,validation=validation,
                 noisy_rf=noisy_rf_result,iq_noise_monte_carlo=noise,cloud_monte_carlo=cloud_rows,
                 limitations=['Single beam only; no three-dimensional reconstruction.',
                              'Point reflectors, constant c and radial velocity; amplitudes fixed.',
                              'No propagation loss, transducer beam pattern, multipath or hardware ringing.',
                              'Cloud trials use direct baseband evaluation of the same moving-reflector model.',
                              'Monte Carlo realizations use independent scatterer configurations, not a real ping-rate model.',
                              'Pair coherence is normalized over overlap; it is not a vendor correlation-count field.'])
    (OUT/'results.json').write_text(json.dumps(results,indent=2,ensure_ascii=False),encoding='utf-8')
    np.savez_compressed(OUT/'signals.npz',t_rf_s=baseline[0],tx_rf=baseline[1],
                        rx_rf=baseline[2],t_iq_s=baseline[3],iq=baseline[4],
                        gate_time_s=baseline[5],gate_iq=baseline[6],lag_products=baseline[7],
                        receiver_noise_rf=noisy_rf)
    with (OUT/'cloud_lag_comparison.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=cloud_rows[0].keys()); w.writeheader(); w.writerows(cloud_rows)
    write_plots(P,baseline,cloud_rows)
    print(json.dumps(dict(baseline=baseline_result,validation=validation,noise=noise,cloud=cloud_rows),indent=2))


if __name__=='__main__':
    main()

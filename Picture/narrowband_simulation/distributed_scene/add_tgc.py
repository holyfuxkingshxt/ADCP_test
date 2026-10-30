"""Apply receive-time gain to the complete RF record, including receiver noise.

Run on an existing distributed_signals.npz, then export_distributed_mat.py.
The full scene generator also calls these functions before saving its outputs.
"""
from pathlib import Path
import json
import sys

import numpy as np
from reportlab.pdfgen import canvas

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from narrowband_demo import Parameters, demodulate, estimate_velocity, panel, page_title


def gain_curve(t, cfg):
    # A long echo's representative range uses the TX pulse centre, not its front.
    # Clamping prevents singular near-field extrapolation and unlimited late gain.
    nominal_range = cfg['c'] * (t - cfg['pulse_duration'] / 2) / 2
    r = np.clip(nominal_range, cfg['tgc_min_range_m'], cfg['tgc_max_range_m'])
    ref = cfg['tgc_reference_range_m']
    # Inverse of the point-scatterer PRESSURE loss, normalized at ref.
    gain_db = (40 * np.log10(r / ref)
               + (40 / np.log(10)) * cfg['absorption_np_per_m'] * (r - ref))
    gain_db = np.clip(gain_db, 0, cfg['tgc_max_gain_db'])
    return r, 10 ** (gain_db / 20), gain_db


def add_tgc_data(arrays, stats, cfg):
    cfg.setdefault('tgc_reference_range_m', float(cfg['zone_edges_m'][0]))
    cfg.setdefault('tgc_min_range_m', float(cfg['zone_edges_m'][0]))
    cfg.setdefault('tgc_max_range_m', float(cfg['zone_edges_m'][-1]))
    cfg.setdefault('tgc_max_gain_db', 40.0)
    r, gain, gain_db = gain_curve(arrays['t_rf_s'], cfg)
    p = Parameters(c=cfg['c'], f0=cfg['f0'], fs_rf=cfg['fs_rf'],
                   fir_taps=cfg['fir_taps'], lowpass_hz=cfg['lowpass_hz'],
                   decimation=cfg['decimation'])
    clean = gain * arrays['rx_rf_clean']
    noise = gain * arrays['noise_rf']
    noisy = gain * arrays['rx_rf_noisy']
    # Gain is applied BEFORE the same mixer, FIR and decimator as the raw path.
    _, iq = demodulate(arrays['t_rf_s'], clean, p)
    _, iq_noisy = demodulate(arrays['t_rf_s'], noisy, p)
    _, iq_noise = demodulate(arrays['t_rf_s'], noise, p)
    iq_zones = np.column_stack([
        demodulate(arrays['t_rf_s'], gain * arrays['rx_rf_zones'][:, z], p)[1]
        for z in range(arrays['rx_rf_zones'].shape[1])
    ])
    g = arrays['gate_indices_python']
    fs_iq = cfg['fs_rf'] / cfg['decimation']
    clean_est, q = estimate_velocity(iq[g], fs_iq, cfg['lag_seconds'], p)
    noisy_est, q_noisy = estimate_velocity(iq_noisy[g], fs_iq, cfg['lag_seconds'], p)
    unity = cfg.copy()
    unity['tgc_max_gain_db'] = 0.0
    assert np.all(gain_curve(arrays['t_rf_s'], unity)[1] == 1)
    assert np.isfinite(noisy).all() and np.all(gain >= 1)
    assert np.all(np.diff(gain) >= 0)
    assert np.max(gain_db) <= cfg['tgc_max_gain_db'] + 1e-12
    recovery = float(np.max(abs(noisy / gain - arrays['rx_rf_noisy'])))
    np.testing.assert_allclose(noisy / gain, arrays['rx_rf_noisy'], rtol=2e-15, atol=1e-15)
    np.testing.assert_allclose(noisy, clean + noise, rtol=2e-15, atol=1e-15)
    np.testing.assert_allclose(iq_noisy, iq + iq_noise, rtol=1e-12, atol=1e-13)
    np.testing.assert_allclose(iq, iq_zones.sum(axis=1), rtol=1e-12, atol=1e-13)
    # The positive real gain multiplies signal and noise equally at EACH RF sample.
    # A ratio of window-summed powers can change because TGC also changes weights.
    active = (abs(arrays['rx_rf_clean']) > 1e-10) & (abs(arrays['noise_rf']) > 1e-10)
    gain_signal = clean[active] / arrays['rx_rf_clean'][active]
    gain_noise = noise[active] / arrays['noise_rf'][active]
    ratio_error = float(np.max(abs(gain_signal / gain_noise - 1)))
    assert ratio_error < 1e-12
    arrays.update(
        tgc_range_rf_m=r, tgc_gain_rf=gain, tgc_gain_db_rf=gain_db,
        rx_rf_clean_tgc=clean, rx_rf_noisy_tgc=noisy, noise_rf_tgc=noise,
        iq_total_tgc=iq, iq_noisy_tgc=iq_noisy, iq_noise_tgc=iq_noise,
        iq_zones_tgc=iq_zones, iq_gate_tgc=iq[g], iq_gate_noisy_tgc=iq_noisy[g],
        iq_gate_zones_tgc=iq_zones[g], lag_products_tgc=q,
        lag_products_noisy_tgc=q_noisy,
    )
    stats['configuration'] = cfg
    stats['time_gain_compensation'] = dict(
        model='Positive real RF gain; inverse point-scatterer pressure attenuation at representative range',
        representative_range='clip(c*(t-Ttx/2)/2, min_range, max_range); t starts at transmit start',
        normalization='Gain = 1 at 3 m for the default scene; max configured gain 40 dB',
        processing_order='Sum all RF echoes, add receiver noise, apply gain, mix, FIR, decimate, select gate',
        applied_max_gain=float(gain.max()), applied_max_gain_db=float(gain_db.max()),
        rf_recovery_max_error=recovery, signal_noise_gain_ratio_max_error=ratio_error,
        iq_zone_sum_max_error=float(np.max(abs(iq - iq_zones.sum(axis=1)))),
        estimates=dict(clean_before=stats['estimates']['uniform_total'], clean_after=clean_est,
                       noisy_before=stats['estimates']['noisy_total'], noisy_after=noisy_est),
        limitations=[
            'Digital gain after RF sampling: no analogue receiver, clipping or quantization model.',
            'Amplifies receiver noise too; cannot improve samplewise input SNR.',
            'A finite TX pulse mixes ranges at one receive time; one time gain cannot undo every scatterer loss.',
            'Does not remove speckle or correct the simplified constant point density per metre.',
            'Gain reweights lag products and filtering; a finite-gate speed estimate may change.',
        ],
    )
    return arrays, stats


def make_tgc_plots(arrays, stats, cfg):
    pdf = canvas.Canvas(str(HERE / 'tgc_comparison.pdf'), pagesize=(1040, 780))
    pdf.setTitle('Time gain compensation of the full distributed RF record')
    boxes = [(25,390,490,320), (525,390,490,320), (25,50,490,320), (525,50,490,320)]
    t = arrays['t_rf_s']
    block = 240  # 0.1 ms RMS, not carrier subsampling.
    usable = len(t) // block * block
    times = t[:usable].reshape(-1, block).mean(axis=1) * 1000
    def rms(name):
        return np.sqrt(np.mean(arrays[name][:usable].reshape(-1, block)**2, axis=1))
    peak_db = stats['time_gain_compensation']['applied_max_gain_db']
    page_title(pdf, 'Receive-time gain compensation',
               f'Full 3-21 m scene | Gain referenced to 3 m | Applied maximum {peak_db:.2f} dB | Before I/Q filtering', 1)
    panel(pdf, boxes[0], [(t*1000, arrays['tgc_gain_db_rf'], 'TGC curve')],
          'A. Positive real gain increases with receive time', 'Receive time (ms)', 'Amplitude gain (dB)', (0,34,0,40))
    panel(pdf, boxes[1], [(times,rms('rx_rf_clean'),'Raw'), (times,rms('rx_rf_clean_tgc'),'With TGC')],
          'B. Clean echo: 0.1 ms RMS', 'Receive time (ms)', 'RF RMS amplitude')
    panel(pdf, boxes[2], [(times,rms('rx_rf_noisy'),'Raw'), (times,rms('rx_rf_noisy_tgc'),'With TGC')],
          'C. Full noisy received record: 0.1 ms RMS', 'Receive time (ms)', 'RF RMS amplitude')
    panel(pdf, boxes[3], [(times,rms('noise_rf'),'Raw noise'), (times,rms('noise_rf_tgc'),'With TGC')],
          'D. Later receiver noise is amplified too', 'Receive time (ms)', 'Noise RMS amplitude')
    pdf.save()


def main():
    stats = json.loads((HERE/'distributed_results.json').read_text(encoding='utf-8'))
    with np.load(HERE/'distributed_signals.npz') as data:
        arrays = {key: data[key].copy() for key in data.files}
    cfg = stats['configuration']
    arrays, stats = add_tgc_data(arrays, stats, cfg)
    np.savez_compressed(HERE/'distributed_signals.npz', **arrays)
    (HERE/'distributed_results.json').write_text(json.dumps(stats, indent=2), encoding='utf-8')
    make_tgc_plots(arrays, stats, cfg)
    print(json.dumps(stats['time_gain_compensation'], indent=2))


if __name__ == '__main__':
    main()

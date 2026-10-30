"""Known-answer narrowband tests: controlled layers + continuous-scatterer ensemble.

The controlled layers test the RF algorithm against a stable, physical echo.
The ensemble separately tests the expected correlation in continuous scatter.
No velocity is corrected using its truth, and no random seeds are searched.
Run with the bundled NumPy/ReportLab Python, then verify_and_export.py with SciPy.
"""
from pathlib import Path
import sys
import json
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
sys.path.insert(0, str(BASE))
sys.path.insert(0, str(BASE/'distributed_scene'))
from narrowband_demo import Parameters, demodulate, lowpass_kernel, estimate_velocity
from distributed_demo import CONFIG, rf_scene
from add_tgc import gain_curve

CENTRES = np.array([4.5, 7.5, 10.5, 13.5, 16.5, 19.5])
VELOCITIES = np.array([0., .1, .2, -.2, .4, -.4])
CASES = [('narrowband_validation', .2),
         ('narrowband_boundary_valid', 2.4),
         ('narrowband_boundary_alias', 3.)]


def params(cfg):
    return Parameters(c=cfg['c'], f0=cfg['f0'], fs_rf=cfg['fs_rf'],
                      fir_taps=cfg['fir_taps'], lowpass_hz=cfg['lowpass_hz'],
                      decimation=cfg['decimation'])


def save(name, arrays, meta):
    np.savez_compressed(HERE/(name+'.npz'), **arrays)
    (HERE/(name+'.json')).write_text(json.dumps(meta, indent=2), encoding='utf-8')


def layer_oracle(gates, ranges, zones, amp, velocities, gain, cfg):
    """Closed-form, in-gate IQ, independent of the full RF scene summation.

    All particles in this layer remain present over the entire FIR support.
    Explicitly includes the small negative-frequency image from real RF.
    """
    h = lowpass_kernel(params(cfg)); k = np.arange(len(h))-(len(h)-1)//2
    raw = np.zeros(gates.shape, complex); tgc = np.zeros_like(raw)
    for z in range(6):
        own = zones == z; v = velocities[z]
        alpha = (cfg['c']-v)/(cfg['c']+v)
        beta = 2*ranges[own]/(cfg['c']+v)
        B = np.sum(amp[own]*np.exp(-2j*np.pi*cfg['f0']*beta))
        n = gates[:, z]*cfg['decimation']
        for freq, a in [(cfg['f0']*(alpha-1), B), (-cfg['f0']*(alpha+1), np.conj(B))]:
            kernel = h*np.exp(-2j*np.pi*freq*k/cfg['fs_rf'])
            tone = a*np.exp(2j*np.pi*freq*n/cfg['fs_rf'])
            raw[:, z] += tone*kernel.sum()
            tgc[:, z] += tone*np.sum(gain[n[:, None]-k[None, :]]*kernel[None, :], axis=1)
    return raw, tgc


def controlled_cases():
    cfg = CONFIG.copy()
    cfg.update(tgc_reference_range_m=3., tgc_min_range_m=3., tgc_max_range_m=21.,
               tgc_max_gain_db=40., layer_half_width_m=.08,
               clean_tolerance_mps=1e-5, noisy_tolerance_mps=1e-3,
               input_snr_weakest_gate_db=20.)
    rng = np.random.default_rng(20261029)
    count = cfg['scatterers_per_zone']
    ranges = (CENTRES[:, None]+rng.uniform(-.08, .08, (6, count))).ravel()
    zones = np.repeat(np.arange(6), count)
    amp = (rng.normal(size=len(ranges))+1j*rng.normal(size=len(ranges)))/np.sqrt(2*count)
    amp *= cfg['amplitude_scale']*(10/ranges)**2*np.exp(-2*cfg['absorption_np_per_m']*ranges)
    t = np.arange(round(cfg['record_duration']*cfg['fs_rf']))/cfg['fs_rf']
    _, gain, gain_db = gain_curve(t, cfg)
    fs_iq = cfg['fs_rf']/cfg['decimation']
    starts = 2*CENTRES/cfg['c']+.0005
    gates = np.column_stack([np.arange(round(s*fs_iq), round(s*fs_iq)+180) for s in starts])
    selected = gates[:, 2]; L = round(cfg['lag_seconds']*fs_iq)
    h = lowpass_kernel(params(cfg)); half = (len(h)-1)/(2*cfg['fs_rf'])
    tx = (t < cfg['pulse_duration'])*np.cos(2*np.pi*cfg['f0']*t)
    for case_index, (name, target_v) in enumerate(CASES):
        velocities = VELOCITIES.copy(); velocities[2] = target_v
        rxz, onset, end = rf_scene(t, ranges, zones, amp, velocities[zones], cfg)
        rx = rxz.sum(axis=1)
        # Check the actual filter-support interval, not just the gate midpoint.
        geometry = []
        for z in range(6):
            a = gates[0, z]/fs_iq-half; b = gates[-1, z]/fs_iq+half
            own = zones == z
            stable = bool(np.all(onset[own] <= a) and np.all(end[own] > b))
            others_absent = bool(np.all((end[~own] <= a) | (onset[~own] > b)))
            assert stable and others_absent, (name, z)
            geometry.append(dict(zone_matlab=z+1, all_own_returns_present=stable,
                                 other_layers_absent_over_fir=others_absent))
        rf_masks = [(t >= s) & (t < s+.003) for s in starts]
        powers = np.array([np.mean(rx[m]**2) for m in rf_masks])
        sigma = np.sqrt(powers.min()/10**(cfg['input_snr_weakest_gate_db']/10))
        noise = np.random.default_rng(20261030+case_index).normal(0, sigma, len(t))
        noisy = rx+noise
        tb, iq = demodulate(t, rx, params(cfg))
        _, iqn = demodulate(t, noisy, params(cfg))
        _, iqt = demodulate(t, gain*rx, params(cfg))
        _, iqnt = demodulate(t, gain*noisy, params(cfg))
        oracle, oracle_tgc = layer_oracle(gates, ranges, zones, amp, velocities, gain, cfg)
        err = float(np.max(abs(iq[gates]-oracle)))
        err_tgc = float(np.max(abs(iqt[gates]-oracle_tgc)))
        assert err < 1e-10 and err_tgc < 1e-10
        fd = -2*cfg['f0']*velocities/(cfg['c']+velocities)
        phase = np.angle(np.exp(2j*np.pi*fd*cfg['lag_seconds']))
        wrapped = phase/(2*np.pi*cfg['lag_seconds'])
        ambiguous = abs(fd) >= 1/(2*cfg['lag_seconds'])
        expected_small = -cfg['c']*wrapped/(2*cfg['f0'])
        expected_exact = -cfg['c']*wrapped/(2*cfg['f0']+wrapped)
        reference = {}
        for label, x in [('clean', iq), ('noisy', iqn), ('clean_tgc', iqt), ('noisy_tgc', iqnt)]:
            reference[label] = [estimate_velocity(x[gates[:, z]], fs_iq, cfg['lag_seconds'], params(cfg))[0]
                                for z in range(6)]
        arrays = dict(t_rf_s=t, tx_rf=tx, rx_rf_clean=rx, rx_rf_noisy=noisy,
                      noise_rf=noise, rx_rf_zones=rxz, tgc_gain_rf=gain, tgc_gain_db_rf=gain_db,
                      rx_rf_clean_tgc=gain*rx, rx_rf_noisy_tgc=gain*noisy, noise_rf_tgc=gain*noise,
                      fir_coefficients=h, t_iq_s=tb, iq_total=iq, iq_noisy=iqn,
                      iq_total_tgc=iqt, iq_noisy_tgc=iqnt,
                      gate_indices_all_matlab=gates+1, gate_indices_matlab=selected+1,
                      t_gates_s=tb[gates], t_gate_s=tb[selected],
                      iq_gates=iq[gates], iq_gates_noisy=iqn[gates],
                      iq_gates_tgc=iqt[gates], iq_gates_noisy_tgc=iqnt[gates],
                      iq_gate=iq[selected], iq_gate_noisy=iqn[selected],
                      iq_gate_tgc=iqt[selected], iq_gate_noisy_tgc=iqnt[selected],
                      lag_products=iq[selected[L:]]*np.conj(iq[selected[:-L]]),
                      lag_products_noisy_tgc=iqnt[selected[L:]]*np.conj(iqnt[selected[:-L]]),
                      oracle_iq_gates=oracle, oracle_iq_gates_tgc=oracle_tgc,
                      scatterer_range0_m=ranges, scatterer_zone_matlab=zones+1,
                      scatterer_amplitude_complex=amp, scatterer_velocity_mps=velocities[zones],
                      scatterer_echo_start_s=onset, scatterer_echo_end_s=end)
        truth = dict(velocity_mps=velocities.tolist(), fd_true_hz=fd.tolist(),
                     fd_wrapped_expected_hz=wrapped.tolist(), phase_wrapped_expected_rad=phase.tolist(),
                     velocity_small_v_expected_mps=expected_small.tolist(),
                     velocity_exact_from_wrapped_expected_mps=expected_exact.tolist(),
                     is_ambiguous=ambiguous.tolist(), selected_zone_matlab=3,
                     selected_velocity_mps=float(target_v), range_centres_m=CENTRES.tolist())
        meta = dict(params=cfg, truth=truth, reference=reference,
                    model='Controlled thin scattering layers, not a continuous water-column accuracy claim',
                    geometry_checks=geometry, oracle_iq_max_error=err, oracle_iq_tgc_max_error=err_tgc,
                    raw_gate_snr_db=(10*np.log10(powers/sigma**2)).tolist(),
                    noise_seed=20261030+case_index,
                    notes=dict(positive_velocity='Away from the fixed transducer',
                               time_units='seconds, measured from transmit start',
                               gate_indices='MATLAB indices are 1-based',
                               gate_columns='Six gates, one for each layer; selected gate is column 3',
                               ambiguity='Single-lag data cannot identify the true ambiguity interval on its own',
                               exact_vs_small='Compare the usual small-v output with truth.velocity_small_v_expected_mps',
                               rf_model='Each scattered REAL RF echo is generated from exact moving-target propagation, then summed',
                               tgc='Applied to summed RF including noise, before FIR and decimation'))
        save(name, arrays, meta)
        print(name, 'main noisy TGC speed:', reference['noisy_tgc'][2]['velocity_exact'])


def tgc_basis(gate, tr, ranges, velocities, onset, end, gain, cfg):
    """Exact receive-gain/FIR map; matches the real-RF path before Monte Carlo."""
    h = lowpass_kernel(params(cfg)); middle = (len(h)-1)//2
    n = (gate*cfg['decimation'])[:, None]; k = np.arange(len(h))-middle
    n0 = np.searchsorted(tr, onset)[None, :]; n1 = np.searchsorted(tr, end)[None, :]
    lo = np.clip(n-n1+1+middle, 0, len(h)); hi = np.clip(n-n0+middle+1, 0, len(h))
    pos = np.zeros((len(gate), len(ranges)), complex); neg = np.zeros_like(pos)
    for v in np.unique(velocities):
        cols = np.flatnonzero(velocities == v)
        alpha = (cfg['c']-v)/(cfg['c']+v); beta = 2*ranges[cols]/(cfg['c']+v)
        for freq, out, sgn in [(cfg['f0']*(alpha-1), pos, -1), (-cfg['f0']*(alpha+1), neg, 1)]:
            w = h[None, :]*gain[n-k[None, :]]*np.exp(-2j*np.pi*freq*k/cfg['fs_rf'])[None, :]
            prefix = np.column_stack([np.zeros(len(gate), complex), np.cumsum(w, axis=1)])
            rows = np.arange(len(gate))[:, None]
            overlap = prefix[rows, hi[:, cols]]-prefix[rows, lo[:, cols]]
            out[:, cols] = overlap*np.exp(2j*np.pi*freq*n/cfg['fs_rf'])*np.exp(sgn*2j*np.pi*cfg['f0']*beta)[None, :]
    return pos, neg


def continuous_ensemble():
    source = BASE/'distributed_scene'
    cfg = json.loads((source/'distributed_results.json').read_text(encoding='utf-8'))['configuration']
    with np.load(source/'distributed_signals.npz') as f:
        d = {k: f[k].copy() for k in f.files}
    gate = d['gate_indices_python']; L = round(cfg['lag_seconds']*cfg['fs_rf']/cfg['decimation'])
    pos, neg = tgc_basis(gate, d['t_rf_s'], d['scatterer_range0_m'], d['scatterer_velocity_mps'],
                         d['scatterer_echo_start_s'], d['scatterer_echo_end_s'], d['tgc_gain_rf'], cfg)
    amp = d['scatterer_amplitude_complex']
    error = float(np.max(abs(pos@amp+neg@np.conj(amp)-d['iq_gate_tgc'])))
    assert error < 1e-10
    r = d['scatterer_range0_m']
    std = cfg['amplitude_scale']*(cfg['reference_range_m']/r)**2*np.exp(-2*cfg['absorption_np_per_m']*r)/np.sqrt(cfg['scatterers_per_zone'])
    trials = 8192
    rng = np.random.default_rng(20261031)
    x = np.empty((len(gate), trials), complex)
    for start in range(0, trials, 256):
        count = min(256, trials-start)
        a = (rng.normal(size=(len(r), count))+1j*rng.normal(size=(len(r), count)))*std[:, None]/np.sqrt(2)
        x[:, start:start+count] = pos@a+neg@np.conj(a)
    corr = np.mean(x[L:]*np.conj(x[:-L]), axis=0)
    phase_mean = float(np.angle(np.mean(corr)))
    fd_mean = phase_mean/(2*np.pi*cfg['lag_seconds'])
    v_mean = -cfg['c']*fd_mean/(2*cfg['f0']+fd_mean)
    q_expected = (pos[L:]*np.conj(pos[:-L])+neg[L:]*np.conj(neg[:-L]))@std**2
    expected_corr = np.mean(q_expected)
    fd_expected = np.angle(expected_corr)/(2*np.pi*cfg['lag_seconds'])
    v_expected = -cfg['c']*fd_expected/(2*cfg['f0']+fd_expected)
    v_single = -cfg['c']*(np.angle(corr)/(2*np.pi*cfg['lag_seconds']))/(2*cfg['f0']+np.angle(corr)/(2*np.pi*cfg['lag_seconds']))
    tolerance = .01  # Fixed before inspecting this ensemble's result.
    assert abs(v_mean-cfg['uniform_velocity_away']) < tolerance
    arrays = dict(iq_gates_tgc=x, t_gate_s=d['t_gate_s'], ensemble_correlations=corr,
                  ensemble_velocity_exact_mps=v_single,
                  expected_correlation=np.asarray([expected_corr]),
                  example_t_rf_s=d['t_rf_s'], example_rx_rf_clean=d['rx_rf_clean'],
                  example_rx_rf_clean_tgc=d['rx_rf_clean_tgc'], example_iq_gate_tgc=d['iq_gate_tgc'],
                  example_gate_indices_matlab=gate+1, tgc_gain_rf=d['tgc_gain_rf'],
                  fir_coefficients=lowpass_kernel(params(cfg)), scatterer_range0_m=r,
                  scatterer_amplitude_std=std)
    meta = dict(params=cfg, truth=dict(velocity_mps=cfg['uniform_velocity_away'],
                                      tolerance_mps=tolerance, realizations=trials),
                reference=dict(ensemble_velocity_exact_mps=float(v_mean),
                               expected_velocity_exact_mps=float(v_expected),
                               single_estimate_std_mps=float(v_single.std(ddof=1)),
                               exact_basis_vs_full_rf_max_error=error),
                notes=dict(input='180 x 8192 I/Q gates, columns are independent scattering realizations',
                           use='Correlate each column, average COMPLEX correlations, then convert phase to velocity',
                           model='Continuous scatterers from 3 to 21 m, TGC and same FIR as full RF example',
                           independence='Artificial independent scattering realizations; no real ping-to-ping evolution model',
                           rf_scope='One full-RF example validates the exact linear map; not 8192 stored RF records',
                           noise='No added receiver noise; tests random scattering error and ensemble averaging'),
                seed=20261031)
    save('continuous_ensemble_validation', arrays, meta)
    print('continuous ensemble:', trials, 'estimated:', v_mean, 'expected:', v_expected)


if __name__ == '__main__':
    HERE.mkdir(parents=True, exist_ok=True)
    controlled_cases()
    continuous_ensemble()

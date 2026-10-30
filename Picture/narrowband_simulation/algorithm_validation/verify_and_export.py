"""Independent SciPy RF processing, ground-truth checks, and verified MAT export.

Requires NumPy + SciPy; does NOT import the simulation's demodulator/estimator.
"""
from pathlib import Path
import json
import csv
import numpy as np
from scipy.io import savemat, loadmat
from scipy.signal import firwin, fftconvolve

HERE = Path(__file__).resolve().parent
NAMES = ['narrowband_validation', 'narrowband_boundary_valid', 'narrowband_boundary_alias']


def matlabify(x):
    if isinstance(x, dict):
        return {k: matlabify(v) for k, v in x.items()}
    if isinstance(x, list):
        if all(isinstance(v, (int, float, bool)) for v in x):
            return np.asarray(x)
        return np.asarray([matlabify(v) for v in x], dtype=object)
    return x


def export(name, data, meta, verification):
    payload = dict(data)
    payload.update({k: matlabify(v) for k, v in meta.items()})
    payload['verification'] = matlabify(verification)
    target = HERE/(name+'.mat')
    savemat(target, payload, format='5', long_field_names=True, do_compression=True, oned_as='column')
    check = loadmat(target)
    for key, value in data.items():
        wanted = value.reshape(-1, 1) if value.ndim == 1 else value
        np.testing.assert_array_equal(check[key], wanted)
    print('Exported, every array read back:', target.name, target.stat().st_size, 'bytes')


def load(name):
    meta = json.loads((HERE/(name+'.json')).read_text(encoding='utf-8'))
    with np.load(HERE/(name+'.npz')) as f:
        data = {k: f[k].copy() for k in f.files}
    return data, meta


def main():
    rows = []; summaries = {}
    for name in NAMES:
        d, meta = load(name); p = meta['params']; truth = meta['truth']
        h = firwin(p['fir_taps'], p['lowpass_hz'], window='blackman', fs=p['fs_rf'])
        np.testing.assert_allclose(h, d['fir_coefficients'], rtol=1e-12, atol=1e-15)
        np.testing.assert_allclose(d['rx_rf_zones'].sum(axis=1), d['rx_rf_clean'], rtol=1e-14, atol=1e-15)
        np.testing.assert_allclose(d['rx_rf_noisy_tgc'], d['rx_rf_noisy']*d['tgc_gain_rf'], rtol=1e-14, atol=1e-15)
        gates = d['gate_indices_all_matlab']-1
        L = round(p['lag_seconds']*p['fs_rf']/p['decimation'])
        paths = [('clean', 'rx_rf_clean', 'iq_total'),
                 ('noisy', 'rx_rf_noisy', 'iq_noisy'),
                 ('clean_tgc', 'rx_rf_clean_tgc', 'iq_total_tgc'),
                 ('noisy_tgc', 'rx_rf_noisy_tgc', 'iq_noisy_tgc')]
        iq_error = 0.; case_rows = []
        for label, rf_key, iq_key in paths:
            mixed = 2*d[rf_key]*np.exp(-2j*np.pi*p['f0']*d['t_rf_s'])
            iq = fftconvolve(mixed, h, mode='same')[::p['decimation']]
            iq_error = max(iq_error, float(np.max(abs(iq-d[iq_key]))))
            np.testing.assert_allclose(iq, d[iq_key], rtol=1e-10, atol=1e-12)
            x = iq[gates]
            corr = np.sum(x[L:]*np.conj(x[:-L]), axis=0)/(len(x)-L)
            phase = np.arctan2(corr.imag, corr.real)
            fd = phase/(2*np.pi*p['lag_seconds'])
            v_small = -p['c']*fd/(2*p['f0'])
            v_exact = -p['c']*fd/(2*p['f0']+fd)
            tol = p['noisy_tolerance_mps'] if 'noisy' in label else p['clean_tolerance_mps']
            err_small = abs(v_small-np.asarray(truth['velocity_small_v_expected_mps']))
            err_exact = abs(v_exact-np.asarray(truth['velocity_exact_from_wrapped_expected_mps']))
            assert np.all(err_small < tol) and np.all(err_exact < tol), (name, label, err_exact)
            ambiguous = np.asarray(truth['is_ambiguous'])
            # Physical truth is checked only when the single-lag result is unambiguous.
            assert np.all(abs(v_exact[~ambiguous]-np.asarray(truth['velocity_mps'])[~ambiguous]) < tol)
            if label == 'clean':
                # Independent phase-slope method also checks sign, frequency and time units.
                for z in range(6):
                    slope = np.polyfit(d['t_gates_s'][:, z], np.unwrap(np.angle(x[:, z])), 1)[0]
                    assert abs(slope/(2*np.pi)-truth['fd_true_hz'][z]) < .002
            for z in range(6):
                row = dict(dataset=name, path=label, zone=z+1,
                           truth_velocity_mps=truth['velocity_mps'][z],
                           estimated_velocity_exact_mps=float(v_exact[z]),
                           estimated_velocity_small_mps=float(v_small[z]),
                           expected_single_lag_exact_mps=truth['velocity_exact_from_wrapped_expected_mps'][z],
                           model_error_mps=float(err_exact[z]), tolerance_mps=tol,
                           status='EXPECTED_ALIAS' if ambiguous[z] else 'PASS')
                rows.append(row); case_rows.append(row)
        verification = dict(independent_scipy_iq_max_error=iq_error,
                            clean_max_model_error_mps=max(r['model_error_mps'] for r in case_rows if 'noisy' not in r['path']),
                            noisy_max_model_error_mps=max(r['model_error_mps'] for r in case_rows if 'noisy' in r['path']),
                            numeric_checks_passed=len(case_rows), results=case_rows,
                            note='EXPECTED_ALIAS checks the folded result; it does not validate recovery of true ambiguous velocity')
        export(name, d, meta, verification)
        summaries[name] = verification
    name = 'continuous_ensemble_validation'
    d, meta = load(name); p = meta['params']
    L = round(p['lag_seconds']*p['fs_rf']/p['decimation'])
    x = d['iq_gates_tgc']
    corr = np.sum(x[L:]*np.conj(x[:-L]), axis=0)/(len(x)-L)
    np.testing.assert_allclose(corr, d['ensemble_correlations'], rtol=1e-12, atol=1e-14)
    h = firwin(p['fir_taps'], p['lowpass_hz'], window='blackman', fs=p['fs_rf'])
    mixed = 2*d['example_rx_rf_clean_tgc']*np.exp(-2j*np.pi*p['f0']*d['example_t_rf_s'])
    example_iq = fftconvolve(mixed, h, mode='same')[::p['decimation']]
    gate = d['example_gate_indices_matlab']-1
    np.testing.assert_allclose(example_iq[gate], d['example_iq_gate_tgc'], rtol=1e-10, atol=1e-12)
    mean_corr = corr.mean()
    fd = np.arctan2(mean_corr.imag, mean_corr.real)/(2*np.pi*p['lag_seconds'])
    speed = -p['c']*fd/(2*p['f0']+fd)
    error = abs(speed-meta['truth']['velocity_mps'])
    assert error < meta['truth']['tolerance_mps']
    verification = dict(status='PASS', velocity_exact_mps=float(speed), absolute_error_mps=float(error),
                        tolerance_mps=meta['truth']['tolerance_mps'], realizations=x.shape[1],
                        example_full_rf_iq_max_error=float(np.max(abs(example_iq[gate]-d['example_iq_gate_tgc']))))
    export(name, d, meta, verification)
    summaries[name] = verification
    report = dict(status='PASS', controlled_numeric_checks=len(rows), datasets=summaries,
                  matlab_executed=False, independent_processor='SciPy firwin + fftconvolve; direct complex-product estimator')
    (HERE/'validation_report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    with (HERE/'validation_results.csv').open('w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    print('All 72 controlled comparisons and continuous-ensemble test passed.')
    print('Main profile (noisy + TGC):')
    for row in rows:
        if row['dataset'] == NAMES[0] and row['path'] == 'noisy_tgc':
            print(row['zone'], row['truth_velocity_mps'], row['estimated_velocity_exact_mps'], row['status'])
    print('Continuous ensemble:', verification)


if __name__ == '__main__':
    main()

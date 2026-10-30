"""Export saved simulation arrays and metadata to MATLAB MAT format.

Run after narrowband_demo.py. Requires numpy and scipy.
All signal/time vectors are MATLAB column vectors.
"""
from pathlib import Path
import json

import numpy as np
from scipy.io import savemat, loadmat


def export():
    folder = Path(__file__).resolve().parent
    results = json.loads((folder / "results.json").read_text(encoding="utf-8"))
    with np.load(folder / "signals.npz", allow_pickle=False) as z:
        data = {name: z[name].copy() for name in z.files}

    est = results["baseline"]
    n = len(data["gate_iq"])
    lag = est["lag_samples"]
    cloud_columns = ["lag_ms", "pairs", "mean_velocity", "std_velocity", "mean_coherence"]
    validation_columns = ["true_velocity", "true_doppler_hz", "measured_doppler_hz",
                          "measured_velocity", "ambiguity_expected"]
    payload = {
        "params": results["parameters"],
        "truth": results["truth"],
        "estimate": est,
        "fs_iq_hz": results["parameters"]["fs_rf"] / results["parameters"]["decimation"],
        "t_rf_s": data["t_rf_s"],
        "tx_rf": data["tx_rf"],
        "rx_rf": data["rx_rf"],
        # The original NPZ name receiver_noise_rf contains signal + noise.
        "rx_rf_noisy": data["receiver_noise_rf"],
        "noise_rf": data["receiver_noise_rf"] - data["rx_rf"],
        "t_iq_s": data["t_iq_s"],
        "iq_full": data["iq"],
        "t_gate_s": data["gate_time_s"],
        "iq_gate": data["gate_iq"],
        "lag_products": data["lag_products"],
        "gate_indices_matlab": np.arange(est["gate_first_index"] + 1,
                                          est["gate_last_exclusive"] + 1, dtype=np.int64),
        "pair_indices_matlab": np.column_stack((np.arange(1, n-lag+1),
                                                 np.arange(lag+1, n+1))).astype(np.int64),
        "noisy_rf_estimate": results["noisy_rf"],
        "iq_noise_stats": results["iq_noise_monte_carlo"],
        "cloud_lag_table": np.array([[row[k] for k in cloud_columns]
                                      for row in results["cloud_monte_carlo"]]),
        "cloud_lag_columns": np.array(cloud_columns, dtype=object),
        "validation_table": np.array([[row[k] for k in validation_columns]
                                       for row in results["validation"]], dtype=float),
        "validation_columns": np.array(validation_columns, dtype=object),
        "notes": {
            "positive_velocity": "Away from fixed transducer; m/s",
            "time_unit": "s; all time arrays are seconds",
            "signals": "RF is real; I/Q and lag_products are complex; amplitudes arbitrary",
            "indices": "gate_indices_matlab and pair_indices_matlab use 1-based MATLAB indices",
            "estimate_indices": "estimate.gate_first_index and gate_last_exclusive preserve Python 0-based/end-exclusive convention",
            "filter": "Offline symmetric FIR delay compensated; 0.5 ms guard before gate",
            "model": "Single-beam teaching simulation, not manufacturer firmware or hardware validation",
        },
    }
    target = folder / "narrowband_simulation.mat"
    savemat(target, payload, format="5", long_field_names=True,
            do_compression=True, oned_as="column")

    # Verify the saved file, including complex components and MATLAB dimensions.
    loaded = loadmat(target)
    for name in ["t_rf_s", "tx_rf", "rx_rf", "rx_rf_noisy", "noise_rf", "t_iq_s",
                 "iq_full", "t_gate_s", "iq_gate", "lag_products", "gate_indices_matlab"]:
        np.testing.assert_array_equal(loaded[name], np.asarray(payload[name]).reshape(-1, 1))
    np.testing.assert_array_equal(loaded["pair_indices_matlab"], payload["pair_indices_matlab"])
    assert loaded["iq_gate"].dtype.kind == "c"
    meta = loadmat(target, simplify_cells=True)
    assert meta["params"]["velocity_away"] == results["parameters"]["velocity_away"]
    assert meta["estimate"]["phase_rad"] == est["phase_rad"]
    print(f"Saved and verified: {target}")
    for name in ["tx_rf", "rx_rf", "iq_full", "iq_gate", "lag_products", "pair_indices_matlab"]:
        print(name, loaded[name].shape, str(loaded[name].dtype))
    print("Bytes:", target.stat().st_size)


if __name__ == "__main__":
    export()

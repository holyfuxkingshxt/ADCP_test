"""Export all distributed-scene RF, I/Q, scatterers and diagnostics to MAT."""
from pathlib import Path
import json
import numpy as np
from scipy.io import savemat,loadmat


def matlab_struct(x):
    if isinstance(x,dict):
        return {k:matlab_struct(v) for k,v in x.items()}
    if isinstance(x,list):
        if all(isinstance(v,(int,float,bool)) for v in x):
            return np.asarray(x)
        return np.asarray(x,dtype=object)
    return x


def main():
    folder=Path(__file__).resolve().parent
    stats=json.loads((folder/'distributed_results.json').read_text(encoding='utf-8'))
    with np.load(folder/'distributed_signals.npz') as z:
        data={k:z[k].copy() for k in z.files}
    payload=dict(data)
    payload['params']=matlab_struct(stats['configuration'])
    payload['diagnostics']=matlab_struct({k:v for k,v in stats.items() if k!='configuration'})
    payload['notes']={
        'rf_dimensions':'time x range zone; six zone columns',
        'iq_dimensions':'time x range zone; six zone columns',
        'column_vectors':'All 1D arrays are saved as MATLAB column vectors',
        'zone_ids':'MATLAB 1-based zone IDs; edges are params.zone_edges_m',
        'velocity_sign':'Positive away from fixed transducer; units m/s',
        'time_units':'All time vectors and overlaps use seconds',
        'gate_indices':'Use gate_indices_matlab directly; *_python uses zero-based indices',
        'spatial_support':'A receive gate contains finite-pulse echoes from several physical range zones',
        'tgc_rf':'rx_rf_noisy is raw; rx_rf_noisy_tgc = tgc_gain_rf .* rx_rf_noisy',
        'tgc_iq':'*_tgc I/Q is recomputed from gain-compensated RF BEFORE FIR, not scaled after decimation',
        'tgc_noise':'noise_rf_tgc = tgc_gain_rf .* noise_rf; TGC amplifies signal and noise together',
    }
    target=folder/'distributed_narrowband_simulation.mat'
    savemat(target,payload,format='5',long_field_names=True,do_compression=True,oned_as='column')
    check=loadmat(target)
    for key,original in data.items():
        wanted=original.reshape(-1,1) if original.ndim==1 else original
        np.testing.assert_array_equal(check[key],wanted)
    assert check['rx_rf_zones'].shape==(81600,6)
    assert check['iq_gate_zones'].shape==(180,6)
    assert check['iq_gate'].dtype.kind=='c'
    assert check['rx_rf_noisy_tgc'].shape==(81600,1)
    assert check['iq_gate_noisy_tgc'].shape==(180,1)
    print('Saved and all arrays verified:',target)
    print('MAT file bytes:',target.stat().st_size)
    for key in ['rx_rf_clean','rx_rf_zones','iq_total','iq_zones','iq_gate_zones','lag_products',
                'tgc_gain_rf','rx_rf_noisy_tgc','iq_noisy_tgc','iq_gate_noisy_tgc']:
        print(key,check[key].shape,check[key].dtype)


if __name__=='__main__':
    main()

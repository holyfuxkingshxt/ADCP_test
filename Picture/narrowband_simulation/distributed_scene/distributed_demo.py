"""Full real-RF distributed narrowband scene, then estimate ONE receive gate.

Run with the bundled NumPy/ReportLab Python. Export to MAT with the companion
export_distributed_mat.py using a Python containing NumPy and SciPy.
"""
from pathlib import Path
import sys
import json
import csv

import numpy as np
from reportlab.pdfgen import canvas

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from narrowband_demo import Parameters, demodulate, estimate_velocity, lowpass_kernel, panel, page_title
from add_tgc import add_tgc_data, make_tgc_plots


CONFIG = dict(
    c=1500.0, f0=300000.0, fs_rf=2400000.0, pulse_duration=0.004,
    record_duration=0.034, fir_taps=1601, lowpass_hz=6000.0, decimation=40,
    lag_seconds=0.0005, gate_start_s=0.0145, gate_duration_s=0.003,
    zone_edges_m=[3., 6., 9., 12., 15., 18., 21.],
    scatterers_per_zone=240, uniform_velocity_away=0.2,
    target_zone_matlab=3, adjacent_deeper_zone_matlab=4,
    far_zone_matlab=6, far_amplitude_multiplier=10.0,
    adjacent_deeper_shear_velocity=0.6,
    reference_range_m=10.0, amplitude_scale=0.04,
    absorption_np_per_m=0.005, input_snr_gate_db=20.0, seed=20261029,
)


def generate_scatterers(cfg):
    rng = np.random.default_rng(cfg['seed'])
    edges = np.asarray(cfg['zone_edges_m'])
    count = cfg['scatterers_per_zone']
    ranges = np.concatenate([rng.uniform(a,b,count) for a,b in zip(edges[:-1],edges[1:])])
    zones = np.repeat(np.arange(len(edges)-1),count)
    # Random complex scattering amplitude: every draw remains fixed during ping.
    scatter = (rng.normal(size=len(ranges))+1j*rng.normal(size=len(ranges))) / np.sqrt(2*count)
    # Simplified TWO-WAY point-scatterer pressure loss, evaluated at initial range.
    gain = ((cfg['reference_range_m']/ranges)**2
            * np.exp(-2*cfg['absorption_np_per_m']*ranges))
    amplitudes = cfg['amplitude_scale'] * scatter * gain
    return ranges, zones, amplitudes


def rf_scene(t, ranges, zones, amplitudes, velocities, cfg):
    nz = len(cfg['zone_edges_m'])-1
    rx = np.zeros((len(t),nz))
    onset = np.zeros(len(ranges)); end = np.zeros_like(onset)
    for i,(r0,z,a,v) in enumerate(zip(ranges,zones,amplitudes,velocities)):
        alpha=(cfg['c']-v)/(cfg['c']+v)
        beta=2*r0/(cfg['c']+v)
        onset[i]=beta/alpha
        end[i]=onset[i]+cfg['pulse_duration']/alpha
        n0=int(np.searchsorted(t,onset[i])); n1=int(np.searchsorted(t,end[i]))
        phase=2*np.pi*cfg['f0']*(alpha*t[n0:n1]-beta)
        # Re{a * exp(j*phase)}. Sum every scatterer's REAL waveform BEFORE I/Q.
        rx[n0:n1,z] += a.real*np.cos(phase)-a.imag*np.sin(phase)
    return rx,onset,end


def geometric_weights(onset,end,gate_start,gate_end,lag):
    # Echo-energy gate overlap and same-scatterer lag-product overlap.
    # These geometric expectations exclude FIR spreading and coherent cross terms.
    w_echo=np.maximum(0.,np.minimum(end,gate_end)-np.maximum(onset,gate_start))
    w_pair=np.maximum(0.,np.minimum(end-lag,gate_end-lag)-np.maximum(onset,gate_start))
    return w_echo,w_pair


def estimate(x, cfg):
    p=Parameters(c=cfg['c'],f0=cfg['f0'])
    return estimate_velocity(x,cfg['fs_rf']/cfg['decimation'],cfg['lag_seconds'],p)


def expected_phase_velocity(amplitudes,velocities,pair_weights,cfg):
    fd=-2*cfg['f0']*velocities/(cfg['c']+velocities)
    value=np.sum(abs(amplitudes)**2*pair_weights*np.exp(2j*np.pi*fd*cfg['lag_seconds']))
    phase=float(np.angle(value))
    speed=-cfg['c']*phase/(4*np.pi*cfg['f0']*cfg['lag_seconds'])
    return dict(phase_rad=phase,velocity_small_v=float(speed),
                interpretation='Ensemble autocovariance expectation without random cross terms; geometric, unfiltered')


def exact_gate_bases(gate,tr,ranges,velocities,onset,end,p,cfg):
    """Exact linear map from scatterer amplitudes to selected filtered I/Q.

    Uses cumulative FIR sums for rectangular, time-scaled echoes. Includes both
    the desired baseband component and the residual real-RF image. This avoids
    repeating all carrier calculations for independent scattering realizations.
    """
    n=(gate*cfg['decimation'])[:,None]
    n0=np.searchsorted(tr,onset)[None,:]; n1=np.searchsorted(tr,end)[None,:]
    h=lowpass_kernel(p); middle=(len(h)-1)//2
    lo=np.clip(n-n1+1+middle,0,len(h))
    hi=np.clip(n-n0+middle+1,0,len(h))
    pos=np.zeros((len(gate),len(ranges)),complex)
    neg=np.zeros_like(pos)
    signed_k=np.arange(len(h))-middle
    for velocity in np.unique(velocities):
        cols=np.flatnonzero(velocities==velocity)
        alpha=(cfg['c']-velocity)/(cfg['c']+velocity)
        beta=2*ranges[cols]/(cfg['c']+velocity)
        for freq,output,phase_sign in [(cfg['f0']*(alpha-1),pos,-1),
                                       (-cfg['f0']*(alpha+1),neg,1)]:
            prefix=np.concatenate(([0j],np.cumsum(h*np.exp(-2j*np.pi*freq*signed_k/cfg['fs_rf']))))
            overlap=prefix[hi[:,cols]]-prefix[lo[:,cols]]
            output[:,cols]=(overlap*np.exp(2j*np.pi*freq*(n/cfg['fs_rf']))
                            *np.exp(phase_sign*2j*np.pi*cfg['f0']*beta)[None,:])
    return pos,neg


def ensemble_check(pos,neg,ranges,amp,gate_iq,cfg):
    mapped=pos@amp+neg@np.conj(amp)
    error=float(np.max(abs(mapped-gate_iq)))
    assert error<max(float(np.max(abs(gate_iq)))*1e-9,1e-13)
    gain=((cfg['reference_range_m']/ranges)**2
          *np.exp(-2*cfg['absorption_np_per_m']*ranges))
    scatter_std=cfg['amplitude_scale']*gain/np.sqrt(cfg['scatterers_per_zone'])
    L=round(cfg['lag_seconds']*cfg['fs_rf']/cfg['decimation'])
    # Zero-mean circular independent scatterer amplitudes remove all i!=j terms.
    expected_q=((pos[L:]*np.conj(pos[:-L])+neg[L:]*np.conj(neg[:-L]))@scatter_std**2)
    expected_corr=np.mean(expected_q)
    factor=-cfg['c']/(4*np.pi*cfg['f0']*cfg['lag_seconds'])
    expected_speed=float(factor*np.angle(expected_corr))
    rng=np.random.default_rng(cfg['seed']+1900)
    trials=1000
    random_amp=(rng.normal(size=(len(ranges),trials))+1j*rng.normal(size=(len(ranges),trials)))
    random_amp*=scatter_std[:,None]/np.sqrt(2)
    x=pos@random_amp+neg@np.conj(random_amp)
    correlations=np.mean(x[L:]*np.conj(x[:-L]),axis=0)
    speeds=factor*np.angle(correlations)
    rows=[]
    for count in [10,100,1000]:
        rows.append(dict(realizations=count,
                         velocity_from_mean_correlation=float(factor*np.angle(np.mean(correlations[:count])))))
    return dict(exact_basis_vs_rf_max_error=error,
                expected_velocity_from_full_filter=float(expected_speed),
                monte_carlo_realizations=trials,
                single_ping_velocity_mean=float(speeds.mean()),
                single_ping_velocity_std=float(speeds.std(ddof=1)),
                convergence=rows,
                note='Independent scattering realizations at the same ranges; not a physical ping-time evolution model'),speeds,correlations


def make_plots(arrays,stats,cfg):
    pdf=canvas.Canvas(str(HERE/'distributed_scene_figures.pdf'),pagesize=(1040,780))
    pdf.setTitle('Distributed narrowband echoes: all ranges, one receive gate')
    boxes=[(25,390,490,320),(525,390,490,320),(25,50,490,320),(525,50,490,320)]
    tr=arrays['t_rf_s']; tb=arrays['t_iq_s']; gate=arrays['gate_indices_python']
    iq=arrays['iq_total']; iz=arrays['iq_zones']; target=cfg['target_zone_matlab']-1
    block=240 # 0.1 ms, computes RMS instead of alias-prone waveform subsampling.
    usable=len(tr)//block*block
    times=tr[:usable].reshape(-1,block).mean(axis=1)*1000
    clean_rms=np.sqrt(np.mean(arrays['rx_rf_clean'][:usable].reshape(-1,block)**2,axis=1))
    noisy_rms=np.sqrt(np.mean(arrays['rx_rf_noisy'][:usable].reshape(-1,block)**2,axis=1))
    page_title(pdf,'1 / Full received record, then one gate',
               '1440 scatterers in 3-21 m | 4 ms narrowband pulse | Real RF is summed before demodulation',1)
    panel(pdf,boxes[0],[(times,clean_rms,'Clean RF'),(times,noisy_rms,'Noisy RF')],
          'A. Whole-record RF amplitude (0.1 ms RMS)','Receive time from TX start (ms)','RF RMS amplitude')
    panel(pdf,boxes[1],[(tb*1000,abs(iq),'All ranges'),(tb[gate]*1000,abs(iq[gate]),'Selected gate')],
          'B. Full I/Q envelope and selected gate','Receive time (ms)','I/Q magnitude')
    series=[]
    for z in [target-1,target,target+1]:
        a=cfg['zone_edges_m'][z]; b=cfg['zone_edges_m'][z+1]
        series.append((tb[gate]*1000,abs(iz[gate,z]),f'{a:g}-{b:g} m'))
    series.append((tb[gate]*1000,abs(iq[gate]),'Total'))
    panel(pdf,boxes[2],series,'C. Several range zones enter ONE gate','Receive time (ms)','I/Q magnitude')
    panel(pdf,boxes[3],[(tb[gate]*1000,abs(iq[gate]),'Coherent total'),
                       (tb[gate]*1000,np.sqrt(np.sum(abs(iz[gate])**2,axis=1)),'Incoherent RSS'),
                       (tb[gate]*1000,abs(iz[gate,target]),'Target zone only')],
          'D. Contributions can reinforce or cancel','Receive time (ms)','I/Q magnitude')
    pdf.showPage()
    page_title(pdf,'2 / Which ranges influence this gate?',
               'Gate = 14.5-17.5 ms | Nominal range center = 10.5 m | No multipath or previous-ping echoes',2)
    r=np.linspace(3,21,600); v=cfg['uniform_velocity_away']; c=cfg['c']
    on=2*r/(c-v); dur=cfg['pulse_duration']*(c+v)/(c-v)
    we,wp=geometric_weights(on,on+dur,cfg['gate_start_s'],cfg['gate_start_s']+cfg['gate_duration_s'],cfg['lag_seconds'])
    panel(pdf,boxes[0],[(r,we*1000,'Gate overlap'),(r,wp*1000,'Same-particle pair overlap')],
          'A. Geometric contribution vs initial range','Initial scatterer range (m)','Overlap duration (ms)')
    energy=np.asarray(stats['zone_energy_in_gate'])
    series=[]
    for z,e in enumerate(energy):
        series.append((np.array([z+1-.32,z+1-.32,z+1+.32,z+1+.32]),np.array([0,e,e,0]),f'Z{z+1}'))
    panel(pdf,boxes[1],series,'B. Per-zone gate energy (before cross terms)','Zone index (see zone table)','Sum of |I/Q| squared')
    near=iq-iz[:,cfg['adjacent_deeper_zone_matlab']-1]
    panel(pdf,boxes[2],[(tb[gate]*1000,abs(iq[gate]),'All zones'),
                       (tb[gate]*1000,abs(near[gate]),'Remove adjacent deeper')],
          'C. Remove the 12-15 m zone','Receive time (ms)','I/Q magnitude')
    panel(pdf,boxes[3],[(tb[gate]*1000,abs(iq[gate]),'Original'),
                       (tb[gate]*1000,abs(arrays['iq_far_boosted'][gate]),'18-21 m amplitude x10')],
          'D. A distant direct echo arrives too late','Receive time (ms)','I/Q magnitude')
    pdf.save()


def main():
    cfg=CONFIG.copy(); HERE.mkdir(parents=True,exist_ok=True)
    ranges,zones,amp=generate_scatterers(cfg)
    velocities=np.full(len(ranges),cfg['uniform_velocity_away'])
    tr=np.arange(round(cfg['record_duration']*cfg['fs_rf']))/cfg['fs_rf']
    rxz,onset,end=rf_scene(tr,ranges,zones,amp,velocities,cfg)
    rx=rxz.sum(axis=1)
    tx=(tr<cfg['pulse_duration'])*np.cos(2*np.pi*cfg['f0']*tr)
    p=Parameters(c=cfg['c'],f0=cfg['f0'],fs_rf=cfg['fs_rf'],
                 fir_taps=cfg['fir_taps'],lowpass_hz=cfg['lowpass_hz'],decimation=cfg['decimation'])
    tb,iq=demodulate(tr,rx,p)
    iz=np.column_stack([demodulate(tr,rxz[:,z],p)[1] for z in range(rxz.shape[1])])
    fsbb=cfg['fs_rf']/cfg['decimation']
    first=int(round(cfg['gate_start_s']*fsbb)); count=int(round(cfg['gate_duration_s']*fsbb))
    gate=np.arange(first,first+count)
    assert tb[gate[-1]]<cfg['gate_start_s']+cfg['gate_duration_s']
    rng=np.random.default_rng(cfg['seed']+900)
    rf_gate=(tr>=cfg['gate_start_s'])&(tr<cfg['gate_start_s']+cfg['gate_duration_s'])
    noise_sigma=np.sqrt(np.mean(rx[rf_gate]**2)/10**(cfg['input_snr_gate_db']/10))
    noise=rng.normal(0,noise_sigma,len(tr)); noisy=rx+noise
    _,iqnoisy=demodulate(tr,noisy,p)
    far=cfg['far_zone_matlab']-1; deeper=cfg['adjacent_deeper_zone_matlab']-1; target=cfg['target_zone_matlab']-1
    boosted_rf=rx+(cfg['far_amplitude_multiplier']-1)*rxz[:,far]
    _,boosted_iq=demodulate(tr,boosted_rf,p)
    changed_velocity=velocities.copy()
    changed_velocity[zones==deeper]=cfg['adjacent_deeper_shear_velocity']
    shear_z,_,_=rf_scene(tr,ranges[zones==deeper],np.full(np.sum(zones==deeper),deeper),
                         amp[zones==deeper],changed_velocity[zones==deeper],cfg)
    shear_rf=rx-rxz[:,deeper]+shear_z[:,deeper]
    _,shear_iq=demodulate(tr,shear_rf,p)
    we,wp=geometric_weights(onset,end,cfg['gate_start_s'],
                            cfg['gate_start_s']+cfg['gate_duration_s'],cfg['lag_seconds'])
    uniform_est,q=estimate(iq[gate],cfg); noisy_est,_=estimate(iqnoisy[gate],cfg)
    target_est,_=estimate(iz[gate,target],cfg)
    removed_est,_=estimate((iq-iz[:,deeper])[gate],cfg)
    shear_est,_=estimate(shear_iq[gate],cfg)
    boosted_est,_=estimate(boosted_iq[gate],cfg)
    energy=np.sum(abs(iz[gate])**2,axis=0)
    total_energy=float(np.sum(abs(iq[gate])**2))
    linearity_error=float(np.max(abs(iq-iz.sum(axis=1))))
    far_difference=float(np.max(abs(boosted_iq[gate]-iq[gate])))
    typical=float(np.max(abs(iq[gate])))
    assert np.max(abs(rx-rxz.sum(axis=1)))==0
    assert linearity_error<max(typical*1e-10,1e-13)
    assert far_difference<max(typical*1e-10,1e-13)
    assert float(np.sum(we[zones==far]))==0.
    assert np.sum(we[zones==deeper])>0
    c=cfg['c']; v=cfg['uniform_velocity_away']; start=cfg['gate_start_s']; finish=start+cfg['gate_duration_s']
    support_lo=((c-v)*start-(c+v)*cfg['pulse_duration'])/2
    support_hi=(c-v)*finish/2
    stats=dict(
        configuration=cfg,rf_samples=len(tr),iq_samples=len(tb),scatterer_count=len(ranges),
        selected_gate_samples=count,selected_gate_pairs=len(q),gate_start_index_python=first,
        gate_end_exclusive_python=first+count,
        nominal_gate_center_range_m=c*((start+finish)/2-cfg['pulse_duration']/2)/2,
        unfiltered_echo_support_initial_range_m=[support_lo,support_hi],
        fir_half_support_s=(cfg['fir_taps']-1)/(2*cfg['fs_rf']),
        zone_energy_in_gate=energy.tolist(),total_gate_energy=total_energy,
        incoherent_sum_zone_energy=float(energy.sum()),
        cross_zone_interference_energy=total_energy-float(energy.sum()),
        fraction_gate_samples_total_above_target_only=float(np.mean(abs(iq[gate])>abs(iz[gate,target]))),
        fraction_deeper_zone_increases_magnitude=float(np.mean(abs(iq[gate])>abs((iq-iz[:,deeper])[gate]))),
        iq_linearity_max_error=linearity_error,far_boost_gate_max_error=far_difference,
        estimates=dict(uniform_total=uniform_est,noisy_total=noisy_est,target_zone_only=target_est,
                       without_adjacent_deeper=removed_est,adjacent_deeper_velocity_changed=shear_est,
                       far_zone_amplitude_x10=boosted_est),
        expected_uniform=expected_phase_velocity(amp,velocities,wp,cfg),
        expected_shear=expected_phase_velocity(amp,changed_velocity,wp,cfg),
        limitations=[
            'Single beam with point scatterers and fixed c; no hardware blanking.',
            'No multipath, bottom boundary or previous-ping echoes.',
            'Point-scatterer amplitudes include assumed r^-2 two-way spreading and absorption, fixed at initial range.',
            'Uniform number of points per metre; no full 3-D beam-volume or calibrated backscatter model.',
            'One speckle realization; estimates can fluctuate around ensemble expectations.',
            'Gate is a time window, not an isolated 9-12 m slab; finite pulse and FIR spread its spatial sensitivity.',
            'Coherence is overlap-normalized, not a vendor correlation-count field.',
        ])
    basis_pos,basis_neg=exact_gate_bases(gate,tr,ranges,velocities,onset,end,p,cfg)
    ens_stats,ens_speeds,ens_corr=ensemble_check(basis_pos,basis_neg,ranges,amp,iq[gate],cfg)
    stats['ensemble_validation']=ens_stats
    arrays=dict(t_rf_s=tr,tx_rf=tx,rx_rf_clean=rx,rx_rf_noisy=noisy,noise_rf=noise,
                rx_rf_zones=rxz,rx_rf_shear=shear_rf,t_iq_s=tb,iq_total=iq,iq_noisy=iqnoisy,
                iq_zones=iz,iq_shear=shear_iq,iq_far_boosted=boosted_iq,
                gate_indices_python=gate,gate_indices_matlab=gate+1,t_gate_s=tb[gate],
                iq_gate=iq[gate],iq_gate_noisy=iqnoisy[gate],iq_gate_zones=iz[gate],
                lag_products=q,scatterer_range0_m=ranges,scatterer_zone_matlab=zones+1,
                scatterer_amplitude_complex=amp,scatterer_velocity_mps=velocities,
                scatterer_shear_velocity_mps=changed_velocity,scatterer_echo_start_s=onset,
                scatterer_echo_end_s=end,scatterer_gate_overlap_s=we,scatterer_pair_overlap_s=wp,
                ensemble_velocity_mps=ens_speeds,ensemble_correlations=ens_corr)
    arrays,stats=add_tgc_data(arrays,stats,cfg)
    np.savez_compressed(HERE/'distributed_signals.npz',**arrays)
    (HERE/'distributed_results.json').write_text(json.dumps(stats,indent=2),encoding='utf-8')
    with (HERE/'range_zone_contributions.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.writer(f);writer.writerow(['zone','range_start_m','range_end_m','points','direct_gate_points','direct_pair_points','filtered_iq_gate_energy'])
        for z in range(len(energy)):
            writer.writerow([z+1,cfg['zone_edges_m'][z],cfg['zone_edges_m'][z+1],int(np.sum(zones==z)),
                             int(np.sum((zones==z)&(we>0))),int(np.sum((zones==z)&(wp>0))),energy[z]])
    make_plots(arrays,stats,cfg)
    make_tgc_plots(arrays,stats,cfg)
    print(json.dumps({k:stats[k] for k in ['rf_samples','iq_samples','scatterer_count',
          'unfiltered_echo_support_initial_range_m','zone_energy_in_gate',
          'cross_zone_interference_energy','fraction_deeper_zone_increases_magnitude',
          'iq_linearity_max_error','far_boost_gate_max_error','estimates','expected_uniform','expected_shear','ensemble_validation']},indent=2))


if __name__=='__main__':
    main()

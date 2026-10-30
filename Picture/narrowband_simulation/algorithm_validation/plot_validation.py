"""Scientific verification figures made with the existing vector plot helpers."""
from pathlib import Path
import sys
import json
import numpy as np
from reportlab.pdfgen import canvas

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from narrowband_demo import panel,page_title


def main():
    with np.load(HERE/'narrowband_validation.npz') as f:
        d={k:f[k].copy() for k in f.files}
    meta=json.loads((HERE/'narrowband_validation.json').read_text(encoding='utf-8'))
    report=json.loads((HERE/'validation_report.json').read_text(encoding='utf-8'))
    pdf=canvas.Canvas(str(HERE/'validation_figures.pdf'),pagesize=(1040,780))
    pdf.setTitle('Narrowband algorithm validation: controlled truth and statistical truth')
    boxes=[(25,390,490,320),(525,390,490,320),(25,50,490,320),(525,50,490,320)]
    page_title(pdf,'1 / Known-answer full-RF validation',
               '6 controlled scattering layers | 1440 scatterers | Same 300 kHz carrier, 4 ms TX and 0.5 ms lag',1)
    block=240
    times=d['t_rf_s'].reshape(-1,block).mean(axis=1)*1000
    def rms(key):
        return np.sqrt(np.mean(d[key].reshape(-1,block)**2,axis=1))
    panel(pdf,boxes[0],[(times,rms('rx_rf_noisy'),'Raw noisy RF'),(times,rms('rx_rf_noisy_tgc'),'With TGC')],
          'A. All layers are synthesized before selecting a gate','Receive time (ms)','RF RMS amplitude')
    tt=d['t_gate_s']-d['t_gate_s'][0]
    phase=np.unwrap(np.angle(d['iq_gate_noisy_tgc'])); phase-=phase[0]
    expected=2*np.pi*meta['truth']['fd_true_hz'][2]*tt
    panel(pdf,boxes[1],[(tt*1000,expected,'Truth'),(tt*1000,phase,'Noisy + TGC')],
          'B. Selected gate: phase follows the known Doppler','Time within gate (ms)','Relative phase (rad)')
    rows=[r for r in report['datasets']['narrowband_validation']['results'] if r['path']=='noisy_tgc']
    truth=np.asarray(meta['truth']['velocity_mps'])
    estimate=np.array([r['estimated_velocity_exact_mps'] for r in rows])
    zones=np.arange(1,7)
    panel(pdf,boxes[2],[(zones,truth,'Truth'),(zones,estimate,'Noisy + TGC')],
          'C. Zero, positive and negative velocities all match','Layer index','Velocity (m/s)')
    error=np.abs(estimate-truth)*1000
    panel(pdf,boxes[3],[(zones,error,'Absolute error'),(zones,np.ones(6),'Tolerance')],
          'D. All six noisy gates pass the 1 mm/s criterion','Layer index','Absolute error (mm/s)',(1,6,0,1.1))
    pdf.showPage()
    with np.load(HERE/'continuous_ensemble_validation.npz') as f:
        x=f['iq_gates_tgc']; single=f['ensemble_velocity_exact_mps']; corr=f['ensemble_correlations']
    m=json.loads((HERE/'continuous_ensemble_validation.json').read_text(encoding='utf-8'))
    p=m['params']; vtrue=m['truth']['velocity_mps']
    page_title(pdf,'2 / Continuous scattering: verify the ensemble',
               'Original 3-21 m continuous scene | 8192 independent configurations | TGC | No added receiver noise',2)
    hist,edges=np.histogram(single,bins=70,density=True)
    xx=np.repeat(edges,2)[1:-1]; yy=np.repeat(hist,2)
    panel(pdf,boxes[0],[(xx,yy,'Single realization'),(np.array([vtrue,vtrue]),np.array([0,hist.max()]),'Truth')],
          'A. Individual estimates retain physical random error','Single-estimate velocity (m/s)','Probability density')
    counts=np.unique(np.round(np.geomspace(1,len(corr),70)).astype(int))
    def speed(r):
        fd=np.angle(r)/(2*np.pi*p['lag_seconds'])
        return -p['c']*fd/(2*p['f0']+fd)
    estimates=np.array([speed(corr[:n].mean()) for n in counts])
    panel(pdf,boxes[1],[(np.log10(counts),estimates,'Complex-mean estimate'),
                       (np.log10(counts),np.full(len(counts),vtrue),'Truth')],
          'B. Average correlations, then compute velocity','log10(number of independent configurations)','Ensemble velocity (m/s)')
    panel(pdf,boxes[2],[(np.log10(counts),abs(estimates-vtrue),'Absolute error'),
                       (np.log10(counts),np.full(len(counts),.01),'Tolerance')],
          'C. Final ensemble error is 0.00201 m/s','log10(number of independent configurations)','Absolute error (m/s)')
    left,bottom,width,height=boxes[3]
    pdf.setFont('Helvetica-Bold',13); pdf.drawString(left+35,bottom+height-22,'What each dataset validates')
    lines=[
        'Controlled layers: RF mixing, FIR, delay alignment, gates,',
        'lag, velocity sign, scale and receive-time gain.',
        '',
        'Continuous ensemble: random scattering and averaging.',
        'Fixed seeds; no selection of unusually accurate realizations.',
        '',
        'Boundary +2.4 m/s: physical velocity is recovered.',
        'Boundary +3.0 m/s: expected alias is verified and labelled.',
        'A single lag cannot identify the true ambiguity interval.',
        '',
        'The ensemble configurations are independent synthetic draws.',
        'They do not model real ping-to-ping time evolution.',
    ]
    pdf.setFont('Helvetica',10)
    for i,line in enumerate(lines):
        pdf.drawString(left+35,bottom+height-55-i*17,line)
    pdf.save()


if __name__=='__main__':
    main()

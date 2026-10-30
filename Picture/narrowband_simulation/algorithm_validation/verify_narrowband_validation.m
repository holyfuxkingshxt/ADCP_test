function results = verify_narrowband_validation(matPath)
% Recompute from REAL RF with MATLAB core functions; no stored IQ is needed.
% Exact inversion is used to compare with physical velocity truth.
% The usual small-v inversion is separately compared with its expected value.
if nargin < 1
    matPath = fullfile(fileparts(mfilename('fullpath')), 'narrowband_validation.mat');
end
S = load(matPath); p = S.params;
t = S.t_rf_s; M = double(p.fir_taps); fs = p.fs_rf;
k = (-(M-1)/2:(M-1)/2).';
a = 2*p.lowpass_hz*k/fs;
si = ones(size(a)); nz = a~=0; si(nz) = sin(pi*a(nz))./(pi*a(nz));
j = (0:M-1).';
w = 0.42-0.5*cos(2*pi*j/(M-1))+0.08*cos(4*pi*j/(M-1));
h = (2*p.lowpass_hz/fs)*si.*w; h = h/sum(h);
assert(max(abs(h-S.fir_coefficients)) < 1e-14, 'FIR coefficients differ');
g = S.gate_indices_all_matlab;
L = round(p.lag_seconds*fs/p.decimation);
rf = {S.rx_rf_clean, S.rx_rf_noisy, S.rx_rf_clean_tgc, S.rx_rf_noisy_tgc};
iqReference = {S.iq_total, S.iq_noisy, S.iq_total_tgc, S.iq_noisy_tgc};
pathNames = {'clean','noisy','clean_tgc','noisy_tgc'};
datasetCol = cell(24,1); zoneCol=zeros(24,1); truthCol=zeros(24,1);
exactCol=zeros(24,1); smallCol=zeros(24,1); errorCol=zeros(24,1);
statusCol=cell(24,1); nrow=0;
for cidx=1:4
    mixed = 2*rf{cidx}.*exp(-1j*2*pi*p.f0*t);
    filtered = conv(mixed,h,'same'); % Odd FIR: offline delay already compensated.
    iq = filtered(1:p.decimation:end);
    assert(max(abs(iq-iqReference{cidx})) < 1e-10, 'Full RF-to-IQ mismatch');
    x=iq(g);
    q=x(1+L:end,:).*conj(x(1:end-L,:));
    R=mean(q,1).'; fd=angle(R)/(2*pi*p.lag_seconds);
    vSmall=-p.c*fd/(2*p.f0); vExact=-p.c*fd./(2*p.f0+fd);
    if mod(cidx,2)==0
        tol=p.noisy_tolerance_mps;
    else
        tol=p.clean_tolerance_mps;
    end
    modelError=abs(vExact-S.truth.velocity_exact_from_wrapped_expected_mps);
    assert(all(modelError<tol), 'Exact-inversion result outside tolerance');
    assert(all(abs(vSmall-S.truth.velocity_small_v_expected_mps)<tol), 'Small-v result outside tolerance');
    valid=~logical(S.truth.is_ambiguous);
    assert(all(abs(vExact(valid)-S.truth.velocity_mps(valid))<tol), 'Physical speed outside tolerance');
    for z=1:6
        nrow=nrow+1;
        datasetCol{nrow}=pathNames{cidx}; zoneCol(nrow)=z;
        truthCol(nrow)=S.truth.velocity_mps(z); exactCol(nrow)=vExact(z);
        smallCol(nrow)=vSmall(z); errorCol(nrow)=modelError(z);
        if S.truth.is_ambiguous(z)
            statusCol{nrow}='EXPECTED_ALIAS';
        else
            statusCol{nrow}='PASS';
        end
    end
end
results=table(datasetCol,zoneCol,truthCol,exactCol,smallCol,errorCol,statusCol,...
    'VariableNames',{'Path','Zone','Truth_mps','Exact_mps','SmallV_mps','ModelError_mps','Status'});
disp(results);
fprintf('PASS: %s\n',matPath);
fprintf('EXPECTED_ALIAS validates the wrapped result, not true-velocity recovery.\n');
end

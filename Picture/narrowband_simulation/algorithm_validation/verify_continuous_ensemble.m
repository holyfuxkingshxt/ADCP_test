function v = verify_continuous_ensemble(matPath)
% Columns are artificial independent scattering realizations, not timed pings.
if nargin < 1
    matPath=fullfile(fileparts(mfilename('fullpath')),'continuous_ensemble_validation.mat');
end
S=load(matPath); p=S.params; x=S.iq_gates_tgc;
L=round(p.lag_seconds*p.fs_rf/p.decimation);
R=mean(x(1+L:end,:).*conj(x(1:end-L,:)),1);
Rmean=mean(R); % Average COMPLEX correlations before taking the angle.
fd=angle(Rmean)/(2*pi*p.lag_seconds);
v=-p.c*fd/(2*p.f0+fd);
assert(abs(v-S.truth.velocity_mps)<S.truth.tolerance_mps,'Ensemble outside tolerance');
fprintf('PASS: %d independent realizations, truth %.8f, estimate %.8f m/s\n',...
    size(x,2),S.truth.velocity_mps,v);
end

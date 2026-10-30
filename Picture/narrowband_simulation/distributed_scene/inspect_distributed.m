% Inspect all ranges first, then the single selected receive gate.
% This script locates the MAT file relative to itself.
S = load(fullfile(fileparts(mfilename('fullpath')), 'distributed_narrowband_simulation.mat'));
g = S.gate_indices_matlab;

figure('Name','Full record and selected gate');
subplot(2,1,1);
plot(S.t_rf_s*1000,S.rx_rf_noisy); hold on;
xline(S.params.gate_start_s*1000,'r--');
xline((S.params.gate_start_s+S.params.gate_duration_s)*1000,'r--');
xlabel('Receive time from transmit start (ms)'); ylabel('Real RF amplitude');
title('Full 3-21 m distributed echo; dashed lines select one gate'); grid on;
subplot(2,1,2);
plot(S.t_gate_s*1000,abs(S.iq_gate_zones)); hold on;
plot(S.t_gate_s*1000,abs(S.iq_gate),'k','LineWidth',2);
xlabel('Receive time (ms)'); ylabel('Complex I/Q magnitude');
legend('3-6 m','6-9 m','9-12 m','12-15 m','15-18 m','18-21 m','Coherent total');
title('Contributions from each physical range zone'); grid on;

figure('Name','Far-range control');
plot(S.t_gate_s*1000,abs(S.iq_gate),'LineWidth',1.5); hold on;
plot(S.t_gate_s*1000,abs(S.iq_far_boosted(g)),'--','LineWidth',1.5);
legend('Original','18-21 m amplitude x10');
xlabel('Receive time (ms)'); ylabel('I/Q magnitude'); grid on;

L=round(S.params.lag_seconds*S.params.fs_rf/S.params.decimation);
q=S.iq_gate(1+L:end).*conj(S.iq_gate(1:end-L));
phi=angle(mean(q));
fd=phi/(2*pi*S.params.lag_seconds);
v=-S.params.c*fd/(2*S.params.f0);
fprintf('Single realization: phase=%.6f rad, Doppler=%.6f Hz, speed=%.6f m/s\n',phi,fd,v);
disp(S.diagnostics.estimates);

% Time gain is applied to the complete RF record, before I/Q and its FIR.
% Keep both paths visible; *_tgc names always identify compensated data.
figure('Name','Receive-time gain compensation');
subplot(3,1,1);
plot(S.t_rf_s*1000,S.tgc_gain_db_rf,'LineWidth',1.5);
xlabel('Receive time (ms)'); ylabel('Gain (dB)'); grid on;
title('TGC: inverse assumed pressure loss, referenced to 3 m');
subplot(3,1,2);
plot(S.t_rf_s*1000,S.rx_rf_noisy); hold on;
plot(S.t_rf_s*1000,S.rx_rf_noisy_tgc);
xlabel('Receive time (ms)'); ylabel('Real RF amplitude'); grid on;
legend('Raw','With TGC'); title('Signal and receiver noise both receive gain');
subplot(3,1,3);
plot(S.t_iq_s*1000,abs(S.iq_noisy)); hold on;
plot(S.t_iq_s*1000,abs(S.iq_noisy_tgc));
xlabel('Receive time (ms)'); ylabel('I/Q magnitude'); grid on;
legend('Raw','With TGC');

q_tgc=S.iq_gate_noisy_tgc(1+L:end).*conj(S.iq_gate_noisy_tgc(1:end-L));
v_tgc=-S.params.c*angle(mean(q_tgc))/(4*pi*S.params.f0*S.params.lag_seconds);
fprintf('After TGC (noisy gate): speed=%.6f m/s\n',v_tgc);
disp(S.diagnostics.time_gain_compensation);

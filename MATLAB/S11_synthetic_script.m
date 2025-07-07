% This creates a synthetic multi lorentzian peak S11 curve to test the
% capability of the fitting algorithm

% frequency range and number of points
f = linspace(3.93e9, 3.936e9, 4001);

% p: [f0_1, k_int_1, k_ext_1, ..., f0_N, k_int_N, k_ext_N]
p = [3.93235e9,90e3, 10e3,...
     3.93258e9,40e3, 40e3,...
     3.93263e9,90e3, 1e3,...
     3.93275e9,40e3, 1e3,...
     3.93285e9,40e3, 40e3,...
     3.93301e9,40e3, 10e3];

% p = [3.93235e9,10e3, 10e3,...
%      3.93236e9,10e3, 10e3,...
%      3.93237e9,10e3, 10e3,...
%      3.93238e9,10e3, 10e3,...
%      3.93239e9,10e3, 10e3,...
%      3.93236e9,10e3, 10e3,...
%      3.93237e9,10e3, 10e3,...
%      3.93238e9,10e3, 10e3];

num_peaks = length(p)/3;

S11 = aafunc_S11_complex_MPL(p,f);
% Magnitude in linear units
mag_lin = abs(S11);             % |S11|
% Magnitude in dB
mag_db = 20 * log10(mag_lin);   % 20*log10(|S11|)
% Phase in radians
phase_rad = angle(S11);         % phase(S11) in radian
% Phase in degrees
phase_deg = rad2deg(phase_rad);

figure;
subplot(2,2,1);
plot(f/1e9, real(S11), 'b');
xlabel('Frequency (GHz)'); ylabel('Re(S_{11})');

subplot(2,2,3);
plot(f/1e9, imag(S11), 'b');
xlabel('Frequency (GHz)'); ylabel('Im(S_{11})');

subplot(2,2,[2 4]);
plot(real(S11), imag(S11), 'b');
axis equal; grid on; xlabel('Re'); ylabel('Im'); title('S_{11} in Complex Plane');

date = datetime('now');
date.Format = 'yyyy_MM_dd';
datestr = string(date);

filename = sprintf('../Data/synthetic_S11_%d_peaks_%s',num_peaks,datestr);

saveas(gcf,append(filename,'.png'))
savefig(append(filename,'.fig'))
save(append(filename,'.mat'),'f','S11','mag_lin','mag_db','phase_deg','phase_rad');
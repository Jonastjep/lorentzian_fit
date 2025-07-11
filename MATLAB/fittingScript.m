clear

%% Loading data
filename = 'synthetic_S11_6_peaks_2025_07_07';

% full .mat datafile filepath
filepath = append('..\Data\', filename,'.mat');

% Synthetic data loading in case of need
load(filepath);

% Real NEMS data loading
% [f,mag_db,mag_lin,phase_deg,phase_rad] = aafunc_real_data_import(filepath);

%% Data processing
% The phase is not perfect from the measurement, so we need to unwrap it 
% (to remove discontinuities) and and to fit a line to it (there is a 
% phase shift due to cable propagation delay)
phase_rad = unwrap(phase_rad); % unwrapping

% The large numbers in the frequencies (GHz, ~1e9) cause the polynomial 
% fitting matrix to be ill-conditioned. The easiest and most stable 
% solution is to center the frequency vector. (We could also just scale the
% x axis to GHz scale e.g. f./1e9)
f0 = mean(f);  
pf = polyfit(f - f0, phase_rad, 1); % line fitting
phase_rad = phase_rad - polyval(pf, f - f0); % substract the fitted line

% recreate the complete complex signal from the magnitude and phase
S11_measured = mag_lin .* exp(1j * phase_rad);


%% Manual peak selection (f0 initial guesses)
% take the derivative of the phase and try to find possible peaks of interest
thresh_val = 0.5; %minima threshold values
smooth_val = 20;  %derivative smoothing window
[dphi_df,res_freqs,f_mid,final_locs] = aafunc_peaks_from_S11phase(f, phase_rad, thresh_val, smooth_val);

% manual peak selection process for initial guesses of peak positions
f0_init = aafunc_peak_selection(f,mag_db,phase_rad, dphi_df,res_freqs,f_mid,final_locs);

%% Initial guesses for ki and k1
% Default initial guesses for k1 and ki
num_peaks = length(f0_init);
k1_init = 5e4 * ones(1, num_peaks);
ki_init = 6e4 * ones(1, num_peaks);

% Build initial parameter vector p0: [x0_1, ki_1, k1_1, x0_2, ki_2, k1_2, ...]
p0 = zeros(1, 3 * num_peaks);
for i = 1:num_peaks
    p0(3*i - 2) = f0_init(i);
    p0(3*i - 1) = ki_init(i);
    p0(3*i) = k1_init(i);
end

%% Set upper and lower bounds for parameter search
%Upper and lower bounds on parameters
lb = zeros(1, 3*num_peaks);        % Lower bound: 0
ub = inf(1, 3*num_peaks);          % Upper bound: infinity
for i = 1:num_peaks
    %upper and lower bounds of frequency
    lb(3*i - 2) = f0_init(i)-1e6;
    ub(3*i - 2) = f0_init(i)+1e6;
end

%% Fitting real+imag components
opts = optimoptions('lsqnonlin', 'Display', 'iter', 'TolFun', 1e-12,'MaxIterations',1000);
[p_fit, resnorm] = lsqnonlin(@(p) aafunc_S11_residual(p, f, S11_measured),...
                               p0, lb, ub, opts);

%% Plotting Re and Im parts (also the combined complex plane plot)
% Reconstruct the fitted S11 complex function with resulting fitting
% parameters p_fit
S11_fit = aafunc_S11_complex_MPL(p_fit, f);

figure;
subplot(2,2,1);
plot(f/1e9, real(S11_measured), 'b', f/1e9, real(S11_fit), 'r--');
xlabel('Frequency (GHz)'); ylabel('Re(S_{11})'); legend('Measured','Fit');

subplot(2,2,3);
plot(f/1e9, imag(S11_measured), 'b', f/1e9, imag(S11_fit), 'r--');
xlabel('Frequency (GHz)'); ylabel('Im(S_{11})');

subplot(2,2,[2 4]);
plot(real(S11_measured), imag(S11_measured), 'b', real(S11_fit), imag(S11_fit), 'r--');
axis equal; grid on; xlabel('Re'); ylabel('Im'); title('S_{11} in Complex Plane');

%% Plotting data and fitted magnitude and phase
% Initialize for plotting individual peaks
S11_indPeaks = zeros(num_peaks, length(f));
% Loop to extract individual peaks
for i = 1:num_peaks
    p = p_fit(3*i - 2:3*i);

    % Evaluate single Lorentzian
    S11_indPeaks(i, :) = aafunc_S11_complex_SPL(p,f);
end

colors = lines(num_peaks);  % distinct colors for each peak
figure;
subplot(2,1,1);
hold on
plot(f/1e9, 20*log10(abs(S11_measured)), 'b', f/1e9, 20*log10(abs(S11_fit)), 'r--');
xlabel('Frequency (GHz)'); ylabel('|S_{11}| dB'); legend('Measured','Fit');
for i = 1:num_peaks
    plot(f/1e9, 20*log10(abs(S11_indPeaks(i, :))),'Color', colors(i,:),'DisplayName', sprintf('Peak %d fit', i));
end
legend;
hold off

subplot(2,1,2);
hold on
plot(f/1e9, angle(S11_measured), 'b', f/1e9, angle(S11_fit),'r--','DisplayName', sprintf('Peak %d fit', i));
xlabel('Frequency (GHz)'); ylabel('phase S_{11}');
for i = 1:num_peaks
    plot(f/1e9, unwrap(angle(S11_indPeaks(i, :))),'Color', colors(i,:),'DisplayName', sprintf('Peak %d fit', i));
end
legend;
hold off

figure;
subplot(2,1,1);
hold on
plot(f/1e9, 20*log10(abs(S11_measured)), 'b', f/1e9, 20*log10(abs(S11_fit)), 'r--');
xlabel('Frequency (GHz)'); ylabel('|S_{11}| dB'); legend('Measured','Fit');
legend;
hold off

subplot(2,1,2);
hold on
plot(f/1e9, angle(S11_measured), 'b', f/1e9, angle(S11_fit),'r--','DisplayName', sprintf('Peak %d fit', i));
xlabel('Frequency (GHz)'); ylabel('Phase S_{11}'); legend('Measured','Fit');
legend;
hold off

date = datetime('now');
date.Format = 'yyyy_MM_dd_hhmmss';
datestr = string(date);

filename = sprintf('../Fits/%s_S11_fit_%d_peaks_%s',filename,num_peaks,datestr);

savefig(append(filename,'.fig'))
save(append(filename,'.mat'),'S11','mag_lin','mag_db','phase_deg','phase_rad');

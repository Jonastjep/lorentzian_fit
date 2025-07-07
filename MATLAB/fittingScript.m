% Load your data
clear
load("Data\2025-06-05-103509_HBAR_coupling_J3_50mK.mat",'dfreq','ddata','ddatamag','ddatamaglin');
% dfreq -> the frequency range of the measurement (in Hz)
% ddata -> the phase of the S11 measurement (in deg)
% ddatamag -> the magnitude (dB) of the S11 measurement
% ddatamaglin -> the magnitude (linear units) of the S11 measurement
% ddatamag -> the magnitude (linear units) of the S11 measurement

%%% PROOF DATAMAG IS IN DB %%%
% figure;
% hold on;
% plot(f/1e9, 20*log10(ddatamaglin), 'r-');
% plot(f/1e9, ddatamag, 'b--');
% legend;

% dfreq is a row vector instead of a column vector like the other variables
f = transpose(dfreq);

% convert the magnitude data (|S11| dB) back to linear magnitude
mag_db = ddatamag;
mag_lin = ddatamaglin;

% the phase is messy from the measurement, so we need to unwrap it (to
% remove discontinuities), convert it to radiansand and to fit it to a line
%(there is a phase shift due to propagation delay)
phase_rad = deg2rad(ddata);

phase_rad = unwrap(phase_rad);

f0 = mean(f);  % center frequency
pf = polyfit(f - f0, phase_rad, 1);

phase_rad = phase_rad - polyval(pf, f - f0);

% recreate the complete complex signal from the magnitude and phase
S11_measured = ddatamaglin .* exp(1j * phase_rad);

[dphi_df,res_freqs,f_mid,final_locs] = peaks_from_S11phase(f, phase_rad, 0.5, 20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FAKE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This creates a synthetic multi lorentzian peak S11 curve to test the
% capability of the fitting algorithm
f = linspace(3.93e9, 3.936e9, 4001);
% p = [3.93235e9,90e3, 10e3,...
%      3.93258e9,40e3, 40e3,...
%      3.93263e9,90e3, 1e3,...
%      3.93275e9,40e3, 1e3,...
%      3.93285e9,40e3, 40e3,...
%      3.93301e9,40e3, 10e3];
p = [3.93235e9,10e3, 10e3,...
     3.93236e9,10e3, 10e3,...
     3.93237e9,10e3, 10e3,...
     3.93238e9,10e3, 10e3,...
     3.93239e9,10e3, 10e3,...
     3.93236e9,10e3, 10e3,...
     3.93237e9,10e3, 10e3,...
     3.93238e9,10e3, 10e3];

S11_fake = S11_complex_MPL(p,f);
% Magnitude in linear units
mag_lin = abs(S11_fake);             % |S11|
% Magnitude in dB
mag_db = 20 * log10(mag_lin);   % 20*log10(|S11|)
% Phase in radians
phase_rad = angle(S11_fake);         % phase(S11) in radian
% Phase in degrees
phase_deg = rad2deg(phase_rad);

[dphi_df,res_freqs,f_mid,final_locs] = peaks_from_S11phase(f, phase_deg, 0.5, 20);

S11_measured = S11_fake;

figure;
hold on;
plot(f/1e9, mag_db, 'r-');
plot(f/1e9, ddatamag, 'b--');
legend;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%% PEAK SELECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

selectFig = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5 0 0.5 1]);
ax1 = subplot(3,1,1);
p1 = plot(f, mag_db, 'b');
hold on;
s1 = scatter(res_freqs, interp1(f, mag_db, res_freqs), ...
    20, 'ro');
xlabel('Frequency (GHz)');
ylabel('|S_{11}|');
title('Magnitude with Estimated Resonances');
legend('boxoff');
legend('','Estimated resonances');
grid on;

ax2 = subplot(3,1,2);
plot(f, phase_rad, 'b');
hold on;
scatter(res_freqs, interp1(f, phase_rad, res_freqs), ...
    20, 'ro');
xlabel('Frequency (GHz)');
ylabel('Phase');
title('Phase with Estimated Resonances');
grid on;

ax3 = subplot(3,1,3);
plot(f_mid, dphi_df, 'k');
hold on;
scatter(res_freqs, dphi_df(final_locs), ...
    20, 'ro');
xlabel('Frequency (GHz)');
ylabel('dPhase/df');
title('Phase Derivative');
grid on;
linkaxes([ax1,ax2,ax3],'x');

num_peaks = input('Enter the number of peaks you want to fit to this spectrum:');
fprintf('Please select the %d peaks from the plot (a crosshair will appear)\n',num_peaks)

[x0,y0] = ginput(1);
Lx = plot(x0,y0,'rx','LineWidth',2,'MarkerSize',12);
fprintf('Peak 1/%d selected.\n',num_peaks);
for k = 1:num_peaks-1
    [x0,y0] = ginput(1);
    XData = get(Lx, 'XData');
    YData = get(Lx, 'YData');
    set(Lx, 'Xdata', [x0, XData], 'YData', [y0, YData]);
    fprintf('Peak %d/%d selected.\n',k+1,num_peaks);
end
XData = get(Lx, 'XData');
%YData = get(Lx, 'YData');
hold off
close(selectFig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial guesses for peak positions
x0_init = XData;

% Default initial guesses for k1 and ki
k1_init = 5e4 * ones(1, num_peaks);
ki_init = 6e4 * ones(1, num_peaks);

%Upper and lower bounds on parameters
lb = zeros(1, 3*num_peaks);        % Lower bound: 0
ub = inf(1, 3*num_peaks);          % Upper bound: infinity

% Build initial parameter vector: p contains [x0_1, ki_1, k1_1, x0_2, ki_2, k1_2, ...]
p0 = zeros(1, 3 * num_peaks);
for i = 1:num_peaks
    p0(3*i - 2) = x0_init(i);
    p0(3*i - 1) = ki_init(i);
    p0(3*i) = k1_init(i);

    %upper and lower bounds of frequency
    lb(3*i - 2) = x0_init(i)-1e6;
    ub(3*i - 2) = x0_init(i)+1e6;
end

% Use lsqnonlin to fit real+imag components
opts = optimoptions('lsqnonlin', 'Display', 'iter', 'TolFun', 1e-12,'MaxIterations',1000);
[p_fit, resnorm] = lsqnonlin(@(p) S11_residual(p, f, S11_measured),...
                               p0, lb, ub, opts);

%%%%%%%%%% PLOTTING THE DATA AND FITTED Re AND Im PARTS %%%%%%%%%%%%%%%%%%%
S11_fit = S11_complex_MPL(p_fit, f);

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

%%%%%%%%%%%% PLOTTING THE DATA AND FITTED MAG AND PHASE %%%%%%%%%%%%%%%%%%%
%%%%%% OF COMPLETE S11 CURVE AND INDIVIDUAL FITTED S11 CURVES %%%%%%%%%%%%%
% Initialize for plotting individual peaks
S11_indPeaks = zeros(num_peaks, length(f));
% Loop to extract individual peaks
for i = 1:num_peaks
    p = p_fit(3*i - 2:3*i);

    % Evaluate single Lorentzian
    S11_indPeaks(i, :) = S11_complex_SPL(p,f);
end

colors = lines(num_peaks);  % distinct colors for each peak
figure;
subplot(2,2,1);
hold on
plot(f/1e9, 20*log10(abs(S11_measured)), 'b', f/1e9, 20*log10(abs(S11_fit)), 'r--');
xlabel('Frequency (GHz)'); ylabel('|S_{11}| dB'); legend('Measured','Fit');
for i = 1:num_peaks
    plot(f/1e9, 20*log10(abs(S11_indPeaks(i, :))),'Color', colors(i,:),'DisplayName', sprintf('Peak %d fit', i));
end
legend;
hold off

subplot(2,2,3);
hold on
plot(f/1e9, imag(S11_measured), 'b', f/1e9, imag(S11_fit),'Color', colors(i,:),'DisplayName', sprintf('Peak %d fit', i));
xlabel('Frequency (GHz)'); ylabel('Im(S_{11})');
for i = 1:num_peaks
    plot(f/1e9, unwrap(angle(S11_indPeaks(i, :))),'Color', colors(i,:),'DisplayName', sprintf('Peak %d fit', i));
end
legend;
hold off

subplot(2,2,2);
hold on
plot(f/1e9, 20*log10(abs(S11_measured)), 'b', f/1e9, 20*log10(abs(S11_fit)), 'r--');
xlabel('Frequency (GHz)'); ylabel('|S_{11}| dB'); legend('Measured','Fit');
legend;
hold off

subplot(2,2,4);
hold on
plot(f/1e9, imag(S11_measured), 'b', f/1e9, imag(S11_fit),'Color', colors(i,:),'DisplayName', sprintf('Peak %d fit', i));
xlabel('Frequency (GHz)'); ylabel('Phase S_{11}');
legend;
hold off
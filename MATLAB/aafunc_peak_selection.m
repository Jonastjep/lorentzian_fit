function XData = aafunc_peak_selection(f,mag_db,phase_rad, dphi_df,res_freqs,f_mid,final_locs)
%AAFUNC_PEAK_SELECTION Summary of this function goes here
%   Detailed explanation goes here
input('Welcome to the lorentzian fitting script. Press ENTER to start, and follow the instructions in the terminal.');

selectFig = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5 0 0.5 1]);
ax1 = subplot(3,1,1);
plot(f, mag_db, 'b');
hold on;
scatter(res_freqs, interp1(f, mag_db, res_freqs), 20, 'ro');
xlabel('Frequency (GHz)');
ylabel('|S_{11}|');
title('Magnitude with Estimated Resonances');
legend('boxoff');
legend('','Estimated resonances');
grid on;

ax2 = subplot(3,1,2);
plot(f, phase_rad, 'b');
hold on;
scatter(res_freqs, interp1(f, phase_rad, res_freqs), 20, 'ro');
xlabel('Frequency (GHz)');
ylabel('Phase');
title('Phase with Estimated Resonances');
grid on;

ax3 = subplot(3,1,3);
plot(f_mid, dphi_df, 'k');
hold on;
scatter(res_freqs, dphi_df(final_locs), 20, 'ro');
xlabel('Frequency (GHz)');
ylabel('dPhase/df');
title('Phase Derivative');
grid on;

linkaxes([ax1,ax2,ax3],'x'); %make axes same scale and move together 

num_peaks = input(['The interactive plot shown is to help you identify peaks.' ...
    ' Once you enter the number of peaks you will be brought to an interactive' ...
    ' graphical peak selection mode. At that stage you won''t be able to pan and' ...
    ' zoom the plots anymore.\n\n' ...
    'Enter the number of peaks you want to fit to this spectrum:']);
fprintf('Please select the %d peaks from the plot (a crosshair will appear)\n',num_peaks)

[x0,y0] = ginput(1);
Lx = plot(x0,y0,'rx','LineWidth',2,'MarkerSize',8);
fprintf('Peak 1/%d selected.\n',num_peaks);
for k = 1:num_peaks-1
    [x0,y0] = ginput(1);
    XData = get(Lx, 'XData');
    YData = get(Lx, 'YData');
    set(Lx, 'Xdata', [x0, XData], 'YData', [y0, YData]);
    fprintf('Peak %d/%d selected.\n',k+1,num_peaks);
end
XData = get(Lx, 'XData');
% YData = get(Lx, 'YData');
hold off
close(selectFig)
end


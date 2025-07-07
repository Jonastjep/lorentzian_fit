function [dphi_df,res_freqs,freq_mid,final_locs] = aafunc_peaks_from_S11phase(freq, phase_unwrapped, thresh_val, smooth_val)
% FIND_RESONANCES_FROM_PHASE Detects resonance peaks from the unwrapped S11 phase.
%
% INPUTS:
%   freq            - frequency vector (Hz)
%   phase_unwrapped - unwrapped phase of S11 (radians or degrees)
%   thresh_val      - minima threshold factor (0.5)
%   smooth_val      - smoothing window factor (5)
%
% OUTPUTS:
%   dphi_df         - phase derivative with respect to frequency
%   res_freqs       - estimated resonance frequencies (Hz)
%   freq_mid        - midpoints between adjacent frequency values

    % Ensure column vectors
    freq = freq(:);
    phase_unwrapped = phase_unwrapped(:);

    % Compute derivative of phase w.r.t. frequency
    dphi_df = diff(phase_unwrapped) ./ diff(freq);
    freq_mid = (freq(1:end-1) + freq(2:end)) / 2;

    % Smooth derivative to suppress noise (optional)
    dphi_df_smooth = smooth(dphi_df, smooth_val);  % adjust window size as needed

    % Find local maxima and minima in derivative magnitude
    % [~, locs_max] = findpeaks(dphi_df_smooth);
    [~, locs_min] = findpeaks(-dphi_df_smooth);
    all_locs = sort(locs_min);
    % all_locs = sort([locs_max; locs_min]);

    % Apply threshold to avoid false positives (adjust as needed)
    threshold = thresh_val * max(abs(dphi_df_smooth));
    significant = abs(dphi_df_smooth(all_locs)) > threshold;
    final_locs = all_locs(significant);

    % Estimate resonance frequencies
    res_freqs = freq_mid(final_locs);
end
function [f,mag_db,mag_lin,phase_deg,phase_rad] = aafunc_real_data_import(filepath)
%AAFUNC_REAL_DATA_IMPORT Summary of this function goes here
%   Detailed explanation goes here
load(filepath,'dfreq','ddata','ddatamag','ddatamaglin');
% dfreq -> the frequency range of the measurement (in Hz)
% ddata -> the phase of the S11 measurement (in deg)
% ddatamag -> the magnitude (dB) of the S11 measurement
% ddatamaglin -> the magnitude (linear units) of the S11 measurement

% dfreq is a row vector instead of a column vector like the other variables
if size(dfreq) ~= size(ddatamaglin)
    f = transpose(dfreq);
else
    f = dfreq;
end

% PHASE
phase_deg = ddata;
phase_rad = deg2rad(ddata);

% MAGNITUDE
mag_db = ddatamag;
mag_lin = ddatamaglin;

end


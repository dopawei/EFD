function [inst_freq, inst_amp] = IFIA(x,fs)
% =========================================================================
% This function is used to obtain instantaneous frequencies and amplitudes
% of signal component by Hilber transform
% Inputs:
%   -x: signal component
%   -fs:sampling frequency
%
% Outputs:
%   -inst_freq: instantaneous frequencies
%   -inst_amp: instantaneous amplitudes
%
% Author: Wei Zhou
% Institution: Department of Mechanical and Materials Engineering,
% University of Cincinnati, Cincinnati, OH 45221, USA
% Year: 2022
% Version: 2.0
% Reference: Empirical Fourier decomposition: An accurate signal decomposition method
% for nonlinear and non-stationary time series analysis
% https://doi.org/10.1016/j.ymssp.2021.108155
% =========================================================================
[m,n] = size(x);
% convert row data to coloumn data
if (m<n)
    x = x';
end
[m,n] = size(x);
% initialize variables
inst_freq0 = zeros(m,n);
inst_amp0 = inst_freq0;
% compute instantaneous frequencies and amplitudes
for i = 1:n
    [inst_amp, inst_freq] = comp_inst_fre_amp(x(:,i),fs);
    inst_amp0(:,i) = inst_amp;
    inst_freq0(:,i) = inst_freq;
end

function [inst_amp, inst_freq] = comp_inst_fre_amp(x,fs)
z = hilbert(x);
inst_amp = abs(z);
phi = unwrap(angle(z));
tmp_phi = phi;
diff_phi = ((tmp_phi(3:end)-tmp_phi(1:end-2))/2);
diff_phi = [diff_phi(1);diff_phi;diff_phi(end)];
inst_freq = diff_phi*(fs/(2*pi));
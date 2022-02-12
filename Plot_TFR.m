function [nt,tscale,fscale] = Plot_TFR(inst_freq,inst_amp,T)
% =========================================================================
% This function is used to plot TFR
% Inputs:
%   -inst_freq: instantaneous frequencies
%   -inst_amp: instantaneous amplitudes
%   -T: total time
%
% Outputs:
%   -nt: amplitude for TFR
%   -tscale: gird for time
%   -fscale: gird for frequency
%
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
t0 = 0;
t1 = T;
multfactor = 4;
if(length(inst_freq(:,1)) >= 100*multfactor)
    fres = 100*multfactor;
    tres = 100*multfactor;
else
    fres = length(inst_freq(:,1));
    tres = fres;
end

fw0 = min(min(inst_freq));
fw1 = max(max(inst_freq));

if (fw0 < 0)
    fw0 = 0;
end
tw0 = t0;
tw1 = t1;

%  nspplote.m to plot the time-frequency spectrum
lscale = 0;
[nt,tscale,fscale] = nspplote(inst_freq,inst_amp,t0,t1,fres,tres,fw0,fw1,tw0,tw1,lscale);
end
function plotbounds(f,boundaries)
% =========================================================================
% This function is used to plot the detected boundaries under the Fourier spectrum 
% Inputs:
%   -f: the Fourier spectrum of input signal
%   -boundaries: the set of boundaries corresponding are detected by  Boundaries_Dect
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
figure
magf=abs(fft(f));
freq=2*pi*[0:length(magf)-1]/length(magf);
R=round(length(magf)/2);
plot(freq(1:R),magf(1:R));
xlim([0,pi])
ylim([0,max(magf)])
for i=1:length(boundaries)
    hold on
    line([boundaries(i) boundaries(i)],[0 max(magf)],'Marker','o','MarkerSize',3,'LineWidth',1,'LineStyle','--','Color',[1 0 1]);
end
end


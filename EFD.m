function [efd,cerf,bounds] = EFD(x,N)

% =========================================================================
% This function is used to construct a zero-phase filter bank
% function [efd,cerf,bounds] = EFD(x,N)

% Inputs:
%   -x: the input signal, a vector 
%   -N: maximum number of segments
%
% Outputs:
%   -efd: cell containing first the successives frequency subbands
%   -cerf: vector containing the central frequency of each band, [0,Pi]
%   -bounds: vector containing the set of boundaries corresponding
%                to the Fourier line segmentation (normalized between
%                0 and Pi)
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

%% Boundary detection

% We compute the Fourier transform of f
[a,b] = size(x);
if b > a
    x = x';
end
ff = fft(x);

% We extract the boundaries of Fourier segments
[bounds,cerf] = Segm_tec(abs(ff(1:round(length(ff)/2))),N);

% We trun the boundaries to [0,pi]
bounds = bounds*pi/round(length(ff)/2);
%% Filtering

% We extend the signal by miroring to deal with the boundaries
l = round(length(x)/2);
x = [x(l-1:-1:1);x;x(end:-1:end-l+1)];
ff = fft(x);

% We obtain the boundaries in the extend f
bound2 = ceil(bounds*round(length(ff)/2)/pi);

% We get the core of filtering
efd = cell((length(bound2)-1),1);
ft = zeros(length(efd),length(ff));

% We define an ideal functions and extract components
for k = 1:length(efd)
    if bound2(k)==0
       ft(k,1:bound2(k+1)) = ff(1:bound2(k+1));
       ft(k,length(ff)+2-bound2(k+1):(length(ff))) = ff(length(ff)+2-bound2(k+1):length(ff));
    else
        ft(k,bound2(k):bound2(k+1)) = ff(bound2(k):bound2(k+1));
        ft(k,length(ff)+2-bound2(k+1):length(ff)+2-bound2(k)) = ff(length(ff)+2-bound2(k+1):length(ff)+2-bound2(k));
    end
    efd{k} = real(ifft(ft(k,:)));
    efd{k} = efd{k}(l:end-l);
end

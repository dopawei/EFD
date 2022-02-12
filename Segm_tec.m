function [bounds,cerf] = Segm_tec(f,N)
% =========================================================================
% This function is used to implement the improved segmentation technique 
% function [boundaries,cerf] = Boundaries_Dect(ff,N)
% Inputs:
%   -f: the Fourier spectrum of input signal
%   -N: maximum number of segments
%
% Outputs:
%
%   -cerf: vector containing the central frequency of each band,[0,pi]
%   -bound: vector containing the set of boundaries corresponding
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

% 1. detect the local maxima and minina
locmax = zeros(size(f));
locmin = max(f)*ones(size(f));
for i=2:length(f)-1
    if ((f(i-1)<f(i)) && (f(i)>f(i+1)))
        locmax(i)=f(i);
    end
    if ((f(i-1)>f(i)) && (f(i)<f(i+1)))
        locmin(i)=f(i);
    end
end
locmax(1) = f(1);
locmax(end) = f(end);

% 2. keep the N-th highest maxima and their index
if N~=0
    [lmax,Imax] = sort(locmax,1,'descend');
    if length(lmax)>N
        Imax = sort(Imax(1:N));
    else
        Imax = sort(Imax);
        N = length(lmax);
    end
    % 3. detect the lowest minima between two consecutive maxima
    M = N+1;% numbers of the boundaries
    omega = [1;Imax;length(f)]; %location
    bounds = zeros(1,M);
    for i=1:M
        if (i == 1 || i == M) && (omega(i) == omega(i+1))
            bounds(i) = omega(i)-1;
        else
            [lmin,ind] = min(f(omega(i):omega(i+1)));
             bounds(i) = omega(i)+ind-2;
        end
        cerf = Imax*pi/round(length(f));
    end
end


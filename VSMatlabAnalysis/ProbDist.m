%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8/1/07 ProbDist.mat: Creates a normalized histogram,
%  then overlays a plot of the Probability distribution
%  based from the scattering theory calculation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs:
%   histdata-histogram data from the ScatteringHistogram.cpp
%       code.
%   x0-min value of histogram domain
%   xf-max value of histogram domain
%   binsize-size of histogram bin
%   probdata - probability distribution data from the 
%      ScatteringHistogram.cpp code.


function [h,hn,xspan] = ProbDist(histdata, x0, binsize, xf, probdata)

xspan=[x0:binsize:xf];
h=hist(histdata,xspan); % Generate histogram
hn = h/(length(histdata)*binsize); % Normalize histogram to 1
bar(xspan, hn); % Plot histogram

hold on;

plot( probdata(:,1), probdata(:,2), 'k+' ) % Plot overlayed Prob Density
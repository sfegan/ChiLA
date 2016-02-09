%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8/1/07 histn.m: Creates a normalized histogram.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h,hn,xspan] = histn(data, x0, binsize, xf)

xspan=[x0:binsize:xf];
h=hist(data,xspan);            % Generate histogram
hn = h/(length(data)*binsize); % Normalize histogram to 1
bar(xspan, hn);                % Plot histogram

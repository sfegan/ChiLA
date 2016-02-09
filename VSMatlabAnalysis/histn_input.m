%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8/1/07 histn_input.mat: Creates a normalized histogram
%    using two input files, 1) histogram data 2) Probability
%    density plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hn,xspan] = histn(x0, binsize, xf)

xspan=[x0:binsize:xf];
data =load('HistogramData_z0.1_E0.05.txt');
ProbDist=load('EBLPhotonProbDens_z0.1_E0.05.txt');
h=hist(data,xspan); % Generate histogram
hn = h/(length(data)*binsize); % Normalize histogram to 1
bar(xspan, hn); % Plot histogram
hold on;
plot(ProbDist(:,1),ProbDist(:,2), 'xr');
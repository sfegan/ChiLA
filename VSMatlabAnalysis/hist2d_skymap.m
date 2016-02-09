function [] = hist2d_skymap(hist_data)
% hist2d_skymap displays the 2-dimensional histogram
%
% EXAMPLE
%   
%   
%   

if ~isstruct(hist_data) || ...
    ~(isfield(hist_data,'contents') && isfield(hist_data,'variance'))
    error('hist2d_skymap:incorrectInput', ...
        'hist2d_skymap expects a struct with fields contents and variance');
end


x = hist_data.x;
y = hist_data.y;
z = hist_data.z;

imagesc(x,y,z)
colorbar
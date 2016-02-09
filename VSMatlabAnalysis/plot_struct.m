function [] = plot_struct(data,logpl)
% PLOT_STRUCT(data)  Displays stage3 histograms.
%
%   INPUTS:
%
%     data: histogram struct from an .h5 file, using loadh5.
%
%     logpl: On a 1D error bar plot, if logpl == 0, then plot will
%     be standard 'plot'. Otherwise, it is semilogy.
%
%   EXAMPLE:
%
%     plot_struct(mrk421.results.significance_hist)
%
%     - OR - 
%
%     plot_struct(mrk421.results.significance_hist,0)
%
%   AUTHOR:
%
%       Timothy C. Arlen - timothyarlen@gmail.com - 2008-09-15
%

if ~exist('logpl')
   logpl = 1;
end

if isstruct(data) && ...
    ~(isfield(data,'contents') || isfield(data,'variance'))
    %----------------------------------
    %          2D Histogram     
    %----------------------------------
    x = data.x;
    y = data.y;
    z = data.z;

    imagesc(x,y,z)
    colorbar
    
elseif isstruct(data) && ...
    isfield(data,'contents') && isfield(data,'variance')
    %----------------------------------
    %      1DPlot with error bars
    %----------------------------------
    x = data.contents.x;
    y = data.contents.y;
    x0 = data.contents.lo_limit;
    xf = data.contents.hi_limit;
    E = sqrt(data.variance.y);
    
    if logpl == 0
        plot(x,y,'.r')
    else
        semilogy(x,y,'.r')
    end
    xlim([x0 xf]);
    hold on;
    errorbar(x,y,E,'.')
    hold off;

else
    error('incorrect input format for plot_struct')
end

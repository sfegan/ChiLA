%STAGE3_ANALYSIS   Plot VERITAS stage3 analysis results
%
%
%   INPUTS:
%
%       file_or_data: filename to load and plot, or variable containing
%                     data previously loaded
%
%   EXAMPLE:
%
%       stage3_analysis('crab_s3.h5');
%
%       - OR -
%
%       data=loadh5('crab_s3.h5');
%       stage3_analysis(data)
%
%   AUTHOR:
%
%       Matthew Wood - mdwood@astro.ucla.edu - 2009-05-11
%
%   VERSION:
%
%       $Id: stage3_analysis.m,v 1.13 2010/03/17 22:16:52 matthew Exp $
%
function stage3_analysis(file_or_data,varargin)

    global analyze_revision_str analyze_revision
    analyze_revision_str = '1.00';
    analyze_revision = 1000;
    
    global sheet printed_sheet
    sheet = 0;
    printed_sheet = 0;

    global mx my mw mh
    mx=0.015;
    my=0.030;
    mw=0.97;
    mh=0.95;   
    
    global run_info
    run_info = struct;
    
    if(ischar(file_or_data))
        fprintf(1,'Loading %s...\n',file_or_data);
        data = loadh5(file_or_data);
    else
        data = file_or_data;
    end
    
    if(isfield(file_or_data,'results'))
        data = data.stage2;
    end

    sheets = [];

    src_name = data.stage3.src_name;    
    src_name(src_name==' ') = '_';
    analysis_name = '';
    if(nargin>3)
      src_name = varargin{3};
    end
      
    if(nargin>4)
      analysis_name = varargin{4};
    end
      
    
    printfile = '';
    if(nargin>1)
      printfile = varargin{1};
    else
      printfile=sprintf('%s_analysis.ps',src_name);
    end
    
    if(isfield(data.analyze,'revision'))
        analyze_revision = data.analyze.revision;
        analyze_revision_str = data.analyze.revision_string;
    end
    
    close all
    
    revision = '$Revision: 1.13 $';
    version = revision(12:length(revision)-2);
    sheettitle=sprintf(strcat('VERITAS analyze v%s/results v%s',...
                              ' - %s@%s'),...
                       analyze_revision_str,version,...
                       data.analyze.user,data.analyze.host);
    
    
    % ---------------------------------------------------------------------
    % Sheet 1
    % ---------------------------------------------------------------------

    if(newfigure(sheets))
        zsubplot(4,4,1);  draw_summary(data);
	zsubplot(4,4,2);  draw_obs(data);
	zsubplot(4,4,3);
	draw_hist2d(data,data.stage3.sky_counts_hist,...
		    'Uncorrelated Counts Map',4,4);
	
%	
	zsubplot(4,4,4);  draw_acceptance_summary(data);
	zsubplot(4,4,5);  draw_fov_acceptance(data);
	zsubplot(4,4,6);  draw_acceptance(data);
	
	zsubplot(4,4,7);  draw_projection(data,'x');
	zsubplot(4,4,8);  draw_projection(data,'y');

	zsubplot(4,4,9);  
	draw_skymap(data,data.results.sky_bkgnd_density_rate_hist,...
		    'Bkgnd Rate [min^{-1} deg^{-2}]');
	
	zsubplot(4,4,10);  
	draw_skymap(data,data.results.sky_significance_hist,...
		    'Significance [\sigma]');
	zsubplot(4,4,11);  draw_significance(data);
	
	
	zsubplot(4,4,12);  draw_thetasq(data);

	zsubplot(4,4,13);  
	draw_skymap(data,data.results.sky_excess_rate_hist,...
		    'Excess Rate [min^{-1}]');
	zsubplot(4,4,14);  
	draw_skymap(data,data.results.sky_flux_ul95_hist,...
		    '95% C.L. Flux Upper Limit [m^{-2}s^{-1}TeV^{-1}]','log');
	
	zsubplot(4,4,15);  draw_cumulative_significance(data);
	zsubplot(4,4,16);  draw_cumulative_excess(data);
	
	endfigure(sheettitle,printfile,4,4);
    end
    
    % ---------------------------------------------------------------------
    % Sheet 2
    % ---------------------------------------------------------------------
    if(newfigure(sheets))
      
      zsubplot(4,4,1);  draw_spectrum_summary(data);	
      zsubplot(4,4,2);  draw_spectrum_fit(data);
      zsubplot(4,4,3);  draw_reconstructed_energy(data);            
      zsubplot(4,4,4);  draw_spectrum(data);
      zsubplot(4,4,5);  draw_scope_participation(data);
      
      endfigure(sheettitle,printfile,4,4);
    end
    
    % ---------------------------------------------------------------------
    % Write Run Info
    % ---------------------------------------------------------------------
    
    if(nargin>2)
      write_source_info(data, varargin{2},src_name,analysis_name);
    end
    
    fprintf(1,'\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dosheet = newfigure(sheets)
    global sheet printed_sheet
    sheet = sheet+1;
    dosheet=1;
    if((~isempty(sheets))&&(isempty(find(sheets==sheet,1))))
        dosheet=0;
        return;
    end
    fprintf(1,'[%d]',sheet);
    printed_sheet = printed_sheet+1;
    
    figure
    orient landscape;
    set(gcf,'PaperPosition',[0 0 11 8.5]);
    set(gcf,'DefaultAxesFontSize',5)
    set(gcf,'DefaultTextFontSize',5)        
    set(gcf,'DefaultAxesLineWidth',get(gcf,'DefaultAxesLineWidth')*0.5);
    set(gcf,'DefaultLineLineWidth',get(gcf,'DefaultLineLineWidth')*0.5);
    colormap(my_colormap);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function endfigure(pagetitle,printfile,nx,ny)
    global sheet printed_sheet
    global mx my mw mh

    axes('Position',[0 0 1 1],'Visible','off');
    dx=0.01;
    dy=dx*11.0/8.5;
    text(0.0+dx,1.0-dx,pagetitle,...
         'VerticalAlignment','top','HorizontalAlignment','left',...
         'Interpreter','none');
    text(1.0-dx,1.0-dy,sprintf('Page %d',sheet),...
         'VerticalAlignment','top','HorizontalAlignment','right');
    
    dx=0.015;
    dy=-0.015;
    for ilab=1:nx*ny
        iy = floor((ilab-1)/nx);
        ix = ilab-iy*nx-1;
        x = mx+ix*mw/nx+dx;
        y = my+(ny-iy)*mh/ny+dy;
        text(x,y,sprintf('(%d%s)',sheet,char(96+ilab-32)),'Color','r');
    end
     
    if(printed_sheet==1)
        print('-dpsc2',printfile);
    else
        print('-dpsc2','-append',printfile);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZSUBPLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zsubplot(m,n,i)
    global mx my mw mh

%    subplot(m,n,i);
%    p=get(gca,'Position')
%    factor=0.02;
%    dx=p(3)*factor;
%    dy=p(4)*factor;
%    p=[p(1)-dx p(2)-dy p(3)+2*dx p(4)+2*dy];
%    set(gca,'Position',p);
    
    plotsize=0.8;
    ix=mod((i-1),m);
    iy=floor((i-1)/m);
    mx=0.015;
    my=0.030;
    mw=0.97;
    mh=0.95;   
    px=mx+mw*(ix+(1-plotsize)/1.4)/m;
    py=my+mh*(n-iy-1+(1-plotsize)/2)/n;
    pw=mw*plotsize/m;
    ph=mh*plotsize/n;
    p=[px py pw ph];
    subplot('Position',p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZAXIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zaxis(a)
    if(sum(isnan(a))==0)
        if(a(1)>=a(2))
            a(2)=a(1)*(1+eps)+eps;
        end
        if(a(3)>=a(4))
            a(4)=a(3)*(1+eps)+eps;
        end
        axis(a);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZSTAIRS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xs, ys] = zstairs(varargin)
    if(nargin==1)
        y=varargin{1};
        x=1:numel(y);
        if(size(y,1)>1)
            x=x';
        end
    else
        x=varargin{1};
        y=varargin{2};
    end
    dx=eps;
    if(nargin>=3)
        dx=varargin{3};
    end

    dim=1;

    if(numel(x) > 1)
        dx=median(diff(x));
        [xs,ys] = stairs(x,y);
        xs=cat(dim,xs,xs(numel(xs))+dx);
        ys=cat(dim,ys,ys(numel(ys)));
    elseif(numel(x) == 1)
        if(dx > 0)
            xs = [x x+dx]';
            ys = [y y]';
        else
            xs = [ x x*(1+2*eps) ]';
            ys = [ y y ]';
        end
    else
        xs = 0;
        ys = 0;
    end
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT_SIMPLEHIST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = plot_simplehist(hist, varargin)
    limits=[-inf inf];
    xscale = 1;
    yscale = 1;
    dx = 0;
    plotspec = '';
    if(nargin>=2)
        limits = varargin{1};
    end
    if(nargin>=3)
        plotspec = varargin{2};
    end
    if(nargin>=4)
        yscale = varargin{3};
    end
    if(nargin>=5)
        xscale = varargin{4};
    end
    if(nargin>=6)
        dx = varargin{5};
    end

    xh = double(hist.x);
    yh = double(hist.y);
          
    m = xh>=limits(1) & xh<=limits(2);

    if(~isempty(xh))
        [xs,ys] = zstairs(xh(m)*xscale,yh(m).*yscale,dx);
    else
        xs = 0;
        ys = 0;
    end

    if(nargout==2)
        varargout{1} = xs;    
        varargout{2} = ys;
    else
        plot(xs,ys,plotspec);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT_SIMPLEHIST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h = transpose(z)

  [nx ny] = size(z);
  h = zeros(ny,nx);
  
  for i=1:nx
    for j=1:ny	
      h(j,i) = z(nx-i+1,ny-j+1);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROUNDAXIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lim = roundaxis(val)
    max_val=max(double(val));
    if(max_val == 0)
        lim=1;
    else
        rounder = 10^ceil(log10(max_val/7))/2;
        lim = ceil(max_val/rounder+0.2)*rounder;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEGEND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = zlegend(varargin)
    if((nargin>0)&&(~isempty(varargin{1})))
        h=legend(varargin{:});
        legend boxoff
        if(nargout)
            varargout{1}=h;
        end
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MY COLORMAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cm = my_colormap
%    rcm = [ 0.0000 , 0.0000 , 0.2000 
%            0.0000 , 0.0000 , 1.0000
%            1.0000 , 0.0000 , 0.0000
%            1.0000 , 1.0000 , 0.3000 ];
%    x = (0:size(rcm,1)-1)/(size(rcm,1)-1);
%    cm = interp1(x,rcm,0:1/64:1);
%    cm(cm>1)=1;
%    cm(cm<0)=0;
    cm=jet;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND LIMITS CONTAINING SOME FRACTION OF ALL VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lo hi] = find_limits(data, fraction, extension)
    ds=sort(data(~isnan(data)&~isinf(data)));
    if(numel(ds) > 0)
        dsm=median(ds);
        dspos=ds(ds>dsm);
        dsneg=ds(ds<dsm);
        ihi=floor(numel(dspos)*fraction);
        ilo=floor(numel(dsneg)*(1-fraction));
        if(ihi > 0)
            hi=dspos(ihi);
            hi=extension*(hi-dsm)+dsm;
        else
            hi=dsm;
        end
        if(ilo > 0)
            lo=dsneg(ilo);
            lo=dsm-extension*(dsm-lo);
        else
            lo=dsm;
        end
        if(hi==lo)
            hi=hi+0.5;
            lo=lo-0.5;
        end
    else
        hi=1;
        lo=0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TWO DIMENSIONAL HISTOGRAM - WHY DOES MATLAB NOT HAVE THIS??
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h2 = simple_hist2(data, xlim, ylim, nxbin, nybin)
    h2 = zeros(nybin,nxbin);

    if(size(data,1)*size(data,2)==0)
        return
    end
    
    minx = xlim(1);
    maxx = xlim(2);
    miny = ylim(1);
    maxy = ylim(2);

    if(size(data,2)==2)
        data = data';
    end

    m = ~isnan(data(1,:)) & ~isnan(data(1,:)) ...
        & data(1,:)>minx & data(1,:)<maxx ...
        & data(2,:)>miny & data(2,:)<maxy;
    
    if(sum(m)==0)
        return
    end
    
    xc = data(1,m);
    yc = data(2,m);
    
    index = floor((xc-minx)/(maxx-minx)*nxbin)*nybin...
            + floor((yc-miny)/(maxy-miny)*nybin);

    h1 = hist(index, 0:(nxbin*nybin-1));
    h2(:)= h1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m = zmean(data)
    if(sum(~isnan(data)) > 0)
        m = mean(data(~isnan(data)));
    else
        m = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STANDARD DEVIATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = zstd(data)
    if(sum(~isnan(data)) > 0)
        d = std(data(~isnan(data)));
    else
        d = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = zvar(data)
    if(sum(~isnan(data)) > 0)
        v = var(data(~isnan(data)));
    else
        v = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEDIAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m = zmedian(data)
    if(sum(~isnan(data)) > 0)
        m = median(data(~isnan(data)));
    else
        m = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ELAPSED TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sec = elapsed_sec(data)
    D=data.diagnostics;
    if(isfield(D,'gps_ticks_elapsed'))
        sec = double(D.gps_ticks_elapsed)/1e7;
    else
        sec = double(D.l3_ticks_elapsed)/1e7;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_summary(data)
    D=data.results;

    disp('draw_summary');
    
    set(gca,'Box','on');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    if(0)
        p=get(gca,'Position');
        factor=0.1;
        dx=p(3)*factor;
        dy=p(4)*factor;
        p=[p(1)-2*dx p(2)-dy p(3)+2*dx p(4)+2*dy];
        set(gca,'Position',p)
    end

    nl=0;
    x0=+0.05;
    y0=+0.92;
    dy=-0.06;
    dx=+0.06;

    
    text(x0,y0+nl*dy,sprintf('%-20s %-30s','Source:',data.stage3.src_name));
    nl=nl+1;
    text(x0,y0+nl*dy,sprintf('%-20s %-30s','Analysis:',data.results.method));
    nl=nl+1;    
    text(x0,y0+nl*dy,sprintf('RA/DEC [hms,dms]: %15s %15s',...
			     data.stage3.src_ra_hms,data.stage3.src_dec_dms));

    nl=nl+1;    
    text(x0,y0+nl*dy,sprintf('Exposure [min]:'));
    text(x0+0.5,y0+nl*dy,sprintf('%-6.2f',D.elaptime_min));    
    nl=nl+1;
    text(x0,y0+nl*dy,sprintf('Livetime [min]:'));
    text(x0+0.5,y0+nl*dy,sprintf('%-6.2f',D.livetime_min));    
    nl=nl+1;
    
%    text(x0,y0+nl*dy,sprintf('%-20s %5.2f','Significance:',...
%			     D.significance));
    text(x0,y0+nl*dy,sprintf('Significance:'));
    text(x0+0.5,y0+nl*dy,sprintf('%-6.2f',D.significance));
    nl=nl+1;
%    text(x0,y0+nl*dy,sprintf('%-20s %5.1f \\pm %3.1f','Excess:',...
%			     D.excess,D.excess_err));
        
    text(x0,y0+nl*dy,sprintf('Background Rate [min^{-1} deg^{-2}]:'));
    text(x0+0.5,y0+nl*dy,sprintf('%-6.2f \\pm %4.2f',...
				 D.bkgnd_density_rate,...
				 D.bkgnd_density_rate_err));

    nl=nl+1;
    text(x0,y0+nl*dy,sprintf('Excess:'));
    text(x0+0.5,y0+nl*dy,sprintf('%-6.1f \\pm %3.1f',D.excess,D.excess_err));
    
    nl=nl+1;
%    text(x0,y0+nl*dy,sprintf('Excess Rate [min^{-1}]  %5.2f \\pm %3.2f',...
%			     D.excess_rate,D.excess_rate_err));

    text(x0,y0+nl*dy,sprintf('Excess Rate [min^{-1}]:'));
    text(x0+0.5,y0+nl*dy,sprintf('%-7.3f \\pm %4.3f',...
				 D.excess_rate,D.excess_rate_err));
    
    nl=nl+1;
    text(x0,y0+nl*dy,sprintf('Excess Rate 95%% C.L. UL [min^{-1}]:'));
    text(x0+0.5,y0+nl*dy,sprintf('%-6.2f',D.excess_rate_ul95));
    
    nl=nl+1;
    text(x0,y0+nl*dy,sprintf('Flux [m^{-2} s^{-1} TeV^{-1}]:'));
    text(x0+0.5,y0+nl*dy,sprintf('%6.3d \\pm %6.3d',D.flux,D.flux_err));
    
    nl=nl+1;
    text(x0,y0+nl*dy,sprintf('Flux 95%% C.L. UL [m^{-2} s^{-1} TeV^{-1}]:'));
    text(x0+0.5,y0+nl*dy,sprintf('%6.3d',D.flux_ul95));
    
    title('Summary');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACCEPTANCE SUMMARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_acceptance_summary(data)
  D=data.results;

  disp('draw_acceptance_summary');
  
  set(gca,'Box','on');
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);
  if(0)
    p=get(gca,'Position');
    factor=0.1;
    dx=p(3)*factor;
    dy=p(4)*factor;
    p=[p(1)-2*dx p(2)-dy p(3)+2*dx p(4)+2*dy];
    set(gca,'Position',p)
  end
  
  nl1=0;
  nl2=0;
  x0=+0.05;
  y0=+0.94;
  dy=-0.06;
  dx=+0.06;
  
  
  
  txt1 = {};
  txt2 = {};
  
  text(x0,y0+nl1*dy,...
       sprintf('Acceptance Model:  %s',data.stage3.acceptance.model));
  
  nl1 = nl1 + 2;
%  
%  text(x0,y0+nl1*dy,'   m = 0','FontWeight','bold');
%  nl1 = nl1 + 1;  
  
%  nl2 = nl2 + 2;
  
%  text(x0+0.5,y0+nl2*dy,'   m = 1','FontWeight','bold');
%  nl2 = nl2 + 1;  
  
%  txt1{numel(txt1)+1}=...
%      sprintf('Acceptance Model:  %s',data.stage3.acceptance.model);
%  txt2{numel(txt2)+1}='';
  
%  txt1{numel(txt1)+1}='';
%  txt2{numel(txt2)+1}='';
  
%  txt1{numel(txt1)+1}='m=0';
%  txt2{numel(txt2)+1}='m=1';
  

  if(data.stage3.acceptance.model == 'bessel2')

    
    
    for i = 1:numel(data.stage3.acceptance.m0_param.n)

      n = data.stage3.acceptance.m0_param.n(i);
      
      text(x0,y0+nl1*dy,...
	   sprintf('\\delta c_{0%i}:  %5.3f \\pm %5.3f',...
		   n,data.stage3.acceptance.m0_param.c(i),...
		   data.stage3.acceptance.m0_param.c_err(i)));
      nl1 = nl1 + 1;

    end
    
    
    nl1 = 8;
    
    for i = 1:numel(data.stage3.acceptance.mn_param{1}.n)

      p = data.stage3.acceptance.mn_param{1};      
      
      text(x0,y0+nl1*dy,...
	   sprintf('|\\delta c_{1%i}|:  %5.3f \\pm %5.3f',...
		   p.n(i),p.cr(i),p.cr_err(i)));
      nl1 = nl1 + 1;
      text(x0,y0+nl1*dy,...
	   sprintf('\\phi_{1%i}:  %5.1f \\pm %5.1f',...
		   p.n(i),p.cphi(i),p.cphi_err(i)));
      nl1 = nl1 + 1;
    end
    
    
    nl1 = 8;
    
    for i = 1:numel(data.stage3.acceptance.mn_param{2}.n)

      p = data.stage3.acceptance.mn_param{2};      
      
      text(x0+0.5,y0+nl1*dy,...
	   sprintf('|\\delta c_{2%i}|:  %5.3f \\pm %5.3f',...
		   p.n(i),p.cr(i),p.cr_err(i)));
      nl1 = nl1 + 1;
      text(x0+0.5,y0+nl1*dy,...
	   sprintf('\\phi_{2%i}:  %5.1f \\pm %5.1f',...
		   p.n(i),p.cphi(i),p.cphi_err(i)));
      nl1 = nl1 + 1;
    end
    
    nl2 = 2;
    
    nparm = data.stage3.acceptance.nparam_shape;
    
    for i = nparm+1:numel(data.stage3.acceptance.param)
      text(x0+0.5,y0+nl2*dy,...
	    sprintf('B_{%i}:  %6.1f \\pm %6.1f',...
		    i-nparm-1,data.stage3.acceptance.param(i),...
		    data.stage3.acceptance.param_err(i)));
	nl2 = nl2 + 1;
    end
    
    
  end
    
%  text(x0,y0,txt1,'HorizontalAlignment','left',...
%       'VerticalAlignment','top','Units','Normalized');
  
%  text(x0+0.5,y0,txt2,'HorizontalAlignment','left',...
%       'VerticalAlignment','top','Units','Normalized');
  
  title('Acceptance Summary');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW Az/El of Individual Runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_obs(data)
  tlab = 'Run Az/El';
  
  disp('draw_obs');
  
  wob = [];

  
  maxy = -inf;
  miny = inf;
  maxx = -inf;
  minx = inf;
  
  
  for irun = 1:numel(data.stage3.run)
    
    phi = data.stage3.run{irun}.wobble_phi_deg;    
    [rte,index]=min(abs(wob-phi));
    if(isempty(rte) || rte > 1)
      wob = [wob phi];
    end
  end
  
  wob = sort(wob);
  
  colors = { 'k' 'r' 'g' 'b'};  
  legendtext = {};
  txt = {};
  
  for iwob = 1:numel(wob)
    
    legendtext = { legendtext{:} sprintf('\\phi = %3.0f^\\circ',wob(iwob)) };
    
    x = [];
    y = [];
    nwob = 0;
    for irun = 1:numel(data.stage3.run)
      if(data.stage3.run{irun}.wobble_phi_deg == wob(iwob))    

	el = 90-data.stage3.run{irun}.zn_mean_deg;
	az = data.stage3.run{irun}.az_mean_deg;
	
	% Subtract 180 deg for sources that transit in the north
	if(data.stage3.src_dec_rad > 0.53537 && az>180)
	  az = az-360;
	end
	
	x = [x az];
	y = [y el];    
      
	maxy = max([maxy el]);
	miny = min([miny el]);
	maxx = max([maxx az]);
	minx = min([minx az]);
	
	nwob = nwob+1;
      end
    end
      
    txt{numel(txt)+1}=...
	sprintf('N(\\phi = %3.0f^\\circ): %i',wob(iwob),nwob);
      
    draw = '';
    
    if(iwob <= 4)      
      plot(x,y,strcat(colors{iwob},'o'),'MarkerFaceColor',colors{iwob},...
	   'MarkerSize',3);
    end
    hold on;
  end
  
  dy = maxy-miny;
  dx = max([1 maxx-minx]);
  
  minx = minx-dx*0.08;
  maxx = maxx+dx*0.08;
  miny = miny-dy*0.05;
  maxy = maxy+dy*0.50;
  
  zaxis([minx maxx miny maxy]);
    
  zlegend(legendtext{:});
  
  xlabel('Azimuth [deg]')
  ylabel('Elevation [deg]')

  text(0.02,0.98,txt,'HorizontalAlignment','left',...
       'VerticalAlignment','top','Units','Normalized');
  
  title(tlab)    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW Acceptance Distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_acceptance(data)
  
  disp('draw_acceptance')
  
  tlab = 'FoV \theta^2 Distribution';
  
  fovr2_x = data.stage3.fovr2_counts_hist.contents.x;
  fovr2_y = data.stage3.fovr2_counts_hist.contents.y;
  
  fovr2_bkgnd_x = data.stage3.acceptance.fovr2_bkgnd_hist.contents.x;
  fovr2_bkgnd_y = data.stage3.acceptance.fovr2_bkgnd_hist.contents.y;
  
  fovr2_x = fovr2_x + data.stage3.fovr2_counts_hist.contents.bin_width/2.;
  fovr2_bkgnd_x = fovr2_bkgnd_x + ...
      data.stage3.fovr2_counts_hist.contents.bin_width/2.;
  
  x0 = data.stage3.fovr2_counts_hist.contents.lo_limit;
  xf = data.stage3.fovr2_counts_hist.contents.hi_limit;
  err = sqrt(data.stage3.fovr2_counts_hist.variance.y);
  
  plot(fovr2_x,fovr2_y,'.r')
  hold on;
  plot(fovr2_bkgnd_x,fovr2_bkgnd_y,'r')
%  plot(src_x,src_y,'g')
  
  xlim([x0 xf]);
  ylim([-inf 1.1*max(fovr2_y)]);
  
  errorbar(fovr2_x,fovr2_y,err,'.')
  hold off;
  
% $$$   txt = {};
% $$$   txt{numel(txt)+1}=...
% $$$       sprintf('Acceptance Model:  %s',data.stage3.acceptance.model);
% $$$ 
% $$$   if(data.stage3.acceptance.model == 'bessel2')
% $$$ 
% $$$     for i = 1:numel(data.stage3.acceptance.m0_param)
% $$$ 
% $$$       if(data.stage3.acceptance.m0_param(i) ~= 0)
% $$$ 	txt{numel(txt)+1}=...
% $$$ 	    sprintf('\\delta c_{0%i}:  %5.3f \\pm %5.3f',...
% $$$ 		    i,data.stage3.acceptance.m0_param(i),...
% $$$ 		    data.stage3.acceptance.m0_param_err(i));
% $$$       end
% $$$ 
% $$$     end
% $$$     
% $$$     for i = 1:numel(data.stage3.acceptance.m1r_param)
% $$$ 
% $$$       if(data.stage3.acceptance.m1r_param(i) ~= 0)
% $$$ 	txt{numel(txt)+1}=...
% $$$ 	    sprintf('|\\delta c_{1%i}|:  %5.3f \\pm %5.3f',...
% $$$ 		    i,data.stage3.acceptance.m1r_param(i),...
% $$$ 		    data.stage3.acceptance.m1r_param_err(i));
% $$$ 	txt{numel(txt)+1}=...
% $$$ 	    sprintf('\\phi_{1%i}:  %5.1f \\pm %5.1f',...
% $$$ 		    i,data.stage3.acceptance.m1phi_param(i),...
% $$$ 		    data.stage3.acceptance.m1phi_param_err(i));
% $$$       end
% $$$ 
% $$$     end
% $$$   end
% $$$   
% $$$ 
% $$$   
% $$$   
% $$$   text(0.5,0.98,txt,'HorizontalAlignment','left',...
% $$$        'VerticalAlignment','top','Units','Normalized');
  
  xlabel('\theta^2 [deg^2]')
  title(tlab)    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW FOV Acceptance Distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_fov_acceptance(data)
  
  disp('draw_fov_acceptance')
  
  tlab = 'FoV Acceptance - dP/d\Omega [deg^{-2}]';
  
  draw_hist2d(data,data.stage3.acceptance.fov_acceptance_hist,tlab);
  
  phi = (180-data.stage3.pangle_mean_deg)*pi/180.;
  
  a=axis;
  
  ax = [ 5 -5 -5 5 ];
  ay = [ 0  2 -2 0 ];
  ct = cos(phi);
  st = sin(phi);
  set(line((ax*ct-ay*st)*(a(2)-a(1))*0.003+a(2)-(a(2)-a(1))*0.05,...
	   (ay*ct+ax*st)*(a(2)-a(1))*0.003+a(4)-(a(4)-a(3))*0.15),...
      'Color','k')
  
  txt = sprintf('Mean Parallactic Angle: %5.1f',data.stage3.pangle_mean_deg);
  
  text(0.4,0.98,txt,'HorizontalAlignment','left',...
       'VerticalAlignment','top','Units','Normalized');
  
%  xlabel('\theta^2 [deg^2]')
%  title(tlab)    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW Acceptance Distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_projection(data,paxis)
    
  disp('draw_projection')
  
  dx = data.stage3.sky_counts_hist.xbincalc.bin_width;
  dy = data.stage3.sky_counts_hist.ybincalc.bin_width;
  
  x = data.stage3.sky_counts_hist.x+dx/2.;
  y = data.stage3.sky_counts_hist.y+dy/2.;
  z = data.stage3.sky_counts_hist.z;
  bkgnd_z = data.results.sky_bkgnd_hist.z;
  source_z = data.results.sky_source_hist.z;
  
  if(numel(source_z) == 0)
    source_z = bkgnd_z;
  end
  
  [nx ny] = size(z);
  
  if(paxis == 'x')
  
    px = zeros(size(x));
    bkgnd_px = zeros(size(x));
    source_px = zeros(size(x));
    
    for i=1:nx
      for j=1:ny	
	if(y(j) > -0.11 && y(j) < 0.11)
	  px(i) = px(i) + z(i,j);
	  bkgnd_px(i) = bkgnd_px(i) + bkgnd_z(i,j);
	  source_px(i) = source_px(i) + source_z(i,j);	  
	end
      end
    end
    
    x0 = data.stage3.sky_counts_hist.xlo;
    xf = data.stage3.sky_counts_hist.xhi;    
    
    px = rebin(px,4);
    err = sqrt(px);
    bkgnd_px = rebin(bkgnd_px,4);
    source_px = rebin(source_px,4);

    x = linspace(x0+2*dx,xf-2*dx,numel(px));
    
    errorbar(x,px,err,'.');
    hold on;
    plot(x,bkgnd_px,'r','LineWidth',1.0);
    plot(x,source_px,'g','LineWidth',1.0);

    xlim([x0 xf]);
    
  else
    
    px = zeros(size(y));
    bkgnd_px = zeros(size(y));
    source_px = zeros(size(y));
    
    for i=1:nx
      for j=1:ny	
	if(x(i) > -0.11 && x(i) < 0.11)
	  px(j) = px(j) + z(i,j);
	  bkgnd_px(j) = bkgnd_px(j) + bkgnd_z(i,j);
	  source_px(j) = source_px(j) + source_z(i,j);
	end
      end
    end
        
    x0 = data.stage3.sky_counts_hist.ylo;
    xf = data.stage3.sky_counts_hist.yhi;    

    px = rebin(px,4);
    err = sqrt(px);
    bkgnd_px = rebin(bkgnd_px,4);
    source_px = rebin(source_px,4);
    
    y = linspace(x0+2*dy,xf-2*dy,numel(px));
    
    errorbar(y,px,err,'.');
    hold on;
    plot(y,bkgnd_px,'r','LineWidth',1.0);
    plot(y,source_px,'g','LineWidth',1.0);
    
    xlim([x0 xf]);
    
  end

  if(paxis == 'x')
    title('Counts Projection X');
    xlabel('X [deg]');    
  else
    title('Counts Projection Y');
    xlabel('Y [deg]');    
  end
  

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW Theta Squared Distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_thetasq(data)
  
  disp('draw_thetasq')
  
  tlab = '\theta^2 Distribution';
  
  th2_x = data.stage3.th2_on_counts_hist.contents.x;
  th2_y = data.stage3.th2_on_counts_hist.contents.y;
  bkgnd_x = data.results.th2_bkgnd_hist.contents.x;
  bkgnd_y = data.results.th2_bkgnd_hist.contents.y;
  src_x = data.results.th2_source_hist.contents.x;
  src_y = data.results.th2_source_hist.contents.y;
  
  th2_x = th2_x + data.stage3.th2_on_counts_hist.contents.bin_width/2.;
  bkgnd_x = bkgnd_x + data.stage3.th2_on_counts_hist.contents.bin_width/2.;
  src_x = src_x + data.stage3.th2_on_counts_hist.contents.bin_width/2.;
  
  
  x0 = data.stage3.th2_on_counts_hist.contents.lo_limit;
  xf = data.stage3.th2_on_counts_hist.contents.hi_limit;
  err = sqrt(data.stage3.th2_on_counts_hist.variance.y);
  
  plot(th2_x,th2_y,'.r')
  hold on;
  plot(bkgnd_x,bkgnd_y,'r')
  plot(src_x,src_y,'g')
  
  xlim([x0 xf]);

  errorbar(th2_x,th2_y,err,'.')
  hold off;
  
  xlabel('\theta^2 [deg^2]')
  title(tlab)    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW Significance Distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_significance(data)
  
  disp('draw_significance')
  
  tlab = 'Significance Distribution';

  bin_width = data.results.significance_hist.contents.bin_width;
  
  sig_x = data.results.significance_excluded_hist.contents.x;
  sig_y = data.results.significance_excluded_hist.contents.y;
  sig_err = sqrt(data.results.significance_excluded_hist.variance.y);
  
  sig_x = sig_x + bin_width/2.;
  
  x0 = data.results.significance_excluded_hist.contents.lo_limit;
  xf = data.results.significance_excluded_hist.contents.hi_limit;

  
%  plot(sig_x,sig_y,'.r')
  errorbar(sig_x,sig_y,sig_err,'.')
  xlim([x0 xf]);
  hold on;
  
  x = -10:0.1:10;
  
  gauss_y = bin_width*sum(sig_y)/sqrt(2*pi)*exp(-x.^2/2.);
  plot(x,gauss_y,'r')
  
  
  xlabel('Significance [\sigma]')
  title(tlab)    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRUM SUMMARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_spectrum_summary(data)

  set(gca,'Box','on');
  set(gca,'XTick',[]);
  set(gca,'YTick',[]);
  if(0)
    p=get(gca,'Position');
    factor=0.1;
    dx=p(3)*factor;
    dy=p(4)*factor;
    p=[p(1)-2*dx p(2)-dy p(3)+2*dx p(4)+2*dy];
    set(gca,'Position',p)
  end
  
  nl1=0;
  nl2=0;
  x0=+0.05;
  y0=+0.94;
  dy=-0.06;
  dx=+0.06;
  
  
  txt1 = {};
  txt2 = {};
  
  text(x0,y0+nl1*dy,...
       sprintf('Spectrum Model:  %s',data.spectrum_fit.model));
  
  nl1 = nl1 + 1;
  
  text(x0,y0+nl1*dy-0.02,sprintf('dF/dE @ E_{0} [m^{-2}s^{-1}TeV^{-1}]:'));  
  text(x0+0.45,y0+nl1*dy,...
       sprintf(' %8.2d \\pm %8.2d',...
		data.spectrum_fit.dfde,...		
		data.spectrum_fit.dfde_err));
  nl1 = nl1 + 1;
  
  text(x0,y0+nl1*dy-0.02,...
       sprintf('E^{2}dF/dE @ E_{0} [m^{-2}s^{-1}TeV]:'));  
  text(x0+0.45,y0+nl1*dy,...
       sprintf(' %8.2d \\pm %8.2d',...
		data.spectrum_fit.e2dfde,...		
		data.spectrum_fit.e2dfde_err));
  nl1 = nl1 + 1;
  
  text(x0,y0+nl1*dy-0.01,'E_{0} [TeV]:');  
  text(x0+0.45,y0+nl1*dy,...
       sprintf(' %6.3f',...
	       10^data.spectrum_fit.log10_enorm));
  nl1 = nl1 + 1;
  
  text(x0,y0+nl1*dy,'\Gamma:');  
  text(x0+0.45,y0+nl1*dy,...
       sprintf(' %6.2f \\pm %6.2f',...
		data.spectrum_fit.param(2),...		
		data.spectrum_fit.param_err(2)));
  nl1 = nl1 + 1;
  
  sigma = erfinv(1-data.spectrum_fit.chi2_ml_pval)*sqrt(2);
  
  text(x0,y0+nl1*dy,'\chi^{2}/ndf:');  
  text(x0+0.45,y0+nl1*dy,...
       sprintf('%6.2f/%4d (p = %6.3f \\sigma = %6.2f)',...
	       data.spectrum_fit.chi2_ml,...
	       data.spectrum_fit.ndf,...
	       data.spectrum_fit.chi2_ml_pval,sigma));
  nl1 = nl1 + 1;
  
  text(x0,y0+nl1*dy,'F(>100 GeV) [m^{-2}s^{-1}]:');  
  text(x0+0.45,y0+nl1*dy,sprintf('%6.2d',data.spectrum_fit.flux100));
  nl1 = nl1 + 1;

  text(x0,y0+nl1*dy,'F(>316 GeV) [m^{-2}s^{-1}]:');  
  text(x0+0.45,y0+nl1*dy,sprintf('%6.2d',data.spectrum_fit.flux316));
  nl1 = nl1 + 1;
  
  text(x0,y0+nl1*dy,'F(>1 TeV) [m^{-2}s^{-1}]:');  
  text(x0+0.45,y0+nl1*dy,sprintf('%6.2d',data.spectrum_fit.flux1000));
  nl1 = nl1 + 1;
  
  title('Spectrum Fit Summary');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW Spectrum Fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_spectrum_fit(data)
  tlab = 'Spectrum Fit Likelihood';

  if(numel(data.spectrum_fit.lcont68) == 0 || ...
     numel(data.spectrum_fit.lcont90) == 0)
    return;
  end
  
  x68 = data.spectrum_fit.lcont68{1}.x;
  y68 = data.spectrum_fit.lcont68{1}.y;
  
  x90 = data.spectrum_fit.lcont90{1}.x;
  y90 = data.spectrum_fit.lcont90{1}.y;
  
  legendtext = {};    
  legendtext = { legendtext{:} '68% C.L.' };
  legendtext = { legendtext{:} '90% C.L.' };

  ymax = 1.3*max(y90);

  plot(x68,y68);
  hold on;
  plot(x90,y90,'--');
  hold on;

  ymin = max( [0 0.8*min(y90)] );
  
  ylim([ymin ymax]);

  zlegend(legendtext{:});
  
  xlabel('\Gamma')
  ylabel('dF/dE [m^{-2}s^{-1}TeV^{-1}]')
  
  labels = {'\phi_{0} [m^{-2}s^{-1}TeV^{-1}]' '\Gamma'};
  
  txt = {};
  txt1 = {};
  txt2 = {};
  
% $$$   if(numel(data.spectrum_fit.param) == 2)
% $$$     
% $$$     txt1{numel(txt1)+1}=sprintf('%10s:',labels{1});
% $$$     txt2{numel(txt2)+1}=...
% $$$ 	sprintf(' %6.2d \\pm %6.2d',...
% $$$ 		data.spectrum_fit.param(1),...		
% $$$ 		data.spectrum_fit.param_err(1));
% $$$ 
% $$$     
% $$$     txt1{numel(txt1)+1}=sprintf('%10s:',labels{2});
% $$$     txt2{numel(txt2)+1}=...
% $$$ 	sprintf(' %6.2f \\pm %6.2f',...
% $$$ 		data.spectrum_fit.param(2),...		
% $$$ 		data.spectrum_fit.param_err(2));
% $$$     
% $$$     txt1{numel(txt1)+1}=sprintf('%10s:','\chi^{2}/ndf');
% $$$     txt2{numel(txt2)+1}=sprintf('%6.2f/%4d',data.spectrum_fit.chi2_ml,...
% $$$ 				data.spectrum_fit.ndf);
% $$$     
% $$$   end
  

% $$$   text(0.02,0.98,txt1,'HorizontalAlignment','left',...
% $$$        'VerticalAlignment','top','Units','Normalized');
% $$$   
% $$$   text(0.3,0.98,txt2,'HorizontalAlignment','left',...
% $$$        'VerticalAlignment','top','Units','Normalized');
  
 
  plot(data.spectrum_fit.param(2),data.spectrum_fit.param(1),...
       '+b','LineWidth',0.5,'MarkerSize',5);
  
  plot(data.spectrum_fit.param(2),data.spectrum_fit.param(1),...
       'sb','LineWidth',0.2,'MarkerSize',2,'MarkerFaceColor','b');
  
  title(tlab)    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW Differential Spectrum 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_spectrum(data)
  tlab = 'Differential Spectrum';

  dx = data.spectrum.e2dfde_hist.contents.bin_width;  
  x = data.spectrum.e2dfde_hist.contents.x;  
  x = x + dx/2.;

  dfde = data.spectrum_fit.param(1);
  lambda = data.spectrum_fit.param(2);
  enorm = data.spectrum_fit.log10_enorm;
  
  x2 = -1:0.1:1;
  y2 = dfde*power(10,2*x2-lambda*(x2-enorm));  
  
  y = data.spectrum.e2dfde_hist.contents.y;
  err = sqrt(data.spectrum.e2dfde_hist.variance.y);
  
  errorbar(x,y,err,'.')
  hold on;
  set(gca,'yscale','log');
  
%  errorbare('vlogy',x,y,err,err,'.')
%  semilogy(x,y,'.')
  plot(x2,y2);

  for i = 1:numel(x)
    plot([x(i)-dx/2. ; x(i)+dx/2.],[y(i) y(i)],'-b')
    
    if(err(i) == 0)
      xy1 = [ x(i) y(i) ];
      xy2 = [ x(i) y(i)*0.1 ];
      h=arrow(xy1,xy2,2,'BaseAngle',90);
      set(h, 'FaceColor', 'b'); 
      set(h, 'EdgeColor', 'b'); 
    end
    
  end


  
  xlabel('Log_{10}(E/TeV)')
  ylabel('E^{2}dF/dE [m^{-2}s^{-1}TeV]')
  
  title(tlab)    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW Significance Distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_scope_participation(data)

  disp('draw_scope_participation');
  
  tlab = 'Scope Particpation';

  bin_width = data.stage3.has_image_hist.contents.bin_width;
  
  tel_combo = [ 3, 5, 6, 9, 10, 12, 7, 11, 13, 14, 15 ];
%  tel_combo_name = { 'T1T2', 'T1T3', 'T2T3', 'T1T4', 'T2T4',...
%		     'T3T4', 'T1T2T3', 'T1T2T4', 'T1T3T4', ...
%		     'T2T3T4', 'T1T2T3T4' };
  
  tel_combo_name = { '12', '13', '23', '14', '24',...
		     '34', '123', '124', '134', ...
		     '234', '1234' };
  
  has_image_y = zeros(size(tel_combo));
  used_in_reconstruction_y = zeros(size(tel_combo));
  triggered_y = zeros(size(tel_combo));
  
  ymax = 0;
  
  for i=1:numel(data.stage3.has_image_hist.contents.y)
    
    n = sum(bitget(i,1:16));

    x = data.stage3.has_image_hist.contents.x(i);
    y = data.stage3.has_image_hist.contents.y(i);
    
    ymax = max( [y ymax] );
    
    for j = 1:numel(tel_combo)
      if(x == tel_combo(j))
	has_image_y(j) = y;
      end
    end

    
    
    
%    bitget(i,1:16)
%    if(n > 1)
%      
%      has_image_y = [has_image_y data.stage3.has_image_hist.contents.y(i)];
%      
%    end
    
  end
  
  for i=1:numel(data.stage3.used_in_reconstruction_hist.contents.y)
    
    x = data.stage3.used_in_reconstruction_hist.contents.x(i);
    y = data.stage3.used_in_reconstruction_hist.contents.y(i);
    
    ymax = max( [y ymax] );
    
    for j = 1:numel(tel_combo)
      if(x == tel_combo(j))
	used_in_reconstruction_y(j) = y;
      end
    end
  end
  
  for i=1:numel(data.stage3.triggered_hist.contents.y)
    
    x = data.stage3.triggered_hist.contents.x(i);
    y = data.stage3.triggered_hist.contents.y(i);
    
    ymax = max( [y ymax] );
    
    for j = 1:numel(tel_combo)
      if(x == tel_combo(j))
	triggered_y(j) = y;
      end
    end
  end
  
  
  used_in_reconstruction_y = [used_in_reconstruction_y 0];
  has_image_y = [has_image_y 0];
  triggered_y = [triggered_y 0];
  
    
  stairs( [-0.5:10.5], used_in_reconstruction_y, 'k' );
  hold on;
  stairs( [-0.5:10.5], has_image_y, 'b' );
  stairs( [-0.5:10.5], triggered_y, 'r' );
%  errorbar(sig_x,sig_y,sig_err,'.')
%  xlim([x0 xf]);
  hold on;
  
  
  zaxis([-0.5 10.5 -inf 1.2*ymax ])
  
  
  set(gca,'XTick',[0:10]);  % This automatically sets 
                         % the XTickMode to manual.
			 % Set the XTickLabels so that abbreviations for the
			 % months are used.

  set(gca,'XTickLabel',tel_combo_name)
  
  zlegend({'Used In Reconstruction','Has Image','Triggered'});
  
  
  xlabel('Scope Combination')
  title(tlab)    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW Cumulative Significance vs. Livetime
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_cumulative_significance(data)
  
  x = [];
  y1 = [];
  y2 = [];
  
  for i=1:numel(data.cumulative_run_results)
    
    x = [x data.cumulative_run_results{i}.livetime_min];
    y1 = [y1 data.cumulative_run_results{i}.significance];
    y2 = [y2 data.cumulative_run_results{i}.excess_rate];
  end
  
  plot(x,y1,'-s','MarkerSize',3,'MarkerFaceColor','b');
  
  xlabel('Livetime [min]')
  ylabel('Cumulative Significance [\sigma]');
  title('Cumulative Source Significance');
  
end

function draw_cumulative_excess(data)
  
  x = [];
  y = [];
  yerr = [];
  
  for i=1:numel(data.cumulative_run_results)    
    x = [x data.cumulative_run_results{i}.livetime_min];
    y = [y data.cumulative_run_results{i}.excess];
    yerr = [yerr data.cumulative_run_results{i}.excess_err];
  end
  
  errorbar(x,y,yerr,'-s','MarkerSize',3,'MarkerFaceColor','b');
  
  xlabel('Livetime [min]')
  ylabel('Cumulative Excess [Counts]');
  title('Cumulative Source Excess');
  
end

function draw_reconstructed_energy(data)

  x = data.spectrum_fit.on_hist.contents.x + ...
      data.spectrum_fit.on_hist.contents.bin_width/2.;
  on_y = data.spectrum_fit.on_hist.contents.y;
  on_yerr = sqrt(data.spectrum_fit.on_hist.variance.y);
  
  mu_on_y = data.spectrum_fit.mu_on_hist.contents.y;
  mu_bkgnd_y = data.spectrum_fit.mu_bkgnd_hist.contents.y;
    
  errorbar(x,on_y,on_yerr,'.');
  hold on;
  plot(x,mu_bkgnd_y,'-r');
  plot(x,mu_on_y,'-g');
  
  xlim([-1.5 1.5]);
  
  txt1 = {};
  txt2 = {};
  
  sigma = erfinv(1-data.spectrum_fit.chi2_ml_pval)*sqrt(2);
  
  txt1{numel(txt1)+1}=sprintf('%10s:','\chi^{2}/ndf');
  txt2{numel(txt2)+1}=sprintf('%6.2f/%4d (p = %6.3f \\sigma = %6.2f)',...
			      data.spectrum_fit.chi2_ml,...
			      data.spectrum_fit.ndf,...
			      data.spectrum_fit.chi2_ml_pval,sigma);
    
  text(0.02,0.98,txt1,'HorizontalAlignment','left',...
       'VerticalAlignment','top','Units','Normalized');
  
  text(0.3,0.98,txt2,'HorizontalAlignment','left',...
       'VerticalAlignment','top','Units','Normalized');
  
  title('ON Source Reconstructed Energy Distribution');
  
  xlabel('Reconstructed Energy [log_{10}(E/TeV)]');
  ylabel('Events');
  
  zlegend({'Data','Background','Signal+Background'});
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW Significance Distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = rebin(varargin)

%   if(nargin >= 4)

%     x = varargin{1};
%     y = varargin{2};
%     err = varargin{3};    
%     nbin = varargin{4};
    
%     x2 = [];
%     y2 = [];
%     err2 = [];
    
%     for i = 1:nbin:numel(y)
%       xavg = 0;
%       var = 0;
%       ysum = 0;
%       n = min([i+nbin numel(y)];

%       for j = i:n
% 	ysum = ysum + y(j);
% 	xavg = xavg + x(j);
% 	var = var + err(j)*err(j);
%       end
      
%       xavg = xavg/(n-i);
      
%       y2 = [y2 ysum];
%       x2 = [x2 xavg];
%       err2 = [err2 sqrt(var)];    
%     end
    
%     varargout{1} = x2;
%     varargout{2} = y2;
%     varargout{3} = err2;
%   else
    y = varargin{1};
    nbin = varargin{2};
    
    y2 = [];
    
    for i = 1:nbin:numel(y)
      
      ysum = 0;
      n = min([i+nbin numel(y)]);
      for j = i:n
	ysum = ysum + y(j);
      end
      
      y2 = [y2 ysum];
    end
    
    varargout{1} = y2;
    
%  end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW THIS SPACE INTENTIONALLY BLANK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_this_space_intentionally_blank
    set(gca,'Box','on');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    a=axis;
    text((a(1)+a(2))*0.5,(a(3)+a(4))*0.5,...
         'This space intentionally blank',...
         'HorizontalAlignment','center','VerticalAlignment','middle');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW THIS PLOT NOT AVAILABLE IN THIS VERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_not_available
    global analyze_revision_str
    set(gca,'Box','on');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    a=axis;
    text((a(1)+a(2))*0.5,(a(3)+a(4))*0.5,...
         sprintf('This plot unavailable in version %s',...
                 analyze_revision_str),...
         'HorizontalAlignment','center','VerticalAlignment','middle');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW THIS SPACE INTENTIONALLY BLANK SQUARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_this_space_intentionally_blank_square
    set(gca,'Box','on');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    axis('equal')
    a=axis;
    text((a(1)+a(2))*0.5,(a(3)+a(4))*0.5,...
         'This space intentionally blank',...
         'HorizontalAlignment','center','VerticalAlignment','middle');
    p=get(gca,'Position');
    p(1)=p(1)-0.1*p(3);
    set(gca,'Position',p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW THIS PLOT NOT AVAILABLE IN THIS VERSION SQUARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_not_available_square
    global analyze_revision_str
    set(gca,'Box','on');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    axis('equal')
    a=axis;
    text((a(1)+a(2))*0.5,(a(3)+a(4))*0.5,...
         sprintf('This plot unavailable in version %s',...
                 analyze_revision_str),...
         'HorizontalAlignment','center','VerticalAlignment','middle');
    p=get(gca,'Position');
    p(1)=p(1)-0.1*p(3);
    set(gca,'Position',p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW NO SCOPE PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_no_scope_plot(t,iscope)
    set(gca,'Box','on');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    if(0)
        p=get(gca,'Position');
        factor=0.1;
        dx=p(3)*factor;
        dy=p(4)*factor;
        p=[p(1)-2*dx p(2)-dy p(3)+2*dx p(4)+2*dy];
        set(gca,'Position',p)
    end
    a=axis;
    text((a(1)+a(2))*0.5,(a(3)+a(4))*0.5,...
             sprintf('No data for T%1d available',iscope),...
             'HorizontalAlignment','center',...
             'VerticalAlignment','middle');
    title(t);
end

function draw_no_plot(t)
    set(gca,'Box','on');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    if(0)
        p=get(gca,'Position');
        factor=0.1;
        dx=p(3)*factor;
        dy=p(4)*factor;
        p=[p(1)-2*dx p(2)-dy p(3)+2*dx p(4)+2*dy];
        set(gca,'Position',p)
    end
    a=axis;
    text((a(1)+a(2))*0.5,(a(3)+a(4))*0.5,...
             'No data available',...
             'HorizontalAlignment','center',...
             'VerticalAlignment','middle');
    title(t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE SOURCE INFO ELEMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_source_info(data,file,src_name,analysis_name)
  src_name(src_name==' ') = '_';
  
  fp=fopen(file,'w');
  fprintf(fp,...
	  strcat('%-25s %-25s %-10s %-10s %10.3f %8.3f %8.3f %8.3f ',...
		 '%8.3f %8.3f'),...
	  src_name,...
	  analysis_name,...
	  data.stage3.src_ra_hms,data.stage3.src_dec_dms,...
	  data.stage3.livetime_min,...
	  data.results.significance,...
	  data.results.excess_rate,...
	  data.results.excess_rate_err,...
	  data.results.bkgnd_density_rate,...
	  data.results.bkgnd_density_rate_err);
  
  fprintf(fp,'\n');
  fclose(fp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

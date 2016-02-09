%MAT_DIAGNOSTICS   Plot VERITAS diagnostics information
%
%   Diagnostic plots are displayed in figures and printed to a postscript
%   file and to PNGs. The PS file is names 'x12345_diagnostics.ps' where
%   12345 is the run number. The PNG files are similarly named, but each
%   has a different page number, e.g. 'x12345_diagnostics_1.png' etc.
%
%   INPUTS:
%
%       file_or_data: filename to load and plot, or variable containing
%                     data previously loaded
%
%   EXAMPLE:
%
%       mat_diagnostics('/veritas/data/analyzed/d20060922/x30168.h5');
%
%       - OR -
%
%       data=loadh5('/veritas/data/analyzed/d20060922/x30168.h5');
%       mat_diagnostics(data)
%
%   AUTHOR:
%
%       Stephen Fegan - sfegan@astro.ucla.edu - 2006-09-22
%
%   VERSION:
%
%       $Id: mat_diagnostics.m,v 3.73 2010/10/20 02:13:17 matthew Exp $
%
%   TODO:
%
%       - Add pedestal distribution
%       - Add NSB with cutoff
%       - Add noise plot in RA/Dec
%       - Add plots of our reconstruction parameters
%       - Add use of NSpace cuts
%       - Add autocorrelation plots
%       - Add muon yield
%       - Add plot of channel-by-channel autocorrelation tau parameter
function mat_diagnostics(file_or_data,varargin)

    global l3_label l3_plotspec
    l3_plotspec = 'k';
    l3_label = 'L3';
    
    global tel_label tel_plotspec
    tel_plotspec = { 'b' 'r' 'g' 'm' };
    tel_label = { 'T1', 'T2', 'T3', 'T4' };
    
    global spect_Qcut spect_F_ref spect_gamma_ref
    spect_Qcut = 4;
    spect_F_ref = 1;
    spect_gamma_ref = 1;

    global analyze_revision_str analyze_revision
    analyze_revision_str = '1.51';
    analyze_revision = 1051;
    
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
    fprintf(1,'Drawing plots for %d ... ',data.run_number);
    
    if(isfield(file_or_data,'stage2'))
        data = data.stage2;
    end

    sheets = [];
    if(nargin>1)
        sheets = varargin{1};
    end
    
    run_info_file = '';
    if(nargin>2)
        run_info_file = varargin{2};
    end
    
    if(isfield(data.analyze,'revision'))
        analyze_revision = data.analyze.revision;
        analyze_revision_str = data.analyze.revision_string;
    end
    
    close all
    printfile=sprintf('x%d_diagnostics.ps',data.run_number);
    pngfile=sprintf('x%d_diagnostics_%%d.png',data.run_number);
    
    revision = '$Revision: 3.73 $';
    version = revision(12:length(revision)-2);
    sheettitle=sprintf(strcat('VERITAS analyze v%s/diagnostics v%s',...
                              ' - Run %d on %s/%s from %s - %s@%s'),...
                       analyze_revision_str,version,data.run_number,...
                       data.observation.name,...
                       data.observation.mode_string,...
                       data.run_start_time_string(1:19),...
                       data.analyze.user,data.analyze.host);

    % ---------------------------------------------------------------------
    % Sheet 1
    % ---------------------------------------------------------------------

    if(newfigure(sheets))
        zsubplot(4,4,1);  draw_summary(data);
        zsubplot(4,4,2);  draw_l3(data);
        zsubplot(4,4,3);  draw_l2(data);
        zsubplot(4,4,4);  
        if(analyze_revision < 2022)
            draw_not_available
        else
            draw_multitel_hist(data,'median_l1_rate',0,inf,...
                               'Run time [min]','L1 rate [Hz]',...
                               'Median L1 (CFD) rate');
        end

        zsubplot(4,4,5);  draw_total_signal(data);
        zsubplot(4,4,6);  
        if(analyze_revision < 2022)
            draw_not_available
        else
            a=draw_fir_hist(data);
            Dy=20; dy=(a(4)-a(3));
            if(dy<Dy); 
                a(4)=a(4)+(Dy-dy)/2; a(3)=a(3)-(Dy-dy)/2; zaxis(a); 
            end
        end
        zsubplot(4,4,7);  
        if(analyze_revision < 2022)
            draw_not_available
        else
            a=draw_multitel_hist(data,'median_pedvar',0,inf,...
                                 'Run time [min]',...
                                 'Pedestal varience [DC]',...
                                 'Median start of trace pedestal RMS');        
            Dy=0.5; dy=(a(4)-a(3));
            if(dy<Dy); 
                a(4)=a(4)+(Dy-dy)/2; a(3)=a(3)-(Dy-dy)/2; zaxis(a); 
            end
        end
        zsubplot(4,4,8);  draw_dead(data);
                       
        zsubplot(4,4,9);  draw_nchan_hist(data,'camera_nimage',...
                                          'Channels in image',...
                                          'Number of channels in image');
        zsubplot(4,4,10); draw_nscope(data);
        zsubplot(4,4,11); draw_iscope(data);
        zsubplot(4,4,12); draw_cfd_thresholds(data);
        zsubplot(4,4,13); draw_ped_hist(data);
        zsubplot(4,4,14); draw_gain_hist(data);
        zsubplot(4,4,15); draw_suppressed_channels(data);
        zsubplot(4,4,16); draw_mean_higain(data);
        endfigure(sheettitle,printfile,pngfile,4,4);
    end
    
    % ---------------------------------------------------------------------
    % Sheet 2
    % ---------------------------------------------------------------------

    if(newfigure(sheets))
        zsubplot(4,4,1);  draw_elevation(data);
        zsubplot(4,4,2);  draw_azimuth(data);
        zsubplot(4,4,3);  draw_ra(data);
        zsubplot(4,4,4);  draw_dec(data);
        zsubplot(4,4,5);  draw_delta_t(data);
        zsubplot(4,4,6);  draw_delta_t_log(data);
        zsubplot(4,4,7);  draw_gps_diff(data);
        zsubplot(4,4,8);  draw_mean_logain(data);
        zsubplot(4,4,9);  draw_log_l2(data);
        zsubplot(4,4,10); draw_laser_sample_time(data);
        zsubplot(4,4,11); draw_laser_crate_time(data);
        zsubplot(4,4,12); draw_signal_projected_crate_time(data);
        zsubplot(4,4,13); draw_nimage_ntrig(data);
        zsubplot(4,4,14); draw_linear_total_signal(data);
        zsubplot(4,4,15); 
        draw_nchan_hist_linear(data,...
                               'camera_ntrigger_largest_region',...
                               'Number of channels',...
                      'Largest number of contiguous triggering channels');
        zsubplot(4,4,16); 
        draw_nchan_hist_linear(data,'camera_ntrigger',...
                               'Number of channels',...
                               'Number of triggering channels');
        endfigure(sheettitle,printfile,pngfile,4,4);
    end
        
    % ---------------------------------------------------------------------
    % Sheet 3
    % ---------------------------------------------------------------------

    if(newfigure(sheets))
        zsubplot(4,4,1);  draw_centroid_hist2d(data,1);
        zsubplot(4,4,2);  draw_centroid_hist2d(data,2);
        zsubplot(4,4,3);  draw_centroid_hist2d(data,3);
        zsubplot(4,4,4);  draw_centroid_hist2d(data,4);
        zsubplot(4,4,5);  draw_multitel_hist(data,'camera_length',0,1,...
                                            'Length [deg]','Events',...
                                            'Image length');
        zsubplot(4,4,6);  draw_multitel_hist(data,'camera_width',0,0.5,...
                                            'Width [deg]','Events',...
                                            'Image width');
        zsubplot(4,4,7); 
        draw_multitel_hist(data,'camera_xc',-Inf,Inf,...
                           'Centroid-X [deg]','Events',...
                            'X-coordinate of image centroid in camera');
        zsubplot(4,4,8); 
        draw_multitel_hist(data,'camera_yc',-Inf,Inf,...
                           'Centroid-Y [deg]','Events',...
                           'Y-coordinate of image centroid in camera');
        zsubplot(4,4,9);  draw_reconstructed_hist2d(data,0)
        zsubplot(4,4,10); draw_reconstructed_hist2d_simplecuts(data)
        zsubplot(4,4,11); draw_reconstructed_core_hist2d(data)
        zsubplot(4,4,12); draw_multitel_hist(data,'camera_psi',-180,180,...
                                            'Psi [deg]','Events',...
                                            'Orientation of image');
%        zsubplot(4,4,13);  draw_l2_traces(data,2);
%        zsubplot(4,4,14);  draw_l2_traces(data,2);
%        zsubplot(4,4,15);  draw_l2_traces(data,3);
%        zsubplot(4,4,16);  draw_l2_traces(data,4);
%        zsubplot(4,4,13); draw_l2_fft(data)
        zsubplot(4,4,13); draw_datum_event_number_hist(data)
        zsubplot(4,4,14); draw_trigger_event_number_hist(data)
        zsubplot(4,4,15); draw_missing_event_number_hist(data)
        zsubplot(4,4,16); draw_all_l2_traces(data);
        endfigure(sheettitle,printfile,pngfile,4,4);
    end
    
    % ---------------------------------------------------------------------
    % Sheet 4
    % ---------------------------------------------------------------------

    if(newfigure(sheets))
        zsubplot(4,4,1);  draw_mean_forced_higain(data);
        zsubplot(4,4,2);  draw_mean_forced_logain(data);
        zsubplot(4,4,3);  draw_gps_diff_mean_hist(data);
        zsubplot(4,4,4);  draw_gps_diff_rms_hist(data);
        zsubplot(4,4,5);  draw_this_space_intentionally_blank
        zsubplot(4,4,6);  draw_this_space_intentionally_blank
        zsubplot(4,4,7);  draw_this_space_intentionally_blank
        zsubplot(4,4,8);  draw_this_space_intentionally_blank
        zsubplot(4,4,9);  draw_this_space_intentionally_blank
        zsubplot(4,4,10); draw_this_space_intentionally_blank
        zsubplot(4,4,11); draw_this_space_intentionally_blank
        zsubplot(4,4,12); draw_this_space_intentionally_blank
        zsubplot(4,4,13); draw_crate_mean_higain(data,1);
        zsubplot(4,4,14); draw_crate_mean_higain(data,2);
        zsubplot(4,4,15); draw_crate_mean_higain(data,3);
        zsubplot(4,4,16); draw_crate_mean_higain(data,4);
        endfigure(sheettitle,printfile,pngfile,4,4);
    end

    % ---------------------------------------------------------------------
    % Sheet 5
    % ---------------------------------------------------------------------
    
    if(newfigure(sheets))
        zsubplot(4,4,1);  draw_telescope_signal_plots(data,1);
        zsubplot(4,4,2);  draw_telescope_signal_plots(data,2);
        zsubplot(4,4,3);  draw_telescope_signal_plots(data,3);
        zsubplot(4,4,4);  draw_telescope_signal_plots(data,4);
        zsubplot(4,4,5);  draw_tdc(data,1);
        zsubplot(4,4,6);  draw_tdc(data,2);
        zsubplot(4,4,7);  draw_tdc(data,3);
        zsubplot(4,4,8);  draw_tdc(data,4);
        zsubplot(4,4,9);  draw_single_datum_event_number_hist(data,1)
        zsubplot(4,4,10); draw_single_datum_event_number_hist(data,2)
        zsubplot(4,4,11); draw_single_datum_event_number_hist(data,3)
        zsubplot(4,4,12); draw_single_datum_event_number_hist(data,4)
        zsubplot(4,4,13); draw_gain_corrected_pedestal_rms(data,1)
        zsubplot(4,4,14); draw_gain_corrected_pedestal_rms(data,2)
        zsubplot(4,4,15); draw_gain_corrected_pedestal_rms(data,3)
        zsubplot(4,4,16); draw_gain_corrected_pedestal_rms(data,4)
        endfigure(sheettitle,printfile,pngfile,4,4);
    end
        
    % ---------------------------------------------------------------------
    % Sheet 5.1
    % ---------------------------------------------------------------------

if(0)    
    if(newfigure(sheets))
        zsubplot(4,4,1);  draw_gain_xy(data,1,1);
        zsubplot(4,4,2);  draw_gain_xy(data,2,1);
        zsubplot(4,4,3);  draw_gain_xy(data,3,1);
        zsubplot(4,4,4);  draw_gain_xy(data,4,1);
        zsubplot(4,4,5);  draw_gain_xy(data,1,2);
        zsubplot(4,4,6);  draw_gain_xy(data,2,2);
        zsubplot(4,4,7);  draw_gain_xy(data,3,2);
        zsubplot(4,4,8);  draw_gain_xy(data,4,2);
        zsubplot(4,4,9);  draw_pedvar_xy(data,1,1);
        zsubplot(4,4,10); draw_pedvar_xy(data,2,1);
        zsubplot(4,4,11); draw_pedvar_xy(data,3,1);
        zsubplot(4,4,12); draw_pedvar_xy(data,4,1);
        zsubplot(4,4,13); draw_pedvar_xy(data,1,2);
        zsubplot(4,4,14); draw_pedvar_xy(data,2,2);
        zsubplot(4,4,15); draw_pedvar_xy(data,3,2);
        zsubplot(4,4,16); draw_pedvar_xy(data,4,2);
        endfigure(sheettitle,printfile,pngfile,4,4);
    end
end        
    % ---------------------------------------------------------------------
    % Sheet 6
    % ---------------------------------------------------------------------

    if(newfigure(sheets))
        zsubplot(4,3,1);  draw_gain_camera(data,1);
        zsubplot(4,3,2);  draw_gain_camera(data,2);
        zsubplot(4,3,3);  draw_gain_camera(data,3);
        zsubplot(4,3,4);  draw_gain_camera(data,4);
        zsubplot(4,3,5);  draw_absgain_camera(data,1);
        zsubplot(4,3,6);  draw_absgain_camera(data,2);
        zsubplot(4,3,7);  draw_absgain_camera(data,3);
        zsubplot(4,3,8);  draw_absgain_camera(data,4);
        zsubplot(4,3,9);  draw_collection_eff_camera(data,1);
        zsubplot(4,3,10); draw_collection_eff_camera(data,2);
        zsubplot(4,3,11); draw_collection_eff_camera(data,3);
        zsubplot(4,3,12); draw_collection_eff_camera(data,4);
%        zsubplot(4,3,9);  draw_laser_time_camera(data,1);
%        zsubplot(4,3,10); draw_laser_time_camera(data,2);
%        zsubplot(4,3,11); draw_laser_time_camera(data,3);
%        zsubplot(4,3,12); draw_laser_time_camera(data,4);
        endfigure(sheettitle,printfile,pngfile,4,3);
    end
    
    % ---------------------------------------------------------------------
    % Sheet 7
    % ---------------------------------------------------------------------

    if(newfigure(sheets))
        zsubplot(4,3,1);  draw_pedvar_camera(data,1);
        zsubplot(4,3,2);  draw_pedvar_camera(data,2);
        zsubplot(4,3,3);  draw_pedvar_camera(data,3);
        zsubplot(4,3,4);  draw_pedvar_camera(data,4);
        zsubplot(4,3,5);  draw_median_current_camera(data,1);
        zsubplot(4,3,6);  draw_median_current_camera(data,2);
        zsubplot(4,3,7);  draw_median_current_camera(data,3);
        zsubplot(4,3,8);  draw_median_current_camera(data,4);
        zsubplot(4,3,9);  draw_suppressed_nslice(data,1)
        zsubplot(4,3,10); draw_suppressed_nslice(data,2)
        zsubplot(4,3,11); draw_suppressed_nslice(data,3)
        zsubplot(4,3,12); draw_suppressed_nslice(data,4)
        endfigure(sheettitle,printfile,pngfile,4,3);
    end

    % ---------------------------------------------------------------------
    % Sheet 8
    % ---------------------------------------------------------------------

    if(newfigure(sheets))        
        zsubplot(4,3,1);  draw_cfd_participation_camera(data,1)
        zsubplot(4,3,2);  draw_cfd_participation_camera(data,2)
        zsubplot(4,3,3);  draw_cfd_participation_camera(data,3)
        zsubplot(4,3,4);  draw_cfd_participation_camera(data,4)
        zsubplot(4,3,5);  draw_cfd_rate_camera(data,1)
        zsubplot(4,3,6);  draw_cfd_rate_camera(data,2)
        zsubplot(4,3,7);  draw_cfd_rate_camera(data,3)
        zsubplot(4,3,8);  draw_cfd_rate_camera(data,4)
        zsubplot(4,3,9);  draw_cfd_threshold_camera(data,1)
        zsubplot(4,3,10); draw_cfd_threshold_camera(data,2)
        zsubplot(4,3,11); draw_cfd_threshold_camera(data,3)
        zsubplot(4,3,12); draw_cfd_threshold_camera(data,4)
        endfigure(sheettitle,printfile,pngfile,4,3);
    end
    
    % ---------------------------------------------------------------------
    % Sheet 9
    % ---------------------------------------------------------------------

    if(newfigure(sheets))
        zsubplot(4,3,1);  draw_raw_max1_camera(data,1);
        zsubplot(4,3,2);  draw_raw_max1_camera(data,2);
        zsubplot(4,3,3);  draw_raw_max1_camera(data,3);
        zsubplot(4,3,4);  draw_raw_max1_camera(data,4);
        zsubplot(4,3,5);  draw_raw_top3_camera(data,1)
        zsubplot(4,3,6);  draw_raw_top3_camera(data,2)
        zsubplot(4,3,7);  draw_raw_top3_camera(data,3)
        zsubplot(4,3,8);  draw_raw_top3_camera(data,4)
        zsubplot(4,3,9);  draw_subthreshold_frequency_camera(data,1)
        zsubplot(4,3,10); draw_subthreshold_frequency_camera(data,2)
        zsubplot(4,3,11); draw_subthreshold_frequency_camera(data,3)
        zsubplot(4,3,12); draw_subthreshold_frequency_camera(data,4)
%        zsubplot(4,3,9);  draw_isolated_frequency_camera(data,1)
%        zsubplot(4,3,10); draw_isolated_frequency_camera(data,2)
%        zsubplot(4,3,11); draw_isolated_frequency_camera(data,3)
%        zsubplot(4,3,12); draw_isolated_frequency_camera(data,4)
        endfigure(sheettitle,printfile,pngfile,4,3);
    end
    
    % ---------------------------------------------------------------------
    % Sheet 10
    % ---------------------------------------------------------------------

    if(newfigure(sheets))
        zsubplot(4,3,1);  draw_lo_gain_camera(data,1);
        zsubplot(4,3,2);  draw_lo_gain_camera(data,2);
        zsubplot(4,3,3);  draw_lo_gain_camera(data,3);
        zsubplot(4,3,4);  draw_lo_gain_camera(data,4);
	zsubplot(4,3,5);  draw_laser_time_camera(data,1);
        zsubplot(4,3,6);  draw_laser_time_camera(data,2);
        zsubplot(4,3,7);  draw_laser_time_camera(data,3);
        zsubplot(4,3,8);  draw_laser_time_camera(data,4);
 	zsubplot(4,3,9);  draw_this_space_intentionally_blank
	zsubplot(4,3,10);  draw_this_space_intentionally_blank
	zsubplot(4,3,11);  draw_this_space_intentionally_blank
	zsubplot(4,3,12);  draw_this_space_intentionally_blank
        endfigure(sheettitle,printfile,pngfile,4,3);
    end
        
    % ---------------------------------------------------------------------
    % Write Run Info
    % ---------------------------------------------------------------------
    
    if(~isempty(run_info_file))
        finish_run_info(data);
        write_run_info(run_info_file,run_info);
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

function endfigure(pagetitle,printfile,pngfile,nx,ny)
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
   if(0)
        print('-dpng',sprintf(pngfile,sheet));
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
% RAD2HMSSTRING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = rad2hmsstring(rad)
    hrs = mod(rad/pi*12+24,24);
    ifracsec = mod(round(hrs*3600*10),24*3600*10);
    str = sprintf('%02d:%02d:%02d.%01d',...
                  floor(ifracsec/3600/10),...
                  floor(mod(ifracsec/60/10,60)),...
                  floor(mod(ifracsec/10,60)),...
                  floor(mod(ifracsec,10)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAD2DMSSTRING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = rad2dmsstring(rad)
    deg = mod(rad/pi*180+360,360);
    if(deg>180)
        deg = deg-360;
    end
    ifracsec = mod(round(abs(deg)*3600*10),180*3600*10);
    if(deg>=0)
        str = sprintf('+%02d:%02d:%02d.%01d',...
                      floor(ifracsec/3600/10),...
                      floor(mod(ifracsec/60/10,60)),...
                      floor(mod(ifracsec/10,60)),...
                      floor(mod(ifracsec,10)));
    else
        str = sprintf('-%02d:%02d:%02d.%01d',...
                      floor(ifracsec/3600/10),...
                      floor(mod(ifracsec/60/10,60)),...
                      floor(mod(ifracsec/10,60)),...
                      floor(mod(ifracsec,10)));
    end
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

function m = zmean(varargin)
  
  m = 0;
  
  if(nargin == 1 && (~isempty(varargin{1})))  
    if(sum(~isnan(varargin{1})) > 0)
      m = mean(varargin{1}(~isnan(varargin{1})));
    end
  elseif(nargin == 2 && (~isempty(varargin{1})) ) 
    
    n = 0;
    s = 0;
    x = double(varargin{1});
    y = double(varargin{2});
    
    s = dot(x,y);
    n = sum(y);
    
    if(n > 0)
      m = s/n;
    end
    
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% $$$ function m = zmean(data)
% $$$     if(sum(~isnan(data)) > 0)
% $$$         m = mean(data(~isnan(data)));
% $$$     else
% $$$         m = 0;
% $$$     end
% $$$ end

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
    global run_info

    global analyze_revision
    if(analyze_revision < 1051)
        draw_not_available
        return
    end
    
    D=data.diagnostics;

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

    nscope = 0;
    run_info.nscope = 0;
    run_info.scope_mask = 0;
    mean_zn       = zeros(1,numel(D.t));
    mean_az       = zeros(1,numel(D.t));
    mean_ra       = zeros(1,numel(D.t));
    mean_dec      = zeros(1,numel(D.t));
    dev_ra        = zeros(1,numel(D.t));
    dev_dec       = zeros(1,numel(D.t));
    mean_gal_l    = zeros(1,numel(D.t));
    mean_gal_b    = zeros(1,numel(D.t));
    mean_moon_sep = zeros(1,numel(D.t));

    for iscope = 1:numel(D.t)      
        if(isstruct(D.t{iscope}))
	  run_info.nscope = run_info.nscope+1;
	  run_info.scope_mask = run_info.scope_mask + bitshift(1,iscope-1);
	  run_info.t{iscope}.laser_runno = 0;
	  
            smp = D.t{iscope}.scope_mean_position;
            nscope=nscope+1;
            mean_zn(iscope)    = zmean(smp.zn);
            mean_az(iscope)    = zmean(smp.az);
            mean_ra(iscope)    = zmean(smp.ra);
            mean_dec(iscope)   = zmean(smp.dec);
            dev_ra(iscope)     = ...
                sqrt(zvar(smp.ra) + zmean(smp.ra_dev.^2));
            dev_dec(iscope)    = ...
                sqrt(zvar(smp.dec) + zmean(smp.dec_dev.^2));
            mean_gal_l(iscope) = zmean(smp.l);
            mean_gal_b(iscope) = zmean(smp.b);
            if(isfield(smp,'moon_separation'))
                mean_moon_sep(iscope) = zmean(smp.moon_separation);
            else
                mean_moon_sep(iscope) = -1;
            end
        else
            mean_zn(iscope)    = 0;
            mean_az(iscope)    = 0;
            mean_ra(iscope)    = 0;
            mean_dec(iscope)   = 0;
            dev_ra(iscope)     = 0;
            dev_dec(iscope)    = 0;
            mean_gal_l(iscope) = 0;
            mean_gal_b(iscope) = 0;
            mean_moon_sep(iscope) = 0;
        end
    end
    
    run_info.mean_zn_vec       = mean_zn;
    run_info.mean_az_vec       = mean_az;
    run_info.mean_ra_vec       = mean_ra;
    run_info.mean_dec_vec      = mean_dec;
    run_info.dev_ra_vec        = dev_ra;
    run_info.dev_dec_vec       = dev_dec;
    run_info.mean_moon_sep_vec = mean_moon_sep;
    
    nevent = double(D.events_found);
    if(isfield(D,'gps_ticks_elapsed'))
        gpseltime = double(D.gps_ticks_elapsed)/1e7;
    end
    elapsedtime = double(D.l3_ticks_elapsed)/1e7;
    livetime = elapsedtime-double(D.l3_ticks_veto_both)/1e7;

    text(x0,y0+nl*dy,sprintf('Run number: %d',data.run_number));
    run_info.runno = data.run_number;
        
    laser_run_string = '';
    
    if(isfield(data.stage1,'laser'))
      
      laser_runs = [];
      laser_scopes = {};
      
      for iscope = 1:numel(data.stage1.laser.scope)
	
	if(isstruct(data.stage1.laser.scope{iscope}))	
	  
	  if(isfield(data.stage1.laser.scope{iscope},'runno'))	  
	    runno = data.stage1.laser.scope{iscope}.runno;
	  else
	    runno = data.stage1.laser.m_runno;
	  end
	  
	  index = find(laser_runs == runno );
	  if(numel(index) == 0)
	    laser_runs = [laser_runs runno];
	    laser_scopes{numel(laser_scopes)+1}=sprintf('%d',iscope);
	  else
	    laser_scopes{index(1)}=...
		strcat(laser_scopes{index(1)},sprintf('%d',iscope));	
	  end
	  
	  if(iscope <= numel(D.t))
	    run_info.t{iscope}.laser_runno = runno;
	  end
	end
      end
	
      laser_runs
      laser_scopes
      
      if(numel(laser_runs) == 1)
	laser_run_string = sprintf('%d',laser_runs(1));
      else
	for irun = 1:numel(laser_runs)
	  laser_run_string = ...
	      strcat(laser_run_string,...
		     sprintf(' %s/%d',laser_scopes{irun},...
			     laser_runs(irun)));	
	end
      end
    end
          
    nl=nl+1;
    text(x0,y0+nl*dy,sprintf('Laser Run number: %s',laser_run_string));
    
    nl=nl+1;
    text(x0,y0+nl*dy,sprintf('Run start: %s (%11.5f)',...
         data.run_start_time_string,data.run_start_time));
    run_info.run_start_time_string = data.run_start_time_string;
    run_info.run_start_time = data.run_start_time;

    if(analyze_revision >= 2001)
        nl=nl+1;
        text(x0,y0+nl*dy,sprintf('Observation: %s/%s',...
             data.observation.name,data.observation.mode_string),...
             'Interpreter','none');
        run_info.target = data.observation.name;
        run_info.mode = data.observation.mode_string;
    end

    nl=nl+1;
    text(x0,y0+nl*dy,sprintf('Events found: %d',nevent));

    nl=nl+1;
    text(x0+dx,y0+nl*dy,sprintf('processed: %d (%.2f%%)',...
         D.events_processed,double(D.events_processed)/nevent*100.0));
    run_info.nevents = D.events_processed;

    nl=nl+1;
    text(x0+dx,y0+nl*dy,sprintf('reconstructed: %d  (%.2f%%)',...
         D.events_reconstructed,...
         double(D.events_reconstructed)/nevent*100.0));

    nl=nl+1;
    text(x0+dx,y0+nl*dy,sprintf('written: %d (%.2f%%)',...
         D.events_written,double(D.events_written)/nevent*100.0));

    if(isfield(D,'gps_ticks_elapsed'))
        nl=nl+1;
        text(x0,y0+nl*dy,...
             sprintf('Elapsed time GPS [sec]: %.1f',gpseltime));
        run_info.gpseltime = gpseltime;
    end

    nl=nl+1;
    text(x0,y0+nl*dy,sprintf('Elapsed time L3 [sec]: %.1f',elapsedtime));
    run_info.elapsed = elapsedtime;

    nl=nl+1;
    text(x0,y0+nl*dy,...
         sprintf('Live time L3 [sec]: %.1f (%.2f%%)   Mean rate [Hz]: %.1f',...
                 livetime,livetime/elapsedtime*100,...
                 double(D.events_processed)/livetime));
    run_info.livetime = livetime;
    run_info.rate = double(D.events_processed)/livetime;

    nl=nl+1;
    text(x0,y0+nl*dy,sprintf('Mean elevation [deg]: %.1f',...
         90-mean(mean_zn)));
    if(~isnan(mean_zn))
        run_info.mean_el = 90-mean(mean_zn);
        run_info.mean_az = mean(mean_az);
    else
        run_info.mean_el = 0;
        run_info.mean_az = 0;
    end

    nl=nl+1;
    text(x0,y0+nl*dy,sprintf('Mean RA , DEC [hms,dms]: %s , %s',...
               rad2hmsstring(mean(mean_ra)/180*pi),...
               rad2dmsstring(mean(mean_dec)/180*pi)))
    run_info.mean_ra = mean(mean_ra);
    run_info.mean_dec = mean(mean_dec);

    nl=nl+1;
    text(x0,y0+nl*dy,sprintf('Mean Galactic l, b [deg]: %.2f , %.2f',...
         mean(mean_gal_l),mean(mean_gal_b)));
    run_info.mean_gal_l = mean(mean_gal_l);
    run_info.mean_gal_b = mean(mean_gal_b);

    if(isfield(D,'moon_phase'))
        nl=nl+1;
        if(isfield(D,'moon_dphase_dt'))
            waxwane='\downarrow';
            if(D.moon_dphase_dt>0)
               waxwane='\uparrow';
            end
            text(x0,y0+nl*dy,...
                 sprintf('Moon elevation [deg]: %.2f. Phase: %.0f%% (%s)',...
                 D.moon_el,round(D.moon_phase*100),waxwane));
        else
            text(x0,y0+nl*dy,...
                 sprintf('Moon elevation [deg]: %.2f. Phase: %.0f%%',...
                 D.moon_el,round(D.moon_phase*100)));
        end
        run_info.moon_phase = round(D.moon_phase*100);
        run_info.moon_angle = D.moon_angle;
        run_info.moon_el = D.moon_el;
    end
    
    title('Summary');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L3 RATE HISTOGRAM - LINEAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_l3(data)
    global run_info
    D=data.diagnostics;
    Ni=D.l3_rate_hist;
    % Calculate time scaling
    T=60*ones(size(Ni.x));
    if(isempty(T))
        run_info.rate_per_min_mean = 0;
        run_info.rate_per_min_std = 0;
        run_info.rate_constant_chi2dof = 0;
        run_info.rate_constant_prob = 0;
        return
    end

    Telapsed=elapsed_sec(data)/60;
    T(numel(T)) = (Telapsed-Ni.x(numel(T)))*60;
    [xs,ys]=plot_simplehist(Ni,[-inf inf],'',1./T,1,1);
    plot(xs,ys,'b');
    zaxis([min(xs) max(xs) 0 roundaxis(ys)]);
    xlabel('Event time [min]');
    ylabel('Rate [Hz]');
    title('L3 rate');
    
    % Calculate probability that rate is constant (apart from last bin)
    % by maximum Likelihood and the chi^2 distribution
    ti=double(Ni.x);
    Ni=double(Ni.y(ti>=0));
    if(numel(Ni)>2)
        Ni=Ni(1:numel(Ni)-1);
        Nbar=mean(Ni);
        Nstd=std(Ni);
        chi2=-2*(sum(Ni*log(Nbar) - Nbar - gammaln(Ni+1)) ...
                 + 1/2*log(2*pi*Nbar)*numel(Ni));
        prob=1-chi2cdf(chi2,numel(Ni)-1);
        chi2dof=chi2/(numel(Ni)-1);
%        Nbar_alt=numel(Ni)/sum(1./Ni);
%        Chi2_alt=sum((Ni-Nbar_alt).^2 ./ Ni);
    else
        if(~isempty(Ni))
            Nbar=mean(Ni);
            Nstd=std(Ni);
        else
            Nbar=0;
            Nstd=0;
        end
        chi2dof=0;
        prob=1;
    end
    text(0.02,0.02,sprintf(strcat('Mean rate [Hz]: %.1f. Constant fit:',...
                                 ' (\\chi^2/dof)^{1/2}: %.2f (prob: %.2g)'),...
                           Nbar/60,sqrt(chi2dof),prob),...
        'HorizontalAlignment','left',...
        'VerticalAlignment','bottom','Units','Normalized');

    run_info.rate_per_min_mean = Nbar/60;
    run_info.rate_per_min_std = Nstd/60;
    run_info.rate_constant_chi2dof = chi2dof;
    run_info.rate_constant_prob = prob;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L2 RATE HISTOGRAM - LINEAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_l2(data)
    global tel_label tel_plotspec
    global run_info
    maxy=0;
    isc=[];
    D=data.diagnostics;
    Telapsed=elapsed_sec(data)/60;
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            Ni=D.t{iscope}.l2_rate_hist;
	    if(~isempty(Ni.x))
                % Calculate time scaling
		T=60*ones(size(Ni.x));
		T(numel(T)) = (Telapsed-Ni.x(numel(T)))*60;
		[xs ys]=plot_simplehist(Ni,[-inf inf],'',1./T,1,1);
		maxy = max([maxy max(ys)]);
		plot(xs,ys,tel_plotspec{iscope});
		hold on;
		isc=cat(2,isc,iscope);
		run_info.t{iscope}.l2_rate_per_min_median = ...
                    zmedian(double(Ni.y)./T);
	    else
		run_info.t{iscope}.l2_rate_per_min_median = 0;
	    end
        else
            run_info.t{iscope}.l2_rate_per_min_median = 0;
        end
    end
    zaxis([-inf inf 0 roundaxis(maxy)]);
    xlabel('Event time [min]');
    ylabel('Rate [Hz]');
    title('L2 rates');
    zlegend(tel_label{isc});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIR HISTOGRAM - LINEAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = draw_fir_hist(data)
    global tel_label tel_plotspec
    global run_info
    a=[inf -inf inf -inf];
    D=data.diagnostics;

    fir_label = { };
    txt = {'FIR Mean/RMS'};
    
    run_info.fir0_mean = 0;
    run_info.fir0_dev = 0;
    
    if(numel(D.fir) >= 5 && isstruct(D.fir{5}) && (~isempty(D.fir{5}.time)))
      h.x=D.fir{5}.time;
      h.y=D.fir{5}.sky;
      [xs ys]=plot_simplehist(h);
      a(1) = min([a(1) min(xs)]);
      a(2) = max([a(2) max(xs)]);
      a(3) = min([a(3) min(ys)]);
      a(4) = max([a(4) max(ys)]);            
      plot(xs,ys,'k-');
      hold on;	
      fir_label = { 'FIR0' };
      txt{numel(txt)+1}=...
	  sprintf('FIR0:  %.1f/%.2f',zmean(h.y),zstd(h.y));
      
      run_info.fir0_mean = zmean(h.y);
      run_info.fir0_dev = zstd(h.y);
    end
    
    for iscope = 1:numel(D.t)
      if(isstruct(D.t{iscope}))
	run_info.t{iscope}.fir_mean = 0;
	run_info.t{iscope}.fir_dev = 0;
      end
    end

    if(isfield(D,'fir'))

      for iscope = 1:numel(D.t)
	  
	  if(numel(D.fir) >= iscope && ~isempty(D.fir{iscope}.time))
            h.x=D.fir{iscope}.time;	   	    
            h.y=D.fir{iscope}.sky;

	    txt{numel(txt)+1}=...
		sprintf('FIR%1d:  %.1f/%.2f',iscope,zmean(h.y),zstd(h.y));
	    
            if(iscope == 3)
	      h.y = h.y + 40;
	      fir_label = { fir_label{:} sprintf('FIR%d + 40',iscope) };
	    else
	      fir_label = { fir_label{:} sprintf('FIR%d',iscope) };
	    end
	    
	    run_info.t{iscope}.fir_mean = zmean(h.y);
	    run_info.t{iscope}.fir_dev = zstd(h.y);
	    
	    [xs ys]=plot_simplehist(h);
            a(1) = min([a(1) min(xs)]);
            a(2) = max([a(2) max(xs)]);
            a(3) = min([a(3) min(ys)]);
            a(4) = max([a(4) max(ys)]);            
            plot(xs,ys,tel_plotspec{iscope});
            hold on;
	  end
	end
    else
        for iscope = 1:numel(D.t)
            if((isstruct(D.t{iscope}))&&(isfield(D.t{iscope},'fir_sky')))
                y=D.t{iscope}.fir_sky;
                [xs ys]=plot_simplehist(y);
                a(1) = min([a(1) min(xs)]);
                a(2) = max([a(2) max(xs)]);
                a(3) = min([a(3) min(ys)]);
                a(4) = max([a(4) max(ys)]);            
                plot(xs,ys,tel_plotspec{iscope});
                hold on;
		fir_label = { fir_label{:} sprintf('FIR%d',iscope) };
            end
        end
    end
    maxy=a(4);
    maxy=roundaxis(maxy);
%    if(maxy>0)
%        zaxis([0 inf 0 maxy]);
%    end
    xlabel('Run time [min]');
    ylabel('Temperature [^\circC]');
    title('FIR sky temperature');
    zlegend(fir_label{:});
    text(0.02,0.98,txt,'HorizontalAlignment','left',...
         'VerticalAlignment','top','Units','Normalized');
    if(nargout==1)
        varargout{1}=a;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L2 RATE HISTOGRAM - LOG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_log_l2(data)
    global tel_label tel_plotspec
    maxx=0;
    maxy=0;
    isc=[];
    D=data.diagnostics;
    Telapsed=elapsed_sec(data)/60;
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            Ni=D.t{iscope}.l2_rate_hist;
	    if(~isempty(Ni.x))
                % Calculate time scaling
		T=60*ones(size(Ni.x));
		T(numel(T)) = (Telapsed-Ni.x(numel(T)))*60;
		[xs ys]=plot_simplehist(Ni,[-inf inf],'',1./T,1,1);
		maxx = max([maxx max(xs)]);
		maxy = max([maxy max(ys)]);
		semilogy(xs,ys,tel_plotspec{iscope});
		hold on;
		isc=cat(2,isc,iscope);
	    end
        end
    end
    zaxis([-inf inf 1 10^roundaxis(log10(maxy))]);
    xlabel('Event time [min]');
    ylabel('Rate [Hz]');
    title('L2 rates');
    zlegend(tel_label{isc});
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEAD TIME HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_dead(data)
    global tel_label tel_plotspec
    D=data.diagnostics;
    Ni=D.l3_ticks_veto_both_hist;
    Telapsed=elapsed_sec(data)/60;
    T=6e8*ones(size(Ni.x));
    if(~isempty(T))
        T(numel(T)) = (Telapsed-Ni.x(numel(T)))*6e8;
    end
    [xs1,ys1]=plot_simplehist(Ni,[-inf inf],'',100./T,1,1);
%    [xs2,ys2]=plot_simplehist(D.l3_ticks_veto_vdaq_hist,...
%                              [-inf inf],'',1/6000000,1,1);
%    [xs3,ys3]=plot_simplehist(D.l3_ticks_veto_lev3_hist,...
%                              [-inf inf],'',1/6000000,1,1);
%    plot(xs1,ys1,'r',xs2,ys2,'b',xs3,ys3,'g');
%    zaxis([min(xs1) max(xs1) 0 roundaxis([max(ys1) max(ys2) max(ys3)])]); 
%    plot(xs1,ys1,'r',xs2,ys2,'b');
%    zaxis([min(xs1) max(xs1) 0 roundaxis([max(ys1) max(ys2)])]); 
    plot(xs1,ys1,'k-');
    hold on;
    maxx = max(xs1);
    maxy = max(ys1);
    legendtext= { 'Effective' };
    D=data.diagnostics;
    for iscope = 1:numel(D.t)
        if((isstruct(D.t{iscope}))...
           &&(isfield(D.t{iscope},'scope_l3_ticks_veto_vdaq_hist'))...
           &&(~isempty(D.t{iscope}.scope_l3_ticks_veto_vdaq_hist.x)))
            Ni=D.t{iscope}.scope_l3_ticks_veto_vdaq_hist;
            T=6e8*ones(size(Ni.x));
            T(numel(T)) = (Telapsed-Ni.x(numel(T)))*6e8;
            [xs,ys]=plot_simplehist(Ni,[-inf inf],'',100./T,1,1);
            plot(xs,ys,tel_plotspec{iscope});
            maxx = max([maxx max(xs)]);
            maxy = max([maxy max(ys)]);
            legendtext = { legendtext{:} tel_label{iscope} };
        end
    end
    t='Deadtime rate';
    if(~isempty(legendtext))
        zaxis([0 maxx -roundaxis(maxy)*0.05 roundaxis(maxy)]);
        xlabel('Event time [min]');
        ylabel('Deadtime [%]');
        title(t);
%        zlegend('Effective','VDAQ','L3');
%        zlegend('Effective','VDAQ');
        zlegend(legendtext{:});
    else
        draw_no_plot(t)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOTAL SIGNAL HISTOGRAM - LOG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_total_signal(data)
    global run_info
    global spect_Qcut
    global tel_label tel_plotspec
    global analyze_revision

    maxx=0;
    maxy=0;
    isc=[];
    D=data.diagnostics;
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            [xs ys]=plot_simplehist(D.t{iscope}.camera_log10_charge);
            maxx = max([maxx max(xs)]);
            maxy = max([maxy max(ys)]);
            loglog(10.^xs,ys,tel_plotspec{iscope});
            hold on;
            isc=cat(2,isc,iscope);
        end
    end
    zaxis([1 10^maxx 0.5 10^roundaxis(log10(maxy))]);
    xlabel('Total signal [DC]');
    ylabel('Events');
    if(analyze_revision < 3032)
        title('Total signal in image');
    else
        title('Total signal in triggered images');
    end

    zlegend(tel_label{isc});

    elapsedtime = elapsed_sec(data);
    livetime = elapsedtime-double(D.l3_ticks_veto_both)/1e7;
    nl=1;

    t=sprintf(strcat('dF/dQ = F (Q/%gkDC)^{-\\gamma}',...
                     ' [event/s/kDC]'),10^(spect_Qcut-3));
    txt={ t };
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            Qx=double(D.t{iscope}.camera_log10_charge.x);
            Qy=double(D.t{iscope}.camera_log10_charge.y);
            Qm=Qx>=spect_Qcut;
            Qn=sum(Qy(Qm));
            gamma=Qn/sum(Qy(Qm).*(Qx(Qm)-spect_Qcut)*log(10))+1;
            gmo=gamma-1;
            F=Qn*gmo/10^(spect_Qcut-3)/livetime;
            nl=nl+1;
            t=sprintf('T%1d, F=%.3f, \\gamma=%.2f, F^{1/(\\gamma-1)}=%.3f',...
                      iscope,F,gamma,F^(1/(gamma-1)));
            txt{nl}=t;
            if(~isnan(F))
                run_info.t{iscope}.spect_F = F;
                run_info.t{iscope}.spect_gamma = gamma;
                run_info.t{iscope}.spect_throughput = F^(1/(gamma-1));
            else
                run_info.t{iscope}.spect_F = 0;
                run_info.t{iscope}.spect_gamma = 0;
                run_info.t{iscope}.spect_throughput = 0;
            end
        else
            run_info.t{iscope}.spect_F = 0;
            run_info.t{iscope}.spect_gamma = 0;
            run_info.t{iscope}.spect_throughput = 0;
        end
    end
    text(0.02,0.02,txt,'Units','normalized',...
        'HorizontalAlignment','left','VerticalAlignment','bottom');
%        'FontSize',get(gcf,'defaultTextFontSize')*0.8);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOTAL SIGNAL HISTOGRAM - LINEAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_linear_total_signal(data)
    global tel_label tel_plotspec
    global analyze_revision

    maxx=0;
    maxy=0;
    isc=[];
    D=data.diagnostics;
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            [xs ys]=plot_simplehist(D.t{iscope}.camera_log10_charge);
            maxx = max([maxx max(xs)]);
            maxy = max([maxy max(ys)]);
            semilogx(10.^xs,ys,tel_plotspec{iscope});
            hold on;
            isc=cat(2,isc,iscope);
        end
    end
    zaxis([1 10^maxx 0 roundaxis(maxy)]);
    xlabel('Total signal [DC]');
    ylabel('Events');
    if(analyze_revision < 3032)
        title('Total signal in image');
    else
        title('Total signal in triggered images');
    end
    zlegend(tel_label{isc});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N-SCOPE HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_nscope(data)
    D=data.diagnostics;
    x=(0:4)'; y=zeros(size(x))+0.01;
    d1=struct('x',x,'y',y);
    d2=struct('x',x,'y',y);
    d3=struct('x',x,'y',y);
    d4=struct('x',x,'y',y);
    d1.y(D.scope_nimage.x+1)=double(D.scope_nimage.y)+0.01;
    d2.y(D.scope_ntrigger.x+1)=double(D.scope_ntrigger.y)+0.01;
    d3.y(D.scope_nsent_l3.x+1)=double(D.scope_nsent_l3.y)+0.01;
    d4.y(D.scope_nhas_event.x+1)=double(D.scope_nhas_event.y)+0.01;
    [xs1 ys1]=plot_simplehist(d1,[-inf inf],'',1,1,1);
    [xs2 ys2]=plot_simplehist(d2,[-inf inf],'',1,1,1);
    %[xs3 ys3]=plot_simplehist(d3,[-inf inf],'',1,1,1);
    [xs4 ys4]=plot_simplehist(d4,[-inf inf],'',1,1,1);
    semilogy(xs2-0.5,ys2,'r',xs4-0.5+0.01,ys4,'b',xs1-0.5+0.02,ys1,'g');
    maxy=max([max(ys1) max(ys2) max(ys4)]);
    zaxis([-0.5 4.5 0.5 10^roundaxis(log10(maxy))]);
    set(gca,'XTick',0:4);
    xlabel('Number of telescopes');
    ylabel('Events');
    title('Histograms of number of telescopes');
    zlegend('Triggered','Has event','In image');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I-SCOPE HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_iscope(data)
    D=data.diagnostics;
    x=(1:4)'-0.5; y=zeros(size(x));
    d1=struct('x',x,'y',y);
    d2=struct('x',x,'y',y);
    d3=struct('x',x,'y',y);
    d4=struct('x',x,'y',y);
    for iscope=1:4
        if((numel(D.t)>=iscope)&&(isstruct(D.t{iscope})))
            d1.y(iscope) = double(D.t{iscope}.scope_image);
            d2.y(iscope) = double(D.t{iscope}.scope_trigger);
            d3.y(iscope) = double(D.t{iscope}.scope_sent_l3);
            d4.y(iscope) = double(D.t{iscope}.scope_has_event);
        end
    end
    [xs1 ys1]=plot_simplehist(d1,[-inf inf],'',1,1,1);
    [xs2 ys2]=plot_simplehist(d2,[-inf inf],'',1,1,1);
    %[xs3 ys3]=plot_simplehist(d3,[-inf inf],'',1,1,1);
    [xs4 ys4]=plot_simplehist(d4,[-inf inf],'',1,1,1);
    plot(xs2,ys2,'r',xs4+0.01,ys4,'b',xs1+0.02,ys1,'g');
    maxy=max([max(ys1) max(ys2) max(ys4)]);
    zaxis([0.5 4.5 0.0 roundaxis(maxy)]);
    set(gca,'XTick',1:4);
    xlabel('Telescope number');
    ylabel('Events');
    title('Histograms of telescope participation in events');
    zlegend('Triggered','Has event','In image');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N-CHANNEL HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_nchan_hist(data,thehist,thexlabel,thetitle)
    global tel_label tel_plotspec
    maxx=0;
    maxy=0;
    isc=[];
    D=data.diagnostics;
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            eval(sprintf('nchan.x=double(D.t{iscope}.%s.x);',thehist));
            eval(sprintf('nchan.y=double(D.t{iscope}.%s.y);',thehist));
            if(~isempty(nchan.x))
                if(nchan.x(1)==0)
                    nchan.x(1)=0.01;
                end
                [xs ys]=plot_simplehist(nchan);
                maxx = max([maxx max(xs)]);
                maxy = max([maxy max(ys)]);
                loglog(xs,ys,tel_plotspec{iscope});
                hold on;
                isc=cat(2,isc,iscope);
            end
        end
    end
    zaxis([0.5 maxx 0.5 10^roundaxis(log10(maxy))]);
    xlabel(thexlabel);
    ylabel('Events');
    title(thetitle);
    zlegend(tel_label{isc});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N-CHANNEL HISTOGRAM LINEAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_nchan_hist_linear(data,thehist,thexlabel,thetitle)
    global tel_label tel_plotspec
    maxx=0;
    maxy=0;
    isc=[];
    D=data.diagnostics;
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            eval(sprintf('nchan.x=double(D.t{iscope}.%s.x);',thehist));
            eval(sprintf('nchan.y=double(D.t{iscope}.%s.y);',thehist));
            if(~isempty(nchan.x))
                if(nchan.x(1)==0)
                    nchan.x(1)=0.01;
                end
                [xs ys]=plot_simplehist(nchan);
                maxx = max([maxx max(xs)]);
                maxy = max([maxy max(ys)]);
                semilogx(xs,ys,tel_plotspec{iscope});
                hold on;
                isc=cat(2,isc,iscope);
            end
        end
    end
    zaxis([0.5 maxx 0 roundaxis(maxy)]);
    xlabel(thexlabel);
    ylabel('Events');
    title(thetitle);
    zlegend(tel_label{isc});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN CRATE HIGH-GAIN PULSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_crate_mean_higain(data,iscope)
    global tel_plotspec
    global analyze_revision
    if(analyze_revision < 3063)
        draw_not_available
        return
    end
    t=sprintf('Mean signal per crate in each HIGH gain trace - T%d',iscope);
    if((iscope>numel(data.channel_info))...
        ||(~isstruct(data.channel_info{iscope}))...
        ||(iscope>numel(data.diagnostics.t))...
        ||(~isstruct(data.diagnostics.t{iscope})))
        draw_no_scope_plot(t,iscope);        
        return
    end
    crate_label={'Crate 1','Crate 2','Crate 3','Crate 4'};    
    D=data.diagnostics.t{iscope};
    C=data.channel_info{iscope};

    ncrate=max(C.crate)-min(C.crate)+1;
    maxx=0;
    maxy=0;
    icr=[];
    for icrate = 0:(ncrate-1);
        cmhgs=D.channel_mean_higain_sample;
        m1=C.crate==icrate;
        m2=~C.suppress_all_events;
        m3=~isnan(sum(cmhgs,2))';
        m=m1&m2&m3;
        s=sum(cmhgs(m,:));
        c=sum(m);
        [xs ys]=zstairs((s-min(s))/c);
        maxx = max([maxx max(xs)]);
        maxy = max([maxy max(ys)]);
        plot(xs-1,ys,tel_plotspec{icrate+1});
        hold on;
        icr=cat(2,icr,icrate+1);
    end
    zaxis([0 maxx-1 0 roundaxis(maxy)]);
    xlabel('Sample number');
    ylabel('Mean signal');
    title(t);
    zlegend(crate_label{icr});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN HIGH-GAIN PULSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_mean_higain(data)
    global tel_label tel_plotspec
    maxx=0;
    maxy=0;
    isc=[];
    D=data.diagnostics;
    for iscope = 1:numel(D.t)
        if((isstruct(D.t{iscope}))...
           &&(~isempty(D.t{iscope}.channel_mean_higain_sample)))
            cmhgs=D.t{iscope}.channel_mean_higain_sample;
            m1=~data.channel_info{iscope}.suppress_all_events;
            m2=~isnan(sum(cmhgs,2))';
            m=m1&m2;
            s=sum(cmhgs(m,:));
            c=sum(m);
            [xs ys]=zstairs((s-min(s))/c);
            maxx = max([maxx max(xs)]);
            maxy = max([maxy max(ys)]);
            plot(xs-1,ys,tel_plotspec{iscope});
            hold on;
            isc=cat(2,isc,iscope);
        end
    end
    zaxis([0 maxx-1 0 roundaxis(maxy)]);
    xlabel('Sample number');
    ylabel('Mean signal');
    title('Mean signal in each HIGH gain trace');
    zlegend(tel_label{isc});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN LOW-GAIN PULSE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_mean_logain(data)
    global tel_label tel_plotspec
    maxx=0;
    maxy=0;
    isc=[];
    D=data.diagnostics;
    for iscope = 1:numel(D.t)
        if((isstruct(D.t{iscope}))...
           &&(~isempty( D.t{iscope}.channel_mean_logain_sample)))
            cmlgs=D.t{iscope}.channel_mean_logain_sample;
            m1=~data.channel_info{iscope}.suppress_all_events;
            m2=~isnan(sum(cmlgs,2))';
            m=m1&m2;
            s=sum(cmlgs(m,:));
            c=sum(m);
            [xs ys]=zstairs((s-min(s))/c);
            maxx = max([maxx max(xs)]);
            maxy = max([maxy max(ys)]);
            plot(xs-1,ys,tel_plotspec{iscope});
            hold on;
            isc=cat(2,isc,iscope);
        end
    end
    zaxis([0 maxx-1 0 roundaxis(maxy)]);
    xlabel('Sample number');
    ylabel('Mean signal');
    title('Mean signal in each LOW gain trace');
    zlegend(tel_label{isc});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN HIGH-GAIN PULSE - FORCED EVENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_mean_forced_higain(data)
    global tel_label tel_plotspec
    maxx=0;
    maxy=0;
    isc=[];
    D=data.diagnostics;
    for iscope = 1:numel(D.t)
        if((isstruct(D.t{iscope}))...
           &&(~isempty(D.t{iscope}.channel_mean_higain_sample_no_l2)))
            cmhgs=D.t{iscope}.channel_mean_higain_sample_no_l2;
            m1=~data.channel_info{iscope}.suppress_all_events;
            m2=~isnan(sum(cmhgs,2))';
            m=m1&m2;
            s=sum(cmhgs(m,:));
            c=sum(m);
            [xs ys]=zstairs((s-min(s))/c);
            maxx = max([maxx max(xs)]);
            maxy = max([maxy max(ys)]);
            plot(xs-1,ys,tel_plotspec{iscope});
            hold on;
            isc=cat(2,isc,iscope);
        end
    end
    zaxis([0 maxx-1 0 roundaxis(maxy)]);
    xlabel('Sample number');
    ylabel('Mean signal');
    title('Mean signal in each HIGH gain trace - FORCED events');
    zlegend(tel_label{isc});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN LOW-GAIN PULSE - FORCED EVENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_mean_forced_logain(data)
    global tel_label tel_plotspec
    maxx=0;
    maxy=0;
    isc=[];
    D=data.diagnostics;
    for iscope = 1:numel(D.t)
        if((isstruct(D.t{iscope}))...
           &&(~isempty( D.t{iscope}.channel_mean_logain_sample_no_l2)))
            cmlgs=D.t{iscope}.channel_mean_logain_sample_no_l2;
            m1=~data.channel_info{iscope}.suppress_all_events;
            m2=~isnan(sum(cmlgs,2))';
            m=m1&m2;
            s=sum(cmlgs(m,:));
            c=sum(m);
            [xs ys]=zstairs((s-min(s))/c);
            maxx = max([maxx max(xs)]);
            maxy = max([maxy max(ys)]);
            plot(xs-1,ys,tel_plotspec{iscope});
            hold on;
            isc=cat(2,isc,iscope);
        end
    end
    zaxis([0 maxx-1 0 roundaxis(maxy)]);
    xlabel('Sample number');
    ylabel('Mean signal');
    title('Mean signal in each LOW gain trace - FORCED events');
    zlegend(tel_label{isc});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERAL MULTI-TELESCOPE DIAGNOSTICS HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = draw_multitel_hist(data,dataset,xmin,xmax,xlab,...
                                        ylab,tlab)
    global tel_label tel_plotspec
    a=[inf -inf inf -inf];
    isc=[];
    D=data.diagnostics;
    for iscope = 1:numel(D.t)
        if((isstruct(D.t{iscope}))&&(isfield(D.t{iscope},dataset)))
            eval(sprintf('y=D.t{%d}.%s;',iscope,dataset));
            [xs ys]=plot_simplehist(y);
            a(1) = min([a(1) min(xs)]);
            a(2) = max([a(2) max(xs)]);
            a(3) = min([a(3) min(ys)]);
            a(4) = max([a(4) max(ys)]);            
            plot(xs,ys,tel_plotspec{iscope});
            hold on;
            isc=cat(2,isc,iscope);
        end
    end
    maxy=a(4);
    maxy=roundaxis(maxy);
    if(maxy>0)
        zaxis([xmin xmax 0 maxy]);
    end
    xlabel(xlab);
    ylabel(ylab);
    title(tlab);
    zlegend(tel_label{isc});
    if(nargout==1)
        varargout{1}=a;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PEDESTAL VARIANCE HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_ped_hist(data)
    global tel_label tel_plotspec
    global run_info
    maxy=0;
    minx=-Inf;
    maxx=Inf;
    isc=[];
    
    txt = { 'Median Pedvar' };
    run_info.mean_pedvar = 0;
    run_info.mean_pedvar_pe = 0;
    
    for iscope = 1:numel(data.diagnostics.t)      
      run_info.t{iscope}.median_pedvar_pe = 0;
      run_info.t{iscope}.median_pedvar = 0;
    end
      
    for iscope = 1:numel(data.channel_info)
      if(isstruct(data.channel_info{iscope}))
	C=data.channel_info{iscope};
	m=data.channel_info{iscope}.has_pmt...
	  & ~(data.channel_info{iscope}.suppress_reason&1);
	
	eff  = C.pmt_efficiency(m);
	absgain = C.pmt_multiplier_gain(m);
	uvar = C.median_dev(m).^2;
	udev = sqrt(uvar(uvar>0));
	pevar = uvar./(eff.*absgain.^2);
	pedev = sqrt(pevar(pevar>0));
	
	run_info.t{iscope}.median_pedvar_pe = ...
	    median(pedev(~isnan(pedev)));
	run_info.t{iscope}.median_pedvar = median(udev(~isnan(udev)));
	  
	if(isnan(run_info.t{iscope}.median_pedvar))
	  run_info.t{iscope}.median_pedvar = 0;
	end
	
	if(isnan(run_info.t{iscope}.median_pedvar_pe))
	  run_info.t{iscope}.median_pedvar_pe = 0;
	end
	
	run_info.mean_pedvar = ...
	    run_info.mean_pedvar + run_info.t{iscope}.median_pedvar;
	run_info.mean_pedvar_pe = ...
	    run_info.mean_pedvar_pe + run_info.t{iscope}.median_pedvar_pe;
		
	txt{numel(txt)+1}=...
	    sprintf('T%1d:  %.2f [DC]/%.3f [PE]',iscope,...
		    run_info.t{iscope}.median_pedvar,...
		    run_info.t{iscope}.median_pedvar_pe);
		
	val=C.median_dev;
	val(val==0)=NaN;
	[xlo xhi]=find_limits(val,0.95,1.2);
	[y,x]=hist(val,min(val):.1:max(val));
	if(~isempty(x))
	  [xs,ys]=stairs(x,y);
	  minx = min([minx xlo]);
	  maxx = max([maxx xhi]);
	  maxy = max([maxy max(ys)]);
	  plot(xs,ys,tel_plotspec{iscope});
	  hold on;
	  isc=cat(2,isc,iscope);
	end
      end
    end
    
    run_info.mean_pedvar = run_info.mean_pedvar/run_info.nscope;
    run_info.mean_pedvar_pe = run_info.mean_pedvar_pe/run_info.nscope;
    
    zaxis([0 maxx 0 roundaxis(maxy)]);
    xlabel('Pedestal RMS [DC]');
    ylabel('Channels');
    title('Pedestal RMS');
    if(~isempty(isc))
        zlegend(tel_label{isc});
    end
    
    text(0.02,0.98,txt,'HorizontalAlignment','left',...
         'VerticalAlignment','top','Units','Normalized');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAIN HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_gain_hist(data)
    global tel_label tel_plotspec
    maxy=0;
    isc=[];
    
    txt = { 'Gain RMS' };
    
    for iscope = 1:numel(data.channel_info)
        if(isstruct(data.channel_info{iscope}))
            m=data.channel_info{iscope}.gain>0;
	    m2=data.channel_info{iscope}.gain>0 ...
	       & ~(bitand(data.channel_info{iscope}.suppress_reason,1));
            [y,x]=hist(-log10(data.channel_info{iscope}.gain(m)),...
                       -1:.0025:1);
	    	 
	    rms = zstd(data.channel_info{iscope}.gain(m2));
	    txt{numel(txt)+1} = sprintf('T%1d:  %.3f',iscope,rms);	    
	    
            xlo=max([1 find(y>0,1)-1]);
            xhi=min([numel(y) find(y>0,1,'last')+1]);
            [xs,ys]=zstairs(x(xlo:xhi),y(xlo:xhi));
            maxy = max([maxy max(ys)]);
            semilogx(10.^xs,ys,tel_plotspec{iscope});
            hold on;
            isc=cat(2,isc,iscope);
        end
    end
    zaxis([1/1.5 1.5 0 roundaxis(maxy)]);
    set(gca,'XTick',[0.5 0.6 0.7 0.8 0.9 1.0 1.2 1.4 1.6 1.8 2.0]);
    xlabel('Gain [DC]');
    ylabel('Channels');
    title('Channel gain');
    zlegend(tel_label{isc});
    
    text(0.02,0.98,txt,'HorizontalAlignment','left',...
         'VerticalAlignment','top','Units','Normalized');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CFD THRESHOLDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_cfd_thresholds(data)
    global run_info
    global tel_label tel_plotspec
    maxy=0;
    isc=[];
    D=data.diagnostics;
    txt = { };
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            cfd=D.t{iscope}.channel_cfd_threshold;
            ntr=D.t{iscope}.channel_is_triggered;
            mtr=ntr>10 & ~isinf(cfd) & ~isnan(cfd);
            cfdmed=median(cfd(mtr));
            cfdstd=std(cfd(mtr & cfd>cfdmed*0.5 & cfd<cfdmed*2.0));
            cfdwid=4;
            mcfd=(cfd<(cfdmed+cfdwid*cfdstd))&(cfd>(cfdmed-cfdwid*cfdstd));
            [y,x]=hist(cfd(mcfd & mtr),0:2:200);
            xlo=max([1 find(y>0,1)-1]);
            xhi=min([numel(y) find(y>0,1,'last')+1]);
            [xs,ys]=stairs(x(xlo:xhi),y(xlo:xhi));
            maxy = max([maxy max(ys)]);
            plot(xs,ys,tel_plotspec{iscope});
            hold on;
            isc=cat(2,isc,iscope);
            cfdmedrel=sort(abs(cfd(mtr)-cfdmed));
            if(numel(cfdmedrel)>10)
                cfdmedrel80=cfdmedrel(floor(numel(cfdmedrel)*0.80));
            else
                cfdmedrel80=nan;
            end
            txt{numel(txt)+1}=...
                sprintf('T%1d - %.1f +/- %.1f',iscope,cfdmed,cfdmedrel80);
            if(~isnan(cfdmed))
                run_info.t{iscope}.cfd_threshold = cfdmed;
            else
                run_info.t{iscope}.cfd_threshold = 0;
            end
            if(~isnan(cfdmedrel80))
                run_info.t{iscope}.cfd_threshold_err = cfdmedrel80;
            else
                run_info.t{iscope}.cfd_threshold_err = 0;
            end
        else
            run_info.t{iscope}.cfd_threshold = 0;
            run_info.t{iscope}.cfd_threshold_err = 0;
        end
    end
    zaxis([-Inf Inf 0 roundaxis(maxy)]);
    xlabel('Threshold [mV]');
    ylabel('Channels');
    title('CFD 50% trigger threshold');
    zlegend(tel_label{isc});
    text(0.02,0.98,txt,'HorizontalAlignment','left',...
         'VerticalAlignment','top','Units','Normalized');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUPPRESSED CHANNELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_suppressed_channels(data)
    global run_info

    xcirc=cos((0:5:360)/180*pi);
    ycirc=sin((0:5:360)/180*pi);
    r=0.6;
    r2=1.7;
    r3=1.2;
    dy=-r*(sin(90/180*pi)+sin(210/180*pi))/2;

    c090=cos(90/180*pi);
    s090=sin(90/180*pi);
    c210=cos(210/180*pi);
    s210=sin(210/180*pi);
    c330=cos(330/180*pi);
    s330=sin(330/180*pi);

    plot(xcirc+r*cos(90/180*pi),ycirc+r*sin(90/180*pi)+dy,'r',...
         xcirc+r*cos(210/180*pi),ycirc+r*sin(210/180*pi)+dy,'b',...
         xcirc+r*cos(330/180*pi),ycirc+r*sin(330/180*pi)+dy,'g');
    axis('equal');
    set(gca,'YLim',[-1.75 1.75]);
    title('Suppressed channel count');
    zlegend('Laser','Peds','Pad peds');

    set(gca,'Box','on');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    if(0)
        p=get(gca,'Position');
        factor=0.08;
        dx2=p(3)*factor;
        dy2=p(4)*factor;
        p=[p(1)-2*dx2 p(2)-2*dy2 p(3)+2*dx2 p(4)+2*dy2];
        set(gca,'Position',p)
    end
    
    sep=0.13;
    nscope=0;
    for iscope = 1:numel(data.channel_info)
        if(isstruct(data.channel_info{iscope}))
            nscope=nscope+1;
        end
    end

    a=axis; 
    xleft=0.95*a(1)+0.05*a(2);
    xtop=0.95*a(4)+0.05*a(3);
    c=0;
    for iscope = 1:numel(data.diagnostics.t)
        if((iscope<=numel(data.channel_info))...
           &&(isstruct(data.channel_info{iscope})))
            csc=data.channel_info{iscope}.suppressed_count;
            yoff = dy-(c-nscope/2)*sep;
            text(0,yoff,sprintf('%d',csc.all));
            text(c090*r*r2,s090*r*r2+yoff,sprintf('%d',csc.laser));
            text(c210*r*r2,s210*r*r2+yoff,sprintf('%d',csc.ped))
            text(c330*r*r2,s330*r*r2+yoff,sprintf('%d',csc.pad));
            text(-c090*r*r3,-s090*r*r3+yoff,sprintf('%d',csc.ped_pad));
            text(-c210*r*r3,-s210*r*r3+yoff,sprintf('%d',csc.pad_laser));
            text(-c330*r*r3,-s330*r*r3+yoff,sprintf('%d',csc.ped_laser));
            text(xleft,xtop-c*sep,sprintf('T%1d - %3d',iscope,csc.any));
            c=c+1;
            run_info.t{iscope}.nsuppressed = csc.any;
        else
            run_info.t{iscope}.nsuppressed = 0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ELEVATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_elevation(data)
    global tel_label tel_plotspec
    minx=-Inf;
    maxx=Inf;
    isc=[];
    D=data.diagnostics;
    for iscope = 1:numel(D.t)
        if((isstruct(D.t{iscope}))...
           &&(~isempty(D.t{iscope}.scope_mean_position.t)))
            x=D.t{iscope}.scope_mean_position.t/60;
            minx=min([minx min(x)]);
            maxx=max([maxx max(x)]);
            plot(x,90-D.t{iscope}.scope_mean_position.zn,...
                 tel_plotspec{iscope});
            hold on;
            isc=cat(2,isc,iscope);
        end
    end
    zaxis([minx maxx -Inf Inf]);
    xlabel('Event time [min]');
    ylabel('Elevation [deg]');
    title('Telescope elevation angle');
    if(~isempty(isc))
        zlegend(tel_label{isc});
    end

    for iscope = 1:numel(D.t)
        for jscope = iscope+1:numel(D.t)
            x = 0.03 + (iscope-1)*0.10;
            y = 0.02 + (numel(D.t)-jscope+1.25)*0.06;
            t = '-';
            if((isstruct(D.t{iscope}))...
                &&(isstruct(D.t{jscope}))...
                &&(~isempty(D.t{iscope}.scope_mean_position.t))...
                &&(~isempty(D.t{jscope}.scope_mean_position.t)))
                [t,ii,ij]=intersect(D.t{iscope}.scope_mean_position.t,...
                                    D.t{jscope}.scope_mean_position.t);
                azi=D.t{iscope}.scope_mean_position.az(ii)/180*pi;
                azj=D.t{jscope}.scope_mean_position.az(ij)/180*pi;
                zni=D.t{iscope}.scope_mean_position.zn(ii)/180*pi;
                znj=D.t{jscope}.scope_mean_position.zn(ij)/180*pi;
                ctheta=cos(zni).*cos(znj)+sin(zni).*sin(znj).*cos(azi-azj);
                theta=mean(acos(ctheta(~isnan(ctheta))));
            t=sprintf('%.2f^\\circ',theta/pi*180);
            end
            text(x,y,t,'HorizontalAlignment','left',...
                 'VerticalAlignment','top','Units','Normalized');                
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AZIMUTH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_azimuth(data)
    global tel_label tel_plotspec
    minx=-Inf;
    maxx=Inf;
    isc=[];
    D=data.diagnostics;
    for iscope = 1:numel(D.t)
        if((isstruct(D.t{iscope}))...
           &&(~isempty(D.t{iscope}.scope_mean_position.t)))
            x=D.t{iscope}.scope_mean_position.t/60;
            minx=min([minx min(x)]);
            maxx=max([maxx max(x)]);
            plot(x,D.t{iscope}.scope_mean_position.az,...
                 tel_plotspec{iscope});
            hold on;
            isc=cat(2,isc,iscope);
        end
    end
    zaxis([minx maxx -Inf Inf]);
    xlabel('Event time [min]');
    ylabel('Azimuth [deg]');
    title('Telescope azimuth angle');
    if(~isempty(isc))
        zlegend(tel_label{isc});
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIGHT ASCENSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_ra(data)
    global tel_label tel_plotspec
    minx=-Inf;
    maxx=Inf;
    miny=-Inf;
    maxy=Inf;
    mint=Inf;
    maxt=-Inf;
    isc=[];
    h=[];    
    txt= {};
    D=data.diagnostics;
    
    obs_ra=0;
    if(isfield(data.observation,'obs_ra_rad'))
      txt={ 'Mean Error/RMS [deg]' };
      obs_ra = data.observation.obs_ra_rad*180./pi;
      maxy=obs_ra;
      miny=obs_ra;
    end
    
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            t = D.t{iscope}.scope_mean_position.t/60;
            minx=min([minx min(t)]);
            maxx=max([maxx max(t)]);
	    mint=min([mint min(t)]);
            maxt=max([maxt max(t)]);
	    
            ra = D.t{iscope}.scope_mean_position.ra;
            radev = D.t{iscope}.scope_mean_position.ra_dev;
	    ramin = ra-radev;
	    ramax = ra+radev;
	    
	    miny=min([miny ramin]);
	    maxy=max([maxy ramax]);
	    
	    hh=plot(t,ra+radev,tel_plotspec{iscope},...
                    t,ra-radev,tel_plotspec{iscope});
            hold on;
            if(~isempty(hh))
                isc=cat(2,isc,iscope);
                h=cat(2,h,hh(1));
            end

	    if(isfield(data.observation,'obs_ra_rad'))
	      dra = zmean(ra)-obs_ra;
	      txt={txt{:} sprintf('T%d: %+07.5f %07.5f',iscope,dra,zmean(radev))};
	    end
	end
    end
    
    maxy = maxy + 0.5*(maxy-miny);
    miny = miny - 0.1*(maxy-miny);
    
    zaxis([minx maxx miny maxy]);
    xlabel('Event time [min]');
    ylabel('RA [deg]');
    title('Apparent telescope right ascension - 1sigma range');
    zlegend(h,tel_label{isc});
    
    if(isfield(data.observation,'obs_ra_rad'))
      plot([mint ; maxt],[obs_ra ; obs_ra],'--k','LineWidth',1)
    end
      
    text(0.02,0.98,txt,'Units','normalized',...
            'HorizontalAlignment','left','VerticalAlignment','top');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DECLINATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_dec(data)
    global tel_label tel_plotspec
    minx=-Inf;
    maxx=Inf;
    miny=-Inf;
    maxy=Inf;
    mint=Inf;
    maxt=-Inf;
    isc=[];
    h=[];    
    txt={};
    D=data.diagnostics;

    obs_dec=0;
    if(isfield(data.observation,'obs_dec_rad'))
      txt={ 'Mean Error/RMS [deg]' };
      obs_dec = data.observation.obs_dec_rad*180./pi;
      miny=obs_dec;
      maxy=obs_dec;
    end
    
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            t = D.t{iscope}.scope_mean_position.t/60;
            minx=min([minx min(t)]);
            maxx=max([maxx max(t)]);
	    mint=min([mint min(t)]);
            maxt=max([maxt max(t)]);
	    
            dec = D.t{iscope}.scope_mean_position.dec;
            decdev = D.t{iscope}.scope_mean_position.dec_dev;
	    decmin = dec-decdev;
	    decmax = dec+decdev;
	    
	    miny=min([miny decmin]);
	    maxy=max([maxy decmax]);
	    
            hh=plot(t,dec+decdev,tel_plotspec{iscope},...
                    t,dec-decdev,tel_plotspec{iscope});
            hold on;
            if(~isempty(hh))
                isc=cat(2,isc,iscope);
                h=cat(2,h,hh(1));
            end
	    
	    if(isfield(data.observation,'obs_dec_rad'))
	      ddec = zmean(dec)-obs_dec;
	      txt={txt{:} sprintf('T%d: %+07.5f %07.5f',iscope,ddec,zmean(decdev))};
	    end
        end
    end
    
    maxy = maxy + 0.5*(maxy-miny);
    miny = miny - 0.1*(maxy-miny);
    
    zaxis([minx maxx miny maxy]);
    xlabel('Event time [min]');
    ylabel('Dec [deg]');
    title('Apparent telescope declination - 1sigma range');
    zlegend(h,tel_label{isc});
    
    if(isfield(data.observation,'obs_dec_rad'))
      plot([mint ; maxt],[obs_dec ; obs_dec],'--k','LineWidth',1)
    end
      
    text(0.02,0.98,txt,'Units','normalized',...
            'HorizontalAlignment','left','VerticalAlignment','top');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELTA T HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_delta_t(data)
    global run_info
%    global analyze_revision
    D=data.diagnostics;

    x=D.l3_dt.x;
    y=D.l3_dt.y;
    [xs,ys] = zstairs(x,y);

    m = x>0.001 & x<0.025;
    if(sum(m) > 2)
        xf = x(m)+D.l3_dt.bin_width/2;
        yf = log(double(y(m)));
        p = polyfit(xf,yf,1);
        run_info.tau = -1/p(1);
        run_info.delta_t_fit=[p(1) exp(p(2)-log(D.l3_dt.bin_width))];
        semilogy(x,exp(polyval(p,x)),'--k',xs,ys,'b')
        zlegend(sprintf('exp(-t/%.2f ms)',run_info.tau*1000),'Data');
    else
        semilogy(xs,ys,'b')
    end

    maxx = 1;
    if(sum(D.l3_dt.y>2)>0)
       maxx = D.l3_dt.x(find(D.l3_dt.y>2,1,'last'))*1.05;
    end
    zaxis([0 maxx 0.5 Inf]);
    xlabel('Delta T [sec]');
    ylabel('Events');

    if((isfield(run_info,'delta_t_fit'))&&(~isnan(run_info.tau)))
%        if(analyze_revision < 3031)
%            run_info.delta_t_live = double(D.events_processed)*run_info.tau;
%        else
%            run_info.delta_t_live = double(D.events_l2+D.events_ped)*run_info.tau;
%        end
        run_info.delta_t_live = double(D.events_found)*run_info.tau;
        txt=sprintf('Livetime estimate [sec]: %.1f',run_info.delta_t_live);
        if(isfield(D,'gps_ticks_elapsed'))
            gpseltime = double(D.gps_ticks_elapsed)/1e7;
            livepct = run_info.delta_t_live/gpseltime*100.0;
            txt = strcat(txt,sprintf(' (%.1f%%)',livepct));
        end

        text(0.02,0.02,txt,'Units','normalized',...
            'HorizontalAlignment','left','VerticalAlignment','bottom');
    else
        run_info.delta_t_live = 0;
    end

    title('Event separation');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELTA T HISTOGRAM - LOG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_delta_t_log(data)
    global run_info
    D=data.diagnostics;
    x=D.l3_log10_dt.x;
    y=double(D.l3_log10_dt.y);
    m=y>0;
    [xs,ys] = zstairs(x(m),y(m));

    if(isfield(run_info,'delta_t_fit'))
        xx = 10.^(x+D.l3_log10_dt.bin_width/2);
        xl = 10.^x;
        xh = 10.^(x+D.l3_log10_dt.bin_width);
        p0 = run_info.delta_t_fit(1);
        p1 = run_info.delta_t_fit(2);
        yy = p1/p0 * (exp(xh*p0)-exp(xl*p0));
        loglog(xx,yy,'--k',10.^xs,ys,'b')
        zlegend(sprintf('exp(-t/%.2f ms)',run_info.tau*1000),'Data',4);
    else
        loglog(10.^xs,ys,'b')
    end

    if(sum(m)>0)
        zaxis([10^min(x(m)) 10^max(x(m)) 1 10^roundaxis(log10(y(m)))]);
    end
    xlabel('Delta T [sec]');
    ylabel('Events');
    title('Event separation');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPS CLOCK TIME DIFFERENCE HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_gps_diff(data)
    global l3_label l3_plotspec
    global tel_label tel_plotspec
    maxy=0;
    isc=[];
    txt={ 'Under/overflow' };
    D=data.diagnostics;
    legendtext = { };

    if(isfield(D,'gps_l3_event_dt'))
        [xs ys]=plot_simplehist(D.gps_l3_event_dt,[-inf inf],'b',1,1/1000);
        maxy = max([maxy max(ys)]);
        plot(xs,ys,l3_plotspec);
        hold on;
        legendtext=cat(2,legendtext,l3_label);
        underflow = 0;
        overflow = 0;
        if(isfield(D.gps_l3_event_dt,'underflow'))
            underflow = D.gps_l3_event_dt.underflow;
        end
        if(isfield(D.gps_l3_event_dt,'overflow'))
            overflow = D.gps_l3_event_dt.overflow;
        end
        txt={txt{:} sprintf('L3: %d/%d', underflow, overflow)};
    end
    
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            if(isfield(D.t{iscope},'gps_scope_event_dt'))
                gps=D.t{iscope}.gps_scope_event_dt;
            else
                gps=D.t{iscope}.gps_scope_l3_dt;
            end
            [xs ys]=plot_simplehist(gps,[-inf inf],'b',1,1/1000);
            maxy = max([maxy max(ys)]);
            plot(xs,ys,tel_plotspec{iscope});
            hold on;
            legendtext=cat(2,legendtext,tel_label{iscope});
            underflow = 0;
            overflow = 0;
            if(isfield(gps,'underflow'))
                underflow = gps.underflow;
            end
            if(isfield(gps,'overflow'))
                overflow = gps.overflow;
            end
            txt={txt{:} sprintf('T%d: %d/%d',iscope,underflow,overflow)};
        end
    end
    t='GPS time diffence';
    if(~isempty(legendtext))
        zaxis([-Inf Inf 0 roundaxis(maxy)]);
        xlabel('GPS time relative to designated event time [microsecond]');
        ylabel('Events');
        title(t);
        zlegend(legendtext);
        text(0.02,0.98,txt,'Units','normalized',...
            'HorizontalAlignment','left','VerticalAlignment','top');
%            'FontSize',get(gcf,'defaultTextFontSize')*0.8);
    else
        draw_no_plot(t)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPS CLOCK TIME MEAN DIFFERENCE EVENT NUMBER HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_gps_diff_mean_hist(data)
    global l3_label l3_plotspec
    global tel_label tel_plotspec

    t='Mean GPS time difference vs event number';
    D=data.diagnostics;
    legendtext = { };
    maxy=[];
    miny=[];
    
    if((isfield(D,'gps_l3_diff_mean_ev_hist'))...
       &&(~isempty(D.gps_l3_diff_mean_ev_hist)))
        gps=D.gps_l3_diff_mean_ev_hist;
        ev=(0:numel(gps)-1)*1000;
        [xs,ys]=zstairs(ev,gps);
        maxy = max(ys);
        miny = min(ys);
        plot(xs,ys,l3_plotspec);
        hold on;
        legendtext=cat(2,legendtext,l3_label);
    end
    
    for iscope = 1:numel(D.t)
        if((isstruct(D.t{iscope}))...
            &&(isfield(D.t{iscope},'gps_scope_diff_mean_ev_hist'))...
            &&(~isempty(D.t{iscope}.gps_scope_diff_mean_ev_hist)))
            gps=D.t{iscope}.gps_scope_diff_mean_ev_hist;
            ev=(0:numel(gps)-1)*1000;
            [xs,ys]=zstairs(ev,gps);
            maxy = max([maxy max(ys)]);
            miny = min([miny min(ys)]);
            plot(xs,ys,tel_plotspec{iscope});
            hold on;
            legendtext=cat(2,legendtext,tel_label{iscope});
        end
    end
    if(~isempty(legendtext))
        maxy=min([max([roundaxis(maxy) roundaxis(abs(miny))]) 4000]);
        zaxis([-Inf Inf -maxy maxy]);
        xlabel('Event number');
        ylabel('Mean GPS time difference [ns]');
        title(t);
        zlegend(legendtext);
    else
        draw_no_plot(t);
    end
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPS CLOCK TIME RMS DIFFERENCE EVENT NUMBER HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_gps_diff_rms_hist(data)
    global l3_label l3_plotspec
    global tel_label tel_plotspec

    t='RMS GPS time diffence vs event number';
    D=data.diagnostics;
    legendtext = { };
    maxy=[];
    
    if((isfield(D,'gps_l3_diff_rms_ev_hist'))...
       &&(~isempty(D.gps_l3_diff_rms_ev_hist)))
        gps=D.gps_l3_diff_rms_ev_hist;
        ev=(0:numel(gps)-1)*1000;
        [xs,ys]=zstairs(ev,gps);
        maxy = max(ys);
        semilogy(xs,ys+eps,l3_plotspec);
        hold on;
        legendtext=cat(2,legendtext,l3_label);
    end
    
    for iscope = 1:numel(D.t)
        if((isstruct(D.t{iscope}))...
            &&(isfield(D.t{iscope},'gps_scope_diff_rms_ev_hist'))...
            &&(~isempty(D.t{iscope}.gps_scope_diff_rms_ev_hist)))
            gps=D.t{iscope}.gps_scope_diff_rms_ev_hist;
            ev=(0:numel(gps)-1)*1000;
            [xs,ys]=zstairs(ev,gps);
            maxy = max([maxy max(ys)]);
            semilogy(xs,ys+eps,tel_plotspec{iscope});
            hold on;
            legendtext=cat(2,legendtext,tel_label{iscope});
        end
    end
    if(~isempty(legendtext))
        zaxis([-Inf Inf 1 10^roundaxis(log10(maxy))]);
        xlabel('Event number');
        ylabel('RMS GPS time difference [ns]');
        title(t);
        zlegend(legendtext);
    else
        draw_no_plot(t);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LASER TIME FROM SAMPLE START HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_laser_sample_time(data)
    global tel_label tel_plotspec
    maxy=0;
    isc=[];
    CI=data.channel_info;
    for iscope = 1:numel(CI)
        if(isstruct(CI{iscope}))
            m=~isnan(CI{iscope}.chantime);
            [y,x]=hist((CI{iscope}.chantime(m)+CI{iscope}.l2time(m))*2,...
                       -40:0.2:100);
            xlo=max([1 find(y>0,1)-1]);
            xhi=min([numel(y) find(y>0,1,'last')+1]);
            [xs,ys]=zstairs(x(xlo:xhi),y(xlo:xhi));
            maxy=max([maxy max(ys)]);
            plot(xs,ys,tel_plotspec{iscope});
            hold on;
            isc=cat(2,isc,iscope);
        end
    end
    zaxis([-Inf Inf 0 roundaxis(maxy)]);
    xlabel('Sample time [ns]');
    ylabel('Channels');
    title('Mean LASER time with respect to start of sample');
    zlegend(tel_label{isc});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LASER TIME FROM CRATE TIME HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_laser_crate_time(data)
    global tel_label tel_plotspec
    maxy=0;
    isc=[];
    CI=data.channel_info;
    for iscope = 1:numel(CI)
        if(isstruct(CI{iscope}))
            m=~isnan(CI{iscope}.chantime);
            [y,x]=hist((CI{iscope}.chantime(m)...
                        -CI{iscope}.cratetime(m))*2,...
                        -40:0.2:40);
            xlo=max([1 find(y>0,1)-1]);
            xhi=min([numel(y) find(y>0,1,'last')+1]);
            [xs,ys]=zstairs(x(xlo:xhi),y(xlo:xhi));
            maxy=max([maxy max(ys)]);
            plot(xs,ys,tel_plotspec{iscope});
            hold on;
            isc=cat(2,isc,iscope);
        end
    end
    zaxis([-Inf Inf 0 roundaxis(maxy)]);
    xlabel('Relative time [ns]');
    ylabel('Channels');
    title('Mean LASER time with respect mean crate time');
    zlegend(tel_label{isc});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGNAL TIME FROM PROJECTED CRATE TIME HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_signal_projected_crate_time(data)
    global tel_label tel_plotspec
    maxy=0;
    isc=[];
    D=data.diagnostics;
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            ns = -15:0.5:30;
            c = uint32(zeros(size(ns)));
            for jchan = 1:numel(D.t{iscope}.channel_tij)
                if(~isempty(D.t{iscope}.channel_tij{jchan}.x))
                    m=(ns>=min(D.t{iscope}.channel_tij{jchan}.x))...
                       &(ns<=max(D.t{iscope}.channel_tij{jchan}.x));
                    c(m) = c(m) + D.t{iscope}.channel_tij{jchan}.y;
                end
            end
            [xs,ys]=zstairs(ns,c);
            maxy=max([maxy max(ys)]);
            plot(xs,ys,tel_plotspec{iscope});
            hold on;
            isc=cat(2,isc,iscope);
        end
    end
    zaxis([-Inf Inf 0 roundaxis(maxy)]);
    xlabel('Relative time [ns]');
    ylabel('Events times channels');
    title('Signal time with respect to projected crate time');
    zlegend(tel_label{isc});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L2 RATE FFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_l2_fft(data)
    global analyze_revision
    if(analyze_revision < 1051)
        draw_not_available
        return
    end
    
    global tel_label tel_plotspec
    isc=[];
    D=data.diagnostics;
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            x=D.t{iscope}.l2_high_res_rate_hist.x;
            y=double(D.t{iscope}.l2_high_res_rate_hist.y);
            if(numel(x)>1)
                dt=x(2)-x(1);
                n=floor(1/dt);
                N=floor(numel(y)/n);
                yy=zeros(n,N);
                yy(:)=y(1:n*N);
                YY=fft(yy);
                PP=YY.*conj(YY);
                if(n==1)
                    power=PP;
                else
                    power=sum(PP,2);
                end
                frequency=(0:numel(power)-1)/(numel(power)-1)/dt;
                m=2:numel(power)/2;
                semilogy(frequency(m),power(m)/power(1),tel_plotspec{iscope});
                hold on;
                isc=cat(2,isc,iscope);
            end
        end
    end
    a=axis;
    rate=double(D.events_processed)/elapsed_sec(data);
    l3nyquist=rate/2;
    if(l3nyquist<0.9*a(2))
        line([l3nyquist l3nyquist],[a(3) a(4)],'LineStyle','--','Color','k');
        text(l3nyquist*1.02,sqrt(a(3)*a(4)),...
             'This region affected by low L3 rate',...
             'VerticalAlignment','top','HorizontalAlignment','center',...
             'Rotation',90);
    end
    xlabel('Frequency [Hz]');
    ylabel('Power relative to DC [events^2]');
    title('Fourier transform of L2 event rate');
    zlegend(tel_label{isc});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L2 HISTOGRAM OF NUMBER OF IMAGE MINUS TRIGGERED CHANNELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_nimage_ntrig(data)
    global analyze_revision
    if(analyze_revision < 2000)
        draw_not_available
        return
    end
    
    global tel_label tel_plotspec
    minx=0;
    maxx=0;
    maxy=0;
    isc=[];
    D=data.diagnostics;
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            [xs ys]=...
 plot_simplehist(D.t{iscope}.camera_nimage_minus_ntrigger_largest_region);
            minx = min([minx min(xs)]);
            maxx = max([maxx max(xs)]);
            maxy = max([maxy max(ys)]);
            plot(xs,ys,tel_plotspec{iscope});
            hold on;
            isc=cat(2,isc,iscope);
        end
    end
    zaxis([-30 30 0 roundaxis(maxy)]);
    a=axis;
    set(line([0 0],[a(3) a(4)]),...
        'LineStyle','--','LineWidth',0.01,'Color','k');
    xlabel('Difference in number of channels');
    ylabel('Events');
    title('Num. image channels minus num. in largest triggering region');
    zlegend(tel_label{isc});
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
% DRAW TELESCOPE TOTAL AND CLEANED SIGNAL PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_telescope_signal_plots(data,iscope)
    global analyze_revision
    t=sprintf('Total and cleaned signal - T%1d',iscope);
    D=data.diagnostics;
    if((iscope<=numel(D.t))&&(isstruct(D.t{iscope})))
        [xs ys]=plot_simplehist(D.t{iscope}.camera_log10_charge);
        maxx=ceil(max(D.t{iscope}.camera_log10_charge.x));
        if(isempty(maxx))
            maxx = 0;
        end
        maxy=max(ys);
        if((iscope<=numel(data.events.t))...
           &&(isstruct(data.events.t{iscope}))...
           &&(~isempty(data.events.t{iscope}.fp_N)))
            m=data.events.t{iscope}.fp_N>0;
            N=data.events.t{iscope}.fp_N(m);
            maxx2=ceil(max(log10(N)));
            [yh2 xh2]=hist(log10(N),-2:.1:maxx2);
            [xs2 ys2]=zstairs(xh2,yh2);
            
            if(analyze_revision < 3032)
                maxx=max([maxx maxx2]);
                maxy=max([maxy max(ys)]);
                loglog(10.^xs,ys,'b',10.^xs2,ys2,'r');
                zlegend('Total','Cleaned');
            else
               [xs1 ys1]=...
                  plot_simplehist(D.t{iscope}.camera_log10_charge_no_l2);
                maxx1=ceil(max(D.t{iscope}.camera_log10_charge_no_l2.x));
                if(isempty(maxx1))
                  maxx1 = 0;
                end
                maxy1=max(ys1);
                maxx=max([maxx maxx1 maxx2]);
                maxy=max([maxy maxy1 max(ys)]);
                loglog(10.^xs,ys,'b',10.^xs1,ys1,'g',10.^xs2,ys2,'r');
                zlegend('Total (L2)','Total (no L2)','Cleaned');
            end
        else
            if(analyze_revision < 3032)
                loglog(10.^xs,ys,'b');
                zlegend('Total');
            else
               [xs1 ys1]=...
                  plot_simplehist(D.t{iscope}.camera_log10_charge_no_l2);
                maxx1=ceil(max(D.t{iscope}.camera_log10_charge_no_l2.x));
                if(isempty(maxx1))
                  maxx1 = 0;
                end
                maxy1=max(ys1);
                maxx=max([maxx maxx1]);
                maxy=max([maxy maxy1]);
                loglog(10.^xs,ys,'b',10.^xs1,ys1,'g');
                zlegend('Total (L2)','Total (no L2)');
            end
        end
        zaxis([1 10^maxx 0.5 10^roundaxis(log10(maxy))]);
        title(t);
    else
        draw_no_scope_plot(t,iscope);
    end
    xlabel('Total signal [DC]');
    ylabel('Events');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW L2 CHANNEL TRACES (HI GAIN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_l2_traces(data,iscope)
    global tel_plotspec
    t=sprintf('Mean L2 channel traces - T%1d',iscope);
    D=data.diagnostics;
    if((iscope<=numel(D.t))&&(isstruct(D.t{iscope})))
        l2chan=unique(data.channel_info{iscope}.l2chan);
        l2chan=l2chan(l2chan<numel(data.channel_info{iscope}.l2chan));
        if(~isempty(l2chan))
            ltext={};
            for ichan=1:numel(l2chan)
                jchan=l2chan(ichan)+1;
                plot(D.t{iscope}.channel_mean_higain_sample(jchan,:),...
                     tel_plotspec{mod(ichan-1,numel(tel_plotspec))+1});
                hold on
                ltext{ichan}=sprintf('Ch. %d',jchan);
            end
            zlegend(ltext);
        else
            draw_no_scope_plot(t,iscope);
        end
    else
        draw_no_scope_plot(t,iscope);
    end
    title(t);
    xlabel('Sample number');
    ylabel('Mean sample [DC]');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW ALL L2 CHANNEL TRACES (HI GAIN) - SAVE REAL ESTATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_all_l2_traces(data)
    global tel_label tel_plotspec
    isc=[];
    iha=[];
    D=data.diagnostics;
    CI=data.channel_info;
    for iscope = 1:numel(D.t)
        if((isstruct(D.t{iscope}))&&(iscope<=numel(CI))...
           &&(isstruct(CI{iscope})))
            l2chan=unique(data.channel_info{iscope}.l2chan);
            l2chan=l2chan(l2chan<numel(data.channel_info{iscope}.l2chan));
            if(~isempty(l2chan))
                for ichan=1:numel(l2chan)
                    jchan=l2chan(ichan)+1;
                    h=plot(D.t{iscope}.channel_mean_higain_sample(jchan,:),...
                           tel_plotspec{mod(iscope-1,numel(tel_plotspec))+1});
                    hold on
                    if(ichan==1)
                        iha=cat(2,iha,h);
                    end
                end
                isc=cat(2,isc,iscope);
            end
        end
    end
    xlabel('Sample number');
    ylabel('Signal [dc]');
    title('Mean L2 channel traces');
    zlegend(iha,tel_label{isc});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW TDC DIFFERENCES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_tdc(data,iscope)
    global analyze_revision
    if(analyze_revision < 2022)
        draw_not_available
        return
    end

    global tel_label tel_plotspec
    
    txt = { '\Delta t Mean/RMS [ns]' };
    
    t=sprintf('TDC L2 arrival time differences - T%1d',iscope);
    isc=[];
    D=data.diagnostics;
    if((iscope<=numel(D.t))&&(isstruct(D.t{iscope})))
        for jscope = 1:numel(D.t{iscope}.scope_tdc_diff)
            if((iscope~=jscope)...
               &&(isstruct(D.t{iscope}.scope_tdc_diff{jscope})))
                [xs,ys]=zstairs(D.t{iscope}.scope_tdc_diff{jscope}.x,...
                                D.t{iscope}.scope_tdc_diff{jscope}.y);
                isc=cat(2,isc,jscope);
                plot(xs,ys,tel_plotspec{jscope});
                hold on
		
		mean = zmean(D.t{iscope}.scope_tdc_diff{jscope}.x,...
			     D.t{iscope}.scope_tdc_diff{jscope}.y);
		
		dev = zstd(D.t{iscope}.scope_tdc_diff{jscope}.x,...
			   D.t{iscope}.scope_tdc_diff{jscope}.y,-100,100);
		
		txt{numel(txt)+1} = ...
		    sprintf('T%1d:  %+6.2f/%.2f',jscope,mean,dev);
            end
        end
        zaxis([-100 100 0 Inf]);
        title(t);
        zlegend(tel_label{isc});
    else
        draw_no_scope_plot(t,iscope);
    end
    xlabel(sprintf('L2 arrival time difference T%d-Tn [ns]',iscope));
    ylabel('Events');
    
    text(0.02,0.98,txt,'HorizontalAlignment','left',...
         'VerticalAlignment','top','Units','Normalized');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW HAS DATUM EVENT FRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_datum_event_number_hist(data)
    global analyze_revision
    if(analyze_revision < 2022)
        draw_not_available
        return
    end
    
    global tel_label tel_plotspec
    t='Fraction of events for which telescope has datum';
    legendtext={};
    D=data.diagnostics;
    
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            nev=double(data.stage1.run_info.hi_event_num);
            dnev = double(D.t{iscope}.scope_has_event_ev_num_hist.bin_width);

            x = 0:dnev:nev;
            y = x*0;

            Nev=[];
            if(nev>(dnev-1))
                Nev = dnev-1;
                nev = nev-(dnev-1);
            end
            while(nev>dnev)
                Nev = cat(2,Nev,dnev);
                nev = nev-dnev;
            end
            Nev = cat(2,Nev,nev);

            X = double(D.t{iscope}.scope_has_event_ev_num_hist.x);
            if(~isempty(X))
                m = (x>=X(1))&(x<=X(numel(X)));
                y(m)=double(D.t{iscope}.scope_has_event_ev_num_hist.y);
            end
            
            [xs,ys]=zstairs(x,y./Nev);

            legendtext={ legendtext{:} tel_label{iscope} };
            plot(xs,ys,tel_plotspec{iscope});
            hold on
        end
    end
    if(~isempty(legendtext))
        zaxis([0 ceil(double(data.stage1.run_info.hi_event_num)/1000)*1000 ...
               -0.05 1.05]);
        xlabel(sprintf('Event number'));
        ylabel('Fraction of events with datum');
        title(t);
        zlegend(legendtext{:},4);
    else
        draw_no_plot(t)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW HAS DATUM EVENT FRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_single_datum_event_number_hist(data,iscope)
    global analyze_revision
    if(analyze_revision < 2022)
        draw_not_available
        return
    end
    
    t=sprintf('Fraction of events for which telescope has datum - T%d',...
              iscope);
    D=data.diagnostics;
    
    if((iscope<=numel(D.t))&&(isstruct(D.t{iscope})))
        nev=double(data.stage1.run_info.hi_event_num);
        dnev = double(D.t{iscope}.scope_has_event_ev_num_hist.bin_width);

        x = 0:dnev:nev;
        y = x*0;

        Nev=[];
        if(nev>(dnev-1))
            Nev = dnev-1;
            nev = nev-(dnev-1);
        end
        while(nev>dnev)
            Nev = cat(2,Nev,dnev);
            nev = nev-dnev;
        end
        Nev = cat(2,Nev,nev);

        X = double(D.t{iscope}.scope_has_event_ev_num_hist.x);
        if(~isempty(X))
            m = (x>=X(1))&(x<=X(numel(X)));
            y(m)=double(D.t{iscope}.scope_has_event_ev_num_hist.y);
        end
            
        [xs,ys]=zstairs(x,y./Nev);

        plot(xs,ys,'b');
        zaxis([0 ceil(double(data.stage1.run_info.hi_event_num)/1000)*1000 ...
               -0.05 1.05]);
        title(t);
    else
        draw_no_scope_plot(t,iscope);
    end
    xlabel(sprintf('Event number'));
    ylabel('Fraction of events with datum');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW GAIN CORRECTED PEDESTAL RMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_gain_corrected_pedestal_rms(data,iscope)
    global analyze_revision
    global run_info
    if(analyze_revision < 3005)
        draw_not_available
        return
    end
    
    t=sprintf('Gain corrected pedestal RMS - T%d',iscope);
    if((iscope<=numel(data.channel_info))&&...
       (isstruct(data.channel_info{iscope})))
        C=data.channel_info{iscope};

        m=data.channel_info{iscope}.has_pmt...
          & ~(data.channel_info{iscope}.suppress_reason&1);
        uvar = C.median_dev(m).^2;
        gain = C.pmt_multiplier_gain(m)/median(C.pmt_multiplier_gain(m));
        eff  = C.pmt_efficiency(m);
        zvar = 0;
        cvar = (uvar-zvar)./(eff.*gain.^2);
	
        udev = sqrt(uvar(uvar>0));
        cdev = sqrt(cvar(cvar>0));
        [xlo xhi]=find_limits(cat(2,udev,cdev),0.99,1.2);
	if(xhi > (xlo+0.25))
	  [uy,ux]=hist(udev(~isnan(udev)),xlo:.25:xhi);
          [cy,cx]=hist(cdev(~isnan(cdev)),xlo:.25:xhi);
          [usx,usy]=zstairs(ux,uy);
          [csx,csy]=zstairs(cx,cy);
        else
	  usx = [];
          usy = [];
	  csx = [];
          csy = [];
        end
        plot(usx,usy,'b',csx,csy,'r');
        legend('Unequalized','Equalized');
        title(t);
    else
        draw_no_scope_plot(t,iscope);
    end
    xlabel(sprintf('RMS [DC]'));
    ylabel('Number of channels');    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW TRIGGER EVENT FRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_trigger_event_number_hist(data)
    global analyze_revision
    if(analyze_revision < 2022)
        draw_not_available
        return
    end
    
    global tel_label tel_plotspec
    t='Fraction of events for which telescopes trigger';
    legendtext={};
    D=data.diagnostics;
    
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            nev=double(data.stage1.run_info.hi_event_num);
            dnev = double(D.t{iscope}.scope_trigger_ev_num_hist.bin_width);            

            x = 0:dnev:nev;
            y = x*0;

            Nev=[];
            if(nev>(dnev-1))
                Nev = dnev-1;
                nev = nev-(dnev-1);
            end
            while(nev>dnev)
                Nev = cat(2,Nev,dnev);
                nev = nev-dnev;
            end
            Nev = cat(2,Nev,nev);

            X = double(D.t{iscope}.scope_trigger_ev_num_hist.x);
            if(~isempty(X))
                m = (x>=X(1))&(x<=X(numel(X)));
                y(m)=double(D.t{iscope}.scope_trigger_ev_num_hist.y);
            end
            
            [xs,ys]=zstairs(x,y./Nev);

            legendtext={ legendtext{:} tel_label{iscope} };
            plot(xs,ys,tel_plotspec{iscope});
            hold on
        end
    end
    if(~isempty(legendtext))
        zaxis([0 ceil(double(data.stage1.run_info.hi_event_num)/1000)*1000 ...
               -0.05 1.05]);
        xlabel(sprintf('Event number'));
        ylabel('Fraction of events triggered');
        title(t);
        zlegend(legendtext{:},4);
    else
        draw_no_plot(t)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW MISSING FRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_missing_event_number_hist(data)
    global analyze_revision
    if(analyze_revision < 2022)
        draw_not_available
        return
    end
    
    global l3_label l3_plotspec
    global tel_label tel_plotspec
    t='Number of DATUMS missing per 1000 event numbers';
    D=data.diagnostics;
    
    nev=double(data.stage1.run_info.hi_event_num);
    dnev = double(D.has_l3_ev_num_hist.bin_width);            

    tr_x = 0:dnev:nev;
    tr_y = tr_x*0;

    Nev=[];
    if(nev>(dnev-1))
        Nev = dnev-1;
        nev = nev-(dnev-1);
    end
    while(nev>dnev)
        Nev = cat(2,Nev,dnev);
        nev = nev-dnev;
    end
    Nev = cat(2,Nev,nev);

    X = double(D.has_l3_ev_num_hist.x);
    if(~isempty(X))
        m = (tr_x>=X(1))&(tr_x<=X(numel(X)));
        tr_y(m)=double(D.has_l3_ev_num_hist.y);
    end
    
    [xs,ys]=zstairs(tr_x,Nev-tr_y+0.001);
    
    semilogy(xs,ys,strcat(l3_plotspec,'-'));
    hold on
    legendtext={ l3_label };
    
    for iscope = 1:numel(D.t)
        if(isstruct(D.t{iscope}))
            Nev = double(D.t{iscope}.scope_sent_l3_ev_num_hist.y);
            x = double(D.t{iscope}.scope_sent_l3_ev_num_hist.x);
            y = Nev*0;
            X = double(D.t{iscope}.scope_has_event_ev_num_hist.x);
            if((~isempty(X))&&(~isempty(x)))
                xmlo=max([min(X) min(x)]);
                xmhi=min([max(X) max(x)]);
                mx = (x>=xmlo)&(x<=xmhi);
                mX = (X>=xmlo)&(X<=xmhi);
                y(mx)=double(D.t{iscope}.scope_has_event_ev_num_hist.y(mX));
            end
            [xs,ys]=zstairs(x,Nev-y+0.001);
            legendtext={ legendtext{:} tel_label{iscope} };
            semilogy(xs,ys,tel_plotspec{iscope});
            hold on
        end
    end
    if(~isempty(legendtext))
        zaxis([0 double(D.events_found) 0.5 Inf]);
        xlabel(sprintf('Event number'));
        ylabel('Number of missing DATUMS');
        title(t);
        zlegend(legendtext{:},4);
    else
        draw_no_plot(t)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW IMAGE CENTROIDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_centroid_hist2d(data,iscope)
    global analyze_revision
    t=sprintf('Centroid, camera coordinates [events] - T%1d',iscope);
    E=data.events;
    if((iscope<=numel(E.t))&&(isstruct(E.t{iscope})))
        nbins=20;
        rmax=1.6;
        if(analyze_revision < 3033)
            m=E.t{iscope}.nimage>0 & E.t{iscope}.fp_N>0;
        else
	    m=E.t{iscope}.has_image;
        end
        h2=simple_hist2([E.t{iscope}.fp_xc(m); E.t{iscope}.fp_yc(m)],...
                        [-rmax rmax], [-rmax rmax], nbins, nbins);
        h2ext=zeros(size(h2)+1);
        h2ext(1:size(h2,1),1:size(h2,2))=h2;
        x=(((0:nbins)+0.5)/nbins*2.0-1.0)*rmax;
        y=(((0:nbins)+0.5)/nbins*2.0-1.0)*rmax;
        [X Y]=meshgrid(x,y);
        R=sqrt(X.^2 + Y.^2);
        h2ext(R>rmax)=nan;
        x=((0:nbins)/nbins*2.0-1.0)*rmax;
        y=((0:nbins)/nbins*2.0-1.0)*rmax;
        pcolor(x,y,h2ext);
        shading('flat');
        axis('square');
        title(t);
        p=get(gca,'Position');
        p(1)=p(1)-0.1*p(3);
        set(gca,'Position',p);
        p(1)=p(1)+1.02*p(3);
        p(2)=p(2)+0.58*p(4);
        p(3)=p(3)*0.05;
        p(4)=p(4)*0.4;
        colorbar('Position',p);
    else
        draw_no_scope_camera(t,iscope);
    end
    xlabel('X [deg]')
    ylabel('Y [deg]')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW RECONSTRUCTED POSITION - SIZE CUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_reconstructed_hist2d(data,sizecut)
    global analyze_revision
    if(analyze_revision < 2000)
        draw_not_available
        return
    end
    E=data.events;
    m=E.nscope_image>=2;
    maxsize=zeros(size(m));
    for iscope=1:numel(E.t)
	if(isstruct(E.t{iscope}))
	        maxsize=max(maxsize,E.t{iscope}.fp_N);
	end
    end
    m=m&(maxsize>=sizecut);
    if(sizecut==0)
        t='Reconstructed direction in derotated camera [events]';
    else
        t=sprintf(strcat('Reconstructed direction in derotated',...
                         ' camera [events,size>%d]'),sizecut);
    end
    do_draw_reconstructed_hist2d(data,m,t)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW RECONSTRUCTED POSITION - SIMPLE SHAPE CUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_reconstructed_hist2d_simplecuts(data)
    global analyze_revision
    if(analyze_revision < 2000)
        draw_not_available
        return
    end
    E=data.events;
    %m=E.nscope_image>=2;

    pass_count = zeros(size(E.nscope_image));
    if((isfield(data,'cuts'))&&(isstruct(data.cuts))...
       &&(isfield(data.cuts,'npassed_scope_cut'))...
       &&(numel(data.cuts.npassed_scope_cut)==numel(pass_count)))
        pass_count = data.cuts.npassed_scope_cut;
        tt='[events,nspace cuts]';
    else
        for iscope=1:numel(E.t)
	    if(isstruct(E.t{iscope}))
                M = ones(size(E.t{iscope}.fp_N));
	        M = M & log10(E.t{iscope}.fp_N)>2.5; % 2.75
	        M = M & E.t{iscope}.delta2l<150;
	        M = M & E.t{iscope}.d1<-5000; % -4500
	        M = M & E.t{iscope}.d2>750; % 800
	        %sum(M)
	        pass_count(M) = pass_count(M)+1;
	    end
        end
        tt='[events,simple cuts]';
    end
    m=pass_count >= 2;
    t=strcat('Reconstructed direction in derotated camera ',tt);
    %sum(m)
    do_draw_reconstructed_hist2d(data,m,t)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMMON DRAW RECONSTRUCTED POSITION CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function do_draw_reconstructed_hist2d(data,m,t)
    global analyze_revision
    if(analyze_revision < 2000)
        draw_not_available
        return
    end
    E=data.events;
    nbins=20;
    rmax=1.6;
    fov_x=E.mean_derotated_fov_x(m);
    fov_y=E.mean_derotated_fov_y(m);
    h2=simple_hist2([fov_x; fov_y],...
                    [-rmax rmax], [-rmax rmax], nbins, nbins);
    h2ext=zeros(size(h2)+1);
    h2ext(1:size(h2,1),1:size(h2,2))=h2;
    x=(((0:nbins)+0.5)/nbins*2.0-1.0)*rmax;
    y=(((0:nbins)+0.5)/nbins*2.0-1.0)*rmax;
    [X Y]=meshgrid(x,y);
    R=sqrt(X.^2 + Y.^2);
    h2ext(R>rmax)=nan;
    x=((0:nbins)/nbins*2.0-1.0)*rmax;
    y=((0:nbins)/nbins*2.0-1.0)*rmax;
    pcolor(x,y,h2ext);
    shading('flat');
    axis('square');
    title(t);
    p=get(gca,'Position');
    p(1)=p(1)-0.1*p(3);
    set(gca,'Position',p);
    p(1)=p(1)+1.02*p(3);
    p(2)=p(2)+0.58*p(4);
    p(3)=p(3)*0.05;
    p(4)=p(4)*0.4;
    colorbar('Position',p);
    xlabel('X [deg]')
    ylabel('Y [deg]')
    a=axis;
    dx=a(2)-a(1);
    dy=a(4)-a(3);
    set(line([0.050 0.050]*dx+a(1),[0.025 0.100]*dy+a(3)),'Color','k');
    set(line([0.035 0.050 0.065]*dx+a(1),[0.085 0.100 0.085]*dy+a(3)),...
        'Color','k');
    text(0.050*dx+a(1),0.110*dy+a(3),'N','HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
    set(line([0.025 0.100]*dx+a(1),[0.050 0.050]*dy+a(3)),'Color','k');
    set(line([0.085 0.100 0.085]*dx+a(1),[0.035 0.050 0.065]*dy+a(3)),...
        'Color','k');
    text(0.110*dx+a(1),0.050*dy+a(3),'W','HorizontalAlignment','left',...
        'VerticalAlignment','middle');    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW RECONSTRUCTED CORE POSITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_reconstructed_core_hist2d(data)
    global analyze_revision
    if(analyze_revision < 2000)
        draw_not_available
        return
    end
    E=data.events;
    nbins=20;
    rmax=150;
    fov_x=E.mean_derotated_fov_x;
    fov_y=E.mean_derotated_fov_y;
    m=(E.nscope_image>=2) & (sqrt(fov_x.^2+fov_y.^2)<1.5);
    h2=simple_hist2([E.Rx(m); E.Ry(m)],...
                    [-rmax rmax], [-rmax rmax], nbins, nbins);
    h2ext=zeros(size(h2)+1);
    h2(h2<=0)=nan;
    h2ext(1:size(h2,1),1:size(h2,2))=log10(h2);
    x=(((0:nbins)+0.5)/nbins*2.0-1.0)*rmax;
    y=(((0:nbins)+0.5)/nbins*2.0-1.0)*rmax;
    [X Y]=meshgrid(x,y);
    R=sqrt(X.^2 + Y.^2);
    h2ext(R>rmax)=nan;
    x=((0:nbins)/nbins*2.0-1.0)*rmax;
    y=((0:nbins)/nbins*2.0-1.0)*rmax;
    pcolor(x,y,h2ext);
    shading('flat');
    axis('square');
    title('Reconstructed core position [log events]');
    p=get(gca,'Position');
    p(1)=p(1)-0.1*p(3);
    set(gca,'Position',p);
    p(1)=p(1)+1.02*p(3);
    p(2)=p(2)+0.58*p(4);
    p(3)=p(3)*0.05;
    p(4)=p(4)*0.4;
    colorbar('Position',p);
    xlabel('X [m]')
    ylabel('Y [m]')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW CAMERA FACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function camera = get_camera_500
    camera.x = [ ...  
 0.000,  0.148,  0.074, -0.074, -0.148, -0.074,  0.074,  0.296,  0.222, ...
 0.148,  0.000, -0.148, -0.222, -0.296, -0.222, -0.148,  0.000,  0.148, ...
 0.222,  0.444,  0.370,  0.296,  0.222,  0.074, -0.074, -0.222, -0.296, ...
-0.370, -0.444, -0.370, -0.296, -0.222, -0.074,  0.074,  0.222,  0.296, ...
 0.370,  0.592,  0.518,  0.444,  0.370,  0.296,  0.148,  0.000, -0.148, ...
-0.296, -0.370, -0.444, -0.518, -0.592, -0.518, -0.444, -0.370, -0.296, ...
-0.148,  0.000,  0.148,  0.296,  0.370,  0.444,  0.518,  0.740,  0.666, ...
 0.592,  0.518,  0.444,  0.370,  0.222,  0.074, -0.074, -0.222, -0.370, ...
-0.444, -0.518, -0.592, -0.666, -0.740, -0.666, -0.592, -0.518, -0.444, ...
-0.370, -0.222, -0.074,  0.074,  0.222,  0.370,  0.444,  0.518,  0.592, ...
 0.666,  0.888,  0.814,  0.740,  0.666,  0.592,  0.518,  0.444,  0.296, ...
 0.148,  0.000, -0.148, -0.296, -0.444, -0.518, -0.592, -0.666, -0.740, ...
-0.814, -0.888, -0.814, -0.740, -0.666, -0.592, -0.518, -0.444, -0.296, ...
-0.148,  0.000,  0.148,  0.296,  0.444,  0.518,  0.592,  0.666,  0.740, ...
 0.814,  1.036,  0.962,  0.888,  0.814,  0.740,  0.666,  0.592,  0.518, ...
 0.370,  0.222,  0.074, -0.074, -0.222, -0.370, -0.518, -0.592, -0.666, ...
-0.740, -0.814, -0.888, -0.962, -1.036, -0.962, -0.888, -0.814, -0.740, ...
-0.666, -0.592, -0.518, -0.370, -0.222, -0.074,  0.074,  0.222,  0.370, ...
 0.518,  0.592,  0.666,  0.740,  0.814,  0.888,  0.962,  1.184,  1.110, ...
 1.036,  0.962,  0.888,  0.814,  0.740,  0.666,  0.592,  0.444,  0.296, ...
 0.148,  0.000, -0.148, -0.296, -0.444, -0.592, -0.666, -0.740, -0.814, ...
-0.888, -0.962, -1.036, -1.110, -1.184, -1.110, -1.036, -0.962, -0.888, ...
-0.814, -0.740, -0.666, -0.592, -0.444, -0.296, -0.148,  0.000,  0.148, ...
 0.296,  0.444,  0.592,  0.666,  0.740,  0.814,  0.888,  0.962,  1.036, ...
 1.110,  1.332,  1.258,  1.184,  1.110,  1.036,  0.962,  0.888,  0.814, ...
 0.740,  0.666,  0.518,  0.370,  0.222,  0.074, -0.074, -0.222, -0.370, ...
-0.518, -0.666, -0.740, -0.814, -0.888, -0.962, -1.036, -1.110, -1.184, ...
-1.258, -1.332, -1.258, -1.184, -1.110, -1.036, -0.962, -0.888, -0.814, ...
-0.740, -0.666, -0.518, -0.370, -0.222, -0.074,  0.074,  0.222,  0.370, ...
 0.518,  0.666,  0.740,  0.814,  0.888,  0.962,  1.036,  1.110,  1.184, ...
 1.258,  1.480,  1.406,  1.332,  1.258,  1.184,  1.110,  1.036,  0.962, ...
 0.888,  0.814,  0.740,  0.592,  0.444,  0.296,  0.148,  0.000, -0.148, ...
-0.296, -0.444, -0.592, -0.740, -0.814, -0.888, -0.962, -1.036, -1.110, ...
-1.184, -1.258, -1.332, -1.406, -1.480, -1.406, -1.332, -1.258, -1.184, ...
-1.110, -1.036, -0.962, -0.888, -0.814, -0.740, -0.592, -0.444, -0.296, ...
-0.148,  0.000,  0.148,  0.296,  0.444,  0.592,  0.740,  0.814,  0.888, ...
 0.962,  1.036,  1.110,  1.184,  1.258,  1.332,  1.406,  1.628,  1.554, ...
 1.480,  1.406,  1.332,  1.258,  1.184,  1.110,  1.036,  0.962,  0.888, ...
 0.814,  0.666,  0.518,  0.370,  0.222,  0.074, -0.074, -0.222, -0.370, ...
-0.518, -0.666, -0.814, -0.888, -0.962, -1.036, -1.110, -1.184, -1.258, ...
-1.332, -1.406, -1.480, -1.554, -1.628, -1.554, -1.480, -1.406, -1.332, ...
-1.258, -1.184, -1.110, -1.036, -0.962, -0.888, -0.814, -0.666, -0.518, ...
-0.370, -0.222, -0.074,  0.074,  0.222,  0.370,  0.518,  0.666,  0.814, ...
 0.888,  0.962,  1.036,  1.110,  1.184,  1.258,  1.332,  1.406,  1.480, ...
 1.554,  1.702,  1.628,  1.554,  1.480,  1.406,  1.332,  1.258,  1.184, ...
 1.110,  1.036,  0.962,  0.740,  0.592,  0.444,  0.296,  0.148,  0.000, ...
-0.148, -0.296, -0.444, -0.592, -0.740, -0.962, -1.036, -1.110, -1.184, ...
-1.258, -1.332, -1.406, -1.480, -1.554, -1.628, -1.702, -1.702, -1.628, ...
-1.554, -1.480, -1.406, -1.332, -1.258, -1.184, -1.110, -1.036, -0.962, ...
-0.740, -0.592, -0.444, -0.296, -0.148,  0.000,  0.148,  0.296,  0.444, ...
 0.592,  0.740,  0.962,  1.036,  1.110,  1.184,  1.258,  1.332,  1.406, ...
 1.480,  1.554,  1.628,  1.702,  1.628,  1.554,  1.480,  1.406,  1.332, ...
 1.258,  0.370,  0.222,  0.074, -0.074, -0.222, -0.370, -1.258, -1.332, ...
-1.406, -1.480, -1.554, -1.628, -1.628, -1.554, -1.480, -1.406, -1.332, ...
-1.258, -0.370, -0.222, -0.074,  0.074,  0.222,  0.370,  1.258,  1.332, ...
 1.406,  1.480,  1.554,  1.628 ];
    camera.y = [ ...
 0.000,  0.000, -0.128, -0.128,  0.000,  0.128,  0.128,  0.000, -0.128, ...
-0.256, -0.256, -0.256, -0.128,  0.000,  0.128,  0.256,  0.256,  0.256, ...
 0.128,  0.000, -0.128, -0.256, -0.384, -0.384, -0.384, -0.384, -0.256, ...
-0.128,  0.000,  0.128,  0.256,  0.384,  0.384,  0.384,  0.384,  0.256, ...
 0.128,  0.000, -0.128, -0.256, -0.384, -0.513, -0.513, -0.513, -0.513, ...
-0.513, -0.384, -0.256, -0.128,  0.000,  0.128,  0.256,  0.384,  0.513, ...
 0.513,  0.513,  0.513,  0.513,  0.384,  0.256,  0.128,  0.000, -0.128, ...
-0.256, -0.384, -0.513, -0.641, -0.641, -0.641, -0.641, -0.641, -0.641, ...
-0.513, -0.384, -0.256, -0.128,  0.000,  0.128,  0.256,  0.384,  0.513, ...
 0.641,  0.641,  0.641,  0.641,  0.641,  0.641,  0.513,  0.384,  0.256, ...
 0.128,  0.000, -0.128, -0.256, -0.384, -0.513, -0.641, -0.769, -0.769, ...
-0.769, -0.769, -0.769, -0.769, -0.769, -0.641, -0.513, -0.384, -0.256, ...
-0.128,  0.000,  0.128,  0.256,  0.384,  0.513,  0.641,  0.769,  0.769, ...
 0.769,  0.769,  0.769,  0.769,  0.769,  0.641,  0.513,  0.384,  0.256, ...
 0.128,  0.000, -0.128, -0.256, -0.384, -0.513, -0.641, -0.769, -0.897, ...
-0.897, -0.897, -0.897, -0.897, -0.897, -0.897, -0.897, -0.769, -0.641, ...
-0.513, -0.384, -0.256, -0.128,  0.000,  0.128,  0.256,  0.384,  0.513, ...
 0.641,  0.769,  0.897,  0.897,  0.897,  0.897,  0.897,  0.897,  0.897, ...
 0.897,  0.769,  0.641,  0.513,  0.384,  0.256,  0.128,  0.000, -0.128, ...
-0.256, -0.384, -0.513, -0.641, -0.769, -0.897, -1.025, -1.025, -1.025, ...
-1.025, -1.025, -1.025, -1.025, -1.025, -1.025, -0.897, -0.769, -0.641, ...
-0.513, -0.384, -0.256, -0.128,  0.000,  0.128,  0.256,  0.384,  0.513, ...
 0.641,  0.769,  0.897,  1.025,  1.025,  1.025,  1.025,  1.025,  1.025, ...
 1.025,  1.025,  1.025,  0.897,  0.769,  0.641,  0.513,  0.384,  0.256, ...
 0.128,  0.000, -0.128, -0.256, -0.384, -0.513, -0.641, -0.769, -0.897, ...
-1.025, -1.153, -1.153, -1.153, -1.153, -1.153, -1.153, -1.153, -1.153, ...
-1.153, -1.153, -1.025, -0.897, -0.769, -0.641, -0.513, -0.384, -0.256, ...
-0.128,  0.000,  0.128,  0.256,  0.384,  0.513,  0.641,  0.769,  0.897, ...
 1.025,  1.153,  1.153,  1.153,  1.153,  1.153,  1.153,  1.153,  1.153, ...
 1.153,  1.153,  1.025,  0.897,  0.769,  0.641,  0.513,  0.384,  0.256, ...
 0.128,  0.000, -0.128, -0.256, -0.384, -0.513, -0.641, -0.769, -0.897, ...
-1.025, -1.153, -1.281, -1.281, -1.281, -1.281, -1.281, -1.281, -1.281, ...
-1.281, -1.281, -1.281, -1.281, -1.153, -1.025, -0.897, -0.769, -0.641, ...
-0.513, -0.384, -0.256, -0.128,  0.000,  0.128,  0.256,  0.384,  0.513, ...
 0.641,  0.769,  0.897,  1.025,  1.153,  1.281,  1.281,  1.281,  1.281, ...
 1.281,  1.281,  1.281,  1.281,  1.281,  1.281,  1.281,  1.153,  1.025, ...
 0.897,  0.769,  0.641,  0.513,  0.384,  0.256,  0.128,  0.000, -0.128, ...
-0.256, -0.384, -0.513, -0.641, -0.769, -0.897, -1.025, -1.153, -1.281, ...
-1.409, -1.409, -1.409, -1.409, -1.409, -1.409, -1.409, -1.409, -1.409, ...
-1.409, -1.409, -1.409, -1.281, -1.153, -1.025, -0.897, -0.769, -0.641, ...
-0.513, -0.384, -0.256, -0.128,  0.000,  0.128,  0.256,  0.384,  0.513, ...
 0.641,  0.769,  0.897,  1.025,  1.153,  1.281,  1.409,  1.409,  1.409, ...
 1.409,  1.409,  1.409,  1.409,  1.409,  1.409,  1.409,  1.409,  1.409, ...
 1.281,  1.153,  1.025,  0.897,  0.769,  0.641,  0.513,  0.384,  0.256, ...
 0.128, -0.128, -0.256, -0.384, -0.513, -0.641, -0.769, -0.897, -1.025, ...
-1.153, -1.281, -1.409, -1.538, -1.538, -1.538, -1.538, -1.538, -1.538, ...
-1.538, -1.538, -1.538, -1.538, -1.538, -1.409, -1.281, -1.153, -1.025, ...
-0.897, -0.769, -0.641, -0.513, -0.384, -0.256, -0.128,  0.128,  0.256, ...
 0.384,  0.513,  0.641,  0.769,  0.897,  1.025,  1.153,  1.281,  1.409, ...
 1.538,  1.538,  1.538,  1.538,  1.538,  1.538,  1.538,  1.538,  1.538, ...
 1.538,  1.538,  1.409,  1.281,  1.153,  1.025,  0.897,  0.769,  0.641, ...
 0.513,  0.384,  0.256,  0.128, -0.513, -0.641, -0.769, -0.897, -1.025, ...
-1.153, -1.666, -1.666, -1.666, -1.666, -1.666, -1.666, -1.153, -1.025, ...
-0.897, -0.769, -0.641, -0.513,  0.513,  0.641,  0.769,  0.897,  1.025, ...
 1.153,  1.666,  1.666,  1.666,  1.666,  1.666,  1.666,  1.153,  1.025, ...
 0.897,  0.769,  0.641,  0.513 ];
    camera.r = ones(size(camera.x))*0.074;
end

function draw_camera(vals,vlo,vhi,t,camera,fontsize,top10,gradient,varargin)
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

    if(nargin>8)
        units=varargin{1};
    else
        units='';
    end
    
    if(units == '1')
      units = '';
    end
    
    if(nargin>9)
      plotlog = 1;
    else
      plotlog = 0;
    end
    
    if(numel(vals)>numel(camera.x))
        vals=vals(1:numel(camera.x));
    end
    v=vals;


    
    if(plotlog == 1)
            
      v = log10(v);
      vpos = v(~isnan(v)&~isinf(v)&v>0);
      
      if(numel(vpos) == 0)
	vlo = -0.5;
	vhi = 0.5;
      else
      
	if(vlo > 0)
	  vlo = log10(vlo);
	elseif(vlo <= 0 && numel(vpos) > 0)
	  vlo = log10(min(vpos));
	else 
	  vlo = -0.5;
	end
	
	if(vhi > 0)
	  vhi = log10(vhi);
	elseif(vhi <= 0 && numel(vpos) > 0)
	  vhi = log10(max(vpos));
	else 
	  vhi = 0.5;
	end
	
	if(vlo >= vhi)
	  vhi = vlo+1;
	end
	
      end
    end
    

    
    index=1:numel(v); %find(~isnan(v));
    xchan=camera.x(index);
    ychan=camera.y(index);
    rchan=camera.r(index);
    if(isinf(vlo))
        vlo=min(v(~isnan(v)&~isinf(v)));
    end
    if(isinf(vhi))
        vhi=max(v(~isnan(v)&~isinf(v)));
    end
    v(v>vhi)=vhi;
    v(v<vlo)=vlo;
    cchan=v(index);
    
    circ=0;
    if(circ==1)
        xcirc=cos((0:5:360)/180*pi);
        ycirc=sin((0:5:360)/180*pi);
    else
        xcirc=cos((30:60:360)/180*pi)/cos(30/180*pi);
        ycirc=sin((30:60:360)/180*pi)/cos(30/180*pi);        
    end
    [xa,xb]=meshgrid(xchan,xcirc);
    [ya,yb]=meshgrid(ychan,ycirc);
    ra=meshgrid(rchan,zeros(size(ycirc)));
    x=xa+ra.*xb;
    y=ya+ra.*yb;
    caxis([vlo,vhi]);
    patch(x,y,cchan);
    axis('square');

    for ichan = 1:numel(index)
        text(xchan(ichan),ychan(ichan),sprintf('%d',index(ichan)),...
             'HorizontalAlignment','center',...
             'VerticalAlignment','middle',...
             'FontSize',fontsize);
    end
        
    title(t);
    p=get(gca,'Position');
    p(1)=p(1)-0.1*p(3);
    set(gca,'Position',p);
    p(1)=p(1)+1.02*p(3);
    p(2)=p(2)+0.58*p(4);
    p(3)=p(3)*0.05;
    p(4)=p(4)*0.4;
    colorbar('Position',p);

    if ( top10 == 1 )
        mgn=0.03;
        [sv si]=sort(vals);
        si=si(~isnan(sv));
        si=fliplr(flipud(si));
        if(numel(si)<12)
            si = [ si zeros(1,12-numel(si)) ];
        end
        txt = { 'Top 12'
                sprintf('%03d %03d %03d %03d',si(1:4))
                sprintf('%03d %03d %03d %03d',si(5:8))
                sprintf('%03d %03d %03d %03d',si(9:12)) };
        a=axis;
        set(text(a(1)+(a(2)-a(1))*mgn,a(3)+(a(4)-a(3))*mgn,txt),...
            'HorizontalAlignment','left','VerticalAlignment','bottom',...
            'FontName','fixed','FontSize',fontsize*1.25);
        si=fliplr(flipud(si));
        txt = { 'Bottom 12'
                sprintf('%03d %03d %03d %03d',si(1:4))
                sprintf('%03d %03d %03d %03d',si(5:8))
                sprintf('%03d %03d %03d %03d',si(9:12)) };
        a=axis;
        set(text(a(1)+(a(2)-a(1))*(1-mgn),a(3)+(a(4)-a(3))*mgn,txt),...
            'HorizontalAlignment','right','VerticalAlignment','bottom',...
            'FontName','fixed','FontSize',fontsize*1.25);
        if((~isempty(vals))&&(sum(~isnan(vals))>0))
	  mv = median(vals(~isnan(vals)));	  
	  txt=sprintf('Median: %g %s',mv,units);
	  
	  set(text(a(1)+(a(2)-a(1))*mgn,a(3)+(a(4)-a(3))*(1-mgn),txt),...
	      'HorizontalAlignment','left','VerticalAlignment','top',...
	      'FontName','fixed','FontSize',fontsize*1.25);
        end
    end
    
    if ( ( gradient == 1 ) && ( sum(~isnan(vals)) > 5 ) )
        sv=sort(vals(~isnan(vals)));
	sfact = 1.5;
        s80 = sv(floor(numel(sv)*0.8)+1);
        s20 = sv(floor(numel(sv)*0.2)+1);
        m=~isnan(vals) & vals<=s80*sfact & vals>=s20/sfact;

        % Linear gradient
        n=sum(m);
        x=sum(camera.x(m));
        y=sum(camera.y(m));
        z=sum(vals(m));
        xx=sum(camera.x(m).*camera.x(m));
        yy=sum(camera.y(m).*camera.y(m));
        xy=sum(camera.x(m).*camera.y(m));
        xz=sum(camera.x(m).*vals(m));
        yz=sum(camera.y(m).*vals(m));
        M=[n n n; x xx xy; y xy yy];
        V=[z; xz; yz];
        Minv=inv(M);
        ABC=Minv*V;
        B=ABC(2);
        C=ABC(3);
        
        % Radial gradient
        r=sqrt(camera.x.^2 + camera.y.^2);
        f=polyfit(r(m),vals(m),1);
        
        txt=sprintf('Gradient [%s/deg]\nLin: %.2g\nRad: %.2g',...
                    units,sqrt(B.^2+C^2),abs(f(1)));
        set(text(a(2)-(a(2)-a(1))*mgn,a(3)+(a(4)-a(3))*(1-mgn),txt),...
            'HorizontalAlignment','right','VerticalAlignment','top',...
            'FontName','fixed','FontSize',fontsize*1.25);

	ax = [ 5 -5 -5 5 ];
        ay = [ 0  2 -2 0 ];
	ct = cos(atan2(C,B));
	st = sin(atan2(C,B));
	set(line((ax*ct-ay*st)*(a(2)-a(1))*0.003+a(2)-(a(2)-a(1))*0.05,...
                 (ay*ct+ax*st)*(a(2)-a(1))*0.003+a(4)-(a(4)-a(3))*0.15),...
            'Color','k')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW NO DATA CAMERA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_no_scope_camera(t,iscope)
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
    axis('square');
    a=axis;
    text((a(1)+a(2))*0.5,(a(3)+a(4))*0.5,...
             sprintf('No data for T%1d available',iscope),...
             'HorizontalAlignment','center',...
             'VerticalAlignment','middle');
    title(t);
    p=get(gca,'Position');
    p(1)=p(1)-0.1*p(3);
    set(gca,'Position',p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PEDESTAL VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_pedvar_camera(data,iscope)
    t=sprintf('Pedestal RMS [DC] - T%1d',iscope);
    if((iscope<=numel(data.channel_info))...
       &&(isstruct(data.channel_info{iscope})))
        camera = get_camera_500;
        alldev=[];
        for jscope=1:numel(data.channel_info)
            if(isstruct(data.channel_info{jscope}))
                m=data.channel_info{jscope}.has_pmt;
                alldev=cat(2,alldev,...
                           data.channel_info{jscope}.median_dev(m));
            end
        end
        alldev=alldev(alldev>0);
        [xlo xhi]=find_limits(alldev,0.95,1.1);
        val=data.channel_info{iscope}.median_dev;
        m=data.channel_info{iscope}.has_pmt & val>0;
        val(~m)=nan;
        draw_camera(val,xlo,xhi,t,camera,3,1,1,'DC');
    else
        draw_no_scope_camera(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEDIAN CURRENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_median_current_camera(data,iscope)
    t=sprintf('Median anode current [\\muA] - T%1d',iscope);
    if((iscope<=numel(data.channel_info))...
       &&(isstruct(data.channel_info{iscope}))...
       &&(isfield(data.channel_info{iscope},'median_current')))
        camera = get_camera_500;
        val=data.channel_info{iscope}.median_current;
        [xlo xhi]=find_limits(val,0.97,1.2);
        draw_camera(val,xlo,xhi,t,camera,3,1,0,'uA');
    else
        draw_no_scope_camera(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_gain_camera(data,iscope)
    t=sprintf('Gains - T%1d',iscope);
    if((iscope<=numel(data.channel_info))...
       &&(isstruct(data.channel_info{iscope})))
        camera = get_camera_500;
        val=1.0./data.channel_info{iscope}.gain;
        m=data.channel_info{iscope}.has_pmt...
          & ~(data.channel_info{iscope}.suppress_reason&1);
        val(~m)=nan;
        draw_camera(val,0.85,1/0.85,t,camera,3,1,1,'1');
    else
        draw_no_scope_camera(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABSOLUTE GAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_absgain_camera(data,iscope)
    t=sprintf('Absolute gains [DC/PE] - T%1d',iscope);
    if((iscope<=numel(data.channel_info))...
       &&(isstruct(data.channel_info{iscope}))...
       &&(isfield(data.channel_info{iscope},'pmt_multiplier_gain')))
        camera = get_camera_500;
        allabsgain=[];
        for jscope=1:numel(data.channel_info)
            if(isstruct(data.channel_info{jscope}))
                m=data.channel_info{jscope}.has_pmt...
                  & ~(data.channel_info{jscope}.suppress_reason&1);
                val=data.channel_info{jscope}.pmt_multiplier_gain(m);
                allabsgain=cat(2,allabsgain,val);
            end
        end
        val=data.channel_info{iscope}.pmt_multiplier_gain;
        m=data.channel_info{iscope}.has_pmt...
          & ~(data.channel_info{iscope}.suppress_reason&1);
        val(~m)=nan;
        [xlo xhi]=find_limits(allabsgain,0.95,1.2);
        xlo = max([xlo 0]);
        draw_camera(val,xlo,xhi,t,camera,3,1,1,'DC/PE');
    else
        draw_no_scope_camera(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PMT EFFICIENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_collection_eff_camera(data,iscope)
    t=sprintf('Relative collection efficiency - T%1d',iscope);
    if((iscope<=numel(data.channel_info))...
       &&(isstruct(data.channel_info{iscope}))...
       &&(isfield(data.channel_info{iscope},'pmt_efficiency')))
        camera = get_camera_500;
        alleff=[];
        for jscope=1:numel(data.channel_info)
            if(isstruct(data.channel_info{jscope}))
                m=data.channel_info{jscope}.has_pmt...
                  & ~(data.channel_info{jscope}.suppress_reason&1);
                val=data.channel_info{jscope}.pmt_efficiency(m);
                alleff=cat(2,alleff,val);
            end
        end
        val=data.channel_info{iscope}.pmt_efficiency;
        m=data.channel_info{iscope}.has_pmt...
          & ~(data.channel_info{iscope}.suppress_reason&1);
        val(~m)=nan;
        [xlo xhi]=find_limits(alleff,0.95,1.2);
        xlo = max([xlo 0]);
        draw_camera(val,xlo,xhi,t,camera,3,1,0,'1');
    else
        draw_no_scope_camera(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LASER TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_laser_time_camera(data,iscope)
    t=sprintf(strcat('Mean LASER time with respect mean crate time ',...
                     '[ns] - T%1d'),iscope);
    if((iscope<=numel(data.channel_info))...
       &&(isstruct(data.channel_info{iscope})))
        camera = get_camera_500;
        allt0=[];
        for jscope=1:numel(data.channel_info)
            if(isstruct(data.channel_info{jscope}))
                m=data.channel_info{jscope}.has_pmt...
                  & ~(bitand(data.channel_info{iscope}.suppress_reason,1));
                val=(data.channel_info{jscope}.chantime(m)...
                     -data.channel_info{jscope}.cratetime(m))*2.0;
                allt0=cat(2,allt0,val);
            end
        end
        val=(data.channel_info{iscope}.chantime...
             -data.channel_info{iscope}.cratetime)*2.0;
        m=data.channel_info{iscope}.has_pmt...
          & ~(bitand(data.channel_info{iscope}.suppress_reason,1));
        val(~m)=nan;
        [xlo xhi]=find_limits(allt0,0.95,1.2);
        draw_camera(val,xlo,xhi,t,camera,3,1,0,'ns');
    else
        draw_no_scope_camera(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CFD THRESHOLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_cfd_threshold_camera(data,iscope)
    t=sprintf('CFD trigger threshold [mV] - T%1d',iscope);
    D=data.diagnostics;
    CI=data.channel_info;
    if((iscope<=numel(D.t))&&(isstruct(D.t{iscope}))...
       &&(iscope<=numel(CI))&&(isstruct(CI{iscope})))
        camera = get_camera_500;
        allcfd=[];
        for jscope=1:numel(D.t)
            if((isstruct(D.t{jscope}))...
               &&(jscope<=numel(CI))&&(isstruct(CI{jscope})))
                ntr=D.t{iscope}.channel_is_triggered;
                mtr=ntr>10;
                m=CI{jscope}.has_pmt & mtr;
                allcfd=cat(2,allcfd,D.t{jscope}.channel_cfd_threshold(m));
            end
        end
        cfd=D.t{iscope}.channel_cfd_threshold;
        ntr=D.t{iscope}.channel_is_triggered;
        [xlo xhi]=find_limits(allcfd,0.90,1.2);
        cfd(ntr<=10)=nan;
        m=CI{iscope}.has_pmt;
        cfd(~m)=nan;
        draw_camera(cfd,xlo,xhi,t,camera,3,1,0,'mV');
    else
        draw_no_scope_camera(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CFD PARTICIPATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_cfd_participation_camera(data,iscope)
    t=sprintf('CFD trigger participation [events] - T%1d',iscope);
    D=data.diagnostics;
    if((iscope<=numel(D.t))&&(isstruct(D.t{iscope})))
        camera = get_camera_500;
        val=D.t{iscope}.channel_is_triggered;
        [xlo xhi]=find_limits(val,0.97,1.2);
        xlo = 0;
        draw_camera(val,xlo,xhi,t,camera,3,1,0,'events');
    else
        draw_no_scope_camera(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CFD L1 RATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_cfd_rate_camera(data,iscope)
    t=sprintf('Median CFD trigger rate [log_{10}(Hz)] - T%1d',iscope);
    if((iscope<=numel(data.channel_info))...
       &&(isstruct(data.channel_info{iscope}))...
       &&(isfield(data.channel_info{iscope},'median_l1_rate')))
        camera = get_camera_500;
        val=data.channel_info{iscope}.median_l1_rate;
        [xlo xhi]=find_limits(val(val>100),0.97,1.0);
		
        draw_camera(val,xlo,xhi,t,camera,3,1,0,'Hz',1);
    else
        draw_no_scope_camera(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ISOLATED CFD FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_isolated_frequency_camera(data,iscope)
    t=sprintf('Isolated CFD trigger frequency [events] - T%1d',iscope);
    D=data.diagnostics;
    if((iscope<=numel(D.t))&&(isstruct(D.t{iscope})))
        camera = get_camera_500;
        val=D.t{iscope}.channel_is_triggered_isolated;
        [xlo xhi]=find_limits(val,0.97,1.2);
        xlo = 0;
        draw_camera(val,xlo,xhi,t,camera,3,1,0,'events');
    else
        draw_no_scope_camera(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUB THRESHOLD CFD FREQUENCY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_subthreshold_frequency_camera(data,iscope)
    global analyze_revision
    if(analyze_revision < 1052)
        draw_not_available_square
        return
    end

    t=sprintf('CFD trigger frequency in sub-threshold events [events] - T%1d',iscope);
    D=data.diagnostics;
    if((iscope<=numel(D.t))&&(isstruct(D.t{iscope})))
        camera = get_camera_500;
        val=D.t{iscope}.channel_is_triggered_in_sub_threshold_event;
        [xlo xhi]=find_limits(val,0.97,1.2);
        xlo = 0;
        draw_camera(val,xlo,xhi,t,camera,3,1,0,'events');
    else
        draw_no_scope_camera(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOP-3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_raw_top3_camera(data,iscope)
    t=sprintf('Raw TOP-3 frequency [events] - T%1d',iscope);
    D=data.diagnostics;
    CI=data.channel_info;
    if((iscope<=numel(D.t))&&(isstruct(D.t{iscope}))...
       &&(iscope<=numel(CI))&&(isstruct(CI{iscope})))
        camera = get_camera_500;
        val=D.t{iscope}.channel_is_raw_top3;
        m=CI{iscope}.has_pmt;
        val(~m)=nan;
        [xlo xhi]=find_limits(val,0.97,1.2);
        draw_camera(val,xlo,xhi,t,camera,3,1,0,'events');
    else
        draw_no_scope_camera(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAX-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_raw_max1_camera(data,iscope)
    t=sprintf('Raw MAX-1 frequency [events] - T%1d',iscope);
    D=data.diagnostics;
    CI=data.channel_info;
    if((iscope<=numel(D.t))&&(isstruct(D.t{iscope}))...
       &&(iscope<=numel(CI))&&(isstruct(CI{iscope})))
        camera = get_camera_500;
        val=D.t{iscope}.channel_is_raw_max1;
        m=CI{iscope}.has_pmt;
        val(~m)=nan;
        [xlo xhi]=find_limits(val,0.97,1.2);
        draw_camera(val,xlo,xhi,t,camera,3,1,0,'events');
    else
        draw_no_scope_camera(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LO GAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_lo_gain_camera(data,iscope)
    t=sprintf('LOW Gain frequency [events] - T%1d',iscope);
    D=data.diagnostics;
    CI=data.channel_info;
    if((iscope<=numel(D.t))&&(isstruct(D.t{iscope}))...
       &&(iscope<=numel(CI))&&(isstruct(CI{iscope})))
        camera = get_camera_500;
        val=D.t{iscope}.channel_is_lo_gain;
        m=CI{iscope}.has_pmt;
        val(~m)=nan;
        [xlo xhi]=find_limits(val,0.97,1.2);
        draw_camera(val,xlo,xhi,t,camera,3,1,0,'events');
    else
        draw_no_scope_camera(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XY PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_xy(x,vals,units,style)
    sv=sort(vals(~isnan(vals)));
    sfact = 1.5;
    s80 = sv(floor(numel(sv)*0.8)+1);
    s20 = sv(floor(numel(sv)*0.2)+1);
    m=~isnan(vals) & vals<=s80*sfact & vals>=s20/sfact;
    plot(x(m),vals(m),style)
    P1=polyfit(x(m),vals(m),1);
    P2=polyfit(x(m),vals(m),2);
    a=axis;
    mgn=0.03;
    txt=sprintf(['Lin coeff: %.2g [%s/deg]\n',...
                 'Quad coeff: %.2g [%s/deg]'],...
                P1(1),units,P2(1),units);
    set(text(a(1)+(a(2)-a(1))*mgn,a(3)+(a(4)-a(3))*(1-mgn),txt),...
        'Units','Normalized',...
        'HorizontalAlignment','left','VerticalAlignment','top',...
        'FontName','fixed'); %,'FontSize',fontsize*1.25);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAIN XY PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_gain_xy(data,iscope,ixy)
    t=sprintf('Gains - T%1d',iscope);
    if((iscope<=numel(data.channel_info))...
       &&(isstruct(data.channel_info{iscope})))
        camera = get_camera_500;
	style={'bs','bx'};
	xlab={'Tube X-coordinate [deg]','Tube Y-coordinate [deg]'};
	XY={camera.x,camera.y};
        val=1.0./data.channel_info{iscope}.gain;
        m=data.channel_info{iscope}.has_pmt...
          & ~(data.channel_info{iscope}.suppress_reason&1);
        val(~m)=nan;
        draw_xy(XY{ixy},val,'1',style{ixy});
	xlabel(xlab{ixy})
	ylabel('Gain [1]');
	title(t);
    else
        draw_no_scope(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABSGAIN XY PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_absgain_xy(data,iscope,ixy)
    t=sprintf('Absolute Gain - T%1d',iscope);
    if((iscope<=numel(data.channel_info))...
       &&(isstruct(data.channel_info{iscope})))
        camera = get_camera_500;
	style={'bs','bx'};
	xlab={'Tube X-coordinate [deg]','Tube Y-coordinate [deg]'};
	XY={camera.x,camera.y};
        val=data.channel_info{iscope}.pmt_multiplier_gain;
        m=data.channel_info{iscope}.has_pmt...
          & ~(data.channel_info{iscope}.suppress_reason&1);
        val(~m)=nan;
        draw_xy(XY{ixy},val,'DC/PE',style{ixy});
	xlabel(xlab{ixy})
	ylabel('Absolute Gain [DC/PE]');
	title(t);
    else
        draw_no_scope(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PEDVAR XY PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_pedvar_xy(data,iscope,ixy)
    t=sprintf('Pedestal RMS - T%1d',iscope);
    if((iscope<=numel(data.channel_info))...
       &&(isstruct(data.channel_info{iscope})))
        camera = get_camera_500;
	style={'bs','bx'};
	xlab={'Tube X-coordinate [deg]','Tube Y-coordinate [deg]'};
	XY={camera.x,camera.y};
        val=data.channel_info{iscope}.median_dev;
        m=data.channel_info{iscope}.has_pmt...
          & ~(data.channel_info{iscope}.suppress_reason&1);
        val(~m)=nan;
        draw_xy(XY{ixy},val,'DC',style{ixy});
	xlabel(xlab{ixy})
	ylabel('Pedestal RMS [DC]');
	title(t);
    else
        draw_no_scope(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUPPRESSED NSLICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_suppressed_nslice(data,iscope)
    t=sprintf('Channels suppressed - T%1d',iscope);
    if((iscope<=numel(data.channel_info))...
       &&(isstruct(data.channel_info{iscope})))
        nslice=0;
        for jscope=1:numel(data.channel_info)
            if(isstruct(data.channel_info{jscope}))
                val=double(data.channel_info{jscope}.suppress_nslice);
                nslice=max([nslice max(val)]);
            end
        end
        camera = get_camera_500;
        val=double(data.channel_info{iscope}.suppress_nslice);
        if(nslice>1)
            t=sprintf('Number of slices channel suppressed - T%1d',iscope);
        end
        val(val==0)=nan;
        draw_camera(val,0,nslice,t,camera,3,0,0,'slice');
    else
        draw_no_scope_camera(t,iscope);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT MUON PAREMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nmuon, u0_mean, u0_dev] = muon_stats(muons)
    edge_dist=sqrt(muons.x0.^2+muons.y0.^2)+muons.r0;
    u0=muons.r_U0_r0_corr./(muons.r0/pi*180);
    mask=(edge_dist<1.72/180*pi)...
         &(muons.r_rms<0.07/180*pi)...
         &(muons.r_xi<0.9)...
         &(u0>2e3)...
         &(u0<7e3);

    nmuon = sum(mask);
    if(nmuon>0)
        u0_mean = mean(u0(mask));
        u0_dev = std(u0(mask));
    else
        u0_mean = 0;
        u0_dev = 0;
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE FINAL RUN INFO ELEMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function finish_run_info(data)
    global run_info;
    D=data.diagnostics;
    C=data.channel_info;
    
    itrl3 = double(D.scope_ntrigger.x);
    ntrl3 = double(D.scope_ntrigger.y);
    npeak = max(ntrl3);
    ipeak = find(ntrl3==npeak,1,'first');
    ithrl3 = find(ntrl3>npeak*0.5);
    if(~isempty(ithrl3))
        ithrl3 = min(ithrl3(ithrl3<=ipeak));
    else
        ithrl3 = ipeak;
    end
    run_info.l3_threshold=(itrl3(ithrl3));

    run_info.laser = 0;
    if(isfield(data,'stage1') && isfield(data.stage1,'laser'))
        run_info.laser = data.stage1.laser.m_runno;
    end
    
    for iscope=1:numel(D.t)
        if(~isstruct(D.t{iscope}))
	  run_info.t{iscope}.laser_runno = 0;
	  run_info.t{iscope}.median_pedvar = 0;
	  run_info.t{iscope}.median_pedvar_pe = 0;
	  run_info.t{iscope}.fir_mean = 0;
	  run_info.t{iscope}.fir_dev = 0;
	end
    end
	  
    for iscope=1:numel(D.t)
        if(isstruct(D.t{iscope}))
            Dt=D.t{iscope};
            itlr = double(Dt.camera_ntrigger_largest_region.x);
            ntlr = double(Dt.camera_ntrigger_largest_region.y);
            npeak = max(ntlr);
            ipeak = find(ntlr==npeak,1,'first');
            ithr = find(ntlr>npeak*0.5);
            if(~isempty(ithr))
                ithr = min(ithr(ithr<=ipeak));
            else
                ithr = ipeak;
            end
            if(~isempty(ithr))
                run_info.t{iscope}.l2_threshold = itlr(ithr);
            else
                run_info.t{iscope}.l2_threshold = 0;
            end
            run_info.t{iscope}.nmissing_events = ...
                double(Dt.scope_sent_l3)-double(Dt.scope_has_event);
        else
            run_info.t{iscope}.l2_threshold = 0;
            run_info.t{iscope}.nmissing_events = 0;
        end
    end

    for iscope=1:numel(D.t);
        if((iscope<=numel(C))&&(isstruct(C{iscope})))
            mpmt=C{iscope}.has_pmt;

            if(isfield(C{iscope},'pmt_multiplier_gain'))
                val=C{iscope}.pmt_multiplier_gain;
                m = mpmt & ~(C{iscope}.suppress_reason&1);
                val(~m)=nan;
                run_info.t{iscope}.median_absgain = zmedian(val);
            else
                run_info.t{iscope}.median_absgain = 0;
            end

            if(isfield(C{iscope},'median_current'))
                val=C{iscope}.median_current;
                val(~mpmt)=nan;
                run_info.t{iscope}.median_current = zmedian(val);
            else
                run_info.t{iscope}.median_current = 0;
            end

            if(isfield(C{iscope},'median_l1_rate'))
                val=C{iscope}.median_l1_rate;
                val(~mpmt)=nan;
                run_info.t{iscope}.median_l1_rate = zmedian(val);
            else
                run_info.t{iscope}.median_l1_rate = 0;
            end
        else
            run_info.t{iscope}.median_absgain = 0;
            run_info.t{iscope}.median_current = 0;
            run_info.t{iscope}.median_l1_rate = 0;          
        end
    end
    
    if(isfield(data,'muon_analysis'))
        for iscope=1:numel(data.muon_analysis.scope)
            if(isstruct(data.muon_analysis.scope{iscope}))
                [ run_info.t{iscope}.muon_n,...
                  run_info.t{iscope}.muon_u0r0_mean,...
                  run_info.t{iscope}.muon_u0r0_dev ] = ...
                  muon_stats(data.muon_analysis.scope{iscope});
            else
                run_info.t{iscope}.muon_n = 0;
                run_info.t{iscope}.muon_u0r0_mean = 0;
                run_info.t{iscope}.muon_u0r0_dev = 0;
            end
        end
    else
        for iscope=1:numel(run_info.t)
            if(isstruct(run_info.t{iscope}))
                run_info.t{iscope}.muon_n = 0;
                run_info.t{iscope}.muon_u0r0_mean = 0;
                run_info.t{iscope}.muon_u0r0_dev = 0;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE RUN INFO ELEMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_run_info(file,info)
    target = info.target;
    target(target==' ') = '_';
    fp=fopen(file,'a');
    fprintf(fp,...
            strcat('%-5d 2 %s %11.5f %-20s %-13s %5.2f %+7.2f',...
                   ' %5.2f %5.2f %6.2f %6.2f %6.2f %5.2f %1d'),...
            info.runno,...
            info.run_start_time_string,...
            info.run_start_time,...
            target,...
            info.mode,...
            info.mean_el,...
            info.mean_az,...
            info.elapsed/60,...
            info.livetime/60,...
            info.rate,...
            info.rate_per_min_mean,...
            info.rate_per_min_std,...
            sqrt(info.rate_constant_chi2dof),...
            info.l3_threshold);

    for iscope=1:4
        if((iscope<=numel(info.t))&&(isstruct(info.t{iscope})))
	  
	  info.t{iscope}
	  
            fprintf(fp,...
                    ' %3d %7.2f %5.2f %5.3f %4.2f %5.3f',...
                    info.t{iscope}.l2_threshold,...
                    info.t{iscope}.cfd_threshold,...
                    info.t{iscope}.cfd_threshold_err,...
                    info.t{iscope}.spect_F,...
                    info.t{iscope}.spect_gamma,...
                    info.t{iscope}.spect_throughput);
        else
            fprintf(fp,...
                    ' %3d %7.2f %5.2f %5.3f %4.2f %5.3f',...
                    0,0,0,0,0,0);
        end
    end

    fprintf(fp,...
            strcat(' %5.2f %5.2f %-5d %+6.2f %3d'),...
            info.gpseltime/60,...
            info.delta_t_live/60,...
            info.laser,...
            info.moon_el,...
            info.moon_phase);

    for iscope=1:4
        if((iscope<=numel(info.t))&&(isstruct(info.t{iscope})))
	  	  	  
            fprintf(fp,...
                    strcat(' %+8.3f %6.3f %5.3f %5.3f %6.1f %4.1f',...
                           ' %7.1f %-6d %-3d %6.1f %6.1f %5.2f %5.1f',...
                           ' %3d'),...
                    info.mean_ra_vec(iscope),...
                    info.mean_dec_vec(iscope),...
                    info.dev_ra_vec(iscope),...
                    info.dev_dec_vec(iscope),...
                    info.t{iscope}.median_l1_rate/1000,...
                    info.t{iscope}.median_current,...
                    info.t{iscope}.l2_rate_per_min_median,...
                    info.t{iscope}.nmissing_events,...
                    info.t{iscope}.muon_n,...
                    info.t{iscope}.muon_u0r0_mean,...
                    info.t{iscope}.muon_u0r0_dev,...
                    info.t{iscope}.median_absgain,...
                    info.mean_moon_sep_vec(iscope),...
                    info.t{iscope}.nsuppressed);
        else
            fprintf(fp,...
                    strcat(' %+8.3f %6.3f %5.3f %5.3f %6.1f %4.1f',...
                           ' %7.1f %-6d %-3d %6.1f %6.1f %5.2f %5.1f',...
                           ' %3d'),...
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0);
        end
    end

    fprintf(fp,...
            strcat(' %3d %3d %6.3f %6.3f %6.3f %6.2f '),...
	    info.nscope,...
	    info.scope_mask,...
	    info.mean_pedvar,...
	    info.mean_pedvar_pe,...
            info.fir0_mean,...
            info.fir0_dev);
    
    for iscope=1:4
        if((iscope<=numel(info.t))&&(isstruct(info.t{iscope})))
	  iscope
	  
            fprintf(fp,...
                    strcat(' %-5d %6.3f %6.3f %6.3f %6.2f'),...
		    info.t{iscope}.laser_runno,...
		    info.t{iscope}.median_pedvar,...
		    info.t{iscope}.median_pedvar_pe,...
		    info.t{iscope}.fir_mean,...
		    info.t{iscope}.fir_dev);
        else
            fprintf(fp,...
                    strcat(' %-5d %6.3f %6.3f %6.3f %6.2f'),0,0,0,0,0);
        end
    end
    
    
    fprintf(fp,'\n');
    fclose(fp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

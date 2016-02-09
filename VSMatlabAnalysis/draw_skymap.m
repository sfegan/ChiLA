function draw_skymap(file_or_data,hist,tlab,varargin)
  
  if(ischar(file_or_data) && ischar(hist))
    data = loadh5(file_or_data);
    h = loadh5(file_or_data,hist);
  else
    data = file_or_data;
    h = hist;
  end
  
  x = h.x;
  y = h.y;
  z = h.z;
    
  if(numel(z) == 0)
    return;
  end
    
  zmin = min(min(z));
  zmax = max(max(z));  
    
  if(numel(varargin) == 1)
    
    zmin = inf;
    zmax = -inf;
    
    for i=1:numel(z)
      
      if(z(i) > 0)
	z(i) = log10(z(i));
	zmin = min( [z(i) zmin] );
	zmax = max( [z(i) zmax] );
      end
      
    end
    
    if(zmin == inf || zmax == -inf)
      return
    end

  end
    
    zexp = data.results.sky_exposure_hist.z;
    
    colormap([1 1 1; jet(128)])
            
    if(zmin == zmax)
      return
    end
    
    zmin = zmin - (zmax-zmin)/127.;    
    [nx ny] = size(z);
    for i=1:nx
      for j=1:ny
	
	if(zexp(i,j) == 0)
	  z(i,j) = -999;
	end
	  
      end
    end
        
    [x1 y1] = radec_to_xy(data.stage3.origin_ra_rad, ...
			  pi/2.-data.stage3.origin_dec_rad,...
			  data.stage3.src_ra_rad,...
			  pi/2.-data.stage3.src_dec_rad);
    
    
    z2 = z';
    clims = [zmin,zmax];

    imagesc(x,y,z2,clims)
    hold on;
    axis xy;
    %    colorbar;
    p=get(gca,'Position');
    p
    p(1)=p(1)-0.1*p(3);
%    set(gca,'Position',p);
    p(1)=p(1)+1.02*p(3);
    p(2)=p(2)+0.58*p(4);
    p(3)=p(3)*0.05;
    p(4)=p(4)*0.4;
%    p(1)=p(1)+1.02*p(3);
%    p(2)=0.1;
%    p(3)=p(1)+0.1;
%    p(4)=0.9;
%    p
%    colorbar('Position',p);
    colorbar('location','EastOutside');
    
%    colorbar('location','east');
    axis('image');
    
    draw_radec_grid(data);
    
    plot(x1,y1,'+k','LineWidth',0.2,'MarkerSize',5);
    
    draw_exclusion_regions(data);
      
    xlabel('{\Delta}RA [deg]')
    ylabel('{\Delta}Dec [deg]')
    title(tlab,'FontWeight','bold','FontSize',14)    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ELAPSED TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x y] = radec_to_xy(origin_ra, origin_dec, ra, dec)
  
  ra = ra - origin_ra;
  
  sth1 = sin(dec);
  x1 = sth1*cos(ra);
  y1 = sth1*sin(ra);
  z1 = cos(dec);
  
  sth2 = sin(-origin_dec);
  cth2 = cos(-origin_dec);
  
  xx = x1*cth2+z1*sth2;
  zz = z1*cth2-x1*sth2;
  
  th = atan2(sqrt(y1*y1+xx*xx),zz);
  phi = atan2(y1,xx);
  
  x = -th*180/pi*sin(phi);
  y = -th*180/pi*cos(phi);
  
end

function draw_radec_grid(data)
    
  ora = data.stage3.origin_ra_rad;
  odec = data.stage3.origin_dec_rad;

  histxlo = 1.1*data.stage3.sky_counts_hist.xlo;
  histxhi = 1.1*data.stage3.sky_counts_hist.xhi;
  histylo = 1.1*data.stage3.sky_counts_hist.ylo;
  histyhi = 1.1*data.stage3.sky_counts_hist.yhi;
  
  
  ra_step = (1.0/cos(odec))*pi/180.;  
%  ra_step = 0.001454441*floor(ra_step/0.001454441);
  ra_step = 0.00436332*floor(ra_step/0.00436332);
  
  ralo = ra_step*floor((ora - 3*ra_step)/ra_step);
  rahi = ra_step*floor((ora + 3*ra_step)/ra_step);
  
      
  src_dec = data.stage3.src_dec_rad;
  
  
  cx = xlim;
  cy = ylim;
      
  for ra = ralo:ra_step:rahi
      
    decx = [];
    decy = [];
    
    xlab = 0;
    ylab = cy(1) + 0.04*(cy(2)-cy(1));
    
    declo = odec+histylo*pi/180.;
    dechi = odec+histyhi*pi/180.;
    
    for dec = declo:0.001:dechi
      
      [x y] = radec_to_xy(ora, pi/2.-odec, ra, pi/2.-dec);
    
      if(y > ylab && xlab == 0)
	xlab = x;
      end
      
      decx = [decx x];
      decy = [decy y];
      
    end
  
    dx = decx(numel(decx)) - decx(1);
    dy = max(decy)-min(decy);
    rot_angle = 180+atan(dy/dx)*180./pi;
    
    if(dx > 0)
      rot_angle = rot_angle - 180;
    end
    
        
    xlab = xlab - 0.1;
    
    if(xlab > cx(1) + 0.1 && xlab < cx(2) - 0.1)
      
      label = rad2hmsstring(ra);
      
      label
      
      label = label(1:5);      
      
      text(xlab,ylab,sprintf('%s',label),...
	   'Rotation',rot_angle,'FontSize',4);
    end
      
    plot(decx,decy,'--k','LineWidth',0.1);
    hold on;
  end
  
  dec_step = 1.0*pi/180.;  
  dec_step = 0.001454441*floor(dec_step/0.001454441);
  
  dec_step*180/pi
  
  declo = dec_step*floor((odec - 3*dec_step)/dec_step);
  dechi = dec_step*floor((odec + 3*dec_step)/dec_step);
    
  for dec = declo:dec_step:dechi
      
    rax = [];
    ray = [];
    
    xlab = cx(1) + 0.04*(cx(2)-cx(1));
    ylab = 0;      
    dydx = 0;
    

    
    ralo = ora-histxhi/cos(odec)*pi/180.;
    rahi = ora-histxlo/cos(odec)*pi/180.;
    
    for ra = ralo:0.001:rahi
      
      [x y] = radec_to_xy(ora, pi/2.-odec, ra, pi/2.-dec);  
       
      if(x < xlab && ylab == 0 && numel(ray) > 0)
	ylab = y;
	dydx = (y - ray(numel(ray)))/(x - rax(numel(rax)));
      end
      
      rax = [rax x];
      ray = [ray y];
      
    end
  
    rot_angle = atan(dydx)*180./pi;
    ylab = ylab + 0.1;
    
    if(ylab > cy(1) + 0.1 && ylab < cy(2) - 0.1)
      
      label = rad2dmsstring(dec);
      
      label
      
      label = label(1:6);      
      text(xlab,ylab,sprintf('%s',label),...
	   'Rotation',rot_angle,'FontSize',4);
    end
      
    plot(rax,ray,'--k','LineWidth',0.1);
    hold on;
  end
  
end

function draw_exclusion_regions(data)
  
  for i=1:numel(data.exclusion_regions)
    
    x = data.exclusion_regions{i}.xy(1);
    y = data.exclusion_regions{i}.xy(2);    
    
         
    if(strcmp(data.exclusion_regions{i}.type,'star'))
    
      circle(x,y,data.exclusion_regions{i}.radius_deg);      
      
      marker_size = 1.0;
      
      if(data.exclusion_regions{i}.vmag < 6)
      
	marker_size = marker_size + (6 - data.exclusion_regions{i}.vmag)/3.;
      
      end
      
      plot(x,y,'ko','LineWidth',0.1,'MarkerSize',marker_size,...
	   'MarkerFaceColor','k');
      
    end
    
  end
  
end

function circle(x0, y0, r)
  
  t = linspace(0,2*pi,50);

  x = r*cos(t)+x0;                  
  y = r*sin(t)+y0;                  
  plot(x,y,'k','LineWidth',0.5);           

  hold on;
  
end

function load_contours(file,origin_ra_rad,origin_dec_rad)
  
  fid = fopen(file);
  
  C = textscan(fid,'%f %f %f %f %f','commentStyle', '#');
  
  fclose(fid);
  
  ra = C{1,1};
  dec = C{1,2};
  r = C{1,3};
  e = C{1,4};
  phi = C{1,5};
  
  nx = size(ra);
  
  for i=1:nx
    
    ra(i) = ra(i)*pi/180.;
    dec(i) = dec(i)*pi/180.;
    
    [x0 y0] = radec_to_xy(origin_ra_rad,pi/2.-origin_dec_rad,...
			  ra(i),pi/2.-dec(i));
        
    a = r(i)/60.;
    b = r(i)*(1-e(i))/60.;
    
    ellipse(x0,y0,a,b,phi(i)*pi/180.+pi/2.);

  end
  
end

function ellipse(x0, y0, a, b, phi)


  
  t = linspace(0,2*pi,100);

  x = x0+a*cos(t)*cos(phi)-b*sin(t)*sin(phi);                  
  y = y0+a*cos(t)*sin(phi)+b*sin(t)*cos(phi);
  plot(x,y,'-k','LineWidth',1.0);           

  hold on;
  
end
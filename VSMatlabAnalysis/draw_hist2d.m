%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAW 2D HIST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_hist2d(file_or_data,hist,tlab,varargin)
  
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
  
  colormap([1 1 1; jet(128)])
  
  zt = z';
  
  if(numel(varargin) == 2)    
    
    dx = varargin{1};
    dy = varargin{2};
    
    ztmp = zeros(int32(floor(numel(y)/dy)),int32(floor(numel(x)/dx)));
    
    [rows cols] = size(ztmp);
    
    x2 = zeros(cols,1);
    y2 = zeros(rows,1);
    
    for ix1=1:numel(x2)      
      nx = min([(ix1-1)*dx+dx numel(x)]);
      for ix2=(ix1-1)*dx+1:nx    
	x2(ix1) = x2(ix1) + x(ix2)/dx;
      end      
    end
    
    for iy1=1:numel(y2) 
      ny = min([(iy1-1)*dy+dy numel(y)]);
      for iy2=(iy1-1)*dy+1:ny    
	y2(iy1) = y2(iy1) + y(iy2)/dy;
      end      
    end
    
    for ix1=1:numel(x2) 
      for iy1=1:numel(y2)   
	
	nx = min([(ix1-1)*dx+dx numel(x)]);
	ny = min([(iy1-1)*dy+dy numel(y)]);
	
	for ix2=(ix1-1)*dx+1:nx      
	  for iy2=(iy1-1)*dy+1:ny      	    	    
	    ztmp(iy1,ix1) = ztmp(iy1,ix1) + zt(iy2,ix2);
	  end
	end
      end
    end
    
    zt = ztmp;
    x = x2;
    y = y2;
  end
    
    
  imagesc(x,y,zt)
%  surfc(y,x,zt)
  axis xy;
  
  p=get(gca,'Position');
  p(1)=p(1)-0.1*p(3);
  set(gca,'Position',p);
  p(1)=p(1)+1.02*p(3);
  p(2)=p(2)+0.58*p(4);
  p(3)=p(3)*0.05;
  p(4)=p(4)*0.4;
  colorbar('Position',p);
  
  axis('image');
  
  xlabel('X [deg]')
  ylabel('Y [deg]')
  title(tlab)    
end
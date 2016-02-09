%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STANDARD DEVIATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = zstd(varargin)
  
  d = 0;
  
  if(nargin == 1 && (~isempty(varargin{1})) )  
    if(sum(~isnan(varargin{1})) > 0)
      d = std(varargin{1}(~isnan(varargin{1})));
    end
  elseif(nargin >= 2 && (~isempty(varargin{1})) && (~isempty(varargin{2}))) 
    
    n = 0;
    s = 0;
        
    x = double(varargin{1});
    y = double(varargin{2});
    
    if(nargin == 4)      
      xlo = varargin{3};
      xhi = varargin{4};
      [rte,ix1]=min(abs(xlo-x));
      [rte,ix2]=min(abs(xhi-x));
      
       x = x(ix1:ix2);
       y = y(ix1:ix2);
    end
    
    x2 = x.^2;    
    s = dot(x,y);
    s2 = dot(x2,y);
    n = sum(y);
    
    if(n > 0)
      d = sqrt(1/(n-1)*(s2 - s*s/n));
    end
    
  end

end
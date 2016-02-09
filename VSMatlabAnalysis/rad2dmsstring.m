%
%   rad2dmsstring
%
%   Convert angle in radians to DMS.
%
%   INPUTS:
%
%       rad: Angle in radians.
%
%   AUTHOR:
%
%       Matthew Wood - mdwood@astro.ucla.edu - 2009-05-11
%
%   VERSION:
%
%       $Id: rad2dmsstring.m,v 3.1 2009/12/16 00:39:27 matthew Exp $
%

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
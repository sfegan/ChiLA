%
%   rad2hmsstring
%
%   Convert angle in radians to HMS.
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
%       $Id: rad2hmsstring.m,v 3.1 2009/12/16 00:39:27 matthew Exp $
%

function str = rad2hmsstring(rad)
    hrs = mod(rad/pi*12+24,24);
    ifracsec = mod(round(hrs*3600*10),24*3600*10);
    str = sprintf('%02d:%02d:%02d.%01d',...
                  floor(ifracsec/3600/10),...
                  floor(mod(ifracsec/60/10,60)),...
                  floor(mod(ifracsec/10,60)),...
                  floor(mod(ifracsec,10)));
end
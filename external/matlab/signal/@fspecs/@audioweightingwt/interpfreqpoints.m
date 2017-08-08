function [Fv Av] = interpfreqpoints(this,Flimits,Alimits,interpFactor,type,varargin) %#ok<INUSL>
%INTERPFREQPOINTS  Perform linear or log-linear interpolation between two
%frequency values

%   Copyright 2009 The MathWorks, Inc.

Fv = linspace(Flimits(1),Flimits(2),interpFactor);

if strcmp(type,'linear')
    if nargin < 6
        Av = interp1(Flimits,Alimits,Fv,'linear');
    else
        p = polyfit(Flimits,Alimits,1);
        Av = polyval(p,varargin{1});
        Fv = varargin{1};
    end
else
    numDecades = log10(Flimits(2)/Flimits(1));
    Ad = (Alimits(2) - Alimits(1))/numDecades;
    
    if nargin < 6
        Av = log10(Fv/Flimits(1))* Ad + Alimits(1);
    else
        Av = log10(varargin{1}/Flimits(1))* Ad + Alimits(1);
        Fv = varargin{1};
    end
end
% [EOF]
function this = fracdelaymeas(hfilter, hfdesign, varargin)
%FRACDELAYMEAS   Construct a FRACDELAYMEAS object.

%   Author(s): V. Pellissier
%   Copyright 2005-2006 The MathWorks, Inc.

error(nargchk(1,inf,nargin,'struct'));

this = fdesign.fracdelaymeas;

% Parse the inputs to get the specification and the measurements list.
minfo = parseconstructorinputs(this, hfilter, hfdesign, varargin{:});

if this.NormalizedFrequency, Fs = 2;
else,                        Fs = this.Fs; end

[gpd w]= grpdelay(hfilter,1024,Fs);

% Measure the fractional delay filter.
if isempty(minfo.Apass),
    this.Fpass1 = findfpass(this, hfilter, minfo.Fpass1, minfo.Apass, 'up',...
        [],{gpd w minfo.FracDelayError 'up' Fs});
    this.Fpass2 = findfpass(this, hfilter, minfo.Fpass2, minfo.Apass, 'down',...
        [],{gpd w minfo.FracDelayError 'down' Fs});
else
    this.Fpass1 = findfpass(this, hfilter, minfo.Fpass1, minfo.Apass, 'up');
    this.Fpass2 = findfpass(this, hfilter, minfo.Fpass2, minfo.Apass, 'down');
end
   
% Measure Apass always after Fpass1 and Fpass2
this.Apass  = measureripple(this, hfilter, minfo.Fpass1, minfo.Fpass2, minfo.Apass);

if this.NormalizedFrequency,
    this.NomGrpDelay = gpd(1)-hfilter.FracDelay;
else
    this.NomGrpDelay = (gpd(1)-hfilter.FracDelay)/Fs;
end
this.FracDelayError = measurefracdelayerr(this, hfilter, minfo.Fpass1, minfo.Fpass2, ...
    minfo.Apass, minfo.FracDelayError);

% [EOF]

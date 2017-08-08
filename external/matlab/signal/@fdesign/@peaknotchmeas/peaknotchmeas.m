function this = peaknotchmeas(hfilter, varargin)
%PEAKNOTCHMEAS   Construct a PEAKNOTCHMEAS object.

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

error(nargchk(1,inf,nargin,'struct'));

this = fdesign.peaknotchmeas;

% Parse the inputs.
minfo = parseconstructorinputs(this, hfilter, varargin{:});

if this.NormalizedFrequency, Fs = 2;
else,                        Fs = this.Fs; end

% Copy known values
f=fieldnames(minfo);

for k = 1:length(f)
    if ~isempty(minfo.(f{k})),
        this.(f{k}) = minfo.(f{k});
    end
end


[H,w] = freqz(hfilter,8192,Fs);

% Get freq. range in which to search
Frange = getrange(this,H,w,minfo);

% Measure the filter.
if strcmpi(class(varargin{1}),'fdesign.notch'),
    this.BWpass = findbwstop(this, hfilter, minfo.BWpass, this.Apass,Frange);
    this.BWstop = findbwstop(this, hfilter,minfo.BWstop, this.Astop,Frange);
else
    this.BWpass = findbwpass(this, hfilter, minfo.BWpass, this.Apass,Frange);
    this.BWstop = findbwpass(this, hfilter,minfo.BWstop, this.Astop,Frange);
end

% Set lower and upper transition widths
if ~isempty(this.BWpass) && ~isempty(this.BWstop),
    [Fpl,Fph] = parameqbandedge(pi*this.F0/(Fs/2),pi*this.BWpass/(Fs/2),0);
    [Fsl,Fsh] = parameqbandedge(pi*this.F0/(Fs/2),pi*this.BWstop/(Fs/2),0);
    this.LowTransitionWidth = abs(Fpl-Fsl)*Fs/(2*pi);
    this.HighTransitionWidth = abs(Fsh-Fph)*Fs/(2*pi);    
end
%%
function BWpass = findbwpass(this, hfilter, BWpass, Apass,Frange)

if isempty(BWpass),
    if isempty(Apass),
        BWpass = [];
    else
        Fph = findfpass(this, hfilter, [], Apass, 'down',Frange);
        Fpl = findfpass(this, hfilter, [], Apass, 'up',Frange);
        BWpass = abs(Fph-Fpl);
    end
end  
%%
function BWstop = findbwstop(this, hfilter, BWstop, Astop,Frange)

if isempty(BWstop),
    if isempty(Astop),
        BWstop = [];
    else
        Fph = findfstop(this, hfilter, [], Astop, 'up',Frange);
        Fpl = findfstop(this, hfilter, [], Astop, 'down',Frange);
        BWstop = abs(Fph-Fpl);
    end
end  


% [EOF]

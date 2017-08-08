function this = parameqmeas(hfilter, varargin)
%PARAMEQMEAS   Construct a PARAMEQMEAS object.

%   Copyright 2008 The MathWorks, Inc.

error(nargchk(1,inf,nargin,'struct'));

% Construct an "empty" object.
this = fdesign.parameqmeas;

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

[H,w] = freqz(hfilter,2^14,Fs);

% Get freq. range in which to search
Frange = getrange(this,H,w,minfo);

fd = varargin{1};
%%
%Special measurements for variable slope shelving filters
if strcmp(fd.Specification,'N,F0,Fc,S,G0') || strcmp(fd.Specification,'N,F0,Fc,Qa,G0') 
    %shelving filters are always referenced to 0 dB and do not have a Gref
    %property. Measurements for these filters are slightly different from
    %the common parametric equalizer measurements.
    Gref = 0; 
    Hm = 20*log10(abs(H));
    F0 = minfo.F0/(Fs/2);
    
    %Measure reference gain, and boost/cut gain
    [this.G0,this.Gref]  = measshelfgogref(Hm,w,F0,fd.Fc,fd.G0,Gref);
          
    %If S>1 the filter has ripple so we can measure ripple related
    %parameters, otherwise, we leave these measurements empty
    if isprop(fd,'Qa')
        Ag=10^(fd.G0/40);
        S = inv((((1/fd.Qa^2)-2)/(Ag+1/Ag) ) + 1);
    else
        S = fd.S;
    end
    if S>1        
        if minfo.F0 %highpass shelf
            [this.BWpass,this.BWstop,this.Gpass,this.Gstop,...
            this.LowTransitionWidth]  = measshelfripple(Hm,w,F0,fd.Fc,Fs,fd.G0,Gref);
            this.HighTransitionWidth  = 0;
        else %lowpass shelf
            this.LowTransitionWidth  = 0;
            [this.BWpass,this.BWstop,this.Gpass,this.Gstop,...
            this.HighTransitionWidth]  = measshelfripple(Hm,w,F0,fd.Fc,Fs,fd.G0,Gref);
        end
    else
        this.BWpass = [];
        this.BWstop = [];
        this.Gpass = [];
        this.Gstop = [];
        this.HighTransitionWidth = [];
        this.LowTransitionWidth = [];
    end
%%   
%Measurements for any other type of filters
else 
    Gref = fd.Gref;
    
    if fd.G0 > Gref,
        this.G0   = 20*log10(max(abs(H)));
        this.Gref = 20*log10(min(abs(H)));
        this.BWpass = findbwpass(this, hfilter, minfo.BWpass, this.G0,minfo.Gpass,Frange);
        this.BWstop = findbwpass(this, hfilter,minfo.BWstop, this.G0,minfo.Gstop,Frange);
    else
        this.Gref = 20*log10(max(abs(H)));
        this.G0   = 20*log10(min(abs(H)));
        this.BWpass = findbwstop(this, hfilter, minfo.BWpass, this.G0,minfo.Gpass,Frange);
        this.BWstop = findbwstop(this, hfilter,minfo.BWstop, this.G0,minfo.Gstop,Frange);
    end
    
    
    % Measure the filter.
    this.Gpass  = findgpass(this, hfilter, minfo.BWpass, minfo.F0, minfo.Gpass,Fs);
    this.Gstop  = findgpass(this, hfilter, minfo.BWstop, minfo.F0, minfo.Gstop,Fs);
    
    % Set lower and upper transition widths
    if ~isempty(this.BWpass) && ~isempty(this.BWstop),
        [Fpl,Fph] = parameqbandedge(pi*this.F0/(Fs/2),pi*this.BWpass/(Fs/2),0);
        [Fsl,Fsh] = parameqbandedge(pi*this.F0/(Fs/2),pi*this.BWstop/(Fs/2),0);
        this.LowTransitionWidth = abs(Fpl-Fsl)*Fs/(2*pi);
        this.HighTransitionWidth = abs(Fsh-Fph)*Fs/(2*pi);
    end
    
end
%Give correct measurements in all-pass cases that arise when G0=Gref
if fd.G0 == Gref    
    this.BW = Fs/2;
    this.BWpass = [];
    this.BWstop = [];
    this.Flow = 0;
    this.Fhigh = Fs/2;
    this.GBW = this.G0;
    this.LowTransitionWidth = [];
    this.HighTransitionWidth = [];
    this.Gpass = [];
    this.Gstop = [];        
end
%%
%Create dynamic private properties depending on whether we have Qa or S specs
g = get(varargin{:});
if isfield(g,'Qa') || isfield(g,'S')
    p = schema.prop(this, 'Qa', 'mxArray');
    if g.FilterOrder == 2        
        A = hfilter.sosMatrix(4:end);
        this.Qa = sqrt(((1+A(3))^2) - (A(2)^2))/(2*(1-A(3)));
        if strcmp(fd.Specification,'N,F0,Qa,Gref,G0')
            %Un-normalize Qa so that the measurement equals the Qa 
            %specified by the user
            this.Qa = this.Qa*10^((this.Gref-this.G0)/40); 
        end
    else
        this.Qa = [];
    end
    set(p, 'AccessFlags.PublicSet', 'Off');    
end

end

%%
function BWpass = findbwpass(this, hfilter, BWpass, G0, Gpass,Frange)

if isempty(BWpass),
    if isempty(Gpass),
        BWpass = [];
    else
        Fpl = findfpass(this, hfilter, [], G0-Gpass, 'down',Frange);
        Fph = findfpass(this, hfilter, [], G0-Gpass, 'up',Frange);
        BWpass = abs(Fph-Fpl);
    end
end  
end
%%
function BWstop = findbwstop(this, hfilter, BWstop, G0,Gstop,Frange)

if isempty(BWstop),
    if isempty(Gstop),
        BWstop = [];
    else
        Fpl = findfstop(this, hfilter, [], -Gstop, 'down',Frange);
        Fph = findfstop(this, hfilter, [], -Gstop, 'up',Frange);
        BWstop = abs(Fph-Fpl);
    end
end  
end
%%
function Gpass  = findgpass(this, hfilter,BWpass,F0,Gpass,Fs)

if isempty(Gpass),
    if isempty(BWpass),
        Gpass = [];
    else
        [Flow,Fhigh] = parameqbandedge(pi*F0/(Fs/2),pi*BWpass/(Fs/2),0);
        Gpass = 20*log10(abs(freqz(hfilter,[Flow,Fhigh])));
        Gpass = Gpass(1); % Both should be the same
    end
end

end

%%
function [G0,Gref]  = measshelfgogref(H,w,F0,Fc,g0,gref)
    idxL = w<(Fc);
    idxH = ~idxL;
    HL = flipud(H(idxL));
    HH = H(idxH);    
    if (F0 && g0>gref) %Highpass boost
        idx1 = find(HL<=gref);               
        idx2 = find(HH >= g0);
        if isempty(idx1), idx1 = length(HL); end
        if isempty(idx2), idx2 = length(HH); end        
        Gref = max(HL(idx1));
        G0 = min(HH(idx2));
    elseif (F0 && g0<gref) %highpass cut        
        idx1 = find(HL>=gref);
        idx2 = find(HH<=g0);
        if isempty(idx1), idx1 = length(HL); end
        if isempty(idx2), idx2 = length(HH); end        
        Gref = min(HL(idx1));
        G0 = max(HH(idx2));     
    elseif (~F0 && g0>gref) %lowpass boost
        idx1 = find(HL>=g0);
        idx2 = find(HH<=gref);  
        if isempty(idx1), idx1 = length(HL); end
        if isempty(idx2), idx2 = length(HH); end        
        G0 = min(HL(idx1));
        Gref = max(HH(idx2)); 
    else %(~F0 &&  g0<gref)lowpass cut
        idx1 = find(HL<=g0);
        idx2 = find(HH >= gref); 
        if isempty(idx1), idx1 = length(HL); end
        if isempty(idx2), idx2 = length(HH); end        
        G0 = max(HL(idx1));
        Gref = min(HH(idx2)); 
    end    
end

%%
function [BWpass,BWstop,Gp,Gst,TB] = measshelfripple(H,w,F0,Fc,Fs,g0,gref)
%MEASSHELFRIPPLE measure the shelving filter parameters
%Input F0 should be normalized so that it can be used as a flag to
%determine if filter design is lowpass or highpass.

%Measure Gp, Gst
    idxL = w<(Fc);
    idxH = ~idxL;
    HL = flipud(H(idxL));
    wL = flipud(w(idxL));
    HH = H(idxH);
    wH = w(idxH);
    
    if (F0 && g0>gref) %Highpass boost
        idx1 = find(HL<=gref);               
        idx2 = find(HH >= g0);
        Gp = max(HH);
        Gst = min(HL); 
    elseif (F0 && g0<gref) %highpass cut        
        idx1 = find(HL>=gref);
        idx2 = find(HH<=g0);
        Gp = min(HH);
        Gst = max(HL); 
    elseif (~F0 && g0>gref) %lowpass boost
        idx1 = find(HL>=g0);
        idx2 = find(HH<=gref);  
        Gp = max(HL);
        Gst = min(HH);  
    else %(~F0 &&  g0<gref)lowpass cut
        idx1 = find(HL<=g0);
        idx2 = find(HH >= gref); 
        Gst = max(HH);
        Gp = min(HL);  
    end  
    
%Measure BWpass, BWstop, transition BW
    if isempty(idx1), idx1 = length(HL); end
    if isempty(idx2), idx2 = length(HH); end
    w1 = wL(idx1(1)); 
    w2 = wH(idx2(1));
    
    if F0
        BWpass = (Fs/2)-w2;
        BWstop = (Fs/2)-w1;
    else
        BWpass = w1;
        BWstop = w2;
    end  
    TB = (w2-w1);
end

% [EOF]

function this = combmeas(hfilter, varargin)
%PEAKNOTCHMEAS   Construct a COMBMEAS object.

%   Copyright 2008 The MathWorks, Inc.

error(nargchk(1,inf,nargin,'struct'));

this = fdesign.combmeas;

% Parse the inputs.
minfo = parseconstructorinputs(this, hfilter, varargin{:});
N = minfo.FilterOrder;

if this.NormalizedFrequency
    Fs = 2;
else
    Fs = this.Fs;
end

[H,w] = freqz(hfilter,32768,Fs);
HdB = 20*log10(abs(H));

%Measure BW
% Get freq. range in which to search (first peak or notch centered at zero)
idx = w < (1/N)*(Fs/2);
wL = w(idx);
HbwL = HdB(idx);
[dummy,idxL] = min(abs(HdB(idx)-minfo.GBW));
idxL = idxL(1);
if idxL > 1
    if isequal(lower(minfo.CombType),'peak')
        if HbwL(idxL) > minfo.GBW
           X = [HbwL(idxL) HbwL(idxL+1)];
           Y = [wL(idxL) wL(idxL+1)];
        else
           X = [HbwL(idxL-1) HbwL(idxL)];
           Y = [wL(idxL-1) wL(idxL)];
        end
    else %notch
        if HbwL(idxL) > minfo.GBW
           X = [HbwL(idxL-1) HbwL(idxL)];
           Y = [wL(idxL-1) wL(idxL)];            
        else
           X = [HbwL(idxL) HbwL(idxL+1)];
           Y = [wL(idxL) wL(idxL+1)];            
        end
    end  
    X(X==-inf) = -3000; %to avoid NaN errors in tostring functions
    p = polyfit(X,Y,1);            
    this.BW = 2*polyval(p,minfo.GBW);
    this.GBW = minfo.GBW;        
else
    this.BW = 2*wL(idxL);
    this.GBW = HbwL(idxL);
end

if N>1
    this.Q = (2/N)*(Fs/2)/this.BW;
else
    this.Q = [];
end

%Measure peak/notch positions and their respective magnitudes
%Find windows for each peak/notch
%Measure frequency position and magnitude value at each window
%Measure reference gain and peak/notch gain
if isequal(lower(minfo.CombType),'notch')    
    this.Gref = 20*log10(max(abs(H)));
    this = measureFreqs(this,HdB,w,N,Fs,'notch');   
else
    this.Gref = 20*log10(min(abs(H)));   
    this = measureFreqs(this,HdB,w,N,Fs,'peak');   
end

end

%-------------------------------------------------------------------------p
%,
function this = measureFreqs(this,HdB,w,N,Fs,ctype)
%Measure position and magnitude of peak or notch frequencies

if ~rem(N,2)
    K = (N/2)+1;
else
    K = (N+1)/2;
end

mGains = zeros(1,K);
mFreqs = zeros(1,K);

for k = 1:K
    idx = (w>=0)&(w>(2*k-3)*(Fs/2)/N) & (w < (2*k-1)*(Fs/2)/N) & (w<=Fs/2);
    W = w(idx);
    Hw = HdB(idx);
    if isequal(ctype,'notch')
        [mag, idx2] = min(Hw);
    else
        [mag, idx2] = max(Hw);
    end
    mGains(k) = mag(1);
    mFreqs(k) = W(idx2(1));
    clear W Hw
end

%Create dynamic private properties depending on whether we have a peak or
%notch comb filter. Set the corresponding properties with the measured
%gains and frequency values. 
if isequal(ctype,'notch')
    p = schema.prop(this, 'Gnotch', 'mxArray');
    this.Gnotch = mGains;
    set(p, 'AccessFlags.PublicSet', 'Off');
    
    p = schema.prop(this, 'NotchFrequencies', 'mxArray');
    this.NotchFrequencies = mFreqs;
    set(p, 'AccessFlags.PublicSet', 'Off');
else
    p = schema.prop(this, 'Gpeak', 'mxArray');
    this.Gpeak = mGains;
    set(p, 'AccessFlags.PublicSet', 'Off');
    
    p = schema.prop(this, 'PeakFrequencies', 'mxArray');
    this.PeakFrequencies = mFreqs;
    set(p, 'AccessFlags.PublicSet', 'Off');
end

    
end


% [EOF]

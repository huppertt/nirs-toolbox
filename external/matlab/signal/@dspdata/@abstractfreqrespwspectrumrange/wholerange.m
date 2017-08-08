function wholerange(this)
%WHOLERANGE   Power spectrum calculated over the whole Nyquist interval.

%   Author(s): P. Pacheco
%   Copyright 1988-2004 The MathWorks, Inc.

newSpectrumRange = 'whole';
if strcmpi(this.SpectrumRange,newSpectrumRange),
    % Already a 'whole' spectrum.
    return;
end

if this.NormalizedFrequency,
    fnyq = pi; 
else
    fnyq = this.getfs/2;
end

% Convert a spectrum calculated over the 'half' the Nyquist interval to
% spectrum calculated over the 'whole' Nquist interval.
Sxx = this.Data;
[Nfft,nchans] = size(Sxx);
startIdx = Nfft+1;
W = this.Frequencies;
if isevenwholenfft(this,Nfft,W),      % EVEN "whole" NFFT
    endIdx = (Nfft-1)*2;
    Sxx(startIdx:endIdx,:) = Sxx(Nfft-1:-1:2,:);  % Add positive half.
    W(startIdx:endIdx) = W(2:Nfft-1)+fnyq;
    
else                                  % ODD "whole" NFFT
    endIdx = (Nfft*2)-1;
    Sxx(startIdx:endIdx,:) = Sxx(Nfft:-1:2,:);    % Add positive half.
    W(startIdx:endIdx) = W(2:Nfft)+fnyq;
end

set(this,...
    'Data',Sxx,...
    'Frequencies',W);

setspectrumtype(this,newSpectrumRange);

% [EOF]

function twosided(this)
%TWOSIDED   Convert a one-sided spectrum to a two-sided spectrum.

%   Copyright 1988-2012 The MathWorks, Inc.

newSpectrumType = 'twosided';
if strcmpi(this.SpectrumType,newSpectrumType),
    return;    % Spectrum already two-sided.
end

if this.NormalizedFrequency,
    fnyq = pi; 
else
    fnyq = this.getfs/2;
end

Pxx = this.Data;
W   = this.Frequencies;
CI  = this.ConfInterval;
[Nfft,~] = size(Pxx);

% Rebuild the 'twosided' PSD from the 'onesided' PSD.
startIdx = Nfft+1;
if isevenwholenfft(this,Nfft,W),      % EVEN "whole" NFFT
    endIdx = (Nfft-1)*2;
    Pxx(2:Nfft-1,:) = Pxx(2:Nfft-1,:)/2;
    Pxx(startIdx:endIdx,:) = Pxx(Nfft-1:-1:2,:);  % Add positive half.
    W(startIdx:endIdx) = W(2:Nfft-1)+fnyq;
    if ~isempty(CI)
        CI(2:Nfft-1,:) = CI(2:Nfft-1,:)/2;
        CI(startIdx:endIdx,:) = CI(Nfft-1:-1:2,:);
    end    
else                                  % ODD "whole" NFFT
    endIdx = (Nfft*2)-1;
    Pxx(2:Nfft,:) = Pxx(2:Nfft,:)/2;
    Pxx(startIdx:endIdx,:) = Pxx(Nfft:-1:2,:);    % Add positive half.
    W(startIdx:endIdx) = W(2:Nfft)+fnyq;
    if ~isempty(CI)
        CI(2:Nfft,:) = CI(2:Nfft,:)/2;
        CI(startIdx:endIdx,:) = CI(Nfft:-1:2,:);    % Add positive half.
    end
end

this.Data = Pxx;
this.Frequencies = W;
this.ConfInterval = CI;
setspectrumtype(this,newSpectrumType); % Uses priv property to produce better error msg.

% [EOF]

function fderr = measurefracdelayerr(this, hfilter, Fstart, Fend, Apass, FracDelayError)
%MEASUREFRACDELAYERR   

%   Author(s): V. Pellissier
%   Copyright 2006 The MathWorks, Inc.

error(nargchk(4,6,nargin,'struct'));

if nargin < 5
    FracDelayError = [];
end
fderr = [];

if isempty(Fstart) || isempty(Fend)
    
    if isempty(Apass) && ~isempty(FracDelayError),
        % If we are missing either the start of the passband or the end, we
        % cannot measure the error.  Return the spec.
        fderr = FracDelayError;
    else
        % Read measured Fpass1 Fpass2 computed from the FracDelayError
        Fstart = this.Fpass1;
        Fend = this.Fpass2;
        if isempty(Fstart) || isempty(Fend)
            fderr = FracDelayError;
        else
            fderr = measurefracdelayerr(this, hfilter, Fstart, Fend);
        end
    end
else
    % Measure fractional delay error between Fstart and Fend
    nomgdp = floor(order(hfilter)/2);
    fd = hfilter.FracDelay;
    idealgpd = nomgdp + fd;
    
    if this.NormalizedFrequency, Fs = 2;
    else,                        Fs = this.Fs; end
    
    % Refine search
    N = 2^10;
    gpd= grpdelay(hfilter, linspace(Fstart, Fend, N), Fs);
    
    fderr = max(abs(gpd-idealgpd));
    
    if ~this.NormalizedFrequency,
        fderr = fderr/Fs; 
    end

end

% [EOF]

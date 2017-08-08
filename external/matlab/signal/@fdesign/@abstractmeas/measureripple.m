function rip = measureripple(this, hfilter, Fstart, Fend, Apass, idealfcn)
%MEASURERIPPLE   Return the ripple in the passband.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

error(nargchk(4,6,nargin,'struct'));

if nargin < 5
    Apass = [];
end

if isempty(Fstart) || isempty(Fend)
    
    % If we are missing either the start of the passband or the end, we
    % cannot measure the ripple.  Return either the spec or [].
    rip = Apass;
else

    N = 2^10;

    if this.NormalizedFrequency, Fs = 2;
    else,                        Fs = this.Fs; end

    [h, w] = freqz(hfilter, linspace(Fstart, Fend, N), Fs);
    h      = abs(h);
    
    if nargin > 5
        if iscell(idealfcn)
            idealh = feval(idealfcn{1}, w/(Fs/2), idealfcn{2:end});
        else
            idealh = feval(idealfcn, w/(Fs/2))
        end
        h = h-idealh+1;
    end
    
    % The ripple is defined as the amplitude (dB) variation between the two
    % specified frequency points.
    rip = db(max(h))-db(min(h));
end

% [EOF]

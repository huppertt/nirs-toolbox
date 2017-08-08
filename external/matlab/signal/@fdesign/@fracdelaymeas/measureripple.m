function rip = measureripple(this, hfilter, Fstart, Fend, Apass)
%MEASURERIPPLE   Return the ripple in the passband.

%   Author(s): J. Schickler
%   Copyright 2005-2006 The MathWorks, Inc.

error(nargchk(4,5,nargin,'struct'));

if nargin < 5
    Apass = [];
end

if isempty(Fstart) || isempty(Fend)
    
    if ~isempty(Apass),
        % If we are missing either the start of the passband or the end, we
        % cannot measure the ripple.  Return the spec.
        rip = Apass;
    else
        % Read measured Fpass1 Fpass2 computed from the FracDelayError
        Fstart = this.Fpass1;
        Fend = this.Fpass2;
        if isempty(Fstart) || isempty(Fend)
            rip = Apass;
        else
            rip = measureripple(this, hfilter, Fstart, Fend);
        end
    end
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
